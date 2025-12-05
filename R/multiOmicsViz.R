utils::globalVariables(c("x", "y", "value"))
#' Multi-Omics Correlation Visualization Along Chromosomes
#'
#' This function calculates multi-omics correlations and visualizes them 
#' along genomic coordinates using ggplot2. It supports one or multiple 
#' target omics layers and produces integrated heatmap + barplot panels.
#'
#' @param sourceOmics A source omics matrix or `SummarizedExperiment`.
#' @param sourceOmicsName Character string, name of the source omics.
#' @param chrome_sourceOmics Chromosomes to include for source omics.
#' @param targetOmicsList A list containing one or multiple target omics matrices
#'   or `SummarizedExperiment` objects.
#' @param targetOmicsName Character vector naming each element in `targetOmicsList`.
#' @param chrome_targetOmics Chromosomes to include for target omics.
#' @param fdrThr Numeric FDR threshold.
#' @param outputfile File name prefix for saving the PNG output.
#' @param nThreads Integer number of parallel threads (required if multiple 
#'   target omics are supplied).
#' @param legend Logical, whether to show the color legend.
#' @param point_size Numeric, point size in heatmaps.
#' @param axis_label_style Character, one of `"plain"` or `"zigzag"`, controlling
#'   chromosome axis label spacing.
#'
#' @return A `ggplot` object (single omics) or a list of `ggplot` objects 
#'   (multi-omics).
#'
#' @import ggplot2
#' @import ggpubr
#' @import patchwork
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom SummarizedExperiment assays
#' @importFrom stats cor p.adjust pt
#' @importFrom grDevices png dev.off
#'
#' @export
#'
#' @examples
#' \dontrun{
#' multiOmicsViz(
#'   sourceOmics = cna_matrix,
#'   sourceOmicsName = "CNA",
#'   chrome_sourceOmics = "1",
#'   targetOmicsList = list(Proteome = prot_matrix),
#'   targetOmicsName = "Proteome",
#'   chrome_targetOmics = "1",
#'   fdrThr = 0.05,
#'   outputfile = "CNA_Protein"
#' )
#' }


multiOmicsViz <- function(
    sourceOmics,
    sourceOmicsName,
    chrome_sourceOmics,
    targetOmicsList,
    targetOmicsName,
    chrome_targetOmics,
    fdrThr,
    outputfile,
    nThreads = NULL,
    legend = TRUE,
    point_size = 0.5,             # 新增参数
    axis_label_style = "zigzag"   # 新增参数
){
    
    outputfile <- paste(outputfile,".png",sep="")
    
    if(inherits(sourceOmics, "SummarizedExperiment")){     #第一处
      sourceOmics <- assays(sourceOmics, n=1)
    }else{
      if(!is.matrix(sourceOmics) && !is.data.frame(sourceOmics)){
        stop("Source omics data (e.g. CNA data) should be a R matrix, 
        data.frame or SummarizedExperiment object.")
      }
    }
    
    if(!is.character(sourceOmicsName)){
      stop("sourceOmicsName is the name of the source omics data, which 
      should be a R character object.") 
    }
    
    if(!is.list(targetOmicsList)){
      stop("Multipl target omics data (e.g. mRNA or protein data) 
      should be saved in the list object.")
    }
    
    if(length(targetOmicsList)>5){
      stop("targetOmicsList can only contain at most five omics data.")
    }
    
    for(i in seq_len(length(targetOmicsList))){
      if(inherits(targetOmicsList[[i]], "SummarizedExperiment")){
        targetOmicsList[[i]] <- assays(targetOmicsList[[i]], 1)
      } else {
        if(!is.matrix(targetOmicsList[[i]]) && 
           !is.data.frame(targetOmicsList[[i]])){
          stop("Each of all target omics data in the list (e.g. mRNA or 
           protein data) should be a R matrix, data.frame or 
           SummarizedExperiment object.")
        }
      }
    }  #修改4
    
    
    if(!is.character(targetOmicsName)){
      stop("targetOmicsName should be a R vector object 
      containing the name for each omics data in the targetOmicsList.")
    }
    
    if(length(targetOmicsList)!=length(targetOmicsName)){
      stop("targetOmicsName should have the same length 
      with targetOmicsList.")
    }   
    
    if(length(targetOmicsList)>1 && is.null(nThreads)){
      stop("Please input nThreads for the parallel computing.")
    }
    
    if(length(targetOmicsList)>1 && !is.null(nThreads)){
    
      if(nThreads>length(targetOmicsList)){
        stop("nThreads should be at most the length of targetOmicsList.")
      }
    }
    
        
    ########find the intersect genes among all omics data######
    intG <- c()
    for(i in seq_len(length(targetOmicsList))){
      if(i==1){
        intG <- rownames(targetOmicsList[[i]])
      }else{
        intG <- intersect(intG,rownames(targetOmicsList[[i]]))
      }
    }
    
    if(length(intG)==0){
      stop("The ID types of all omics data in the targetOmicsList 
      should be the same.")
    }
    
    #####process chrome location###
    
    chromeList <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
    "14","15","16","17","18","19","20","21","22","X","Y","All")
    
    x <- setdiff(chrome_sourceOmics,chromeList)
    if(length(x)>0){
      stop('The input chrome infomation for source omics data contains the 
      invalid information. Please only input the chromosome names from 
      the following list: "1","2","3","4","5","6","7","8","9","10","11","12",
      "13","14","15","16","17","18","19","20","21","22","X","Y" and "All".')
    }
    
    x <- setdiff(chrome_targetOmics,chromeList)
    if(length(x)>0){
      stop('The input chrome infomation for target omics data contains the 
      invalid information. Please only input the chromosome names from 
      the following list: "1","2","3","4","5","6","7","8","9","10","11","12",
      "13","14","15","16","17","18","19","20","21","22","X","Y" and "All".')
    }

    if((length(chrome_sourceOmics)>1 && which(chrome_sourceOmics=="All")>1) 
    || (length(chrome_sourceOmics)==1 && chrome_sourceOmics=="All")){
      chrome_sourceOmics <- "All"
    }
    
    if((length(chrome_targetOmics)>1 && which(chrome_targetOmics=="All")>1) 
    || (length(chrome_targetOmics)==1 && chrome_targetOmics=="All")){
      chrome_targetOmics <- "All"
    }
    
    if(chrome_sourceOmics=="All"){
      chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
      "12","13","14","15","16","17","18","19","20","21","22","X","Y")
    }
      
    if(chrome_targetOmics=="All"){
      chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
      "12","13","14","15","16","17","18","19","20","21","22","X","Y")
    }
      
    #######Extract sub list#########
    genelocate <- load_genelocate()
    
    genelocate_sourceOmics <- genelocate[genelocate[,2] %in%
chrome_sourceOmics,]
    genelocate_targetOmics <- genelocate[genelocate[,2] %in%
chrome_targetOmics,]
  
    intG <- intersect(intG,genelocate_targetOmics[,1])
    
    if(length(intG)==0){
      stop("The ID types for all omics data in the targetOmicsList should be 
      gene symbol or all genes in the target omics data are not in the 
      selected chromosomal location chrome_targetOmics.")
    }
    
    for(i in seq_len(length(targetOmicsList))){
      targetOmicsList[[i]] <- targetOmicsList[[i]][intG,]
    }
        
    source_gene <- rownames(sourceOmics)
    source_gene_locate <-
intersect(unique(genelocate_sourceOmics[,1]),source_gene)
    if(length(source_gene_locate)==0){
      stop("The ID type in the source omics data should be gene symbol or all 
      genes in the source omics data are not in the selected chromosomal 
      location chrome_sourceOmics.")
    }
    source_gene <- sourceOmics[source_gene_locate,]
    
    genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[,1] 
    %in% source_gene_locate,]
    genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[,1] 
    %in% intG,]
  
    
    ###Calculate the correlation between cna and other omics data######
    cat("Identify the significant correlations...\n")
    if(length(targetOmicsList)==1){
      resultList <-
calculateCorForTwoMatrices(source_gene,targetOmicsList[[1]],fdrThr)
    }else{
      cl <- parallel::makeCluster(nThreads)
      doParallel::registerDoParallel(cl)
      
      resultList <- foreach::foreach(
        i = seq_len(length(targetOmicsList)),
        .packages="multiOmicsVizEnhanced"
      ) %dopar% {
        calculateCorForTwoMatrices(source_gene, targetOmicsList[[i]], fdrThr)
      }
      
      parallel::stopCluster(cl)
    }
        
    
    ##Calculate the location of genes in the heatmap
    chromLength <- load_chromLength()
    
    re <-
.calculateChromLength(chromLength,chrome_sourceOmics,genelocate_sourceOmics)
    genelocate_sourceOmics <- re$genelocate
    chromLength_sourceOmics <- re$chromLength
    
    re <-
.calculateChromLength(chromLength,chrome_targetOmics,genelocate_targetOmics)
    genelocate_targetOmics <- re$genelocate
    chromLength_targetOmics <- re$chromLength
    
    ##########Plot Figure############
    
    cat("Plot figure...\n")
      

    # 单组学模式
    if (length(targetOmicsList) == 1) {
      
      layout <- autoLayout(nrow(genelocate_sourceOmics), 1)
      
      # 1) heatmap
      p <- .plotHeatMap(
        resultList,
        genelocate_sourceOmics,
        chromLength_sourceOmics,
        genelocate_targetOmics,
        chromLength_targetOmics,
        sourceOmicsName,
        targetOmicsName,
        point_size       = point_size,
        axis_label_style = axis_label_style
      )
      if (!legend) p <- p + theme(legend.position = "none")
      
      # 2) barplot
      bar_p <- .plotSummaryBar_singleOmics(
        resultList,
        chromLength_sourceOmics,
        genelocate_sourceOmics,
        sourceOmicsName,
        axis_label_style = axis_label_style
      )
      
      # 3) assemble
      final_plot <- p / bar_p +
        patchwork::plot_layout(heights = c(3, 1)) &
        ggplot2::theme(
          plot.margin   = ggplot2::margin(4, 4, 4, 4),
          panel.spacing = grid::unit(1, "mm")
        )
      
      png(outputfile, width = layout$width, height = layout$height, res = 600)
      print(final_plot)
      dev.off()
      
      # 返回 ggplot
      return(final_plot)
    }else {
      
      # 1) Heatmaps
      heat_list <- .plotHeatMap_multi(
        resultList,
        genelocate_sourceOmics,
        chromLength_sourceOmics,
        genelocate_targetOmics,
        chromLength_targetOmics,
        sourceOmicsName,
        targetOmicsName,
        point_size       = point_size,
        axis_label_style = axis_label_style
      )
      
      if (!legend) {
        heat_list <- lapply(heat_list, function(pp) {
          pp + theme(legend.position = "none")
        })
      }
      
      # 2) Barplots
      bar_list <- .plotSummaryBar_multi(
        resultList,
        chromLength_sourceOmics,
        genelocate_sourceOmics,
        sourceOmicsName,
        axis_label_style = axis_label_style
      )
      
      # 3) Combine per panel
      columns <- lapply(seq_along(heat_list), function(i){
        heat_list[[i]] / bar_list[[i]] +
          patchwork::plot_layout(heights = c(3, 1))
      })
      
      
      final_plot <- patchwork::wrap_plots(columns, ncol = length(columns)) &
        ggplot2::theme(
          plot.margin   = ggplot2::margin(4, 4, 4, 4),
          panel.spacing = grid::unit(1, "mm")
        )
      
      layout <- autoLayout(
        n_genes  = nrow(genelocate_sourceOmics),
        n_panels = length(heat_list)
      )
      
      png(outputfile, width = layout$width, height = layout$height, res = 600)
      print(final_plot)
      dev.off()
      
      # ⭐⭐⭐ 关键：返回所有 ggplot 对象
      names(columns) <- targetOmicsName
      return(columns)
    }
    
    
    
}

#' @keywords internal
.calculateChromLength <- function(chromLength,selectedChrom,genelocate){
    chromLength <- chromLength[chromLength[,1] %in% selectedChrom,,drop=FALSE]
    
    if(length(selectedChrom)==1){
      x <- 0
    }else{
      x <- c(0,chromLength[1:(nrow(chromLength)-1),2])
    }
    chromLength[,3] <- cumsum(as.numeric(x))
    chromLength[,4] <- cumsum(as.numeric(chromLength[,2]))
    
    genelocate <- cbind(genelocate,0,0)
    
    colnames(genelocate)[5:6] <- c("finalstart","finalend")
        
    for(i in c(1:nrow(genelocate))){
        chr <- genelocate[i,2]
        s <- genelocate[i,3]
        e <- genelocate[i,4]
        cs <- chromLength[chromLength[,1]==chr,3]
        genelocate[i,5] <- s+cs
        genelocate[i,6] <- e+cs
    }
    re <- list(chromLength=chromLength,genelocate=genelocate)
    return(re)
}



#统一的刻度排布-----------------------------------------------------
#' @keywords internal
.formatChromLabels <- function(chrom_labels, style = "plain"){
  
  style <- match.arg(style, c("plain", "zigzag"))
  
  # ─────────────────────────────
  # 1) plain → 不做任何换行
  # ─────────────────────────────
  if(style == "plain"){
    return(chrom_labels)
  }
  
  # ─────────────────────────────
  # 2) zigzag → 9–Y 交错上下分布
  #     上排：9, 11, 13, 15, 17, 19, 21
  #     下排：10, 12, 14, 16, 18, 20, 22, X, Y
  # ─────────────────────────────
  labels2 <- chrom_labels
  
  odd_chr  <- c("9","11","13","15","17","19","21","X")
  even_chr <- c("10","12","14","16","18","20","22","Y")
  
  # 偶数与 X/Y 换到第二行（不挤）
  labels2[chrom_labels %in% even_chr] <-
    paste0("\n", chrom_labels[chrom_labels %in% even_chr])
  
  return(labels2)
}

#' @keywords internal
.formatChromLabelsY <- function(chrom_labels, style = "zigzag"){
  
  style <- match.arg(style, c("plain", "zigzag"))
  if (style == "plain") return(chrom_labels)
  
  odd_chr  <- c("1","3","5","7","9","11","13","15","17","19","21","X")
  even_chr <- c("2","4","6","8","10","12","14","16","18","20","22","Y")
  
  labels2 <- chrom_labels
  
  # EM SPACE：宽度远大于普通空格，旋转后仍能产生横向偏移
  offset <- "\u2003\u2003"   # 两个 EM SPACE ＝ 一个汉字的宽度
  
  labels2[chrom_labels %in% even_chr] <-
    paste0(offset, chrom_labels[chrom_labels %in% even_chr])
  
  return(labels2)
}



# 全局自动布局函数：自动计算 Heatmap 与 Barplot 的宽与高-----------------------------
#' @keywords internal
autoLayout <- function(
    n_genes,
    n_panels = 1,
    base_side = 3200,   # 每列基础宽度（heatmap 这一块），你可以改大一点
    gene_unit = 2500
){
  # 基因多的话可以稍微放大一点；也可以直接写 scale_factor <- 1 固定大小
  scale_factor <- max(1, sqrt(n_genes / gene_unit))
  
  side <- base_side * scale_factor  # 每一列中 heatmap 的“目标边长（像素）”
  
  width  <- side * n_panels         # 有几列就乘几列
  height <- side * 4/3              # 4/3 * side → 上面 3/4 是 heatmap，下面 1/4 是 barplot
  
  list(width = width, height = height)
}




#画图函数-----------------------------------------------------
#' @keywords internal
.plotHeatMap <- function(
    corrArray,
    genelocate_sourceOmics,
    chromLength_sourceOmics,
    genelocate_targetOmics,
    chromLength_targetOmics,
    sourceOmicsName,
    targetOmicsName,
    point_size = 0.5,
    axis_label_style = "zigzag"   # ← X 轴仍然 zigzag
){
 
  
  # ---------- 1 构建数据 ----------
  idx <- which(corrArray != 0, arr.ind = TRUE)
  
  df <- data.frame(
    cna_gene  = rownames(corrArray)[idx[,1]],
    targ_gene = colnames(corrArray)[idx[,2]],
    corr      = corrArray[idx]
  )
  
  df$x <- genelocate_sourceOmics$finalstart[
    match(df$cna_gene, genelocate_sourceOmics[,1])
  ]
  df$y <- genelocate_targetOmics$finalstart[
    match(df$targ_gene, genelocate_targetOmics[,1])
  ]
  
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  df$sign <- ifelse(df$corr > 0, "pos", "neg")
  
  col_map <- c(pos="#d62728", neg="#1f77b4")
  
  # ---------- 2 breaks ----------
  x_breaks <- chromLength_sourceOmics[,4] - chromLength_sourceOmics[,2]/2
  y_breaks <- chromLength_targetOmics[,4] - chromLength_targetOmics[,2]/2
  
  # X 轴保持 zigzag，Y 轴保持正常
  x_labels <- .formatChromLabels(chromLength_sourceOmics[,1], axis_label_style)
  y_labels <- chromLength_targetOmics[,1]  # ← 恢复正常 Y 轴
  
  # ---------- 3 主图 ----------
  p <- ggplot2::ggplot(df, ggplot2::aes(x=x, y=y, color=sign)) +
    ggplot2::geom_point(size=point_size, alpha=0.9) +
    ggplot2::scale_color_manual(values = col_map) +
    
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0,0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = c(0,0)
    ) +
    
    ggplot2::geom_vline(xintercept = chromLength_sourceOmics[,4], color="grey85", size=0.3) +
    ggplot2::geom_hline(yintercept = chromLength_targetOmics[,4], color="grey85", size=0.3) +
    
    ggplot2::labs(
      x = "",
      y = paste(targetOmicsName, "chromosomal position")
    ) +
    
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size=9),
      axis.text.y = ggplot2::element_text(size=7),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}

  

#' @keywords internal
.plotHeatMap_multi <- function(
    resultList,
    genelocate_sourceOmics,
    chromLength_sourceOmics,
    genelocate_targetOmics,
    chromLength_targetOmics,
    sourceOmicsName,
    targetOmicsName,
    point_size       = 0.5,
    axis_label_style = "zigzag"
){
  # 返回的多个热图
  plot_list <- list()
  
  # 多组学：resultList 的每个元素 = 一个 target omics 的相关性矩阵
  for (i in seq_along(resultList)) {
    
    corr_mat <- resultList[[i]]
    this_target <- targetOmicsName[i]
    
    # -------------------------------
    # 构成单张 heatmap 图：调用 .plotHeatMap
    # -------------------------------
    p <- .plotHeatMap(
      corrArray                 = corr_mat,
      genelocate_sourceOmics    = genelocate_sourceOmics,
      chromLength_sourceOmics   = chromLength_sourceOmics,
      genelocate_targetOmics    = genelocate_targetOmics,
      chromLength_targetOmics   = chromLength_targetOmics,
      sourceOmicsName           = sourceOmicsName,
      targetOmicsName           = this_target,
      point_size                = point_size,
      axis_label_style          = axis_label_style
    )
    
    # 保存到 list
    plot_list[[i]] <- p
  }
  
  names(plot_list) <- targetOmicsName
  return(plot_list)
}




#' @keywords internal
.plotSummaryBar_singleOmics <- function(
    corrArray,
    chromLength_sourceOmics,
    genelocate_sourceOmics,
    sourceOmicsName,
    axis_label_style = "zigzag"
){

  
  spe <- apply(abs(corrArray), 1, sum)
  
  df <- data.frame(
    gene = names(spe),
    value = spe,
    x = genelocate_sourceOmics$finalstart[
      match(names(spe), genelocate_sourceOmics[,1])
    ]
  )
  df <- df[!is.na(df$x), ]
  
  xmax <- chromLength_sourceOmics[nrow(chromLength_sourceOmics), 4]
  
  chrom_labels <- .formatChromLabels(chromLength_sourceOmics[,1], axis_label_style)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value)) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0), color = "blue", size = 0.3) +
    ggplot2::geom_vline(xintercept = chromLength_sourceOmics[,4], color = "grey80", size = 0.3) +
    ggplot2::scale_x_continuous(
      limits = c(0, xmax),
      breaks = chromLength_sourceOmics[,4] - chromLength_sourceOmics[,2]/2,
      labels = chrom_labels,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = paste(sourceOmicsName, "chromosomal location"),
      y = "Number of significant correlations"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    )
  return(p)
}



#' @keywords internal
.plotSummaryBar_multi <- function(
    resultList,
    chromLength_sourceOmics,
    genelocate_sourceOmics,
    sourceOmicsName,
    axis_label_style = "zigzag"
){
  out <- list()
  
  for(i in seq_along(resultList)){
    out[[i]] <- .plotSummaryBar_singleOmics(
      corrArray = resultList[[i]],
      chromLength_sourceOmics = chromLength_sourceOmics,
      genelocate_sourceOmics = genelocate_sourceOmics,
      sourceOmicsName = sourceOmicsName,
      axis_label_style = axis_label_style
    )
  }
  
  names(out) <- names(resultList)
  return(out)
}

