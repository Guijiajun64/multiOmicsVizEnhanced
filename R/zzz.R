# 新的 datacache 环境
datacache <- new.env(parent = emptyenv())

# 懒加载函数：在真正需要时才读文件，而不是在加载包时读
load_genelocate <- function() {
  if (!exists("genelocate", envir = datacache)) {
    path <- system.file("extdata", "genelocate.txt",
                        package = "multiOmicsVizEnhanced")
    datacache$genelocate <- utils::read.table(
      path, header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
  }
  return(datacache$genelocate)
}

load_chromLength <- function() {
  if (!exists("chromLength", envir = datacache)) {
    path <- system.file("extdata", "chromLength.txt",
                        package = "multiOmicsVizEnhanced")
    datacache$chromLength <- utils::read.table(
      path, header = FALSE, sep = "\t", stringsAsFactors = FALSE
    )
  }
  return(datacache$chromLength)
}

utils::globalVariables(c("x", "y", "value"))
