### **Multi-Omics Correlation Visualization Along Chromosomes Using ggplot2**

**multiOmicsVizEnhanced** 提供一个基于 **ggplot2** 的可视化框架，用于展示源组学数据（如 CNA）与多个目标组学数据（如 Proteomics、RNA-seq）沿染色体位置的相关性。

相较于原始版本 **multiOmicsViz**，本包提供：

- **更高质量的 ggplot2 图形输出**
- **兼容 DEP2 / SummarizedExperiment**
- **可自定义美学：点大小、坐标轴风格、图例等**
- **支持 1–5 层组学同时绘制**
- **自动多面板布局**
- **并行计算加速 Spearman 相关系数计算**

---

## **📘 Input Data Requirements**

### **1️⃣ Source Omics（例如 CNA）**
必须是以下之一：
- matrix
- data.frame
- SummarizedExperiment

结构要求：
- rownames：基因符号  
- colnames：样本 ID  

---

### **2️⃣ Target Omics List**
传入一个包含 **1–5 个组学矩阵**的列表，例如：

list(
  Protein = proteome_matrix,
  Phospho = phospho_matrix
)

每个元素必须满足：
- rownames 为基因符号  
- colnames 与 source omics 完全一致  
- 不同 omics 之间样本顺序一致  

---

## **🧬 Chromosome Annotation Files**

以下基因注释文件已内置于包中（inst/extdata/）：

- **genelocate.txt**：基因所在染色体、起始位置、终止位置  
- **chromLength.txt**：各染色体长度（用于绘图定位）

这些文件将由内部函数自动加载，无需用户手动指定。

---

## **📦 Installation**

### From CRAN
~~install.packages("multiOmicsVizEnhanced")~~  还没上传

### From GitHub (development version)
remotes::install_github("Guijiajun64/multiOmicsVizEnhanced")

---

## **🔍 Example Usage**

result <- multiOmicsViz(  sourceOmics        = cna_matrix,
                          
                          sourceOmicsName    = "CNA",
                          
                          chrome_sourceOmics = "All",
                        
                          targetOmicsList    = list(
                            Protein = protein_matrix,
                            Phospho = phospho_matrix
                          ),
                          
                          targetOmicsName    = c("Protein", "Phospho"),
                          
                          chrome_targetOmics = "All",
                        
                          fdrThr       = 0.05,
                          
                          outputfile   = "CNA_multiOmics",
                          
                          nThreads     = 2,
                          
                          legend       = TRUE,
                          
                          point_size   = 0.5
                       )

---

## **📌 Why Use multiOmicsVizEnhanced?**

- 专为 **多组学整合研究（CNA–表达–蛋白）** 设计  
- 基于染色体的可视化方式清晰展示跨组学一致性与差异  
- 输出图形适合 **发表文章（publication-ready）**  
- 兼容分析流程（如 DEP2）  
- 可灵活定制绘图参数以适应各种科研需求  

---

## **📚 Citation**

Gui J. *multiOmicsVizEnhanced: Multi-Omics Correlation Visualization Along Chromosomes Using ggplot2*. R package version 0.1.0.

---

## **🐞 Bug Reports & Development**

欢迎提交 issue 或功能请求：  
https://github.com/Guijiajun64/multiOmicsVizEnhanced

-------------------------------------------------------------------------------------------------------------------------------------

Multi-Omics Correlation Visualization Along Chromosomes Using ggplot2

multiOmicsVizEnhanced provides a ggplot2-based visualization framework for displaying correlations between a source omics layer (e.g., CNA) and multiple target omics layers (e.g., proteomics, RNA-seq) along chromosomal positions.

Compared with the original multiOmicsViz, this package offers:

Higher-quality ggplot2 visual outputs

Compatibility with DEP2 / SummarizedExperiment

Customizable aesthetics: point size, axis label style, legend control, etc.

Support for plotting 1–5 omics layers simultaneously

Automatic multi-panel layout

Parallelized Spearman correlation computation

## 📘 **Input Data Requirements**
### 1️⃣ **Source Omics (e.g., CNA)**

Accepted formats:

matrix

data.frame

SummarizedExperiment

Requirements:

rownames: gene symbols

colnames: sample IDs

### 2️⃣ **Target Omics List**

Provide a list containing 1–5 omics matrices, for example:

list(
  Protein = proteome_matrix,
  Phospho = phospho_matrix
)


Each matrix must satisfy:

rownames are gene symbols

colnames are identical to the source omics

sample order is consistent across all omics layers

## 🧬 **Chromosome Annotation Files**

The package includes the following annotation files under inst/extdata/:

genelocate.txt — gene symbol, chromosome, start, end

chromLength.txt — chromosome lengths (for coordinate scaling)

These files are automatically loaded internally; no user input is required.

## 📦 **Installation**
From CRAN

install.packages("multiOmicsVizEnhanced") (Not yet released)

From GitHub (development version)
remotes::install_github("Guijiajun64/multiOmicsVizEnhanced")



```r
result <- multiOmicsViz(
  sourceOmics        = cna_matrix,
  sourceOmicsName    = "CNA",
  chrome_sourceOmics = "All",

  targetOmicsList    = list(
    Protein = protein_matrix,
    Phospho = phospho_matrix
  ),

  targetOmicsName    = c("Protein", "Phospho"),
  chrome_targetOmics = "All",

  fdrThr       = 0.05,
  outputfile   = "CNA_multiOmics",
  nThreads     = 2,
  legend       = TRUE,
  point_size   = 0.5
)
```


## 📌 **Why Use multiOmicsVizEnhanced?**

Designed specifically for multi-omics integration (CNA → expression → protein)

Chromosome-based visualization clearly highlights cross-omics concordance and divergence

Produces publication-ready figures suitable for manuscripts

Compatible with existing analysis pipelines (e.g., DEP2)

Highly customizable aesthetics for flexible scientific visualization

## 📚 **Citation**

Gui J. multiOmicsVizEnhanced: Multi-Omics Correlation Visualization Along Chromosomes Using ggplot2. R package version 0.1.0.

## 🐞 **Bug Reports & Development**

Issues and feature requests are welcome at:
https://github.com/Guijiajun64/multiOmicsVizEnhanced
