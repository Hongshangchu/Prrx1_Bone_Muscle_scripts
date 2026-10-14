# ============================================================================
# 00_setup.R —— 环境准备、依赖安装 (R 4.4, Windows/Linux 通用, v2 修正安装源)
# 依赖安装源分类 (关键!):
#   CRAN          : install.packages()
#   Bioconductor  : BiocManager::install()
#   GitHub        : devtools::install_github()  (Windows 需 Rtools)
# ============================================================================

install_pkg <- function(pkg) {
  if (pkg %in% rownames(installed.packages())) return(invisible(TRUE))
  message("Installing (CRAN): ", pkg)
  tryCatch(install.packages(pkg), error = function(e) message("FAILED: ", pkg))
}
bioc_install <- function(pkg) {
  if (pkg %in% rownames(installed.packages())) return(invisible(TRUE))
  message("Installing (Bioconductor): ", pkg)
  tryCatch(BiocManager::install(pkg, ask = FALSE, update = FALSE),
           error = function(e) message("FAILED: ", pkg))
}
gh_install <- function(repo) {
  pkg <- strsplit(repo, "/")[[1]][2]
  if (pkg %in% rownames(installed.packages())) return(invisible(TRUE))
  message("Installing (GitHub): ", repo)
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  tryCatch(devtools::install_github(repo),
           error = function(e) message("FAILED: ", repo, " - ", e$message))
}

## ---- CRAN ----
cran_pkgs <- c("Seurat","harmony","Matrix","tidyverse","patchwork","future",
               "cowplot","RColorBrewer","WGCNA","glmnet","randomForest",
               "xgboost","caret","pROC","msigdbr","monocle3","BisqueRNA","SoupX")
invisible(lapply(cran_pkgs, install_pkg))

## ---- Bioconductor ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("GEOquery","SingleR","celldex","clusterProfiler","org.Mm.eg.db",
               "fgsea","enrichplot","slingshot","tradeSeq","apeglm","celda")
invisible(lapply(bioc_pkgs, bioc_install))

## ---- GitHub (不在 CRAN/Bioc 的包) ----
gh_install("sqjin/CellChat")      # CellChat 仅 GitHub
gh_install("saeyslab/nichenetr")  # nichenetr 仅 GitHub

## ---- 验证 ----
all_pkgs <- c(cran_pkgs, bioc_pkgs, "CellChat", "nichenetr")
missing  <- all_pkgs[!all_pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  warning("以下包未安装成功, 请单独处理:\n  ", paste(missing, collapse = ", "))
} else message("All packages installed successfully.")

## ---- 加载 ----
to_load <- all_pkgs[all_pkgs %in% rownames(installed.packages())]
suppressPackageStartupMessages(lapply(to_load, library, character.only = TRUE))
plan("multisession", workers = 8)
options(future.globals.maxSize = 16 * 1024^3)
## 注意: 对超大对象 (如 02 整合后的 ~27GB 对象), future 并行需序列化
## 整个对象发给工作线程, 必然超限 —— 大对象脚本内已改为 plan("sequential")
set.seed(1234)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)

## ---- 网络代理 (国内用户如 R 下载失败, 取消注释并改为实际端口) ----
# Sys.setenv(https_proxy = "http://127.0.0.1:7890",
#            http_proxy  = "http://127.0.0.1:7890")

## 注意: GEO 数据下载已移至 scripts/00b_download_geo.R
##       辅助数据库下载已移至 scripts/00c_download_resources.R
