# ============================================================================
# 01b_decontx.R —— DecontX 环境 RNA 校正 (在 01 之后、02 之前运行)
# 背景: 损伤组织环境 RNA 污染极强 (可扩散转录本如 Ptprc/S100a8/Lyz2),
#       不校正会污染聚类本身, 必须在聚类前去除
# 方法: DecontX (celda, Bioconductor) —— 仅需 表达矩阵+聚类标签,
#       用贝叶斯模型估计并去除每个细胞的环境 RNA; 不依赖空液滴 (适合本数据)
# 输出: results/decontx_objects.rds (02 脚本检测到该文件会优先使用)
# ============================================================================
library(celda); library(Seurat)

## ---- 样本配置 (与 01_load_QC.R 对应) ----
samples <- list(
  list(dir = "data/GSE195940/raw/ss/", project = "Periosteum_SS",
       meta = list(tissue="periosteum", time="d0", injury="none", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE195940/raw/d3/", project = "Periosteum_d3",
       meta = list(tissue="periosteum", time="d3", injury="fracture", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE164573/raw/d0/",  project = "Muscle_d0_none",
       meta = list(tissue="muscle", time="d0", injury="none", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE164573/raw/fd3/", project = "Muscle_d3_fracture",
       meta = list(tissue="muscle", time="d3", injury="fracture", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE164573/raw/fd5/", project = "Muscle_d5_fracture",
       meta = list(tissue="muscle", time="d5", injury="fracture", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE164573/raw/pd3/", project = "Muscle_d3_polytrauma",
       meta = list(tissue="muscle", time="d3", injury="polytrauma", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE164573/raw/pd5/", project = "Muscle_d5_polytrauma",
       meta = list(tissue="muscle", time="d5", injury="polytrauma", age="adult", cre="Prx1Cre")),
  list(dir = "data/GSE138689/raw/ss/", project = "BM_SS",
       meta = list(tissue="BM", time="d0", injury="none", age="adult", cre="LeprCre")),
  list(dir = "data/GSE138689/raw/d7/", project = "BM_d7",
       meta = list(tissue="BM", time="d7", injury="fracture", age="adult", cre="LeprCre"))
)

find_dir <- function(d) {   # 兼容 raw/<cond>/ 与 raw/<cond>/<子目录>/ 两种结构
  fs <- list.files(d, full.names = TRUE)
  sub <- fs[file.info(fs)$isdir]
  if (length(sub) > 0) sub[1] else d
}

## ---- 单样本 DecontX 校正 ----
decontx_one <- function(dir = NULL, h5 = NULL, project = "X", meta = list()) {
  message("==== DecontX: ", project, " ====")
  counts <- if (!is.null(h5)) {
    m <- Read10X_h5(h5)
    m[, colSums(m > 0) >= 200, drop = FALSE]          # 粗定细胞
  } else {
    Read10X(find_dir(dir))
  }
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]

  ## 快速聚类 (供 DecontX 的 z 参数)
  tmp <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  tmp <- NormalizeData(tmp, verbose = FALSE) |>
         FindVariableFeatures(verbose = FALSE) |>
         ScaleData(verbose = FALSE) |> RunPCA(verbose = FALSE)
  tmp <- FindNeighbors(tmp, dims = 1:20, verbose = FALSE) |>
         FindClusters(resolution = 0.3, verbose = FALSE)

  ## DecontX: 估计并去除环境 RNA
  dx <- decontX(x = counts, z = as.integer(tmp$seurat_clusters),
                verbose = FALSE)
  message(project, " 估计污染比例中位数: ",
          round(stats::median(dx$contamination), 3))

  corrected <- round(dx$decontXcounts)                   # 取整为 counts
  corrected[is.na(corrected)] <- 0

  obj <- CreateSeuratObject(corrected, project = project,
                            min.cells = 3, min.features = 200)
  obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
  for (k in names(meta)) obj[[k]] <- meta[[k]]
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                          percent.mt < 15 & nCount_RNA < 50000)
  obj
}

objs <- lapply(samples, function(s) decontx_one(dir = s$dir, h5 = s$h5,
                                                project = s$project, meta = s$meta))
stopifnot(length(objs) == 9)    # 9 个 MTX 样本
saveRDS(objs, "results/decontx_objects.rds")
message("完成: results/decontx_objects.rds (9 个样本) -> 重新运行 02_integrate_harmony.R")
