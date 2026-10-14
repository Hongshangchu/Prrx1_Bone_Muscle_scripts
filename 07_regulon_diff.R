# ============================================================================
# 07_regulon_diff.R —— Aim2c: Regulon AUC 组间差异 (pySCENIC v0.12 loom 适配)
# 读取: col_attrs/RegulonsAUC (复合属性: 每细胞一行, 每 regulon 一列)
# 比较: fracture vs polytrauma (loom 含 none/fracture/polytrauma 三组,
#       必须先过滤成两组, 否则 wilcox.test 报"分组因子只能有两层")
# 前置: 06a 不过滤 injury + 06_scenic.py 重建 loom/aucell
# ============================================================================
library(hdf5r); library(tidyverse); library(Seurat)

h5 <- H5File$new("results/scenic_aucell.loom", mode = "r")
cell_ids <- h5[["col_attrs/CellID"]]$read()
auc_df   <- h5[["col_attrs/RegulonsAUC"]]$read()   # data.frame: 细胞 x regulon
h5$close_all()

stopifnot(nrow(auc_df) == length(cell_ids))
AUC <- as.matrix(auc_df)
rownames(AUC) <- cell_ids
colnames(AUC) <- sub("\\(.*\\)$", "", colnames(AUC))   # 去 "(+)" 后缀

obj <- readRDS("results/integrated_atlas_scored.rds")
common <- intersect(rownames(AUC), colnames(obj))
AUC  <- AUC[common, , drop = FALSE]
meta <- obj@meta.data[common, ]
cat("匹配细胞:", length(common), "/", nrow(AUC), "\n")

## 关键: 先过滤出 fracture + polytrauma 两组再比较 (loom 含三组)
use <- meta$injury %in% c("fracture", "polytrauma")
x_all <- AUC[use, , drop = FALSE]
g_all <- droplevels(factor(meta$injury[use]))

res <- lapply(colnames(AUC), function(r) {
  x <- x_all[, r]
  f <- suppressWarnings(wilcox.test(x ~ g_all))
  tibble(regulon = r, p = f$p.value,
         mean_frac = round(mean(x[g_all == "fracture"]), 4),
         mean_poly = round(mean(x[g_all == "polytrauma"]), 4))
}) |> bind_rows() |> mutate(padj = p.adjust(p, "BH"))

arrange(res, padj) |>
  write.csv("results/tables/regulon_diff.csv", row.names = FALSE)
slice_head(arrange(res, padj), n = 20) |> print()
message("07 完成 -> results/tables/regulon_diff.csv")
# 预期: Stat3 / Nfkb1 / Hif1a 的 mean_poly > mean_frac 且 padj < 0.05
