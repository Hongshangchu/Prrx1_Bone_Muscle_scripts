# ============================================================================
# 07b_per_celltype_regulon.R —— 逐细胞类型 regulon 活性差异 (polytrauma vs fracture)
# 目的: 分解混合池比较中的"组成效应 vs 单细胞活性效应", 堵住审稿人
# 输出: regulon_diff_per_celltype.csv (长表) + regulon_keyregulons_summary.csv
#       + FigS6_regulon_heatmap.pdf (关键 regulon x 细胞类型 log2FC 热图)
# ============================================================================
library(hdf5r); library(tidyverse); library(Seurat); library(ggplot2)
source("scripts/fig_utils.R")

## ---- 读取 AUC ----
h5 <- H5File$new("results/scenic_aucell.loom", mode = "r")
cell_ids <- h5[["col_attrs/CellID"]]$read()
auc_df   <- h5[["col_attrs/RegulonsAUC"]]$read()
h5$close_all()
AUC <- as.matrix(auc_df)
rownames(AUC) <- cell_ids
colnames(AUC) <- sub("\\(.*\\)$", "", colnames(AUC))

obj <- readRDS("results/integrated_atlas_scored.rds")
common <- intersect(rownames(AUC), colnames(obj))
AUC  <- AUC[common, , drop = FALSE]
meta <- obj@meta.data[common, ]

## ---- 逐细胞类型比较 ----
MIN_N <- 30
cts <- intersect(c("FCP","Fibro","SSPC_a","SSPC_b","Cycling",
                   "Chondro","Pericyte","Osteo"), unique(meta$celltype))

res_all <- bind_rows(lapply(cts, function(ct) {
  use <- meta$celltype == ct & meta$injury %in% c("fracture", "polytrauma")
  g <- droplevels(factor(meta$injury[use]))
  if (length(g) < 2 * MIN_N || min(table(g)) < MIN_N) {
    message("跳过 ", ct, " (fracture=", sum(g == "fracture"),
            ", polytrauma=", sum(g == "polytrauma"), ")")
    return(NULL)
  }
  x <- AUC[use, , drop = FALSE]
  lapply(colnames(AUC), function(r) {
    f <- suppressWarnings(wilcox.test(x[, r] ~ g))
    tibble(regulon = r, p = f$p.value,
           mean_frac = round(mean(x[g == "fracture", r]), 4),
           mean_poly = round(mean(x[g == "polytrauma", r]), 4))
  }) |> bind_rows() |>
    mutate(padj = p.adjust(p, "BH"),
           log2FC = round(log2((mean_poly + 1e-4) / (mean_frac + 1e-4)), 3),
           celltype = ct)
}))

write.csv(res_all, "results/tables/regulon_diff_per_celltype.csv", row.names = FALSE)

## ---- 关键 regulon 汇总 + 热图 ----
key <- c("Jun","Fos","Jund","Fosb","Egr1","Atf3","Nfkb1","Nfkb2","Stat3",
         "Twist1","Tead1","Sox9","Sox6","Sp7","Runx2","Prrx2","Gli1","Cebpd")
key <- intersect(key, unique(res_all$regulon))
summ <- res_all |> filter(regulon %in% key) |>
  select(celltype, regulon, log2FC, padj) |>
  mutate(regulon = factor(regulon, levels = rev(key)))
write.csv(summ, "results/tables/regulon_keyregulons_summary.csv", row.names = FALSE)

p <- ggplot(summ, aes(celltype, regulon, fill = log2FC)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "log2FC\n(poly/frac)") +
  geom_text(data = summ |> filter(padj < 0.05), aes(label = "*"), size = 5) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Regulon activity: polytrauma vs fracture")
save_fig(p, "results/figures/FigS6_regulon_heatmap.pdf", width = 8, height = 6.5)

## ---- 各细胞类型 top10 打印 ----
for (ct in unique(res_all$celltype)) {
  cat("\n=== ", ct, " top10 (padj 排序) ===\n", sep = "")
  res_all |> filter(celltype == ct) |>
    arrange(padj) |> slice(1:10) |> print()
}
message("07b 完成")
