# ============================================================================
# 02b_annotate.R —— 最终注释定稿 + 终版 Fig1
# 前置: 01b_decontx.R + 02_integrate_harmony.R 已完成
# 依据: results/tables/singler_crosstab_v3.csv 与 cluster_marker_pct_v3.csv
# ============================================================================
library(Seurat); library(ggplot2)
source("scripts/fig_utils.R")
obj <- readRDS("results/integrated_atlas.rds")
stopifnot(length(levels(obj$seurat_clusters)) == 23)

## ============ 最终映射 (9样本终版, 54,606细胞/23群) ============
ann <- c("0"="Immune","1"="SSPC_a","2"="FCP","3"="Cycling","4"="SSPC_b",
         "5"="FCP","6"="Cycling","7"="SSPC_a","8"="SSPC_a","9"="FCP",
         "10"="Fibro","11"="Immune","12"="Chondro","13"="Immune","14"="Pericyte",
         "15"="FCP","16"="Fibro","17"="Chondro","18"="Fibro","19"="Pericyte",
         "20"="Cycling","21"="Osteo","22"="Eryth")
obj$celltype <- unname(ann[as.character(obj$seurat_clusters)])
if (any(is.na(obj$celltype))) stop("有 cluster 未定义: ",
  paste(unique(obj$seurat_clusters[is.na(obj$celltype)]), collapse = ","))
obj <- subset(obj, celltype != "Eryth")
saveRDS(obj, "results/integrated_atlas.rds")
print(table(obj$celltype))
print(round(prop.table(table(obj$orig.ident, obj$celltype), 1), 2))

## ---- 终版 Fig1 ----
p1 <- DimPlot(obj, group.by = "celltype", label = TRUE) + NoLegend()
p2 <- DimPlot(obj, group.by = "celltype", split.by = "orig.ident", ncol = 3, pt.size = 0.05) + NoLegend()
p3 <- DimPlot(obj, group.by = "injury") |
     (DimPlot(obj, group.by = "age") + NoLegend())   # age 仅一个水平(adult), 图例冗余去掉
save_fig(p1 / p3, "results/figures/Fig1_atlas.pdf", width = 12, height = 10)
save_fig(p2, "results/figures/Fig1b_sample_split.pdf", width = 14, height = 8)

comp <- as.data.frame(prop.table(table(obj$orig.ident, obj$celltype), margin = 1))
colnames(comp) <- c("sample", "celltype", "freq")
p4 <- ggplot(comp, aes(sample, freq, fill = celltype)) + geom_col(width = 0.8) +
  coord_flip() + scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion") + theme_classic(base_size = 11)
save_fig(p4, "results/figures/Fig1c_sample_composition.pdf", width = 9, height = 4.5)
message("注释完成 -> results/integrated_atlas.rds + Fig1 已更新, 可进 03")
