# ============================================================================
# 04_trajectory.R —— Aim2b: Monocle3 拟时轨迹 (FCP->纤维软骨->成熟软骨->成骨)
# 构建方式: new_cell_data_set 原生构建 (不依赖 SeuratWrappers 的转换器)
# ============================================================================
library(monocle3); library(Seurat)
## dplyr/tidyr 必须显式加载: 单独运行本脚本时若未加载,
## count 等函数会被其他包同名函数遮蔽报错 (全流程跑时因前置脚本
## 已加载 tidyverse 而侥幸通过)。dplyr 放在后面加载以保证遮蔽优先。
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
source("scripts/fig_utils.R")

obj <- readRDS("results/integrated_atlas.rds")
## 注意: 不过滤 injury —— polytrauma 细胞必须进入轨迹才能做组间拟时比较
msc <- subset(obj, celltype %in% c("SSPC_a","SSPC_b","Osteo","Chondro","Fibro",
                                   "Cycling","Pericyte","FCP"))

## 每类抽样至 2500, 控制 monocle3 运行时间 (~17k 细胞)
set.seed(1234)
keep <- unlist(lapply(unique(msc$celltype), function(ct) {
  cs <- colnames(msc)[msc$celltype == ct]
  if (length(cs) > 2500) sample(cs, 2500) else cs
}))
msc <- subset(msc, cells = keep)

## ---- 原生构建 CDS ----
expr <- GetAssayData(msc, assay = "SCT", layer = "counts")
meta <- msc@meta.data
fd   <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr))

cds <- new_cell_data_set(expression_data = expr,
                         cell_metadata = meta,
                         gene_metadata = fd)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)          # monocle3 自己的 UMAP
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

## 根 = 稳态 SSPC_a (d0); 不足则放宽到全部 d0 基质
root_cells <- colnames(cds)[colData(cds)$time == "d0" &
                            colData(cds)$celltype == "SSPC_a"]
if (length(root_cells) == 0)
  root_cells <- colnames(cds)[colData(cds)$time == "d0"]
cds <- order_cells(cds, root_cells = root_cells[seq_len(min(100, length(root_cells)))])
saveRDS(cds, "results/monocle_cds.rds")

save_fig_multi(list(
  plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE),
  plot_cells(cds, color_cells_by = "celltype", label_cell_groups = TRUE),
  plot_cells(cds, color_cells_by = "time"),
  plot_cells(cds, color_cells_by = "injury")),        # 复合伤/骨折细胞沿拟时分布
  "results/figures/Fig3_trajectory.pdf", width = 18, height = 5)

## 拟时相关基因 (轨迹决策点)
deg_pt <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
write.csv(deg_pt, "results/tables/pseudotime_DEGs.csv")
message("04 完成 -> Fig3_trajectory.pdf + pseudotime_DEGs.csv")


## ============================================================================
## 重要: 多组织图谱在 monocle3 主图上碎成多岛, ~90% 损伤细胞拟时为 Inf,
## 图基拟时不可用于组间比较。改用 "分化模块得分(DiffScore)" 作为正文量化:
## 基于轨迹 DEG (Moran's I) 的 late-early 基因模块, 在全细胞上计算, 无 Inf。
## ============================================================================
pt <- pseudotime(cds)
inf_pct <- dplyr::tibble(injury = colData(cds)$injury, inf = is.infinite(pt)) |>
  dplyr::count(injury, inf) |> dplyr::group_by(injury) |>
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
write.csv(inf_pct, "results/tables/pseudotime_inf_fraction.csv", row.names = FALSE)
print(inf_pct)

## late = ECM/软骨成熟程序 (轨迹 top Moran's I), early = 祖细胞程序
late_genes  <- c("Col9a1","Col9a3","Serpine2","Pcp4l1","Serping1","Fn1","Fmod",
                 "Dpt","Bgn","Timp1","Col3a1","Col5a2","Igfbp5","Ucma","Thbs1")
early_genes <- c("Prrx1","Prrx2","Grem1","Fgf7","Twist2")
## 用 02b 注释后的 obj (含 celltype/injury), 重新读入保证完整
obj <- readRDS("results/integrated_atlas.rds")
late_genes  <- intersect(late_genes,  rownames(obj))
early_genes <- intersect(early_genes, rownames(obj))
obj <- AddModuleScore(obj, features = list(late_genes, early_genes),
                      name = c("Late", "Early"))
obj$DiffScore <- obj$Late1 - obj$Early2   # AddModuleScore 按列表序号命名: Late1/Early2
saveRDS(obj, "results/integrated_atlas_scored.rds")

df_ds <- obj@meta.data |>
  dplyr::filter(celltype %in% c("FCP","Cycling","Chondro","Fibro","SSPC_a"),
                injury %in% c("fracture","polytrauma"))
tab <- df_ds |> dplyr::group_by(celltype, injury) |>
  dplyr::summarise(m = round(mean(DiffScore), 3), .groups = "drop") |>
  tidyr::pivot_wider(names_from = injury, values_from = m)
write.csv(tab, "results/tables/diffscore_by_celltype.csv", row.names = FALSE)
print(tab)

pval <- df_ds |> dplyr::group_by(celltype) |>
  dplyr::summarise(p = tryCatch(wilcox.test(DiffScore[injury=="fracture"],
                                     DiffScore[injury=="polytrauma"])$p.value,
                         error = function(e) NA), .groups = "drop")
write.csv(pval, "results/tables/diffscore_pvalues.csv", row.names = FALSE)
print(pval)

p_ds <- ggplot(df_ds, aes(celltype, DiffScore, fill = injury)) +
  geom_boxplot(outlier.shape = NA) + theme_classic()
save_fig(p_ds, "results/figures/FigS4c_diffscore_by_ct.pdf", width = 8, height = 4)
message("04 完成: Fig3(参考) + DiffScore(正文) 已输出")
