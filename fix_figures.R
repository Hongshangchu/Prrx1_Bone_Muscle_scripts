# ============================================================================
# fix_figures.R —— 重建损坏的 Fig2 / Fig3
# v2 修复: Fig3 弃用易损坏的 base pdf() 设备, 改用 ggsave(plot_cells 返回 ggplot)
# ============================================================================
library(Seurat); library(tidyverse); library(ggplot2); library(patchwork)
library(CellChat); library(monocle3)
source("scripts/fig_utils.R")

## ---- Fig2: CellChat 差异 (ggsave, 稳定) ----
cc <- readRDS("results/cellchat_objects.rds")
merged <- cc$merged

w_frac <- merged@net[["fracture"]]$weight
w_poly <- merged@net[["polytrauma"]]$weight
dd <- as.data.frame(as.table(w_poly - w_frac))
colnames(dd) <- c("sender", "receiver", "delta")
pA <- ggplot(dd, aes(receiver, sender, fill = delta)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Delta communication (polytrauma - fracture)", fill = "delta")

np_frac <- subsetCommunication(cc$frac, slot.name = "netP") |> rename(pathway = pathway_name)
np_poly <- subsetCommunication(cc$poly, slot.name = "netP") |> rename(pathway = pathway_name)
np <- bind_rows(mutate(np_frac, cond = "fracture"),
                mutate(np_poly, cond = "polytrauma")) |>
  group_by(pathway, cond) |> summarise(total = sum(prob), .groups = "drop")
focus <- c("IL6","TNF","CCL","CXCL","COMPLEMENT","OSM","TGFb","BMP","IL1","CSF")
pB <- ggplot(np |> filter(pathway %in% focus),
             aes(reorder(pathway, total), total, fill = cond)) +
  geom_col(position = "dodge") + coord_flip() + theme_classic() +
  labs(x = NULL, y = "Total communication prob", fill = NULL) +
  theme(legend.position = "top")

save_fig(pA | pB, "results/figures/Fig2_cellchat_comparison.pdf", width = 15, height = 6.5)

## ---- Fig3: 轨迹四联图 —— 关键修复: plot_cells 返回 ggplot, 用 ggsave 不用 base pdf() ----
cds <- readRDS("results/monocle_cds.rds")
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE)
p2 <- plot_cells(cds, color_cells_by = "celltype", label_cell_groups = TRUE)
p3 <- plot_cells(cds, color_cells_by = "time")
p4 <- plot_cells(cds, color_cells_by = "injury")
save_fig(p1 | p2 | p3 | p4,
         "results/figures/Fig3_trajectory.pdf", width = 20, height = 5.5,
         limitsize = FALSE)

## ---- 600 dpi TIFF 已由 save_fig 同步输出 (替代旧 PNG 兜底) ----
message("Fig2 / Fig3 重建完成 (PDF + 600 dpi TIFF)")
