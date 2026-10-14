# ============================================================================
# 03_cellchat_nichenet.R —— Aim2a: 免疫-间充质互作 (v2 修正版)
# NicheNet 使用带日期后缀的真实文件名 (来自 00c 第 4 步)
# ============================================================================
library(CellChat); library(nichenetr); library(Seurat); library(tidyverse)
library(patchwork)
source("scripts/fig_utils.R")

obj <- readRDS("results/integrated_atlas.rds")

## ---- CellChat: fracture vs polytrauma (条件内按 celltype 分组, 官方比较流程) ----
run_cellchat <- function(seu, group.by = "celltype", max.cells = 2000) {
  seu <- subset(seu, celltype %in% c("Immune","SSPC_a","SSPC_b","Osteo","Fibro",
                                     "Chondro","Cycling","FCP","Pericyte"))
  ## 每类抽样加速 (CellChat 在 >2万细胞上很慢; 抽样不改变通讯方向性结论)
  set.seed(1234)
  keep <- unlist(lapply(unique(seu$celltype), function(ct) {
    cs <- colnames(seu)[seu$celltype == ct]
    if (length(cs) > max.cells) sample(cs, max.cells) else cs
  }))
  seu <- subset(seu, cells = keep)
  ## 用 RNA assay 的 log-normalized 数据 (DecontX 校正后 counts -> NormalizeData)
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  message("CellChat 输入: ", ncol(seu), " 细胞, ", length(unique(seu$celltype)), " 群")

  ## 官方推荐: 直接从 Seurat 对象创建 (避免手动抽矩阵的版本兼容问题)
  cc <- createCellChat(seu, group.by = group.by, assay = "RNA")
  cc@DB <- CellChatDB.mouse
  message("DB interactions: ", nrow(cc@DB$interaction),
          " | signaling genes: ", length(intersect(rownames(cc@data), unique(cc@DB$interaction$ligand))))  # 诊断
  cc <- subsetData(cc)
  message("data.signaling dim: ", paste(dim(cc@data.signaling), collapse = " x "))  # 诊断: 若为 0 行说明基因名不匹配
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc
}
cc_frac <- run_cellchat(subset(obj, orig.ident %in% c("Muscle_d3_fracture","Muscle_d5_fracture")))
cc_poly <- run_cellchat(subset(obj, orig.ident %in% c("Muscle_d3_polytrauma","Muscle_d5_polytrauma")))

## 两条件群数不同 (fracture 缺 Osteo): 先 lift 到同一群空间再合并 (CellChat 官方方案)
group.new <- union(levels(cc_frac@idents), levels(cc_poly@idents))
cc_frac <- liftCellChat(cc_frac, group.new = group.new)
cc_poly <- liftCellChat(cc_poly, group.new = group.new)
object.list <- list(fracture = cc_frac, polytrauma = cc_poly)
merged <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(list(frac = cc_frac, poly = cc_poly, merged = merged),
        "results/cellchat_objects.rds")

sig_path <- c("IL6","TNF","CCL","CXCL","COMPLEMENT","OSM","TGFb","BMP")
df_frac <- subsetCommunication(cc_frac, signaling = sig_path)
df_poly <- subsetCommunication(cc_poly, signaling = sig_path)
write.csv(rbind(df_frac %>% mutate(cond = "fracture"),
                df_poly %>% mutate(cond = "polytrauma")),
          "results/tables/ligand_receptor_pairs.csv", row.names = FALSE)

## ---- Fig2: 自绘差异图 (netVisual_diffInteraction/rankNet 与新版 ggplot2 不兼容) ----
## 图A: 差异互作热图 (polytrauma - fracture)
w_frac <- merged@net[["fracture"]]$weight
w_poly <- merged@net[["polytrauma"]]$weight
dd <- as.data.frame(as.table(w_poly - w_frac))
colnames(dd) <- c("sender", "receiver", "delta")
pA <- ggplot(dd, aes(receiver, sender, fill = delta)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Delta communication (polytrauma - fracture)", fill = "delta")

## 图B: 关键炎症通路通讯强度对比 ( dodge 条形图 )
np_frac <- subsetCommunication(cc_frac, slot.name = "netP") |> rename(pathway = pathway_name)
np_poly <- subsetCommunication(cc_poly, slot.name = "netP") |> rename(pathway = pathway_name)
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
write.csv(np, "results/tables/pathway_communication_strength.csv", row.names = FALSE)

## ---- NicheNet v2: 免疫配体 -> 间充质修复程序 ----
lr_network           <- readRDS("nichenetr/lr_network_mouse_21122021.rds")
ligand_target_matrix <- readRDS("nichenetr/ligand_target_matrix_nsga2r_final_mouse.rds")

geneset_oi <- c("Sox9","Col2a1","Acan","Runx2","Sp7","Col1a1","Alpl","Fn1","Cthrc1")
receiver <- c("SSPC_a","SSPC_b","Chondro","Osteo","FCP"); sender <- "Immune"

expressed_genes_receiver <- rownames(GetAssayData(subset(obj, celltype %in% receiver), assay = "SCT"))
expressed_genes_sender   <- rownames(GetAssayData(subset(obj, celltype == sender),  assay = "SCT"))

## 与 NicheNet 基因宇宙取交集, 防止下标越界 (部分目标基因不在矩阵基因集内)
geneset_oi <- intersect(geneset_oi, colnames(ligand_target_matrix))
ligands <- unique(lr_network$from)
expressed_ligands <- Reduce(intersect, list(ligands, expressed_genes_sender,
                                            rownames(ligand_target_matrix)))
## 仅保留对目标基因集有至少一条预测链接的配体 (>0, 而非>0.05; 0.05会把配体滤光)
potential_ligands <- expressed_ligands[
  rowSums(ligand_target_matrix[expressed_ligands, geneset_oi, drop = FALSE] > 0) > 0]
background_expressed_genes <- expressed_genes_receiver   # v2 必填参数

nn <- predict_ligand_activities(geneset = geneset_oi,
                                background_expressed_genes = background_expressed_genes,
                                ligand_target_matrix = ligand_target_matrix,
                                potential_ligands = potential_ligands)
write.csv(nn, "results/tables/nichenet_ligand_activity.csv", row.names = FALSE)

## top 配体下游靶基因 (配体->靶点网络图数据)
## v1/v2 列名自适应: 校正AUPR 与 配体名列
corr_col  <- intersect(c("test_aupr_corrected","aupr_corrected"), colnames(nn))[1]
lig_col   <- intersect(c("test_ligand","ligand"), colnames(nn))[1]
stopifnot(!is.na(corr_col), !is.na(lig_col))
top_ligands <- nn |> arrange(desc(.data[[corr_col]])) |> slice(1:10) |> pull(.data[[lig_col]])
active_links <- lapply(top_ligands, function(l) {
  w <- ligand_target_matrix[l, ]
  tibble(ligand = l, target = names(w[w > 0.10]), weight = w[w > 0.10])
}) |> bind_rows() |> filter(target %in% geneset_oi)
write.csv(active_links, "results/tables/nichenet_ligand_target_links.csv",
          row.names = FALSE)
