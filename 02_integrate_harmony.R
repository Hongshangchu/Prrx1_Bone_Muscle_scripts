# ============================================================================
# 02_integrate_harmony.R —— Aim1: Harmony 整合、聚类、注释、Fig1 图谱
# 注意: 运行后须根据 FigS1_marker_dotplot.pdf 人工注释 celltype (见文件末尾)
# ============================================================================
library(Seurat); library(harmony); library(patchwork); library(ggplot2)
library(SingleR); library(celldex)
source("scripts/fig_utils.R")

## 整合对象为全细胞合并版 (~27GB), future 多会话并行会把整个对象
## 序列化发给每个工作线程 -> 触发 globals.maxSize 上限且极慢。
## 内存才是本步瓶颈, 全程串行最稳 (00_setup.R 的 plan 对小对象仍有效)
future::plan("sequential")

## 若存在 DecontX 校正结果 (01b 输出), 优先使用; 否则用 01 的 QC 对象
objs <- if (file.exists("results/decontx_objects.rds")) {
  message("使用 DecontX 9样本对象 (results/decontx_objects.rds)")
  readRDS("results/decontx_objects.rds")
} else {
  message("未找到校正结果, 使用原始 QC 对象 (注意: 环境 RNA 可能污染聚类)")
  readRDS("results/qc_filtered_objects.rds")
}
obj  <- merge(objs[[1]], y = objs[-1])
obj  <- JoinLayers(obj)

obj <- SCTransform(obj, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
obj <- RunPCA(obj, npcs = 50)
## 注: 不用 IntegrateLayers —— 在 JoinLayers+SCT 的对象上它会报
## "contrasts 只适用于两个或多于两个水平的因子" (Seurat v5 已知问题)。
## 直接 RunHarmony 等价且更稳。
obj <- harmony::RunHarmony(obj, group.by.vars = "orig.ident",
                           reduction = "pca", reduction.save = "harmony",
                           verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.6)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)

## SingleR 自动注释初筛
ref  <- celldex::MouseRNAseqData()
pred <- SingleR(GetAssayData(obj, layer = "data", assay = "SCT"),
                ref = ref, labels = ref$label.main)
obj$singler_label <- pred$labels[match(colnames(obj), rownames(pred))]

## 文献 marker (Colnot / Xing 体系)
markers <- list(
  SSPC_Pdgfra  = c("Pdgfra","Ly6a","Cd34"),
  SSPC_Pdgfrb  = c("Pdgfrb","Cspg4","Rgs5"),
  Periosteal   = c("S100a4","Thy1","Col1a1","Fn1"),
  FibroChondro = c("Sox9","Col2a1","Acan"),
  Osteo        = c("Runx2","Sp7","Bglap","Alpl"),
  Chondro      = c("Col10a1","Mmp13"),
  Adipo        = c("Adipoq","Pparg","Fabp4"),
  Pericyte     = c("Notch3","Kcnj8"),
  Endo         = c("Cdh5","Pecam1","Kdr"),
  Immune       = c("Ptprc","Lyz2","Cd3e","Nkg7","Cd79a","Ms4a1"))
## 去重 (Sox9/Col2a1 等重复会导致 factor level duplicated 错误) + 过滤缺失基因
feats <- unlist(markers)
feats <- feats[!duplicated(feats)]
feats <- feats[feats %in% rownames(obj)]
save_fig(DotPlot(obj, features = feats) + RotatedAxis(),
         "results/figures/FigS1_marker_dotplot.pdf", width = 14, height = 8)

## >>> 最终注释 (经 marker dotplot x2 + SingleR + 分选策略 + 定量阳性率四重验证) <<<
## 注意: 骨痂/损伤组织环境 RNA 极强 (S100a8/Lyz2/Ptprc/Col1a1 均受污染),
## 注释依据: 细胞内在 marker (C1qa/Ms4a7/Cspg4/Runx2) + 样本分选策略
## 用命名向量映射 (base R, 不依赖 dplyr::recode)
## DecontX 后的最终注释 (29 clusters, 依据 results/tables/cluster_marker_pct.csv)
## 历史教训: 聚类编号随输入数据改变, 自动套旧映射曾造成错标签。
## 现改为占位保存, 注释在外部 02b_annotate.R 定稿 (核对 singler/pct 表后填入)。
obj$celltype <- paste0("C", obj$seurat_clusters)
saveRDS(obj, "results/integrated_atlas.rds")

p1 <- DimPlot(obj, group.by = "celltype", label = TRUE) + NoLegend()
p3 <- DimPlot(obj, group.by = "injury") | DimPlot(obj, group.by = "age")
save_fig(p1 / p3, "results/figures/Fig1_atlas.pdf", width = 12, height = 10)

## Fig1b: 样本贡献分面图 (9 样本各一格, celltype 着色, 避免拥挤)
p2 <- DimPlot(obj, group.by = "celltype", split.by = "orig.ident",
              ncol = 4, pt.size = 0.05) + NoLegend()
save_fig(p2, "results/figures/Fig1b_sample_split.pdf", width = 16, height = 9)

## Fig1c: 样本 x 细胞类型组成柱状图 (比 UMAP 更适合展示样本贡献)
comp <- as.data.frame(prop.table(table(obj$orig.ident, obj$celltype), margin = 1))
colnames(comp) <- c("sample", "celltype", "freq")
comp$sample <- factor(comp$sample, levels = names(sort(tapply(
  as.numeric(obj$celltype == "Immune"), obj$orig.ident, mean))))
p4 <- ggplot(comp, aes(sample, freq, fill = celltype)) +
  geom_col(width = 0.8) + coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")
save_fig(p4, "results/figures/Fig1c_sample_composition.pdf", width = 9, height = 5)
cat("Integration done.\n")
