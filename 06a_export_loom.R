# ============================================================================
# 06a_export_loom.R —— 为 pySCENIC 导出表达矩阵 (抽样版)
# 全量 3-4 万细胞对 GRNBoost2 过重; 按 MAX_PER_CT 抽样 (set.seed 固定, 可复现)
# ============================================================================
library(Seurat)

obj <- readRDS("results/integrated_atlas_scored.rds")
## 注意: 不过滤 injury —— SCENIC 需含 polytrauma 细胞才能做组间比较
msc <- subset(obj, celltype %in% c("SSPC_a","SSPC_b","FCP","Osteo","Chondro",
                                   "Fibro","Cycling","Pericyte"))

MAX_PER_CT <- 1000  # 每类抽样数 (8类约7-8k细胞; 卡调度就调小)
set.seed(1234)
keep <- unlist(lapply(unique(msc$celltype), function(ct) {
  cs <- colnames(msc)[msc$celltype == ct]
  if (length(cs) > MAX_PER_CT) sample(cs, MAX_PER_CT) else cs
}))
msc <- subset(msc, cells = keep)
message("SCENIC 输入: ", ncol(msc), " 细胞 x ",
        length(unique(msc$celltype)), " 群")

mat <- as.matrix(GetAssayData(msc, assay = "SCT", layer = "counts"))
write.csv(mat, "results/msc_counts.csv")
cat("导出完成: results/msc_counts.csv\n")
