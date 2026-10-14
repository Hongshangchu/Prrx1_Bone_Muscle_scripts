# ============================================================================
# 09_bulk_rnaseq.R —— Aim3: 自有 Bulk RNA-seq + 单细胞反卷积对接
# 注: 自有 Bulk 仅稳态 Prrx1+ (n=3, 单组), 无法做组间差异检验 ->
#     不再做 DESeq 差异分析 (design ~ condition 会因单水平因子报错),
#     只做: 标准化表达量 + BisqueRNA 反卷积 + 与单细胞骨折激活 DEG 交叉
# ============================================================================
library(DESeq2); library(Seurat); library(BisqueRNA); library(tidyverse)

## ---- 1. 稳态骨髓 Prrx1+ FACS, n=3 (单组: design ~ 1, 只取标准化表达) ----
counts  <- read.csv("data/bulk/prrx1_bulk_counts.csv", row.names = 1, check.names = FALSE)
message("Bulk 样本: ", ncol(counts), " 个 -> ", paste(colnames(counts), collapse = ", "))
coldata <- data.frame(condition = factor(rep("steady", ncol(counts))),
                      row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ 1)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
res_df <- data.frame(baseMean = rowMeans(norm_counts), check.names = FALSE)
write.csv(res_df, "results/tables/bulk_norm_expression.csv")  # 标准化表达 (替代原 bulk_DEG_steady.csv)
message("Bulk 标准化表达: ", nrow(res_df), " 基因 -> results/tables/bulk_norm_expression.csv")

## ---- 2. BisqueRNA 反卷积 (Seurat -> ExpressionSet 手动转换, 避免函数版本依赖) ----
obj <- readRDS("results/integrated_atlas_scored.rds")
seu2eset <- function(seu) {
  em <- GetAssayData(seu, assay = "SCT", layer = "data")
  ## BisqueRNA 要求 phenoData 含 cellType (细胞类型) 与 SubjectName (个体/样本)
  ## 两个标签, 原脚本均缺失, 必报错; 本数据每细胞一个样本, 故用 orig.ident 充当个体
  pd <- data.frame(cellType = as.character(seu$celltype),
                   SubjectName = as.character(seu$orig.ident),
                   row.names = colnames(seu))
  Biobase::ExpressionSet(assayData = as.matrix(em),
                         phenoData = Biobase::AnnotatedDataFrame(pd))
}
ref_cells <- subset(obj, orig.ident %in% c("Periosteum_SS","Muscle_d0_none","BM_SS"))
ref.eset  <- seu2eset(ref_cells)
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(counts))
## use.overlap = FALSE: 自有 Bulk 与公共 sc 参考无重叠个体,
## 默认 TRUE 会报 "No overlapping samples"; FALSE 即标准跨队列参考模式
dec <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, ref.eset, use.overlap = FALSE)
print(dec$bulk.props)
write.csv(dec$bulk.props, "results/tables/bulk_deconvolution.csv")

## ---- 3. 交叉: 稳态 Prrx1+ (Bulk 表达) x 骨折激活 (GSE138689 BM d7 vs d0) ----
bm <- subset(obj, tissue == "BM")
Idents(bm) <- "time"
bm_deg <- FindMarkers(bm, ident.1 = "d7", ident.2 = "d0", logfc.threshold = 0.25)
cross <- bm_deg |>
  rownames_to_column("gene") |>
  inner_join(res_df |> rownames_to_column("gene"), by = "gene") |>
  mutate(score = avg_log2FC * log1p(baseMean)) |>   # sc 激活倍数 x bulk 表达量
  arrange(desc(score))
write.csv(cross, "results/tables/steady_vs_activated.csv", row.names = FALSE)

## ---- 4. 切片 IF 验证目标清单 = AI signature ∩ Bulk 高表达 ----
sig_file <- "results/tables/AI_core_signature.csv"
if (file.exists(sig_file)) {
  sig <- read.csv(sig_file)$gene
} else {
  message("!! 未找到 AI_core_signature.csv (08 未完成), 暂用 cross top 30 基因")
  sig <- head(cross$gene, 30)
}
bulk_mean <- res_df$baseMean; names(bulk_mean) <- rownames(res_df)
candidates <- intersect(sig, names(bulk_mean)[bulk_mean > median(bulk_mean)])
writeLines(candidates, "results/IF_validation_targets.txt")
# 切片染色: tdTomato + Cxcl12/Lepr (身份) + p-STAT3/NF-kB (炎症) + top 候选基因
