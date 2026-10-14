# ============================================================================
# 09b_bulk_identity_tpm.R —— Aim3 (TPM版): Prrx1+ 骨髓 Bulk 的身份鉴定
# 输入: data/bulk/TPM_matrix.csv-your_file_withSymbol.csv (只用 NativePrrx101-103)
# 方法: 与图谱各细胞类型平均表达谱相关性 + AI signature 表达检查
# ============================================================================
library(Seurat); library(tidyverse)

## ---- 1. 清洗: 只取 NativePrrx 三列 ----
raw <- read.csv("data/bulk/TPM_matrix.csv-your_file_withSymbol.csv",
                check.names = FALSE)
bulk <- raw |>
  dplyr::select(gene = GeneSymbol,   # 显式命名空间, 防 AnnotationDbi 遮蔽
         Prrx1_a = `NativePrrx101.bam`,
         Prrx1_b = `NativePrrx102.bam`,
         Prrx1_c = `NativePrrx103.bam`) |>
  filter(!is.na(gene), gene != "") |>
  mutate(across(-gene, as.numeric)) |>
  group_by(gene) |>
  summarise(across(everything(), mean), .groups = "drop")
write.csv(bulk, "data/bulk/prrx1_bulk_tpm.csv", row.names = FALSE)
cat("基因数:", nrow(bulk), "\n")

## log2 转换
bmat <- bulk |> column_to_rownames("gene") |> as.matrix()
bmat <- log2(bmat + 1)

## ---- 2. 与图谱细胞类型相关性 ----
obj <- readRDS("results/integrated_atlas.rds")
ref <- AverageExpression(obj, group.by = "celltype", assays = "SCT")$SCT

common <- intersect(rownames(ref), rownames(bmat))
message("交集基因: ", length(common))
cor_mat <- cor(bmat[common, ], as.matrix(ref[common, ]), method = "spearman")
write.csv(cor_mat, "results/tables/bulk_celltype_correlation.csv")
print(round(cor_mat, 3))

## ---- 3. AI 核心特征基因在 Bulk 中的表达 ----
sig <- read.csv("results/tables/AI_core_signature.csv")$gene
sig <- intersect(sig, rownames(bmat))
sig_expr <- data.frame(gene = sig,
                       mean_log2TPM = rowMeans(bmat[sig, , drop = FALSE]))
write.csv(sig_expr, "results/tables/bulk_signature_expression.csv", row.names = FALSE)
print(arrange(sig_expr, desc(mean_log2TPM)))

## ---- 4. 关键身份基因检查 (应高表达) ----
id_genes <- c("Prrx1","Prrx2","Cxcl12","Lepr","Pdgfra","Ly6a","Col1a1","Dcn",
              "Ptprc","Lyz2","S100a8","C1qa","Hbb-bs")
chk <- data.frame(gene = id_genes,
                  log2TPM = round(rowMeans(bmat[intersect(id_genes, rownames(bmat)), , drop = FALSE]), 2))
write.csv(chk, "results/tables/bulk_identity_check.csv", row.names = FALSE)
print(chk)

## ---- 5. 切片 IF 验证目标 = AI signature ∩ Bulk高表达 ----
high <- rownames(bmat)[rowMeans(bmat) > median(rowMeans(bmat))]
if_targets <- intersect(sig, high)
writeLines(if_targets, "results/IF_validation_targets.txt")
cat("\nIF 验证目标(", length(if_targets), "): ", paste(if_targets, collapse = ", "), "\n")
message("09b 完成")
