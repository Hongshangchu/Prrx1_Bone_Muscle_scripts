# ============================================================================
# 08_ML_signature.R —— Aim2e: AI 核心特征基因筛选 (v4 最终版)
# 演进: v1 WGCNA(样本太少报错) -> v2 LASSO(类样本<2被拒) -> v3 XGBoost
#       (n=9 LOOCV全NA) -> v4 定稿: RFx10重复 ∩ 单变量t过滤, 均适配 n=9
# ============================================================================
library(randomForest); library(Seurat); library(tidyverse)
source("scripts/fig_utils.R")

obj <- readRDS("results/integrated_atlas.rds")

## ---- 1. 特征池: injury DEG (top100/类) ∪ 轨迹 DEG ----
injury_deg <- readRDS("results/injury_DEGs.rds")
pool <- bind_rows(lapply(names(injury_deg), function(ct) {
  x <- injury_deg[[ct]]
  if (is.null(x)) return(NULL)
  x |> rownames_to_column("gene") |> mutate(celltype = ct)
})) |>
  filter(p_val_adj < 0.05) |>
  group_by(celltype) |>
  slice_max(order_by = abs(avg_log2FC), n = 100) |>
  pull(gene) |> unique()
pt <- read.csv("results/tables/pseudotime_DEGs.csv")
pt_genes <- pt |> filter(status == "OK", q_value < 0.01) |>
  slice_max(morans_I, n = 100) |> pull(gene_short_name)
pool <- union(pool, pt_genes)
pool <- intersect(pool, rownames(obj))
cat("特征池大小:", length(pool), "\n")

## ---- 2. 样本级表达矩阵 ----
expr <- AverageExpression(obj, features = pool, group.by = "orig.ident",
                          assays = "SCT")$SCT
X <- t(scale(t(log1p(as.matrix(expr)))))
X <- t(X)
rn <- gsub("-", "_", rownames(X))          # Seurat 曾把 "_" 换成 "-"
ii <- match(rn, obj$orig.ident)
stopifnot("仍有未匹配样本" = !any(is.na(ii)))
y <- factor(obj$injury[ii]); names(y) <- rownames(X)
print(table(y))

## ---- 3. 随机森林 x10 重复 (稳定性排名) ----
set.seed(1234)
rf_list <- lapply(1:10, function(i) randomForest(X, y, ntree = 2000, importance = TRUE))
print(rf_list[[10]])                       # 摘要 (OOB error rate)
imp_mean <- sort(rowMeans(sapply(rf_list, function(r)
  importance(r)[, "MeanDecreaseAccuracy"])), decreasing = TRUE)
rf_genes <- names(imp_mean)[1:30]

## ---- 4. 单变量过滤: 三类两两比较最大 |t| (与 RF 互补) ----
tt <- sapply(colnames(X), function(g) {
  s <- split(X[, g], y)
  if (any(sapply(s, sd, na.rm = TRUE) == 0)) return(0)   # 组内恒量 -> 记0, 跳过t检验
  suppressWarnings(max(abs(c(
    t.test(s[[1]], s[[2]])$statistic,
    t.test(s[[1]], s[[3]])$statistic,
    t.test(s[[2]], s[[3]])$statistic))))
})
tt_genes <- names(sort(tt, decreasing = TRUE))[1:30]

## ---- 5. 集成定稿: RF-top30 ∩ 过滤-top30 ----
final_genes <- intersect(rf_genes, tt_genes)
if (length(final_genes) < 5) final_genes <- rf_genes[1:20]
final_genes <- final_genes[1:min(30, length(final_genes))]
write.csv(data.frame(gene = final_genes),
          "results/tables/AI_core_signature.csv", row.names = FALSE)
cat("最终核心特征:", length(final_genes), "\n"); print(final_genes)

## ---- 6. 池基因重聚类验证 ----
sub <- subset(obj, features = pool)
sub <- ScaleData(sub) |> RunPCA(npcs = 20) |> RunUMAP(dims = 1:20)
save_fig(DimPlot(sub, group.by = "injury") | DimPlot(sub, group.by = "celltype"),
         "results/figures/FigS2_poolgenes_umap.pdf", width = 10, height = 5)
message("08 完成")
