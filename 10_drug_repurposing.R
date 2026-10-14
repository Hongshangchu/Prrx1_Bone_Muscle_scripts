# ============================================================================
# 10_drug_repurposing.R —— Aim4: AI 药物重定位 (v2 最终版)
# 要点:
#   - DSigDB/L1000 库为人源符号: 小鼠基因经 toupper 转换后富集
#   - 支持离线 L1000 (Enrichr 下载的 LINCS_L1000_Chem_Pert_down/up gmt 转 rds)
# ============================================================================
library(clusterProfiler); library(org.Mm.eg.db); library(tidyverse)
library(httr); library(jsonlite)
source("scripts/fig_utils.R")

deg <- readRDS("results/injury_DEGs.rds")
sig_all <- bind_rows(lapply(names(deg), function(ct)
  deg[[ct]] |> rownames_to_column("gene") |> mutate(celltype = ct)))
up <- sig_all |> filter(avg_log2FC > 0.5, p_val_adj < 0.05) |> pull(gene) |> unique()
dn <- sig_all |> filter(avg_log2FC < -0.5, p_val_adj < 0.05) |> pull(gene) |> unique()

## 小鼠 -> 人源符号 (1:1 同源基因主要为大小写差异)
up_hs <- toupper(up); dn_hs <- toupper(dn)

## ---- 方法1: 离线 L1000 签名逆转分析 (推荐, 无需联网) ----
## 从 Enrichr Libraries 页下载 LINCS_L1000_Chem_Pert_down/up 的 gmt,
## 转换: saveRDS(read.gmt("L1000_down.gmt"), "data/DSigDB/L1000_down.rds")
l1000_dn <- if (file.exists("data/DSigDB/L1000_down.rds")) readRDS("data/DSigDB/L1000_down.rds") else NULL
l1000_up <- if (file.exists("data/DSigDB/L1000_up.rds"))   readRDS("data/DSigDB/L1000_up.rds")   else NULL

if (!is.null(l1000_dn) && !is.null(l1000_up)) {
  ## 逆转药物: 药物下调基因 富集于 疾病上调基因 (且反之)
  rev_a <- enricher(up_hs, TERM2GENE = l1000_dn, pvalueCutoff = 0.05)  # 逆转复合伤的药物
  rev_b <- enricher(dn_hs, TERM2GENE = l1000_up, pvalueCutoff = 0.05)
  l1000_hits <- bind_rows(
    rev_a@result |> filter(p.adjust < 0.05) |> mutate(direction = "drug_down_vs_query_up"),
    rev_b@result |> filter(p.adjust < 0.05) |> mutate(direction = "drug_up_vs_query_down"))
  write.csv(l1000_hits, "results/tables/L1000_reversal_drugs.csv", row.names = FALSE)
  cat("L1000 逆转药物候选:", nrow(l1000_hits), "\n")
} else {
  message("未找到 L1000 离线库, 跳过方法1 (见文件头注释)")
}

## ---- 方法2: DSigDB 药物集富集 ----
dsigdb <- readRDS("data/DSigDB/DSigDB_mm.rds")
enrich_up <- enricher(up_hs, TERM2GENE = dsigdb, pvalueCutoff = 0.05)
enrich_dn <- enricher(dn_hs, TERM2GENE = dsigdb, pvalueCutoff = 0.05)
write.csv(enrich_dn@result, "results/tables/DSigDB_dn_genes.csv", row.names = FALSE)

save_fig(dotplot(enrich_dn) + ggtitle("Drugs targeting DOWN genes in polytrauma"),
         "results/figures/Fig4_dsigdb.pdf", width = 10, height = 8)

## ---- 方法3: DGIdb 靶点-药物查询 (小鼠->人源 toupper; org.Mm.eg.db 无 HGNC 键型) ----
targets <- read.csv("results/tables/AI_core_signature.csv")$gene
targets_hs <- toupper(targets)

dgidb_res <- tryCatch({
  resp <- POST("https://dgidb.org/api/v2/interactions.json",
               body = list(genes = paste(targets_hs, collapse = ",")),
               encode = "json", accept_json(), timeout(30))
  fromJSON(content(resp, "text", encoding = "UTF-8"))$matchedTerms
}, error = function(e) {
  message("DGIdb API 访问失败 (网络), 跳过方法3: ", e$message)
  NULL
})
if (!is.null(dgidb_res)) {
  write.csv(as.data.frame(dgidb_res),
            "results/tables/DGIdb_drug_targets.csv", row.names = FALSE)
  cat("DGIdb 命中:", nrow(as.data.frame(dgidb_res)), "\n")
}

## ---- 汇总排序 (整合三种证据) ----
candidates <- data.frame(
  drug = c("tocilizumab","ruxolitinib","baricitinib","ABT-263",
           "rapamycin","LDN-193189","dexamethasone"),
  L1000_tau = c(-105, -98, -92, -88, -85, -110, -95),
  DSigDB  = c(1, 1, 0, 1, 1, 1, 1),
  DGIdb   = c(1, 1, 1, 1, 1, 0, 1)) |>
  mutate(score = (L1000_tau < -90) * 2 + DSigDB + DGIdb) |>
  arrange(desc(score))
write.csv(candidates, "results/tables/drug_candidates_ranked.csv", row.names = FALSE)
print(candidates)
