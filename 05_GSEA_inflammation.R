# ============================================================================
# 05_GSEA_inflammation.R —— Aim2c/d: 炎症通路 GSEA + 模块得分
# ============================================================================
library(clusterProfiler); library(msigdbr); library(fgsea)
library(Seurat); library(tidyverse)

obj <- readRDS("results/integrated_atlas.rds")

Idents(obj) <- "celltype"
injury_deg <- lapply(unique(obj$celltype), function(ct) {
  sub <- subset(obj, celltype == ct)
  if (length(unique(sub$injury)) < 2) return(NULL)
  Idents(sub) <- "injury"
  tryCatch(FindMarkers(sub, ident.1 = "polytrauma", ident.2 = "fracture",
                       logfc.threshold = 0.25, min.pct = 0.1),
           error = function(e) NULL)
})
names(injury_deg) <- unique(obj$celltype)
saveRDS(injury_deg, "results/injury_DEGs.rds")

hall <- msigdbr(species = "Mus musculus", category = "H")
hall_list <- split(hall$gene_symbol, hall$gs_name)

gsea_run <- function(deg, label) {
  ranks <- deg$avg_log2FC; names(ranks) <- rownames(deg)
  ranks <- sort(ranks, decreasing = TRUE)
  fg <- fgsea(pathways = hall_list, stats = ranks, minSize = 10, maxSize = 500)
  fg$comparison <- label; fg
}
gsea_res <- bind_rows(lapply(names(injury_deg), function(ct)
  if (!is.null(injury_deg[[ct]])) gsea_run(injury_deg[[ct]], ct) else NULL))

focus <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_COMPLEMENT","HALLMARK_INFLAMMATORY_RESPONSE",
           "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_KRAS_SIGNALING_UP")
gsea_out <- gsea_res |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))  # list列转字符串
filter(gsea_out, pathway %in% focus) |>
  write.csv("results/tables/GSEA_inflammation.csv", row.names = FALSE)

inflam_genes <- list(IL6_JAK_STAT3 = c("Il6","Il6st","Jak1","Stat3","Socs3"),
                     NFkB          = c("Rela","Nfkb1","Tnfaip3","Nfkbia"),
                     SASP          = c("Cdkn1a","Cdkn2a","Il6","Cxcl1","Mmp13"))
obj <- AddModuleScore(obj, features = inflam_genes, name = c("IL6","NFkB","SASP"))
saveRDS(obj, "results/integrated_atlas_scored.rds")
cat("GSEA & module scoring done.\n")
