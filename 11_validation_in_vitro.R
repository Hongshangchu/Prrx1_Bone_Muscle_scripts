# ============================================================================
# 11_validation_in_vitro.R —— Aim4: 外植体 rescue 实验分析模板
# 设计: 骨膜外植体 + IL-6/TNF 处理, 候选药物 (Toci/Ruxo) 共处理 rescue
# ============================================================================
library(ggplot2); library(ggpubr); library(patchwork); library(tidyverse)
source("scripts/fig_utils.R")

df <- read.csv("data/functional/qPCR_rescue.csv")   # 列: gene, group, value
genes <- c("Sox9","Col2a1","Runx2","Alpl","Cxcl1","Il6")

plist <- lapply(genes, function(g) {
  d <- subset(df, gene == g)
  ggplot(d, aes(group, value, fill = group)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 1) +
    stat_compare_means(comparisons = list(c("IL6","IL6+Toci"),
                                          c("IL6","IL6+Ruxo"))) +
    ggtitle(g) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
})
save_fig(wrap_plots(plist, ncol = 3),
         "results/figures/Fig6_rescue.pdf", width = 12, height = 8)
