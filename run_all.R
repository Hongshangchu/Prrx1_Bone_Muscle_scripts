# ============================================================================
# run_all.R —— 一键主控 (数据/依赖就绪后运行)
# ============================================================================
steps <- c(
  "scripts/00d_organize_geo.R",         # 解包+整理 (首次运行)
  "scripts/01b_decontx.R",             # 环境RNA校正 (Fig1c 验证后启用)
  "scripts/02b_annotate.R",            # 注释定稿 (02之后)
  "scripts/01_load_QC.R",
  "scripts/02_integrate_harmony.R",   # 运行后先人工注释 celltype 再继续!
  "scripts/03_cellchat_nichenet.R",
  "scripts/04_trajectory.R",
  "scripts/05_GSEA_inflammation.R",
  "scripts/06a_export_loom.R",
  # Python 环境运行: python scripts/06_scenic.py
  "scripts/07_regulon_diff.R",
  "scripts/08_ML_signature.R",
  "scripts/09_bulk_rnaseq.R",
  "scripts/10_drug_repurposing.R",
  "scripts/11_validation_in_vitro.R")

for (s in steps) {
  message("========== Running: ", s, " ==========")
  source(s)
}
