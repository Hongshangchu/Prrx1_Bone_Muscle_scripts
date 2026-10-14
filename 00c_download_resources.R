# ============================================================================
# 00c_download_resources.R —— SCENIC / NicheNet / DSigDB 自动下载 (v3 最终版)
# 要点:
#   - NicheNet 真实文件名带日期后缀 (lr_network_mouse_21122021.rds 等)
#   - DSigDB 优先 Enrichr API, 备选 tanlab 官方站
#   - 已存在的文件自动跳过; 失败打印浏览器直链
# ============================================================================

robust_dl <- function(urls, dest, min_size = 1e5, tries = 3) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(dest) && file.size(dest) > min_size) {
    message("已存在, 跳过: ", dest); return(invisible(TRUE))
  }
  for (u in urls) {
    for (i in 1:tries) {
      message("下载 (", i, "/", tries, "): ", basename(dest), " <- ", u)
      ok <- tryCatch({
        download.file(u, dest, method = "libcurl", mode = "wb", quiet = TRUE)
        file.exists(dest) && file.size(dest) > min_size
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) {
        message("成功: ", dest, " (", round(file.size(dest)/1e6, 1), " MB)")
        return(invisible(TRUE))
      }
      Sys.sleep(3)
    }
  }
  message("!! 失败, 请浏览器手动下载: ", paste(urls, collapse = " 或 "))
  invisible(FALSE)
}

## ---- 1. SCENIC: mm10 cisTarget 数据库 (约 1 GB) ----
## 本地文件名: ...mc9nr.genes_vs_motifs.rankings.feather
robust_dl("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
          "scenic/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather", min_size = 5e8)

## ---- 2. SCENIC: motif 注释 ----
robust_dl("https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
          "scenic/motifs-v9-nr.mgi-m0.001-o0.0.tbl")

## ---- 3. 小鼠 TF 列表 (jsdelivr CDN 优先, 国内可达性好) ----
## 本地文件名: mm_mgi_tfs.txt
robust_dl(c("https://cdn.jsdelivr.net/gh/aertslab/pySCENIC@master/resources/mm_tfs.txt",
            "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/mm_tfs.txt"),
          "scenic/mm_mgi_tfs.txt", min_size = 1e3)

## ---- 4. NicheNet v2 网络 (Zenodo 7074291; 文件名带日期后缀!) ----
for (f in c("lr_network_mouse_21122021.rds",
            "signaling_network_mouse_21122021.rds",
            "ligand_target_matrix_nsga2r_final_mouse.rds",
            "weighted_networks_nsga2r_final_mouse.rds")) {
  robust_dl(paste0("https://zenodo.org/records/7074291/files/", f, "?download=1"),
            file.path("nichenetr", f), min_size = 1e4)
}

## ---- 5. DSigDB (Enrichr API; 失败则 tanlab 官方站) ----
dsig_ok <- robust_dl(
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DSigDB_2024",
  "data/DSigDB/DSigDB_2024.gmt", min_size = 1e5)
if (!dsig_ok) {
  robust_dl("https://dsigdb.tanlab.org/DSigDBv1.0/downloads/DSigDB_All.gmt",
            "data/DSigDB/DSigDB_All.gmt", min_size = 1e5)
}

## ---- 6. GMT -> R 对象 ----
library(clusterProfiler)

## DSigDB (优先用浏览器手动下载的 DSigDB.gmt)
gmt_candidates <- c("data/DSigDB/DSigDB.gmt",
                    "data/DSigDB/DSigDB_2024.gmt",
                    "data/DSigDB/DSigDB_All.gmt")
gmt_file <- gmt_candidates[file.exists(gmt_candidates)][1]
if (!is.na(gmt_file)) {
  saveRDS(clusterProfiler::read.gmt(gmt_file), "data/DSigDB/DSigDB_mm.rds")
  message("DSigDB_mm.rds 已生成 <- ", gmt_file)
}

## LINCS L1000 化学扰动签名 (供 10_drug_repurposing.R 离线逆转分析)
## 显式写法, 避免命名向量与函数遮蔽问题
if (file.exists("data/DSigDB/LINCS_L1000_Chem_Pert_down.gmt")) {
  saveRDS(clusterProfiler::read.gmt("data/DSigDB/LINCS_L1000_Chem_Pert_down.gmt"),
          "data/DSigDB/L1000_down.rds")
  message("L1000_down.rds 已生成")
}
if (file.exists("data/DSigDB/LINCS_L1000_Chem_Pert_up.gmt")) {
  saveRDS(clusterProfiler::read.gmt("data/DSigDB/LINCS_L1000_Chem_Pert_up.gmt"),
          "data/DSigDB/L1000_up.rds")
  message("L1000_up.rds 已生成")
}
message("全部完成。失败的文件请按上方直链用浏览器手动下载。")
