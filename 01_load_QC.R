# ============================================================================
# 01_load_QC.R —— Aim1: 三个 scRNA-seq 数据集读取、质控
# 前置: 按 data/GSM_sample_map.csv 解压整理 RAW.tar 到 raw/ 子目录
# ============================================================================
library(Seurat); library(tidyverse)

## raw/<cond>/ 若含子目录(多样本条件), 自动进入第一个子目录
find_sample_dir <- function(cond_dir) {
  fs <- list.files(cond_dir, full.names = TRUE)
  sub <- fs[file.info(fs)$isdir]
  if (length(sub) > 0) sub[1] else cond_dir
}

load10x <- function(data.dir, project, meta = list()) {
  obj <- CreateSeuratObject(Read10X(find_sample_dir(data.dir)), project = project,
                            min.cells = 3, min.features = 200)
  obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
  for (k in names(meta)) obj[[k]] <- meta[[k]]
  obj
}

## GSE195940: 骨膜 Prx1Cre;mTmG, 稳态 vs 骨折 d3
perio_ss <- load10x("data/GSE195940/raw/ss/", "Periosteum_SS",
  meta = list(tissue="periosteum", time="d0", injury="none", age="adult", cre="Prx1Cre"))
perio_d3 <- load10x("data/GSE195940/raw/d3/", "Periosteum_d3",
  meta = list(tissue="periosteum", time="d3", injury="fracture", age="adult", cre="Prx1Cre"))

## GSE164573: 骨骼肌 Prx1-GFP+, d0 / 骨折 d3,d5 / 复合伤 d3,d5
load_muscle <- function(dir, time, injury) {
  load10x(dir, paste0("Muscle_", time, "_", injury),
    meta = list(tissue="muscle", time=time, injury=injury, age="adult", cre="Prx1Cre"))
}
mus_d0 <- load_muscle("data/GSE164573/raw/d0/",  "d0", "none")
mus_f3 <- load_muscle("data/GSE164573/raw/fd3/", "d3", "fracture")
mus_f5 <- load_muscle("data/GSE164573/raw/fd5/", "d5", "fracture")
mus_p3 <- load_muscle("data/GSE164573/raw/pd3/", "d3", "polytrauma")
mus_p5 <- load_muscle("data/GSE164573/raw/pd5/", "d5", "polytrauma")

## GSE138689: 骨髓 Lepr+ SSPC, 稳态 vs 骨折 d7
bm_ss <- load10x("data/GSE138689/raw/ss/", "BM_SS",
  meta = list(tissue="BM", time="d0", injury="none", age="adult", cre="LeprCre"))
bm_d7 <- load10x("data/GSE138689/raw/d7/", "BM_d7",
  meta = list(tissue="BM", time="d7", injury="fracture", age="adult", cre="LeprCre"))

## QC
objs <- list(perio_ss, perio_d3, mus_d0, mus_f3, mus_f5, mus_p3, mus_p5,
             bm_ss, bm_d7)
objs <- lapply(objs, function(x)
  subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                  percent.mt < 15 & nCount_RNA < 50000))
saveRDS(objs, "results/qc_filtered_objects.rds")
cat("QC done:", length(objs), "objects\n")
