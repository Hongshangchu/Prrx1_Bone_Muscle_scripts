# ============================================================================
# 00d_organize_geo.R —— 自动解包 RAW.tar 并按条件整理到 raw/<cond>/ 目录
# 输出: 每个样本一个文件夹, 内含标准 10X 三件套
#       (barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz)
# 同时生成 results/tables/geo_sample_organization.csv 供核对
# 无需联网。整理完成后可手动删除各 GSE 文件夹下的 extracted/ 省空间
# ============================================================================

## ---------- 工具函数 ----------
untar_once <- function(tarfile, exdir) {
  if (!file.exists(tarfile)) stop("找不到: ", tarfile)
  if (!dir.exists(exdir)) {
    message("解压: ", tarfile)
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    untar(tarfile, exdir = exdir)
  } else message("已解压, 跳过: ", exdir)
  invisible(exdir)
}

## 按样本名关键词分类 -> raw 子目录; 返回 NA=未知, "SKIP"=本分析不需要
classify_sample <- function(gse, prefix) {
  p <- tolower(prefix)
  if (gse == "GSE195940") {                    # 骨膜: 稳态 vs 骨折 d3
    if (grepl("ss|steady|p0|uninj", p)) return("ss")
    if (grepl("frac|d3|day3|injur", p)) return("d3")
  }
  if (gse == "GSE164573") {                    # 骨骼肌: d0/骨折/复合伤
    if (grepl("poly", p)) {
      if (grepl("d3|day3", p)) return("pd3")
      if (grepl("d5|day5", p)) return("pd5")
    } else if (grepl("frac", p)) {
      if (grepl("d3|day3", p)) return("fd3")
      if (grepl("d5|day5", p)) return("fd5")
    } else if (grepl("d0|uninj|naive|steady|ss|ctrl|control", p)) return("d0")
  }
  if (gse == "GSE138689") {                    # 骨髓 LeprCre (MCY=Lepr-Cre;tdTomato):
    if (grepl("frac", p))                     return("d7")   # MCY_Fracture -> bm_d7
    if (grepl("adult2", p))                   return("SKIP") # 第二重复, 可选加入
    if (grepl("adult", p) && !grepl("aging", p)) return("ss") # MCY_Adult -> bm_ss
    return("SKIP")                             # Aging/Rosi/Irradiation/IRRA/Prx1_BM 不用
  }
  NA
}

## 把解压出的散文件整理成 raw/<cond>/<样本>/ 标准三件套
organize <- function(gse) {
  exdir <- file.path("data", gse, "extracted")
  untar_once(list.files(file.path("data", gse), pattern = "_RAW\\.tar$",
                        full.names = TRUE)[1], exdir)

  files <- list.files(exdir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("(barcodes|features|genes|matrix)\\.(tsv|mtx)(\\.gz)?$",
                       basename(files))]
  if (length(files) == 0) stop(gse, ": 解压后未找到 10X 文件, 请检查 tar 内容")

  ## 按前缀分组 (去掉 _barcodes / _features / _genes / _matrix 后缀)
  prefix <- sub("_(barcodes|features|genes|matrix)\\.[^.]+(\\.gz)?$", "",
                basename(files))
  groups <- split(files, prefix)
  log <- list()

  for (pf in names(groups)) {
    cond <- classify_sample(gse, pf)
    if (is.na(cond)) {
      message("!! 无法分类: ", gse, " / ", pf, "  -> 请人工检查并加入 classify_sample()")
      next
    }
    if (cond == "SKIP") { message("跳过(不需要): ", pf); next }
    ## 单样本条件: 三件套直接放 raw/<cond>/; 多样本条件: raw/<cond>/<prefix>/
    cond_dir <- file.path("data", gse, "raw", cond)
    dest <- if (dir.exists(cond_dir) &&
                length(list.files(cond_dir, pattern = "\\.(tsv|mtx)\\.gz$")) > 0)
              file.path(cond_dir, pf) else cond_dir
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    for (f in groups[[pf]]) {
      b <- basename(f)
      newname <- if (grepl("barcodes", b)) "barcodes.tsv.gz"
                 else if (grepl("features|genes", b)) "features.tsv.gz"
                 else "matrix.mtx.gz"
      file.copy(f, file.path(dest, newname), overwrite = TRUE)
    }
    log[[length(log) + 1]] <- data.frame(GSE = gse, sample = pf, condition = cond,
                                         dir = dest, n_files = length(groups[[pf]]))
    message("整理: ", pf, " -> ", gse, "/raw/", cond, "/")
  }
  bind_rows(log)
}

## ---------- 主流程 ----------
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

log_all <- bind_rows(lapply(c("GSE195940", "GSE164573", "GSE138689"), organize))

write.csv(log_all, "results/tables/geo_sample_organization.csv", row.names = FALSE)
message("全部完成! 映射表: results/tables/geo_sample_organization.csv")
message("提示: 核对无误后可删除各 GSE 文件夹下的 extracted/ 与 *_RAW.tar 释放空间")
