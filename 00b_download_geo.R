# ============================================================================
# 00b_download_geo.R —— GEO 补充文件下载 (https + 重试; 失败时打印浏览器直链)
# 国内建议: 优先浏览器/迅雷手动下载, 本脚本检测到已存在文件会自动跳过
# ============================================================================

dl_geo <- function(gse, pattern, dest = "data", tries = 3) {
  dir.create(file.path(dest, gse), showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(dest, gse, paste0(gse, "_", pattern))
  if (file.exists(outfile) && file.size(outfile) > 1e6) {
    message("已存在, 跳过: ", outfile); return(invisible(TRUE))
  }
  subdir <- paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", subdir, "/", gse,
                "/suppl/", gse, "_", pattern)
  for (i in 1:tries) {
    message("尝试 ", i, "/", tries, ": ", url)
    ok <- tryCatch({
      download.file(url, outfile, method = "libcurl", mode = "wb", quiet = TRUE)
      file.exists(outfile) && file.size(outfile) > 1e6
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) { message("成功: ", outfile); return(invisible(TRUE)) }
    Sys.sleep(5)
  }
  message("!! 下载失败, 请浏览器手动下载:\n   ", url)
  invisible(FALSE)
}

dl_geo("GSE195940", "RAW.tar")
dl_geo("GSE164573", "RAW.tar")
dl_geo("GSE138689", "RAW.tar")

## 解压指引:
##   tar -xf data/GSE195940/GSE195940_RAW.tar -C data/GSE195940/
##   再按 data/GSM_sample_map.csv 将各 GSM 文件夹移入 raw/ss, raw/d3 等
