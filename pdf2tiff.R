# ============================================================================
# pdf2tiff.R —— 批量把 results/figures 下已有 PDF 转成 600 dpi TIFF (可选)
# 依赖: 系统 poppler-utils 的 pdftoppm (Linux: apt/yum install poppler-utils)
# 注意: 此脚本仅用于补历史文件; 新生成的图请走 fig_utils.R 的 save_fig
# ============================================================================
fig_dir <- "results/figures"
pdfs <- list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE)
if (length(pdfs) == 0) stop("未找到 PDF: ", fig_dir)
stopifnot("未安装 pdftoppm (poppler-utils)" = nzchar(Sys.which("pdftoppm")))
for (f in pdfs) {
  stem <- sub("\\.pdf$", "", f)
  ## -r 600: 600 dpi; -tiff: 输出 TIFF; 多页 PDF 自动编号 -1, -2 ...
  system2("pdftoppm", c("-r", "600", "-tiff", shQuote(f), shQuote(stem)))
  message("转换完成: ", f)
}
