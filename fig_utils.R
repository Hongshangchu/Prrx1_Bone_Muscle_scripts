# ============================================================================
# fig_utils.R —— 统一图形导出: PDF (矢量) + 600 dpi TIFF (期刊投稿要求)
# ----------------------------------------------------------------------------
# 用法 (所有绘图脚本开头先 source):
#   source("scripts/fig_utils.R")
#
# 单页图:
#   save_fig(ggplot对象, "results/figures/FigX.pdf", width = 8, height = 6)
#     -> FigX.pdf + FigX.tiff (600 dpi, LZW 无损压缩)
#
# 多页图 (一个 PDF 含多页):
#   save_fig_multi(图列表, "results/figures/FigX.pdf", width, height)
#     -> FigX.pdf (多页) + FigX_p1.tiff / FigX_p2.tiff ... (每页一张 600 dpi TIFF)
#
# 规则: 全项目图形一律经 save_fig / save_fig_multi 导出,
#       不再直接调用 pdf() / ggsave() 输出到 figures 目录,
#       保证每个 PDF 都有对应 600 dpi TIFF。
# ============================================================================

FIG_DPI <- 600  # 期刊要求的分辨率

## TIFF 设备工厂: cairo 渲染 (无 X11 的服务器也能出图), LZW 无损压缩
.tiff_device <- function(dpi = FIG_DPI) {
  force(dpi)
  function(filename, width, height, ...) {   # ... 吸收 ggsave 传入的 bg 等参数
    grDevices::tiff(filename, width = width, height = height, units = "in",
                    res = dpi, compression = "lzw",
                    type = if (capabilities("cairo")) "cairo"
                           else getOption("bitmapType", "Xlib"),
                    ...)
  }
}

## ---- 单页图: PDF + TIFF ----
save_fig <- function(plot, filename, width, height, dpi = FIG_DPI, ...) {
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE))
    stop("save_fig: 文件名必须以 .pdf 结尾 -> ", filename)
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  tiff_name <- sub("\\.pdf$", ".tiff", filename, ignore.case = TRUE)
  ## 矢量 PDF (cairo 嵌入字体)
  ggplot2::ggsave(filename, plot, width = width, height = height,
                  device = grDevices::cairo_pdf, ...)
  ## 600 dpi TIFF
  ggplot2::ggsave(tiff_name, plot, width = width, height = height,
                  device = .tiff_device(dpi), dpi = dpi, ...)
  invisible(c(pdf = filename, tiff = tiff_name))
}

## ---- 多页图: 多页 PDF + 每页一张 TIFF ----
save_fig_multi <- function(plots, filename, width, height, dpi = FIG_DPI, ...) {
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE))
    stop("save_fig_multi: 文件名必须以 .pdf 结尾 -> ", filename)
  if (!is.list(plots) || length(plots) == 0)
    stop("save_fig_multi: plots 必须是非空列表")
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  stem <- sub("\\.pdf$", "", filename, ignore.case = TRUE)

  ## 多页 PDF (与旧 pdf() 行为一致: 每个元素一页)
  grDevices::cairo_pdf(filename, width = width, height = height, onefile = TRUE)
  for (p in plots) print(p)
  grDevices::dev.off()

  ## 每页一张 TIFF: FigX_p1.tiff, FigX_p2.tiff, ...
  tiff_names <- sprintf("%s_p%d.tiff", stem, seq_along(plots))
  for (i in seq_along(plots)) {
    ggplot2::ggsave(tiff_names[i], plots[[i]], width = width, height = height,
                    device = .tiff_device(dpi), dpi = dpi, ...)
  }
  invisible(c(pdf = filename, tiff = tiff_names))
}
