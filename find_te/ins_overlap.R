
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  print("Usage: Rscript overlap2.R <overlap.csv>")
  print("")
  print("Option:")
  print("    overlap.csv    The overlapped csv file within TE folder")
  stop("arguments failed")
}

f <- args[1]
file_png <- gsub(".csv$", ".png", f)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(hiseqr))
suppressPackageStartupMessages(library(patchwork))

df <- read.csv(f)
if(nrow(df) > 3) {
  message("more than 2 rows found, choose top 2")
  df <- df[1:3, ]
}

p_list <- lapply(df$fileB, function(i) {
  di <- dplyr::filter(df, fileB == i)
  hiseqr::overlap_p2(n1 = di$numA, n2 = di$numTE, n12 = di$BinA,
                     name1 = di$fileA, name2 = di$fileB,
                     height = 4, width = 4)
})
# patchwork
p <- patchwork::wrap_plots(p_list, ncol = 2) +
  patchwork::plot_annotation(
    title = "Overlap between TEseq and ONT,Illumina",
    tag_levels = "A"
  )
message(glue::glue("Save overlap to file: {file_png}"))
ggplot2::ggsave(file_png, p, width = 5, height = 4, units = "in", dpi = 300)
