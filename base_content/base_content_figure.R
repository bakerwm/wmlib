  #!/usr/bin/env Rscript
  ## create base-content figure
  ## Author: Ming Wang
  ## Date: 2018-05-20

  args = commandArgs(trailingOnly = TRUE)

  # test arguments
  if (length(args) != 2) {
    stop("Rscript base_content_figure.R <input.txt> <out.png>", call. = FALSE)
  }
  indir <- args[1] # statdir
  outname <- args[2] # name

  f1 <- file.path(indir, "base_content.txt")
  f2 <- file.path(indir, "length_distribution.txt")

  if (! file.exists(f1) & file.exists(f2)) {
    stop("file not exists: base_content.txt, length_distribution.txt")
  }

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))

  ##----------------------------------------------------------------------------##
  df <- read.table(f1, header = TRUE, sep = " ")
  df[is.na(df)] <- 0 # NA to 0
  df <- tidyr::gather(df, base, count, -position)
  df <- filter(df, base %in% c("A", "C", "G", "T")) %>%
    mutate(base = factor(base, levels = c("A", "C", "G", "T"))) %>%
    group_by(position) %>%
    mutate(freq = count / sum(count))

  p1 <- ggplot(df, aes(position, freq, fill = base)) +
    geom_bar(position = "fill", stat = "identity") +
    xlab("Position in read (bp)") +
    ylab("Per base sequence content (%)") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    scale_fill_manual(values = c("green3", "blue3", "grey30", "red3")) +
    theme_classic()

  p2 <- ggplot(df, aes(position, freq, color = base)) +
    geom_line(size = .5) +
    #geom_bar(position = "fill", stat = "identity") +
    xlab("Position in read (bp)") +
    ylab("Per base sequence content (%)") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    scale_color_manual(values = c("green3", "blue3", "grey30", "red3")) +
    theme_classic()

  ##----------------------------------------------------------------------------##
  df2 <- read.table(f2, header = FALSE, sep = " ")
  names(df2) <- c("length", "count")
  # if only one record in df2
  if (nrow(df2) == 1) {
    n <- df2[1, "length"]
    df2x <- data.frame(length = c(n - 1, n + 1),
                       count = 0)
    df2 <- rbind(df2, df2x)
  }

  p3 <- ggplot(df2, aes(length, count)) +
    geom_line(size = .5, color = "red3") +
    xlab("Length of reads (nt)") +
    ylab("Number of reads") +
    theme_classic()

  ##----------------------------------------------------------------------------##
  png1 <- file.path(indir, paste0(outname, ".base_content.bar.png"))
  png2 <- file.path(indir, paste0(outname, ".base_content.line.png"))
  png3 <- file.path(indir, paste0(outname, ".length_distribution.line.png"))
  ggsave(png1, p1, width = 5, height = 3, units = "in")
  ggsave(png2, p2, width = 5, height = 3, units = "in")
  ggsave(png3, p3, width = 5, height = 3, units = "in")
  ##----------------------------------------------------------------------------##

  ## EOF
