library(tidyverse)
library(pheatmap)

depth_files <- snakemake@input[["depth_files"]]

mat_list <- lapply(depth_files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  parts <- str_match(basename(dirname(f)), "^([^.]+)\\.([^.]+)\\.([^.]+)\\.([^.]+)$")
  lib <- parts[,2]
  setNames(df$copy.norm, paste0(df$chr, ":", df$start)) -> vals
  tibble(bin = names(vals), logR = vals, library_id = lib)
})

combined <- bind_rows(mat_list)
mat <- combined %>%
  pivot_wider(names_from = bin, values_from = logR) %>%
  column_to_rownames("library_id") %>%
  as.matrix()

pdf(snakemake@output[["pdf"]], width = 20, height = max(4, nrow(mat)))
pheatmap(mat, cluster_cols = FALSE, show_colnames = FALSE,
         main = "LogR Heatmap", color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
dev.off()

png(snakemake@output[["png"]], width = 2000, height = max(400, nrow(mat) * 80), res = 150)
pheatmap(mat, cluster_cols = FALSE, show_colnames = FALSE,
         main = "LogR Heatmap", color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
dev.off()
