library(tidyverse)

df <- read_tsv(snakemake@input[["tsv"]], show_col_types = FALSE)

p <- ggplot(df, aes(x = library_id, y = as.numeric(tumor_fraction),
                     color = method, shape = preset)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  facet_wrap(~window, scales = "free_x") +
  labs(x = "Library", y = "Tumor Fraction", title = "TF Comparison: ichorCNA vs Fragle") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(snakemake@output[["pdf"]], p, width = 12, height = 6)
ggsave(snakemake@output[["png"]], p, width = 12, height = 6, dpi = 150)
