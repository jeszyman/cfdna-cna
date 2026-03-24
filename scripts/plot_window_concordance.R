library(tidyverse)

df <- read_tsv(snakemake@input[["tsv"]], show_col_types = FALSE) %>%
  mutate(tumor_fraction = as.numeric(tumor_fraction))

p <- ggplot(df, aes(x = window, y = tumor_fraction,
                     group = library_id, color = library_id)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~method, scales = "free_y") +
  labs(x = "Window", y = "Tumor Fraction",
       title = "Window Concordance: TF across fragment length windows") +
  theme_bw()

ggsave(snakemake@output[["pdf"]], p, width = 10, height = 6)
ggsave(snakemake@output[["png"]], p, width = 10, height = 6, dpi = 150)
