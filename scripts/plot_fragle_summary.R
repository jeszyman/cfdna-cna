library(tidyverse)

fragle_files <- snakemake@input[["fragle_csvs"]]

df <- map_dfr(fragle_files, function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  parts <- str_match(f, "/fragle/([^/]+)\\.([^/]+)/Fragle\\.csv")
  d %>% mutate(library_id = parts[,2], window = parts[,3])
})

p <- ggplot(df, aes(x = library_id, y = Fragle, fill = window)) +
  geom_col(position = "dodge") +
  labs(x = "Library", y = "Fragle ctDNA Fraction",
       title = "Fragle Predictions by Window") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(snakemake@output[["pdf"]], p, width = 10, height = 6)
ggsave(snakemake@output[["png"]], p, width = 10, height = 6, dpi = 150)
