library(tidyverse)

# Collect all ichorCNA TF files
ichor_files <- snakemake@input[["ichor_tfs"]]
fragle_files <- snakemake@input[["fragle_csvs"]]

# Read and bind ichorCNA TFs
ichor_df <- map_dfr(ichor_files, read_tsv, show_col_types = FALSE)

# Read and bind Fragle TFs
fragle_df <- map_dfr(fragle_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  # Extract library and window from path
  parts <- str_match(f, "/fragle/([^/]+)\\.([^/]+)/[^/]+/Fragle\\.csv")
  tibble(
    library_id = parts[,2],
    window = parts[,3],
    pon = NA_character_,
    preset = NA_character_,
    method = "Fragle",
    tumor_fraction = df$Fragle[1],
    ploidy = NA_real_
  )
})

# Combine and write
bind_rows(ichor_df, fragle_df) %>%
  write_tsv(snakemake@output[["tsv"]])
