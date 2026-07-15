# =============================================================================
# 00_download_gbif.R
#
# REQUIRED CITATION
#   GBIF.org (20 February 2026) GBIF Occurrence Download
#   https://doi.org/10.15468/dl.n295k4
# =============================================================================

library(rgbif)
library(here)

# --- Download identifiers -----------------------------------------------------
GBIF_DOI <- "10.15468/dl.n295k4"
GBIF_KEY <- "0031909-260208012135463"

# Filter used in the original query, for reference:
#   TaxonKey           = Squamata
#   HasCoordinate      = TRUE
#   HasGeospatialIssue = FALSE
#   Geometry           = POLYGON((-165.09376 -59.05944, -29.42038 -59.05944,
#                                 -29.42038 75.94475, -170.19808 76.39084,
#                                 -165.09376 -59.05944))

dest_dir <- here::here("..", "data", "raw")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
dest_csv <- file.path(dest_dir, "gbif_occurrences.csv")

if (file.exists(dest_csv)) {
  message("File already present at: ", dest_csv, "\nNothing to do.")
} else {

  message("Downloading GBIF request ", GBIF_KEY, " (~479 MB zipped).")
  message("This may take several minutes.")

  # occ_download_get fetches an existing download by key. No credentials needed
  # for a download that already exists -- only for creating a new one.
  zipfile <- rgbif::occ_download_get(
    key       = GBIF_KEY,
    path      = tempdir(),
    overwrite = TRUE
  )

  message("Unzipping...")
  utils::unzip(zipfile, exdir = tempdir())

  tsv <- file.path(tempdir(), paste0(GBIF_KEY, ".csv"))
  if (!file.exists(tsv)) {
    # internal name may vary; take the single extracted .csv
    cands <- list.files(tempdir(), pattern = "\\.csv$", full.names = TRUE)
    stopifnot(length(cands) == 1)
    tsv <- cands
  }

  file.copy(tsv, dest_csv, overwrite = TRUE)
  message("Done: ", dest_csv)
}

# --- Provenance record --------------------------------------------------------
writeLines(c(
  "GBIF Occurrence Download",
  paste0("DOI:  https://doi.org/", GBIF_DOI),
  paste0("Key:  ", GBIF_KEY),
  "Date: 20 February 2026",
  "Records: 3,780,267 from 1,225 published datasets",
  "Format: simple tab-separated values (TSV)",
  "",
  "Cite as:",
  paste0("  GBIF.org (20 February 2026) GBIF Occurrence Download ",
         "https://doi.org/", GBIF_DOI)
), file.path(dest_dir, "gbif_download_doi.txt"))

# --- Retention note -----------------------------------------------------------
# GBIF keeps the TSV for six months (until 20 August 2026). Downloads cited by
# DOI are detected and retained indefinitely. If this script fails with a 404,
# the download expired: rebuild the query from the filter above with
# rgbif::occ_download(), and record the NEW DOI here. Results will not be
# identical.
