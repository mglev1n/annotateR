#' annotate_rsids
#'
#' Annotates genome wide association study (GWAS) results (or other tabular data including chromosome, position, and allele data) with rsids from a reference dataset.
#' @param df data frame or tibble containing data for annotation
#' @param chr_col column containing chromosome data
#' @param pos_col column containing position data
#' @param effect_allele_col column containing effect alleles (used only for rsid matching, no allele harmonization is performed)
#' @param other_allele_col column containing effect alleles (used only for rsid matching, no allele harmonization is performed)
#' @param multi_thread enable multi-threaded processing. Use with caution if resource-limited (default = TRUE)
#' @param progress_bar enable progress bar (default = TRUE)
#' @param ref_dataset path to folder containing reference dataset in .parquet format, which will be rapidly loaded using the arrow::open_dataset() function
#' @return Tibble containing original dataset + additional column including rsids
#' @export

annotate_rsids <- function(df, chr_col = chr, pos_col = pos, effect_allele_col = EFFECT_ALLELE, other_allele_col = OTHER_ALLELE, multi_thread = TRUE, progress_bar = TRUE, ref_dataset = NA) {

  # Check for reference dataset
  if(is.na(ref_dataset)) {
    stop("Path to reference dataset is required")
  }

  # Rename columns for standardized processing
  df <- df %>%
    rename(chr = {{ chr_col }}, pos = {{ pos_col }}, EFFECT_ALLELE = {{ effect_allele_col }}, OTHER_ALLELE = {{ other_allele_col }})

  # Convert chromosome column to character to ensure compatibility with 1kg data
  if(class(df$chr) != "character") {
    df <- df %>% mutate(chr = as.character(chr))
  }



  # Annotation function
  .annotate <- function(dat, chrom) {

    p(message = glue("Reading chromosome {chrom}"), class = "sticky")
    # cli::cli_alert_info("Reading chromosome {chrom}")

    annotation_df <- open_dataset("../RawData/1kg_all_phase3_snps_markers/") %>%
      filter(pos %in% unique(dat$pos)) %>%
      filter(chr == chrom) %>%
      select(pos, ref, alt, rsid)  %>%
      # head() %>%
      collect()

    p(message = glue("Annotating chromosome {chrom}"), class = "sticky")
    # cli::cli_alert_info("Annotating chromosome {chrom}")

    dat %>%
      dtplyr::lazy_dt() %>%
      # dtplyr::lazy_dt(key_by = c(pos)) %>%
      left_join(annotation_df, by = c("pos" = "pos")) %>%
      as_tibble() %>%
      mutate(rsid = case_when(
        (EFFECT_ALLELE == ref & OTHER_ALLELE == alt) | (OTHER_ALLELE == ref & EFFECT_ALLELE == alt) ~ rsid,
        TRUE ~ NA_character_
      )) %>%
      select(-ref, -alt) %>%
      unique()
  }

  # Set multithreading mode; need to keep number of processes low to ensure enough memory available
  if(multi_thread) {
    plan(multisession, workers = availableCores()/4)
  } else {
    plan(sequential)
  }

  with_progress({
    p <- progressor(along = 1:(length(unique(df$chr))*2+1))
    df %>%
      group_nest(chr) %>%
      arrange(as.numeric(chr)) %>%
      mutate(annotated = future_map2(data, chr, ~.annotate(dat = .x, chrom = .y), .options = furrr_options(seed = TRUE))) %>%
      select(-data) %>%
      unnest(annotated) %>%
      rename("{{ chr_col }}" := chr, "{{ pos_col }}" := pos, "{{ effect_allele_col }}" := EFFECT_ALLELE, "{{ other_allele_col }}" := OTHER_ALLELE)}, enable = progress_bar, handlers = handler_progress())

}
