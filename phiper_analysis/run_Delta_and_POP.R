# load required packages
library(pak)
library(phiperio)
library(phiper)
library(svglite)
library(rlang)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(purrr)
library(locfdr)
library(rsvg)
library(openxlsx)
library(purrr)
library(locfdr)
require(extrafont)
# phip_use_montserrat()
# RUn this for installing phiper from source
# install.packages("/Users/gabrielinnocenti/Desktop/phd_tmp/phiper/phiper_new_version/phiper/", repos = NULL, type = "source")

setwd("/Users/gabrielinnocenti/Desktop/phd_tmp/Projects/CORSA/")

# Exist matrix
exist_path <- file.path(getwd(), "Final_output/sampleset_bin_more_informative.csv")
# Counts matrix
hit_path <- file.path(getwd(), "Final_output/sampleset_count_more_informative.csv")
# Input counts
input_path <- file.path(getwd(), "input_counts/inputATC2_L75_M2.csv")
# read.csv(input_path)
# Samples metadata
samples_path <- file.path(getwd(), "samples.csv")

ps <- convert_legacy(
  exist_file = exist_path,
  input_file = input_path,
  hit_file = hit_path,
  samples_file = samples_path,
  extra_cols = c("age", "sex", "group"),
  peptide_library  = FALSE,
  materialise_table = FALSE
)

# Pull into memory
ps_data_long_df <- ps$data_long %>% collect()
ps_data_long_df$subject_id <- ps_data_long_df$sample_id
ps$data_long <- ps_data_long_df
# 
# # Rename columns
# colnames(ps_data_long_df) <- c("peptide_id", "sample_id", "exist", "inut_count", "hit_count", "group_char", "sex", "age")



plot_enrichment_counts(
  ps,
  group_cols        = c("group"),
  group_interaction = FALSE,
  interaction_only  = FALSE,
  annotation_size   = 3
)

# Prepare pop and delta frameworkwar
set.seed(632961)

# parse command-line arguments for optional parameters
N_CORES <- 14
LOG <- TRUE
LOG_FILE <- NULL
MAX_GB <- 20
CMP_INDEX <- NA_integer_
CMP_KEY <- NULL
args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("=", arg)) {
    parts <- strsplit(arg, "=")[[1]]
    key <- parts[1]
    value <- parts[2]
    if (key == "N_CORES") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) N_CORES <- val_num
    } else if (key == "MAX_GB") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) MAX_GB <- val_num
    } else if (key == "LOG") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        LOG <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        LOG <- FALSE
      }
    } else if (key == "LOG_FILE") {
      LOG_FILE <- sub("^['\\\"]|['\\\"]$", "", value)
    } else if (key %in% c("CMP_INDEX", "COMPARISON_INDEX")) {
      val_num <- suppressWarnings(as.integer(value))
      if (!is.na(val_num)) CMP_INDEX <- val_num
    } else if (key %in% c("CMP_KEY", "COMPARISON_KEY")) {
      CMP_KEY <- sub("^['\\\"]|['\\\"]$", "", value)
    }
  }
}


# create base results directory and save peptide library
# dir.create("results_maxmean_strat_bin_2perc", recursive = TRUE, showWarnings = FALSE)
phiperio::get_peptide_library() %>%
  collect() %>%
  as.data.frame() -> peplib
# Relabel in a single label all the CagA antigens
# peplib$Description[
#   peplib$Description %in% c(
#     "CagA",
#     "Cytotoxicity-associated immunodominant antigen",
#     "Cag pathogenicity island protein 7"
#   )
# ] <- "Cytotoxicity-associated immunodominant antigen (CagA)"

# Load improved annotations with VFDB naming
vfdb_annotations <- read.csv('/Users/gabrielinnocenti/Desktop/phd_tmp/Projects/CORSA/Current_output_all/vfdb_custom_name_mapping.csv')
colnames(vfdb_annotations)
peplib <- peplib %>%
  left_join(vfdb_annotations, by = c("peptide_id" = "peptide_name")) %>%
  mutate(
    Description = ifelse(
      !is.na(custom_name),
      custom_name,
      Description
    )
  ) 


# FUll info on proteins
prot_to_pept <- readRDS('/Users/gabrielinnocenti/Downloads/protein-peptide_seq_agilent-twist.rds')
saveRDS(peplib, file = "results/peptide_library.rds")

# define list of comparisons and longitudinal flags

source("/Users/gabrielinnocenti/Downloads/zzz.R")

if (is.na(CMP_INDEX)) {
  env_idx <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (nzchar(env_idx)) {
    val_num <- suppressWarnings(as.integer(env_idx))
    if (!is.na(val_num)) CMP_INDEX <- val_num
  }
}

if (!is.null(CMP_KEY) && !is.na(CMP_INDEX)) {
  stop("Provide only one of CMP_KEY or CMP_INDEX (or SLURM_ARRAY_TASK_ID).")
}

if (!is.null(CMP_KEY)) {
  cmp_keys <- vapply(comparisons, paste, collapse = "|", FUN.VALUE = character(1))
  match_idx <- match(CMP_KEY, cmp_keys)
  if (is.na(match_idx)) {
    stop("CMP_KEY not found in comparisons list: ", CMP_KEY)
  }
  comparisons <- list(comparisons[[match_idx]])
  longitudinal <- longitudinal[match_idx]
} else if (!is.na(CMP_INDEX)) {
  if (CMP_INDEX < 1L || CMP_INDEX > length(comparisons)) {
    stop("CMP_INDEX out of range: ", CMP_INDEX, " (1..", length(comparisons), ")")
  }
  comparisons <- list(comparisons[[CMP_INDEX]])
  longitudinal <- longitudinal[CMP_INDEX]
}

# columns always to keep in saved data
base_cols <- c("sample_id", "peptide_id", "group", "exist", "fold_change")

# helper function to safely save an RDS file after creating directories
save_rds_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
}

# helper function: select base columns (and any extras), drop NA rows, collect and save to RDS
make_and_save <- function(data, out_path, extra_vars = NULL, drop_na = TRUE) {
  cols_to_select <- unique(c(base_cols, extra_vars %||% character(0)))
  avail <- names(data)
  missing <- setdiff(cols_to_select, avail)
  if (length(missing)) {
    message("Skipping missing columns: ", paste(missing, collapse = ", "))
  }
  df <- data %>% dplyr::select(dplyr::any_of(cols_to_select))
  if (isTRUE(drop_na)) {
    df <- df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(intersect(cols_to_select, names(.))), ~ !is.na(.)))
  }
  df <- df %>% dplyr::collect()
  save_rds_safe(df, out_path)
  invisible(df)
}

# set up parallel processing for DELTA framework analysis
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
options(future.globals.maxSize = MAX_GB * 1024^3, future.scheduling = 1)
original_plan <- future::plan()
if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_CORES)
} else {
  if (N_CORES > 1L) {
    future::plan(future::multicore, workers = N_CORES)
  } else {
    future::plan(future::sequential)
  }
}

# loop through each comparison and perform the same analysis
for (cmp_idx in seq_along(comparisons)) {
  cmp <- comparisons[[cmp_idx]]
  print(cmp)
  var1 <- cmp[1]
  var2 <- cmp[2]
  is_longitudinal <- isTRUE(longitudinal[cmp_idx])
  message("Running comparison: ", var1, " vs ", var2)
  label_dir <- paste(var1, "vs", var2, sep = "_")
  paired_col <- if (is_longitudinal) "subject_id" else NULL
  # create output directory for this comparison
  out_dir <- file.path("results_maxmean_cagA_twist_reannotated_filtered_peptides", label_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # filter data for the two groups (always OR), then, if paired, keep only IDs with both groups
  # ps_cmp <- ps %>%
  #   filter((!!sym(var1) == 1L) | (!!sym(var2) == 1L)) %>%
  #   mutate(
  #     group_char  = if_else(!!sym(var1) == 1L, var1, var2),
  #     group_dummy = if_else(!!sym(var1) == 1L, 1L, 0L)
  #   )
  # # Fix
  ps_cmp <- ps %>%
    filter(group %in% c(var1, var2)) %>%  # <- this is the fix
    mutate(
      group_char  = group,                   # just copy the existing group
      group_dummy = if_else(group == var1, 1L, 0L)
    )
  
  if (!is.null(paired_col)) {
    ps_cmp <- ps_cmp %>%
      group_by(!!sym(paired_col)) %>%
      filter(any(!!sym(var1) == 1L) & any(!!sym(var2) == 1L)) %>%
      ungroup()
  }
  
  # save the filtered dataset for this comparison (keep both pairing candidates)
  make_and_save(
    data       = ps_cmp,
    out_path   = file.path(out_dir, paste0(label_dir, "_data.rds")),
    extra_vars = unique(c(paired_col, "subject_id"))
  )
  # ------------------ POP framework ----------------------
  data_frameworks <- readRDS(file.path(out_dir, paste0(label_dir, "_data.rds")))

  extract_tbl <- function(obj) {
    if (is.data.frame(obj)) {
      return(tibble::as_tibble(obj))
    }
    for (nm in c("data", "table", "tbl", "df", "result", "results")) {
      if (!is.null(obj[[nm]])) {
        return(tibble::as_tibble(obj[[nm]]))
      }
    }
    out <- try(tibble::as_tibble(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
    out <- try(as.data.frame(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(tibble::as_tibble(out))
    }
    stop("Cannot extract a data table from the POP result object.")
  }

  downsample_for_static <- function(df, prop = 0.1, seed = 1L) {
    if (is.null(df) || !nrow(df)) {
      return(df)
    }
    n <- nrow(df)
    size <- max(1L, floor(n * prop))
    if (size >= n) {
      return(df)
    }
    set.seed(seed)
    dplyr::slice_sample(df, n = size)
  }

  dir.create(file.path(out_dir, "POP_framework"), recursive = TRUE, showWarnings = FALSE)

  paired_prev <- if (is_longitudinal) "subject_id" else FALSE

  prev_res_pep <- phiper::compute_pop(
    x                 = data_frameworks,
    group_cols        = "group",
    rank_cols         = "peptide_id",
    paired            = paired_prev,
  )
  pep_tbl <- extract_tbl(prev_res_pep)
  write.csv(pep_tbl, file.path(out_dir, "POP_framework", "single_peptide.csv"))

  ranks_tax <- c("phylum", "class", "order", "family", "genus", "species")
  prev_res_rank <- phiper::compute_pop(
    x                 = data_frameworks,
    group_cols        = "group",
    rank_cols         = ranks_tax,
    paired            = paired_prev,
    peptide_library   = peplib,
  )
  rank_tbl <- extract_tbl(prev_res_rank)
  write.csv(rank_tbl, file.path(out_dir, "POP_framework", "taxa_ranks.csv"))

  ranks_combined <- c(ranks_tax, "peptide_id")
  plots_dir <- file.path(out_dir, "POP_framework", "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  for (rank_name in ranks_combined) {
    rank_chr <- as.character(rank_name)
    out_name <- file.path(plots_dir, rank_chr)
    df_rank <- if (rank_chr == "peptide_id") {
      pep_tbl
    } else {
      rank_tbl %>% filter(rank == rank_chr)
    }
    df_rank_static <- downsample_for_static(df_rank, prop = 1, seed = 1L)
    p_static <- scatter_static(
      df = df_rank_static,
      rank = rank_chr,
      xlab = df_rank$group1[1],
      ylab = df_rank$group2[1],
      point_size = 2,
      jitter_width_pp = 0.15,
      jitter_height_pp = 0.15,
      point_alpha = 0.85,
      font_size = 12
    ) +
      ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
        text        = ggplot2::element_text(family = "Montserrat")
      )


    p_inter <- scatter_interactive(
      df = df_rank,
      rank = rank_chr,
      xlab = df_rank$group1[1],
      ylab = df_rank$group2[1],
      peplib = peplib,
      point_size = 10,
      jitter_width_pp = 0.25,
      jitter_height_pp = 0.25,
      point_alpha = 0.85,
      font_size = 12
    )

    p_inter <- plotly::layout(
      p_inter,
      autosize = TRUE,
      margin   = list(l = 70, r = 30, t = 10, b = 70),
      xaxis    = list(range = c(-2, 102), automargin = TRUE),
      yaxis    = list(range = c(-2, 102), automargin = TRUE)
    )
    htmlwidgets::saveWidget(
      p_inter,
      file = paste0(out_name, "_interactive.html"),
      selfcontained = FALSE
    )  #   p_static <- scatter_static(
  #     df = df_rank_static,
  #     rank = rank_chr,
  #     xlab = df_rank$group1[1],
  #     ylab = df_rank$group2[1],
  #     point_size = 2,
  #     jitter_width_pp = 0.15,
  #     jitter_height_pp = 0.15,
  #     point_alpha = 0.85,
  #     font_size = 12
  #   ) +
  #     ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
  #     ggplot2::theme(
  #       plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
  #       text        = ggplot2::element_text(family = "Montserrat")
  #     )
  # 
  # 
  #   p_inter <- scatter_interactive(
  #     df = df_rank,
  #     rank = rank_chr,
  #     xlab = df_rank$group1[1],
  #     ylab = df_rank$group2[1],
  #     peplib = peplib,
  #     point_size = 10,
  #     jitter_width_pp = 0.25,
  #     jitter_height_pp = 0.25,
  #     point_alpha = 0.85,
  #     font_size = 12
  #   )
  # 
  #   p_inter <- plotly::layout(
  #     p_inter,
  #     autosize = TRUE,
  #     margin   = list(l = 70, r = 30, t = 10, b = 70),
  #     xaxis    = list(range = c(-2, 102), automargin = TRUE),
  #     yaxis    = list(range = c(-2, 102), automargin = TRUE)
  #   )
  #   htmlwidgets::saveWidget(
  #     p_inter,
  #     file = paste0(out_name, "_interactive.html"),
  #     selfcontained = FALSE
  #   )
  }
  # ------------------ DELTA framework ---------------------
  dir.create(file.path(out_dir, "DELTA_framework"), recursive = TRUE, showWarnings = FALSE)
  
  peplib[] <- lapply(peplib, as.character)
  log_file_current <- if (is.null(LOG_FILE)) {
    file.path(out_dir, "DELTA_framework", "log.txt")
  } else {
    LOG_FILE
  }
  
  stat_mode <- if (is_longitudinal) {
    "srlr_paired"
  } else {
    "srlr"
  }
  paired_by <- if (is_longitudinal) paired_col else NULL
  
  res <- phiper::compute_delta(
    x = data_frameworks,
    exist_col = "exist",
    rank_cols = c("genus", "species", "Description"
      # "is_auto", "is_infect", "is_EBV", "is_toxin", "is_PNP", "is_EM",
      # "is_MPA", "is_patho", "is_probio", "is_IgA", "is_flagellum", "is_allergens",
      # "anno_other", "anno_eucaryotic_pathogens", "anno_animals",
      # "anno_fly_worm_bee_cocroach_mite_mosquito", "anno_plant_grass",
      # "anno_is_fungi", "anno_fungi_type", "anno_is_food",
      # "anno_food_subcategory", "anno_food_item", "anno_bacteria_archaea",
      # "anno_is_homo_sapiens", "anno_viruses_bacteriophage", "anno_is_lacto_phage"
    ),
    group_cols = "group",
    peptide_library = peplib,
    B_permutations = 150000L,
    weight_mode = "equal",
    stat_mode = stat_mode,
    aggregate_stat = "maxmean",
    winsor_z = Inf,
    # strat_bins = 0.02,
    rank_feature_keep = list(
      phylum = NULL, Description = NULL, genus = NULL, species = NULL,
      is_auto = "TRUE", is_infect = "TRUE", is_EBV = "TRUE", is_toxin = "TRUE", is_PNP = "TRUE", is_EM = "TRUE",
      is_MPA = "TRUE", is_patho = "TRUE", is_probio = "TRUE", is_IgA = "TRUE", is_flagellum = "TRUE", is_allergens = "TRUE",
      anno_is_fungi = "TRUE", anno_is_food = "TRUE", anno_is_homo_sapiens = "TRUE", anno_is_lacto_phage = "TRUE"
    ),
    log = LOG,
    log_file = log_file_current,
    paired_by = paired_by
  )
  res <- as.data.frame(res)
  
  write.csv(res, file = file.path(out_dir, "DELTA_framework", "delta_table.csv"))
  tax_ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  std_idx <- res$rank %in% tax_ranks
  res$feature[!std_idx] <- as.character(res$rank[!std_idx])
  res$rank <- "all"
  
  # Open SVG device
  svglite(
    filename = file.path(out_dir, "DELTA_framework", "uncorrected_significant_static_forestplot.svg"),
    width = 30/2.54,   # convert cm to inches
    height = 30/2.54,
    bg = "white"
  )
  
  p_forest_unc <- phiper::forestplot(
    results_tbl = res,
    rank_of_interest = "all",
    use_diverging_colors = TRUE,
    filter_significant = "p_perm",
    left_label = paste0("More in ", df_rank$group1[1]),
    right_label = paste0("More in ", df_rank$group2[1]),
    label_vjust = -0.9,
    y_pad = 0.3,
    label_x_gap_frac = -0.3,
    statistic_to_plot = "T_stand"
  )
  print(p_forest_unc)
  dev.off()
  
  p_inter <- phiper::forestplot_interactive(
    results_tbl = res,
    rank_of_interest = "all",
    statistic_to_plot = "T_stand",
    use_diverging_colors = TRUE,
    filter_significant = "p_perm",
    left_label = paste0("More in ", df_rank$group1[1]),
    right_label = paste0("More in ", df_rank$group2[1]),
    arrow_length_frac = 0.35,
    label_x_gap_frac = -0.3,
    label_y_offset = -0.9
  )$plot
  htmlwidgets::saveWidget(
    p_inter,
    file = file.path(out_dir, "DELTA_framework", "uncorrected_significant_interactive_forestplot.html"),
    selfcontained = FALSE
  )
  
  # ------------------ plot interesting features ---------
  # Always include these features in the "interesting" set (even if not significant)
  always_keep <- c("")
  
  res_filtered <- res %>%
    dplyr::mutate(
      .force_keep = (.data$feature %in% always_keep) | (.data$rank %in% always_keep)
    ) %>%
    dplyr::filter(
      .force_keep | (.data$p_perm < 0.05)
    ) %>%
    dplyr::arrange(
      dplyr::desc(.force_keep), # optional: keep forced ones at the top
      dplyr::desc(.data$T_obs_stand)
    ) %>%
    dplyr::select(-.force_keep)
  
  write.csv(
    res_filtered,
    file      = file.path(out_dir, "DELTA_framework", "delta_table_interesting.csv"),
    row.names = FALSE
  )
  
  tax_cols <- intersect(tax_ranks, names(peplib))
  
  special_features <- c(
    "is_IEDB_or_cntrl", "is_auto", "is_infect", "is_EBV",
    "is_toxin", "is_PNP", "is_EM", "is_MPA", "is_patho",
    "is_probio", "is_IgA", "is_flagellum", "signalp6_slow",
    "is_topgraph_new", "is_allergens",
    "anno_is_fungi", "anno_is_food", "anno_is_homo_sapiens", "anno_is_lacto_phage"
  )
  
  get_binary_and_ids <- function(feature, peplib, tax_cols,
                                 peptide_col = "peptide_id") {
    if (feature %in% special_features) {
      vals <- peplib[[feature]]
      present <- !is.na(vals) & as.logical(vals)
    } else if (feature %in% names(peplib)) {
      vals <- peplib[[feature]]
      present <- as.logical(vals)
      present[is.na(present)] <- FALSE
    } else {
      if (length(tax_cols) == 0L) {
        stop("No taxonomic columns found in peptide library.")
      }
      matches <- lapply(tax_cols, function(col) {
        vals <- peplib[[col]]
        !is.na(vals) & vals == feature
      })
      present <- Reduce(`|`, matches)
    }
    
    peptide_ids <- as.character(peplib[[peptide_col]][present])
    peptide_ids <- unique(peptide_ids)
    list(present = present, peptide_ids = peptide_ids)
  }
  
  res_with_pep <- res_filtered %>%
    mutate(
      match_info = map(feature, ~ get_binary_and_ids(.x, peplib, tax_cols)),
      binary_in_peplib = map(match_info, "present"),
      peptide_ids = map(match_info, "peptide_ids"),
      pep_tbl_subset = map(peptide_ids, ~ pep_tbl %>% filter(feature %in% .x))
    ) %>%
    select(-match_info)
  
  dir.create(file.path(out_dir, "DELTA_framework", "interesting_features"), recursive = TRUE, showWarnings = FALSE)
  
  add_background_static <- function(p, bg,
                                    size = 0.8,
                                    alpha = 0.12,
                                    color = "#808080") {
    if (is.null(bg) || !nrow(bg)) {
      return(p)
    }
    
    bg_layer <- ggplot2::geom_point(
      data = bg,
      mapping = ggplot2::aes(x = percent1, y = percent2),
      inherit.aes = FALSE,
      color = color,
      size = size,
      alpha = alpha,
      show.legend = FALSE
    )
    
    p$layers <- c(list(bg_layer), p$layers)
    p
  }
  
  plot_feature_all <- function(feature_name, group1, group2, feature_data, out_dir) {
    if (is.null(feature_data) || nrow(feature_data) == 0L) {
      message("Skipping ", feature_name, " (no peptide data)")
      return(invisible(NULL))
    }
    safe_name <- gsub("[^A-Za-z0-9_-]+", "_", as.character(feature_name))
    base_dir <- file.path(out_dir, "DELTA_framework", "interesting_features")
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
    scatter_dir <- file.path(base_dir, "scatter")
    dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
    file_prefix <- file.path(base_dir, safe_name)
    scatter_prefix <- file.path(scatter_dir, safe_name)
    
    SHOW_BG <- TRUE
    BG_MAX_N_INTERACTIVE <- Inf
    BG_SEED <- 1L
    
    bg_df <- NULL
    if (isTRUE(SHOW_BG)) {
      bg_df <- pep_tbl
      
      if ("feature" %in% names(bg_df) && "feature" %in% names(feature_data)) {
        bg_df <- bg_df %>% dplyr::filter(!(feature %in% feature_data$feature))
      }
      
      # keep only unique (percent1, percent2) combos; drop the rest at random
      if (all(c("percent1", "percent2") %in% names(bg_df))) {
        set.seed(BG_SEED) # makes the random pick reproducible
        bg_df <- bg_df %>%
          dplyr::slice_sample(prop = 1) %>% # shuffle rows
          dplyr::distinct(percent1, percent2, .keep_all = TRUE)
      }
      
      if (nrow(bg_df) > BG_MAX_N_INTERACTIVE) {
        set.seed(BG_SEED)
        bg_df <- bg_df %>% dplyr::slice_sample(n = BG_MAX_N_INTERACTIVE)
      }
    }
    
    feature_data_static <- downsample_for_static(feature_data, prop = 1, seed = 1L)
    bg_df_static <- downsample_for_static(bg_df, prop = 1, seed = BG_SEED)
    
    ## ---------------- SCATTER STATIC ----------------
    svglite(
      filename = paste0(scatter_prefix, "_scatter_static.svg"),
      width = 30/2.54,    # convert cm to inches
      height = 30/2.54,
      bg = "white"
    )
    
    p_scatter <- scatter_static(
      df = feature_data_static,
      xlab = group1,
      ylab = group2,
      point_size = 2,
      point_alpha = 0.85,
      jitter_width_pp = 0.15,
      jitter_height_pp = 0.15,
      font_family = "Montserrat",
      font_size = 12
    ) +
      ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
        text        = ggplot2::element_text(family = "Montserrat")
      ) +
      ggplot2::ggtitle(feature_name)
    
    p_scatter <- add_background_static(
      p_scatter,
      bg = bg_df_static,
      size = 1,
      alpha = 0.35
    )
    
    print(p_scatter)
    dev.off()
    
    ## ---------------- SCATTER INTERACTIVE ----------------
    # Make sure background cannot “inherit” any categories
    if (!is.null(bg_df) && nrow(bg_df)) {
      bg_df <- bg_df %>%
        dplyr::select(-dplyr::any_of(c("category"))) %>% # safety
        dplyr::mutate(category = "all peptides") # force separate legend/trace
    }
    
    cat_cols <- c(
      "significant (wBH, per rank)" = "#FF1744",
      "nominal only"                = "#00E676",
      "not significant"             = "#2979FF",
      "all peptides"                = "#7A7A7A"
    )
    
    p_inter <- scatter_interactive(
      df = feature_data,
      xlab = group1,
      ylab = group2,
      peplib = peplib,
      show_background = TRUE,
      background_df = bg_df,
      background_name = "all peptides",
      background_color = "#808080",
      background_alpha = 0.40,
      background_size = 7,
      background_max_n = Inf,
      category_colors = cat_cols,
      point_size = 16,
      point_alpha = 0.95,
      jitter_width_pp = 0.05,
      jitter_height_pp = 0.05,
      font_size = 12
    )
    
    p_inter <- plotly::layout(
      p_inter,
      autosize = TRUE,
      margin   = list(l = 70, r = 30, t = 10, b = 70),
      xaxis    = list(range = c(-2, 102), automargin = TRUE),
      yaxis    = list(range = c(-2, 102), automargin = TRUE)
    )
    
    htmlwidgets::saveWidget(
      p_inter,
      file = paste0(file_prefix, "_scatter_interactive.html"),
      selfcontained = FALSE
    )
    
    ## ---------------- DELTA PLOT: conditional smooth ----------------
    use_smooth <- nrow(feature_data) >= 7
    smooth_k <- if (use_smooth) 3L else 1L
    
    ## ---------------- DELTA PLOT STATIC ----------------
    # open SVG device
    svglite(
      filename = paste0(file_prefix, "_deltaplot_static.svg"),
      width = 30/2.54,   # convert cm to inches
      height = 30/2.54,
      bg = "white"
    )
    
    p_delta_static <- tryCatch(
      {
        deltaplot(
          prev_tbl              = feature_data,
          group_pair_values     = c(group1, group2),
          group_labels          = c(group1, group2),
          point_jitter_width    = 0.01,
          point_jitter_height   = 0.01,
          point_alpha           = 0.6,
          point_size            = 6,
          add_smooth            = use_smooth,
          smooth_k              = smooth_k,
          arrow_head_length_mm  = 4
        ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat", size = 12))
      },
      error = function(e) {
        message(
          "Delta static plot failed for ", feature_name,
          " (", group1, " vs ", group2, "): ", conditionMessage(e)
        )
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(
            paste0("No delta plot for ", feature_name, "\n(", group1, " vs ", group2, ")")
          ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat"))
      }
    )
    print(p_delta_static)
    dev.off()
    
    ## ---------------- DELTA PLOT INTERACTIVE ----------------
    # p_delta <- tryCatch(
    #   deltaplot_interactive(
    #     prev_tbl            = feature_data,
    #     group_pair_values   = c(group1, group2),
    #     group_labels        = c(group1, group2),
    #     point_alpha         = 0.6,
    #     point_size          = 6,
    #     add_smooth          = use_smooth,
    #     smooth_k            = smooth_k,
    #     arrow_length_frac   = 0.35, # old arrow_frac_h
    #     point_jitter_width  = 0.01,
    #     point_jitter_height = 0.01
    #   ),
    #   error = function(e) {
    #     message(
    #       "Delta interactive plot failed for ", feature_name,
    #       " (", group1, " vs ", group2, "): ", conditionMessage(e)
    #     )
    #     NULL
    #   }
    # )
    # if (!is.null(p_delta)) {
    #   htmlwidgets::saveWidget(
    #     p_delta,
    #     file = paste0(file_prefix, "_deltaplot_interactive.html"),
    #     selfcontained = FALSE
    #   )
    # }
    
    ## ---------------- ECDF STATIC ----------------
    svglite(
      filename = paste0(file_prefix, "_ecdfplot_static.svg"),
      width = 30/2.54,
      height = 30/2.54,
      bg = "white"
    )
    
    p_ecdf_static <- tryCatch(
      {
        ecdf_plot(
          prev_tbl            = feature_data,
          group_pair_values   = c(group1, group2),
          group_labels        = c(group1, group2),
          line_width_pt       = 1,
          line_alpha          = 1,
          show_median_lines   = TRUE,
          show_ks_test        = TRUE
        ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat", size = 12))
      },
      error = function(e) {
        message(
          "ECDF static plot failed for ", feature_name,
          " (", group1, " vs ", group2, "): ", conditionMessage(e)
        )
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(
            paste0("No ECDF plot for ", feature_name, "\n(", group1, " vs ", group2, ")")
          ) +
          ggplot2::theme(text = ggplot2::element_text(family = "Montserrat"))
      }
    )
    print(p_ecdf_static)
    dev.off()
    
    ## ---------------- ECDF INTERACTIVE ----------------
    # p_ecdf <- tryCatch(
    #   ecdf_plot_interactive(
    #     prev_tbl            = feature_data,
    #     group_pair_values   = c(group1, group2),
    #     group_labels        = c(group1, group2),
    #     line_width_px       = 2,
    #     line_alpha          = 1,
    #     show_median_lines   = TRUE,
    #     show_ks_test        = TRUE
    #   ),
    #   error = function(e) {
    #     message(
    #       "ECDF interactive plot failed for ", feature_name,
    #       " (", group1, " vs ", group2, "): ", conditionMessage(e)
    #     )
    #     NULL
    #   }
    # )
    # if (!is.null(p_ecdf)) {
    #   htmlwidgets::saveWidget(
    #     p_ecdf,
    #     file = paste0(file_prefix, "_ecdfplot_interactive.html"),
    #     selfcontained = FALSE
    #   )
    # }
    
    invisible(NULL)
  }
  
  n_features <- nrow(res_with_pep)
  for (i in seq_len(n_features)) {
    feature_name <- res_with_pep$feature[i]
    group1 <- res_with_pep$pep_tbl_subset[[i]]$group1[1]
    group2 <- res_with_pep$pep_tbl_subset[[i]]$group2[1]
    feature_data <- res_with_pep$pep_tbl_subset[[i]]
    message("Plotting [", i, "/", n_features, "]: ", feature_name)
    plot_feature_all(feature_name, group1, group2, feature_data, out_dir)
  }
}

# restore original future plan
future::plan(original_plan)
# 
# 
# 
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# # 1) Keep only the two groups
# df_sub <- ps$data_long %>%
#   filter(group %in% c("Controls", "Adenoma"))
# 
# # 2) Total number of unique samples in each group
# n_samples_group <- df_sub %>%
#   distinct(sample_id, group) %>%
#   count(group, name = "n_samples_total")
# 
# # 3) Keep only present peptides
# df_present <- df_sub %>%
#   filter(exist == 1) %>%
#   distinct(sample_id, group, peptide_id)
# 
# # 4) Number of unique samples in which each peptide is present, by group
# prev_long <- df_present %>%
#   count(peptide_id, group, name = "n_present") %>%
#   left_join(n_samples_group, by = "group") %>%
#   mutate(prevalence = n_present / n_samples_total)
# 
# # 5) Make one row per peptide
# prev_df <- prev_long %>%
#   select(peptide_id, group, prevalence) %>%
#   pivot_wider(
#     names_from = group,
#     values_from = prevalence,
#     values_fill = 0
#   )
# 
# # Check: should be one row per peptide
# stopifnot(!anyDuplicated(prev_df$peptide_id))
# 
# # 6) Merge with peplib
# plot_df <- prev_df %>%
#   left_join(peplib, by = "peptide_id")
# 
# # 7) Flag Subdoligranulum
# # Replace 'genus' with the correct column in peplib if needed
# plot_df <- plot_df %>%
#   mutate(
#     is_subdoligranulum = ifelse(
#       grepl("Subdoligranulum", genus, ignore.case = TRUE),
#       "Subdoligranulum",
#       "Other"
#     )
#   )
# 
# # Optional check
# stopifnot(!anyDuplicated(plot_df$peptide_id))
# 
# # 8) Scatterplot
# ggplot(plot_df, aes(x = Controls, y = Adenoma)) +
#   geom_point(
#     data = subset(plot_df, is_subdoligranulum == "Other"),
#     color = "grey70", alpha = 0.5, size = 2
#   ) +
#   geom_point(
#     data = subset(plot_df, is_subdoligranulum == "Subdoligranulum"),
#     color = "red", alpha = 0.9, size = 2.5
#   ) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   labs(
#     x = "Prevalence in Controls",
#     y = "Prevalence in Adenoma",
#     title = "Peptide prevalence: Adenoma vs Controls"
#   ) +
#   theme_classic()