suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
  library(meta)
})

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------
table_s11_file <- "Table S11 Data from Studies Included in the Meta-analysis.xlsx"
table_s12_file <- "Table S12 List of Studies and Reasons for Exclusion.xlsx"
asian_file <- "Asian Melanoma.xlsx"
results_dir <- "results"

key_genes <- c("BRAF", "NRAS", "TP53", "TERT", "NF1")
table2_method_order <- c(
  "NGS",
  "PCR",
  "PCR and Sanger sequencing",
  "Pyrosequencing",
  "Target NGS",
  "Sanger",
  "Sequenom Massarray",
  "WES",
  "WGS"
)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
clean_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "[\r\n]+", " ")
  x <- str_squish(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

sanitize_name_one <- function(x) {
  x <- clean_text(x)
  if (is.na(x)) {
    return("blank")
  }

  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  if (identical(x, "")) {
    x <- "blank"
  }
  x
}

sanitize_names_unique <- function(x) {
  make.unique(vapply(x, sanitize_name_one, character(1)), sep = "_")
}

extract_first_numeric <- function(x) {
  x <- clean_text(x)
  x <- ifelse(is.na(x), NA_character_, str_replace_all(x, ",", "."))
  suppressWarnings(as.numeric(str_extract(x, "-?[0-9]+\\.?[0-9]*")))
}

collapse_methodology <- function(x) {
  y <- str_to_lower(clean_text(x))

  case_when(
    is.na(y) ~ NA_character_,
    y %in% c("target ngs", "targeted ngs") ~ "Target NGS",
    str_detect(y, "sequenom|massarray|mela.?carta") ~ "Sequenom Massarray",
    y == "pcr and sanger sequencing" ~ "PCR and Sanger sequencing",
    str_detect(y, "pyro") ~ "Pyrosequencing",
    y %in% c(
      "sanger sequencing",
      "conventional dna sequencing",
      "direct sequencing",
      "direct sequencing of pcr products"
    ) ~ "Sanger",
    y == "ngs" ~ "NGS",
    str_detect(y, "\\bwes\\b") ~ "WES",
    str_detect(y, "\\bwgs\\b") ~ "WGS",
    str_detect(y, "pcr|rtpcr|rt pcr|qrt pcr|prt pcr|castpcr|cast pcr|ms pcr|snapshot|droplet digital pcr") ~ "PCR",
    TRUE ~ clean_text(x)
  )
}

read_gene_mutation_sheet <- function(path) {
  raw <- read_excel(path, sheet = "Gene mutation", col_names = FALSE)

  if (nrow(raw) < 3) {
    stop("The 'Gene mutation' sheet must contain two header rows plus data.")
  }

  if (ncol(raw) < 162) {
    stop("Unexpected 'Gene mutation' sheet format.")
  }

  header_row_2 <- as.character(unlist(raw[2, ]))
  data <- raw[-c(1, 2), ]

  metadata_headers <- c(
    "TT",
    "Author",
    "Method",
    "Type_tumor",
    "Multiple_sample",
    "Acral_melanoma",
    "Sample_N",
    "Patient_N",
    "Age",
    "Male",
    "Female",
    "Sample_primary",
    "Patient_primary",
    "Tumor_burden_score",
    "Fraction_genome_altered",
    "Ethnic",
    "Methodology"
  )

  primary_gene_headers <- sanitize_names_unique(header_row_2[18:155])

  metastatic_metadata_headers <- c(
    "Method_metastatic",
    "Type_tumor_metastatic",
    "Multiple_sample_metastatic",
    "Metastatic_site",
    "Metastatic_samples_N",
    "Patient_metastatic_N"
  )

  metastatic_gene_headers <- paste0(
    sanitize_names_unique(header_row_2[162:ncol(raw)]),
    "_metastatic"
  )

  names(data) <- c(
    metadata_headers,
    primary_gene_headers,
    metastatic_metadata_headers,
    metastatic_gene_headers
  )

  data %>%
    mutate(
      TT = as.integer(extract_first_numeric(TT)),
      Author = clean_text(Author),
      Method = clean_text(Method),
      Methodology = clean_text(Methodology),
      Method_group = collapse_methodology(Methodology),
      across(
        c(
          Sample_N,
          Patient_N,
          Sample_primary,
          Patient_primary,
          Metastatic_samples_N,
          Patient_metastatic_N
        ),
        extract_first_numeric
      )
    )
}

prepare_gene_data <- function(data, gene, tumor = c("primary", "metastatic")) {
  tumor <- match.arg(tumor)

  gene_column <- if (tumor == "primary") {
    gene
  } else {
    paste0(gene, "_metastatic")
  }

  n_column <- if (tumor == "primary") "Sample_primary" else "Metastatic_samples_N"

  if (!gene_column %in% names(data)) {
    return(tibble())
  }

  data %>%
    transmute(
      TT = TT,
      Author = Author,
      Methodology = Methodology,
      Method_group = Method_group,
      n = extract_first_numeric(.data[[n_column]]),
      event = extract_first_numeric(.data[[gene_column]])
    ) %>%
    filter(!is.na(TT), !is.na(event), !is.na(n), n > 0) %>%
    mutate(studylab = TT)
}

fit_prevalence_model <- function(d) {
  if (is.null(d) || nrow(d) == 0) {
    return(NULL)
  }

  metaprop(
    event = event,
    n = n,
    studlab = studylab,
    data = d,
    sm = "PFT",
    method = "Inverse",
    method.tau = "DL",
    common = FALSE,
    random = TRUE,
    backtransf = TRUE
  )
}

transf_ipft_hm <- function(theta, n_vec) {
  n_vec <- as.numeric(n_vec)
  nhm <- length(n_vec) / sum(1 / n_vec)
  s <- sin(theta)
  inside <- 1 - (s + (s - 1 / s) / nhm)^2
  inside[inside < 0] <- 0
  0.5 * (1 - sign(cos(theta)) * sqrt(inside))
}

format_p_het <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x < 0.01, "<0.01", sprintf("%.2f", x)))
}

summarise_metaprop <- function(model, d, gene, tumor_label) {
  if (is.null(model) || nrow(d) == 0) {
    return(tibble())
  }

  tibble(
    Gene = gene,
    `Tumor type` = tumor_label,
    Studies = dplyr::n_distinct(d$TT),
    `Sample size` = sum(d$n, na.rm = TRUE),
    `Total mutated` = sum(d$event, na.rm = TRUE),
    `TT IDs` = paste(sort(unique(d$TT)), collapse = ", "),
    `Pooled prevalence (%)` = round(100 * transf_ipft_hm(model$TE.random, d$n), 1),
    `95% CI` = sprintf(
      "%.1f–%.1f",
      100 * transf_ipft_hm(model$lower.random, d$n),
      100 * transf_ipft_hm(model$upper.random, d$n)
    ),
    `I2 (%)` = ifelse(is.null(model$I2), NA_real_, round(as.numeric(model$I2), 1)),
    `Heterogeneity P` = format_p_het(model$pval.Q)
  )
}

build_method_summary <- function(data, genes) {
  purrr::map_dfr(genes, function(gene) {
    d_gene <- prepare_gene_data(data, gene = gene, tumor = "primary") %>%
      filter(!is.na(Method_group))

    if (nrow(d_gene) == 0) {
      return(tibble())
    }

    method_levels <- intersect(table2_method_order, unique(d_gene$Method_group))

    purrr::map_dfr(method_levels, function(method_name) {
      d_sub <- d_gene %>%
        filter(Method_group == method_name)

      model <- fit_prevalence_model(d_sub)

      if (is.null(model)) {
        return(tibble())
      }

      tibble(
        Gene = gene,
        Subgroup = method_name,
        Studies = dplyr::n_distinct(d_sub$TT),
        `Sample size` = sum(d_sub$n, na.rm = TRUE),
        `Total mutated` = sum(d_sub$event, na.rm = TRUE),
        `TT IDs` = paste(sort(unique(d_sub$TT)), collapse = ", "),
        `Pooled prevalence (%)` = round(100 * transf_ipft_hm(model$TE.random, d_sub$n), 1),
        `95% CI` = sprintf(
          "%.1f–%.1f",
          100 * transf_ipft_hm(model$lower.random, d_sub$n),
          100 * transf_ipft_hm(model$upper.random, d_sub$n)
        ),
        `I2 (%)` = ifelse(is.null(model$I2), NA_real_, round(as.numeric(model$I2), 1)),
        `Heterogeneity P` = format_p_het(model$pval.Q)
      )
    })
  })
}

build_asian_summary <- function(data, genes = key_genes) {
  purrr::map_dfr(genes, function(gene) {
    d_primary <- prepare_gene_data(data, gene = gene, tumor = "primary")
    d_metastatic <- prepare_gene_data(data, gene = gene, tumor = "metastatic")

    model_primary <- fit_prevalence_model(d_primary)
    model_metastatic <- fit_prevalence_model(d_metastatic)

    bind_rows(
      summarise_metaprop(model_primary, d_primary, gene, "Primary"),
      summarise_metaprop(model_metastatic, d_metastatic, gene, "Metastatic")
    )
  })
}

build_all_gene_summary <- function(data) {
  gene_columns <- names(data)[18:155]
  gene_columns <- gene_columns[gene_columns != "blank"]

  purrr::map_dfr(gene_columns, function(gene) {
    d_primary <- prepare_gene_data(data, gene = gene, tumor = "primary")
    d_metastatic <- prepare_gene_data(data, gene = gene, tumor = "metastatic")

    model_primary <- fit_prevalence_model(d_primary)
    model_metastatic <- fit_prevalence_model(d_metastatic)

    bind_rows(
      summarise_metaprop(model_primary, d_primary, gene, "Primary"),
      summarise_metaprop(model_metastatic, d_metastatic, gene, "Metastatic")
    )
  })
}

save_model_summary <- function(model, file_path) {
  if (is.null(model)) {
    return(invisible(NULL))
  }

  txt <- capture.output({
    print(summary(model))
  })

  writeLines(txt, con = file_path)
}

write_methodology_map <- function(data, file_path) {
  data %>%
    distinct(TT, Author, Methodology, Method_group) %>%
    arrange(TT) %>%
    write_csv(file_path, na = "")
}

# -------------------------------------------------------------------
# Main analysis
# -------------------------------------------------------------------
if (!file.exists(table_s11_file)) {
  stop(sprintf("Input file not found: %s", table_s11_file))
}

main_data <- read_gene_mutation_sheet(table_s11_file)

write_csv(main_data, file.path(results_dir, "tableS11_gene_mutation_clean.csv"), na = "")
write_methodology_map(main_data, file.path(results_dir, "methodology_mapping.csv"))

available_primary_genes <- names(main_data)[18:155]
analysis_genes <- key_genes[key_genes %in% available_primary_genes]

# Overall primary models for key genes
primary_models <- setNames(
  lapply(analysis_genes, function(gene) fit_prevalence_model(prepare_gene_data(main_data, gene, "primary"))),
  analysis_genes
)

saveRDS(primary_models, file.path(results_dir, "primary_overall_models.rds"))
iwalk(
  primary_models,
  function(model, gene) {
    save_model_summary(model, file.path(results_dir, paste0("primary_", gene, "_summary.txt")))
  }
)

# Overall metastatic models for key genes
metastatic_analysis_genes <- analysis_genes[analysis_genes %in% sub("_metastatic$", "", names(main_data)[162:ncol(main_data)])]
metastatic_models <- setNames(
  lapply(metastatic_analysis_genes, function(gene) fit_prevalence_model(prepare_gene_data(main_data, gene, "metastatic"))),
  metastatic_analysis_genes
)

saveRDS(metastatic_models, file.path(results_dir, "metastatic_overall_models.rds"))
iwalk(
  metastatic_models,
  function(model, gene) {
    save_model_summary(model, file.path(results_dir, paste0("metastatic_", gene, "_summary.txt")))
  }
)

# Table 2 method summary
table2_method_summary <- build_method_summary(main_data, analysis_genes)
write_csv(table2_method_summary, file.path(results_dir, "table2_method_summary.csv"), na = "")

# Optional Asian pooled prevalence analysis
if (file.exists(asian_file)) {
  asian_data <- read_gene_mutation_sheet(asian_file)

  write_csv(asian_data, file.path(results_dir, "asian_gene_mutation_clean.csv"), na = "")
  write_methodology_map(asian_data, file.path(results_dir, "asian_methodology_mapping.csv"))

  asian_key_summary <- build_asian_summary(asian_data, genes = analysis_genes)
  asian_all_summary <- build_all_gene_summary(asian_data)

  write_csv(asian_key_summary, file.path(results_dir, "asian_key_genes.csv"), na = "")
  write_csv(asian_all_summary, file.path(results_dir, "asian_all_genes.csv"), na = "")

  asian_primary_models <- setNames(
    lapply(analysis_genes, function(gene) fit_prevalence_model(prepare_gene_data(asian_data, gene, "primary"))),
    analysis_genes
  )

  asian_metastatic_models <- setNames(
    lapply(metastatic_analysis_genes, function(gene) fit_prevalence_model(prepare_gene_data(asian_data, gene, "metastatic"))),
    metastatic_analysis_genes
  )

  saveRDS(
    list(primary = asian_primary_models, metastatic = asian_metastatic_models),
    file.path(results_dir, "asian_models.rds")
  )

  iwalk(
    asian_primary_models,
    function(model, gene) {
      save_model_summary(model, file.path(results_dir, paste0("asian_primary_", gene, "_summary.txt")))
    }
  )

  iwalk(
    asian_metastatic_models,
    function(model, gene) {
      save_model_summary(model, file.path(results_dir, paste0("asian_metastatic_", gene, "_summary.txt")))
    }
  )
} else {
  message("Asian Melanoma.xlsx not found. Asia-specific analysis was skipped.")
}

capture.output(sessionInfo(), file = file.path(results_dir, "sessionInfo.txt"))

message("Analysis complete. Output files are available in: ", results_dir)
