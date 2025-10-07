#' Pivot laboratory sample metadata to wide format
#'
#' @param si_long A data frame or lazy tibble containing `MetagenNumber`,
#'   `variable`, and `value` columns.
#'
#' @return A data frame with one row per `MetagenNumber`.
#' @export
create_sampleInfo_table <- function(si_long) {
  data <- ensure_long_format(
    si_long,
    id_col = "MetagenNumber",
    key_col = "variable",
    value_col = "value",
    context = "sample information"
  )

  tidyr::pivot_wider(data, names_from = "variable", values_from = "value") %>%
    dplyr::arrange(rlang::.data$MetagenNumber) %>%
    tibble::as_tibble()
}

#' Retrieve and merge sample and CMS metadata
#'
#' @param flist Optional logical expression evaluated with [subset()] to filter
#'   returned samples.
#' @param con An existing database connection. If `NULL`, a cached connection is
#'   reused.
#' @param cms_format One of "long" (default) to collect the `cmsdatalong`
#'   records or "wide" to use `cmsdata` directly.
#'
#' @return A data frame keyed by `MetagenNumber` containing laboratory and CMS
#'   metadata.
#' @export
fetch_sampleInfo <- function(flist = NULL, con = NULL, cms_format = c("long", "wide")) {
  con_obj <- resolve_connection(con)
  cms_format <- match.arg(cms_format)

  lab_tbl <- dplyr::tbl(con_obj, "labdata") %>%
    dplyr::select("MetagenNumber", "variable", "value") %>%
    dplyr::collect()
  sample_info <- create_sampleInfo_table(lab_tbl)

  cms_info <- switch(
    cms_format,
    long = {
      cms_long <- dplyr::tbl(con_obj, "cmsdatalong") %>%
        dplyr::select("MetagenNumber", "Factor", "Level") %>%
        dplyr::collect()
      create_cms_table(cms_long)
    },
    wide = {
      ensure_dataframe(
        dplyr::tbl(con_obj, "cmsdata") %>% dplyr::collect(),
        required = c(
          "MetagenNumber", "PropertyName", "SurveyID", "BlockName",
          "CropName", "SurveyDate", "Location", "AgronomistName",
          "TargetCropTypeName"
        ),
        context = "CMS metadata"
      )
    }
  )

  sample_info <- dplyr::full_join(sample_info, cms_info, by = "MetagenNumber")

  if (!is.null(flist)) {
    filter_expr <- substitute(flist)
    sample_info <- subset(sample_info, eval(filter_expr, sample_info, parent.frame()))
  }

  rownames(sample_info) <- sample_info$MetagenNumber
  sample_info
}

#' Retrieve sparse ASV abundance matrix
#'
#' @param con Database connection or pool. Defaults to the cached connection.
#' @param database Name of the ASV table to query.
#' @param phylo Unused but retained for backwards compatibility.
#' @param whichSamples Optional character vector of metagen numbers to restrict
#'   the result.
#'
#' @return A sparse matrix (`dgCMatrix`) of ASV abundances.
#' @export
fetch_asv_table_sparse <- function(con = NULL, database = "eukaryota_sv", phylo = FALSE, whichSamples = NULL) { # nolint
  build_asv_sparse(con, database, whichSamples)
}

#' Retrieve ASV abundance matrix
#'
#' @inheritParams fetch_asv_table_sparse
#'
#' @return A sparse `dgCMatrix` abundance matrix.
#' @export
fetch_asv_table <- function(con = NULL, database = "eukaryota_sv", phylo = FALSE, whichSamples = NULL) { # nolint
  build_asv_sparse(con, database, whichSamples)
}

#' Retrieve taxonomy records
#'
#' @param con Database connection or pool. Defaults to cached connection.
#' @param database Name of the taxonomy table to query.
#' @param whichTaxa Optional character vector restricting the ASVs returned.
#' @param phylo Logical; if `TRUE`, returns a matrix formatted for phyloseq.
#'
#' @return Either a tibble with taxonomy columns or a matrix if `phylo = TRUE`.
#' @export
fetch_taxonomy <- function(con = NULL, database = "eukaryota_tax", whichTaxa = NULL, phylo = FALSE) {
  con_obj <- resolve_connection(con)

  tax_tbl <- dplyr::tbl(con_obj, database) %>%
    dplyr::select(-dplyr::any_of("type"))

  if (!is.null(whichTaxa)) {
    whichTaxa <- unique(as.character(whichTaxa))
    tax_tbl <- dplyr::filter(tax_tbl, rlang::.data$SV %in% !!whichTaxa)
  }

  tax <- dplyr::collect(tax_tbl)

  tax[] <- lapply(tax, function(col) {
    if (is.character(col)) {
      gsub("\r", "", col)
    } else {
      col
    }
  })

  if (!isTRUE(phylo)) {
    attr(tax, "type") <- "taxonomy"
    return(tibble::as_tibble(tax))
  }

  sv_ids <- tax$SV
  tax_matrix <- tax %>%
    dplyr::select(-"SV") %>%
    as.matrix()
  rownames(tax_matrix) <- sv_ids
  tax_matrix
}

#' Pivot CMS metadata to wide format
#'
#' @param cms_long A data frame or lazy tibble containing `MetagenNumber`,
#'   `Factor`, and `Level` columns.
#'
#' @return A tibble with one row per `MetagenNumber` and CMS factors as columns.
#' @export
create_cms_table <- function(cms_long) {
  data <- ensure_long_format(
    cms_long,
    id_col = "MetagenNumber",
    key_col = "Factor",
    value_col = "Level",
    context = "CMS metadata"
  )

  tidyr::pivot_wider(data, names_from = "Factor", values_from = "Level") %>%
    dplyr::arrange(rlang::.data$MetagenNumber) %>%
    tibble::as_tibble()
}

#' Retrieve ASV records for specific samples
#'
#' @param samples Character vector of metagen numbers.
#' @param database Table containing ASV abundances.
#' @param con Database connection.
#' @param col Column name identifying the sample identifier field.
#'
#' @return A tibble with ASV abundance rows.
#' @export
fetch_asvs_by_sample <- function(samples, database = NULL, con = NULL, col = "MetagenNumber") {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when fetching ASVs by sample.", call. = FALSE)
  }

  con_obj <- resolve_connection(con)
  samples <- unique(as.character(samples))

  column <- rlang::sym(col)

  dplyr::tbl(con_obj, database) %>%
    dplyr::filter(!!column %in% !!samples) %>%
    dplyr::collect()
}

#' Retrieve sparse ASV matrix for specified samples
#'
#' @inheritParams fetch_asv_table_sparse
#'
#' @return A sparse `dgCMatrix` abundance matrix limited to the requested samples.
#' @export
fetch_asv_table_sparse_by_sample <- function(con = NULL, database = "eukaryota_sv", phylo = FALSE, whichSamples = NULL) { # nolint
  if (is.null(whichSamples)) {
    stop("`whichSamples` must be supplied when requesting a sample-specific ASV table.", call. = FALSE)
  }

  build_asv_sparse(con, database, unique(as.character(whichSamples)))
}

#' Retrieve taxonomy for specific ASVs
#'
#' @param taxa Character vector of ASV identifiers.
#' @inheritParams fetch_taxonomy
#'
#' @return A tibble of taxonomy records.
#' @export
fetch_taxonomy_by_asv <- function(taxa, database = NULL, con = NULL, col = "SV") {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when fetching taxonomy by ASV.", call. = FALSE)
  }

  con_obj <- resolve_connection(con)
  taxa <- unique(as.character(taxa))

  column <- rlang::sym(col)

  dplyr::tbl(con_obj, database) %>%
    dplyr::filter(!!column %in% !!taxa) %>%
    dplyr::collect() %>%
    dplyr::mutate(dplyr::across(where(is.character), ~ gsub("\r", "", .x))) %>%
    tibble::as_tibble()
}

#' Retrieve distinct ASV identifiers from a table
#'
#' @inheritParams fetch_asv_table_sparse
#'
#' @return A character vector of distinct ASV identifiers.
#' @export
get_svs <- function(database = NULL, con = NULL) {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when retrieving ASV identifiers.", call. = FALSE)
  }

  con_obj <- resolve_connection(con)

  dplyr::tbl(con_obj, database) %>%
    dplyr::distinct(rlang::.data$SV) %>%
    dplyr::arrange(rlang::.data$SV) %>%
    dplyr::pull("SV")
}

#' Add single quotes to character values
#'
#' @param x Character vector of values to quote.
#'
#' @return A character vector wrapped in single quotes.
#' @export
add_quotes <- function(x) {
  sprintf("'%s'", unique(x))
}

#' Construct a sparse ASV matrix from long-format data
#'
#' @param asv_long Data frame containing `MetagenNumber`, `SV`, and `Abundance` columns.
#'
#' @return A sparse `dgCMatrix` matrix.
#' @export
construct_asv_table <- function(asv_long) {
  data <- ensure_dataframe(
    asv_long,
    required = c("MetagenNumber", "SV", "Abundance"),
    context = "ASV abundance"
  )

  if (!nrow(data)) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(0, 0)))
  }

  samples <- unique(data$MetagenNumber)
  taxa <- unique(data$SV)

  i <- match(data$MetagenNumber, samples)
  j <- match(data$SV, taxa)

  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = data$Abundance,
    dims = c(length(samples), length(taxa)),
    dimnames = list(samples, taxa)
  )
}

build_asv_sparse <- function(con, database, samples = NULL) {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when fetching ASV data.", call. = FALSE)
  }

  con_obj <- resolve_connection(con)

  asv_tbl <- dplyr::tbl(con_obj, database) %>%
    dplyr::select("MetagenNumber", "SV", "Abundance")

  if (!is.null(samples)) {
    samples <- unique(as.character(samples))
    asv_tbl <- dplyr::filter(asv_tbl, rlang::.data$MetagenNumber %in% !!samples)
  }

  asv_long <- dplyr::collect(asv_tbl)
  construct_asv_table(asv_long)
}

ensure_dataframe <- function(data, required, context) {
  if (inherits(data, "tbl_lazy")) {
    data <- dplyr::collect(data)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      sprintf(
        "Missing required column(s) for %s: %s",
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  data
}

ensure_long_format <- function(data, id_col, key_col, value_col, context) {
  data <- ensure_dataframe(
    data,
    required = c(id_col, key_col, value_col),
    context = context
  )

  data[[id_col]] <- as.character(data[[id_col]])
  data[[key_col]] <- as.character(data[[key_col]])
  data
}
