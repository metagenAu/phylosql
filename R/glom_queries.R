#' Build a phyloseq object aggregated at a taxonomy rank
#'
#' @param taxalevel Character name of the taxonomy column that defines the
#'   aggregation level (for example, `"Family"`).
#' @param taxa_group Character prefix used to identify the source tables. The
#'   function expects `<taxa_group>_sv` and `<taxa_group>_tax` tables to exist in
#'   the connected schema.
#' @param con Optional database connection, pool, or connection expression. When
#'   omitted the function will reuse the most recently established connection.
#'
#' @return A `phyloseq` object built from the aggregated counts and taxonomy
#'   metadata.
#' @export
sql_phyloseq_by_tax_glom <- function(taxalevel, taxa_group, con = NULL) {
  con_obj <- resolve_connection(con)
  validate_taxa_group(taxa_group)

  tax_fields <- DBI::dbListFields(con_obj, paste0(taxa_group, "_tax"))
  if (is.null(taxalevel) || !nzchar(taxalevel)) {
    stop("`taxalevel` must be a non-empty taxonomy column name.", call. = FALSE)
  }
  if (!taxalevel %in% tax_fields) {
    stop(
      sprintf(
        "`taxalevel` %s is not present in %s_tax.",
        shQuote(taxalevel),
        taxa_group
      ),
      call. = FALSE
    )
  }

  idx <- match(taxalevel, tax_fields)
  tax_cols <- unique(c("SV", tax_fields[seq_len(idx)]))
  taxonomy_group_cols <- setdiff(tax_cols, "SV")
  if (length(taxonomy_group_cols) == 0) {
    stop("`taxalevel` must reference a taxonomy column other than `SV`.", call. = FALSE)
  }

  sv_tax_data <- collect_sv_tax(con_obj, taxa_group)
  if (!nrow(sv_tax_data)) {
    stop("No sequence variants were found for the requested inputs.", call. = FALSE)
  }

  aggregated <- dplyr::group_by(
    sv_tax_data,
    dplyr::across(dplyr::all_of(c(taxonomy_group_cols, "MetagenNumber")))
  )
  aggregated <- dplyr::summarise(aggregated, Abundance = sum(Abundance), .groups = "drop")
  aggregated <- dplyr::filter(aggregated, Abundance > 0)

  if (!nrow(aggregated)) {
    stop("All abundances were zero after glomming at the requested taxonomy level.", call. = FALSE)
  }

  taxonomy_lookup <- collect_taxonomy_lookup(con_obj, taxa_group, tax_cols)
  aggregated <- dplyr::left_join(
    aggregated,
    taxonomy_lookup,
    by = taxonomy_group_cols
  )
  aggregated <- dplyr::filter(aggregated, !is.na(SV))

  if (!nrow(aggregated)) {
    stop("Unable to map glommed abundances to representative SV identifiers.", call. = FALSE)
  }

  taxonomy <- dplyr::distinct(aggregated, SV, .keep_all = TRUE)
  build_phyloseq_object(
    con_obj = con_obj,
    counts = dplyr::select(aggregated, MetagenNumber, SV, Abundance),
    taxonomy = dplyr::select(taxonomy, dplyr::all_of(tax_cols))
  )
}

#' Build a phyloseq object for selected samples
#'
#' @param taxa_group Character prefix used to identify the source tables. The
#'   function expects `<taxa_group>_sv` and `<taxa_group>_tax` tables to exist in
#'   the connected schema.
#' @param samples Optional character vector of sample identifiers to keep.
#' @param con Optional database connection, pool, or connection expression.
#'
#' @return A `phyloseq` object restricted to the requested samples.
#' @export
sql_phyloseq_by_sample <- function(taxa_group, samples, con = NULL) {
  con_obj <- resolve_connection(con)
  validate_taxa_group(taxa_group)

  if (missing(samples)) {
    samples <- NULL
  }

  sv_tax_data <- collect_sv_tax(con_obj, taxa_group, samples)
  sv_tax_data <- dplyr::filter(sv_tax_data, Abundance > 0)

  if (!nrow(sv_tax_data)) {
    stop("No abundances were found for the requested samples.", call. = FALSE)
  }

  tax_fields <- DBI::dbListFields(con_obj, paste0(taxa_group, "_tax"))
  taxonomy <- dplyr::distinct(sv_tax_data, SV, .keep_all = TRUE)

  build_phyloseq_object(
    con_obj = con_obj,
    counts = dplyr::select(sv_tax_data, MetagenNumber, SV, Abundance),
    taxonomy = dplyr::select(taxonomy, dplyr::all_of(tax_fields))
  )
}

#' Build a phyloseq object for selected samples aggregated at a taxonomy rank
#'
#' @inheritParams sql_phyloseq_by_tax_glom
#' @inheritParams sql_phyloseq_by_sample
#'
#' @return A `phyloseq` object restricted to the supplied samples and aggregated
#'   at the requested taxonomy level.
#' @export
sql_phyloseq_by_sample_and_tax_glom <- function(taxa_group, taxalevel, samples, con = NULL) {
  con_obj <- resolve_connection(con)
  validate_taxa_group(taxa_group)

  if (missing(samples) || is.null(samples) || !length(samples)) {
    stop("`samples` must contain at least one sample identifier.", call. = FALSE)
  }

  tax_fields <- DBI::dbListFields(con_obj, paste0(taxa_group, "_tax"))
  if (!taxalevel %in% tax_fields) {
    stop(
      sprintf(
        "`taxalevel` %s is not present in %s_tax.",
        shQuote(taxalevel),
        taxa_group
      ),
      call. = FALSE
    )
  }

  idx <- match(taxalevel, tax_fields)
  tax_cols <- unique(c("SV", tax_fields[seq_len(idx)]))
  taxonomy_group_cols <- setdiff(tax_cols, "SV")
  if (length(taxonomy_group_cols) == 0) {
    stop("`taxalevel` must reference a taxonomy column other than `SV`.", call. = FALSE)
  }

  sv_tax_data <- collect_sv_tax(con_obj, taxa_group, samples)
  if (!nrow(sv_tax_data)) {
    stop("No sequence variants were found for the requested inputs.", call. = FALSE)
  }

  aggregated <- dplyr::group_by(
    sv_tax_data,
    dplyr::across(dplyr::all_of(c(taxonomy_group_cols, "MetagenNumber")))
  )
  aggregated <- dplyr::summarise(aggregated, Abundance = sum(Abundance), .groups = "drop")
  aggregated <- dplyr::filter(aggregated, Abundance > 0)

  if (!nrow(aggregated)) {
    stop("All abundances were zero after glomming at the requested taxonomy level.", call. = FALSE)
  }

  taxonomy_lookup <- collect_taxonomy_lookup(con_obj, taxa_group, tax_cols)
  aggregated <- dplyr::left_join(
    aggregated,
    taxonomy_lookup,
    by = taxonomy_group_cols
  )
  aggregated <- dplyr::filter(aggregated, !is.na(SV))

  if (!nrow(aggregated)) {
    stop("Unable to map glommed abundances to representative SV identifiers.", call. = FALSE)
  }

  taxonomy <- dplyr::distinct(aggregated, SV, .keep_all = TRUE)
  build_phyloseq_object(
    con_obj = con_obj,
    counts = dplyr::select(aggregated, MetagenNumber, SV, Abundance),
    taxonomy = dplyr::select(taxonomy, dplyr::all_of(tax_cols))
  )
}

resolve_connection <- function(con) {
  if (is.null(con)) {
    con <- try_fetch_connection()
    if (is.null(con)) {
      stop("No connection available. Call `get_mtgn_connection()` or supply `con`.", call. = FALSE)
    }
    return(con)
  }

  if (is.character(con)) {
    con <- eval_con(con)
  }

  if (inherits(con, "Pool")) {
    if (pool::poolClosed(con)) {
      stop("The supplied connection pool has been closed.", call. = FALSE)
    }
    return(con)
  }

  if (inherits(con, "DBIConnection")) {
    if (!DBI::dbIsValid(con)) {
      stop("The supplied connection is no longer valid.", call. = FALSE)
    }
    return(con)
  }

  stop("`con` must be a DBI connection, pool, or connection expression string.", call. = FALSE)
}

collect_sv_tax <- function(con, taxa_group, samples = NULL) {
  sv_table <- paste0(taxa_group, "_sv")
  tax_table <- paste0(taxa_group, "_tax")

  sv_tbl <- dplyr::tbl(con, sv_table)
  tax_tbl <- dplyr::tbl(con, tax_table)

  joined <- dplyr::inner_join(sv_tbl, tax_tbl, by = "SV")

  if (!is.null(samples)) {
    samples <- unique(samples)
    if (!length(samples)) {
      stop("`samples` must contain at least one sample identifier.", call. = FALSE)
    }
    joined <- dplyr::filter(joined, MetagenNumber %in% samples)
  }

  dplyr::collect(joined)
}

collect_taxonomy_lookup <- function(con, taxa_group, tax_cols) {
  tax_table <- paste0(taxa_group, "_tax")

  lookup <- dplyr::tbl(con, tax_table)
  lookup <- dplyr::select(lookup, dplyr::all_of(tax_cols))
  lookup <- dplyr::collect(lookup)

  lookup <- lookup[order(lookup$SV), , drop = FALSE]
  if (length(tax_cols) > 1) {
    dup_cols <- setdiff(tax_cols, "SV")
    lookup <- lookup[!duplicated(lookup[dup_cols]), , drop = FALSE]
  }

  lookup
}

build_phyloseq_object <- function(con_obj, counts, taxonomy) {
  counts <- dplyr::mutate(counts,
    MetagenNumber = as.character(MetagenNumber),
    SV = as.character(SV)
  )
  counts <- dplyr::arrange(counts, MetagenNumber, SV)

  if (!nrow(counts)) {
    stop("No abundance records are available to build a phyloseq object.", call. = FALSE)
  }

  taxonomy <- dplyr::mutate(taxonomy, SV = as.character(SV))
  taxonomy <- taxonomy[match(unique(counts$SV), taxonomy$SV), , drop = FALSE]
  taxonomy <- taxonomy[!is.na(taxonomy$SV), , drop = FALSE]
  taxonomy_unique <- taxonomy
  if (nrow(taxonomy_unique) == 0) {
    stop("No taxonomy records were available for the supplied SV identifiers.", call. = FALSE)
  }

  tax_matrix <- as.matrix(taxonomy_unique[, setdiff(colnames(taxonomy_unique), "SV"), drop = FALSE])
  rownames(tax_matrix) <- taxonomy_unique$SV

  asv_table <- construct_asv_table(
    dplyr::select(counts, SV, MetagenNumber, Abundance)
  )

  sample_info <- fetch_sampleInfo(con = con_obj)
  keep_samples <- unique(counts$MetagenNumber)
  missing_samples <- setdiff(keep_samples, rownames(sample_info))
  if (length(missing_samples) > 0) {
    stop(
      sprintf(
        "Sample metadata is missing for: %s",
        paste(missing_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sample_info <- sample_info[keep_samples, , drop = FALSE]

  phyloseqSparse::phyloseq(
    phyloseqSparse::otu_table(asv_table, taxa_are_rows = FALSE),
    phyloseqSparse::tax_table(tax_matrix),
    phyloseqSparse::sample_data(sample_info)
  )
}

validate_taxa_group <- function(taxa_group) {
  if (missing(taxa_group) || is.null(taxa_group) || length(taxa_group) != 1 || !is.character(taxa_group) || !nzchar(taxa_group)) {
    stop("`taxa_group` must be a non-empty character string.", call. = FALSE)
  }
  invisible(NULL)
}
