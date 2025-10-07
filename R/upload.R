resolve_connection <- function(con = NULL) {
  if (!is.null(con)) {
    if (!DBI::dbIsValid(con)) {
      stop("`con` must be a valid DBI connection or pool object.", call. = FALSE)
    }
    return(con)
  }

  cached <- try_fetch_connection()
  if (is.null(cached)) {
    stop("No cached database connection available. Provide `con` or establish one with `get_mtgn_connection()`.", call. = FALSE)
  }

  if (!DBI::dbIsValid(cached)) {
    stop("Cached database connection is no longer valid. Reconnect with `get_mtgn_connection()`.", call. = FALSE)
  }

  cached
}

validate_upload_data <- function(data, required_cols, context) {
  if (is.null(data)) {
    stop(sprintf("`data` must not be NULL for %s uploads.", context), call. = FALSE)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s) for %s upload: %s",
        context,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  data
}

combine_key <- function(df, columns) {
  if (length(columns) == 1L) {
    return(as.character(df[[columns]]))
  }

  rows <- df[, columns, drop = FALSE]
  apply(rows, 1, function(row) {
    values <- ifelse(is.na(row), "<NA>", as.character(row))
    paste(values, collapse = "\r")
  })
}

collect_existing_keys <- function(con, table, key_cols) {
  tbl <- dplyr::tbl(con, table)
  selected <- dplyr::select(tbl, dplyr::all_of(key_cols))
  dplyr::collect(selected)
}

compute_new_indices <- function(existing, candidate, key_cols) {
  if (!nrow(candidate)) {
    return(integer())
  }

  existing_keys <- if (nrow(existing)) combine_key(existing, key_cols) else character()
  candidate_keys <- combine_key(candidate, key_cols)

  which(!candidate_keys %in% existing_keys)
}

perform_append_upload <- function(con, table, data, key_cols, context, required_cols = key_cols, preprocess = identity) {
  con <- resolve_connection(con)
  data <- validate_upload_data(data, required_cols, context)
  data <- preprocess(data)

  if (!nrow(data)) {
    message("Input data has 0 rows; nothing to upload.")
    return(invisible(data))
  }

  existing <- collect_existing_keys(con, table, unique(key_cols))
  indices <- compute_new_indices(existing, data, unique(key_cols))

  if (!length(indices)) {
    message("No new records to upload.")
    return(invisible(data[0, , drop = FALSE]))
  }

  rows_to_upload <- data[indices, , drop = FALSE]
  message(sprintf("Uploading %d record(s) to %s...", nrow(rows_to_upload), table))
  DBI::dbAppendTable(con, table, value = rows_to_upload)
  message("Upload complete.")

  invisible(rows_to_upload)
}

remove_carriage_returns <- function(df) {
  for (col in names(df)) {
    if (is.character(df[[col]])) {
      df[[col]] <- gsub("\r", "", df[[col]])
    }
  }
  df
}

#' Upload laboratory data to a MySQL database
#'
#' @param data A data frame containing laboratory records.
#' @param database The name of the target table. Defaults to `"labdata"`.
#' @param con An existing database connection. If `NULL`, a cached connection will be used.
#'
#' @return Invisibly returns the rows uploaded.
#' @export
#' @importFrom DBI dbAppendTable dbIsValid
upload_lab_data <- function(data, database = "labdata", con = NULL) {
  perform_append_upload(
    con = con,
    table = database,
    data = data,
    key_cols = c("MetagenNumber", "variable"),
    context = "lab data"
  )
}

#' Upload sequence variant data to a MySQL database
#'
#' @inheritParams upload_lab_data
#'
#' @return Invisibly returns the rows uploaded.
#' @export
upload_sv <- function(data, database = NULL, con = NULL) {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when uploading sequence variant data.", call. = FALSE)
  }

  perform_append_upload(
    con = con,
    table = database,
    data = data,
    key_cols = c("MetagenNumber", "SV"),
    context = "sequence variant data"
  )
}

#' Upload taxonomy data to a MySQL database
#'
#' @inheritParams upload_lab_data
#'
#' @return Invisibly returns the rows uploaded.
#' @export
upload_taxonomy <- function(data, database = NULL, con = NULL) {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when uploading taxonomy data.", call. = FALSE)
  }

  perform_append_upload(
    con = con,
    table = database,
    data = data,
    key_cols = "SV",
    context = "taxonomy data",
    preprocess = remove_carriage_returns
  )
}

#' Upload CMS data to a MySQL database
#'
#' @inheritParams upload_lab_data
#'
#' @return Invisibly returns the rows uploaded.
#' @export
upload_cms_data <- function(data, database = "cmsdata", con = NULL) {
  perform_append_upload(
    con = con,
    table = database,
    data = data,
    key_cols = "MetagenNumber",
    context = "CMS data"
  )
}

#' Upload CMS long-format data to a MySQL database
#'
#' @inheritParams upload_lab_data
#'
#' @return Invisibly returns the rows uploaded.
#' @export
upload_cms_data_Long <- function(data, database = "cmsdatalong", con = NULL) {
  con <- resolve_connection(con)
  data <- validate_upload_data(data, c("MetagenNumber", "Factor", "Level"), "CMS long data")

  if (!nrow(data)) {
    message("Input data has 0 rows; nothing to upload.")
    return(invisible(data))
  }

  if (any(is.na(data$Level))) {
    stop("`Level` column contains missing values. Remove them before uploading.", call. = FALSE)
  }

  existing <- collect_existing_keys(con, database, c("MetagenNumber", "Factor"))
  indices <- compute_new_indices(existing, data, c("MetagenNumber", "Factor"))

  if (!length(indices)) {
    message("No new records to upload.")
    return(invisible(data[0, , drop = FALSE]))
  }

  rows_to_upload <- data[indices, , drop = FALSE]
  message(sprintf("Uploading %d record(s) to %s...", nrow(rows_to_upload), database))
  DBI::dbAppendTable(con, database, value = rows_to_upload)
  message("Upload complete.")

  invisible(rows_to_upload)
}

perform_bulk_upload <- function(data, database, con, key_cols, context, preprocess = identity) {
  con <- resolve_connection(con)
  data <- validate_upload_data(data, unique(c(key_cols)), context)
  data <- preprocess(data)

  if (!nrow(data)) {
    message("Input data has 0 rows; nothing to upload.")
    return(invisible(data))
  }

  existing <- collect_existing_keys(con, database, key_cols)
  indices <- compute_new_indices(existing, data, key_cols)

  if (!length(indices)) {
    message("No new records to upload.")
    return(invisible(data[0, , drop = FALSE]))
  }

  rows_to_upload <- data[indices, , drop = FALSE]
  message(sprintf("Uploading %d record(s) to %s via bulk loader...", nrow(rows_to_upload), database))
  uploadData(rows_to_upload, database, con = con)
  message("Upload complete.")

  invisible(rows_to_upload)
}

#' Bulk upload sequence variant data using `LOAD DATA`
#'
#' @inheritParams upload_lab_data
#'
#' @return Invisibly returns the rows uploaded.
#' @export
upload_bulk_sv <- function(data, database = NULL, con = NULL) {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when bulk uploading sequence variant data.", call. = FALSE)
  }

  perform_bulk_upload(
    data = data,
    database = database,
    con = con,
    key_cols = c("MetagenNumber", "SV"),
    context = "sequence variant bulk upload"
  )
}

#' Bulk upload taxonomy data using `LOAD DATA`
#'
#' @inheritParams upload_lab_data
#'
#' @return Invisibly returns the rows uploaded.
#' @export
upload_bulk_tax <- function(data, database = NULL, con = NULL) {
  if (is.null(database) || !nzchar(database)) {
    stop("`database` must be provided when bulk uploading taxonomy data.", call. = FALSE)
  }

  perform_bulk_upload(
    data = data,
    database = database,
    con = con,
    key_cols = "SV",
    context = "taxonomy bulk upload",
    preprocess = remove_carriage_returns
  )
}

#' Bulk upload arbitrary data using `LOAD DATA`
#'
#' @param data A data frame of rows to upload.
#' @param tableName The fully qualified table name.
#' @param con An existing database connection. If `NULL`, a cached connection will be used.
#' @param use_transaction Whether to wrap the bulk upload in a transaction.
#'
#' @return Invisibly returns the uploaded rows.
#' @export
#' @importFrom DBI dbExecute dbQuoteIdentifier dbQuoteString dbWithTransaction
#' @importFrom utils write.csv
uploadData <- function(data, tableName, con = NULL, use_transaction = FALSE) {
  con <- resolve_connection(con)

  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if (!nrow(data)) {
    message("Input data has 0 rows; nothing to upload.")
    return(invisible(data))
  }

  temp_file <- tempfile(pattern = "phylosql-upload-", fileext = ".csv")
  on.exit(unlink(temp_file), add = TRUE)

  utils::write.csv(data, temp_file, row.names = FALSE, quote = TRUE, na = "")

  query <- paste(
    "LOAD DATA LOCAL INFILE",
    DBI::dbQuoteString(con, normalizePath(temp_file, winslash = "/")),
    "INTO TABLE",
    DBI::dbQuoteIdentifier(con, tableName),
    "FIELDS TERMINATED BY ','",
    "ENCLOSED BY '""'",
    "LINES TERMINATED BY '\\n'",
    "IGNORE 1 LINES"
  )

  if (isTRUE(use_transaction)) {
    DBI::dbWithTransaction(con, DBI::dbExecute(con, query))
  } else {
    DBI::dbExecute(con, query)
  }

  invisible(data)
}
