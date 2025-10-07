.phylosql_state <- new.env(parent = emptyenv())

load_connection_credentials <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    stop("`path` must be a non-empty string to a credentials CSV file.", call. = FALSE)
  }

  if (!file.exists(path)) {
    stop(sprintf("Credential file does not exist: %s", path), call. = FALSE)
  }

  cols <- c(host = "character", dbname = "character", port = "integer", user = "character")

  creds <- tryCatch(
    utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      nrows = 1,
      colClasses = cols
    ),
    error = function(err) {
      stop(sprintf("Unable to read credential file %s: %s", path, conditionMessage(err)), call. = FALSE)
    }
  )

  missing_cols <- setdiff(names(cols), names(creds))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Credential file %s is missing required column(s): %s",
        path,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  creds[1, names(cols), drop = FALSE]
}

#' Establish and cache a connection to the MTGN database
#'
#' @param path Path to a CSV containing connection credentials with columns `host`, `dbname`, `port`, and `user`.
#' @param key Secret corresponding to the credentials in `path` (e.g. password).
#' @param refresh Whether to force reconnection even if a cached connection exists.
#'
#' @return A live `DBIConnection` object.
#' @export
#' @importFrom DBI dbConnect dbDisconnect dbIsValid
#' @importFrom RMariaDB MariaDB
get_mtgn_connection <- function(path = NULL, key = NULL, refresh = FALSE) {
  creds_cached <- .phylosql_state$credentials
  key_cached <- .phylosql_state$key
  con_cached <- .phylosql_state$connection

  if (!isTRUE(refresh) && !is.null(con_cached) && DBI::dbIsValid(con_cached)) {
    return(con_cached)
  }

  if (is.null(path)) {
    if (is.null(creds_cached)) {
      stop("No cached connection available; provide `path` and `key`.", call. = FALSE)
    }
    path <- creds_cached$path
  }

  if (is.null(key)) {
    if (is.null(key_cached)) {
      stop("No cached credentials available; provide `key`.", call. = FALSE)
    }
    key <- key_cached
  }

  if (!nzchar(key)) {
    stop("`key` must be a non-empty password or access token.", call. = FALSE)
  }

  creds_df <- load_connection_credentials(path)

  if (!is.null(con_cached) && DBI::dbIsValid(con_cached)) {
    try(DBI::dbDisconnect(con_cached), silent = TRUE)
  }

  message("Establishing database connection...")
  con <- DBI::dbConnect(
    drv = RMariaDB::MariaDB(),
    host = creds_df$host,
    dbname = creds_df$dbname,
    port = creds_df$port,
    user = creds_df$user,
    password = key
  )
  message("Connection established.")

  .phylosql_state$connection <- con
  .phylosql_state$credentials <- list(path = path, details = creds_df)
  .phylosql_state$key <- key

  con
}

#' Establish a pooled connection to the MTGN database
#'
#' @inheritParams get_mtgn_connection
#' @param size Optional pool size configuration passed to [pool::dbPool()].
#'
#' @return A live `pool::Pool` object.
#' @export
#' @importFrom pool dbPool poolClose poolClosed
set_pool <- function(path, key, size = NULL) {
  creds_df <- load_connection_credentials(path)

  if (!nzchar(key)) {
    stop("`key` must be a non-empty password or access token.", call. = FALSE)
  }

  existing_pool <- .phylosql_state$pool
  if (!is.null(existing_pool) && !pool::poolClosed(existing_pool)) {
    pool::poolClose(existing_pool)
  }

  pool_args <- list(
    drv = RMariaDB::MariaDB(),
    host = creds_df$host,
    dbname = creds_df$dbname,
    port = creds_df$port,
    user = creds_df$user,
    password = key
  )

  if (!is.null(size)) {
    pool_args <- c(pool_args, size)
  }

  pool <- do.call(pool::dbPool, pool_args)

  .phylosql_state$pool <- pool
  .phylosql_state$credentials <- list(path = path, details = creds_df)
  .phylosql_state$key <- key

  pool
}

#' Attempt to reuse a cached connection or pool
#'
#' @return A live database connection object or `NULL` if none are cached.
#' @export
#' @importFrom DBI dbIsValid
try_fetch_connection <- function() {
  pool <- .phylosql_state$pool
  if (!is.null(pool) && !pool::poolClosed(pool)) {
    return(pool)
  }

  con <- .phylosql_state$connection
  if (!is.null(con) && DBI::dbIsValid(con)) {
    return(con)
  }

  NULL
}

#' Return the expression used to lazily fetch a connection
#'
#' @return A character string representing a call to [get_mtgn_connection()].
#' @export
try_fetch_connection_string <- function() {
  "get_mtgn_connection()"
}

#' Evaluate a connection expression in the caller environment
#'
#' @param expr Character string of R code that resolves to a connection when evaluated.
#'
#' @return The evaluated object (typically a connection).
#' @export
eval_con <- function(expr) {
  if (missing(expr) || !is.character(expr) || length(expr) != 1L) {
    stop("`expr` must be a single character string.", call. = FALSE)
  }

  eval(parse(text = expr), envir = parent.frame())
}
