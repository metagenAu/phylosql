#' Delete rows for specific samples
#'
#' Removes records matching the supplied sample identifiers from a target
#' table. The `MetagenNumber` column is used to identify matching rows.
#'
#' @param samples Character vector of sample identifiers to delete.
#' @param database Name of the table to modify.
#' @param con Database connection or pool. When `NULL`, the cached connection is
#'   reused.
#'
#' @return Invisibly returns the number of affected rows.
#' @export

delete_data_by_sample <-
  function(samples, database = NULL, con = NULL) {
    if (is.null(database) || !nzchar(database)) {
      stop("`database` must be provided when deleting data by sample.", call. = FALSE)
    }

    if (is.null(con)) {
      con <- try_fetch_connection()
    }

    if (is.null(con)) {
      stop("No connection to database.", call. = FALSE)
    }

    samples <- unique(as.character(samples))
    samples <- samples[!is.na(samples)]

    if (!length(samples)) {
      message("No samples provided; nothing to delete.")
      return(invisible(integer()))
    }

    table_sql <- as.character(DBI::dbQuoteIdentifier(con, database))
    column_sql <- as.character(DBI::dbQuoteIdentifier(con, "MetagenNumber"))
    value_sql <- paste(DBI::dbQuoteString(con, samples), collapse = ", ")

    query <- sprintf(
      "DELETE FROM %s WHERE %s IN (%s)",
      table_sql,
      column_sql,
      value_sql
    )

    affected <- DBI::dbExecute(con, query)
    message("Complete.")

    invisible(affected)
  }





#' A phylosql Function
#'
#' function to update CMS data
#' @param newdata data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
update_cms<-

  function(newdata,database=NULL,con=NULL){

    if(is.null(database)){
      stop("You need to specify a database connection and a database")
    }
    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(class(con)!= 'character'){
      'Stop'
    }
    match_idx = match(colnames(newdata),c('MetagenNumber','Factor','Level'))

    if(length(match_idx)==3 & sum(is.na(match_idx))==0){

    delete_data_by_sample(con=eval_con(con), database=database,samples=unique(newdata$MetagenNumber))

    phylosql::upload_cms_data_Long(con=eval_con(con), data=newdata)

    }else{

      print('No upload as columns did not match database requirements')
    }



  }



#' A phylosql Function
#'
#' function to update SV data
#' @param newdata data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
update_sv<-

  function(newdata,database=NULL,con=NULL){

    if(is.null(database)){
      stop("You need to specify a database connection and a database")
    }
    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(class(con)!= 'character'){
      stop('Error with conneciton')
    }
    match_idx = match(colnames(newdata),c('MetagenNumber','SV','Abundance'))

    if(length(match_idx)==3 & sum(is.na(match_idx))==0){


    try({
      delete_data_by_sample(con=eval_con(con),
                            database=database,
                            samples=unique(newdata$MetagenNumber))
      message('Deleting existing data.')
      })

    upload_bulk_sv(con=eval_con(con),database= database, data=newdata)

     }else{

      print('No upload as columns did not match database requirements')
    }


  }





#' A phylosql Function
#'
#' function to update lab data
#' @param newdata data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
update_labdata<-

  function(newdata,database=NULL,con=NULL){

    if(is.null(database)){
      stop("You need to specify a database connection and a database")
    }

    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(class(con)!= 'character'){
      'Stop'
    }

    match_idx = match(colnames(newdata),c('MetagenNumber','value','variable'))

    if(length(match_idx)==3 & sum(is.na(match_idx))==0){

    try(
      delete_labdata_by_sample_and_var(
      con=eval_con(con),
      database=database,
      samples=newdata$MetagenNumber,
      vars= newdata$variable)
    )

    phylosql::upload_lab_data(con=eval_con(con), data=newdata)

    }else{

      print('No upload as columns did not match database requirements')
    }


  }

#' Update metagen numbers across related tables
#'
#' Iterates through one or more tables, replacing occurrences of
#' `MetagenNumber` values with new identifiers. Each table is updated in turn
#' and any execution errors are caught and reported via `try()`.
#'
#' @param old_names Character vector of existing `MetagenNumber` values.
#' @param new_names Character vector of replacement identifiers. Must be the
#'   same length as `old_names`.
#' @param con Database connection or pool. When `NULL`, the cached connection is
#'   reused.
#' @param databases Character vector of table names to update.
#'
#' @return Invisibly returns `NULL`.
#' @export
change_metagen_number <-
  function(old_names,
           new_names,
           con = NULL,
           databases = list(
             'bacteria_sv',
             'eukaryota_sv',
             'labdata',
             'cmsdata'
           )) {

    if (is.null(con)) {

      con <-  try_fetch_connection()

    }


    if (is.null(con)) {

      stop('No connection to database.', call. = FALSE)

    }

    if (length(old_names) == length(new_names)) {

      new_values <- as.character(new_names)
      old_values <- as.character(old_names)

      for (i in seq_along(databases)) {

        table_sql <- as.character(DBI::dbQuoteIdentifier(con, databases[[i]]))
        metagen_sql <- as.character(DBI::dbQuoteIdentifier(con, "MetagenNumber"))

        for (j in seq_along(new_values)) {

          query  <-  sprintf(
            "UPDATE %s SET %s = %s WHERE %s = %s;",
            table_sql,
            metagen_sql,
            DBI::dbQuoteString(con, new_values[j]),
            metagen_sql,
            DBI::dbQuoteString(con, old_values[j])
          )

          try(DBI::dbExecute(con, query))

        }

        message("Complete.")

      }


    }

  }

#' Delete rows using a custom column
#'
#' Removes records from the specified table where the provided column matches
#' any of the supplied values.
#'
#' @param samples Character vector of values to match.
#' @param database Name of the table to modify.
#' @param con Database connection or pool. When `NULL`, the cached connection is
#'   reused.
#' @param col Column name used to identify matching records.
#'
#' @return Invisibly returns the number of affected rows.
#' @export

delete_data_by_sample_custom <-
  function(samples, database = NULL, con = NULL, col = NULL) {
    if (is.null(database) || !nzchar(database) || is.null(col) || !nzchar(col)) {
      stop("`database` and `col` must be provided when deleting data.", call. = FALSE)
    }

    if (is.null(con)) {
      con <- try_fetch_connection()
    }

    if (is.null(con)) {
      stop("No connection to database.", call. = FALSE)
    }

    samples <- unique(as.character(samples))
    samples <- samples[!is.na(samples)]

    if (!length(samples)) {
      message("No values supplied; nothing to delete.")
      return(invisible(integer()))
    }

    table_sql <- as.character(DBI::dbQuoteIdentifier(con, database))
    column_sql <- as.character(DBI::dbQuoteIdentifier(con, col))
    value_sql <- paste(DBI::dbQuoteString(con, samples), collapse = ", ")

    query <- sprintf(
      "DELETE FROM %s WHERE %s IN (%s)",
      table_sql,
      column_sql,
      value_sql
    )

    affected <- DBI::dbExecute(con, query)
    message("Complete.")

    invisible(affected)
  }
#' Delete lab data for specific sample/variable pairs
#'
#' Removes rows from `labdata` (or another table) where both the
#' `MetagenNumber` and `variable` columns match the supplied vectors.
#'
#' @param samples Character vector of sample identifiers.
#' @param vars Character vector of variable names corresponding to `samples`.
#' @param database Name of the table to modify. Defaults to `"labdata"`.
#' @param con Database connection or pool. When `NULL`, the cached connection is
#'   reused.
#'
#' @return Invisibly returns the total number of affected rows.
#' @export

delete_labdata_by_sample_and_var <-
  function(samples, vars, database = "labdata", con = NULL) {
    if (is.null(con)) {
      con <- try_fetch_connection()
    }

    if (is.null(con)) {
      stop('No connection to database.', call. = FALSE)
    }

    sample_values <- as.character(samples)
    var_values <- as.character(vars)

    keep <- !(is.na(sample_values) | is.na(var_values))
    if (!any(keep)) {
      message("No sample-variable pairs supplied; nothing to delete.")
      return(invisible(integer()))
    }

    pairs <- unique(data.frame(
      MetagenNumber = sample_values[keep],
      variable = var_values[keep],
      stringsAsFactors = FALSE
    ))

    table_sql <- as.character(DBI::dbQuoteIdentifier(con, database))
    metagen_sql <- as.character(DBI::dbQuoteIdentifier(con, "MetagenNumber"))
    variable_sql <- as.character(DBI::dbQuoteIdentifier(con, "variable"))

    affected <- vapply(seq_len(nrow(pairs)), function(i) {
      query <- sprintf(
        "DELETE FROM %s WHERE %s = %s AND %s = %s",
        table_sql,
        metagen_sql,
        DBI::dbQuoteString(con, pairs$MetagenNumber[i]),
        variable_sql,
        DBI::dbQuoteString(con, pairs$variable[i])
      )
      DBI::dbExecute(con, query)
    }, integer(1))

    message("Complete.")

    invisible(sum(affected))
  }
