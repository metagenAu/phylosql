
#' A phylosql Function
#'
#' function to upload lab data to mysql database
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'


upload_lab_data<-
  function(data,database="labdata", con=NULL){
    if(is.null(con)){
      stop("You need to specify a database connection")
    }

    si<- dplyr::as_tibble(
      dplyr::tbl(con,database))

    existingID<- paste0(si$MetagenNumber,si$variable)
    newID<- paste0(data$MetagenNumber,data$variable)

    upload<- which(!newID %in% existingID)

    stopifnot(length(upload)>0)
    message(paste0("Uploading ",length(upload)," samples."))

    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")

  }


#' A phylosql Function
#'
#' function to upload sv table to mysql database
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

upload_sv<-
  function(data,database=NULL,con=NULL){

    if(is.null(con)){
      stop("You need to specify a database connection")
    }

    # Preprocess data for sql here

    sv<- dplyr::as_tibble(
      dplyr::tbl(con,database))

    existingID<- paste0(sv$MetagenNumber,sv$SV)
    newID<- paste0(data$MetagenNumber,data$SV)

    upload<- which(!newID %in% existingID)
    stopifnot(length(upload)>0)
    message(paste0("Uploading ",length(upload)," samples."))
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }

#' A phylosql Function
#'
#' function to upload taxonomy table to mysql database
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

upload_taxonomy<-
  function(data,database=NULL,con=NULL){
    if(is.null(con)){
      stop("You need to specify a database connection")
    }

    data<- gsub("\\r","",data)

    # Preprocess data for sql here

    tax<- dplyr::as_tibble(
      dplyr::tbl(con,database))

    existingID<- paste0(tax$SV)
    newID<- paste0(data$SV)

    upload<- which(!newID %in% existingID)
    stopifnot(length(upload)>0)
    message(paste0("Uploading ",length(upload)," samples."))
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }

#' A phylosql Function
#'
#' function to upload cms data to mysql database
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

upload_cms_data<-
  function(data,database="cmsdata",con=NULL){

    if(is.null(con)){
      stop("You need to specify a database connection")
    }


    si<- dplyr::as_tibble(
      dplyr::tbl(con,database))

    existingID<- paste0(si$MetagenNumber)
    newID<- paste0(data$MetagenNumber)

    upload<- which(!newID %in% existingID)
    stopifnot(length(upload)>0)
    message(paste0("Uploading ",length(upload)," samples."))
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }


#' A phylosql Function
#'
#' function to upload cms data to mysql database
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

upload_cms_data_Long<-
  function(data,database="cmsdatalong",con=NULL){

    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }
    if(ncol(data)!=3){
      stop("This data is not the correct format")
    }
    if(any(is.na(data$Level))){
      stop("Some cells contain NAs. Delete these and reattempt upload.")
    }

    si<- dplyr::as_tibble(
      dplyr::tbl(con,database))


    existingID<- paste0(si$MetagenNumber,si$Factor)
    newID<- paste0(data$MetagenNumber,data$Factor)
    upload<- which(!newID %in% existingID)
    stopifnot(length(upload)>0)
    message(paste0("Uploading ",length(upload)," samples."))
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }



#' A phylosql Function
#'
#' function to upload a long format SV table to mysql database (quickly)
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

upload_bulk_sv<-
  function (data, database = NULL, con = NULL)
  {
    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }
    if (is.null(database)) {
      stop("You need to specify a database")
    }
    sv <- dplyr::as_tibble(dplyr::tbl(con, database))
    existingID <- paste0(sv$MetagenNumber, sv$SV)
    newID <- paste0(data$MetagenNumber, data$SV)
    upload <- which(!newID %in% existingID)
    if(length(upload) > 0){
      message(paste0("Uploading ", length(upload), " samples."))
    uploadData(data=data[upload,],database,con=con)
    message("Complete.")
   # dbDisconnect(con)
    }

  }

#' A phylosql Function
#'
#' function to upload a taxonomy table to mysql database (quickly)
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

upload_bulk_tax<-
  function (data, database = NULL, con = NULL)
  {
    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }
    if (is.null(database)) {
      stop("You need to specify a database")
    }
    tax <- dplyr::as_tibble(dplyr::tbl(con, database))
    existingID <- paste0(tax$SV)
    newID <- paste0(data$SV)
    upload <- which(!newID %in% existingID)
    if(length(upload) > 0){
    message(paste0("Uploading ", length(upload), " samples."))
    data<- gsub("\\\r","",as.matrix(data))

    uploadData(data=data[upload,],database,con=con)
    message("Complete.")
   # dbDisconnect(con)
    }

  }



#' A phylosql Function
#'
#' Push rows from a `flagged_sv` table out to their respective SV databases.
#'
#' Reads the flagged_sv table, groups rows by their target database (determined
#' either from a routing column in flagged_sv or from a user-supplied mapping),
#' and appends only the rows whose MetagenNumber+SV combination is not already
#' present in the target. The flagged_sv table itself is not modified.
#'
#' @param source the name of the flagged source table. Defaults to "flagged_sv".
#' @param route_col name of the column in `source` that identifies the target
#'   database for each row (e.g. "target_database", "Kingdom"). The values in
#'   this column are either table names directly, or keys to look up in `mapping`.
#' @param mapping an optional named character vector translating values found in
#'   `route_col` to target table names, e.g.
#'   `c(Bacteria = "bacteria_sv", Eukaryota = "eukaryota_sv")`. If NULL, the
#'   values in `route_col` are used as table names verbatim.
#' @param targets optional character vector restricting which target tables to
#'   push to. If NULL, all targets present in the routed data are pushed.
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
push_flagged_svs <-
  function(source = "flagged_sv",
           route_col = NULL,
           mapping = NULL,
           targets = NULL,
           con = NULL){

    if(is.null(con)){
      con <- try_fetch_connection()
    }

    if(any(class(con) == 'logical')){
      stop('No connection to database.')
    }

    flagged <- dplyr::as_tibble(dplyr::tbl(con, source))

    if(nrow(flagged) == 0){
      message("flagged_sv is empty; nothing to push.")
      return(invisible(NULL))
    }

    required <- c("MetagenNumber", "SV", "Abundance")
    missing_cols <- setdiff(required, colnames(flagged))
    if(length(missing_cols) > 0){
      stop(paste0("`", source, "` is missing required columns: ",
                  paste(missing_cols, collapse = ", ")))
    }

    if(is.null(route_col)){
      stop(paste0("You need to specify `route_col`: the column in `", source,
                  "` that identifies the target database for each row. ",
                  "Available columns: ",
                  paste(colnames(flagged), collapse = ", ")))
    }

    if(!route_col %in% colnames(flagged)){
      stop(paste0("`route_col` '", route_col, "' is not a column of `",
                  source, "`."))
    }

    route_values <- as.character(flagged[[route_col]])

    if(is.null(mapping)){
      target_tbl <- route_values
    }else{
      unmapped <- setdiff(unique(route_values), names(mapping))
      if(length(unmapped) > 0){
        stop(paste0("No mapping provided for value(s): ",
                    paste(unmapped, collapse = ", ")))
      }
      target_tbl <- unname(mapping[route_values])
    }

    flagged$.target <- target_tbl
    flagged <- flagged[!is.na(flagged$.target) & flagged$.target != "", ]

    if(!is.null(targets)){
      flagged <- flagged[flagged$.target %in% targets, ]
    }

    if(nrow(flagged) == 0){
      message("No flagged rows match the requested targets.")
      return(invisible(NULL))
    }

    groups <- split(flagged, flagged$.target)

    for(tbl_name in names(groups)){
      data <- groups[[tbl_name]][, required, drop = FALSE]
      message(paste0("Pushing ", nrow(data), " row(s) to `", tbl_name, "`."))
      try(upload_bulk_sv(data = data, database = tbl_name, con = con))
    }

    message("Complete.")
    invisible(names(groups))
  }


#' A phylosql Function
#'
#'  A backend function for bulk uploading data to a mysql database
#' @param data data to upload
#' @param tableName database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
uploadData <-
  function(data, # a data frame
           tableName, # table name, possibly qualified (e.g. "my_db.customers")
           con=NULL) # arguments to DBI::dbConnect
  {
    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }
   # TEMPFILE  <-  write.csv(fileext='.csv')
   # TEMPFILE<- normalizePath(TEMPFILE, winslash = "/")
    TEMPFILE = 'bulk_upload1.csv'
    query  <-  sprintf("LOAD DATA LOCAL INFILE '%s'
INTO TABLE %s
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\\n'
IGNORE 1 LINES;" , TEMPFILE,tableName)

    write.csv(data,TEMPFILE, row.names = FALSE,quote = FALSE)
    #
    # CONNECT TO THE DATABASE
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    RMariaDB::dbExecute(con, query)
    #dbDisconnect(con)
    on.exit(file.remove(TEMPFILE))
  }

