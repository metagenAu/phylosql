
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
      stop("You need to specify a database connection")
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
    if (is.null(con)) {
      stop("You need to specify a database connection")
    }
    if (is.null(database)) {
      stop("You need to specify a database")
    }
    print(con)
    sv <- dplyr::as_tibble(dplyr::tbl(con, database))
    existingID <- paste0(sv$MetagenNumber, sv$SV)
    newID <- paste0(data$MetagenNumber, data$SV)
    upload <- which(!newID %in% existingID)
    stopifnot(length(upload) > 0)
    message(paste0("Uploading ", length(upload), " samples."))
    uploadData(data=data[upload,],database,con=con)
    message("Complete.")
    dbDisconnect(con)

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
    if (is.null(con)) {
      stop("You need to specify a database connection")
    }
    if (is.null(database)) {
      stop("You need to specify a database")
    }
    tax <- dplyr::as_tibble(dplyr::tbl(con, database))
    existingID <- paste0(tax$SV)
    newID <- paste0(data$SV)
    upload <- which(!newID %in% existingID)
    stopifnot(length(upload) > 0)
    message(paste0("Uploading ", length(upload), " samples."))
    print(paste0('Con file is: ',con))
    uploadData(data=data[upload,],database,con=con)
    message("Complete.")
    dbDisconnect(con)

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
           con) # arguments to DBI::dbConnect
  {
   # TEMPFILE  <-  write.csv(fileext='.csv')
   # TEMPFILE<- normalizePath(TEMPFILE, winslash = "/")
    TEMPFILE = 'bulk_upload1.csv'
    query  <-  sprintf("LOAD DATA LOCAL INFILE '%s'
INTO TABLE %s
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\\n'
IGNORE 1 LINES;" , TEMPFILE,tableName)
    print(query)

    write.csv(data,TEMPFILE, row.names = FALSE,quote = FALSE)
    #on.exit(file.remove(TEMPFILE))

    # CONNECT TO THE DATABASE
    print(con)

    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    RMariaDB::dbExecute(con, query)
    dbDisconnect(con)
    on.exit(file.remove(TEMPFILE))
  }

