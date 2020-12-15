
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

    existingID<- paste0(si$MetagenNumber)
    newID<- paste0(data$MetagenNumber)

    upload<- which(!newID %in% existingID)
    stopifnot(length(upload)>0)
    message(paste0("Uploading ",length(upload)," samples."))
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }
