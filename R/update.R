


#' A phylosql Function
#'
#' function to delete to mysql database
#' @param data data to upload
#' @param database database to send data
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

delete_data_by_sample<-
  function(samples,database="labdata", con=NULL){
    if(is.null(con)){
      stop("You need to specify a database connection")
    }

    si<- dplyr::as_tibble(
      dplyr::tbl(con,database))

    query  <-  sprintf("DELETE FROM %s WHERE MetagenNumber IN (%s)",  database, paste0(samples,collapse=', '))
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    RMariaDB::dbExecute(con, query)
    dbDisconnect(con)
    message("Complete.")

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

  function(newdata,database,con=NULL){

    if(is.null(con)|is.null(database)){
      stop("You need to specify a database connection and a database")
    }

    match_idx = match(colnames(newdata),c('MetagenNumber','Factor','Level'))

    if(length(match_idx)==3 & sum(is.na(match_idx))==0){

    delete_data_by_sample(con=con, database=database,samples=unique(newdata$MetagenNumber))

    phylosql::upload_cms_data_Long(con=con, data=newdata)

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

  function(newdata,database,con){

    if(is.null(con)|is.null(database)){
      stop("You need to specify a database connection and a database")
    }
    match_idx = match(colnames(newdata),c('MetagenNumber','SV','Abundance'))

    if(length(match_idx)==3 & sum(is.na(match_idx))==0){


    delete_data_by_sample(con=con, database=database,samples=unique(newdata$MetagenNumber))

    phylosql:::upload_bulk_sv(con=con, data=newdata)

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

  function(newdata,database,con=NULL){

    if(is.null(con)|is.null(database)){
      stop("You need to specify a database connection and a database")
    }

    match_idx = match(colnames(newdata),c('MetagenNumber','value','variable'))

    if(length(match_idx)==3 & sum(is.na(match_idx))==0){

    delete_data_by_sample(con=con, database=database,samples=unique(newdata$MetagenNumber))

    phylosql::upload_lab_data(con=con, data=newdata)

    }else{

      print('No upload as columns did not match database requirements')
    }


  }
