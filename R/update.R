


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
  function(samples,database=NULL, con=NULL){
    if(is.null(con)){
      stop("You need to specify a database connection")
    }

   # si<- dplyr::as_tibble(
   #   dplyr::tbl(con,database))

    samples = unique(samples)

    query  <-  sprintf("DELETE FROM %s WHERE MetagenNumber IN (%s)",  database, paste0(add_quotes(samples),collapse=', '))
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    RMariaDB::dbExecute(con, query)
   # dbDisconnect(con)
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

    delete_data_by_sample(con=eval(parse(text = paste0(con))), database=database,samples=unique(newdata$MetagenNumber))

    phylosql::upload_cms_data_Long(con=eval(parse(text = paste0(con))), data=newdata)

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


    try({
      delete_data_by_sample(con=eval(parse(text = paste0(con))),
                            database=database,
                            samples=unique(newdata$MetagenNumber))
      message('Deleting existing data.')
      })

    upload_bulk_sv(con=eval(parse(text = paste0(con))),database= database, data=newdata)

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

    try(delete_data_by_sample(con=eval(parse(text = paste0(con))), database=database,samples=unique(newdata$MetagenNumber)))

    phylosql::upload_lab_data(con=eval(parse(text = paste0(con))), data=newdata)

    }else{

      print('No upload as columns did not match database requirements')
    }


  }

#' A phylosql Function
#'
#' function to update lab data
#' @param new_names data to upload
#' @param old_names data to upload
#' @param con connection
#' @param database database to send data
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
change_metagen_number<-
  function(old_names,
           new_names,
           con,
           databases=list('bacteria_sv',
                          'eukaryota_sv',
                          'labdata',
                          'cmsdata')){

    if(length(old_names)==length(new_names)){

      for( i in seq_along(databases)){

        database<- rep(databases[i],length(new_names))
        queries<-vector('list')
        for(j in seq_along(database)){


        query  <-  sprintf("UPDATE `%s` SET `MetagenNumber` = '%s' WHERE (`MetagenNumber` = '%s') ;",
                           database[j],
                           paste0(new_names[j],collapse=', '),
                           paste0(old_names[j],collapse=', ')
        )

        try(RMariaDB::dbExecute(con, query))

        }

        # SUBMIT THE UPDATE QUERY AND DISCONNECT
      #  dbDisconnect(con)
        message("Complete.")

      }


    }

  }

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

delete_data_by_sample_custom<-
  function(samples,database=NULL, con=NULL,col=NULL){
    if(is.null(con)){
      stop("You need to specify a database connection")
    }

    # si<- dplyr::as_tibble(
    #   dplyr::tbl(con,database))
    samples = unique(samples)


    query  <-  sprintf("DELETE FROM %s WHERE %s IN (%s)",  database,col, paste0(add_quotes(samples),collapse=', '))
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    RMariaDB::dbExecute(con, query)
  #  dbDisconnect(con)
    message("Complete.")

  }

