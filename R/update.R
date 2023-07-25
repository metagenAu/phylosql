


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

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

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
           con=NULL,
           databases=list('bacteria_sv',
                          'eukaryota_sv',
                          'labdata',
                          'cmsdata')){

    if(is.null(con)){

      con <-  try_fetch_connection()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

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

      con <-  try_fetch_connection()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(is.null(database)|is.null(col)){
      stop('missing database or column')
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
#' A phylosql Function
#'
#' function to delete to mysql database
#' @param samples data to upload
#' @param database database to send data
#' @param con connection
#' @param vars
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

delete_labdata_by_sample_and_var<-
  function(samples,vars, database=NULL, con=NULL){
    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    samples = unique(samples)

    queries =
      sprintf("DELETE FROM %s WHERE MetagenNumber = `%s` AND variable = `%s` ",
            'labdata',
            samples,
            vars

    )
    for( query in queries){

      RMariaDB::dbExecute(con, query)

    }
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    # dbDisconnect(con)
    message("Complete.")

  }
