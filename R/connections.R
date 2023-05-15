


#' A phylosql Function
#'
#'
#' @param path
#' @param key
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import DBI
#' @export
#'
get_mtgn_connection<-
  function(path=NULL,key=NULL){

    if(exists('sql_creds')){
      print('Fetching Con From Pool')
      return(sql_creds$creds)

    }else{

      file<-
        tryCatch({
          read.csv( path, header=T)
        },
        error=
          function(e){
            return(NA)
          })

      if(class(file)=='data.frame'){

        message("Fetching connection...")
        con<-
          DBI::dbConnect(
            RMariaDB::MariaDB(),
            host=file$host,
            dbname=file$dbname,
            port=file$port,
            user=file$user,
            password=key)
        message("Complete ;)")
        return(con)

      }else{

        stop("Oops! No secret key found.")

      }

    }

  }



#' A phylosql Function
#'
#'
#' @param path
#' @param key
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import pool
#' @import DBI
#' @export
#'

set_pool<- function(path,key){

  sql_creds <- new.env()
  file = read.csv(path)

  pool <- dbPool(
    drv= RMariaDB::MariaDB(),
    host = file$host,
    dbname = file$dbname,
    port = file$port,
    user = file$user,
    password = key
  )

  sql_creds$creds <- pool
  assign("sql_creds", sql_creds, .GlobalEnv)
  #as.environment('package:phylosql')
  # assign("sql_creds", sql_creds, 'package:phylosql')


}


#' A phylosql Function
#'
#'
#' @param path
#' @param key
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import pool
#' @import DBI
#' @export
#'

try_fetch_connection<-

  function(){

  get_mtgn_connection()

  }

#' A phylosql Function
#'
#'
#' @param path
#' @param key
#' @keywords
#' @export
#'

try_fetch_connection_string<-

  function(){

    'get_mtgn_connection()'

  }

#' A phylosql Function
#'
#'
#' @param path
#' @param key
#' @keywords
#' @export
#'
eval_con<- function(x){

  eval(parse(text = paste0(x)))
}
