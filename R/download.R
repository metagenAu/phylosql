
#' A phylosql Function
#'
#'
#' @param si_long data to upload
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import Matrix
#' @export
#'

create_sampleInfo_table<- function(si_long, ... ){

  vars<- unique(si_long$variable)
  samples<- unique(si_long$MetagenNumber)
  rows<- length(samples)
  cols<- length(vars)
  sample_data<- matrix(NA,ncol=cols, nrow=rows)

  sampleList<- split(si_long,si_long$MetagenNumber)

  for(i in seq_along(sampleList)){
    sample_data[i,match(sampleList[[i]]$variable, vars)  ]<- sampleList[[i]]$value
  }

  rownames(sample_data)<- names(sampleList)
  colnames(sample_data)<- vars
  sample_data<- data.frame(sample_data)
  sample_data$MetagenNumber<- samples

  return(sample_data)
}

#' A phylosql Function
#'
#'
#' @param flist
#' @param con connection
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

  fetch_sampleInfo<-
  function(flist=NULL, con=NULL,cms_format="long"){
    if(is.null(con)){

      con <-  try_fetch_connection()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    si<- dplyr::as_tibble(dplyr::tbl(con,"labdata"))
    sample_info<- create_sampleInfo_table(si_long=si)

    if(cms_format=="long"){

      cmslong<- dplyr::as_tibble(dplyr::tbl(con,"cmsdatalong"))
      cms<- create_cms_table(cms_long=cmslong)

    }else{

      cms<- data.frame(dplyr::tbl(con,"cmsdata"))

    }

    sampleInfo <- merge(sample_info,cms,"MetagenNumber")

    missed<- which(!cms$MetagenNumber %in% sampleInfo$MetagenNumber)
    cols<- match(colnames(cms), colnames(sampleInfo))

    nsamps<- nrow(sampleInfo)

    sampleInfo[(nsamps+1):(nsamps+length(missed)),cols] <- cms[missed, ]

    if(!is.null(flist)){
      sampleInfo<- subset(sampleInfo,flist)
    }

   # class(sampleInfo)<- c(class(sampleInfo), "sampledata")
    rownames(sampleInfo)<- sampleInfo$MetagenNumber
    return(sampleInfo)

  }



#' A phylosql Function
#'
#'
#' @param phylo logical
#' @param database database to send data
#' @param con connection
#' @param whichSamples select specific samples to access
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import Matrix
#' @import magrittr
#' @export
#'

fetch_asv_table_sparse<- function(con=NULL,database="eukaryota_sv",phylo=FALSE, whichSamples=NULL ){

  if(is.null(con)){

    con <-  try_fetch_connection()

  }


  if(any(class(con)=='logical')){

    stop('No connection to database.')

  }
# fetch data
asv_long<- dplyr::as_tibble(dplyr::tbl(con,database))

if(!is.null(whichSamples)){
  asv_long<- asv_long %>%
    filter(!!asv_long$MetagenNumber %in% whichSamples )
}

asv_long$MetagenNumber<- as.factor(asv_long$MetagenNumber)
asv_long$SV<- as.factor(asv_long$SV)

asv_table = Matrix::sparseMatrix(i = asv_long$MetagenNumber %>% as.integer,
                                 j = asv_long$SV %>% as.integer,
                                 x = asv_long$Abundance)
rownames(asv_table) = levels(asv_long$MetagenNumber)
colnames(asv_table) = levels(asv_long$SV)

rm(asv_long)

return(asv_table)

}



#' A phylosql Function
#'
#'
#' @param phylo logical
#' @param database database to send data
#' @param con connection
#' @param whichSamples select specific samples to access
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import Matrix
#' @import magrittr
#' @export
#'

fetch_asv_table<- function(con=NULL,database="eukaryota_sv",phylo=FALSE, whichSamples=NULL ){

  if(is.null(con)){

    con <-  try_fetch_connection()

  }


  if(any(class(con)=='logical')){

    stop('No connection to database.')

  }
  # fetch data
  asv_long<- dplyr::as_tibble(dplyr::tbl(con,database))

  if(!is.null(whichSamples)){
    asv_long<- asv_long %>%
      filter(!!asv_long$MetagenNumber %in% whichSamples )
  }

  # Construct table
  asvs<- unique(asv_long$SV)
  samples<- unique(asv_long$MetagenNumber)
  rows<- length(samples)
  cols<- length(asvs)
  asv_table<- matrix(0L,ncol=cols, nrow=rows)

  sampleList<- split(asv_long,asv_long$MetagenNumber)

  for(i in seq_along(sampleList)){
    asv_table[i,match(sampleList[[i]]$SV, asvs)  ]<- sampleList[[i]]$Abundance
  }

  rownames(asv_table) <- names(sampleList)
  colnames(asv_table)<- asvs

  rm(asvlong)
  rm(sampleList)


  # Set class (or not)
 # if(phylo==FALSE){
  asv_table<-  as(asv_table,"dgCMatrix")
  gc()
  #attr(asv_table, "type")<- "abundance"
 # }

  return(asv_table)

}

#' A phylosql Function
#'
#' @param phylo logical. Whether to format data for phyloseq or sparseHDD.
#' @param database database to send data
#' @param con connection
#' @param whichTaxa select specific taxa
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'


fetch_taxonomy<- function(con=NULL, database="eukaryota_tax",whichTaxa=NULL, phylo=FALSE){

  if(is.null(con)){

    con <-  try_fetch_connection()

  }


  if(any(class(con)=='logical')){

    stop('No connection to database.')

  }


  tax<- dplyr::as_tibble(dplyr::tbl(con,database))
  tax <- dplyr::mutate_if(tax,
                    is.character,
                    stringr::str_replace_all, pattern = "\\r", replacement = "")

  if(!is.null(whichTaxa)){
    tax<- tax %>%
      dplyr::filter(!!tax$SV %in% whichTaxa )
  }

  if(phylo==FALSE){

    attr(tax,"type")<- c( "taxonomy")

  }else{
    SV<- tax$SV
    tax<- tax[,-1]
    tax<- as.matrix(tax)
    tax<- gsub("NA",NA,tax)

    rownames(tax)<- SV
  }
  gc()
   return(tax)

}


#' A phylosql Function
#'
#'
#' @param cms_long data to upload
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

create_cms_table<- function(cms_long, ...){

  vars<- unique(cms_long$Factor)
  samples<- unique(cms_long$MetagenNumber)
  rows<- length(samples)
  cols<- length(vars)
  sample_data<- matrix(NA,ncol=cols, nrow=rows)

  sampleList<- split(cms_long,cms_long$MetagenNumber)

  for(i in seq_along(sampleList)){
    sample_data[i,match(sampleList[[i]]$Factor, vars)  ]<- sampleList[[i]]$Level
  }

  rownames(sample_data)<-  names(sampleList)
  colnames(sample_data)<- vars
  sample_data<- data.frame(sample_data)
  sample_data$MetagenNumber<- samples

  return(sample_data)

}






#' A phylosql Function
#'
#'
#' @param samples
#' @param database
#' @param con
#' @param col
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
fetch_asvs_by_sample<-
  function(samples,database=NULL, con=NULL,col='MetagenNumber'){
    if(is.null(con)){
      stop("You need to specify a database connection")
    }
    samples = unique(samples)

    query  <-  sprintf("SELECT * FROM %s WHERE %s IN (%s)",  database,col, paste0(add_quotes(samples),collapse=', '))
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    df <- dbGetQuery(con, query)
    #dbDisconnect(con)
    message("Complete.")
    df

  }

#' A phylosql Function
#'
#'
#' @param phylo logical
#' @param database database to send data
#' @param con connection
#' @param whichSamples select specific samples to access
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import Matrix
#' @import magrittr
#' @export
#'
fetch_asv_table_sparse_by_sample<-

  function(con=NULL,database="eukaryota_sv",phylo=FALSE, whichSamples=NULL ){


    if(is.null(database)){
      stop("You need to specify a database connection and a database")
    }

    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }


  if(is.null(whichSamples)){
    stop("You need to specify samples.")
  }
  asv_long<- fetch_asvs_by_sample(samples=whichSamples,con=eval_con(con),database=database)

  gc()
  asv_long$MetagenNumber<- as.factor(asv_long$MetagenNumber)
  asv_long$SV<- as.factor(asv_long$SV)

  asv_table = Matrix::sparseMatrix(i = asv_long$MetagenNumber %>% as.integer,
                              j = asv_long$SV %>% as.integer,
                              x = asv_long$Abundance)
  rownames(asv_table) = levels(asv_long$MetagenNumber)
  colnames(asv_table) = levels(asv_long$SV)

  rm(asv_long)

  gc()

  return(asv_table)

}

#' A phylosql Function
#'
#'
#' @param taxa
#' @param database
#' @param con
#' @param col
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
fetch_taxonomy_by_asv<-
  function(taxa,database=NULL, con=NULL,col='SV'){
    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    taxa<- unique(taxa)

    query  <-  sprintf("SELECT * FROM %s WHERE %s IN (%s)",  database,col, paste0(add_quotes(taxa),collapse=', '))
    # SUBMIT THE UPDATE QUERY AND DISCONNECT
    res <- dbGetQuery(con, query)
    #dbDisconnect(con)
    message("Complete.")
    tax <- dplyr::mutate_if(df,
                            is.character,
                            stringr::str_replace_all, pattern = "\\r", replacement = "")


    if(phylo==FALSE){

      attr(tax,"type")<- c( "taxonomy")

    }else{
      SV<- tax$SV
      tax<- tax[,-1]
      tax<- as.matrix(tax)
      tax<- gsub("NA",NA,tax)

      rownames(tax)<- SV
    }
    gc()
    return(tax)

  }


#' A phylosql Function
#'
#'
#' @param database
#' @param con
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'
get_svs<-
  function(database=NULL, con=NULL){

    if(is.null(con)){

      con <-  try_fetch_connection()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    query  <-  sprintf("SELECT `SV` FROM %s",  database)
    df<- dbGetQuery(con, query)
    unique(df$SV)

}


#' A phylosql Function to add quotes within an SQL query
#'
#'
#' @param x
#' @keywords
#' @export
#'
add_quotes<- function(x){

  x <- unique(x)
  vec<- vector('list',length=length(x))
  for(i in seq_along(x)){
    x0<-  paste0("'",x[i],"'")
    vec[[i]]<- x0
  }
  unlist(vec)

}


#' A phylosql Function to construct an asv table
#'
#'
#' @param asv_long
#' @keywords
#' @export
#'
construct_asv_table<- function(asv_long){

  asvs <- unique(asv_long$SV)
  samples <- unique(asv_long$MetagenNumber)
  rows <- length(samples)
  cols <- length(asvs)
  asv_table <- matrix(0L, ncol = cols, nrow = rows)
  sampleList <- split(asv_long, asv_long$MetagenNumber)
  for (i in seq_along(sampleList)) {
    asv_table[i, match(sampleList[[i]]$SV, asvs)] <- sampleList[[i]]$Abundance
  }
  rownames(asv_table) <- names(sampleList)
  colnames(asv_table) <- asvs
  rm(asv_long)
  rm(sampleList)
  asv_table <- as(asv_table, "dgCMatrix")
  gc()
  return(asv_table)

}

#' A phylosql Function
#'
#'
#' @param phylo logical
#' @param database database to send data
#' @param con connection
#' @param whichSamples select specific samples to access
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import Matrix
#' @import magrittr
#' @export
#'
fetch_asv_table_by_sample<-

  function(con=NULL,database="eukaryota_sv",phylo=FALSE, whichSamples=NULL ){


    if(is.null(database)){
      stop("You need to specify a database connection and a database")
    }

    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    # fetch data

    if(is.null(whichSamples)){
      stop("You need to specify samples.")
    }
    asv_long<- fetch_asvs_by_sample(samples=whichSamples,con=eval_con(con),database=database)


    asv_table <- construct_asv_table(asv_long)

    gc()

    return(asv_table)

  }

