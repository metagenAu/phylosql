


#' A phylosql Function
#'
#'
#' @param taxalevel
#' @param taxa_group
#' @param con
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import phyloseqSparse
#' @export
#'
sql_phyloseq_by_tax_glom<-
  function(taxalevel,taxa_group, con=NULL ){

    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(class(con)!= 'character'){
      'Stop'
    }

    sv = paste0(taxa_group,'_sv')
    tax = paste0(taxa_group,'_tax')

    sv_cols = c('SV','MetagenNumber') #,'Abundance')
    sv_cols  = paste0(sv,'.',sv_cols )

    res <- DBI::dbGetQuery(
      con=eval_con(con),
      paste0(
        "SHOW COLUMNS FROM ",
        tax))

    id = match(taxalevel,res$Field)
    tax_base =  res$Field[1:id]
    tax_cols = paste0(tax,'.',tax_base)
    tax_cols = tax_cols[!grepl('SV',tax_cols)]
    print(tax_cols)

    grouping_cols =  c(tax_cols,'SV_rep') %>% paste(collapse=', ')
    select_cols = c(tax_cols,paste0(sv,'.MetagenNumber')) %>% paste(collapse=', ')

  # alt way of writing queries. Easier to read the query but harder to track the inputs.
  #  query =
  #    sprintf('SELECT %s , SUM( %s.Abundance) \n FROM %s \n LEFT JOIN %s ON %s.SV = %s.SV \n GROUP BY %s ;',
  #            select_cols , sv, sv, tax, sv,tax, grouping_cols)

    #thesv = paste0(tax,'.SV')
   # therep = sprintf(', MIN(%s) AS SV_rep ', thesv)
   # SELECT t1.column1, SUM(t1.column2) AS sum_column2, t2.column5, MAX(t2.column4) AS max_column4
  #  query2 =
  #    paste0(c(
  #      paste0('SELECT ',select_cols, therep, ', SUM( ',sv , '.Abundance) AS `Abundance, `',collapse=' '),
  #      paste0('FROM ',sv,collapse=' ') ,
 #       paste0('LEFT JOIN ',tax, ' ON ',sv,'.SV =',tax,'.SV',collapse=' '),
 #       paste0('GROUP BY ',grouping_cols,';',collapse=' ')) ,
 #       collapse='\n'
 #     )

    query =
      paste0(c(
        paste0('SELECT ',select_cols, ', SUM( ',sv , '.Abundance) AS `Abundance`',collapse=' '),
        paste0('FROM ',sv,collapse=' ') ,
        paste0('LEFT JOIN ',tax, ' ON ',sv,'.SV =',tax,'.SV',collapse=' '),
        paste0('GROUP BY ',select_cols,';',collapse=' ')) ,
        collapse='\n'
      )


    results <- DBI::dbGetQuery(
      con=eval_con(con),
      query)

    results = results %>% dplyr::filter(Abundance>0)

    print(paste0('results: ',dim(results))) 
    print(head(results))

    select_cols_tax = c(paste0(tax,'.',tax_base)) %>% paste(collapse=', ')
    tax_query =
      paste0(c(
        paste0('SELECT ', select_cols_tax,collapse=' '),
        paste0('FROM ',tax,collapse=' ') ),
        collapse='\n')


    tax_results <-
      DBI::dbGetQuery(
      con=eval_con(con),
      tax_query)
          
    print(paste0('results: ',dim(tax_results))) 
    print(head(tax_results))

    query_id = apply(results %>% dplyr::select(-MetagenNumber,-Abundance),1,function(x)paste0(x,collapse=';'))

    tax_id = apply(tax_results %>% dplyr::select(-SV),1,function(x)paste0(x,collapse=';'))
    print(paste0('match query: ',length(match(query_id,tax_id))))
    results$SV<- tax_results$SV[match(query_id,tax_id)]
    rm(tax_results)
    gc()
    results = results[-which(is.na(results$SV)),]

    sv_keep = c('SV','MetagenNumber','Abundance')

    asv_long = results[ , sv_keep ]
    print(paste0('asv dim: ',dim(asv_long)))

    asv_long$MetagenNumber<- as.factor(asv_long$MetagenNumber)
    asv_long$SV<- as.factor(asv_long$SV)

    asv_table = Matrix::sparseMatrix(i = asv_long$MetagenNumber %>% as.integer,
                                     j = asv_long$SV %>% as.integer,
                                     x = asv_long$Abundance)

    print(paste0('asv dim: ' ,dim(asv_table)))
               
    rownames(asv_table) = levels(asv_long$MetagenNumber)
    colnames(asv_table) = levels(asv_long$SV)

    rm(asv_long)
    gc()

    tax<- results[,  tax_base ]
    tax = tax[!duplicated(tax$SV),]
    print(paste0('tax dim: ' ,dim(tax)))

    rownames(tax)<- tax$SV
    tax = tax %>% dplyr::select(-SV)
    tax = gsub('\\\r','',as.matrix(tax))
    print(paste0('tax dim: ' ,dim(tax)))

    si<- fetch_sampleInfo(con=eval_con(con))
    rm(results)
    gc()
    phyloseqSparse::phyloseq(
      phyloseqSparse::otu_table(
        asv_table,
        taxa_are_rows=FALSE),
      phyloseqSparse::tax_table(
        as.matrix(tax)
      ),
      phyloseqSparse::sample_data(si)
    ) %>% return()


  }

#' A phylosql Function.
#'
#'
#' @param samples
#' @param taxa_group
#' @param con
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import phyloseqSparse
#' @export
#'

sql_phyloseq_by_sample<-
  function(taxa_group, samples, con =NULL ){

    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }

    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(class(con)!= 'character'){
      stop('Connection is not a character a string')
    }

    sv = paste0(taxa_group,'_sv')
    tax = paste0(taxa_group,'_tax')

    res <- dbGetQuery(
      con=eval_con(con),
      paste0(
        "SHOW COLUMNS FROM ",
        tax))

    tax_base =  res$Field

    sv_cols = c('SV','MetagenNumber') #,'Abundance')
    sv_cols  = paste0(sv,'.',sv_cols )



    query =
      paste0(
        c(
          paste0('SELECT *',collapse=' ') ,
          paste0('FROM ',sv,collapse=' ') ,
          paste0('LEFT JOIN ',tax, ' ON ',sv,'.SV = ',tax,'.SV ',collapse=' '),
          paste0( 'WHERE MetagenNumber IN (',sprintf('%s', paste0(add_quotes(samples),collapse=', ')), ')',collapse=', ')),

        collapse='\n'
      )


    results <- dbGetQuery(
      con=eval_con(con), query)

    results = results %>% dplyr::filter(Abundance>0)


    sv_keep = c('SV','MetagenNumber','Abundance')

    asv_long = results[ , sv_keep ]

    asv_long$MetagenNumber<- as.factor(asv_long$MetagenNumber)
    asv_long$SV<- as.factor(asv_long$SV)

    asv_table = Matrix::sparseMatrix(i = asv_long$MetagenNumber %>% as.integer,
                                     j = asv_long$SV %>% as.integer,
                                     x = asv_long$Abundance)

    rownames(asv_table) = levels(asv_long$MetagenNumber)
    colnames(asv_table) = levels(asv_long$SV)

    rm(asv_long)
    gc()

    tax<- results[,  tax_base ]
    tax = tax[!duplicated(tax$SV),]
    rownames(tax)<- tax$SV
    tax = tax %>% dplyr::select(-SV)
    tax = gsub('\\\r','',as.matrix(tax))

    si<- fetch_sampleInfo(con=eval_con(con))

    phyloseqSparse::phyloseq(
      phyloseqSparse::otu_table(
        asv_table,
        taxa_are_rows=FALSE),
      phyloseqSparse::tax_table(
        as.matrix(tax)
      ),
      phyloseqSparse::sample_data(si)
    ) %>% return()

  }



#' A phylosql Function. this function is unfinished. Do not use.
#'
#'
#' @param taxalevel
#' @param taxa_group
#' @param samples
#' @param con
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import phyloseqSparse
#' @export
#'

sql_phyloseq_by_sample_and_tax_glom<-
  function(taxa_group, taxalevel ,samples, con =NULL ){

    if(is.null(con)){

      con <-  try_fetch_connection_string()

    }


    if(any(class(con)=='logical')){

      stop('No connection to database.')

    }

    if(class(con)!= 'character'){
      stop('Connection is not a character a string')
    }

    sv = paste0(taxa_group,'_sv')
    tax = paste0(taxa_group,'_tax')

    sv_cols = c('SV','MetagenNumber') #,'Abundance')
    sv_cols  = paste0(sv,'.',sv_cols )

    res <- dbGetQuery(
      con=eval_con(con),
      paste0(
        "SHOW COLUMNS FROM ",
        tax))

    id = match(taxalevel,res$Field)
    tax_base =  res$Field[1:id]
    tax_cols = paste0(tax,'.',tax_base)
    tax_cols = tax_cols[!grepl('SV',tax_cols)]

    grouping_cols =  c(tax_cols) %>% paste(collapse=', ')
    select_cols = c(tax_cols,sv_cols) %>% paste(collapse=', ')

    query =
      paste0(
        c(
          paste0('SELECT ',select_cols, ', SUM( ',sv , '.Abundance) AS `Abundance`',collapse=' '),
          paste0('FROM ( SELECT *',collapse=' ') ,
          paste0('FROM ',sv,collapse=' ') ,
          paste0( 'WHERE MetagenNumber IN (',sprintf('%s', paste0(add_quotes(samples),collapse=', ')), ')',collapse=', '),
          paste0(') as ',sv,collapse=' ') ,
          paste0('LEFT JOIN ',tax, ' ON ',sv,'.SV = ',tax,'.SV ',collapse=' '),
          paste0('GROUP BY ',grouping_cols,';',collapse=' ') ),
        collapse='\n'
      )


    results <- dbGetQuery(
      con=eval_con(con), query)

    results = results %>% dplyr::filter(Abundance>0)


    sv_keep = c('SV','MetagenNumber','Abundance')

    asv_long = results[ , sv_keep ]

    asv_long$MetagenNumber<- as.factor(asv_long$MetagenNumber)
    asv_long$SV<- as.factor(asv_long$SV)

    asv_table = Matrix::sparseMatrix(i = asv_long$MetagenNumber %>% as.integer,
                                     j = asv_long$SV %>% as.integer,
                                     x = asv_long$Abundance)

    rownames(asv_table) = levels(asv_long$MetagenNumber)
    colnames(asv_table) = levels(asv_long$SV)

    rm(asv_long)
    gc()


    tax<- results[,  tax_base ]
    tax = tax[!duplicated(tax$SV),]
    rownames(tax)<- tax$SV
    tax = tax %>% dplyr::select(-SV)
    tax = gsub('\\\r','',as.matrix(tax))

    si<- fetch_sampleInfo(con=eval_con(con))

    phyloseqSparse::phyloseq(
      phyloseqSparse::otu_table(
        asv_table,
        taxa_are_rows=FALSE),
      phyloseqSparse::tax_table(
        as.matrix(tax)
      ),
      phyloseqSparse::sample_data(si)
    ) %>% return()

  }

