

upload_lab_data<-
  function(data,database="labdata"){

    si<- as_tibble(tbl(con,database))

    existingID<- paste0(si$MetagenNumber,si$variable)
    newID<- paste0(data$MetagenNumber,data$variable)

    upload<- which(!newID %in% existingID)
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")

  }

upload_sv<-
  function(data,database=NULL){

    # Preprocess data for sql here

    si<- as_tibble(tbl(con,database))

    existingID<- paste0(si$MetagenNumber,si$SV)
    newID<- paste0(data$MetagenNumber,data$SV)

    upload<- which(!newID %in% existingID)
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }

upload_taxonomy<-
  function(data,database=NULL){

    # Preprocess data for sql here

    si<- as_tibble(tbl(con,database))

    existingID<- paste0(si$SV)
    newID<- paste0(data$SV)

    upload<- which(!newID %in% existingID)
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }


upload_cms_data<-
  function(data,database="cmsdata"){
    si<- as_tibble(tbl(con,database))

    existingID<- paste0(si$MetagenNumber)
    newID<- paste0(data$MetagenNumber)

    upload<- which(!newID %in% existingID)
    RMariaDB::dbAppendTable(con, database,value= data[upload,] )
    message("Complete.")
  }
