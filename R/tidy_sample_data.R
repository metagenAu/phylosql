#' A phylosql Function
#'
#'
#' @param lab_path path to lab data
#' @keywords
#' @import readxl
#' @export
#'


read_masterSheet<-
  function(
    lab_path="C:/Users/Chris/Metagen/Lab - master sample sheet/master_sample_sheet_v3.0.xlsx"
    ){

  cols<- c("text",rep("numeric",7),rep("text",13),rep("numeric",3))
  sdata1<- readxl::read_xlsx(lab_path, sheet="Data formatted for R",col_types=cols)

  sdata1[sdata1=="NA"]<- NA
  remove<- which(sdata1$MetagenNumber==0 | is.na(sdata1$MetagenNumber))

  if(length(remove)>0){
    sdata1<- sdata1[-remove,]
  }

  sdata1<- data.frame(sdata1)
  rownames(sdata1)<- sdata1$MetagenNumber

  return(sdata1)
}



#' A phylosql Function
#'
#' converts mastersheet data to cms data for sql upload
#' @param data data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

ms_as_cmsData<- function(data){

  data$SurveyID<- NA
  cms_cols<-  c( "MetagenNumber","PropertyName", "SurveyID",
                 "BlockName", "CropName", "SurveyDate",
                 "Location", "AgronomistName")
  data_sql<- data[ , cms_cols]
  data_sql$TargetCropTypeName <- NA
  return(data_sql)

}

#' A phylosql Function
#'
#' formats cms data for sql upload
#' @param data cms data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'

cms_as_cmsData<- function(data){

  colnames(data)[which(colnames(data)=="SurveyId")]<- "SurveyID"
  colnames(data)[which(colnames(data)=="BarcodeId")]<- "MetagenNumber"

  cms_cols<-  c( "MetagenNumber","PropertyName", "SurveyID",
                 "BlockName", "CropName", "SurveyDate",
                 "Location", "AgronomistName","TargetCropTypeName")
  data_sql<- data[ , cms_cols]
  return(data_sql)
}


#' A phylosql Function
#'
#' formats lab to for sql upload
#' @param labdata data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export
#'


as_labData<- function(labdata){

  lab_cols<-   c("pH", "MetagenNumber",
                 "ActiveCarbon", "AggregateStability",
                 "Phosphatase", "B_Glucosidase",
                 "DNA.conc...ng.ul." ,"SoilMoisture"
                 # ,"SoilWW","SoilDW"
  )


  lab_sql<- labdata[ , lab_cols]

  stopifnot(ncol(lab_sql)==(length(lab_cols)))

  #Rename
  to_change<- c("DNA.conc...ng.ul.")
  new_names<- c("DNAConc")
  data_to_rename<- match(to_change,colnames(lab_sql))
  colnames(lab_sql)[data_to_rename]<- new_names

  lab_sql$SoilMoisture<- as.numeric(as.character(lab_sql$SoilMoisture))
  lab_sql$ActiveCarbon[which(lab_sql$ActiveCarbon==0)]<- NA
  lab_sql$pH[which(lab_sql$pH==0)]<- NA
  lab_sql$Phosphatase[which(lab_sql$Phosphatase==0)]<- NA
  lab_sql$B_Glucosidase[which(lab_sql$B_Glucosidase==0)]<- NA
  lab_sql$DNAConc[which(lab_sql$DNAConc==0)]<- NA
  lab_sql$SoilMoisture[which(lab_sql$SoilMoisture< 0)]<- NA
  #lab_sql$SoilDW[which(lab_sql$DW< 0)]<- NA
 # lab_sql$SoilWW[which(lab_sql$WW< 0)]<- NA

  lab_long<- reshape2::melt(lab_sql)
  lab_longf<- lab_long %>% dplyr::filter( !is.na(value))

  return(lab_longf)

}


#' A phylosql Function
#'
#' formats lab to for sql upload
#' @param labdata data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @export

fix_soil_moisture<-
  function(labdata){

   soilMoisture <-  labdata$SoilDW/labdata$SoilWW
   labdata$SoilMoisture[which(!is.na(soilMoisture))]<- soilMoisture[which(!is.na(soilMoisture))]
   labdata
}


#' A phylosql Function
#'
#' formats cms data to long for upload. Strict mode.
#' @param data data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import tidyr
#' @export
#'


cms_as_longCms<- function(data){
  data_long <- tidyr::gather(data,Factor, Level, PropertyName:TargetCropTypeName, factor_key=TRUE)
  return(data_long)
}

#' A phylosql Function
#'
#' formats misc data to cms data
#' @param data data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import tidyr
#' @export
#'

misc_as_longCms<- function(data){
  stopifnot(colnames(data)[1]=="metagenNumber")
  cols<- colnames(data)
  data_long <- gather(cmsdata,Factor, Level, cols[2]:cols[length(cols)], factor_key=TRUE)
  return(data_long)
}




#' A phylosql Function
#'
#' formats misc data to cms data
#' @param labdata data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import tidyr
#' @export
#'

fix_soil_moisture<-
  function (labdata)
  {
    not_nas<- which(!is.na(labdata$SoilDW))

    soilMoisture<- 1-(na.omit(labdata$SoilDW)/ na.omit(labdata$SoilWW))

    labdata$SoilMoisture[not_nas] <- soilMoisture
    labdata$SoilMoisture[(labdata$SoilMoisture>1|labdata$SoilMoisture< 0)]<- NA

    labdata
  }



#' A phylosql Function
#'
#' formats misc data to cms data
#' @param labdata data to format
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import tidyr
#' @export
#'
clean_dna_yields<-
  function(dnas){


    above<- which(grepl(">",dnas$DNA))
    below<- which(grepl("<",dnas$DNA))
    dnas$DNA<- gsub('>',"",dnas$DNA)
    dnas$DNA<- gsub('<',"",dnas$DNA)

    dnas$range<- "inrange"

    if(length(above)>0){

      dnas$range[above]<-"above"

    }

    if(length(below)>0){

      dnas$range[below]<-"below"

    }

    dnas
  }


#' A phylosql Function
#'
#' formats misc data to cms data
#' @param dna data to format
#' @param dna_elution_volume
#' @keywords
#' @import dplyr
#' @import RMariaDB
#' @import tidyr
#' @export
#'
transform_dna_yield_to_kg<-
  function(dna,dna_elution_volume=2000){

   dna * dna_elution_volume * 3.33 * 1.55 * 10* 1e-6

  }


