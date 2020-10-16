
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

  rownames(sample_data)<- samples
  colnames(sample_data)<- vars
  sample_data<- data.frame(sample_data)
  sample_data$MetagenNumber<- samples

  return(sample_data)
}


fetch_sampleInfo<-
  function(...){
    si<- as_tibble(tbl(con,"labdata"))
    sample_info<- create_sampleInfo_table(si_long=si)
    cms<- as_tibble(tbl(con,"cmsdata"))
    sampleInfo <- merge(sample_info,cms,"MetagenNumber")
    class(sampleInfo)<- "sampledata"
    return(sampleInfo)
  }


si<- as_tibble(tbl(con,"labdata"))

sampleInfo<- create_sampleInfo_table(si)



create_asv_table<- function(con=NULL ){

  asv_long<- as_tibble(tbl(con))

  asvs<- unique(asv_long$SV)
  samples<- unique(asv_long$MetagenNumber)
  rows<- length(samples)
  cols<- length(asvs)
  asv_table<- matrix(0L,ncol=cols, nrow=rows)

  sampleList<- split(asv_long,asv_long$MetagenNumber)

  for(i in seq_along(sampleList)){
    asv_table[i,match(sampleList[[i]]$SV, asvs)  ]<- sampleList[[i]]$Abundance
  }

  rownames(asv_table)<- samples
  colnames(asv_table)<- asvs
  class(asv_table)<- "abundance"
  return(asv_table)
}



fetch_taxonomy<- function

class(taxonomy)<- "taxonomy"
