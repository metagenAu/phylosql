
library(RMariaDB)
library(dplyr)


con <- dbConnect(RMariaDB::MariaDB(), user='mtgn_one', password="SQLmtgn_pword1", dbname='mtgn',
                 host='localhost')

library(phylosql)

cms<- read.csv("C:/Users/Chris/Documents/cms.csv",stringsAsFactors = FALSE)
ms<- read_masterSheet()

labdata<- as_labData(ms)
cmsdata<- cms_as_cmsData(cms)
cmsdata.ms<- ms_as_cmsData(ms)

upload_cms_data(cmsdata,con=con)
upload_lab_data(labdata,con=con)
upload_cms_data(cmsdata.ms,con=con)


tax18s<- fetch_taxonomy(con=con,phylo=TRUE)
si<- fetch_sampleInfo(con=con)
sv18s<- fetch_asv_table(con=con,phylo=TRUE)

ps18s<- phyloseq(otu_table(sv18s), tax_table(tax18s),sample_data(si))



tax16s<- fetch_taxonomy(con=con,phylo=TRUE)
sv16s<- fetch_asv_table(con=con,phylo=TRUE)

ps16s<- phyloseq(otu_table(sv16s), tax_table(tax16s),sample_data(si))
