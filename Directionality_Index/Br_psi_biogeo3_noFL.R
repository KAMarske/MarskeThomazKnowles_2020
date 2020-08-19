## Estimates Peter & Slatkin's (2013) directionality index using the X-Origin package by Qixin He

## The code below require the following source code for execution: 

# https://github.com/KnowlesLab/X-ORGIN/blob/master/precheck/re_functions_resistance.r

#Acknowledgement: A huge thanks to Qixin He for help getting started!

# Each geographic cluster (identified via STRUCTURE) is run separately.
#snapp files based on stru file used for regional input in Procrutes, with 25% limit on missing data

#To generate snapp file from structure file:
#Open in excel, text to columns, delete second column (population number), save as csv
#Open in Notepad++, save by overwriting another snapp file to get correct extension


rm(list = ls())

library(geosphere)
library(rworldmap)


###############################################################################################
######################  Western SI, excluding Flora and Kaituna  ######################
rm(list = ls())
source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Brachynopus/psi/biogeo3_noFlora/west_noFL")

snp_file <-c("westNoFL.SNAPP")
coords_file <-c("westNoFL_locs.txt") ##must be comma delimited
psi_name<-c("Brachy_westNoFL_biogeo3.txt")

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid
subSample <- 0 #Value of number of individual subsample down to, if no subsample is required, put 0
outgroup_columns <- NULL 
nsnp <- NULL
out_file_id = sprintf( "out_%s", basename( snp_file ) )


if(!is.data.frame(data))data <- c()
read_data_snapp(snp_file=snp_file,
                coords_file=coords_file, nsnp=nsnp,
                outgroup_columns=outgroup_columns
)


pops <- make_pops_snapp( coords,n=subSample)  
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data,pop_ss))


##Calculate pairwise Psi


all_psi <- get_all_psi(pop_data, pop_ss ,n=subSample)
psiVec<-c()
k = 1
for (i in 1: (nrow(all_psi)-1)) {
  for (j in (i+1): nrow(all_psi)) {
    psiVec<-cbind(psiVec, all_psi[i,j])
    colnames(psiVec)[k]<-paste("psi",i,"_",j,sep = "")
    k = k+1
  }
}

write.table(psiVec,psi_name,col.names=T,row.names=F,sep="\t",quote=F)



#user-defined bounding box
#fullbbox<-rbind(c(167,176.5), c(-48,-40))
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#bbox<-get_sample_bbox(pop_coords)

#generates the statistic
westSI_noFL <- find_origin(pop_coords, all_psi, region=c("Arthur","Buller", "Lewis","Karamea","LowerHaast","UpperHaast","Milford","Southland","Pororari","StewartIs","WestCoast","WesternEdge"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Arthur","Buller", "Lewis","Karamea","LowerHaast","UpperHaast","Milford","Southland","Pororari","StewartIs","WestCoast","WesternEdge"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)


westSI_noFL




###############################################################################################
######################  Western SI, including Flora and Kaituna, small localities  ######################
rm(list = ls())
source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Brachynopus/psi/biogeo3_noFlora/west_withFL")

snp_file <-c("west.SNAPP")
coords_file <-c("west_locs.txt") ##must be comma delimited
psi_name<-c("Brachy_biogeo3.txt")

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid
subSample <- 0 #Value of number of individual subsample down to, if no subsample is required, put 0
outgroup_columns <- NULL 
nsnp <- NULL
out_file_id = sprintf( "out_%s", basename( snp_file ) )


if(!is.data.frame(data))data <- c()
read_data_snapp(snp_file=snp_file,
                coords_file=coords_file, nsnp=nsnp,
                outgroup_columns=outgroup_columns
)


pops <- make_pops_snapp( coords,n=subSample)  
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data,pop_ss))


##Calculate pairwise Psi


all_psi <- get_all_psi(pop_data, pop_ss ,n=subSample)
psiVec<-c()
k = 1
for (i in 1: (nrow(all_psi)-1)) {
  for (j in (i+1): nrow(all_psi)) {
    psiVec<-cbind(psiVec, all_psi[i,j])
    colnames(psiVec)[k]<-paste("psi",i,"_",j,sep = "")
    k = k+1
  }
}

write.table(psiVec,psi_name,col.names=T,row.names=F,sep="\t",quote=F)



#user-defined bounding box
#fullbbox<-rbind(c(167,176.5), c(-48,-40))
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#bbox<-get_sample_bbox(pop_coords)

#generates the statistic
westSI_withFL <- find_origin(pop_coords, all_psi, region=c("Arthur","Buller", "Lewis","Karamea","LowerHaast","UpperHaast","Milford","Southland","Pororari","StewartIs","WestCoast","WesternEdge","Nelson"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Arthur","Buller", "Lewis","Karamea","LowerHaast","UpperHaast","Milford","Southland","Pororari","StewartIs","WestCoast","WesternEdge","Nelson"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)


westSI_withFL


###############################################################################################
######################  NW SI, excluding Flora and Kaituna  ######################
rm(list = ls())
source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Brachynopus/psi/biogeo3_noFlora/NW_noFL")

snp_file <-c("NW_noFlora.SNAPP")
coords_file <-c("NW_NoFL_locs.txt") ##must be comma delimited
psi_name<-c("NW_NoFL_biogeo3.txt")

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid
subSample <- 0 #Value of number of individual subsample down to, if no subsample is required, put 0
outgroup_columns <- NULL 
nsnp <- NULL
out_file_id = sprintf( "out_%s", basename( snp_file ) )


if(!is.data.frame(data))data <- c()
read_data_snapp(snp_file=snp_file,
                coords_file=coords_file, nsnp=nsnp,
                outgroup_columns=outgroup_columns
)


pops <- make_pops_snapp( coords,n=subSample)  
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data,pop_ss))


##Calculate pairwise Psi


all_psi <- get_all_psi(pop_data, pop_ss ,n=subSample)
psiVec<-c()
k = 1
for (i in 1: (nrow(all_psi)-1)) {
  for (j in (i+1): nrow(all_psi)) {
    psiVec<-cbind(psiVec, all_psi[i,j])
    colnames(psiVec)[k]<-paste("psi",i,"_",j,sep = "")
    k = k+1
  }
}

write.table(psiVec,psi_name,col.names=T,row.names=F,sep="\t",quote=F)



#user-defined bounding box
#fullbbox<-rbind(c(167,176.5), c(-48,-40))
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#bbox<-get_sample_bbox(pop_coords)

#generates the statistic
NW_noFL <- find_origin(pop_coords, all_psi, region=c("Arthur","Buller", "Lewis","Karamea","Pororari"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Arthur","Buller", "Lewis","Karamea","Pororari"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)


NW_noFL


###############################################################################################
######################  SW SI, excluding Flora and Kaituna  ######################
rm(list = ls())
source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Brachynopus/psi/biogeo3_noFlora/SW_noFL")

snp_file <-c("SW_noFlora.SNAPP")
coords_file <-c("SW_NoFL_locs.txt") ##must be comma delimited
psi_name<-c("SW_NoFL_biogeo3.txt")

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid
subSample <- 0 #Value of number of individual subsample down to, if no subsample is required, put 0
outgroup_columns <- NULL 
nsnp <- NULL
out_file_id = sprintf( "out_%s", basename( snp_file ) )


if(!is.data.frame(data))data <- c()
read_data_snapp(snp_file=snp_file,
                coords_file=coords_file, nsnp=nsnp,
                outgroup_columns=outgroup_columns
)


pops <- make_pops_snapp( coords,n=subSample)  
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data,pop_ss))


##Calculate pairwise Psi


all_psi <- get_all_psi(pop_data, pop_ss ,n=subSample)
psiVec<-c()
k = 1
for (i in 1: (nrow(all_psi)-1)) {
  for (j in (i+1): nrow(all_psi)) {
    psiVec<-cbind(psiVec, all_psi[i,j])
    colnames(psiVec)[k]<-paste("psi",i,"_",j,sep = "")
    k = k+1
  }
}

write.table(psiVec,psi_name,col.names=T,row.names=F,sep="\t",quote=F)



#user-defined bounding box
#fullbbox<-rbind(c(167,176.5), c(-48,-40))
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#bbox<-get_sample_bbox(pop_coords)

#generates the statistic
SW_noFL <- find_origin(pop_coords, all_psi, region=c("LowerHaast","UpperHaast","Milford","Southland","StewartIs","WestCoast","WesternEdge"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("LowerHaast","UpperHaast","Milford","Southland","StewartIs","WestCoast","WesternEdge"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)


SW_noFL

##########################################################################################
##########   North, including NI      ###########################################################

rm(list = ls())
source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Brachynopus/psi/biogeo2/north")

snp_file <-c("north.snapp")
coords_file <-c("north_locs.txt") ##must be comma delimited
psi_name<-c("Brachy_north_Out.txt")

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid
subSample <- 0 #Value of number of individual subsample down to, if no subsample is required, put 0
outgroup_columns <- NULL 
nsnp <- NULL
out_file_id = sprintf( "out_%s", basename( snp_file ) )


if(!is.data.frame(data))data <- c()
read_data_snapp(snp_file=snp_file,
                coords_file=coords_file, nsnp=nsnp,
                outgroup_columns=outgroup_columns
)

#the number this spits out is twice the number of loci in the file, maybe is the number of snp? when I use the manual snp file after removing missing data
#with the python generated infiles it is close to 2x but not exactly.

pops <- make_pops_snapp( coords,n=subSample)  
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data,pop_ss))


##Calculate pairwise Psi


all_psi <- get_all_psi(pop_data, pop_ss ,n=subSample)
psiVec<-c()
k = 1
for (i in 1: (nrow(all_psi)-1)) {
  for (j in (i+1): nrow(all_psi)) {
    psiVec<-cbind(psiVec, all_psi[i,j])
    colnames(psiVec)[k]<-paste("psi",i,"_",j,sep = "")
    k = k+1
  }
}

write.table(psiVec,psi_name,col.names=T,row.names=F,sep="\t",quote=F)



#user-defined bounding box
fullbbox<-rbind(c(167,176.5), c(-48,-40))

#bbox<-get_sample_bbox(pop_coords)
#generates the statistic
north <- find_origin(pop_coords, all_psi, region=c("Nelson","Inland","Marlborough","NI"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = fullbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Nelson","Inland","Marlborough","NI"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=fullbbox, IBRmap = NA)

north


###############################################################################################
######################  Eastern SI  ######################
rm(list = ls())
source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Brachynopus/psi/biogeo2/east")

snp_file <-c("east.snapp")
coords_file <-c("east_locs.txt") ##must be comma delimited
psi_name<-c("Brachy_eastOut.txt")

ploidy <- 2  #set ploidy of individuals. 1=haploid, 2 =diploid
subSample <- 0 #Value of number of individual subsample down to, if no subsample is required, put 0
outgroup_columns <- NULL 
nsnp <- NULL
out_file_id = sprintf( "out_%s", basename( snp_file ) )


if(!is.data.frame(data))data <- c()
read_data_snapp(snp_file=snp_file,
                coords_file=coords_file, nsnp=nsnp,
                outgroup_columns=outgroup_columns
)


pops <- make_pops_snapp( coords,n=subSample)  ##are there too few localities with enough samples?
pop_data <- make_pop_data_from_pops( pops, data )
pop_coords <- make_pop_coords_from_pops( pops, coords)
pop_ss <- make_pop_ss_from_pops( pops, data, ploidy=ploidy)
pop_coords <- cbind( pop_coords, hets=get_heterozygosity(pop_data,pop_ss))


##Calculate pairwise Psi


all_psi <- get_all_psi(pop_data, pop_ss ,n=subSample)
psiVec<-c()
k = 1
for (i in 1: (nrow(all_psi)-1)) {
  for (j in (i+1): nrow(all_psi)) {
    psiVec<-cbind(psiVec, all_psi[i,j])
    colnames(psiVec)[k]<-paste("psi",i,"_",j,sep = "")
    k = k+1
  }
}

write.table(psiVec,psi_name,col.names=T,row.names=F,sep="\t",quote=F)



#user-defined bounding box
#fullbbox<-rbind(c(167,176.5), c(-48,-40))

#bbox<-get_sample_bbox(pop_coords)
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#generates the statistic
eastSI <- find_origin(pop_coords, all_psi, region=c("Rakaia","Haast","South","Midcant","SouthCant","Kaikoura"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Rakaia","Haast","South","Midcant","SouthCant","Kaikoura"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)


eastSI



