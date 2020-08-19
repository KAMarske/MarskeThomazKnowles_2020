
## Estimates Peter & Slatkin's (2013) directionality index using the X-Origin package by Qixin He
## The code below require the following source code for execution: 

# https://github.com/KnowlesLab/X-ORGIN/blob/master/precheck/re_functions_resistance.r

#Acknowledgement: A huge thanks to Qixin He for help getting started!



# Each geographic cluster (identified via STRUCTURE) is run separately.

#snapp files based on stru file used for regional input in Procrutes, with 25% limit on missing data

#To generate snapp file from structure file:
#Open in excel, text to columns, delete second column (population number), save as csv
#Open in Notepad++, save by overwriting another snapp file to get correct extension



########################### 

rm(list = ls())

library(geosphere)
library(rworldmap)



#########  North   ########################################################

source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")
setwd("~/Agyrtodes/psi/biogeo2/north")

snp_file <-c("north_25.SNAPP")
coords_file <-c("north_locs.txt") ##must be comma delimited
psi_name<-c("Agyrt_north_PsiOut.txt")


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

#bbox<-get_sample_bbox(pop_coords)
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#generates the statistic
north_origin <- find_origin(pop_coords, all_psi, region=c("AbelTas", "MtArthur", "Karamea", "Marlborough"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("AbelTas", "MtArthur", "Karamea", "Marlborough"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)

north_origin




######################################################################################################################
#########  East   ########################################################

rm(list = ls())

source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")

setwd("~/Agyrtodes/psi/biogeo2/east")

snp_file <-c("east_25.SNAPP")
coords_file <-c("east_locs.txt") ##must be comma delimited
psi_name<-c("Agyrt_east_PsiOut.txt")


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

#bbox<-get_sample_bbox(pop_coords)
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#generates the statistic
east_origin <- find_origin(pop_coords, all_psi, region=c("Rakaia","Southcant","Kaikoura","MtFyffe","Midcant"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Rakaia","Southcant","Kaikoura","MtFyffe","Midcant"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)

east_origin




######################################################################################################################
######################################################################################################################
#########  South   ########################################################

rm(list = ls())

source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")

setwd("~/Agyrtodes/psi/biogeo2/south")

snp_file <-c("south_25.SNAPP")
coords_file <-c("south_locs.txt") ##must be comma delimited
psi_name<-c("Agyrt_south_PsiOut.txt")


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

#bbox<-get_sample_bbox(pop_coords)
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#generates the statistic
south_origin <- find_origin(pop_coords, all_psi, region=c("Catlins","Southland","StewartIs","Otago","WesternEdge","Haast","Milford"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Catlins","Southland","StewartIs","Otago","WesternEdge","Haast","Milford"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)

south_origin



######################################################################################################################

######################################################################################################################
#########  West   ########################################################

rm(list = ls())

source("~/XOrigin/X-ORGIN-1.0/X-ORGIN-1.0/precheck/re_functions_resistance.r")

setwd("~/Agyrtodes/psi/biogeo2/west")

snp_file <-c("west_25.SNAPP")
coords_file <-c("west_locs.txt") ##must be comma delimited
psi_name<-c("Agyrt_west_PsiOut.txt")


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

#bbox<-get_sample_bbox(pop_coords)
SIbbox<-rbind(c(167,174.5), c(-49,-40))

#generates the statistic
west_origin <- find_origin(pop_coords, all_psi, region=c("Westland","Buller","Pororari","ArthursPass","Lewis","Haast","Governor","SouthWest"), countries=NULL, xlen=50, ylen=50, doPlot=T, doPdf=T,f.dist="haversine", rescale = 1000,bbox = SIbbox, IBRmap = NA)#fdist could also be euclidean

#makes the plot
run_region(region=c("Westland","Buller","Pororari","ArthursPass","Lewis","Haast","Governor","SouthWest"),loc_file_id=out_file_id, xlen=50, ylen=50, f.dist="haversine", rescale=1000, bbox=SIbbox, IBRmap = NA)

west_origin



######################################################################################################################




