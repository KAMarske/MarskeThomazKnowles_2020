

############Procrustes for 98 Agyrtodes labralis individuals from 67 localities, maximum 25% missing data

##Acknowledgements: A huge thanks to Luciana Resende, Melisa Olave, and Mariah Kenney for their help in developing this code!




#############################################################
##############Geographic Procrustes##########################

rm(list = ls())

library(adegenet)
library(ade4)
library(vegan)
library(ape)
library(dismo)
library(raster)
library(rgdal)
library(geosphere)


setwd("~/Agyrtodes/Procrustes final/all")


###Import structure data and do PCA


#import structure data as genind object
getstr<-read.structure('98inds_25miss_4422loci.stru', n.ind=98, n.loc=4422, onerowperind=T,col.lab=1,col.pop=2,row.marknames=0,  ask=F, NA.char=0)


#this makes matrix
getstr_scaled <- scaleGen(getstr, NA.method = c("mean") )
genom_pca <- dudi.pca(getstr_scaled, center = T, scale = T, scannf=F, nf=2)

#s.class(genom_pca$li, getstr$pop, cpoint = 1)##labelled plot
#plot(genom_pca$li) ##basic plot with just dots




###Import locality data for all 98 individuals
#must be in same order as structure data, visually check that StacksPopNo are in sequential order
#colnames IID, StacksPop, StacksPopNo, x, y

locs<-read.csv('98inds.csv')
View(locs)
#head(locs)

# save a population vector to use for mapping loop
pop <- as.integer(locs$StacksPopNo)

locs$IID<-NULL
locs$StacksPop<-NULL #Brachy input
locs$IndNo<-NULL #Agyrt input
locs$StacksPopNo<-NULL


###Protest and Procrustes

geo_prot <- protest(X=locs,Y=genom_pca$li,scale=T,symmetric=T,permutations=10000)
geo_prot
summary(geo_prot)
resid(geo_prot)

geo_proc <- procrustes(X=locs,Y=genom_pca$li,scale=T,symmetric=F)
summary(geo_proc)
resid(geo_proc)


plot(geo_proc, main=NULL, ar.col='dark gray',xlab='',ylab='') #basic plot in pca space
dev.print(pdf, 'pca_proc_98inds.pdf') #will not overwrite, must delete if rerunning

plot(geo_proc, main=NULL, ar.col='black',xlab='',ylab='', kind=2) #impulse plot of residuals: x is individual, y is resid
dev.print(pdf, 'proc_resids_by_ind.pdf')

plot(geo_prot, main=NULL, ar.col='dark gray',xlab='',ylab='')
dev.print(pdf, 'pca_protest_98inds.pdf')

geo_resids <-as.data.frame(residuals(geo_proc))


geo_offset_pts <-fitted(geo_proc)##test plot offset points
as.data.frame(geo_offset_pts)



#####Mapping Results

#Read in locality info and map color code (HEX codes) for the 67 localities
#table columns are popn, x, y, pca.color, loc.color


colorTable <- read.csv('98pops.csv')
#View(colorTable) 


#The following plots the results in a graph with long and lat as axes

plot.new();
par(oma=c(1,1,1,1));
par(mar=c(3,3,1,0));
par(pty='s');
plot.window(xlim=c(167, 176), ylim=c(-47, -38));
axis(1, at=seq(167, 176, by=2),  cex.axis=1);
axis(2, at=seq(-47, -38, by=2),  cex.axis=1,las=3);
#axis labels
mtext(side=1, text='Longitude',line=2.2, cex=1.15)
mtext(side=2, text='Latitude', line=2.5, cex=1.15)

##add arrows on top of points
arrows<-arrows(geo_proc$Y[,1]+geo_proc$translation[1,1],geo_proc$Y[,2]+geo_proc$translation[1,2],geo_proc$X[,1]+geo_proc$translation[1,1],geo_proc$X[,2]+geo_proc$translation[1,2],length=0.05,col='dark grey', cex=2)


ind <- 1:length(pop);
table <- data.frame(pop,ind);
for(i in unique(pop)){
  indVec <- table[table$pop == i,2];
  points(geo_proc$Y[indVec,1]+geo_proc$translation[1,1],geo_proc$Y[indVec,2]+geo_proc$translation[1,2],cex=0.8,pch=21,bg=paste("#", colorTable[colorTable$pop == i,4], sep=""),col='black')
  points(geo_proc$X[indVec[1],1]+geo_proc$translation[1,1],geo_proc$X[indVec[1],2]+geo_proc$translation[1,2],cex=1.0,pch=24,bg=paste("#",colorTable[colorTable$pop == i,5], sep=""),col='black')
}
dev.print(pdf, 'Agyrt_latlong_98inds.pdf')


dim1 <- geo_proc$translation[,1]
dim2 <- geo_proc$translation[,2]

######generates lat-long coords for pca locations based on geographic centroid (given by proc$translation)
#see above, this may yield the same as the fitted(proc) command
xout <- geo_proc$Y
for(i in 1:length(geo_proc$Y[,1])) {
  xout[i,1] <- geo_proc$Y[i,1]+dim1
  xout[i,2] <- geo_proc$Y[i,2]+dim2
}

for(i in 1:length(xout[,1])) {
  cat(xout[i,],file='98_offsetcoords.txt',sep='\t',append=T)
  cat('\n',file='98_offsetcoords.txt',append=T)
}


#The following plots onto a map of New Zealand

nzcoast <- raster(x='nz_coastline.asc')

#full outline
#plot(nzcoast, col=gray(1:10/10, alpha=0.9),axes=F, xlim=c(165, 185), ylim=c(-50, -30),bty='n',legend=F, box=F)
#zoom to South Island
plot(nzcoast, col=gray(1:10/10, alpha=0.9),axes=F, xlim=c(165, 180), ylim=c(-50, -40),bty='n',legend=F, box=F)
#zoom to northern SI
#plot(nzcoast, col=gray(1:10/10, alpha=0.9),axes=F, xlim=c(170, 175), ylim=c(-43, -40),bty='n',legend=F, box=F)


##add arrows underneath points
arrows(xout[,1],xout[,2],locs[,1],locs[,2],length=0,lwd=1,col='black')

##This requires the previous mapping loop (above setting the dimensions) to run, gives exactly the same result as above
for(i in unique(pop)){
  indVec <- table[table$pop == i,2];
  #points(locs[indVec[1],1],locs[indVec[1],2],cex=1.0,pch=17,col=paste("#",colorTable[colorTable$popn == i,5], sep=""))
  points(locs[indVec[1],1],locs[indVec[1],2],cex=1.0,pch=24,bg=paste("#",colorTable[colorTable$popn == i,5], sep=""),col='black')
  points(xout[indVec,1],xout[indVec,2], cex=0.8, pch=21, bg=paste("#", colorTable[colorTable$popn == i,4], sep=""),col='black')
}

dev.print(pdf, 'Agyrt_map_98inds.pdf') #will not overwrite, must delete if rerunning



#####Distribution of Residuals


actual_distances <- distm(locs[,c('x','y')], fun=distVincentyEllipsoid) #result is in meters
maxdist <- max(actual_distances)/1000

geo_resid_distances<- distGeo(locs[,c('x','y')], geo_offset_pts[,c(1,2)], a=6378137, f=1/298.257223563) #result is in meters
#offset_pts are taken directly from the procrustes, not from the manual translation of locations

#View(resid_distances)
#max(resid_distances)
hist(geo_resid_distances/1000, breaks=50, main=paste("Agyrt resids (km)"))
#hist(geo_resid_distances/1000, breaks=10)#for northern SI
hist(geo_resid_distances/1000, breaks=20,main=paste("Agyrt resids (km)"))


hist(geo_resid_distances/1000, breaks=20,xlim=c(0,350), ylim=c(0,20), main=paste("Agyrt resids (km)"))
dev.print(pdf, 'Agyrt_geodist_resid_hist_98inds.pdf')



###################################################################################################################
##############Partial Procrustes###################################################################################


setwd("~/Worldclim data/Current/working subset")

files <- list.files(pattern='asc', full.names=TRUE )
files

predictors<- stack(files)
predictors
names(predictors)

rasterRescale<-function(r){
  ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
}

predictorsRescale<-rasterRescale(predictors) #rescales so all values between 0-1

setwd("~/Agyrtodes/Procrustes final")






######Generate the environmental PCA for global analysis (PCA on climate rasters, then extract data for point localities)

glenv_pca <- prcomp(na.omit(predictorsRescale[]))

summary(glenv_pca)
glenv_pca$rotation  #eigenvectors
glenv_pca$sdev^2  #eigenvalues


pcaRast <- predict(predictorsRescale, glenv_pca, index=1:3)

pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*400
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*400
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*400


plotRGB(pcaRast, r=1, g=2, b=3)
dev.print(pdf, 'curr_environment.pdf')
writeRaster(pcaRast, "curr_environment_pca.tif", format="GTiff", overwrite=T) #Exports raster of PCA


#Extract pca data to localities
as.data.frame(locs)
glenv_pca_res<-extract(pcaRast, locs, method='simple')
str(glenv_pca_res)



##Use procrustes to look at how localities are clustered in environmental space--FOR VISUALIZATION ONLY
##Using x axis information from regular procrustes because this PCA has three axes of variation instead of 2



glenvmap_prot <- protest(X=locs,Y=glenv_pca_res,scale=T,symmetric=T,permutations=10000)
#glenvmap_prot
glenvmap_proc <- procrustes(X=locs,Y=glenv_pca_res,scale=T,symmetric=F)
#summary(glenvmap_proc)
plot(glenvmap_proc, main=NULL, ar.col='dark gray')#basic plot in pca space
dev.print(pdf, 'glenv_pca_98inds.pdf') #will not overwrite, must delete if rerunning


plot.new();
par(oma=c(1,1,1,1));
par(mar=c(3,3,1,0));
par(pty='s');
plot.window(xlim=c(167, 176), ylim=c(-47, -38));
axis(1, at=seq(167, 176, by=2),  cex.axis=1);
axis(2, at=seq(-47, -38, by=2),  cex.axis=1,las=3);
#axis labels
mtext(side=1, text='Longitude',line=2.2, cex=1.15)
mtext(side=2, text='Latitude', line=2.5, cex=1.15)

arrows<-arrows(glenvmap_proc$Y[,1]+glenvmap_proc$translation[1,1],glenvmap_proc$Y[,2]+glenvmap_proc$translation[1,2],geo_proc$X[,1]+geo_proc$translation[1,1],geo_proc$X[,2]+geo_proc$translation[1,2],length=0.05,col='dark gray', cex=2)


ind <- 1:length(pop);
table <- data.frame(pop,ind);
for(i in unique(pop)){
  indVec <- table[table$pop == i,2];
  points(glenvmap_proc$Y[indVec,1]+glenvmap_proc$translation[1,1],glenvmap_proc$Y[indVec,2]+glenvmap_proc$translation[1,2],cex=0.8,pch=21,bg=paste("#", colorTable[colorTable$pop == i,4], sep=""),col='black')
  points(geo_proc$X[indVec[1],1]+geo_proc$translation[1,1],geo_proc$X[indVec[1],2]+geo_proc$translation[1,2],cex=1.0,pch=24,bg=paste("#",colorTable[colorTable$pop == i,5], sep=""),col='black')
}

#dev.print(pdf, 'Gl_env_latlong_dotsonly.pdf')
dev.print(pdf, 'Gl_env_latlong_98inds.pdf')

######Generate the environmental PCA for local analysis (extract environmental variables from raster stacks and then generate PCA for just those locs)
localenv <- extract(predictorsRescale, locs)
pca_localenv <- prcomp(localenv)

pca_localenv
summary(pca_localenv)
pca_localenv$rotation  #eigenvectors
pca_localenv$sdev^2  #eigenvalues



##Partial protest and procrustes on global environmental data

gl_part_prot <- protest(X=glenv_pca_res, Y=geo_resids,scale=T,symmetric=T,permutations=10000)
gl_part_prot
summary(gl_part_prot)
resid(gl_part_prot)

gl_part_proc <- procrustes(X=glenv_pca_res, Y=geo_resids,scale=T,symmetric=F)
summary(gl_part_proc)
resid(gl_part_proc)

plot(gl_part_prot, cex=1.0, main=NULL, ar.col='dark grey',xlab='',ylab='')
dev.print(pdf, 'pca_gl_part_prot.pdf')
plot(gl_part_proc, cex=1.0, main=NULL, ar.col='dark grey',xlab='',ylab='')#basic plot in pca space
dev.print(pdf, 'pca_gl_part_proc.pdf')
#plot(gl_part_proc, cex=1.0, main=NULL, ar.col='dark grey',xlab='',ylab='', kind=2)


#partial protest and procrustes on local environmental data

local_part_prot <- protest(X=pca_localenv, Y=geo_resids,scale=T,symmetric=T,permutations=10000)
local_part_prot
summary(local_part_prot)
resid(local_part_prot)

local_part_proc <- procrustes(X=pca_localenv, Y=geo_resids,scale=T,symmetric=F)
summary(local_part_proc)
resid(local_part_proc)

plot(local_part_prot, cex=1.0, main=NULL, ar.col='dark grey',xlab='',ylab='')#basic plot in pca space
dev.print(pdf, 'pca_local_part_prot.pdf')
plot(local_part_proc, cex=1.0, main=NULL, ar.col='dark grey',xlab='',ylab='')#basic plot in pca space
dev.print(pdf, 'pca_local_part_proc.pdf')
plot(local_part_proc, main=NULL, ar.col='black',xlab='',ylab='', kind=2)


##########################################################################################
##############Environmental Procrustes####################################################

#These use the genomic and environmental PCAs generated above


#####Full Environmental Space

glenv_prot <- protest(X=glenv_pca_res, Y=genom_pca$li,scale=T,symmetric=T,permutations=10000)
glenv_prot
summary(glenv_prot)
resid(glenv_prot)

glenv_proc <- procrustes(X=glenv_pca_res, Y=genom_pca$li,scale=T,symmetric=F)
summary(glenv_proc)

#plot(glenv_proc, main=NULL, ar.col='gray',xlab='',ylab='')#basic plot in pca space
#plot(glenv_proc, main=NULL, ar.col='gray',xlab='',ylab='',kind=2)
glenv_resids <-as.data.frame(residuals(glenv_proc))



#####Local Environmental Space

localenv_prot <- protest(X=pca_localenv, Y=genom_pca$li,scale=T,symmetric=T,permutations=10000)
localenv_prot
summary(localenv_prot)
resid(localenv_prot)

localenv_proc <- procrustes(X=pca_localenv, Y=genom_pca$li,scale=T,symmetric=F)
summary(localenv_proc)
resid(localenv_prot)

plot(localenv_prot, cex=1.0, main=NULL, ar.col='gray',xlab='',ylab='')#basic plot in pca space
dev.print(pdf, 'localenv_prot.pdf')
plot(localenv_proc, cex=1.0, main=NULL, ar.col='gray',xlab='',ylab='')#basic plot in pca space
dev.print(pdf, 'localenv_proc.pdf')
#plot(localenv_proc, cex=1.0, main=NULL, ar.col='gray',xlab='',ylab='', kind=2)
localenv_resids <-as.data.frame(residuals(localenv_proc))


#####Distribution of Residuals from environmental procrustes


#plot(regular$regular_resid, stable$stability_resid)  ###basically no pattern

plot(geo_resids[,1], localenv_resids[,1], main='geo_resids by localenv_resids', xlab='geo resid',ylab='local enviro resid')
dev.print(pdf, 'resids_localenv.pdf')
plot(geo_resids[,1], glenv_resids[,1], main='geo_resids by globalenv_resids',xlab='geo resid',ylab='global enviro resid')
dev.print(pdf, 'resids_globalenv.pdf')



