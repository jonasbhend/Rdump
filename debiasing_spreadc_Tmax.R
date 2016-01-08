# Load Libraries
#library(ncvar)
library(ncdf)
#
#
# directories, specs
#

# initialise with variable name


dirhc <- '/store/msclim/Prediction/Seasonal/ecmwf-system4/daily/global/November_init/t2m/'
hcfilestem<-"01_55_global_daily.nc"
obsdir <-'/store/msclim/Prediction/ERAINTERIM/tas/daily/erainterim_tas_'
obsfilestem<-".nc"
initmonth<-11

dir.out <- '/store/msclim/bhendj/EUPORIAS/Tmean_fcst/'
if (!file.exists(dir.out)) dir.create(dir.out, recurs=TRUE)
#fcspecs
fcdays<-215
tvec<-1:fcdays
fcmembers<-51
fcyears<-32
inityr<-1981
yearstr<-as.character(seq(inityr,inityr+fcyears-1))
parfc<-"mean2t24"
parobs<-"t2m"
parout<-"mn2t24"

file.names<-rep(NA,fcyears)
outfile.names<-file.names
for (i in (1:fcyears)){
  file.names[i]<-paste(dirhc,yearstr[i],as.character(initmonth),hcfilestem,sep="")
  outfile.names[i]<-paste(dir.out,yearstr[i],as.character(initmonth),hcfilestem,sep="")
}

## replace files with the ones in Data_from_APN
## file.names <- system('\\ls /store/msclim/Prediction/Seasonal/ecmwf-system4/daily/global/Data_from_APN/*/seasfc_mn*110100.nc', TRUE)[1:fcyears]

#------------------------------------------
#read model data and make hindcast mean
#------------------------------------------
#-----------------------------------------------------------------------------
#this routine makes it possible to read in only parts of the file. In this case it loops over the latitudes
#-----------------------------------------------------------------------------

connections <- vector("list",length(file.names))
obsfilecons<-connections
# f?Á?r jedes File eine Connection erstellen
for (i in 1:length(file.names)) {
  file.in <- file.names[i]
  # yearstr[i]<-substring(file.names[i],1,4)
  print(file.in)
  connections[[i]] <- open.ncdf(file.in,readunlim=FALSE)
  obs.file<-paste(obsdir,yearstr[i],obsfilestem,sep="")
  obsfilecons[[i]]<-open.ncdf(obs.file,readunlim=FALSE)
}

mean.tmin <-  array(NA,dim=c(180,91,fcdays))
sd.tmin<-mean.tmin
mean.obs<-mean.tmin
sd.obs<-mean.tmin
sd.obsfit<-mean.tmin
mean.obsfit<-mean.tmin
sd.fcfit<-mean.tmin
mean.fcfit<-mean.tmin
tot.ensn<-(fcyears-1)*fcmembers
t.obs<-matrix(NA,nrow=180,ncol=fcdays)

for (lat in 25) {
  tas.all <- array(NA,dim=c(180,tot.ensn,fcdays))
  t.obs.all<-array(NA,dim=c(180,fcyears-1,fcdays))
  t.obs<-matrix(NA,nrow=180,ncol=fcdays)
  
  for (year in (1:(length(file.names)-1))) {
    #read from yearly erainterim
    doyrange<-seq(as.numeric(format(as.Date(paste(yearstr[year],as.character(initmonth),"01",sep=""),"%Y%m%d"),"%j")),as.numeric(format(as.Date(paste(yearstr[year],"1231",sep=""),"%Y%m%d"),"%j")))
    if (length(doyrange)< fcdays) {
      #i.e. FC range spanning two ERA-Interim years
      t.obs[,1:length(doyrange)]<-get.var.ncdf(obsfilecons[[year]],varid=parobs,start=c(1,lat,doyrange[1]),count=c(-1,1,length(doyrange)))      
      t.obs[,(length(doyrange)+1):fcdays]<-get.var.ncdf(obsfilecons[[year+1]],varid=parobs,start=c(1,lat,1),count=c(-1,1,fcdays-length(doyrange)))     
      
    }else{
      t.obs<-get.var.ncdf(obsfilecons[[year]],varid=parobs,start=c(1,lat,min(doyrange)),count=c(-1,1,fcdays))
      
    }
    
    #read forecasts
    tas  <- get.var.ncdf(connections[[year]],varid=parfc,start=c(1,lat,1,1),count=c(-1,1,-1,-1))
    x0<-year*fcmembers-(fcmembers-1)
    ind.t <- seq(x0, x0+(fcmembers-1),1)
    tas.all[,ind.t,] <- tas
    t.obs.all[,year,]<-t.obs
    if (year==1){ print(paste("lat =",lat))}
  } #year loop
  
  # mean.tmin[,lat,] <- apply(tas.all,c(1,3),mean,na.rm=TRUE)
  mean.tmin[,lat,]<-colMeans(aperm(tas.all,c(2,1,3)))
  sd.tmin[,lat,]<-apply(tas.all,c(1,3),sd,na.rm=TRUE)
  #  mean.obs[,lat,]<-apply(t.obs.all,c(1,3),mean,na.rm=TRUE)
  mean.obs[,lat,]<-colMeans(aperm(t.obs.all,c(2,1,3)))
  sd.obs[,lat,]<-apply(t.obs.all,c(1,3),sd,na.rm=TRUE)   
  
  ## generalise the span for shorter/longer forecasts
  span <- 120 / fcdays  
  sdspan <- 151 / fcdays
  for (lon in (1:180)){
    mean.obsfit[lon,lat,]<-loess(mean.obs[lon,lat,]~tvec,span=span)$fitted
    sd.obsfit[lon,lat,]<-loess(sd.obs[lon,lat,]~tvec,span=sdspan)$fitted
    mean.fcfit[lon,lat,]<-mean.tmin[lon,lat,] ## no fit due to weird jumps
    sd.fcfit[lon,lat,]<-loess(sd.tmin[lon,lat,]~tvec,span=sdspan)$fitted
  } #lon loop

  #plot(1:215,sd.obsfit[45,1,],type="l")
  #points(1:215,sd.obs[45,1,])
  #lines(1:215,sd.tmin[45,1,],col=2)
  
} # lat loop

#turn data around
mean.tmin[,1:91,] <- mean.tmin[,91:1,]
sd.tmin[,1:91,] <- sd.tmin[,91:1,]
mean.obs[,1:91,]<- mean.obs[,91:1,]
sd.obs[,1:91,]<- sd.obs[,91:1,]
mean.obsfit[,1:91,]<- mean.obsfit[,91:1,]
sd.obsfit[,1:91,]<- sd.obsfit[,91:1,]
mean.fcfit[,1:91,]<- mean.fcfit[,91:1,]
sd.fcfit[,1:91,]<-sd.fcfit[,91:1,]


#---------------------------------
#bias
#--------------------------------
#ind=seq(1,fcdays,1)

#bias <- mean.tmin-tmin.fit[,,ind]
bias<-mean.fcfit-mean.obsfit
sdrat<-sd.obsfit/sd.fcfit


#-----------------------------------------
#debias and save data
#-----------------------------------------

#get model data



for (f in 1:fcyears){
  #for (f in 31:32){
  # file.in <- paste(dirhc,file.names[f],sep='')
  TAS<-get.var.ncdf(connections[[f]],varid=parfc)      
  lon<-get.var.ncdf(connections[[f]],varid="longitude")       
  lat<-get.var.ncdf(connections[[f]],varid="latitude")       
  # TAS <- var.get.ncv(file.in,parfc,data=TRUE,collapse=TRUE)$data
  #  lon <- var.get.ncv(file.in,"longitude",data=TRUE,collapse=TRUE)$data
  #  lat <- var.get.ncv(file.in,"latitude",data=TRUE,collapse=TRUE)$data
  
  TAS[,1:91,,] <- TAS[,91:1,,]
  
  lat[1:91] <- lat[91:1]
  
  system.time({
    tas.debias <- array(NA,dim=c(180,91,51,215))
    tas.anom <- tas.debias 
    ens.mean<- colMeans(aperm(TAS,c(3,1,2,4)))
    tas.anom<- TAS-aperm(replicate(fcmembers,ens.mean),c(1,2,4,3))  
    #spreadcorr
    tas.anom<-tas.anom*aperm(replicate(fcmembers,sdrat),c(1,2,4,3))
    tas.debias<-aperm(replicate(fcmembers,ens.mean),c(1,2,4,3))-aperm(replicate(fcmembers,bias),c(1,2,4,3))+tas.anom
  })
  
  bias2 <- array(bias, c(180,91,1,215))  
  sdrat2 <- array(sdrat, c(180,91,1,215))
  system.time({
    ens.mean <- array(colMeans(aperm(TAS,c(3,1,2,4))), c(180,91,1,215))
    tas.anom <- TAS - ens.mean[,,rep(1,51),]
    tas.anom <- tas.anom * sdrat2[,,rep(1,51),]
    tas.debias <- (ens.mean - bias2)[,,rep(1,51),]  + tas.anom
  }
  )
  
  #  lonw <- coord.def.ncv('longitude',data=lon,att=list('long_name','longitude','units','degrees_east'))
  #  latw <- coord.def.ncv('latitude',data=lat,att=list('long_name','latitude','units','degrees_north'))
  #  ens <- coord.def.ncv('ens',data=c(0:50),att=list('long_name','ensemble members','units','-'))
  lonw <- dim.def.ncdf(name="longitude",units="degrees_east",vals=lon)
  latw <- dim.def.ncdf(name='latitude',units="degrees_north",vals=lat)
  ens <- dim.def.ncdf(name='ens',units="-",vals=as.integer(c(0:50)))
  time<- dim.def.ncdf(name='time',units="days",vals=as.integer(c(0:214)),unlim=TRUE) 
  
  
  #  time <- coord.def.ncv("time",data=c(0:214),att=list("axis", "T", "calendar", "standard", "long_name", "time","units", "days"))
  data<- var.def.ncdf(name=parout,units="K",dim=list(lonw,latw,ens,time),missval=-999,longname="debiased model data")
  # data <- var.def.ncv(name=parout,data=tas.debias,dim=list(lonw,latw,ens,time),att=list('long_name','debiased model data','units','K','_FillValue', -999))
  ncout<-create.ncdf(outfile.names[f],data) 
  put.var.ncdf(ncout,varid=data,vals=tas.debias)
  close.ncdf(ncout)
  #var.put.ncv(paste(dir.out,file.names[f],sep=''),data,new=TRUE,define=TRUE)
}#fcyears



