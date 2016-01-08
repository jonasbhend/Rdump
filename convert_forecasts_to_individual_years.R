# this script converts all stored forecast files for visualisation
# to files individually for months and years

Rdfiles <- system("\\ls /store/msclim/bhendj/EUPORIAS/ecmwf-system4/*/seasonal/*/*/*.Rdata", TRUE)


for (file in Rdfiles){

  basedir <- gsub(".Rdata$", "", file) 
  fi <- basename(basedir)
  print(basedir)
  dir.create(basedir)
  setwd(basedir)
  
  load(file)
  years <- sort(unique(as.numeric(format(fcst.seastimes, '%Y'))))
  if (length(years) > 30 & length(years) == dim(fcst.seas)[5])
{
    
    if (length(grep("initmon11", file))) years <- years - 1
    leads <- apply(outer(1:nrow(fcst.seas), 0:2, '+'), 1, paste, collapse='')
    for (yi in seq(years)){
      for (li in seq(leads)){
        prob <- fcst.prob[yi,,li,,]
        fcst <- fcst.seas[li,,,,yi]
        fcst.time <- fcst.seastimes[li,yi]
        save(prob, fcst, fcst.time, 
             file=paste0(basedir, '/', fi, '_', years[yi], '_', leads[li], '.Rdata'))
      }
    }
    
    
    file.remove(file)
  }  else {
    warnings(file)
  }
  
}
