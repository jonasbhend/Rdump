## read in ECAnD data and convert to blended

stations <- read.table("/store/msclim/bhendj/tmp/ECAnD/ECA_blend_station_tg.txt", 
                       skip=16, sep=',', header=TRUE)

degfun <- function(x) {
  deg <- as.numeric(x)
  return(sign(deg[1]) * (abs(deg[1]) + deg[2]/60 + deg[3]/3600))
}

stations$lon <- sapply(strsplit(stations$LON, ':'), degfun)
stations$lat <- sapply(strsplit(stations$LAT, ':'), degfun)

stmp <- read.fwf("/store/msclim/bhendj/tmp/ECAnD/ECA_blend_source_tg.txt",
                    widths=c(5,-1,6,-1,40,-1,2,-1,9,-1,10,-1,4,-1,4,-1,8,-1,8,-1,5,-1,51),
                    skip=23)
sources <- as.data.frame(stmp[-(1:2),])
names(sources) <- gsub(" ", "", stmp[1,])
sources$station <- gsub("[ ]*$", "", sources$SOUNAME)
sources$contributor <- gsub("[ ]*$", "", sources$PARNAME)
sources$startdate <- as.Date(sources$START, format='%Y%m%d')
sources$enddate <- as.Date(sources$STOP, format='%Y%m%d')
sources$lon <- sapply(strsplit(sources$LON, ':'), degfun)
sources$lat <- sapply(strsplit(sources$LAT, ':'), degfun)


