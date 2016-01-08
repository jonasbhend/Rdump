library(reshape2)

## read in ghcnd stations
read.ghcnd <- function(file, par='TMAX', 
                       daterange=as.Date(c('1981-01-01', '2014-12-31')), 
                       ...){
  fwidths <- c(11,4,2,4,rep(c(5,1,1,1), 31))
  ftmp <- read.fwf(file, widths=fwidths, na.string='-9999', ...)
  names(ftmp) <- c('id', 'year', 'month', 'element', outer(c('value', 'mflag', 'qflag', 'sflag'), 1:31, paste0))
  ftmp <- subset(ftmp, element == par)  
  if (nrow(ftmp) == 0) return(NULL)
  ## convert to time series
  ftmp2 <- reshape(ftmp, varying=list(value=grep('value', names(ftmp)), 
                                      mflag=grep('mflag', names(ftmp)),
                                      qflag=grep('qflag', names(ftmp)),
                                      sflag=grep('sflag', names(ftmp))),
                   v.names=c('value', 'mflag', 'qflag', 'sflag'),
                   direction='long', idvar='day')
  ftmp2$date <- as.Date(apply(as.matrix(ftmp2[,c('year', 'month', 'time')]), 1, paste, collapse='-'))
  ftmp2 <- ftmp2[!is.na(ftmp2$date),]
  if (nrow(ftmp2) == 0) return(NULL)
  ftmp2 <- subset(ftmp2, date >= daterange[1] & date <= daterange[2])
  if (nrow(ftmp2) == 0) return(NULL)
  ftmp2 <- ftmp2[order(ftmp2$date),]
  return(ftmp2)
}


## read GHCN station data

stat <- read.fwf("/store/msclim/bhendj/tmp/GHCN/ghcnd-station2.txt", 
                 c(11, -1, 8, -1, 9, -1, 6, -1, 2, -1, 30, -1, 3, -1, 3, -1, 5))
names(stat) <- c('id', 'lat', 'lon', 'altitude', 'state', 'name', 'gsn', 'hcn', 'wmo')

inventory <- read.table("/store/msclim/bhendj/tmp/GHCN/ghcnd-inventory.txt")
names(inventory) <- c('id', 'lat', 'lon', 'element', 'start', 'stop')
tmax.inv <- subset(inventory, element == 'TMAX' & start <= 1981 & stop >=2014)
tmin.inv <- subset(inventory, element == 'TMIN' & start <= 1981 & stop >=2014)

subinv <- merge(tmax.inv, tmin.inv, by=c('id', 'lat', 'lon'))

## files in gsn database
gfiles <- list.files('/store/msclim/bhendj/tmp/GHCN/ghcnd_gsn', full.names=TRUE)
gsn.id <- gsub('.dly', '', basename(gfiles))

subinv <- subset(subinv, id %in% gsn.id)


dates <- seq(as.Date('1981-01-01'), as.Date('2015-02-28'), 1)

gtmax <- list()

for (i in 1:nrow(subinv)){
  print(subinv$id[i])
  f <- paste0('/store/msclim/bhendj/tmp/GHCN/ghcnd_all/', subinv$id[i], '.dly')
  ftmp <- read.ghcnd(f, par='TMAX', daterange=range(dates))
  if (sum(!is.na(ftmp$value)) > 0.9*length(dates)){
    gtmax[[subinv$id[i]]] <- ftmp
  }
}

## fraction of non-missing values
sinv <- subset(inventory, id %in% names(gtmax) & element == 'TMAX')
sinv$nna <- NA
sinv$nna[match(names(gtmax), sinv$id)] <- sapply(gtmax, function(x) sum(!is.na(x$value))/length(dates))
rownames(sinv) <- sinv$id
sinv <- sinv[names(gtmax),]

sinv$nna2 <- sapply(gtmax, function(x) sum(!is.na(subset(x, date < as.Date('2015-01-01'))$value)) / length(dates[dates < as.Date('2015-01-01')]))


gtmin <- list()
for (id in names(gtmax)){
  print(paste(which(names(gtmax) == id), id, sep=': '))
  f <- paste0('/store/msclim/bhendj/tmp/GHCN/ghcnd_all/', id, '.dly')
  ftmp <- read.ghcnd(f, par='TMIN', daterange=range(dates))
  if (sum(!is.na(ftmp$value)) > 0.9*length(dates)){
    gtmin[[id]] <- ftmp
  }  
}


inames <- c('id', 'date', 'value')
gmean <- list()
for (id in names(gtmin)){
  dtmp <- merge(gtmax[[id]][,inames], gtmin[[id]][,inames], by=inames[1:2])
  names(dtmp)[3:4] <- c('tmax', 'tmin')
  dtmp$tmean <- (dtmp$tmax + dtmp$tmin) / 2
  gmean[[id]] <- dtmp[,c('date', 'tmean')]
  names(gmean[[id]])[2] <- id 
}

allmean <- Reduce(merge, gmean)

allmean <- gmean[[1]]
for (i in 2:length(gmean)) allmean <- merge(allmean, gmean[[i]],

save(sinv, allmean, file='GHCN_tmp.Rdata')










gfiles <- list.files("/store/msclim/bhendj/tmp/GHCN/ghcnd_gsn", full=TRUE)
gdata <- list()
for (f in gfiles) {
  print(paste(which(gfiles == f), f, sep=': '))
  gdata[[gsub('.dly', '', basename(f))]] <- read.ghcnd(f, par='TMAX')
}

## count number of missing values in period from 1981-2014
nvals <- sapply(gdata, function(x) nrow(subset(x, date >= as.Date('1981-01-01') & date <= as.Date('2014-12-31'))))
nnotna <- sapply(gdata, function(x) sum(!is.na(subset(x, date >= as.Date('1981-01-01') & date <= as.Date('2014-12-31'))$value)))










