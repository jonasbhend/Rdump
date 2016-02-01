## script to compute skill metrics from debiased or raw forecast metrics

## load packages
library(myhelpers)
library(easyVerification)

scores <- c('FairRpss', 'FairCrpss', 'EnsRmsess')

scorelist <- list(Ens2AFC='generalized discrimination score', 
                  EnsCorr='correlation', 
                  EnsCrpss='continuous ranked probability skill score', 
                  FairCrpss='fair continuous ranked probability skill score', 
                  FairCrps='mean fair continuous ranked probability score',
                  EnsMae='mean absolute error',
                  EnsMaess='mean absolute error skill score',
                  EnsMe='mean error (bias)',
                  EnsMse='mean squared error',
                  EnsMsess='mean squared error skill score',
                  EnsRmse='root mean squared error',
                  EnsRmsess='root mean sqaured error skill score',
                  EnsRps='ranked probability score',
                  EnsRpss='ranked probability skill score',
                  FairRpss='fair ranked probability skill score',
                  EnsSprErr='spread error ratio',
                  EnsTrend="mean error in linear trend",
                  EnsCond="conditional error",
                  EnsVarlog="log ratio of variances",
                  EnsStdfrac="fraction of standard deviation",
                  FairSprErr='fair spread error ratio',
                  EnsRocss.cat1='ROC area skill score (lower tercile)',
                  EnsRocss.cat2='ROC area skill score (middle tercile)',
                  EnsRocss.cat3='ROC area skill score (upper tercile)',
                  EnsIgnss="ignorance skill score")


## scorefunction for use with veriApply
scorefun <- function(f, o, score, ref, prob=NULL){
  if (length(grep('Rps', score)) == 1 | score == 'EnsRocss'){
    prob <-  c(1/3,2/3)    
  }
  return(veriApply(score, f, o, fcst.ref=ref, ensdim=4, tdim=5, prob=prob))
}


## initialise
scratchdir <- Sys.getenv('SCRATCH')

## read command line arguments:
args <- commandArgs(trailingOnly=TRUE)
## check if there are any command line arguments
if (length(args) == 0){
  args <- c('DWD-CCLM4-8-21', 
            'SMHI-EC-EARTH',
            'WFDEI',
            'tasmax',
            'EAF-22',
            'fastqqmap_????-????_WFDEI',
            '05',
            TRUE,
            FALSE,
            FALSE)
} else if (length(args) < 7){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
model <- args[1]
reference <- args[2]
obsname <- args[3]
index <- args[4]
grid <- args[5]
method <- args[6]
initmon <- args[7]

if (length(args) >= 9){
  seasonals <- as.logical(args[8])
  ccrs <- as.logical(args[9])
} else {
  seasonals <- c(FALSE, TRUE)
  ccrs <- c(FALSE, TRUE)
}
if (length(args) == 10){
  detrends <- as.logical(args[10])
} else {
  detrends <- c(FALSE, TRUE)  
}


## replace placeholder for method string
dpath <- '/store/msclim/bhendj/EUPORIAS'
if (length(grep('????-????', method)) == 1 & method != 'none-forward'){
  mfiles <- system(paste0('ls ', dpath, '/', model,'/', grid, '/monthly/', index, '/', method, '/', index, '_????', initmon, '??_*.nc'), intern=TRUE)
  method <- sapply(strsplit(mfiles, '/'), function(x) x[length(x) - 1])[1]
}
stopifnot(is.character(method))

## specify the file paths
fpath <- paste(dpath, model, grid, 'monthly', index, method, sep='/')
rpath <- paste(dpath, reference, grid, 'monthly', index, method, sep='/')
opath <- paste(dpath, obsname, grid, 'monthly', index, sep='/')

## get forecast and observation files
fcfiles <- list.files(fpath, pattern=paste0('^', index, '_....', initmon, '.._.*_',grid,'_',method, '.nc'), full.name=TRUE)
reffiles <- list.files(rpath, pattern=paste0('^', index, '_....', initmon, '.._.*_', grid, '_', method, '.nc'), full.name=TRUE)
if (method == "none-forward"){
  fpath <- paste(dpath, model, grid, 'monthly', index, 'none', sep='/')
  fcfiles <- list.files(fpath, pattern=paste0('^', index, '_....', initmon, '.._.*_',grid,'_none.nc'), full.name=TRUE)
  rpath <- paste(dpath, reference, grid, 'monthly', index, 'none', sep='/')
  reffiles <- list.files(rpath, pattern=paste0('^', index, '_....', initmon, '.._.*_', gird, '_none.nc'), full.name=TRUE)
}
obsfiles <- list.files(opath, pattern=paste0('^', index, '_'), full.name=TRUE)

## check whether there is at least one of each
if (length(fcfiles) < 1) stop('no forecast file found')
if (length(reffiles) < 1) stop("no reference forecast found")
if (length(obsfiles) < 1) stop('no observation file found')

## only use most recent obsfile
details <- file.info(obsfiles)
details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)),]
obsfiles <- rownames(details)[1]
rm(details)

## enquire file properties (names and dims from forecast)
nc <- nc_open(fcfiles[1])
## get variable name in forecasts
## parfc <- names(nc$var)[names(nc$var) %in% varlist[[index]]]
parfc <- index
## get dimensions
fdims <- nc$var[[parfc]]$dim
fdimnames <- sapply(fdims, function(x) x$name)
names(fdims) <- fdimnames
if (length(fdimnames) == 3){
  nlon <- fdims[[grep('ncells', fdimnames)]]$len
  nlat <- 1
} else {
  nlon <- fdims[[grep('lon', fdimnames)]]$len
  nlat <- fdims[[grep('lat', fdimnames)]]$len  
}
nens <- fdims[[which(!fdimnames %in% grep('tim|lon|lat|ncells|nb2', fdimnames, value=TRUE))]]$len
## get more information on time (i.e. initialisation month)
time.i <- grep('tim', fdimnames)
ntime <- fdims[[time.i]]$len
fctime <- as.Date(nc_time(nc))
## number of months with all days present
initmn <- as.numeric(format(fctime[1], '%m'))
nc_close(nc)

ncref <- nc_open(reffiles[1])
refdims <- nc$var[[parfc]]$dim
refdimnames <- sapply(refdims, function(x) x$name)
names(refdims) <- refdimnames
nrefens <- refdims[[which(!refdimnames %in% grep('tim|lon|lat|ncells|nb2', refdimnames, value=TRUE))]]$len
reftime <- as.Date(nc_time(ncref))
ncomplete <- min(sum(as.numeric(format(fctime, '%d')) >= 28),
                 sum(as.numeric(format(reftime, '%d')) >= 28))
nc_close(ncref)

## name from observations
nc <- nc_open(obsfiles[1])
## parobs <- names(nc$var)[names(nc$var) %in% varlist[[index]]]
parobs <- index
obstime <- as.Date(nc_time(nc))
nc_close(nc)

## set up file connections
fc.con <- lapply(as.list(fcfiles), nc_open)
ref.con <- lapply(as.list(reffiles), nc_open)
obs.con <- lapply(as.list(obsfiles), nc_open)

## get forecast times
fc.times <- lapply(fc.con, function(x) as.Date(nc_time(x)[1:ncomplete] ))
names(fc.times) <- sapply(fc.times, function(x) format(x[1], '%Y'))
ref.times <- lapply(ref.con, function(x) as.Date(nc_time(x)[1:ncomplete] ))
names(ref.times) <- sapply(ref.times, function(x) format(x[1], '%Y'))

## check whether forecasts are in obs
if (model == 'DWD-CCLM4-8-21') fc.times <- lapply(fc.times, function(x) x + c(0,0,0,0,1))
is.in.obs <- sapply(fc.times, function(x) all(x[1:ncomplete] %in% obstime & x[1:ncomplete] %in% unlist(ref.times)))
if (mean(is.in.obs) < 0.5){
  obstime <- obstime + 1
  is.in.obs <- sapply(fc.times, function(x) all(x[1:ncomplete] %in% obstime & x[1:ncomplete] %in% unlist(ref.times)))
}
ref.in.obs <- sapply(ref.times, function(x) all(x[1:ncomplete] %in% obstime & x[1:ncomplete] %in% unlist(fc.times)))

stopifnot(sum(!is.in.obs) < 3)

## reduce to overlapping set
fcfiles <- fcfiles[is.in.obs]
fc.con <- fc.con[is.in.obs]
fc.times <- fc.times[is.in.obs]
reffiles <- reffiles[ref.in.obs]
ref.con <- ref.con[ref.in.obs]
ref.times <- ref.times[ref.in.obs]

## number of years
nyears <- length(fcfiles)
years <- as.numeric(names(fc.times))
if(length(grep('none', method)) == 1) {
  myears <- years  
} else {
  mtmp <- strsplit(method, '_')[[1]]
  mtmp <- as.numeric(strsplit(mtmp[length(mtmp) - 1], '-')[[1]])
  myears <- mtmp[1]:mtmp[2]
} 

## get the forecast ensemble members
nens <- min(sapply(fc.con, function(x) x$dim[[which(! names(x$dim) %in% grep("tim|lon|lat|ncells|nb2", names(x$dim), value=TRUE))]]$len)) 
nrefens <- min(sapply(ref.con, function(x) x$dim[[which(! names(x$dim) %in% grep("tim|lon|lat|ncells|nb2", names(x$dim), value=TRUE))]]$len)) 

## read in the observations (corresponding time steps)
if (length(obsfiles) == 1){
  obs <- array(NA, c(nlon, nlat, ncomplete, nyears))
  fc.i <- sapply(fc.times, function(x) format(x, '%Y%m') %in% unlist(sapply(obstime, format, '%Y%m')) & seq(along=as.numeric(x)) %in% 1:ncomplete)
  obs.i <- sapply(obstime, function(x) format(x, '%Y%m') %in% unlist(sapply(fc.times, function(x) format(x[1:ncomplete], '%Y%m'))))    
  if (length(obs.con[[1]]$dim) == 2){
    obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,obs.i]
  } else {
    obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,,obs.i]    
  }
  obs[fc.i[rep(1:nrow(fc.i), each=nlon*nlat), ]] <- as.vector(obs.tmp)
  obs <- aperm(obs, c(3,1,2,4))
} else {
  stop('Multiple observation files not implemented yet')
}

## fix units
if (index %in% c('tas', 'tasmax', 'tasmin') & max(obs, na.rm=T) < 200){
  obs <- obs + 273.15
}

## read in all the data
ffcst <- array(NA, c(ncomplete, nlon, nlat, nens, nyears))
rref <- array(NA, c(ncomplete, nlon, nlat, nrefens, nyears))
if (length(fc.con[[1]]$dim) == 3) {
  for (i in 1:nyears){
    ffcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,nens,ncomplete)), c(3,1,2))
    rref[,,,,i] <- aperm(ncvar_get(ref.con[[i]], parfc, count=c(-1,nrefens,ncomplete)), c(3,1,2))
  }
} else {
  for (i in 1:nyears){
    ffcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,-1,nens,ncomplete)), c(4,1,2,3))
    rref[,,,,i] <- aperm(ncvar_get(ref.con[[i]], parfc, count=c(-1,-1,nrefens,ncomplete)), c(4,1,2,3))
  }  
}
rownames(ffcst) <- unique(sapply(fc.times, format, '%m'))[1:ncomplete]
dimnames(ffcst)[[5]] <- names(fc.times)
rownames(rref) <- unique(sapply(ref.times, format, '%m'))[1:ncomplete]
dimnames(rref)[[5]] <- names(ref.times)

## get the grid data for display of forecasts
nc <- fc.con[[1]]
if (length(grep('eobs', grid)) == 1){
  lon <- nc$dim$rlon$vals
  lat <- nc$dim$rlat$vals
  plon <- ncatt_get(nc, 'rotated_pole', 'grid_north_pole_longitude')$val
  plat <- ncatt_get(nc, 'rotated_pole', 'grid_north_pole_latitude')$val
} else if (length(nc$dim) == 3) {
  lon <- ncvar_get(nc, 'lon')
  lat <- ncvar_get(nc, 'lat')
} else {
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals        
  plon <- plat <- NULL
}
rm(nc)

## close file connections (desperate debugging)
sapply(fc.con, nc_close)
sapply(ref.con, nc_close)
sapply(obs.con, nc_close)
rm(fc.con, ref.con, obs.con)

## get crossvalidation string
crossval <- length(grep('crossval', method)) == 1
suppressWarnings(nblock <- as.numeric(gsub('_.*', '', gsub('.*crossval', '', method))))
forward <- length(grep("forward", method)) == 1
stopifnot(sum(forward, crossval) <= 1)

## loop over seasonal or monthly forecast skill
for (seasonal in seasonals){
  
  outdir <- paste(dpath, 'skill_against_reference', grid, if (seasonal) 'seasonal' else 'monthly', index, sep='/')
  
  if (seasonal){
    ## compute three-monthly aggregates
    if (grid == 'EAF-22'){
      nred <- 0
    } else {
      nred <- 2
    }
    fcst.seas <- array(NA, dim(ffcst) - c(nred,0,0,0,0))
    ref.seas <- array(NA, dim(rref) - c(nred,0,0,0,0))
    obs.seas <- array(NA, dim(obs) - c(nred,0,0,0))
    for (i in 1:(nrow(ffcst) - 2)){
      fcst.seas[i,,,,] <- colMeans(ffcst[i+0:2,,,,,drop=F], dims=1)
      ref.seas[i,,,,] <- colMeans(rref[i+0:2,,,,,drop=F], dims=1)
      obs.seas[i,,,] <- colMeans(obs[i + 0:2,,,,drop=F], dims=1)
    }
    if (grid == 'EAF-22'){
      fcst.seas[nrow(fcst.seas),,,,] <- colMeans(ffcst[2:5,,,,,drop=F], dims=1)
      ref.seas[nrow(fcst.seas),,,,] <- colMeans(rref[2:5,,,,,drop=F], dims=1)
      obs.seas[nrow(fcst.seas),,,] <- colMeans(obs[2:5,,,,drop=F], dims=1)
    }
    ## rearrange time for output
    fcst.seastimes <- sapply(fc.times, function(x) format(if (grid == 'EAF-22') x[setdiff(1:ncomplete, 1:2)][c(seq(1,ncomplete-2), NA, 3)] else x[setdiff(1:ncomplete, 1:2)], '%Y-%m-%d'))
    fcst.seastimes <- array(0, dim(fcst.seastimes)) + as.Date(fcst.seastimes)  
  } else {
    fcst.seas <- ffcst[1:ncomplete,,,,,drop=F]
    ref.seas <- rref[1:ncomplete,,,,,drop=F]
    obs.seas <- obs[1:ncomplete,,,,drop=F]
    fcst.seastimes <- sapply(fc.times, function(x) format(x[1:ncomplete], '%Y-%m-%d'))
    fcst.seastimes <- array(0, dim(fcst.seastimes)) + as.Date(fcst.seastimes)  
  }
  
  ## mask forecasts and obs where more than 80% of the years are missing
  fcst.mask <- rowMeans(is.na(fcst.seas)*1,  dims=3) <= .2
  ref.mask <- rowMeans(is.na(ref.seas)*1, dims=3) <= .2
  obs.mask <- rowMeans(is.na(obs.seas)*1,  dims=3) <= .2
  
  ## mask according to observation climatology
  if (index %in% c('HDD', 'CDD', 'HDDch')){
    ## new criterion at least 67% of forecasts in climatology are non-zero
    nclim <- min(30, dim(obs.seas)[4])
    obs.mask <- obs.mask & rowMeans(obs.seas[,,,1:nclim] > 0, dims=3, na.rm=T) > 0.67
    ##obs.mask <- obs.mask & rowMeans(obs.seas, dims=3, na.rm=T) > 0.1
    ## this amounts to ~3 HDD per month, or 9 HDD per season
  } else if (index %in% c('FD', 'SD', 'HD')){
    obs.mn <- rowMeans(obs.seas, dims=3, na.rm=T)
    obs.mask <- obs.mask & obs.mn > 0.05 & obs.mn < 0.95
  } else if (index == 'WDF') {
    nclim <- min(30, dim(obs.seas)[4])
    obs.mask <- obs.mask & rowMeans(obs.seas[,,,1:nclim] > 0, dims=3, na.rm=T) > 0.67 & rowMeans(obs.seas[,,,1:nclim] < 1, dims=3, na.rm=T) > 0.67
  }
  
  ## like this, all years are masked out (and therefore the corresponding scores also)
  obs.seas[rep(!(obs.mask & fcst.mask & ref.mask), length=length(obs.seas))] <- NA
  
  ## set up temporary directory
  tmpdir <- paste0(scratchdir, '/skill_', index, '_', grid, '_', method, '_initmon', initmon, '_', ceiling(runif(1, min=0, max=9999)))
  if (! file.exists(tmpdir)) dir.create(tmpdir)
  
  ## loop through raw / detrended forecasts and obs
  for (detrend in detrends){
    
    if (detrend){
      ## detrend obs
      obs.seas <- aperm(apply(obs.seas, 1:3, function(y){
        if (all(!is.na(y))){
          mean(y) + lm(y ~ seq(along=y))$res
        } else {
          NA*y
        }
      }), c(2,3,4,1))
      ## detrend forecasts
      fcst.seas <- aperm(apply(fcst.seas, 1:4, function(y){
        if(all(!is.na(y))){
          mean(y) + lm(y ~ seq(along=y))$res
        } else {
          NA * y
        }
      }), c(2,3,4,5,1))
      ref.seas <- aperm(apply(ref.seas, 1:4, function(y){
        if(all(!is.na(y))){
          mean(y) + lm(y ~ seq(along=y))$res
        } else {
          NA * y
        }
      }), c(2,3,4,5,1))
    }
    
    ## loop through raw / CCR forecasts
    for (is.ccr in ccrs){
      
      if (is.ccr){
        fcst.ccr <- ref.ccr <- array(NA, dim(fcst.seas))
        for (loi in 1:nlon){
          for (lai in 1:nlat){
            if (any(!is.na(obs.seas[,loi,lai,]))){
              fcst.ccr[,loi,lai,,] <- aperm(debias(
                fcst=aperm(fcst.seas[,loi, lai,,which(years %in% myears)], c(1,3,2)),
                obs=obs.seas[,loi,lai,which(years %in% myears)],
                fcst.out=aperm(fcst.seas[,loi,lai,,], c(1,3,2)),
                method='ccr', 
                crossval=crossval,
                blocklength=nblock, 
                forward=forward, 
                type= if(forward | crossval) "prediction" else "calibration"), c(1,3,2))
              ref.ccr[,loi,lai,,] <- aperm(debias(
                fcst=aperm(ref.seas[,loi, lai,,which(years %in% myears)], c(1,3,2)),
                obs=obs.seas[,loi,lai,which(years %in% myears)],
                fcst.out=aperm(ref.seas[,loi,lai,,], c(1,3,2)),
                method='ccr', 
                crossval=crossval,
                blocklength=nblock, 
                forward=forward, 
                type= if(forward | crossval) "prediction" else "calibration"), c(1,3,2))
            } ## end of if on missing values in obs
          } ## end of loop on latitudes
        } ## end of loop on longitudes
        
        fcst.seas <- fcst.ccr
        ref.seas <- ref.ccr
      }
      
      ## get number of years that are analysed (reduced for forward calibration)
      nyears <- dim(obs.seas)[4]
      ## replace if not succesful
      ## yind <- seq(if (forward) 16 else 1, nyears)
      yind <- seq(1,nyears)
      

      ## set up output file to write to
      nc <- nc_open(obsfiles[1])
      nc2 <- nc_open(fcfiles[1])
      if (length(nc$dim) == 2){
        ncells.nc <- nc$dim$ncells
      } else if (length(grep('eobs', grid)) == 1){
        lon.nc <- nc2$dim$rlon
        lat.nc <- nc2$dim$rlat    
      } else {
        lon.nc <- nc2$dim$lon
        lat.nc <- nc2$dim$lat
      }
      time.nc <- ncdim_def(name='time', 
                           units='days since 1980-01-01', 
                           vals=as.numeric(fcst.seastimes[, ncol(fcst.seastimes)] - as.Date('1980-01-01')),
                           unlim=TRUE, 
                           longname=if (seasonal) 'end date of 3-monthly seasons' else 'end date of month')
      if (length(nc$dim) == 2){
        dims.nc <- list(ncells.nc, time.nc)
      } else {
        dims.nc <- list(lon.nc, lat.nc, time.nc)      
      }
      ## create variables
      vars.nc <- list()
      ## check whether this is a rotated grid
      if (length(grep('eobs', grid)) == 1){
        vars.nc[['lon']] <- ncvar_def('lon', 'degrees_east', dims.nc[1:2], missval=-1e20, longname='longitude')
        vars.nc[['lat']] <- ncvar_def('lat', 'degrees_north', dims.nc[1:2], missval=-1e20, longname='latitude')
        vars.nc[['rotated_pole']] <- ncvar_def('rotated_pole', units='', dim=list(), missval=NULL, prec='char')    
      } else if (length(nc$dim) == 2){
        vars.nc[['lon']] <- nc$var$lon
        vars.nc[['lat']] <- nc$var$lat
      }
      for (score in scores){
        if (score == 'EnsRocss'){
          for (i in 1:3) {
            rscore <- paste0(score, '.cat', i)
            vars.nc[[rscore]] <- ncvar_def(rscore, '1', dims.nc, missval=-1e20, longname=scorelist[[rscore]])
            vars.nc[[paste0(rscore, '.sigma')]] <- ncvar_def(paste0(rscore, '.sigma'), '1',
                                                             dims.nc, missval=-1e20,
                                                             longname=scorelist[[rscore]])
          }
        } else {
          vars.nc[[score]] <- ncvar_def(score, '1', dims.nc, missval=-1e20, longname=scorelist[[score]])
          if (score %in% c(scores[grep('pss$', scores)], 'EnsTrend', 'EnsCond', 'EnsVarlog')){
            vars.nc[[paste(score, 'sigma', sep='.')]] <- ncvar_def(paste(score, 'sigma', sep='.'), '1', dims.nc, missval=-1e20, longname=scorelist[[score]])
          }      
        }
      }
      
      
      ## get additional information on content of files
      fc.name <- sapply(strsplit(fcfiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
      ref.name <- sapply(strsplit(reffiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
      obs.name <- sapply(strsplit(obsfiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
      ofile <- paste0(paste(paste0(index, if (detrend) '_detrend' else '', if(is.ccr) '_CCR' else ''), method, paste0(fc.name, '-ref-', ref.name), 'vs', obs.name, paste(range(names(fc.times)), collapse='-'), sep="_"), '_initmon', initmon, '.nc')
      outfile <- paste(tmpdir, ofile, sep='/')
      desc <- paste0('Skill in ', if (detrend) 'detrended ' else '', if (seasonal) 'seasonal ' else 'monthly ', index, ' from bias-corrected (', method, if (is.ccr) ' + CCR' else '', if (crossval) ' cross-validated' else '', ') ', gsub('-', ' ', toupper(fc.name)), ' forecasts verified against ', obs.name, ' with reference forecast ', ref.name, ' for ', paste(range(names(fc.times)[yind]), collapse='-'))
      
      ## create netcdf file
      ncout <- nc_create(outfile, vars.nc)
      ## write variables
      if (length(nc$dim) == 2){
        ncvar_put(ncout, varid='lon', ncvar_get(nc, 'lon'))
        ncvar_put(ncout, varid='lat', ncvar_get(nc, 'lat'))
      } else if (length(grep('eobs', grid)) == 1){
        if ('Actual_longitude' %in% names(nc$var)){
          ncvar_put(ncout, varid='lon', ncvar_get(nc, 'Actual_longitude'))
          ncvar_put(ncout, varid='lat', ncvar_get(nc, 'Actual_latitude'))
        } else {
          library(geocors)
          plon <- ncatt_get(nc, 'rotated_pole', attname='grid_north_pole_longitude')$val
          plat <- ncatt_get(nc, 'rotated_pole', attname='grid_north_pole_longitude')$val    
          lola <- geocors.trafo(rep(nc$dim$rlon$vals, nc$dim$rlat$len), 
                                rep(nc$dim$rlat$vals, each=nc$dim$rlon$len),
                                from.type='rotpol', 
                                from.pars=list(plon=as.numeric(plon), plat=as.numeric(plat)),
                                to.type='lonlat')
          ncvar_put(ncout, varid='lon', lola$lon)
          ncvar_put(ncout, varid='lat', lola$lat)    
        }
        ## fill in the rotated grid coordinates
        ncatt_put(ncout, varid='rotated_pole', attname='grid_mapping_name', 
                  attval='rotated_latitude_longitude', prec='text')
        ncatt_put(ncout, varid='rotated_pole', attname='grid_north_pole_longitude', 
                  attval=-162, prec='double')
        ncatt_put(ncout, varid='rotated_pole', attname='grid_north_pole_latitude', 
                  attval=39.25, prec='double')
      }
      ## add in a global attribute with description
      ncatt_put(ncout, varid=0, attname='description', attval=desc, prec='text')
      nc_close(ncout)
      
      ## compute scores and write to output file
      for (score in scores){
        ncout <- nc_open(outfile, write=TRUE)
        print(score)
        print(system.time( sfo <- scorefun(fcst.seas[,,,,yind, drop=F], obs.seas[,,,yind,drop=F], score, ref=ref.seas[,,,,yind,drop=F])))
        if (is.list(sfo)){
          sfo <- lapply(sfo, function(x){
            x[x == -Inf] <- -9999
            x[abs(x) == Inf] <- NA
            return(x)
          })
          if (score == 'EnsRocss'){
            for (i in 1:3){
              ncvar_put(ncout, varid=paste0(score, '.cat', i), vals=aperm(sfo[[paste0('cat', i)]], c(2,3,1)))
              ncvar_put(ncout, varid=paste0(score, '.cat', i, '.sigma'), vals=aperm(sfo[[paste0('cat', i, '.sigma')]], c(2,3,1)))
            }
          } else {
            ncvar_put(ncout, varid=score, vals=aperm(sfo[[1]], c(2,3,1)))
            ncvar_put(ncout, varid=paste(score, 'sigma', sep='.'), vals=aperm(sfo[[2]], c(2,3,1)))    
          }
        } else {
          ## set negative infinity (division by zero) to large negative value
          sfo[sfo == -Inf] <- -9999
          sfo[abs(sfo) == Inf] <- NA
          ## write to output file
          if (length(sfo) == length(obs.seas[,,,yind])){
            sfo <- rowMeans(sfo, dims=3)
          }
          ncvar_put(ncout, varid=score, vals=aperm(sfo, c(2,3,1))) 
        }
        nc_close(ncout)      
        
        ## readjust grid specification for variables
        ## this is very important to happen AFTER the file connection to the NetCDF
        ## is closed from R - otherwise the output is shifted latitudinally
        if (length(nc$dim) == 2 | length(grep('eobs', grid)) == 1){
          if (length(nc$dim) == 2){
            atts <- c(gridtype='unstructured', coordinates="lon lat")
          } else {
            atts <- c(grid_mapping='rotated_pole')
          }
          for (attn in names(atts)){
            if (score == 'EnsRocss'){
              for (i in 1:3) system(paste0("ncatted -h -a ", attn, ",", score, ".cat", i, ",c,c,'",atts[attn],"' ", outfile))            
            } else if (is.list(sfo)){
              system(paste0("ncatted -h -a ", attn, ",", score, ".sigma,c,c,'",atts[attn],"' ", outfile))      
              system(paste0("ncatted -h -a ", attn, ",", score, ",c,c,'",atts[attn],"' ", outfile))      
            } else {
              system(paste0("ncatted -h -a ", attn, ",", score, ",c,c,'",atts[attn],"' ", outfile))      
            }                
          }
        }
        
        ## make sure we free as much memory as possible/reasonable
        rm(sfo, ncout)
        gc()
        
        
      } ## end of loop on scores
      
      ## copy temporary file to permanent location
      ## make sure the output directory exists
      if (!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
      system(paste0('mv ', outfile, ' ', outdir, '/', ofile))
      
    } ## end of loop on CCR
    
  } ## end of loop on detrend
  
  ## remove temporary directory
  system(paste0('rm -rf ', tmpdir))
  
} ## end of loop on seasonal/monthly skill

warnings()
print('Finished successfully')

q(save='no')
