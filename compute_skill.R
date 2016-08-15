## script to compute skill metrics from debiased or raw forecast metrics

## load packages
library(methods) ## for processing with Rscript
library(ncdf4) ## ncdf is deprecated but still supported
library(biascorrection)
## library(veri)
library(easyVerification)

nc_time <- function(nc){
  t.i <- grep('tim', names(nc$dim))
  t.units <- nc$dim[[t.i]]$units
  t.mul <- c(second=1, 
             minute=60, 
             hour=3600,
             day=1, 
             month=1, 
             year=1) ## special cases for days-years
  mul.i <- pmatch(substr(t.units,1,3), names(t.mul))
  if (length(mul.i) == 1){
    if (mul.i < 4){
      t.out <- as.POSIXct(gsub('.*since ', '', t.units), tz='UTC') + t.mul[mul.i]*nc$dim[[t.i]]$val
    } else if (mul.i == 4){
      t.out <- as.POSIXct(round(as.Date(gsub('.*since ', '', t.units), tz='UTC')) + nc$dim[[t.i]]$val, tz='UTC')
    }
  } else {
    stop('Time coordinate format not implemented yet')
  }
  return(t.out)
}


# ## gather list of scores to run
# scores <- setdiff(ls(pos='package:veri'), c('read.in.data', 'save.data', 'roca', 'rocss'))
# scores <- c(scores, outer(c('roca', 'rocss'), 1:3, paste0))
# 
# ## lookup table for long names of scores
# scorelist <- list(corr='correlation', 
#                   crps='continuous ranked probability score', 
#                   crpss='continuous ranked probability skill score',
#                   mae='mean absolute error',
#                   maess='mean absolute error skill score',
#                   me='mean error (bias)',
#                   mse='mean squared error',
#                   msess='mean squared error skill score',
#                   rmse='root mean squared error',
#                   rmsess='root mean sqaured error skill score',
#                   rps='ranked probability score',
#                   rpss='ranked probability skill score',
#                   spr_err='spread error ratio',
#                   roca1='ROC area (lower tercile)',
#                   roca2='ROC area (middle tercile)',
#                   roca3='ROC area (upper tercile)',
#                   rocss1='ROC area skill score (lower tercile)',
#                   rocss2='ROC area skill score (middle tercile)',
#                   rocss3='ROC area skill score (upper tercile)')

scores <- c('EnsRpss', 'FairRpss', 'EnsCrpss', 'FairCrpss', 'EnsCorr', 'EnsMe', 'EnsMse', 'EnsMae', 'EnsRmse', 'EnsRmsess', 'EnsSprErr', 'EnsRocss')

scorelist <- list(EnsCorr='correlation', 
                  EnsCrpss='continuous ranked probability skill score', 
                  FairCrpss='fair continuous ranked probability skill score', 
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
                  EnsRocss.cat1='ROC area skill score (lower tercile)',
                  EnsRocss.cat2='ROC area skill score (middle tercile)',
                  EnsRocss.cat3='ROC area skill score (upper tercile)')


## exclude all scores that produce values for individual years
# scores <- setdiff(scores, c('crps', 'rps', paste0('roca', 1:3)))


# ## add defaults to verification scores
# ## the category to be used is inferred from the name (last digit, eg. roca1)
# scorefun <- function(f, o, score, unbias=FALSE, prob=c(1/3, 2/3), debias=TRUE, res_dec=1, verbose=FALSE){
#   ## get category to be used
#   suppressWarnings(category <- as.numeric(substr(score, nchar(score), nchar(score))))
#   if (!is.na(category)) score <- substr(score, 1, nchar(score) - 1)
#   ## input arguments as a list
#   inargs <- list(fcst=f, obsv=o, unbias=unbias, prob=prob, debias=debias, category=category, res_dec=res_dec, verbose=verbose)
#   ## get the function to be executed
#   sfun <- get(score)
#   ## get the argument names of the function
#   fargnames <- names(formals(sfun))
#   ## call function with the corresponding arguments
#   return(do.call(sfun, inargs[names(inargs) %in% fargnames]))
# }

## scorefunction for use with veriApply
scorefun <- function(f, o, score, prob=NULL){
  if (length(grep('Rps', score)) == 1 | score == 'EnsRocss'){
    prob <-  c(1/3,2/3)    
  }
  return(veriApply(score, f, o, ensdim=4, tdim=5, prob=prob))
}


## initialise
##varlist <- list(tasmin=c('MN2T24', 'mn2t24', 'tmin', 'tasmin', 'tn'),
##                tasmax=c('MX2T24', 'mx2t24', 'tmax', 'tasmax', 'tx'),
##                tas=c('MEAN2T24', 'mean2t24', 'tas', 'tg', 't2m'),
##                pr=c('TOT_PREC', 'pr', 'tp', 'rr'),
##                HDD=c('heating_degree_days_per_time_period'),
##                CDD=c('cooling_degree_days_per_time_period'))
scratchdir <- Sys.getenv('SCRATCH')

## read command line arguments:
args <- commandArgs(trailingOnly=TRUE)
## check if there are any command line arguments
if (length(args) == 0){
  args <- c('tas',
            'smooth-1981-2010',
            '11',
            system('\\ls /store/msclim/bhendj/EUPORIAS/ecmwf-system4/eobs0.44/monthly/tas/smooth_1981_2012/tas_????11*.nc', intern=TRUE),
            '/store/msclim/bhendj/EUPORIAS/E-OBS/monthly/tas/tas_tg_0.44deg_rot_v10.0.nc',
            '/store/msclim/bhendj/EUPORIAS/skill_scores/eobs0.44/seasonal/tas')
} else if (length(args) < 5){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
index <- args[1]
method <- args[2]
initmon <- args[3]
files <- args[-c(1:3, length(args))]
outdir <- args[length(args)]

## further disaggregate the input files (into obs and fcst files)
## based on similarity of directory tree
fbits <- strsplit(files, '/')
## common length of directory structure
nmin <- min(sapply(fbits, length))
direqual <- sapply(fbits, function(x) sum(x[1:nmin] == fbits[[1]][1:nmin]))
fcfiles <- files[direqual == max(direqual)]
obsfiles <- files[direqual < max(direqual)]

## check whether there is at least one of each
if (length(fcfiles) < 1) stop('no forecast file found')
if (length(obsfiles) < 1) stop('no observation file found')

## enquire file properties (names and dims from forecast)
nc <- nc_open(fcfiles[1])
## get variable name in forecasts
## parfc <- names(nc$var)[names(nc$var) %in% varlist[[index]]]
parfc <- index
## get dimensions
fdims <- nc$var[[parfc]]$dim
fdimnames <- sapply(fdims, function(x) x$name)
names(fdims) <- fdimnames
nlon <- fdims[[grep('lon', fdimnames)]]$len
nlat <- fdims[[grep('lat', fdimnames)]]$len
nens <- fdims[[which(!fdimnames %in% grep('tim|lon|lat', fdimnames, value=TRUE))]]$len
## get more information on time (i.e. initialisation month)
time.i <- grep('tim', fdimnames)
ntime <- fdims[[time.i]]$len
fctime <- nc_time(nc)
## number of months with all days present
ncomplete <- sum(as.numeric(format(fctime, '%d')) >= 28)
initmn <- as.numeric(format(fctime[1], '%m'))
nyears <- length(fcfiles)
nc_close(nc)

## name from observations
nc <- nc_open(obsfiles[1])
## parobs <- names(nc$var)[names(nc$var) %in% varlist[[index]]]
parobs <- index
obstime <- nc_time(nc)
nc_close(nc)

## set up file connections
fc.con <- lapply(as.list(fcfiles), nc_open)
obs.con <- lapply(as.list(obsfiles), nc_open)

## get time axes
fc.times <- lapply(fc.con, function(x) nc_time(x)[1:ncomplete] )
names(fc.times) <- sapply(fc.times, function(x) format(x[1], '%Y'))

## read in all the data
fcst <- array(NA, c(ncomplete, nlon, nlat, nens, nyears))
for (i in 1:nyears){
  fcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,-1,-1,ncomplete)), c(4,1,2,3))
}

## read in the observations (corresponding time steps)
obs.times <- lapply(obs.con, function(x) nc_time(x))
if (length(obsfiles) == 1){
  obs <- array(NA, c(nlon, nlat, ncomplete, nyears))
  fc.i <- sapply(fc.times, function(x) format(x, '%Y%m%d') %in% unlist(sapply(obs.times, format, '%Y%m%d')))
  obs.i <- sapply(obs.times, function(x) format(x, '%Y%m%d') %in% unlist(sapply(fc.times, format, '%Y%m%d')))
  obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,,obs.i]
  obs[fc.i[rep(1:nrow(fc.i), each=nlon*nlat), ]] <- as.vector(obs.tmp)
  obs <- aperm(obs, c(3,1,2,4))
} else {
  stop('Multiple observation files not implemented yet')
}

## fix units
if (index %in% c('tas', 'tasmax', 'tasmin') & max(obs, na.rm=T) < 200){
  obs <- obs + 273.15
}

## compute three-monthly aggregates
fcst.seas <- array(NA, dim(fcst) - c(2,0,0,0,0))
obs.seas <- array(NA, dim(obs) - c(2,0,0,0))
for (i in 1:nrow(fcst.seas)){
  fcst.seas[i,,,,] <- colMeans(fcst[i+0:2,,,,], dims=1)
  obs.seas[i,,,] <- colMeans(obs[i + 0:2,,,], dims=1)
}

## remove the original forecast and observations to save space
rm(obs, fcst)
gc()

## mask forecasts and obs where more than 80% of the years are missing
fcst.mask <- rowSums(is.na(fcst.seas)*1,  dims=3) <= 2
obs.mask <- rowSums(is.na(obs.seas)*1,  dims=3) <= 2

## mask according to observation climatology (more than 10 HDD/CDD per season)
if (index %in% c('HDD', 'CDD', 'HDDch')){
  obs.mask <- obs.mask & rowMeans(obs.seas, dims=3, na.rm=T) > 10/90
} else if (index == 'FD'){
  obs.mn <- rowMeans(obs.seas, dims=3, na.rm=T)
  obs.mask <- obs.mask & obs.mn > 0.05 & obs.mn < 0.95
}

## like this, all years are masked out (and therefore the corresponding scores also)
obs.seas[rep(!(obs.mask & fcst.mask), length=length(obs.seas))] <- NA

## close file connections (desperate debugging)
sapply(fc.con, nc_close)
sapply(obs.con, nc_close)
rm(fc.con, obs.con)

## get crossvalidation string
strategy <- list(type=ifelse(length(grep('crossval', method)) == 1, 
                             "crossval", "none"))
if (strategy$type == "crossval"){
  suppressWarnings(strategy$blocklength <- as.numeric(gsub('_.*', '', gsub('.*crossval', '', method))))
}
  
## loop through raw / CCR forecasts
for (is.ccr in c(FALSE, TRUE)){
  
  if (is.ccr){
    fcst.ccr <- array(NA, dim(fcst.seas))
    for (loi in 1:nlon){
      for (lai in 1:nlat){
        if (any(!is.na(obs.seas[,loi,lai,]))){
          fcst.ccr[,loi,lai,,] <- aperm(debias(
            fcst=aperm(fcst.seas[,loi, lai,,], c(1,3,2)),
            obs=obs.seas[,loi,lai,],
            method='ccr', 
            crossval=crossval,
            blocklength=nblock), c(1,3,2))
        } ## end of if on missing values in obs
      } ## end of loop on latitudes
    } ## end of loop on longitudes
    
    fcst.seas <- fcst.ccr
  }
  
  
  ## set up output file to write to
  nc <- nc_open(obsfiles[1])
  nc2 <- nc_open(fcfiles[1])
  lon.nc <- nc2$dim$rlon
  lat.nc <- nc2$dim$rlat
  time.nc <- ncdim_def(name='time', 
                       units='days since 1980-01-01', 
                       vals=as.numeric(as.Date(fc.times[[1]][3:length(fc.times[[1]])]) - as.Date('1980-01-01')),
                       unlim=TRUE, 
                       longname='end date of 3-monthly seasons')
  dims.nc <- list(lon.nc, lat.nc, time.nc)
  ## create variables
  vars.nc <- list()
  vars.nc[['lon']] <- ncvar_def('lon', 'degrees_east', dims.nc[1:2], missval=-1e20, longname='longitude')
  vars.nc[['lat']] <- ncvar_def('lat', 'degrees_north', dims.nc[1:2], missval=-1e20, longname='latitude')
  vars.nc[['rotated_pole']] <- ncvar_def('rotated_pole', units='', dim=list(), missval=NULL, prec='char')
  for (score in scores){
    if (score == 'EnsRocss'){
      for (i in 1:3) {
        rscore <- paste0(score, '.cat', i)
        vars.nc[[rscore]] <- ncvar_def(rscore, '1', dims.nc, missval=-1e20, longname=scorelist[[rscore]])
      }
    } else {
      vars.nc[[score]] <- ncvar_def(score, '1', dims.nc, missval=-1e20, longname=scorelist[[score]])
      if (score %in% scores[grep('pss$', scores)]){
        vars.nc[[paste(score, 'sigma', sep='.')]] <- ncvar_def(paste(score, 'sigma', sep='.'), '1', dims.nc, missval=-1e20, longname=scorelist[[score]])
      }      
    }
  }
  
  ## make sure the output directory exists
  if (!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
  
  ## get additional information on content of files
  fc.name <- sapply(strsplit(fcfiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
  obs.name <- sapply(strsplit(obsfiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
  ofile <- paste0(paste(paste0(index, if(is.ccr) '_CCR' else ''), method, fc.name, 'vs', obs.name, paste(range(names(fc.times)), collapse='-'), sep="_"), '_initmon', initmon, '.nc')
  outfile <- paste(outdir, ofile, sep='/')
  desc <- paste0('Skill in seasonal ', index, ' from bias-corrected (', method, if (is.ccr) ' + CCR' else '', if (crossval) ' cross-validated' else '', ') ', gsub('-', ' ', toupper(fc.name)), ' forecasts verified against ', obs.name, ' for ', paste(range(names(fc.times)), collapse='-'))
  
  ## create netcdf file
  ncout <- nc_create(outfile, vars.nc)
  ## write variables
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
  ## add in a global attribute with description
  ncatt_put(ncout, varid=0, attname='description', attval=desc, prec='text')
  nc_close(ncout)
  
  ## compute scores and write to output file
  for (score in scores){
    ncout <- nc_open(outfile, write=TRUE)
    print(score)
    print(system.time( sfo <- scorefun(fcst.seas, obs.seas, score)))
    if (is.list(sfo)){
      sfo <- lapply(sfo, function(x){
        x[x == -Inf] <- -9999
        x[abs(x) == Inf] <- NA
        return(x)
      })
      if (length(grep('pss$', score)) == 1){
        ncvar_put(ncout, varid=score, vals=aperm(sfo[[1]], c(2,3,1)))
        ncvar_put(ncout, varid=paste(score, 'sigma', sep='.'), vals=aperm(sfo[[2]], c(2,3,1)))    
      } else if (score == 'EnsRocss'){
        for (i in 1:3){
          ncvar_put(ncout, varid=paste0(score, '.cat', i), vals=aperm(sfo[[paste0('cat', i)]], c(2,3,1)))
        }
      }
    } else {
      ## set negative infinity (division by zero) to large negative value
      sfo[sfo == -Inf] <- -9999
      sfo[abs(sfo) == Inf] <- NA
      ## write to output file
      ncvar_put(ncout, varid=score, vals=aperm(sfo, c(2,3,1))) 
    }
    nc_close(ncout)      
    
    ## readjust grid specification for variables
    ## this is very important to happen AFTER the file connection to the NetCDF
    ## is closed from R - otherwise the output is shifted latitudinally
    if (score == 'EnsRocss'){
      for (i in 1:3) system(paste0("ncatted -h -a grid_mapping,", score, ".cat", i, ",c,c,'rotated_pole' ", outfile))            
    } else if (is.list(sfo)){
      system(paste0("ncatted -h -a grid_mapping,", score, ".sigma,c,c,'rotated_pole' ", outfile))      
    } else {
      system(paste0("ncatted -h -a grid_mapping,", score, ",c,c,'rotated_pole' ", outfile))      
    }
 
    ## make sure we free as much memory as possible/reasonable
    rm(sfo, ncout)
    gc()
    
    
  } ## end of loop on scores
} ## end of loop on CCR

q(save='no')
