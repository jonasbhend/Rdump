## script to plot skill metrics

## load packages
library(methods) ## for processing with Rscript
library(ncdf4) ## ncdf is deprecated but still supported
library(plotmap)
library(colorspace)
library(rworldmap)
## library(rgdal)
library(maptools)

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
      t.out <- as.POSIXct(as.Date(gsub('.*since ', '', t.units), tz='UTC') + nc$dim[[t.i]]$val, tz='UTC')
    }
  } else {
    stop('Time coordinate format not implemented yet')
  }
  return(t.out)
}

addWatermark <- function(txt=NULL, cex=par('cex.axis')*0.5){
  if (is.null(txt)){
    txt <- paste("\uA9 MeteoSwiss,", format(Sys.time(), '%F %H:%M'))
  }
  text(par('usr')[2] - 0.01*diff(par('usr')[1:2]),
       par('usr')[3],
       txt,
       adj=c(1,-0.5),
       cex=cex)
}


## read command line arguments:
infile <- commandArgs(trailingOnly=TRUE)
## check if there are any command line arguments
if (length(infile) != 1){
  infile <- c('/store/msclim/bhendj/EUPORIAS/skill_scores/eobs0.44/seasonal/tas/tas_smooth_1981_2012_ecmwf-system4_vs_E-OBS_1981-2012_initmon11.nc')
}

## output directory
indir <- dirname(infile)
ifile <- basename(infile)
outdir <- gsub('skill_scores', 'figures', indir)
if (!file.exists(outdir)) dir.create(outdir, recursive=TRUE)

## read in grid specification
nc <- nc_open(infile)
rlon <- nc$dim$rlon$vals
rlat <- nc$dim$rlat$vals
plon <- ncatt_get(nc, 'rotated_pole', 'grid_north_pole_longitude')$value
plat <- ncatt_get(nc, 'rotated_pole', 'grid_north_pole_latitude')$value

## read in variable names
scores <- setdiff(names(nc$var), c('lon', 'lat', 'rotated_pole'))
scores <- scores[-grep('sigma', scores)]

## read in forecast times (lead times, are end of months)
fc.times <- as.Date(nc_time(nc))
initmon <- format(fc.times[1] - 80, '%B')
## assemble three-monthly character strings
monnames <- substr(format(as.Date('1990-01-15')+seq(0,330,30), '%b'), 1, 1)
seasnames <- sapply(seq(monnames), 
                    function(x) paste(rep(monnames, 3)[x + 10:12], collapse=''))


## loop through scores and plot
for (score in scores){
  
  for (se in 2){
  #for (se in 1:length(fc.times)){
    seas <- seasnames[as.numeric(format(fc.times[se], '%m'))]    
    
    ## get additional information on score
    longname <- ncatt_get(nc, score, 'long_name')$value
    desc <- ncatt_get(nc, 0, 'description')$value
    plottxt <- paste0(longname, 
                      gsub('Skill in ', ' for ', 
                           gsub('from ', 'from\n', desc)))
    ## chage seasonal to actual season
    plottxt <- gsub('seasonal', seas, plottxt)
    ## add in initialisation
    plottxt <- gsub('forecasts ', paste0('forecasts\n(', initmon, ' initialisation) '), plottxt)
    ## capitalise first letter
    plottxt <- paste0(toupper(substr(plottxt, 1, 1)), substr(plottxt, 2, nchar(plottxt)))
    
    
    ## set up output file for plot
    index <- strsplit(ifile, '_')[[1]][1]
    ifilestem <- paste(strsplit(ifile, '_')[[1]][-1], collapse='_')
    ofile <- paste(score, index, seas, gsub('\\.nc', '.png', ifilestem), sep='_')
    
    ## get the score and add attributes for grid (only one season so far)
    score.grid <- put.grid.atts(ncvar_get(nc, score)[,,se], 
                                grid.type='rotpol',
                                grid.cors=list(rlon, rlat),
                                grid.pars=list(plon=plon, plat=plat))
    ## get stippling
    if (paste0(score, '.sigma') %in% names(nc$var)){
      score.sigma <- score.grid
      score.sigma[,] <- c(score.grid > 1.64*ncvar_get(nc, paste0(score, '.sigma'))[,,se])
    } else {
      score.sigma <- NULL
    }
    
    ## set the level and colours
    if (length(grep('ss$', score)) == 1 | score %in% c('corr', 'EnsCorr') | length(grep('ss.$', score)) == 1){
      ## colours for skill scores from -1 to 1
      breaks <- seq(-0.9, 0.9, 0.1)
      ## colour scale is blue to red
      ## cols <- c(hcl(h=120 + sqrt(seq(140**2, 0, length=10)),
      ##              c=c(seq(10,60,length=4), seq(60,30,length=7)[-1]),
      ##              l=seq(10,98,length=10)),
      ##          rev(heat_hcl(10, l=c(10,98), c=c(90,30))))
      cols <- diverge_hcl(length(breaks) + 1, h=c(260, 10), c=c(90,10), l=c(30, 99), power=c(1/2,1))
    } else if (score %in% c('me', 'EnsMe')) {
      breaks <- pretty(c(score.grid, -score.grid), 14)
      breaks <- breaks[-c(1,length(breaks))]
      if (all(abs(score.grid[!is.na(score.grid)]) < 1e-5)) breaks <- c(-1, -0.1, -0.01, 0.01, 0.1, 1)
      cols <- diverge_hcl(length(breaks)+1, c=c(90,20), l=c(20,95), power=0.7)
    } else if (score %in% c('mae', 'mse', 'rmse', 'EnsMae', 'EnsMse', 'EnsRmse')){
      breaks <- pretty(c(score.grid,0), 10)
      breaks <- breaks[-c(1,length(breaks))]
      cols <- rev(heat_hcl(length(breaks)+1, h=c(-60,120), l=c(10,150), c=c(20,100), power=c(1,1/2)))
    } else if (score %in% c('spr_err', 'EnsSprErr')) {
      breaks <- c(0.33, 0.5, 0.67, 0.8, 0.9, 1.1,1.25,1.5,2,3)
      cols <- diverge_hcl(length(breaks)+1, h=c(260, 10), c=c(90,10), l=c(30,99), power=c(1/2, 1))
    }
    
    lonlim <- c(-10,35)
    latlim <- c(35,70)
    orientation <- c(45, 0, 7.5)
    
    png(paste(outdir, ofile, sep='/'), width=6, height=6, units='in', res=200)
    par(oma=c(0,0,5,0))
    plot.map(score.grid, 
             proj='stereographic', 
             orientation=orientation, 
             add.grid=FALSE, 
             na.col='white', 
             mapdat='none', 
             lonlim=lonlim, 
             latlim=latlim, 
             breaks=breaks,
             col=cols,
             mar.plot=rep(0.5,4), 
             mar.leg=c(0.5, 0.5, 0.5, 3)) 
    ##add.map.outlines(mapdat='worldHires', 
    ##                 lonlim.map=lonlim + c(-30,30),
    ##                 latlim.map=latlim + c(-30,30),
    ##                 col='grey', 
    ##                 lwd=1)
    ###################################################################
    mapoutline <- getMap(resolution='low')
    mm <- sp2tmap(mapoutline)
    mm.proj <- mapproject(mm[,2], mm[,3], projection='stereographic', 
                          orientation=orientation)
    lines(mm.proj, lwd=1, col='grey')
    ## CRS(paste0("+proj=ob_tran +o_proj=longlat ",
    ##                                 "+o_lon_p=", plon, " +o_lat_p=", plat,    
    ##                                 " +ellps=sphere +no_defs")),
    add.map.outlines(mapdat='continents', 
                     lonlim.map=lonlim + c(-30,30),
                     latlim.map=latlim + c(-30,30),
                     col='grey', 
                     lwd=2)
    if (!is.null(score.sigma)){
      xx <- rep(attr(score.sigma, 'rlon'), length(attr(score.sigma, 'rlat')))
      yy <- rep(attr(score.sigma, 'rlat'), each=length(attr(score.sigma, 'rlon')))
      xy.coord <- geocors.trafo(xx, yy, attr(score.sigma, 'grid.type'), attributes(score.sigma)[c('plon', 'plat')], to.type='lonlat', to.pars=list())
      xy <- mapproject(list(x=xy.coord[['lon']], y=xy.coord[['lat']]))
      points(xy$x[score.sigma == 1], xy$y[score.sigma == 1], pch=16, cex=0.3)
    }
    add.map.grid(col='grey')
    addWatermark()
    box()
    axis(3, at=par('usr')[1], plottxt, hadj=0, tick=F)
    dev.off()
    
  } ## end of loop on seasons
} ## end of loop on skill scores

## close the netcdf file
nc_close(nc)
