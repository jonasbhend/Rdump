#!/usr/bin/env Rscript 

args <- commandArgs(TRUE)
if (length(args) != 2){
  cat(
"
Usage: read_MOSMIX.R <infile> <outfile>
      
This script converts (compressed) MOSMIX fimos files into netcdf. If the 
outfile contains cosmoe, then MOSMIX stations in the COSMO-E domain will 
be processed, else only MOSMIX stations in Switzerland are processed.
"
  )
  q(save = "no")
}



point_in_region <- function(lo, la, region = NULL) {
  if (is.null(region)) return(rep(TRUE, length(lo)))
  if (is.character(region)) region <- sanitizePolygon(maps::map(regions = region, plot = FALSE))
  
  if (all(sapply(region, length) == 2)){
    lo > min(region$x) & lo < max(region$x) & la > min(region$y) & la < max(region$y)
  } else {
    sp::point.in.polygon(lo, la, region$x, region$y) == 1
  }
}



## new strategy to read all of the long-range MOSMIX files
read_mosmix <- function(file, region = NULL) {
  
  rr <- vroom::vroom_lines(file)
  nrows <- 115
  rrind <- seq(1, length(rr), by = nrows + 2)
  rrmeta <- rr[rrind] %>%
    readr::read_table(
      col_names = c("x1", "wmo_ind", "datetime", "x2", "x3", "lat", "lon", "height", "type")
    ) %>%
    dplyr::mutate(lon = deg2dec(lon),
                  lat = deg2dec(lat),
                  rowind = rrind) %>%
    dplyr::filter(point_in_region(lon, lat, region))
  
  rrmeta$data <- lapply(
    rrmeta$rowind,
    function(i) {
      rr_matrix <- rr[seq(i + 1, by = 1, length = nrows)] %>%
        strsplit("\\s+") %>%
        sapply(as.numeric)
      rr_matrix[rr_matrix < -9000] <- NA
      rr_matrix[1,2] <- 0 # set analysis time to zero (hour is present in meta data)
      out <- rr_matrix[-1,] %>%
        tibble::as_tibble() %>%
        stats::setNames(rr_matrix[1,]) %>%
        dplyr::rename(lead = `9997`)
    }
  )
  rrmeta
}








## file <- '/store/msclim/bhendj/tmp/mixfimos_2017063010'
# files <- list.files("/store/msclim/NWP/MOSMIX/burglind", full.names=TRUE)
files <- list.files(
  "/store/msclim/NWP/MOSMIX/hourly_temp",
  "mos_mix_20.*[0-9].gz$",
  full.names = TRUE
)
file <- files[1]
swiss <- list(x=c(5.5, 11)[c(1,2,2,1,1)], y=c(45.5, 48)[c(1,1,2,2,1)])
##cosmo_small <- list(x=c(1.425, 15.825)[c(1,2,2,1,1)], y=c(43.1, 49.44)[c(1,1,2,2,1)])
##cosmo <- list(x=c(0.15, 16.75), y=c(42.65, 49.75))
cosmoe <- list(x = c(5.3, 10.9), y = c(45.49, 48.1))

deg2dec <- function(x){
  x %/% 100 + x %% 100 / 60
}

sanitizePolygon <- function(ll){
  ll <- as.data.frame(ll[c("x", "y")])
  ii <- which(apply(ll, 1, function(x) any(is.na(x))))
  ibnds <- cbind(c(1, ii + 1), c(ii - 1, length(ll$x)))
  ll2 <- lapply(1:nrow(ibnds), function(i) ll[ibnds[i,1]:ibnds[i,2],])
  out <- ll2[[1]]
  ll2 <- ll2[-1]
  for (i in seq_along(ll2)){
    endout <- unlist(out[nrow(out),])
    ll.dist <- lapply(ll2, function(x) colMeans((t(as.matrix(x[c(1, nrow(x)),])) - endout)**2))
    sel.i <- which.min(sapply(ll.dist, min))
    o.append <- ll2[[sel.i]]
    if (diff(ll.dist[[sel.i]]) < 0) o.append <- o.append[nrow(o.append):1,]
    out <- rbind(out, o.append)
    ll2 <- ll2[-sel.i]
  }
  out
}

readHeader <- function(f){
  h1 <- suppressWarnings(readLines(f, n=1))
  if (length(h1) == 0) {
    return(h1)
  } else {
    dd <- as_tibble(t(strsplit(h1, "\\s+")[[1]]))
    names(dd) <- c("tag", "id", "date", "nrow", "ntime", "lat", "lon", "height", "null")
    dd <- dd %>%
      mutate(lon=deg2dec(as.numeric(lon)),
             lat=deg2dec(as.numeric(lat)),
             height = as.numeric(height),
             null = as.numeric(null),
             tag=as.numeric(tag),
             nrow=as.numeric(nrow),
             ntime=as.numeric(ntime))
    return(dd)
  }
}

read_fimosloc <- function(file){
  stopifnot(file.exists(file))
  f <- gzfile(file, 'r')
  on.exit(close(f))
  
  out <- list()
  for (i in 1:1e5){
    h1 <- readHeader(f)
    if (length(h1) == 0) break
    out[[i]] <- h1
    tmp <- suppressWarnings(readLines(f, n=h1[['nrow']][1] + 2))
  }
  out <- Reduce(rbind, out)
  return(out)
}

read_fimos <- function(file, region = 'Switzerland'){
  stopifnot(file.exists(file))
  f <- gzfile(file, "r")
  on.exit(close(f))
  
  if (!is.null(region)){
    if (is.character(region)){
      lola <- sanitizePolygon(maps::map(regions=region, plot=FALSE))
    } else if (is.list(region)){
      lola <- region
    }
  }
  
  out <- list()
  for (i in 1:1e5){
    
    h1 <- readHeader(f)
    if (length(h1) == 0) break
    
    ## check conditions
    if (is.null(region)){
      p.in.poly <- 1
    } else {
      p.in.poly <- sp::point.in.polygon(h1$lon, h1$lat, lola$x, lola$y)
    }
    
    if (p.in.poly > 0) {
      h2 <- scan(f, nlines=1, quiet=TRUE)
      data <- array(scan(f, nlines=h1$nrow[1], quiet=TRUE),
                    h1[1,c('ntime', 'nrow')] + c(2,0))
      data[data < -9000] <- NA
      dd <- as_tibble(data[-(1:2),])
      names(dd) <- data[1,]
      dd <- dd[apply(!is.na(dd), 1, any), ]
      odat <- as_tibble(cbind(h1[,c("id", "date", "lon", "lat", "height")],
                              tibble(lead = h2[-(1:2)]),
                              dd))
      t1 <- suppressWarnings(readLines(f, n=1))
      stopifnot(substr(t1, 1, 4) == "9999")
      
      out[[h1[['id']][1]]] <- odat
    } else {
      tmp <- suppressWarnings(readLines(f, n=h1$nrow[1] + 2))
    }
  }
  out <- Reduce(bind_rows, out)
  out <- as_tibble(out[,apply(!is.na(out), 2, any)])
  return(out)
}




library(maps)
library(ncdf4)
library(rgeos)
library(maptools)

switzerland <- map2SpatialPolygons(sanitizePolygon(maps::map(region = 'Switzerland', plot=FALSE)),
                                   IDs = 'Switzerland')
switzerland_buffer <- gBuffer(switzerland, width=0.044, byid=TRUE)
sp2coord <- function(x) {
  xout <- apply(x@polygons[[1]]@Polygons[[1]]@coords, 2, list)
  names(xout) <- c("x", "y")
  lapply(xout, function(x) x[[1]])
}
switzerland <- sp2coord(switzerland)
switzerland_buffer <- sp2coord(switzerland_buffer)

fimosfile <- '/store/msclim/NWP/MOSMIX/fimos_dict.Rdata'
if (!file.exists(fimosfile)){
  fimos_dict <- read.table("/store/msclim/NWP/MOSMIX/fimos_dict.csv", header=TRUE) %>% as_tibble() %>%
    mutate(name = gsub("_", "", name), unit = gsub("-", "1", unit))
  save(fimos_dict, file = fimosfile)
}




fimos <- function(file, region=NULL){
  ## load the dictionary
  load(fimosfile)
  
  stopifnot(file.exists(file))
  f <- gzfile(file, "r")
  on.exit(close(f))
  
  if (!is.null(region)){
    cfun <- function(lo, la, x, y){
      sp::point.in.polygon(lo, la, x, y)
    }
    if (is.character(region)){
      lola <- sanitizePolygon(maps::map(regions=region, plot=FALSE))
    } else {
      lola <- region
      if (all(sapply(lola, length) == 2)){
        cfun <- function(lo, la, x, y){
          min(diff(sign(x - lo)), diff(sign(y - la)))/2
        }
      }
    }
  }
  
  dd <- hh <- list()
  j <- 0
  for (i in 1:1e5){
    
    h1 <- readHeader(f)
    if (length(h1) == 0) break
    
    ## check lat lon
    if (!is.null(region)){
      p.in.poly <- cfun(h1$lon, h1$lat, lola$x, lola$y)
    } else {
      p.in.poly <- 1
    }
    
    if (p.in.poly > 0) {
      j <- j + 1
      h2 <- scan(f, nlines=1, quiet=TRUE)
      dtmp <- array(scan(f, nlines=h1$nrow[1], quiet=TRUE),
                    h1[1,c('ntime', 'nrow')] + c(2,0))
      data <- dtmp[-1:-2,]
      colnames(data) <- as.character(dtmp[1,])
      rownames(data) <- as.character(h2[-1:-2])
      t1 <- suppressWarnings(readLines(f, n=1))
      
      dd[[j]] <- data
      hh[[j]] <- h1
    } else {
      tmp <- suppressWarnings(readLines(f, n=h1$nrow[1] + 2))
    }
  }
  
  nns <- unique(unlist(sapply(dd, colnames)))
  fimel <- as.character(filter(fimos_dict, element %in% nns)$element)
  oo <- sapply(fimel, function(nn) {
    out <- sapply(dd, function(x) x[,nn], simplify='array')
    out[out < filter(fimos_dict, element == nn)$lower[1]] <- NA
    out[out > filter(fimos_dict, element == nn)$upper[1]] <- NA
    out},
    simplify=FALSE)
  
  return(list(data=oo,
              header = bind_rows(hh),
              dict = filter(fimos_dict, element %in% fimel)))
}



fimos_to_netcdf <- function(infile, outfile, region=NULL){
  ff <- fimos(infile, region=region)
  
  ## fix duplicated variable names
  fn <- names(ff$data)
  names(ff$data)[duplicated(fn)] <- as.character(as.numeric(fn[duplicated(fn)]) + 1000)
  
  fe <- ff$dict$element
  ff$dict$element[duplicated(fe)] <- fe[duplicated(fe)] + 1000
  
  ## check that time is constant across stations
  rr1 <- rownames(ff$data[[1]])
  stopifnot(sapply(ff$data, rownames) == rr1)
  
  tmpfile <- paste(system("echo $SCRATCH/fimos_$$", TRUE),
                   basename(outfile), sep='/')
  dir.create(dirname(tmpfile), showWarnings=FALSE,
             recursive=TRUE, mode='0755')
  
  ## create NetCDF variables
  init <- as.POSIXct(ff$header$date[1], format = '%Y%m%d%H', tz='UTC')
  time.nc <- ncdim_def('leadtime',
                       units = 'hours',
                       vals = as.numeric(rr1),
                       unlim=TRUE)
  ncell.nc <- ncdim_def('station',
                        '',
                        vals = 1:nrow(ff$header),
                        unlim = FALSE,
                        create_dimvar = FALSE)
  nchar.nc <- ncdim_def("nchar",
                        "",
                        vals = 1:max(nchar(ff$header$id)),
                        create_dimvar = FALSE)
  dim.nc <- list(ncell.nc, time.nc)
  
  ## hack for units
  ff$dict <- ff$dict %>%
    mutate(scale = ifelse(grepl('mm/10', unit), 0.1, 1),
           scale = ifelse(grepl('/100m', unit), scale*100, scale),
           scale = ifelse(standard_name == 'air_temperature', scale/10, scale),
           unit = gsub("\\?", "1", unit),
           unit = gsub("mm\\/10", "mm", unit),
           unit = gsub("\\/100m", "m", unit),
           name = gsub(">", "_above_", name),
           name = gsub("<", "_below_", name),
           name = gsub("\\/10.m$", "", name),
           name = gsub("\\/", "_per_", name))
  
  ## set up station variables
  vars.nc <- list(ncvar_def("id", "",
                            list(nchar.nc, ncell.nc),
                            missval=NULL,
                            prec = 'char'),
                  ncvar_def("lon", "degrees_east", ncell.nc, missval=NULL),
                  ncvar_def("lat", "degrees_north", ncell.nc, missval=NULL),
                  ncvar_def("height", "m", ncell.nc, missval=NULL),
                  ncvar_def('time',
                            paste0('hours since ', format(init, "%Y-%m-%d %H:%M:%S")),
                            time.nc,
                            missval=NULL),
                  ncvar_def("reftime",
                            paste("hours since", format(init, '%Y-%m-%d %H:%M:%S')),
                            list(),
                            missval=NULL))
  
  vars.nc <- c(vars.nc, sapply(names(ff$data), function(nn){
    dd <- filter(ff$dict, element == nn)
    ncvar_def(dd$name, dd$unit, dim.nc, missval=-1e20, compression=4)
  }, simplify=FALSE))
  
  ## write data to file
  ncout <- nc_create(tmpfile, vars.nc)
  ncvar_put(ncout, "id", ff$header$id)
  ncvar_put(ncout, "lon", ff$header$lon)
  ncvar_put(ncout, "lat", ff$header$lat)
  ncvar_put(ncout, "height", ff$header$height)
  ncvar_put(ncout, "time", as.numeric(rr1))
  ncvar_put(ncout, "reftime", 0)
  for (id in names(ff$data)){
    fi <- which(ff$dict$element == id)
    sname <- ff$dict[[fi,'standard_name']]
    scale <- ff$dict[[fi, 'scale']]
    ncvar_put(ncout, vars.nc[[id]], t(ff$data[[id]]*scale))
    if (!is.na(sname)) {
      ncatt_put(ncout,
                varid = vars.nc[[id]],
                attname = 'standard_name',
                attval = sname,
                prec = 'text')
    }
    ncatt_put(ncout,
              varid = vars.nc[[id]],
              attname = "coordinates",
              attval = "lon lat height reftime time",
              prec = "text")
  }
  nc_close(ncout)
  
  dir.create(dirname(outfile), showWarnings=FALSE,
             recursive=TRUE, mode = "0755")
  system(paste("mv", tmpfile, outfile))
  file.remove(dirname(tmpfile))
  
}

## execute script with input and output
is.cosmo <- grepl("cosmoe", args[2])
fimos_to_netcdf(args[1], args[2], if (is.cosmo) cosmoe else switzerland_buffer)
q(save = 'no')

