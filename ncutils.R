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
      t.out <- as.POSIXct(as.Date(gsub('.*since ', '', t.units), tz='UTC') + round(nc$dim[[t.i]]$val), tz='UTC')
    }
  } else {
    stop('Time coordinate format not implemented yet')
  }
  return(t.out)
}

## function to quickly write data to a NetCDF file
## with the same dimensions as in 
nc_write <- function(nctempfile, file, varname, data, append=FALSE, ...){
  
  nctemplate <- nc_open(nctempfile, readunlim=FALSE, suppress_dimvals=TRUE)
  on.exit(nc_close(nctemplate))
  
  ## check if an accordingly named variable exists
  if (!any(names(nctemplate$var) == varname)) stop('Variable to write does not exist')
  
  ## open netcdf file (and close on exit of function)
  on.exit(nc_close(ncout))
  if (append & file.exists(file)){
    ncout <- nc_open(filename=file, write=TRUE, readunlim=FALSE, suppress_dimvals=TRUE)
  } else {
    if (nctemplate$dim$tim$len == dim(data)[length(dim(data))]){
      system(paste('cdo -s setrtomiss,-1e20,1e20', nctempfile, file))
    } else {
      fdates <- system(paste('cdo -s showdate', nctempfile), TRUE)
      fdates <- gsub("  ", " ", gsub("  ", " ", fdates))
      fdates <- strsplit(fdates, ' ')[[1]]
      fdates <- fdates[nchar(fdates) == 10]
      system(paste('cdo -s -L setrtomiss,-1e20,1e20', paste0('-seldate,',fdates[1], ',', fdates[dim(data)[length(dim(data))]]), nctempfile, file))
    }
    ncout <- nc_open(filename=file, write=TRUE, readunlim=FALSE, suppress_dimvals=TRUE)      
  }
  
  ## write the data
  ncvar_put(ncout, varid=varname, vals=data, ...)
  
}
