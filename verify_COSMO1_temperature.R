library(myhelpers)
library(tidyverse)
library(parallel)

rfile <- "/store/msclim/bhendj/verification/Rdata/t2m_COSMO-1.Rdata"
if (!file.exists(rfile)){
  readf <- function(f, varname=NULL){
    nc <- nc_open(f)
    if (is.null(varname)) varname <- names(nc$var)[! names(nc$var) %in% c('lon', 'lat')]
    ntime <- with_tz(nc_time(nc), 'UTC')
    t2m <- ncvar_get(nc, varname)
    lon <- ncvar_get(nc, 'lon')
    lat <- ncvar_get(nc, 'lat')
    as_tibble(
      data.frame(
        lon = rep(round(lon, 4), length=length(t2m)),
        lat = rep(round(lat, 4), length=length(t2m)),
        time = rep(ntime, each=length(lon)),
        lead = rep(1:ncol(t2m) - 1, each=length(lon)),
        init = rep(hour(ntime[1]), length=length(t2m)),
        fcst = c(t2m))
    )
  }
  
  reado <- function(f, varname="t2m"){
    nc <- nc_open(f)
    varnames <- setdiff(names(nc$var), varname)
    ntime <- with_tz(nc_time(nc), 'UTC')
    t2m <- ncvar_get(nc, varname)
    out <- as_tibble(
      data.frame(
        time = rep(ntime, each=nrow(t2m)),
        obs = c(t2m))
    )
    for (vn in varnames) {
      if (vn %in% c('lon', 'lat')){
        out[[vn]] <- round(ncvar_get(nc, vn), 4)
      } else {
        out[[vn]] <- ncvar_get(nc, vn)
      }
    }
    out
  }
  
  
  obs <- reado("/store/msclim/bhendj/ml_ce_zhang/full/observations/meteoswiss_t2m_20151001-20171114.nc")
  
  ffs <- list.files("/store/msclim/bhendj/verification/COSMO-1/", pattern=".nc", recursive=TRUE, full.names=TRUE)

  ## finit <- substr(ffs, 44, 47)
  ## fini <- table(finit)
  ## fini <- fini[fini > 200]
  
  ## ftmp <- list()
  ## for (ini in names(fini)){
  ##   print(ini)
  ##   fcst <- do.call("rbind", mclapply(ffs[finit == ini], readf, mc.cores=20))
  ##   ftmp[[ini]] <- inner_join(obs, fcst) %>%
  ##     mutate(error = fcst - obs - 273.15) %>%
  ##     group_by(nat_abbr, lead, init, lat, lon, height) %>%
  ##       summarise(mse = mean((error - mean(error, na.rm=T))**2, na.rm=T),
  ##                 mse2 = mean(error**2, na.rm=T),
  ##                 mae = mean(abs(error - mean(error, na.rm=T)), na.rm=T),
  ##                 mae2 = mean(abs(error), na.rm=T)) %>%
  ##     mutate(yearmon=ini)
  ## }
  
  ## fcst <- do.call("rbind", ftmp)
  
  ## fcst <- inner_join(obs, do.call("rbind", mclapply(ffs, readf, mc.cores=20)))
  
  fcst <- mclapply(ffs, function(x) inner_join(obs, readf(x)), mc.cores = 20)
  
  ftmp <- plyr::rbind.fill(fcst)
  fcst <- as_tibble(ftmp)
  
  save(fcst, file = rfile)
  
} else {
  load(rfile)
}


ftmp <- fcst %>%
  mutate(month = factor(month(time), ordered=TRUE)) %>%
  group_by(nat_ind, nat_abbr, name, lon, lat, height, lead, init, month) %>%
  mutate(error = if (n() > 28*0.8) fcst - obs - 273.15 else fcst * NA,
         error2 = if (n() > 28*0.8) (error - mean(error)) else error * NA) %>%
  ungroup()


## consistency check with published verification
JJA.rmse <- ftmp %>%
  filter(time >= as.POSIXct("2017-06-01 00:00:00 UTC"), 
         time <= as.POSIXct("2017-08-31 23:00:00 UTC")) %>%
  group_by(nat_abbr, lat, lon, height) %>%
  summarise(rmse = sqrt(mean((fcst - obs - 273.15)**2, na.rm=T)))


ggplot(JJA.rmse, aes(x = lon, y=lat, colour=rmse)) +
  geom_point()






ant <- filter(ftmp, nat_abbr == 'ANT') %>%
  mutate(doy = pmin(365,as.numeric(format(time, '%j'))), 
         initdate = time - lead*3600, 
         hour = hour(time))
ant$fmn <- NA
for (doi in 1:365) {
  ind <- (doi + seq(-7,7) - 1) %% 365 + 1
  amn <- ant %>%
    filter(doy %in% ind) %>%
    group_by(hour) %>%
    summarise(fmn = mean(fcst, na.rm=T))
  ind2 <- which(ant$doy == doi)
  ant$fmn[ind2] <- amn$fmn[match(ant$hour[ind2], amn$hour)]
}

ant %>% filter(lead == 24) %>%
  mutate(fres = fcst - fmn) %>%
  ggplot(aes(x=fres, y=error)) + 
  facet_wrap(~ month) + 
  geom_point(aes(colour=factor(init)))


## model the error by lead-time:
fn <- function(train, pred, type=c("ar1", "lm", "lar1"), mintrain=200){
  if (nrow(train) < mintrain) return(rep(NA, nrow(pred)))
  
  type <- match.arg(type)
  if (type %in% c("lm", "lar1")) flm <- lm(error ~ factor(hour)*time, train)
  if (type == 'ar1'){
    par1 <- train %>%
      group_by(initdate) %>%
      summarise(ar1 = if (n() > 10) pacf(error, plot=FALSE)$acf[1] else NA)
    ar1 <- median(par1$ar1, na.rm=T)
    out <- arrange(train, initdate, time)[nrow(train),]$error * ar1 ^ seq(1, nrow(pred))
  } else if (type == 'lm') {
    out <- predict(flm, newdata=pred)
  } else if (type == 'lar1') {
    train$res <- flm$res
    par1 <- train %>%
      group_by(initdate) %>%
      summarise(ar1 = if (n() > 10) pacf(res, plot=FALSE)$acf[1] else NA)
    ar1 <- median(par1$ar1, na.rm=T)
    out <- predict(flm, newdata=pred) + arrange(train, initdate, time)[nrow(train), ]$res*ar1^seq(1, nrow(pred))
  } 
  return(out)  
}

ant$lm <- ant$ar1 <- ant$lar1 <- NA
for (iinit in unique(ant$initdate)){
  pred <- filter(ant, initdate == iinit)
  print(pred$initdate[1])
  train <- filter(ant, time < iinit, time >= iinit - 24*3600*5, !is.na(error))
  ind <- which(ant$initdate == iinit)
  ant$lm[ind] <- fn(train, pred, 'lm')
  ant$ar1[ind] <- fn(train, pred, 'ar1')
  ant$lar1[ind] <- fn(train, pred, 'lar1')
}
ant$dec <- NA
for (iinit in unique(ant$initdate)){
  ind <- which(ant$initdate == iinit)
  ind2 <- which(ant$initdate == iinit - 24*3600)
  ind3 <- which(ant$initdate == iinit - 2*24*3600)
  if (length(ind3) == length(ind) & length(ind2) == length(ind)){
    ant$dec[ind] <- ant$error[ind2]*0.02 + ant$error[ind3]*0.98
  }
}

aa <- ant %>% group_by(month, lead, init) %>%
  summarize(
    rmse = if (n() >= 10) sqrt(mean((fcst - obs - 273.15)**2, na.rm=T)) else NA,
    lm.rmse = if (n() >= 10) sqrt(mean((fcst - lm - obs - 273.15)**2, na.rm=T)) else NA,
    ar1.rmse = if (n() >= 10) sqrt(mean((fcst - ar1 - obs - 273.15)**2, na.rm=T)) else NA,
    lar1.rmse = if (n() >= 10) sqrt(mean((fcst - lar1 - obs - 273.15)**2, na.rm=T)) else NA,
    dec.rmse = if (n() >= 10) sqrt(mean((fcst - dec - obs - 273.15)**2, na.rm=T)) else NA)


par(mfrow=c(3,4))
for (moni in 1:12) boxplot(filter(aa, month == moni)[,-(1:3)])


ggplot(aa, aes(x = lead + init, colour=factor(init))) + 
  geom_line(aes(y=1 - lar1.rmse / rmse)) +
  facet_wrap( ~ month) + 
  geom_hline(yintercept=0)



indi <- which(ant$lead == 0)
ant$bias2 <- ant$bias <- NA
for (i in indi){
  nowtime <- ant$initdate[i] - 10
  irange <- which(ant$time > nowtime - 60*3600*24 & ant$time < nowtime & ant$init == ant$init[i])
  ifuture <- which(ant$initdate < nowtime + 3600 & ant$time > nowtime & ant$init == ant$init[i])
  if (length(irange) > 1000){
    ilm <- lm(error ~ factor(hour(time))*time, ant[irange,])
    ant[ifuture,'bias'] <- predict(ilm, ant[ifuture,])
    ilm <- lm(error ~ 1, ant[irange,])
    ant[ifuture, 'bias2'] <- predict(ilm, ant[ifuture,])
  }
}











## Fit a series of regresion models
flm <- ftmp %>% 
  group_by(nat_abbr) %>%
  do(broom::glance(lm(error ~ lead + init + cos(as.numeric(month)/6*pi), data = .)))







ftmp %>% group_by(lead, month, init) %>%
  summarise(rmse= mean(error**2, na.rm=T),
            rmse2 = mean(error2**2, na.rm=T)) %>%
  ggplot(aes(x = lead + init, colour=factor(init))) + 
  geom_line(aes(y=1 - rmse2**2 / rmse**2)) +
  facet_wrap( ~ month)









## plor lead time dependency
ftmp %>%
  mutate(below.1000 = height < 1000) %>%
  group_by(lead, below.1000) %>%
  filter(lead <= 33) %>%
  summarise(rmse2 = sqrt(mean(error2**2, na.rm=T)),
            rmse = sqrt(mean(error**2, na.rm=T))) %>%
  ggplot(aes(x = lead, y=rmse, lty=below.1000)) +
  geom_line() + 
  geom_line(aes(y=rmse2), col=2) +  
  theme_bw()



ftmp %>% 
  group_by(init, lead) %>%
  summarise(rmse = sqrt(mean(error2**2, na.rm=T))) %>%
  ggplot(aes(x = lead + init, y=rmse)) + 
  geom_line(aes(colour=factor(init)))


ftmp %>% 
  mutate(altitude = height < 1000) %>%
  group_by(init, lead, altitude) %>%
  summarise(rmse = sqrt(mean(error2 ** 2, na.rm=T))) %>%
  ggplot(aes(x = lead + init, y = rmse)) +
  geom_line(aes(colour = factor(init), lty=altitude))
















ant <- fcst %>%
  filter(nat_abbr == 'ANT')

ant %>%
  mutate(valid_hour = hour(time),
         valid_month = month(time)) %>%
  group_by(valid_hour, valid_month) %>%
  summarise(fcst.mn = mean(fcst - 273.15),
            obs.mn = mean(obs)) %>%
  ggplot(aes(x = valid_hour)) +
  geom_line(aes(y=fcst.mn)) +
  geom_line(aes(y=obs.mn), col=2) +
  facet_wrap(~ valid_month)

ant %>%
  mutate(month = month(time),
         valid_hour = hour(time),
         error = fcst - obs - 273.15) %>%
  group_by(valid_hour, lead, month) %>%
  summarise(rmse = if (n() > 10) sqrt(mean(error**2)) else NA,
            rmse2 = if (n() > 10) sqrt(mean((error - mean(error))**2)) else NA) %>%
  ggplot(aes(x = valid_hour, colour=factor(lead))) +
  geom_point(aes(y=rmse2)) +
  facet_wrap(~ month) +
  ylim(0, NA)

mn <- fcst %>% group_by(nat_abbr, lat, lon, height, init, lead) %>%
  summarise(rmse = sqrt(mean(mse, na.rm=T)),
            rmse2 = sqrt(mean(mse2, na.rm=T)),
            rmae = mean(mae, na.rm=T),
            rmae2 = mean(mae2, na.rm=T))

ggplot(filter(mn, lead == 24), aes(x = rmse, y=height)) +
  geom_line(aes(group = nat_abbr)) + geom_point(aes(colour = factor(init)))

fcst %>% group_by(lat, lon) %>%
  summarise(rmse = sqrt(mean(mse))) %>%
  ggplot(aes(x=lon, y=lat, colour=rmse)) +
  geom_point() +
  theme_bw()

fcst %>% group_by(init, lead, nat_abbr) %>%
  summarise(rmse = sqrt(mean(mse, na.rm=T)),
            rmse2 = sqrt(mean(mse2, na.rm=T))) %>%
  filter(rmse > 0,
         nat_abbr %in% c("SMA", "BUS", "ANT", "JUN")) %>%
  ggplot(aes(x = lead + init, colour=factor(init))) +
  geom_line(aes(y=rmse)) +
  geom_line(aes(y=rmse2), lty=2) +
  facet_wrap(~nat_abbr)

fcst %>% filter(nat_abbr == 'ANT') %>%
  group_by(init, lead) %>%
  summarise(rmse = sqrt(mean(mse, na.rm=T)),
            rmse2 = sqrt(mean(mse2, na.rm=T))) %>%
  ggplot(aes(x = init + lead, y=rmse, colour=factor(init))) +
  geom_line() +
  geom_line(aes(y=rmse2), lty=2)




