library(easyVerification)
library(rbenchmark)
library(parallel)

if (file.exists('~/tmp/benchmark_easyverification.Rdata')){
  load('~/tmp/benchmark_easyverification.Rdata')
} else {
  tm <- toyarray(1)
  btmp <- benchmark(tmroc <- veriApply("EnsRocss", tm$fcst, tm$obs, parallel=TRUE, ncpu=1, prob=1:4/5), replications=10)
  bench <- expand.grid(size=10^(0:3), ncpu=2^(0:5))
  bench <- cbind(bench, array(NA, c(nrow(bench), length(btmp))))
  names(bench)[-(1:2)] <- names(btmp)  
}

## loop through configs
if (any(is.na(bench$elapsed))){
  for (i in which(is.na(bench$elapsed))){
    print(paste("Size: ", bench[i,'size']))
    print(paste("nCPU: ", bench[i,'ncpu']))
    tm <- toyarray(bench[i,'size'])
    bench[i,-(1:2)] <- benchmark(tmroc <- veriApply("EnsRocss", tm$fcst, tm$obs, parallel=TRUE, ncpu=bench[i,'ncpu'], prob=1:4/5), replications=10)
    rm(tm)
    gc()
    save(bench, file="~/tmp/benchmark_easyverification.Rdata")
  }  
}
