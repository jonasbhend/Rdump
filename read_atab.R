library(tidyverse)
library(lubridate)

read_atab <- function(file){
  require(tidyverse)
  
  ## deparse first few lines
  ll <- substr(readLines(file, n=30), 1, 5)
  nn <- which(ll  == "Indic") - 1
  mm <- which(ll == " YYYY") - 1

  ## read header
  head <- sapply(strsplit(readLines(file, n=nn)[-1], ":\\s+"),
              function(x) gsub("^ .", "", x))
  header <- data.frame(head[2,,drop=FALSE])
  names(header) <- head[1,]
  header <- as.data.frame(
              lapply(header, function(x) {
                htmp <- suppressWarnings(as.numeric(x))
                if (any(is.na(htmp))) x else htmp
              }))
  names(header) <- gsub("\\.", "_", names(header))

  ## read in meta information
  meta <- as_tibble(t(
            read.fwf(file, skip=nn, n=mm - nn,
                     widths=c(15 + 10 + header$Width_of_text_label_column,
                       rep(c(-1, 13), header$Number_of_data_columns)))))[-1,]
  mtmp <- c("nat_abbr", "long", "lati", "grid_i", "grid_j", "height")
  names(meta) <- if (ncol(meta) == 4) mtmp[c(1,2,3,6)] else mtmp
  meta <- as_tibble(lapply(meta, function(x) {
    x <- gsub(" *", "", x)
    xtmp <- suppressWarnings(as.numeric(x))
    if (any(is.na(xtmp))) x else xtmp
  }))
       
  ## read data columns
  data <- as_tibble(read.table(file, skip=mm, header=TRUE))
  data <- data %>%
    mutate(time = as.POSIXct(paste(YYYY, MM, DD, hh, sep='-'), format="%Y-%m-%d-%H", tz='UTC')) %>%
      gather(key="nat_abbr", value="value", ABO:ZER)
  
  out <- full_join(meta, data) %>%
    select(-YYYY, -MM, -DD, -hh, -mm)
  out[out == header$Missing] <- NA
  
  return(out)
}


