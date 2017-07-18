## Global parameters
R         <- 100
n         <- 100
nsplits   <- 25


library(devtools)
source.dir <- function(path=getwd(), pattern='\\.[RrSsQq]$', trace = TRUE, ...){
   for (i in list.files(path, pattern=pattern)) {
      if(trace) message("Sourcing:......",i)
      source(file.path(path, i), ...)
   }
   if(trace) cat("\n")
}
source.dir('R_functions')

## Load setups
setups    <- read.table("setups.dat", header=TRUE, stringsAsFactors=FALSE)


## Do cases
for (ncase in 1:nrow(setups)) {
   set.seed(ncase * 19051977)
   source('docase.R')
}
