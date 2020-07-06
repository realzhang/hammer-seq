# usage: Rscript --vanilla calc-global-maint-ratio.R events_count.txt
# zhangzhuqiang@ibp.ac.cn
# last update: 2020.07.05
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please specify the events count file (input file, for instance: events_count.txt)", call.=FALSE)
}
if( "readr" %in% rownames(installed.packages()) ){
    library(readr)
    d <- read_tsv(args[1], col_types = cols())
}else{
    d <- read.table(args[1], head=TRUE)
}

up.maint <- with(d, sum(Up.MM) / ( sum(Up.MM)+sum(Up.MU) ))
lo.maint <- with(d, sum(Lo.MM) / ( sum(Lo.MM)+sum(Lo.MU) ))
maint <- with(d, (sum(Up.MM)+sum(Lo.MM)) / ( sum(Up.MM)+sum(Up.MU)+sum(Lo.MM)+sum(Lo.MU) ) )

up.denovo <- with(d, sum(Up.UM) / ( sum(Up.UU)+sum(Up.UM) ))
lo.denovo <- with(d, sum(Lo.UM) / ( sum(Lo.UU)+sum(Lo.UM) ))
denovo <- with(d, (sum(Up.UM)+sum(Lo.UM)) / ( sum(Up.UU)+sum(Up.UM)+sum(Lo.UU)+sum(Lo.UM) ) )

cat("Maintenance ratio of upper arm:\t", up.maint, "\n")
cat("Maintenance ratio of lower arm:\t", lo.maint, "\n")
cat("Maintenance ratio (combined upper & lower):\t", maint, "\n")

cat("De novo ratio of upper arm:\t", up.denovo, "\n")
cat("De novo ratio of lower arm:\t", lo.denovo, "\n")
cat("De novo ratio (combined upper & lower):\t", denovo, "\n")
