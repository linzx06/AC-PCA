library(devtools)
library(roxygen2)

setwd("/Users/zhixianglin/AC-PCA/R_package/")
create("acPCA")

setwd("/Users/zhixianglin/AC-PCA/R_package/acPCA/")
document()
setwd("/Users/zhixianglin/AC-PCA/R_package/")
install("acPCA")
