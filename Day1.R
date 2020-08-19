library(EpiModel)
library(statnetWeb)

EpiModel::epiweb("dcm")

statnetWeb::run_sw()

remotes::install_github("statnet/statnetWeb") #bug fix
