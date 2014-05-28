############################
# R make-like file for data gathering
# OKUM GAS and MUH-1 certification project
# updated 20140519
# Thomas Meisel
############################

# Set working directory
setwd("C:/Daten/projects/R/certification/")

# Gather and clean up csv data files which have been prepared in EXCEL
source("GOMGather1.R")

# Merging data
source("GOMMerge.R")
