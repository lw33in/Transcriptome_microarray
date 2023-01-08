# Affymetrix data preprocessing (low-level analysis) 
# Goal: Get low bias, low variance estimates of gene expression from raw data

# Data ploatform: Affymetrix, probe-level data produced by array scanner software (.CEL)

#=====================================================================================
# Directories 
#=====================================================================================
datadir ="/Users/../data/AffyExample/GSE12502/GSE12502_RAW"

#=====================================================================================
# Load packages
#=====================================================================================
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
library("affy")
library(GEOquery) # import data from GEO

#=====================================================================================
# Load data files 
#=====================================================================================
list.files()
list.celfiles()
rawdata=ReadAffy()
rawdata

#=====================================================================================
# Data exploration
#=====================================================================================
length(sampleNames(rawdata))
sampleNames(rawdata)
geneNames(rawdata)[1:5]
length(probeNames(rawdata))

pm(rawdata,"1582358_s_at")[,1:3]

matplot(pm(rawdata,"1582358_s_at"),type="l",xlab="Probe No.",ylab="PM Probe intensity")
matplot(t(pm(rawdata,"1582358_s_at")),type="l",xlab="Array No.",ylab="PM Probe intensity")
hist(rawdata)
boxplot(rawdata)
hist(log2(probes(rawdata,which="pm")[,1]))
hist(log2(probes(rawdata,which="mm")[,1]))

pm(rawdata)[1:4,1:4]
mm(rawdata)[1:4,1:4]

sampleNames(rawdata)

#=====================================================================================
# Run RMA and return log2 transformed data
#=====================================================================================
jrma=justRMA()
exprs(jrma)[1:4,1:4]

gds = getGEO("GDS3538")
eset = GDS2eSet(gds,do.log2=TRUE) # think log2 has been applied, but these data might have been normalized already
anno = pData(eset)

GEOeset = exprs(eset)
2^GEOeset[1:3,1:4]
exprs(jrma)[1:3,1:4]

