# Data: 48 Illumina RefSeq8 arrays

# Factorial design: 
  # 4 Rx: 2 Control, co-reg X, Y
  # 3 timepoints

# Load packages
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
library("affy")
library(limma)
library("genefilter") # QA plot
# install.packages("latticeExtra")
library("latticeExtra") # QA plot
library("genefilter") #required for rowSds

# Load data files - GenomeStudio Summary files
datadir="/Users/../data/Illumina_RefSeq8_arrays"
setwd(datadir)

sfile="mySample Probe Profile FinalReport.txt" #standard probe profile output from GenomeStudio
cfile="myControl Probe Profile FinalReport.txt" #control probe profile output from GenomeStudio
rawdata = read.ilmn(files=sfile, ctrlfiles=cfile,
                    other.columns=c("Detection","Avg_NBEADS","BEAD_STDERR"))
dim(rawdata$E) #  25209    48
table(rawdata$genes$Status)
targets=read.csv(file="SampleChars.csv",row.names=1) 
head(targets) #contains treatment/time
rawdata$targets=targets[colnames(rawdata),] #store targets w/rawdata

# Study Design
Type = paste(rawdata$targets$treatment,rawdata$targets$time,sep="_")
table(rawdata$targets$treatment,rawdata$targets$time)
table(Type) 

# Boxplot of intensity of sample probes across arrays 
boxplot(log2(rawdata$E[rawdata$genes$Status=="regular",]),range=0,xlab="Arrays",
	ylab="log2 intensities",main="Regular probes")

# make density plot of regular probes
# sample first, then add other samples as lines to plot
plot(density(log2(rawdata$E[rawdata$genes$Status=="regular",1])),
	xlab="log2 intensities", ylab="Density",
	main="Density Plot of Intensities",ylim=c(0,1.3)) 
# add expression across each array
na=length(colnames(rawdata$E))
for (i in 2:na)
	lines(density(log2(rawdata$E[rawdata$genes$Status=="regular",i])),col=i,lty=i)
legend(12,1.2,colnames(rawdata$E),lty=1:5,col=1:6,cex=.9) 

# QA plot
library("genefilter") 
library("latticeExtra") 
dd=dist2(log2(rawdata$E[rawdata$genes$Status=="regular",]))
diag(dd)=0
dd.row=as.dendrogram(hclust(as.dist(dd)))
row.ord=order.dendrogram(dd.row)
legend=list(top=list(fun=dendrogramGrob,
                     args=list(x=dd.row,side="top")))
lp=levelplot(dd[row.ord,row.ord],
             scales=list(x=list(rot=90)),xlab="",
             ylab="",legend=legend)
lp
rm(dd,row.ord,lp)

# Estimate the proportion of genes expressed on each array
library(limma)
proportion=propexpr(rawdata)
names(proportion)=targets$treatment
proportion
# average these over different conditions
tapply(proportion,targets$time,mean)
tapply(proportion,targets$treatment,mean)

# Summarize numbers of control probes
head(unique(rawdata$gene$Status),n=10)
length(unique(rawdata$gene$Status))
sum(rawdata$gene$Status!="regular")
# rawdata$gene$Status[which(substr(rawdata$gene$Status,1,4)=="ERCC")]="ERCC"
table(rawdata$gene$Status)

# make boxplot of negative control probes 
boxplot(log2(rawdata$E[rawdata$genes$Status=="NEGATIVE",]),range=0,xlab="Arrays",
        ylab="log2 intensities",main="Negative control probes")

# Plot the control probes across arrays
# first 664 control probes
sum(rawdata$gene$Status=="NEGATIVE")
n.ncp=sum(rawdata$gene$Status=="NEGATIVE")
l2i=t(log2(rawdata$E[rawdata$genes$Status=="NEGATIVE",]))
dim(l2i) #48 x 664
x=array(rep(1:na,n.ncp),dim=c(na,n.ncp))
dim(x) #48 x 664
matplot(x,l2i,type="l",ylim=c(6,max(log2(rawdata$E))),col=1,lty=1,
        xlab="Array",ylab="log2(raw intensity)",main="Control probes aross arrays",axes=F)
axis(1,1:48,rownames(rawdata$targets))
axis(2)
box()

# Background correction + quantile normalization + log2 transformation
dat = neqc(rawdata)
dim(dat) #24526 x 48 #control probes removed!
# illustrate the effect of log2 transf. on variance
library("genefilter") #required for rowSds
# the "unclass" command removes the class of an object: only remains the underlying type (usually, "list")
batch = unclass(dat$targets$treatment)

# MA plot
# create MA plot of each batch
dat=dat[,order(dat$targets$treatment,dat$targets$time)]
idx=(dat$targets$treatment=="Control2")
A=apply(dat[,idx]$E,1,median)
M=dat[,idx]$E-A
#pairs(dat[,idx]$E,lower.panel=NULL,cex=.5)
plot(A,M[,1],ylab="M",cex=.5)
lines(lowess(A,M[,1]),col=2,lwd=4)
# par(mfrow=c(1,1)) if you want one big plot
par(mfrow=c(2,2))
for (i in 1:4) {
  plot(A,M[,i],ylab="M",main=paste("Batch",i,"vs Overall Avg"),cex=.3)
  lines(lowess(A,M[,i]),col=2,lwd=2)
}
dev.off()
rm(M,A)

# Batch correction
dat$targets$Type=paste(dat$targets$treatment, dat$targets$time, sep="_")
library("sva")
plotMDS(dat,labels=dat$targets$Type,
        col=unclass(as.factor(dat$targets$Type)),
        xlim = c(-1.5,1.5), ylim=c(-1,1)) #color by type

trt=unclass(dat$targets$Type)
mod = model.matrix(~as.factor(trt), data=as.data.frame(trt))
mod0= model.matrix(~1,data=as.data.frame(trt))
# estimate number of surrogate variables
n.sv = num.sv(dat$E,mod,method="leek")
n.sv = num.sv(dat$E,mod,method="be")
# estimate 8 surrogate variables
svobj = sva(dat$E,mod,mod0,n.sv=2)
head(svobj$sv)
dim(svobj$sv) #48 2
modsv=model.matrix(~svobj$sv)
fit=lmFit(dat$E,modsv)
yhat=fit$coef %*% t(modsv)
dat.svresid=dat$E-yhat
plotMDS(dat.svresid,labels=dat$targets$Type,
        col=unclass(as.factor(dat$targets$Type)),xlim=c(-1.5,1.5),
        ylim=c(-1,1.3)) 
# Results after estimating 2 surrogate variables looked better

