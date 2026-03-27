#!/share/apps/R-3.6.1/bin/Rscript

wd="C:/Users/dave/OneDrive/msvc/countDNAVariants/counts"
ppi=600

library(ggplot2)
library(tidyverse)

options(bitmapType='cairo')

setwd(wd)

TrioVars=data.frame(read.table("SBSVariants.txt",header=FALSE,sep="\t"))
# get variants in desired order
colnames(TrioVars)=c("Seq","Sig")
TrioF3Freqs=TrioVars
TrioMutRawFreqs=data.frame(read.table("modelMutation.F3.toPlot.txt",header=FALSE,sep="\t"))
colnames(TrioMutRawFreqs)=c("Seq","Beta","SE","zed")
TrioF3Freqs$Beta=TrioMutRawFreqs[match(TrioF3Freqs$Seq,TrioMutRawFreqs$Seq),2]
TrioF3Freqs$Frequency=exp(TrioF3Freqs$Beta)/(exp(TrioF3Freqs$Beta)+1)
TrioF3Freqs$rn=seq(1,nrow(TrioF3Freqs))
g=ggplot(TrioF3Freqs,aes(y=Frequency,x=reorder(Sig,-rn))) + theme_bw() +
	geom_bar(stat="identity")+ coord_flip() +
	scale_y_continuous(breaks = seq(0,0.01,by =0.001),minor_breaks=NULL,limits=c(0,0.01)) +
	theme(axis.text.y=element_text(size=6,hjust=0),axis.title=element_text(size=14),axis.title.y=element_blank())
png("F3.modelMutation.Freqs.png",width=6*ppi, height=12*ppi, res=ppi)
print(g)
dev.off()

TrioMutRawBetas=data.frame(read.table("modelMutation.fittedF3BetasWithF1s.txt",header=FALSE,sep="\t"))
colnames(TrioMutRawBetas)=c("Seq","Beta","SE","Z","LL")
TrioMutORs=TrioVars
TrioMutORs$Beta=TrioMutRawBetas[match(TrioMutORs$Seq,TrioMutRawBetas$Seq),2]
TrioMutORs$OR=exp(TrioMutORs$Beta)
TrioMutORs$rn=seq(1,nrow(TrioMutORs))
g=ggplot(TrioMutORs,aes(x=OR,y=reorder(Sig,-rn))) + geom_point()+theme_bw() +
	scale_x_continuous(breaks = seq(0,3.5,by =0.5),minor_breaks=NULL,limits=c(0,3.5)) +
	theme(axis.text.y=element_text(size=6,hjust=0),axis.title=element_text(size=14),axis.title.y=element_blank())
png("F3.modelMutation.ORs.png",width=6*ppi, height=12*ppi, res=ppi)
print(g)
dev.off()

TrioCountsToPlot=data.frame(read.table("SBSVariants.counts.txt",header=TRUE,sep="\t"))
colnames(TrioCountsToPlot)=c("Sig","Count")
TrioCountsToPlot$rn=seq(1,nrow(TrioCountsToPlot))
g=ggplot(TrioCountsToPlot,aes(y=Count,x=reorder(Sig,-rn))) +theme_bw() +
	geom_bar(stat="identity")+ coord_flip() +
	scale_y_continuous(breaks = seq(0,1.5e5,by =5e4),minor_breaks=NULL,limits=c(0,1.5e5)) +
	theme(axis.text.y=element_text(size=6,hjust=0),axis.title=element_text(size=14),axis.title.y=element_blank())
png("F3.rawCounts.png",width=6*ppi, height=12*ppi, res=ppi)
print(g)
dev.off()

Signatures=data.frame(read.table("SBS.counts.0.fitted.txt",header=TRUE,sep="\t"))
SBSNames=colnames(Signatures)[2:8]
Signatures=Signatures[,1:(ncol(Signatures)-2)]
Signatures=pivot_longer(Signatures,!Sig,names_to="SBSnumber",values_to="PredictedCount")
Signatures$rn=seq(1,nrow(Signatures))
g=ggplot(Signatures,aes(fill=SBSnumber,x=reorder(Sig,-rn),y=PredictedCount)) + theme_bw() +
	geom_bar(position = "stack", stat = "identity") + coord_flip() +
	scale_fill_grey(start=0, end=0.8,limits=SBSNames) + # this keeps them in the correct order in the key
	scale_y_continuous(breaks = seq(0,1.5e5,by =5e4),minor_breaks=NULL,limits=c(0,1.5e5)) +
	ylab("Predicted count") + labs(fill="") +
	theme(axis.text.y=element_text(size=6,hjust=0),axis.title=element_text(size=14),axis.title.y=element_blank())
png("F3.countsPredictedFromSignatures.png",width=6*ppi, height=12*ppi, res=ppi)
print(g)
dev.off()

