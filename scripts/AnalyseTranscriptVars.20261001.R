#!/share/apps/R-3.6.1/bin/Rscript
library(survival)
library(ggplot2)

wd="/cluster/project9/bipolargenomes/UKBB/WGSvars"

LoadFile="transcriptVarLoad.txt"
LoadFileNonCoding="transcriptVarLoad.template.txt" 
AgeFile="/home/rejudcu/UKBB/cancer.20260928/UKBB.Age.txt"
SexFile="/home/rejudcu/UKBB/UKBB.sex.20230807.txt"
CancerFile="/home/rejudcu/UKBB/cancer.20260928/UKBB.Cancerall.txt"
PCsFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt"
EthFile="/home/rejudcu/UKBB/UKBB.ethnicity.20241011.txt"
DatesFile="/home/rejudcu/UKBB/cancer.20260928/UKBB.Dates.txt"
AllDataFile="AllData.20260930.txt"

setwd(wd)

for (chr in c(1:22,"X")) {
	ChrVars=data.frame(read.table(sprintf("transcriptVars.template.%s.raw",chr),header=TRUE,sep="",stringsAsFactors=FALSE))
	ChrVars=ChrVars[,-c(1,3:6)]
	ChrVars[is.na(ChrVars)]=0 
	ChrVars=ChrVars[,colSums(ChrVars)>1000000 | colSums(ChrVars)==1]
	print(head(ChrVars))
	ChrVars$Load=rowSums(ChrVars[,-1])
	ChrLoad=ChrVars[,c("IID","Load")]
	colnames(ChrLoad)[2]=sprintf("Load.%s",chr)
	if (chr == 1) {
		Load=ChrLoad
	} else {
		Load=merge(Load,ChrLoad,by="IID")
	}
	print(head(Load))
}

Load=Load[Load$IID>1,]
Load$TotalLoad=rowSums(Load[,-1])
print(summary(Load$TotalLoad))
print(mean(Load$TotalLoad))
print(var(Load$TotalLoad))

write.table(Load,LoadFileNonCoding,row.names=FALSE,quote=FALSE,sep="\t")

for (chr in c(1:22,"X")) {
	ChrVars=data.frame(read.table(sprintf("transcriptVars.%s.raw",chr),header=TRUE,sep="",stringsAsFactors=FALSE))
	ChrVars=ChrVars[,-c(1,3:6)]
	ChrVars[is.na(ChrVars)]=0 
	ChrVars=ChrVars[,colSums(ChrVars)>1000000 | colSums(ChrVars)==1]
	print(head(ChrVars))
	ChrVars$Load=rowSums(ChrVars[,-1])
	ChrLoad=ChrVars[,c("IID","Load")]
	colnames(ChrLoad)[2]=sprintf("Load.%s",chr)
	if (chr == 1) {
		Load=ChrLoad
	} else {
		Load=merge(Load,ChrLoad,by="IID")
	}
	print(head(Load))
}

Load=Load[Load$IID>1,]
Load$TotalLoad=rowSums(Load[,-1])
print(summary(Load$TotalLoad))
print(mean(Load$TotalLoad))
print(var(Load$TotalLoad))

write.table(Load,LoadFile,row.names=FALSE,quote=FALSE,sep="\t")

Load=data.frame(read.table(LoadFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
AllData=Load[,c("IID","TotalLoad")]
LoadNonCoding=data.frame(read.table(LoadFileNonCoding,header=TRUE,sep="\t",stringsAsFactors=FALSE))
LoadNonCoding=LoadNonCoding[,c("IID","TotalLoad")]
colnames(LoadNonCoding)=c("IID","TotalLoadNonCoding")
AllData=merge(AllData,LoadNonCoding,by="IID")
Cancer=data.frame(read.table(CancerFile,header=TRUE,sep="",stringsAsFactors=FALSE))
AllData=merge(AllData,Cancer,by="IID")
Age=data.frame(read.table(AgeFile,header=TRUE,sep="",stringsAsFactors=FALSE))
AllData=merge(AllData,Age,by="IID")
Sex=data.frame(read.table(SexFile,header=TRUE,sep="",stringsAsFactors=FALSE))
colnames(Sex)=c("IID","Sex")
AllData=merge(AllData,Sex,by="IID")
PCs=data.frame(read.table(PCsFile,header=TRUE,sep="",stringsAsFactors=FALSE))
AllData=merge(AllData,PCs,by="IID")
Dates=data.frame(read.table(DatesFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))
AllData=merge(AllData,Dates,by="IID")
Eth=data.frame(read.table(EthFile,header=TRUE,sep="",stringsAsFactors=FALSE,fill=TRUE))
AllData=merge(AllData,Eth,by="IID")

write.table(AllData,AllDataFile,row.names=FALSE,quote=FALSE,sep="\t")

AllData=data.frame(read.table(AllDataFile,header=TRUE,sep="\t",stringsAsFactors=FALSE))

covars=""
for (p in 1:20) {
  covars=sprintf("%sPC%d + ",covars,p)
}
covars=sprintf("%s Sex",covars)


formulaString=sprintf("Cancer ~ TotalLoad")
model=glm(as.formula(formulaString), data = AllData, family="binomial")
print(summary(model))

formulaString=sprintf("Cancer ~ TotalLoad + Age")
model=glm(as.formula(formulaString), data = AllData, family="binomial")
print(summary(model))

formulaString=sprintf("Cancer ~ %s + Age + TotalLoad",covars)
model=glm(as.formula(formulaString), data = AllData, family="binomial")
print(summary(model))

formulaString=sprintf("Age ~ %s + TotalLoad",covars)
model=glm(as.formula(formulaString), data = AllData)
print(summary(model))

WhiteCodes=c(1,1001,1002,1003)
White=AllData[AllData$ethnicity.20241011 %in% WhiteCodes,]

formulaString=sprintf("Age ~ %s + TotalLoad",covars)
model=glm(as.formula(formulaString), data = White)
print(summary(model))

formulaString=sprintf("Age ~ %s + TotalLoadNonCoding",covars)
model=glm(as.formula(formulaString), data = White)
print(summary(model))

formulaString=sprintf("Age ~ %s + TotalLoad + TotalLoadNonCoding",covars)
model=glm(as.formula(formulaString), data = White)
print(summary(model))

print(nrow(White))

print(mean(White$TimeToEvent[White$Status==0]))


f=survfit(Surv(TimeToEvent,Status)~1,data=AllData[AllData$TotalLoad<3,])
ff=survfit(Surv(TimeToEvent,Status)~1,data=AllData[AllData$TotalLoad>=3,])
summary(f)
summary(ff)

c=coxph(Surv(TimeToEvent,Status)~TotalLoad,data=White)
summary(c)
cc=coxph(Surv(TimeToEvent,Status)~TotalLoad+Age+Sex+PC1+PC2+PC3+PC7+PC14+PC15+PC19,data=White)
summary(cc)

options(bitmapType='cairo')
png(height=600, width=600, pointsize=25, file="TotalLoad.png")
p=ggplot(AllData, aes(as.factor(TotalLoad), Age))
p=p+geom_boxplot() + labs(x = "Variant load", y = "Age") + theme_minimal()
print(p)
dev.off()

N=nrow(AllData)
l=sum(AllData$TotalLoad)/N
f=c(1,1,2,6,24)
ResultsTable=data.frame(matrix(NA, ncol = 6, nrow = 4))
colnames(ResultsTable)=c("Variant load","0","1","2","3","4")
ResultsTable[,1]=c("N",sprintf("Expected N (lamda = %.3f)",l),"Cancer prevalence (95% CI)","Age mean (sd)")
for (vv in 0:4) {
	n=nrow(AllData[AllData$TotalLoad==vv,])
	ResultsTable[1,vv+2]=sprintf("%d",n)
	k=vv
	e=N*(l^k)*exp(-l)/f[vv+1]
	ResultsTable[2,vv+2]=sprintf("%.2f",e)
	c=nrow(AllData[AllData$TotalLoad==vv & AllData$Cancer==1,])/nrow(AllData[AllData$TotalLoad==vv,])
	se=sqrt(c*(1-c)/n)
	ResultsTable[3,vv+2]=sprintf("%.3f (%.3f - %.3f)",c,c-1.96*se,c+1.96*se)
	m=mean(AllData$Age[AllData$TotalLoad==vv])
	s=sd(AllData$Age[AllData$TotalLoad==vv])
	ResultsTable[4,vv+2]=sprintf("%.2f (%.2f)",m,s)
}
print(ResultsTable)
e=as.numeric(ResultsTable[2,5])+as.numeric(ResultsTable[2,6])
o=as.numeric(ResultsTable[1,5])+as.numeric(ResultsTable[1,6])
ch=(o-e)*(o-e)/e
print(sprintf("Chi-squared = %.2f, 1 df, p = %g",ch,pchisq(ch,df=1,lower.tail=FALSE)))
write.table(ResultsTable,"ResultsTable.txt",row.names=FALSE,quote=FALSE,sep="\t")



