# mGWASqual.txt represent LC-MS profiling data with rows and columns
# representing features and accessions (including all 5 replicates).

trt<-read.table("mGWASqual.txt",header=T,sep="\t",stringsAsFactor=F)
colnam<-as.character(trt[,1])
trt<-trt[,-1]
trt[,][trt[,]==0]<-NA

# generate descriptive statistics for all feature abundances

ave<-apply(trt,1,mean,na.rm=T)
std<-apply(trt,1,sd,na.rm=T)
desc<-data.frame(colnam,ave,std)
colnames(desc)<-as.character(c("trait","average","standard deviation"))
write.table(desc,"descrip.txt",sep="\t",row.names=F)

# Convert feature abundances to qualitative traits.

ec<-read.table("ecotypesALL.txt",header=F,sep="\t")[,1]
trt<-t(trt)
colnames(trt)<-colnam
trt<-data.frame(ec,trt)
trt<-trt[order(trt[,1]),]
trtfin<-array(rep(NA,183*6790),dim=c(183,6790))
i=2
while(i<=dim(trt)[2]){
 # For a particular trait, i.e. column, replace NAs with 1 and values with 0 
 sl<-ifelse(is.na(trt[,i]),1,0)
 # Sum all 1's, i.e. the NAs, per accession and divide by the number of reps
 # If more than half of the accession reps were absent: trait is considered absent
 SumPerAccession<-tapply(sl,as.factor(trt[,1]),sum)	
 AccessionLength<-tapply(sl,as.factor(trt[,1]),length)
 trtfin[,i-1]<-ifelse(SumPerAccession/AccessionLength<=0.5,1,0)
 i=i+1
}
ecUNIQ<-unique(trt[,1])
rownames(trtfin)<-ecUNIQ
colnames(trtfin)<-colnam
# save trtfin as .RData

# load snp, trait and kinship data
# trait values indicate presence/absence

load("snps_sel.RData")		# rows are accessions, columns are snps
load("trtfin.RData")		# rows are accessions, columns are traits
load("kinship_sel.RData")	# rows and columns are accessions

# run trait - snp association via fisher exact tests

trt<-colnames(trtfin)
snp<-colnames(snps_sel)
trt_cha<-as.character()
snp_cha<-as.character()
prob_num<-as.numeric()
odds_cha<-as.character()
trtXkinAVE_num<-as.numeric()
trtXkinSD_num<-as.numeric()

j=1
repeat{
 # skip traits that are either present or absent in all accessions
 if((sum(trtfin[,j])==183)|(sum(trtfin[,j])==0)){
  j=j+1
  next
 }
 i=1
 while(i<=dim(snps_sel)[2]){
  # skip snps that are either present or absent in all accessions
  if((sum(snps_sel[,i])==183)|(sum(snps_sel[,i])==0)){
   i=i+1
   next
  }
  fisht<-fisher.test(table(trtfin[,j],snps_sel[,i]))
  # only retain associations with P < 10E-5
  if(fisht$p.value<0.00001){
   trt_cha<-append(trt_cha,trt[j])
   snp_cha<-append(snp_cha,snp[i])
   prob_num<-append(prob_num,fisht$p.value)
   odds_cha<-append(odds_cha,fisht$estimate)
   # compute kinship average and SD across the cells of the contingency table
   trt0sel<-rownames(trtfin[trtfin[,j]==0,])
   trt1sel<-rownames(trtfin[trtfin[,j]==1,])
   snp0sel<-rownames(snps_sel[snps_sel[,i]==0,])
   snp1sel<-rownames(snps_sel[snps_sel[,i]==1,])
   int00<-intersect(trt0sel,snp0sel)
   int01<-intersect(trt0sel,snp1sel)
   int10<-intersect(trt1sel,snp0sel)
   int11<-intersect(trt1sel,snp1sel)
   kin00<-mean(kinship_sel[int00,int00][lower.tri(kinship_sel[int00,int00])])
   kin01<-mean(kinship_sel[int01,int01][lower.tri(kinship_sel[int01,int01])])
   kin10<-mean(kinship_sel[int10,int10][lower.tri(kinship_sel[int10,int10])])
   kin11<-mean(kinship_sel[int11,int11][lower.tri(kinship_sel[int11,int11])])
   trtXkinAVE<-mean(c(kin00,kin01,kin10,kin11))
   trtXkinSD<-sd(c(kin00,kin01,kin10,kin11))
   trtXkinAVE_num<-append(trtXkinAVE_num,trtXkinAVE)
   trtXkinSD_num<-append(trtXkinSD_num,trtXkinSD)
   l1<-length(snp_cha)
   l2<-length(prob_num)
   l3<-length(odds_cha)
   l4<-length(trtXkinAVE_num)
   l5<-length(trtXkinSD_num)
   if (mean(c(l1,l2,l3,l4,l5))!=length(trt_cha)) print(paste("trt",j,"snp",i,sep=" "))
  }
  i=i+1
 }
 if(j==dim(trtfin)[2])break
 j=j+1
}
associations<-data.frame(trt_cha,snp_cha,prob_num,odds_cha,trtXkinAVE_num,
				trtXkinSD_num)
write.table(associations,"fulldataset.txt",sep="\t",row.names=F)

############################################################################

