#Cluster, factoextra, magrittr,NbClust package for generating heat maps
library(cluster)
library(factoextra)
library(magrittr)
library(NbClust)

ds<-read.csv("sampledata.csv", header=TRUE, sep=",")    # ???
library(tidyr)
wsds<- as.data.frame(spread(ds, aminoacid, freq))     #
View(wsds)
write.table(wsds, file="samplefinaldata.csv", sep=",", row.names=FALSE)
ds2<-read.csv("samplefinaldata.csv",header=TRUE, sep=",", row.names=1)

ala<-ds2[,1:4]       # 4    ==== basically the number of codons that codes for a specific amino acids
arg<-ds2[,5:10]      # 6    ====  argininr (R) = 6 (AGA, AGG, CGA, CGC, CGG, CGT)
asn<-ds2[,11:12]     # 2    ==== 18 amino acids (Met and Trp not included)
asp<-ds2[,13:14]     # 2
cys<-ds2[,15:16]     # 2
gln<-ds2[,20:21]     # 2
glu<-ds2[,22:23]     # 2
gly<-ds2[,24:27]     # 4
his<-ds2[,28:29]     # 2
ile<-ds2[,30:32]     # 3
leu<-ds2[,33:38]     # 6
lys<-ds2[,39:40]     # 2
phe<-ds2[,42:43]     # 2
pro<-ds2[,44:47]     # 4
ser<-ds2[,48:53]     # 6
thr<-ds2[,54:57]     # 4
tyr<-ds2[,59:60]     # 2
val<-ds2[,61:64]     # 4

#function for generating heatmaps for the codon usage
heatmap<-function(data) {
  my_data<-data
res.dist<-get_dist(my_data,stand=TRUE,method="spearman")
#Visualize the dynamics of the distances
fviz_dist(res.dist,gradient=list(low="#00AFBB",mid="white",high="#FC4E07"))
#High "Orange"; LOw "Blue
}

heatmap(ala)
heatmap(cys)
heatmap(asn)
heatmap(glu)
heatmap(phe)
heatmap(gly)
heatmap(his)
heatmap(ile)
heatmap(leu)
heatmap(lys)
heatmap(asp)
heatmap(pro)
heatmap(gln)
heatmap(arg)
heatmap(ser)
heatmap(thr)
heatmap(val)
heatmap(tyr)


#Computing for the relative frequencies for alanine
ala$Rfala1<-ala$Ala1/rowSums(ala[,1:4])
ala$rfse1<-sqrt(ala$Rfala1*(1-ala$Rfala1)/rowSums(ala[,1:4]))
ala$Rfala2<-ala$Ala2/rowSums(ala[,1:4])
ala$rfse2<-sqrt(ala$Rfala2*(1-ala$Rfala2)/rowSums(ala[,1:4]))
ala$Rfala3<-ala$Ala3/rowSums(ala[,1:4])
ala$rfse3<-sqrt(ala$Rfala3*(1-ala$Rfala3)/rowSums(ala[,1:4]))
ala$Rfala4<-ala$Ala4/rowSums(ala[,1:4])
ala$rfse4<-sqrt(ala$Rfala4*(1-ala$Rfala4)/rowSums(ala[,1:4]))

#Computing for the relative frequencies for arginine
arg$Rfarg1<-arg$Arg1/rowSums(arg[,1:6])
arg$rfse1<-sqrt(arg$Rfarg1*(1-arg$Rfarg1)/rowSums(arg[,1:6]))
arg$Rfarg2<-arg$Arg2/rowSums(arg[,1:6])
arg$rfse2<-sqrt(arg$Rfarg2*(1-arg$Rfarg2)/rowSums(arg[,1:6]))
arg$Rfarg3<-arg$Arg3/rowSums(arg[,1:6])
arg$rfse3<-sqrt(arg$Rfarg3*(1-arg$Rfarg3)/rowSums(arg[,1:6]))
arg$Rfarg4<-arg$Arg4/rowSums(arg[,1:6])
arg$rfse4<-sqrt(arg$Rfarg4*(1-arg$Rfarg4)/rowSums(arg[,1:6]))
arg$Rfarg5<-arg$Arg5/rowSums(arg[,1:6])
arg$rfse5<-sqrt(arg$Rfarg5*(1-arg$Rfarg5)/rowSums(arg[,1:6]))
arg$Rfarg6<-arg$Arg6/rowSums(arg[,1:6])
arg$rfse6<-sqrt(arg$Rfarg6*(1-arg$Rfarg6)/rowSums(arg[,1:6]))

#Computation of rf for asparagine
asn$Rfasn1<-asn$Asn1/rowSums(asn[,1:2])
asn$rfse1<-sqrt(asn$Rfasn1*(1-asn$Rfasn1)/rowSums(asn[,1:2]))
asn$Rfasn2<-asn$Asn2/rowSums(asn[,1:2])
asn$rfse2<-sqrt(asn$Rfasn2*(1-asn$Rfasn2)/rowSums(asn[,1:2]))

#Computation of rf for aspartate
asp$Rfasp1<-asp$Asp1/rowSums(asp[,1:2])
asp$rfse1<-sqrt(asp$Rfasp1*(1-asp$Rfasp1)/rowSums(asp[,1:2]))
asp$Rfasp2<-asp$Asp2/rowSums(asp[,1:2])
asp$rfse2<-sqrt(asp$Rfasp2*(1-asp$Rfasp2)/rowSums(asp[,1:2]))

#Computation of rf for cysteine
cys$Rfcys1<-cys$Cys1/rowSums(cys[,1:2])
cys$rfse1<-sqrt(cys$Rfcys1*(1-cys$Rfcys1)/rowSums(cys[,1:2]))
cys$Rfcys2<-cys$Cys2/rowSums(cys[,1:2])
cys$rfse2<-sqrt(cys$Rfcys2*(1-cys$Rfcys2)/rowSums(cys[,1:2]))

#Computation of rf for glutamine
gln$Rfgln1<-gln$Gln1/rowSums(gln[,1:2])
gln$rfse1<-sqrt(gln$Rfgln1*(1-gln$Rfgln1)/rowSums(gln[,1:2]))
gln$Rfgln2<-gln$Gln2/rowSums(gln[,1:2])
gln$rfse2<-sqrt(gln$Rfgln2*(1-gln$Rfgln2)/rowSums(gln[,1:2]))


#Computation of rf for glutamate
glu$Rfglu1<-glu$Glu1/rowSums(glu[,1:2])
glu$rfse1<-sqrt(glu$Rfglu1*(1-glu$Rfglu1)/rowSums(glu[,1:2]))
glu$Rfglu2<-glu$Glu2/rowSums(glu[,1:2])
glu$rfse2<-sqrt(glu$Rfglu2*(1-glu$Rfglu2)/rowSums(glu[,1:2]))

#Computing for the relative frequencies for glycine
gly$Rfgly1<-gly$Gly1/rowSums(gly[,1:4])
gly$rfse1<-sqrt(gly$Rfgly1*(1-gly$Rfgly1)/rowSums(gly[,1:4]))
gly$Rfgly2<-gly$Gly2/rowSums(gly[,1:4])
gly$rfse2<-sqrt(gly$Rfgly2*(1-gly$Rfgly2)/rowSums(gly[,1:4]))
gly$Rfgly3<-gly$Gly3/rowSums(gly[,1:4])
gly$rfse3<-sqrt(gly$Rfgly3*(1-gly$Rfgly3)/rowSums(gly[,1:4]))
gly$Rfgly4<-gly$Gly4/rowSums(gly[,1:4])
gly$rfse4<-sqrt(gly$Rfgly4*(1-gly$Rfgly4)/rowSums(gly[,1:4]))

#Computation of rf for histidine
his$Rfhis1<-his$His1/rowSums(his[,1:2])
his$rfse1<-sqrt(his$Rfhis1*(1-his$Rfhis1)/rowSums(his[,1:2]))
his$Rfhis2<-his$His2/rowSums(his[,1:2])
his$rfse2<-sqrt(his$Rfhis2*(1-his$Rfhis2)/rowSums(his[,1:2]))

#Computation of rf for isoleucine
ile$Rfile1<-ile$Ile1/rowSums(ile[,1:3])
ile$rfse1<-sqrt(ile$Rfile1*(1-ile$Rfile1)/rowSums(ile[,1:3]))
ile$Rfile2<-ile$Ile2/rowSums(ile[,1:3])
ile$rfse2<-sqrt(ile$Rfile2*(1-ile$Rfile2)/rowSums(ile[,1:3]))
ile$Rfile3<-ile$Ile3/rowSums(ile[,1:3])
ile$rfse3<-sqrt(ile$Rfile3*(1-ile$Rfile3)/rowSums(ile[,1:3]))

#Computing for the relative frequencies for leucine
leu$Rfleu1<-leu$Leu1/rowSums(leu[,1:6])
leu$rfse1<-sqrt(leu$Rfleu1*(1-leu$Rfleu1)/rowSums(leu[,1:6]))
leu$Rfleu2<-leu$Leu2/rowSums(leu[,1:6])
leu$rfse2<-sqrt(leu$Rfleu2*(1-leu$Rfleu2)/rowSums(leu[,1:6]))
leu$Rfleu3<-leu$Leu3/rowSums(leu[,1:6])
leu$rfse3<-sqrt(leu$Rfleu3*(1-leu$Rfleu3)/rowSums(leu[,1:6]))
leu$Rfleu4<-leu$Leu4/rowSums(leu[,1:6])
leu$rfse4<-sqrt(leu$Rfleu4*(1-leu$Rfleu4)/rowSums(leu[,1:6]))
leu$Rfleu5<-leu$Leu5/rowSums(leu[,1:6])
leu$rfse5<-sqrt(leu$Rfleu5*(1-leu$Rfleu5)/rowSums(leu[,1:6]))
leu$Rfleu6<-leu$Leu6/rowSums(leu[,1:6])
leu$rfse6<-sqrt(leu$Rfleu6*(1-leu$Rfleu6)/rowSums(leu[,1:6]))

#Computation of rf for lysine
lys$Rflys1<-lys$Lys1/rowSums(lys[,1:2])
lys$rfse1<-sqrt(lys$Rflys1*(1-lys$Rflys1)/rowSums(lys[,1:2]))
lys$Rflys2<-lys$Lys2/rowSums(lys[,1:2])
lys$rfse2<-sqrt(lys$Rflys2*(1-lys$Rflys2)/rowSums(lys[,1:2]))

#Computation of rf for phenylalanine
phe$Rfphe1<-phe$Phe1/rowSums(phe[,1:2])
phe$rfse1<-sqrt(phe$Rfphe1*(1-phe$Rfphe1)/rowSums(phe[,1:2]))
phe$Rfphe2<-phe$Phe2/rowSums(phe[,1:2])
phe$rfse2<-sqrt(phe$Rfphe2*(1-phe$Rfphe2)/rowSums(phe[,1:2]))

#Computing for the relative frequencies for proline
pro$Rfpro1<-pro$Pro1/rowSums(pro[,1:4])
pro$rfse1<-sqrt(pro$Rfpro1*(1-pro$Rfpro1)/rowSums(pro[,1:4]))
pro$Rfpro2<-pro$Pro2/rowSums(pro[,1:4])
pro$rfse2<-sqrt(pro$Rfpro2*(1-pro$Rfpro2)/rowSums(pro[,1:4]))
pro$Rfpro3<-pro$Pro3/rowSums(pro[,1:4])
pro$rfse3<-sqrt(pro$Rfpro3*(1-pro$Rfpro3)/rowSums(pro[,1:4]))
pro$Rfpro4<-pro$Pro4/rowSums(pro[,1:4])
pro$rfse4<-sqrt(pro$Rfpro4*(1-pro$Rfpro4)/rowSums(pro[,1:4]))

#Computing for the relative frequencies for serine
ser$Rfser1<-ser$Ser1/rowSums(ser[,1:6])
ser$rfse1<-sqrt(ser$Rfser1*(1-ser$Rfser1)/rowSums(ser[,1:6]))
ser$Rfser2<-ser$Ser2/rowSums(ser[,1:6])
ser$rfse2<-sqrt(ser$Rfser2*(1-ser$Rfser2)/rowSums(ser[,1:6]))
ser$Rfser3<-ser$Ser3/rowSums(ser[,1:6])
ser$rfse3<-sqrt(ser$Rfser3*(1-ser$Rfser3)/rowSums(ser[,1:6]))
ser$Rfser4<-ser$Ser4/rowSums(ser[,1:6])
ser$rfse4<-sqrt(ser$Rfser4*(1-ser$Rfser4)/rowSums(ser[,1:6]))
ser$Rfser5<-ser$Ser5/rowSums(ser[,1:6])
ser$rfse5<-sqrt(ser$Rfser5*(1-ser$Rfser5)/rowSums(ser[,1:6]))
ser$Rfser6<-ser$Ser6/rowSums(ser[,1:6])
ser$rfse6<-sqrt(ser$Rfser6*(1-ser$Rfser6)/rowSums(ser[,1:6]))

#Computing for the relative frequencies for serine
thr$Rfthr1<-thr$Thr1/rowSums(thr[,1:4])
thr$rfse1<-sqrt(thr$Rfthr1*(1-thr$Rfthr1)/rowSums(thr[,1:4]))
thr$Rfthr2<-thr$Thr2/rowSums(thr[,1:4])
thr$rfse2<-sqrt(thr$Rfthr2*(1-thr$Rfthr2)/rowSums(thr[,1:4]))
thr$Rfthr3<-thr$Thr3/rowSums(thr[,1:4])
thr$rfse3<-sqrt(thr$Rfthr3*(1-thr$Rfthr3)/rowSums(thr[,1:4]))
thr$Rfthr4<-thr$Thr4/rowSums(thr[,1:4])
thr$rfse4<-sqrt(thr$Rfthr4*(1-thr$Rfthr4)/rowSums(thr[,1:4]))


#Computation of rf for tyrosine
tyr$Rftyr1<-tyr$Tyr1/rowSums(tyr[,1:2])
tyr$rfse1<-sqrt(tyr$Rftyr1*(1-tyr$Rftyr1)/rowSums(tyr[,1:2]))
tyr$Rftyr2<-tyr$Tyr2/rowSums(tyr[,1:2])
tyr$rfse2<-sqrt(tyr$Rftyr2*(1-tyr$Rftyr2)/rowSums(tyr[,1:2]))

#Computing for the relative frequencies for valine
val$Rfval1<-val$Val1/rowSums(val[,1:4])
val$rfse1<-sqrt(val$Rfval1*(1-val$Rfval1)/rowSums(val[,1:4]))
val$Rfval2<-val$Val2/rowSums(val[,1:4])
val$rfse2<-sqrt(val$Rfval2*(1-val$Rfval2)/rowSums(val[,1:4]))
val$Rfval3<-val$Val3/rowSums(val[,1:4])
val$rfse3<-sqrt(val$Rfval3*(1-val$Rfval3)/rowSums(val[,1:4]))
val$Rfval4<-val$Val4/rowSums(val[,1:4])
val$rfse4<-sqrt(val$Rfval4*(1-val$Rfval4)/rowSums(val[,1:4]))



library(boot)
mean.rf<-function(data, i) {
  d<-data[i,]
  return(mean(d))
}
res<-matrix(1,60)
est<-function(para,dataset) {
  bootest<-boot((na.omit(dataset[, para, drop = FALSE])), statistic=mean.rf,R=1000)
  plot(bootest)
  return(bootest)
}


est("Rfala1",ala)
est("Rfala2",ala)
est("Rfala3",ala)
est("Rfala4",ala)

est("Rfarg1",arg)
est("Rfarg2",arg)
est("Rfarg3",arg)
est("Rfarg4",arg)
est("Rfarg5",arg)
est("Rfarg6",arg)

est("Rfasn1",asn)
est("Rfasn2",asn)

est("Rfasp1",asp)
est("Rfasp2",asp)

est("Rfcys1",cys)
est("Rfcys2",cys)

est("Rfgln1",gln)
est("Rfgln2",gln)

est("Rfglu1",glu)
est("Rfglu2",glu)

est("Rfgly1",gly)
est("Rfgly2",gly)
est("Rfgly3",gly)
est("Rfgly4",gly)

est("Rfhis1",his)
est("Rfhis2",his)

est("Rfile1",ile)
est("Rfile2",ile)
est("Rfile3",ile)

est("Rfleu1",leu)
est("Rfleu2",leu)
est("Rfleu3",leu)
est("Rfleu4",leu)
est("Rfleu5",leu)
est("Rfleu6",leu)

est("Rflys1",lys)
est("Rflys2",lys)

est("Rfphe1",phe)
est("Rfphe2",phe)

est("Rfpro1",pro)
est("Rfpro2",pro)
est("Rfpro3",pro)
est("Rfpro4",pro)

est("Rfser1",ser)
est("Rfser2",ser)
est("Rfser3",ser)
est("Rfser4",ser)
est("Rfser5",ser)
est("Rfser6",ser)

est("Rfthr1",thr)
est("Rfthr2",thr)
est("Rfthr3",thr)
est("Rfthr4",thr)

est("Rftyr1",tyr)
est("Rftyr2",tyr)

est("Rfval1",val)
est("Rfval2",val)
est("Rfval3",val)
est("Rfval4",val)

#computing for relative synonymous codon usage
meanval<-apply(ala[,1:4],1,mean,na.rm=TRUE)
ala$scala1<-ala$Ala1/meanval
ala$scala2<-ala$Ala2/meanval
ala$scala3<-ala$Ala3/meanval
ala$scala4<-ala$Ala4/meanval

meanval<-apply(arg[,1:6],1,mean,na.rm=TRUE)
arg$scarg1<-arg$Arg1/meanval
arg$scarg2<-arg$Arg2/meanval
arg$scarg3<-arg$Arg3/meanval
arg$scarg4<-arg$Arg4/meanval
arg$scarg5<-arg$Arg5/meanval
arg$scarg6<-arg$Arg6/meanval

meanval<-apply(asn[,1:2],1,mean,na.rm=TRUE)
asn$scasn1<-asn$Asn1/meanval
asn$scasn2<-asn$Asn2/meanval

meanval<-apply(asp[,1:2],1,mean,na.rm=TRUE)
asp$scasp1<-asp$Asp1/meanval
asp$scasp2<-asp$Asp2/meanval

meanval<-apply(cys[,1:2],1,mean,na.rm=TRUE)
cys$sccys1<-cys$Cys1/meanval
cys$sccys2<-cys$Cys2/meanval

meanval<-apply(gln[,1:2],1,mean,na.rm=TRUE)
gln$scgln1<-gln$Gln1/meanval
gln$scgln2<-gln$Gln2/meanval

meanval<-apply(glu[,1:2],1,mean,na.rm=TRUE)
glu$scglu1<-glu$Glu1/meanval
glu$scglu2<-glu$Glu2/meanval

meanval<-apply(gly[,1:4],1,mean,na.rm=TRUE)
gly$scgly1<-gly$Gly1/meanval
gly$scgly2<-gly$Gly2/meanval
gly$scgly3<-gly$Gly3/meanval
gly$scgly4<-gly$Gly4/meanval

meanval<-apply(his[,1:2],1,mean,na.rm=TRUE)
his$schis1<-his$His1/meanval
his$schis2<-his$His2/meanval

meanval<-apply(ile[,1:3],1,mean,na.rm=TRUE)
ile$scile1<-ile$Ile1/meanval
ile$scile2<-ile$Ile2/meanval
ile$scile3<-ile$Ile3/meanval

meanval<-apply(leu[,1:6],1,mean,na.rm=TRUE)
leu$scleu1<-leu$Leu1/meanval
leu$scleu2<-leu$Leu2/meanval
leu$scleu3<-leu$Leu3/meanval
leu$scleu4<-leu$Leu4/meanval
leu$scleu5<-leu$Leu5/meanval
leu$scleu6<-leu$Leu6/meanval

meanval<-apply(lys[,1:2],1,mean,na.rm=TRUE)
lys$sclys1<-lys$Lys1/meanval
lys$sclys2<-lys$Lys2/meanval

meanval<-apply(phe[,1:2],1,mean,na.rm=TRUE)
phe$scphe1<-phe$Phe1/meanval
phe$scphe2<-phe$Phe2/meanval

meanval<-apply(pro[,1:4],1,mean,na.rm=TRUE)
pro$scpro1<-pro$Pro1/meanval
pro$scpro2<-pro$Pro2/meanval
pro$scpro3<-pro$Pro3/meanval
pro$scpro4<-pro$Pro4/meanval

meanval<-apply(ser[,1:6],1,mean,na.rm=TRUE)
ser$scser1<-ser$Ser1/meanval
ser$scser2<-ser$Ser2/meanval
ser$scser3<-ser$Ser3/meanval
ser$scser4<-ser$Ser4/meanval
ser$scser5<-ser$Ser5/meanval
ser$scser6<-ser$Ser6/meanval

meanval<-apply(thr[,1:4],1,mean,na.rm=TRUE)
thr$scthr1<-thr$Thr1/meanval
thr$scthr2<-thr$Thr2/meanval
thr$scthr3<-thr$Thr3/meanval
thr$scthr4<-thr$Thr4/meanval

meanval<-apply(tyr[,1:2],1,mean,na.rm=TRUE)
tyr$sctyr1<-tyr$Tyr1/meanval
tyr$sctyr2<-tyr$Tyr2/meanval

meanval<-apply(val[,1:4],1,mean,na.rm=TRUE)
val$scval1<-val$Val1/meanval
val$scval2<-val$Val2/meanval
val$scval3<-val$Val3/meanval
val$scval4<-val$Val4/meanval


est("scala1",ala)
est("scala2",ala)
est("scala3",ala)
est("scala4",ala)

est("scarg1",arg)
est("scarg2",arg)
est("scarg3",arg)
est("scarg4",arg)
est("scarg5",arg)
est("scarg6",arg)

est("scasn1",asn)
est("scasn2",asn)

est("scasp1",asp)
est("scasp2",asp)

est("sccys1",cys)
est("sccys2",cys)

est("scgln1",gln)
est("scgln2",gln)

est("scglu1",glu)
est("scglu2",glu)

est("scgly1",gly)
est("scgly2",gly)
est("scgly3",gly)
est("scgly4",gly)

est("schis1",his)
est("schis2",his)

est("scile1",ile)
est("scile2",ile)
est("scile3",ile)

est("scleu1",leu)
est("scleu2",leu)
est("scleu3",leu)
est("scleu4",leu)
est("scleu5",leu)
est("scleu6",leu)

est("sclys1",lys)
est("sclys2",lys)

est("scphe1",phe)
est("scphe2",phe)

est("scpro1",pro)
est("scpro2",pro)
est("scpro3",pro)
est("scpro4",pro)

est("scser1",ser)
est("scser2",ser)
est("scser3",ser)
est("scser4",ser)
est("scser5",ser)
est("scser6",ser)

est("scthr1",thr)
est("scthr2",thr)
est("scthr3",thr)
est("scthr4",thr)

est("sctyr1",tyr)
est("sctyr2",tyr)

est("scval1",val)
est("scval2",val)
est("scval3",val)
est("scval4",val)

#computing for CAI weight
maxval<-apply(ala[,1:4],1,max,na.rm=TRUE)
ala$wala1<-ala$Ala1/maxval
ala$wala2<-ala$Ala2/maxval
ala$wala3<-ala$Ala3/maxval
ala$wala4<-ala$Ala4/maxval

est("wala1",ala)
est("wala2",ala)
est("wala3",ala)
est("wala4",ala)

maxval<-apply(arg[,1:6],1,max,na.rm=TRUE)
arg$warg1<-arg$Arg1/maxval
arg$warg2<-arg$Arg2/maxval
arg$warg3<-arg$Arg3/maxval
arg$warg4<-arg$Arg4/maxval
arg$warg5<-arg$Arg5/maxval
arg$warg6<-arg$Arg6/maxval

est("warg1",arg)
est("warg2",arg)
est("warg3",arg)
est("warg4",arg)
est("warg5",arg)
est("warg6",arg)

maxval<-apply(asn[,1:2],1,max,na.rm=TRUE)
asn$wasn1<-asn$Asn1/maxval
asn$wasn2<-asn$Asn2/maxval

est("wasn1",asn)
est("wasn2",asn)

maxval<-apply(asp[,1:2],1,max,na.rm=TRUE)
asp$wasp1<-asp$Asp1/maxval
asp$wasp2<-asp$Asp2/maxval

est("wasp1",asp)
est("wasp2",asp)

maxval<-apply(cys[,1:2],1,max,na.rm=TRUE)
cys$wcys1<-cys$Cys1/maxval
cys$wcys2<-cys$Cys2/maxval

est("wcys1",cys)
est("wcys2",cys)


maxval<-apply(gln[,1:2],1,max,na.rm=TRUE)
gln$wgln1<-gln$Gln1/maxval
gln$wgln2<-gln$Gln2/maxval

est("wgln1",gln)
est("wgln2",gln)

maxval<-apply(glu[,1:2],1,max,na.rm=TRUE)
glu$wglu1<-glu$Glu1/maxval
glu$wglu2<-glu$Glu2/maxval

est("wglu1",glu)
est("wglu2",glu)

maxval<-apply(gly[,1:4],1,max,na.rm=TRUE)
gly$wgly1<-gly$Gly1/maxval
gly$wgly2<-gly$Gly2/maxval
gly$wgly3<-gly$Gly3/maxval
gly$wgly4<-gly$Gly4/maxval

est("wgly1",gly)
est("wgly2",gly)
est("wgly3",gly)
est("wgly4",gly)

maxval<-apply(his[,1:2],1,max,na.rm=TRUE)
his$whis1<-his$His1/maxval
his$whis2<-his$His2/maxval

est("whis1",his)
est("whis2",his)

maxval<-apply(ile[,1:3],1,max,na.rm=TRUE)
ile$wile1<-ile$Ile1/maxval
ile$wile2<-ile$Ile2/maxval
ile$wile3<-ile$Ile3/maxval

est("wile1",ile)
est("wile2",ile)
est("wile3",ile)

maxval<-apply(leu[,1:6],1,max,na.rm=TRUE)
leu$wleu1<-leu$Leu1/maxval
leu$wleu2<-leu$Leu2/maxval
leu$wleu3<-leu$Leu3/maxval
leu$wleu4<-leu$Leu4/maxval
leu$wleu5<-leu$Leu5/maxval
leu$wleu6<-leu$Leu6/maxval

est("wleu1",leu)
est("wleu2",leu)
est("wleu3",leu)
est("wleu4",leu)
est("wleu5",leu)
est("wleu6",leu)

maxval<-apply(lys[,1:2],1,max,na.rm=TRUE)
lys$wlys1<-lys$Lys1/maxval
lys$wlys2<-lys$Lys2/maxval

est("wlys1",lys)
est("wlys2",lys)

maxval<-apply(phe[,1:2],1,max,na.rm=TRUE)
phe$wphe1<-phe$Phe1/maxval
phe$wphe2<-phe$Phe2/maxval

est("wphe1",phe)
est("wphe2",phe)

maxval<-apply(pro[,1:4],1,max,na.rm=TRUE)
pro$wpro1<-pro$Pro1/maxval
pro$wpro2<-pro$Pro2/maxval
pro$wpro3<-pro$Pro3/maxval
pro$wpro4<-pro$Pro4/maxval

est("wpro1",pro)
est("wpro2",pro)
est("wpro3",pro)
est("wpro4",pro)

maxval<-apply(ser[,1:6],1,max,na.rm=TRUE)
ser$wser1<-ser$Ser1/maxval
ser$wser2<-ser$Ser2/maxval
ser$wser3<-ser$Ser3/maxval
ser$wser4<-ser$Ser4/maxval
ser$wser5<-ser$Ser5/maxval
ser$wser6<-ser$Ser6/maxval

est("wser1",ser)
est("wser2",ser)
est("wser3",ser)
est("wser4",ser)
est("wser5",ser)
est("wser6",ser)

maxval<-apply(thr[,1:4],1,max,na.rm=TRUE)
thr$wthr1<-thr$Thr1/maxval
thr$wthr2<-thr$Thr2/maxval
thr$wthr3<-thr$Thr3/maxval
thr$wthr4<-thr$Thr4/maxval

est("wthr1",thr)
est("wthr2",thr)
est("wthr3",thr)
est("wthr4",thr)


maxval<-apply(tyr[,1:2],1,max,na.rm=TRUE)
tyr$wtyr1<-tyr$Tyr1/maxval
tyr$wtyr2<-tyr$Tyr2/maxval

est("wtyr1",tyr)
est("wtyr2",tyr)

maxval<-apply(val[,1:4],1,max,na.rm=TRUE)
val$wval1<-val$Val1/maxval
val$wval2<-val$Val2/maxval
val$wval3<-val$Val3/maxval
val$wval4<-val$Val4/maxval

est("wval1",val)
est("wval2",val)
est("wval3",val)
est("wval4",val)

bootest<-boot((na.omit(ala[, "Ala1", drop = FALSE])), sim="parametric", mle=mean(ala$Ala1),R=1000)



#Maximum Likelihood estimation of parameters of the MLE based on the observed data
mle<-function(x) {
  experiments <- rmultinom(n=1, size=1000, prob=x)
  df=data.frame(experiments)/1000
  x<-df[,1]
  se<-sqrt((x*(1-x))/1000)
  return( cbind(x,se))
}

ala1<-sum(ala$Ala1)/sum(ala[,1:4])
ala2<-sum(ala$Ala2)/sum(ala[,1:4])
ala3<-sum(ala$Ala3)/sum(ala[,1:4])
ala4<-sum(ala$Ala4)/sum(ala[,1:4])
mle(c(ala1,ala2,ala3,ala4))

arg1<-sum(arg$Arg1)/sum(arg[,1:6])
arg2<-sum(arg$Arg2)/sum(arg[,1:6])
arg3<-sum(arg$Arg3)/sum(arg[,1:6])
arg4<-sum(arg$Arg4)/sum(arg[,1:6])
arg5<-sum(arg$Arg5)/sum(arg[,1:6])
arg6<-sum(arg$Arg6)/sum(arg[,1:6])
mle(c(arg1,arg2,arg3,arg4, arg5, arg6))

asn1<-sum(asn$Asn1)/sum(asn[,1:2])
asn2<-sum(asn$Asn2)/sum(asn[,1:2])
mle(c(asn1,asn2))

asp1<-sum(asp$Asp1)/sum(asp[,1:2])
asp2<-sum(asp$Asp2)/sum(asp[,1:2])
mle(c(asp1,asp2))

cys1<-sum(cys$Cys1)/sum(cys[,1:2])
cys2<-sum(cys$Cys2)/sum(cys[,1:2])
mle(c(cys1,cys2))

gln1<-sum(gln$Gln1)/sum(gln[,1:2])
gln2<-sum(gln$Gln2)/sum(gln[,1:2])
mle(c(gln1,gln2))

glu1<-sum(glu$Glu1)/sum(glu[,1:2])
glu2<-sum(glu$Glu2)/sum(glu[,1:2])
mle(c(glu1,glu2))


gly1<-sum(gly$Gly1)/sum(gly[,1:4])
gly2<-sum(gly$Gly2)/sum(gly[,1:4])
gly3<-sum(gly$Gly3)/sum(gly[,1:4])
gly4<-sum(gly$Gly4)/sum(gly[,1:4])
mle(c(gly1,gly2,gly3,gly4))

his1<-sum(his$His1)/sum(his[,1:2])
his2<-sum(his$His2)/sum(his[,1:2])
mle(c(his1,his2))

ile1<-sum(ile$Ile1)/sum(ile[,1:3])
ile2<-sum(ile$Ile2)/sum(ile[,1:3])
ile3<-sum(ile$Ile3)/sum(ile[,1:3])
mle(c(ile1,ile2,ile3))


leu1<-sum(leu$Leu1)/sum(leu[,1:6])
leu2<-sum(leu$Leu2)/sum(leu[,1:6])
leu3<-sum(leu$Leu3)/sum(leu[,1:6])
leu4<-sum(leu$Leu4)/sum(leu[,1:6])
leu5<-sum(leu$Leu5)/sum(leu[,1:6])
leu6<-sum(leu$Leu6)/sum(leu[,1:6])
mle(c(leu1,leu2,leu3,leu4, leu5, leu6))

lys1<-sum(lys$Lys1)/sum(lys[,1:2])
lys2<-sum(lys$Lys2)/sum(lys[,1:2])
mle(c(lys1,lys2))

phe1<-sum(phe$Phe1)/sum(phe[,1:2])
phe2<-sum(phe$Phe2)/sum(phe[,1:2])
mle(c(phe1,phe2))

pro1<-sum(pro$Pro1)/sum(pro[,1:4])
pro2<-sum(pro$Pro2)/sum(pro[,1:4])
pro3<-sum(pro$Pro3)/sum(pro[,1:4])
pro4<-sum(pro$Pro4)/sum(pro[,1:4])
mle(c(pro1,pro2,pro3,pro4))

ser1<-sum(ser$Ser1)/sum(ser[,1:6])
ser2<-sum(ser$Ser2)/sum(ser[,1:6])
ser3<-sum(ser$Ser3)/sum(ser[,1:6])
ser4<-sum(ser$Ser4)/sum(ser[,1:6])
ser5<-sum(ser$Ser5)/sum(ser[,1:6])
ser6<-sum(ser$Ser6)/sum(ser[,1:6])
mle(c(ser1,ser2,ser3,ser4, ser5, ser6))

thr1<-sum(thr$Thr1)/sum(thr[,1:4])
thr2<-sum(thr$Thr2)/sum(thr[,1:4])
thr3<-sum(thr$Thr3)/sum(thr[,1:4])
thr4<-sum(thr$Thr4)/sum(thr[,1:4])
mle(c(thr1,thr2,thr3,thr4))

tyr1<-sum(tyr$Tyr1)/sum(tyr[,1:2])
tyr2<-sum(tyr$Tyr2)/sum(tyr[,1:2])
mle(c(tyr1,tyr2))

val1<-sum(val$Val1)/sum(val[,1:4])
val2<-sum(val$Val2)/sum(val[,1:4])
val3<-sum(val$Val3)/sum(val[,1:4])
val4<-sum(val$Val4)/sum(val[,1:4])
mle(c(val1,val2,val3,val4))
