# Original R analysis script used in the melanoma genomic meta-analysis
# Note: local setwd() lines from the author machine have been commented out for portability.
# The study-level source dataset provided in this repository is Table S11.
# Some file names referenced in this historical script correspond to intermediate working CSVs from the original workflow.

library(tidyverse)
library(meta)
library(metafor)
##setwd("C:/Users/buioa/Dropbox/These/France/Paper Duc/statistical manuscript")
# setwd("~/Dropbox/These/France/Paper Duc/statistical manuscript")

#mut<- read.csv("mutationprimaryver2.csv",sep = ';')
mut<- read.csv("mutationprimaryver287.csv",sep = ';')
dim(mut)
names(mut)
summary(mut)
table(mut$Method)
mut$Fixe<-as.factor(mut$Fixe)
mut$Fixe<-relevel(mut$Fixe,ref='Frozen')
mut_fixe<-subset(mut,is.na(Fixe)==FALSE)
mut_multi<-subset(mut,is.na(Mutilple)==FALSE)
dim(mut_fixe)
dim(mut_multi)
##mutation BRAF##
#mut1<-subset(mut, is.na(Method)==FALSE)
#mut1_braf<-subset(mut1, is.na(BRAF)==FALSE)
mut_braf<-subset(mut, is.na(BRAF)==FALSE)
dim(mut_braf)
names(mut_braf)
m2<- metaprop(BRAF, Sample, studlab=Author, sm="PFT", data=mut_braf, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)









##method subgroup
mut1<-subset(mut, is.na(Method)==FALSE)
mut1_braf<-subset(mut1, is.na(BRAF)==FALSE)
m1<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT", data=mut1_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m1, Method)
mu1<-update(m1,subgroup=Method)
forest(mu1)


##ALK#
mut_alk<-subset(mut, is.na(ALK)==FALSE)
dim(mut_alk)
names(mut_alk)
m2<- metaprop(ALK, Sample, studlab=TT, sm="PFT", data=mut_alk, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##ABL1#
mut_abl1<-subset(mut, is.na(ABL1)==FALSE)
dim(mut_abl1)
names(mut_abl1)
m2<- metaprop(ABL1, Sample, studlab=Author, sm="PFT", data=mut_abl1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##AR1D2##
mut_ARID2 <-subset(mut, is.na(ARID2 )==FALSE)
dim(mut_ARID2 )
names(mut_ARID2 )
m2<- metaprop(ARID2 , Sample, studlab=TT, sm="PFT", data=mut_ARID2 , method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##BRAF V600E##
mut_BRAFV600E<-subset(mut, is.na(BRAF.V600E)==FALSE)
dim(mut_BRAFV600E)
names(mut_BRAFV600E)
m2<- metaprop(BRAF.V600E, Sample, studlab=TT, sm="PFT", data=mut_BRAFV600E, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##CDKN2A##
mut_CDKN2A<-subset(mut, is.na(CDKN2A)==FALSE)
dim(mut_CDKN2A)
names(mut_CDKN2A)
m2<- metaprop(CDKN2A, Sample, studlab=TT, sm="PFT", data=mut_CDKN2A, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##CTNNB1#
mut_CTNNB1<-subset(mut, is.na(CTNNB1)==FALSE)
dim(mut_CTNNB1)
names(mut_CTNNB1)
m2<- metaprop(CTNNB1, Sample, studlab=TT, sm="PFT", data=mut_CTNNB1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##GNAQ##
mut_GNAQ<-subset(mut, is.na(GNAQ)==FALSE)
dim(mut_GNAQ)
names(mut_GNAQ)
m2<- metaprop(GNAQ, Sample, studlab=TT, sm="PFT", data=mut_GNAQ, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##HRAS##
mut_HRAS<-subset(mut, is.na(HRAS)==FALSE)
dim(mut_HRAS)
names(mut_HRAS)
m2<- metaprop(HRAS, Sample, studlab=TT, sm="PFT", data=mut_HRAS, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##KIT##
mut_KIT<-subset(mut, is.na(KIT)==FALSE)
dim(mut_KIT)
names(mut_KIT)
m2<- metaprop(KIT, Sample, studlab=TT, sm="PFT", data=mut_KIT, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##KRAS##
mut_KRAS<-subset(mut, is.na(KRAS)==FALSE)
dim(mut_KRAS)
names(mut_KRAS)
m2<- metaprop(KRAS, Sample, studlab=TT, sm="PFT", data=mut_KRAS, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MAP2K1##
mut_MAP2K1<-subset(mut, is.na(MAP2K1)==FALSE)
dim(mut_MAP2K1)
names(mut_MAP2K1)
m2<- metaprop(MAP2K1, Sample, studlab=TT, sm="PFT", data=mut_MAP2K1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MET##
mut_MET<-subset(mut, is.na(MET)==FALSE)
dim(mut_MET)
names(mut_MET)
m2<- metaprop(MET, Sample, studlab=TT, sm="PFT", data=mut_MET, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NF1##
mut_NF1<-subset(mut, is.na(NF1)==FALSE)
dim(mut_NF1)
names(mut_NF1)
m2<- metaprop(NF1, Sample, studlab=TT, sm="PFT", data=mut_NF1, method="Inverse", NF1hod.tau="DL",is.na=TRUE)
forest(m2)
##NRAS##
mut_NRAS<-subset(mut, is.na(NRAS)==FALSE)
dim(mut_NRAS)
names(mut_NRAS)
m2<- metaprop(NRAS, Sample, studlab=TT, sm="PFT", data=mut_NRAS, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NRASQ61K
mut_NRASQ61K<-subset(mut, is.na(NRAS.Q61K)==FALSE)
dim(mut_NRASQ61K)
names(mut_NRASQ61K)
m2<- metaprop(NRAS.Q61K, Sample, studlab=TT, sm="PFT", data=mut_NRASQ61K, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NRAS Q61R
mut_NRASQ61R<-subset(mut, is.na(NRAS.Q61R)==FALSE)
dim(mut_NRASQ61R)
names(mut_NRASQ61R)
m2<- metaprop(NRAS.Q61R, Sample, studlab=TT, sm="PFT", data=mut_NRASQ61R, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##PIK3CA##
mut_PIK3CA<-subset(mut, is.na(PIK3CA)==FALSE)
dim(mut_PIK3CA)
names(mut_PIK3CA)
m2<- metaprop(PIK3CA, Sample, studlab=TT, sm="PFT", data=mut_PIK3CA, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##PTEN##
mut_PTEN<-subset(mut, is.na(PTEN)==FALSE)
dim(mut_PTEN)
names(mut_PTEN)
m2<- metaprop(PTEN, Sample, studlab=TT, sm="PFT", data=mut_PTEN, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##TERT##
mut_TERT<-subset(mut, is.na(TERT)==FALSE)
dim(mut_TERT)
names(mut_TERT)
m2<- metaprop(TERT, Sample, studlab=TT, sm="PFT", data=mut_TERT, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##TP53##
mut_TP53<-subset(mut, is.na(TP53)==FALSE)
dim(mut_TP53)
names(mut_TP53)
m2<- metaprop(TP53, Sample, studlab=TT, sm="PFT", data=mut_TP53, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##ROS1##
mut_ros1<-subset(mut, is.na(ROS1)==FALSE)
dim(mut_ros1)
names(mut_ros1)
m2<- metaprop(ROS1, Sample, studlab=TT, sm="PFT", data=mut_ros1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


###TABLE2 subgroup analysis##
##method subgroup
mut1<-subset(mut, is.na(Method)==FALSE)
mut1_braf<-subset(mut1, is.na(BRAF)==FALSE)
m1<- metaprop(BRAF, Sample,
              studlab=Author, sm="PFT", data=mut1_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m1, Method)
mu1<-update(m1,subgroup=Method)
forest(mu1)




##preservation subgroup
mut_fixe_braf<-subset(mut_fixe, is.na(BRAF)==FALSE)
dim(mut_fixe_braf)
m2<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu2<-update(m2,subgroup=Fixe,fixed = FALSE)
metareg(m2, Fixe)
forest(mu2)
# Do meta-regression with two covariates
metareg(m2, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_braf<-subset(mut_multi, is.na(BRAF)==FALSE)
dim(mut_multi_braf)
m3<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu3<-update(m3, subgroup=Mutilple,fixed = FALSE)
metareg(m3,Mutilple)
forest(mu3)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_braf<-subset(mut4, is.na(BRAF)==FALSE)
dim(mut4_braf)
m4<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT", subgroup = Quality, data=mut4_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)

###NRAS###
#########
mut4<-subset(mut, is.na(Method)==FALSE)
mut4_NRAS<-subset(mut4, is.na(NRAS)==FALSE)
m4<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT", data=mut4_NRAS, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Method)
mu4<-update(m4,subgroup=Method)
forest(mu4)
##preservation subgroup
mut_fixe_NRAS<-subset(mut_fixe, is.na(NRAS)==FALSE)
dim(mut_fixe_NRAS)
m5<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_NRAS, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu5<-update(m5,subgroup=Fixe,fixed = FALSE)
metareg(m5, Fixe)
forest(mu5)
# Do meta-regression with two covariates
metareg(m5, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_NRAS<-subset(mut_multi, is.na(NRAS)==FALSE)
dim(mut_multi_NRAS)
m6<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_NRAS, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu6<-update(m6, subgroup=Mutilple,fixed = FALSE)
metareg(m6,Mutilple)
forest(mu6)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_nras<-subset(mut4, is.na(NRAS)==FALSE)
m4<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT", data=mut4_nras, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)


###TP53###
#########
mut7<-subset(mut, is.na(Method)==FALSE)
mut7_TP53<-subset(mut7, is.na(TP53)==FALSE)
m7<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT", data=mut7_TP53, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m7, Method)
mu7<-update(m7,subgroup=Method)
forest(mu7)
##preservation subgroup
mut_fixe_TP53<-subset(mut_fixe, is.na(TP53)==FALSE)
dim(mut_fixe_TP53)
m5<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_TP53, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu5<-update(m5,subgroup=Fixe,fixed = FALSE)
metareg(m5, Fixe)
forest(mu5)
# Do meta-regression with two covariates
metareg(m5, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_TP53<-subset(mut_multi, is.na(TP53)==FALSE)
dim(mut_multi_TP53)
m9<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_TP53, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu9<-update(m9, subgroup=Mutilple,fixed = FALSE)
metareg(m9,Mutilple)
forest(mu9)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_tp53<-subset(mut4, is.na(TP53)==FALSE)
m4<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT", data=mut4_tp53, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)


###TERT###
#########
mut11<-subset(mut, is.na(Method)==FALSE)
mut11_TERT<-subset(mut11, is.na(TERT)==FALSE)
m11<- metaprop(TERT, Sample,
               studlab=TT, sm="PFT", data=mut11_TERT, method="Inverse",
               method.tau="DL",is.na=TRUE)
metareg(m11, Method)
mu11<-update(m11,subgroup=Method)
forest(mu11)
##preservation subgroup
mut_fixe_TERT<-subset(mut_fixe, is.na(TERT)==FALSE)
dim(mut_fixe_TERT)
m10<- metaprop(TERT, Sample,
               studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_TERT, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu10<-update(m10,subgroup=Fixe,fixed = FALSE)
metareg(m10, Fixe)
forest(mu10)
# Do meta-regression with two covariates
metareg(m10, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_TERT<-subset(mut_multi, is.na(TERT)==FALSE)
dim(mut_multi_TERT)
m12<- metaprop(TERT, Sample,
               studlab=Author, sm="PFT",subgroup = Mutilple, data=mut_multi_TERT, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu12<-update(m12, subgroup=Mutilple,fixed = FALSE)
metareg(m12,Mutilple)
forest(mu12)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_TERT<-subset(mut4, is.na(TERT)==FALSE)
m4<- metaprop(TERT, Sample,
              studlab=TT, sm="PFT", data=mut4_TERT, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)


###NF1###
#########
mut13<-subset(mut, is.na(Method)==FALSE)
mut13_NF1<-subset(mut13, is.na(NF1)==FALSE)
m13<- metaprop(NF1, Sample,
               studlab=TT, sm="PFT", data=mut13_NF1, method="Inverse",
               method.tau="DL",is.na=TRUE)
metareg(m13, Method)
mu13<-update(m13,subgroup=Method)
forest(mu13)
##preservation subgroup
mut_fixe_NF1<-subset(mut_fixe, is.na(NF1)==FALSE)
dim(mut_fixe_NF1)
m14<- metaprop(NF1, Sample,
               studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_NF1, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu14<-update(m14,subgroup=Fixe,fixed = FALSE)
metareg(m14, Fixe)
forest(mu14)
# Do meta-regression with two covariates
metareg(m14, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_NF1<-subset(mut_multi, is.na(NF1)==FALSE)
dim(mut_multi_NF1)
m15<- metaprop(NF1, Sample,
               studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_NF1, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu15<-update(m15, subgroup=Mutilple,fixed = FALSE)
metareg(m15,Mutilple)
forest(mu15)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_NF1<-subset(mut4, is.na(NF1)==FALSE)
m4<- metaprop(NF1, Sample,
              studlab=TT, sm="PFT", data=mut4_NF1, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)

##TABLE 3#####
#########METASTATICS####

##meta<- read.csv("Bsmutationmetastatic.csv",sep = ';')
meta<- read.csv("mutationmetastatic287.csv",sep = ';')
meta_braf<-subset(meta, is.na(BRAF)==FALSE)
dim(meta_braf)
names(meta_braf)
meta1<- metaprop(BRAF, Sample, studlab=TT, sm="PFT", data=meta_braf, method="Inverse", method.tau="DL",is.na=TRUE)
forest(meta1)
names(meta)

##BRAF V600E##
meta_BRAFV600E<-subset(meta, is.na(BRAF.V600E)==FALSE)
dim(meta_BRAFV600E)
names(meta_BRAFV600E)
m2<- metaprop(BRAF.V600E, Sample, studlab=TT, sm="PFT", data=meta_BRAFV600E, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##CDKN2A##
meta_CDKN2A<-subset(meta, is.na(CDKN2A)==FALSE)
dim(meta_CDKN2A)
names(meta_CDKN2A)
m2<- metaprop(CDKN2A, Sample, studlab=TT, sm="PFT", data=meta_CDKN2A, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##CTNNB1#
meta_CTNNB1<-subset(meta, is.na(CTNNB1)==FALSE)
dim(meta_CTNNB1)
names(meta_CTNNB1)
m2<- metaprop(CTNNB1, Sample, studlab=TT, sm="PFT", data=meta_CTNNB1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##ERBB4##
meta_ERBB4<-subset(meta, is.na(ERBB4)==FALSE)
dim(meta_ERBB4)
names(meta_ERBB4)
m2<- metaprop(ERBB4, Sample, studlab=?..TT, sm="PFT", data=meta_ERBB4, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
####GNA11
meta_GNA11<-subset(meta, is.na(GNA11)==FALSE)
dim(meta_GNA11)
names(meta_GNA11)
m2<- metaprop(GNA11, Sample, studlab=TT, sm="PFT", data=meta_GNA11, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##KIT##
meta_KIT<-subset(meta, is.na(KIT)==FALSE)
dim(meta_KIT)
names(meta_KIT)
m2<- metaprop(KIT, Sample, studlab=TT, sm="PFT", data=meta_KIT, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##KRAS##
meta_KRAS<-subset(meta, is.na(KRAS)==FALSE)
dim(meta_KRAS)
names(meta_KRAS)
m2<- metaprop(KRAS, Sample, studlab=?..TT, sm="PFT", data=meta_KRAS, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MAP2K1##
meta_MAP2K1<-subset(meta, is.na(MAP2K1)==FALSE)
dim(meta_MAP2K1)
names(meta_MAP2K1)
m2<- metaprop(MAP2K1, Sample, studlab=?..TT, sm="PFT", data=meta_MAP2K1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MET##
meta_MET<-subset(meta, is.na(MET)==FALSE)
dim(meta_MET)
names(meta_MET)
m2<- metaprop(MET, Sample, studlab=TT, sm="PFT", data=meta_MET, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NF1##
meta_NF1<-subset(meta, is.na(NF1)==FALSE)
dim(meta_NF1)
names(meta_NF1)
m2<- metaprop(NF1, Sample, studlab=TT, sm="PFT", data=meta_NF1, method="Inverse", NF1hod.tau="DL",is.na=TRUE)
forest(m2)
##NRAS##
meta_NRAS<-subset(meta, is.na(NRAS)==FALSE)
dim(meta_NRAS)
names(meta_NRAS)
m2<- metaprop(NRAS, Sample, studlab=TT, sm="PFT", data=meta_NRAS, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NRASQ61K
meta_NRASQ61K<-subset(meta, is.na(Q61K)==FALSE)
dim(meta_NRASQ61K)
names(meta_NRASQ61K)
m2<- metaprop(Q61K, Sample, studlab=TT, sm="PFT", data=meta_NRASQ61K, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NRAS Q61R
meta_NRASQ61R<-subset(meta, is.na(Q61R)==FALSE)
dim(meta_NRASQ61R)
names(meta_NRASQ61R)
m2<- metaprop(Q61R, Sample, studlab=TT, sm="PFT", data=meta_NRASQ61R, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##PIK3CA##
meta_PIK3CA<-subset(meta, is.na(PIK3CA)==FALSE)
dim(meta_PIK3CA)
names(meta_PIK3CA)
m2<- metaprop(PIK3CA, Sample, studlab=?..TT, sm="PFT", data=meta_PIK3CA, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##PTEN##
meta_PTEN<-subset(meta, is.na(PTEN)==FALSE)
dim(meta_PTEN)
names(meta_PTEN)
m2<- metaprop(PTEN, Sample, studlab=TT, sm="PFT", data=meta_PTEN, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##TERT##
meta_TERT<-subset(meta, is.na(TERT)==FALSE)
dim(meta_TERT)
names(meta_TERT)
m2<- metaprop(TERT, Sample, studlab=TT, sm="PFT", data=meta_TERT, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##TP53##
meta_TP53<-subset(meta, is.na(TP53)==FALSE)
dim(meta_TP53)
names(meta_TP53)
m2<- metaprop(TP53, Sample, studlab=TT, sm="PFT", data=meta_TP53, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##ROS1##
meta_ros1<-subset(meta, is.na(ROS1)==FALSE)
dim(meta_ros1)
names(meta_ros1)
m2<- metaprop(ROS1, Sample, studlab=TT, sm="PFT", data=meta_ros1, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##RB1##
meta_rb1<-subset(meta, is.na(RB1)==FALSE)
dim(meta_rb1)
names(meta_rb1)
m2<- metaprop(RB1, Sample, studlab=TT, sm="PFT", data=meta_rb1, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##SF3B1##
meta_SF3B1<-subset(meta, is.na(SF3B1)==FALSE)
dim(meta_SF3B1)
names(meta_SF3B1)
m2<- metaprop(SF3B1, Sample, studlab=?..TT, sm="PFT", data=meta_SF3B1, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


###################CNA#####################################
###########################################################





#####################2.9#####################
###############################################

scna<- read.csv("CNA primary tumor.csv",sep = ';')
scna<- read.csv("CNA primary tumor plus 30.csv",sep = ';')
dim(scna)
names(scna)
head(scna)
library(ggplot2)
theme_update(text = element_text(size=8))
scna$factor<-ifelse(scna$Prevalence<0, "LOSS", "GAIN")
scna<-scna[order(scna$Prevalence),]
scna$CNA <- factor(scna$CNA, levels = scna$CNA)  # convert to factor to retain sorted order in plot.


##graphique##
ggplot(scna, aes(x=CNA, y=Prevalence, label=Prevalence)) + 
  geom_bar(stat='identity', aes(fill=factor), width=.4)  +
  scale_fill_manual(name="SCNA in primary tumors", 
                    labels = c("GAIN", "LOSS"), 
                    values = c("LOSS"="#00ba38", "GAIN"="#f8766d")) + 
  labs(subtitle="SCNA in primary tumor", 
       title= "SCNA in primary tumor") + 
  coord_flip()


###metastases###

cnameta<- read.csv("CNA metastase.csv",sep = ';')
cnameta<- read.csv("CNA metastases plus 30.csv",sep = ';')
dim(cnameta)
names(cnameta)
head(cnameta)
library(ggplot2)
theme_update(text = element_text(size=8))
cnameta$factor<-ifelse(cnameta$Prevalence<0, "LOSS", "GAIN")
cnameta<-cnameta[order(cnameta$Prevalence),]
cnameta$CNA <- factor(cnameta$CNA, levels = cnameta$CNA)  # convert to factor to retain sorted order in plot.


##graphique##
ggplot(cnameta, aes(x=CNA, y=Prevalence, label=Prevalence)) + 
  geom_bar(stat='identity', aes(fill=factor), width=.4)  +
  scale_fill_manual(name="SCNA in metastases", 
                    labels = c("GAIN", "LOSS"), 
                    values = c("LOSS"="#00ba38", "GAIN"="#f8766d")) + 
  labs(subtitle="SCNA in metastases", 
       title= "SCNA in metastases") + 
  coord_flip()


######Test primary tumor and metastase###
########################################
braf <- prop.test(x = c(3454,1799), n = c(8582,3768))
braf 

braf <- prop.test(x = c(3462,1806), n = c(8648,3807))
braf 



nras <- prop.test(x = c(1156,485), n = c(7217,2365))
nras 

brafv600e<- prop.test(x = c(1593,743), n = c(4415,1767))
brafv600e 

tert<- prop.test(x = c(693,146), n = c(1698,284))
tert

Tp53<- prop.test(x = c(275,86), n = c(1505,577))
Tp53

nf1<-prop.test(x = c(65,56), n = c(432,405))
nf1

kit<-prop.test(x = c(166,24), n = c(2585,818))
kit

pik3ca<-prop.test(x = c(88,4), n = c(1911,172))
pik3ca

cdkn2a<-prop.test(x = c(140,93), n = c(1809,764))
cdkn2a

pten<-prop.test(x = c(66,45), n = c(1323,491))
pten

nrasq61r<-prop.test(x = c(159,46), n = c(2152,544))
nrasq61r

nrasq61k<-prop.test(x = c(76,27), n = c(1812,612))
nrasq61k

ros1<-prop.test(x = c(22,37), n = c(262,193))
ros1






###metastatic site###
##skin###
braf <- prop.test(x = c(3454,55), n = c(8582,113))
braf ##0.08
braf <- prop.test(x = c(3454,55), n = c(8648,113))
braf ##0.08

brafv600e<- prop.test(x = c(1593,33), n = c(4415,78))
brafv600e ##0.3

nras <- prop.test(x = c(1156,14), n = c(7217,74))
nras ##0.6

cdkn2a<-prop.test(x = c(140,4), n = c(1809,38))
cdkn2a ##0.74

nrasq61k<-prop.test(x = c(76,5), n = c(1812,38))
nrasq61k #0.02

nrasq61r<-prop.test(x = c(159,4), n = c(2152,74))
nrasq61r#0.67

pten<-prop.test(x = c(66,2), n = c(1323,18))
pten #0.5


###lymphode##
braf <- prop.test(x = c(3454,298), n = c(8582,522))
braf ##<0.0001

brafv600e<- prop.test(x = c(1593,262), n = c(4415,502))
brafv600e ##<0.0001

nras <- prop.test(x = c(1156,81), n = c(7217,434))
nras ##0.16

cdkn2a<-prop.test(x = c(140,23), n = c(1809,176))
cdkn2a ##0.02

nrasq61k<-prop.test(x = c(76,9), n = c(1812,229))
nrasq61k #0.98

nrasq61r<-prop.test(x = c(159,25), n = c(2152,222))
nrasq61r#0.05

pten<-prop.test(x = c(66,8), n = c(1323,93))
pten #0.2

Tp53<- prop.test(x = c(275,13), n = c(1505,124))
Tp53 #0.04

##central nervous##
braf <- prop.test(x = c(3454,134), n = c(8582,276))
braf ##0.007

brafv600e<- prop.test(x = c(1593,78), n = c(4415,221))
brafv600e ##0.8

nras <- prop.test(x = c(1156,71), n = c(7217,271))
nras ##<0.0001

cdkn2a<-prop.test(x = c(140,8), n = c(1809,50))
cdkn2a ##0.06

nrasq61k<-prop.test(x = c(76,10), n = c(1812,129))
nrasq61k #0.09

nrasq61r<-prop.test(x = c(159,11), n = c(2152,143))
nrasq61r#1

pten<-prop.test(x = c(66,5), n = c(1323,56))
pten #0.3

Tp53<- prop.test(x = c(275,9), n = c(1505,50))
Tp53 #1


###skin with lympho node##
braf <- prop.test(x = c(298,55), n = c(522,113))
braf ##0.12

brafv600e<- prop.test(x = c(262,33), n = c(502,78))
brafv600e ##0.13

nras <- prop.test(x = c(81,14), n = c(434,74))
nras ##1

cdkn2a<-prop.test(x = c(23,4), n = c(176,38))
cdkn2a ##0.88

nrasq61k<-prop.test(x = c(9,5), n = c(229,38))
nrasq61k #0.048

nrasq61r<-prop.test(x = c(25,4), n = c(222,74))
nrasq61r#0.21

pten<-prop.test(x = c(8,2), n = c(93,18))
pten#1


##skin with CNS##

braf <- prop.test(x = c(134,55), n = c(276,113))
braf ##1

brafv600e<- prop.test(x = c(78,33), n = c(221,78))
brafv600e ##0.33

nras <- prop.test(x = c(71,14), n = c(271,74))
nras ##0.25

cdkn2a<-prop.test(x = c(8,4), n = c(50,38))
cdkn2a ##0.67

nrasq61k<-prop.test(x = c(10,5), n = c(129,38))
nrasq61k #0.48

nrasq61r<-prop.test(x = c(11,4), n = c(143,74))
nrasq61r#0.73

pten<-prop.test(x = c(6,2), n = c(56,18))
pten#1

##lympho wiwth CNS

braf <- prop.test(x = c(134,298), n = c(276,522))
braf ##<0.026

brafv600e<- prop.test(x = c(78,262), n = c(221,502))
brafv600e ##<0.0001

nras <- prop.test(x = c(71,81), n = c(721,434))
nras ###<0.0001

cdkn2a<-prop.test(x = c(8,23), n = c(50,176))
cdkn2a ##0.76

nrasq61k<-prop.test(x = c(10,9), n = c(129,229))
nrasq61k #0.19

nrasq61r<-prop.test(x = c(11,25), n = c(143,222))
nrasq61r#0.35

pten<-prop.test(x = c(6,8), n = c(56,93))
pten #0.89

Tp53<- prop.test(x = c(9,13), n = c(50,124))
Tp53 #0.27



########

################
#####Graphic compare mutation prevalence from primary and metastase#######
grap<- read.csv("comparegenmutation.csv",sep = ';')
dim(grap)
names(grap)
summary(grap)
grap$Group1
library(forcats)
library(ggplot2)
library(ggsignif)
library(ggpattern)
grap$Gene1<-fct_relevel(grap$Gene,
                        c("TERT","BRAF", "BRAFV600E", "TP53", "NRAS", "NF1", "ROS1","CDKN2A","NRASQ61R",
                          "KIT", "PIK3CA",  "PTEN","NRASQ61K"))
grap$Group1<-fct_relevel(grap$Group,c("Primary tumors","Metastases"))

#remotes::install_github("coolbutuseless/ggpattern")

###dat <- data.frame(Group = c("S1", "S1", "S2", "S2"),
#Sub   = c("A", "B", "A", "B"),
#Value = c(3,5,7,8))  

##ggplot(dat, aes(Group, Value)) +
## geom_bar(aes(fill = Sub), stat="identity", position="dodge", width=.5) +
## geom_signif(stat="identity",
##           data=data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),
##                           y=c(5.8, 8.5), annotation=c("**", "NS")),
##           aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
##geom_signif(comparisons=list(c("S1", "S2")), annotations="***",
##             y_position = 9.3, tip_length = 0, vjust=0.4) +
## scale_fill_manual(values = c("grey80", "grey20"))


x<-ggplot(grap, aes(Gene1, Mutation.Prevalence)) +
  geom_bar(aes(fill = Group1), stat="identity", position="dodge", width=.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12) )+
  geom_signif(stat="identity",
              data=data.frame(x=c(0.675, 1.675,2.675,3.675,4.675,5.675,6.675,7.675,8.675,9.675,10.675,11.675), 
                              xend=c(1.125, 2.125,3.125,4.125,5.125,6.125,7.125,8.125,9.125,10.125,11.125,12.125),
                              y=c(52, 49), annotation=c("***", "***", "***","***"," "," ","***"," ","***","***"," "," ")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) + theme_classic(base_size = 14)

ggsave("mutation_plot1.png", plot = x, width = 13, height = 8, dpi = 300)

######Test primary tumor and metastase###
########################################
chr5q33.34<- prop.test(x = c(1,7), n = c(18,61))
chr5q33.34 #0.78

chr6q<- prop.test(x = c(2,11), n = c(22,114))
chr6q #1

chr9p21<-prop.test(x = c(19,41), n = c(32,156))
chr9p21 # <0.001

chr9p21.3<- prop.test(x = c(103,50), n = c(247,174))
chr9p21.3 # 0.008

chr10q<- prop.test(x = c(12,7), n = c(34,16))
chr10q  # 0.78





chr16q24.3<- prop.test(x = c(1,3), n = c(15,61))
chr16q24.3 #1

chr1p<-prop.test(x = c(5,2), n = c(16,19))
chr1p #0.27


chr4p16.3<-prop.test(x = c(4,3), n = c(22,17))
chr4p16.3 #1
 

chr4q12<-prop.test(x = c(90,8), n = c(718,165))
chr4q12. #0.007

chr6p<-prop.test(x = c(13,16), n = c(38,18))
chr6p #<0.0001

chr7q34<-prop.test(x = c(16,9), n = c(122,498))
chr7q34 #<0.0001

chr8q<-prop.test(x = c(12,12), n = c(38,19))
chr8q #0.046


chr11q13.3<- prop.test(x = c(91,9), n = c(250,43))
chr11q13.3  # 0.07

chr12q13.32<-prop.test(x = c(3,5), n = c(66,61))
chr12q13.32 #0.0.63

chr12q14.1<-prop.test(x = c(23,10), n = c(184,189))
chr12q14.1 #0.023

chr12q15<-prop.test(x = c(10,5), n = c(66,81))
chr12q15#0.1


chr17q12<-prop.test(x = c(16,1), n = c(103,17))
chr17q12 #0.49












##########################
##########
#########################
##########
library(dplyr)
library(forcats)

#grap<- read.csv("Compare chromosome .csv",sep = ';')
##grap<- read.csv("compare chromosome plus 30.csv",sep = ';')
grap<- read.csv("CNA compare plus 30 25.csv",sep = ';')
#grap$Chro1<-fct_relevel(grap$Chro,c("9p21",      "9p21.3",    "4q12",   "5p15.33" , "7q34", "11q13.3" ,"12q13.32" , " 12q14.1" , "12q15" ))
dim(grap)
summary(grap)
grap$Group

# B1. lọc riêng nhóm Primary tumors để lấy thứ tự
# B1. Lấy thứ tự Chro dựa trên Prevalance của Primary tumors
order_primary <- grap %>%
  filter(Group == "Primary tumors") %>%
  arrange(Prevalance) %>%       # sắp xếp từ thấp đến cao (âm trước, dương sau)
  pull(Chro) %>%
  unique()

# B2. Áp dụng thứ tự này cho toàn bộ dataframe
grap$Chro <- factor(grap$Chro, levels = order_primary)

# B2. áp dụng thứ tự này cho toàn bộ dataframe (cả Primary và Metastases)
grap$Chro <- factor(grap$Chro, levels = order_primary)
grap$Group1<-fct_relevel(grap$Group,c("Primary tumors","Metastases"))



#ggplot(grap, aes(Chro, Prevalance)) +
  #geom_bar(aes(fill = Group1), stat="identity", position="dodge", width=.8) +
  # geom_signif(stat="identity",
  #             data=data.frame(x=c(0.675,1.675,2.675,3.675,4.675,5.675,6.675,7.675,8.675,9.675,10.675,11.675,12.675), 
  #                             xend=c(1.125, 2.125,3.125,4.125,5.125,6.125,7.125,8.125,9.125,10.125,11.125,12.125,13.125),
  #                             y=c(-13,-11,-61,-43,-45,-9,33,20,14,91,65,16,17), 
  #                             annotation=c("***", " ","***","***","***","***","***","***","***","***","***","***","***")),
  #             aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  #theme_classic(base_size = 14)

x<-ggplot(grap, aes(x = Chro, y = Prevalance, fill = Group1)) +
  geom_bar(stat = "identity", position = "dodge", width = .8) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("mutation_plot4.png", plot = x, width = 10, height = 8, dpi = 300)











































##########Acral and non acral#####

# setwd("~/Dropbox/These/France/Paper Duc/statistical manuscript")
mut<- read.csv("mutationprimaryver168acral.csv",sep = ';')

dim(mut)
names(mut)
summary(mut)
table(mut$Method)
mut$Fixe<-as.factor(mut$Fixe)
mut$Fixe<-relevel(mut$Fixe,ref='Frozen')
mut_fixe<-subset(mut,is.na(Fixe)==FALSE)
mut_multi<-subset(mut,is.na(Mutilple)==FALSE)
dim(mut_fixe)
dim(mut_multi)
##mutation BRAF##
#mut1<-subset(mut, is.na(Method)==FALSE)
#mut1_braf<-subset(mut1, is.na(BRAF)==FALSE)
mut_braf<-subset(mut, is.na(BRAF)==FALSE)
dim(mut_braf)
names(mut_braf)
m2<- metaprop(BRAF, Sample, studlab=TT, sm="PFT", data=mut_braf, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##method subgroup
mut1<-subset(mut, is.na(Method)==FALSE)
mut1_braf<-subset(mut1, is.na(BRAF)==FALSE)
m1<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT", data=mut1_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m1, Method)
mu1<-update(m1,subgroup=Method)
forest(mu1)


##ALK#
mut_alk<-subset(mut, is.na(ALK)==FALSE)
dim(mut_alk)
names(mut_alk)
m2<- metaprop(ALK, Sample, studlab=TT, sm="PFT", data=mut_alk, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##ABL1#
mut_abl1<-subset(mut, is.na(ABL1)==FALSE)
dim(mut_abl1)
names(mut_abl1)
m2<- metaprop(ABL1, Sample, studlab=Author, sm="PFT", data=mut_abl1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##AR1D2##
mut_ARID2 <-subset(mut, is.na(ARID2 )==FALSE)
dim(mut_ARID2 )
names(mut_ARID2 )
m2<- metaprop(ARID2 , Sample, studlab=TT, sm="PFT", data=mut_ARID2 , method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##BRAF V600E##
mut_BRAFV600E<-subset(mut, is.na(BRAF.V600E)==FALSE)
dim(mut_BRAFV600E)
names(mut_BRAFV600E)
m2<- metaprop(BRAF.V600E, Sample, studlab=TT, sm="PFT", data=mut_BRAFV600E, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##CDKN2A##
mut_CDKN2A<-subset(mut, is.na(CDKN2A)==FALSE)
dim(mut_CDKN2A)
names(mut_CDKN2A)
m2<- metaprop(CDKN2A, Sample, studlab=TT, sm="PFT", data=mut_CDKN2A, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##CTNNB1#
mut_CTNNB1<-subset(mut, is.na(CTNNB1)==FALSE)
dim(mut_CTNNB1)
names(mut_CTNNB1)
m2<- metaprop(CTNNB1, Sample, studlab=TT, sm="PFT", data=mut_CTNNB1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##GNAQ##
mut_GNAQ<-subset(mut, is.na(GNAQ)==FALSE)
dim(mut_GNAQ)
names(mut_GNAQ)
m2<- metaprop(GNAQ, Sample, studlab=TT, sm="PFT", data=mut_GNAQ, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##HRAS##
mut_HRAS<-subset(mut, is.na(HRAS)==FALSE)
dim(mut_HRAS)
names(mut_HRAS)
m2<- metaprop(HRAS, Sample, studlab=TT, sm="PFT", data=mut_HRAS, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##KIT##
mut_KIT<-subset(mut, is.na(KIT)==FALSE)
dim(mut_KIT)
names(mut_KIT)
m2<- metaprop(KIT, Sample, studlab=TT, sm="PFT", data=mut_KIT, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##KRAS##
mut_KRAS<-subset(mut, is.na(KRAS)==FALSE)
dim(mut_KRAS)
names(mut_KRAS)
m2<- metaprop(KRAS, Sample, studlab=TT, sm="PFT", data=mut_KRAS, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MAP2K1##
mut_MAP2K1<-subset(mut, is.na(MAP2K1)==FALSE)
dim(mut_MAP2K1)
names(mut_MAP2K1)
m2<- metaprop(MAP2K1, Sample, studlab=TT, sm="PFT", data=mut_MAP2K1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MET##
mut_MET<-subset(mut, is.na(MET)==FALSE)
dim(mut_MET)
names(mut_MET)
m2<- metaprop(MET, Sample, studlab=TT, sm="PFT", data=mut_MET, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NF1##
mut_NF1<-subset(mut, is.na(NF1)==FALSE)
dim(mut_NF1)
names(mut_NF1)
m2<- metaprop(NF1, Sample, studlab=TT, sm="PFT", data=mut_NF1, method="Inverse", NF1hod.tau="DL",is.na=TRUE)
forest(m2)
##NRAS##
mut_NRAS<-subset(mut, is.na(NRAS)==FALSE)
dim(mut_NRAS)
names(mut_NRAS)
m2<- metaprop(NRAS, Sample, studlab=TT, sm="PFT", data=mut_NRAS, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NRASQ61K
mut_NRASQ61K<-subset(mut, is.na(NRAS.Q61K)==FALSE)
dim(mut_NRASQ61K)
names(mut_NRASQ61K)
m2<- metaprop(NRAS.Q61K, Sample, studlab=TT, sm="PFT", data=mut_NRASQ61K, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##NRAS Q61R
mut_NRASQ61R<-subset(mut, is.na(NRAS.Q61R)==FALSE)
dim(mut_NRASQ61R)
names(mut_NRASQ61R)
m2<- metaprop(NRAS.Q61R, Sample, studlab=TT, sm="PFT", data=mut_NRASQ61R, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##PIK3CA##
mut_PIK3CA<-subset(mut, is.na(PIK3CA)==FALSE)
dim(mut_PIK3CA)
names(mut_PIK3CA)
m2<- metaprop(PIK3CA, Sample, studlab=TT, sm="PFT", data=mut_PIK3CA, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##PTEN##
mut_PTEN<-subset(mut, is.na(PTEN)==FALSE)
dim(mut_PTEN)
names(mut_PTEN)
m2<- metaprop(PTEN, Sample, studlab=TT, sm="PFT", data=mut_PTEN, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##TERT##
mut_TERT<-subset(mut, is.na(TERT)==FALSE)
dim(mut_TERT)
names(mut_TERT)
m2<- metaprop(TERT, Sample, studlab=TT, sm="PFT", data=mut_TERT, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##TP53##
mut_TP53<-subset(mut, is.na(TP53)==FALSE)
dim(mut_TP53)
names(mut_TP53)
m2<- metaprop(TP53, Sample, studlab=TT, sm="PFT", data=mut_TP53, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##ROS1##
mut_ros1<-subset(mut, is.na(ROS1)==FALSE)
dim(mut_ros1)
names(mut_ros1)
m2<- metaprop(ROS1, Sample, studlab=TT, sm="PFT", data=mut_ros1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


########################################
#######################

######

###TABLE5 subgroup analysis##
##method subgroup
mut<- read.csv("primarymutation168.csv",sep = ';')
dim(mut)
names(mut)
mut1_braf<-subset(mut, is.na(BRAF)==FALSE)
mut1_braf$Sample <- as.numeric(mut1_braf$Sample)
m1<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT", data=mut1_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m1, Arcal)
mu1<-update(m1,subgroup=Arcal)
forest(mu1)




##preservation subgroup
mut_fixe_braf<-subset(mut_fixe, is.na(BRAF)==FALSE)
dim(mut_fixe_braf)
m2<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu2<-update(m2,subgroup=Fixe,fixed = FALSE)
metareg(m2, Fixe)
forest(mu2)
# Do meta-regression with two covariates
metareg(m2, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_braf<-subset(mut_multi, is.na(BRAF)==FALSE)
dim(mut_multi_braf)
m3<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu3<-update(m3, subgroup=Mutilple,fixed = FALSE)
metareg(m3,Mutilple)
forest(mu3)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_braf<-subset(mut4, is.na(BRAF)==FALSE)
dim(mut4_braf)
m4<- metaprop(BRAF, Sample,
              studlab=TT, sm="PFT", subgroup = Quality, data=mut4_braf, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)

###NRAS###
#########
mut4<-subset(mut, is.na(Method)==FALSE)
mut4_NRAS<-subset(mut4, is.na(NRAS)==FALSE)
m4<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT", data=mut4_NRAS, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Arcal)
mu4<-update(m4,subgroup=Arcal)
forest(mu4)
##preservation subgroup
mut_fixe_NRAS<-subset(mut_fixe, is.na(NRAS)==FALSE)
dim(mut_fixe_NRAS)
m5<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_NRAS, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu5<-update(m5,subgroup=Fixe,fixed = FALSE)
metareg(m5, Fixe)
forest(mu5)
# Do meta-regression with two covariates
metareg(m5, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_NRAS<-subset(mut_multi, is.na(NRAS)==FALSE)
dim(mut_multi_NRAS)
m6<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_NRAS, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu6<-update(m6, subgroup=Mutilple,fixed = FALSE)
metareg(m6,Mutilple)
forest(mu6)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_nras<-subset(mut4, is.na(NRAS)==FALSE)
m4<- metaprop(NRAS, Sample,
              studlab=TT, sm="PFT", data=mut4_nras, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)


###TP53###
#########
mut7<-subset(mut, is.na(Method)==FALSE)
mut7_TP53<-subset(mut7, is.na(TP53)==FALSE)
m7<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT", data=mut7_TP53, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m7, Arcal)
mu7<-update(m7,subgroup=Arcal)
forest(mu7)
##preservation subgroup
mut_fixe_TP53<-subset(mut_fixe, is.na(TP53)==FALSE)
dim(mut_fixe_TP53)
m5<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_TP53, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu5<-update(m5,subgroup=Fixe,fixed = FALSE)
metareg(m5, Fixe)
forest(mu5)
# Do meta-regression with two covariates
metareg(m5, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_TP53<-subset(mut_multi, is.na(TP53)==FALSE)
dim(mut_multi_TP53)
m9<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_TP53, method="Inverse",
              method.tau="DL",is.na=TRUE)
mu9<-update(m9, subgroup=Mutilple,fixed = FALSE)
metareg(m9,Mutilple)
forest(mu9)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_tp53<-subset(mut4, is.na(TP53)==FALSE)
m4<- metaprop(TP53, Sample,
              studlab=TT, sm="PFT", data=mut4_tp53, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)


###TERT###
#########
mut11<-subset(mut, is.na(Method)==FALSE)
mut11_TERT<-subset(mut11, is.na(TERT)==FALSE)
m11<- metaprop(TERT, Sample,
               studlab=TT, sm="PFT", data=mut11_TERT, method="Inverse",
               method.tau="DL",is.na=TRUE)
metareg(m11, Arcal)
mu11<-update(m11,subgroup=Arcal)
forest(mu11)
##preservation subgroup
mut_fixe_TERT<-subset(mut_fixe, is.na(TERT)==FALSE)
dim(mut_fixe_TERT)
m10<- metaprop(TERT, Sample,
               studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_TERT, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu10<-update(m10,subgroup=Fixe,fixed = FALSE)
metareg(m10, Fixe)
forest(mu10)
# Do meta-regression with two covariates
metareg(m10, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_TERT<-subset(mut_multi, is.na(TERT)==FALSE)
dim(mut_multi_TERT)
m12<- metaprop(TERT, Sample,
               studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_TERT, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu12<-update(m12, subgroup=Mutilple,fixed = FALSE)
metareg(m12,Mutilple)
forest(mu12)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_TERT<-subset(mut4, is.na(TERT)==FALSE)
m4<- metaprop(TERT, Sample,
              studlab=TT, sm="PFT", data=mut4_TERT, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)


###NF1###
#########
mut13<-subset(mut, is.na(Method)==FALSE)
mut13_NF1<-subset(mut13, is.na(NF1)==FALSE)
m13<- metaprop(NF1, Sample,
               studlab=TT, sm="PFT", data=mut13_NF1, method="Inverse",
               method.tau="DL",is.na=TRUE)
metareg(m13, Arcal)
mu13<-update(m13,subgroup=Arcal)
forest(mu13)
##preservation subgroup
mut_fixe_NF1<-subset(mut_fixe, is.na(NF1)==FALSE)
dim(mut_fixe_NF1)
m14<- metaprop(NF1, Sample,
               studlab=TT, sm="PFT",subgroup = Fixe, data=mut_fixe_NF1, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu14<-update(m14,subgroup=Fixe,fixed = FALSE)
metareg(m14, Fixe)
forest(mu14)
# Do meta-regression with two covariates
metareg(m14, Method + Fixe)
##multiple subgroup
names(mut_multi)
mut_multi_NF1<-subset(mut_multi, is.na(NF1)==FALSE)
dim(mut_multi_NF1)
m15<- metaprop(NF1, Sample,
               studlab=TT, sm="PFT",subgroup = Mutilple, data=mut_multi_NF1, method="Inverse",
               method.tau="DL",is.na=TRUE)
mu15<-update(m15, subgroup=Mutilple,fixed = FALSE)
metareg(m15,Mutilple)
forest(mu15)

##quality assessment
mut4<-subset(mut, is.na(Quality)==FALSE)
mut4_NF1<-subset(mut4, is.na(NF1)==FALSE)
m4<- metaprop(NF1, Sample,
              studlab=TT, sm="PFT", data=mut4_NF1, method="Inverse",
              method.tau="DL",is.na=TRUE)
metareg(m4, Quality)
mu4<-update(m4,subgroup=Quality)
forest(mu4)



#########TABLE 6### 
##MEtastatics


meta<- read.csv("Metastaticarcal.csv",sep = ';')

meta_braf<-subset(meta, is.na(BRAF)==FALSE)
dim(meta_braf)
names(meta_braf)
meta1<- metaprop(BRAF, Sample, studlab=TT, sm="PFT", data=meta_braf, method="Inverse", method.tau="DL",is.na=TRUE)
forest(meta1)
names(meta)
##NRAS##
meta_NRAS<-subset(meta, is.na(NRAS)==FALSE)
dim(meta_NRAS)
names(meta_NRAS)
m2<- metaprop(NRAS, Sample, studlab=TT, sm="PFT", data=meta_NRAS, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##BRAF V600E##
meta_BRAFV600E<-subset(meta, is.na(BRAF.V600E)==FALSE)
dim(meta_BRAFV600E)
names(meta_BRAFV600E)
m2<- metaprop(BRAF.V600E, Sample, studlab=TT, sm="PFT", data=meta_BRAFV600E, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##KIT##
meta_KIT<-subset(meta, is.na(KIT)==FALSE)
dim(meta_KIT)
names(meta_KIT)
m2<- metaprop(KIT, Sample, studlab=TT, sm="PFT", data=meta_KIT, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##CDKN2A##
meta_CDKN2A<-subset(meta, is.na(CDKN2A)==FALSE)
dim(meta_CDKN2A)
names(meta_CDKN2A)
m2<- metaprop(CDKN2A, Sample, studlab=TT, sm="PFT", data=meta_CDKN2A, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##CTNNB1#
meta_CTNNB1<-subset(meta, is.na(CTNNB1)==FALSE)
dim(meta_CTNNB1)
names(meta_CTNNB1)
m2<- metaprop(CTNNB1, Sample, studlab=TT, sm="PFT", data=meta_CTNNB1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##MET##
meta_MET<-subset(meta, is.na(MET)==FALSE)
dim(meta_MET)
names(meta_MET)
m2<- metaprop(MET, Sample, studlab=TT, sm="PFT", data=meta_MET, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)



##ERBB4##
meta_ERBB4<-subset(meta, is.na(ERBB4)==FALSE)
dim(meta_ERBB4)
names(meta_ERBB4)
m2<- metaprop(ERBB4, Sample, studlab=TT, sm="PFT", data=meta_ERBB4, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##NRASQ61K
meta_NRASQ61K<-subset(meta, is.na(Q61K)==FALSE)
dim(meta_NRASQ61K)
names(meta_NRASQ61K)
m2<- metaprop(Q61K, Sample, studlab=TT, sm="PFT", data=meta_NRASQ61K, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##TP53##
meta_TP53<-subset(meta, is.na(TP53)==FALSE)
dim(meta_TP53)
names(meta_TP53)
m2<- metaprop(TP53, Sample, studlab=TT, sm="PFT", data=meta_TP53, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##NRAS Q61R
meta_NRASQ61R<-subset(meta, is.na(Q61R)==FALSE)
dim(meta_NRASQ61R)
names(meta_NRASQ61R)
m2<- metaprop(Q61R, Sample, studlab=TT, sm="PFT", data=meta_NRASQ61R, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##PTEN##
meta_PTEN<-subset(meta, is.na(PTEN)==FALSE)
dim(meta_PTEN)
names(meta_PTEN)
m2<- metaprop(PTEN, Sample, studlab=TT, sm="PFT", data=meta_PTEN, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)



##NF1##
meta_NF1<-subset(meta, is.na(NF1)==FALSE)
dim(meta_NF1)
names(meta_NF1)
m2<- metaprop(NF1, Sample, studlab=TT, sm="PFT", data=meta_NF1, method="Inverse", NF1hod.tau="DL",is.na=TRUE)
forest(m2)

##RB1##
meta_rb1<-subset(meta, is.na(RB1)==FALSE)
dim(meta_rb1)
names(meta_rb1)
m2<- metaprop(RB1, Sample, studlab=TT, sm="PFT", data=meta_rb1, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)





####GNA11
meta_GNA11<-subset(meta, is.na(GNA11)==FALSE)
dim(meta_GNA11)
names(meta_GNA11)
m2<- metaprop(GNA11, Sample, studlab=TT, sm="PFT", data=meta_GNA11, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##KRAS##
meta_KRAS<-subset(meta, is.na(KRAS)==FALSE)
dim(meta_KRAS)
names(meta_KRAS)
m2<- metaprop(KRAS, Sample, studlab=?..TT, sm="PFT", data=meta_KRAS, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)
##MAP2K1##
meta_MAP2K1<-subset(meta, is.na(MAP2K1)==FALSE)
dim(meta_MAP2K1)
names(meta_MAP2K1)
m2<- metaprop(MAP2K1, Sample, studlab=?..TT, sm="PFT", data=meta_MAP2K1, method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)




##TERT##
meta_TERT<-subset(meta, is.na(TERT)==FALSE)
dim(meta_TERT)
names(meta_TERT)
m2<- metaprop(TERT, Sample, studlab=TT, sm="PFT", data=meta_TERT, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##ROS1##
meta_ros1<-subset(meta, is.na(ROS1)==FALSE)
dim(meta_ros1)
names(meta_ros1)
m2<- metaprop(ROS1, Sample, studlab=TT, sm="PFT", data=meta_ros1, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)

##PIK3CA##
meta_PIK3CA<-subset(meta, is.na(PIK3CA)==FALSE)
dim(meta_PIK3CA)
names(meta_PIK3CA)
m2<- metaprop(PIK3CA, Sample, studlab=TT, sm="PFT", data=meta_PIK3CA, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


##SF3B1##
meta_SF3B1<-subset(meta, is.na(SF3B1)==FALSE)
dim(meta_SF3B1)
names(meta_SF3B1)
m2<- metaprop(SF3B1, Sample, studlab=TT, sm="PFT", data=meta_SF3B1, Method="Inverse", method.tau="DL",is.na=TRUE)
forest(m2)


######test arcal##

braf <- prop.test(x = c(22,121), n = c(105,952))
braf ##0.03


nras <- prop.test(x = c(96,14), n = c(690,86))
nras # 0.06

brafv600e<- prop.test(x = c(54,2), n = c(489,13))
brafv600e #0.9

tert<- prop.test(x = c(11,8), n = c(125,42))
tert#0.04


nf1<-prop.test(x = c(15,7), n = c(139,42))
nf1 #0.08

kit<-prop.test(x = c(61,12), n = c(498,56))
kit #0.02

pik3ca<-prop.test(x = c(3,2), n = c(84,172))
pik3ca ##0.07



ros1<-prop.test(x = c(5,5), n = c(64,56))
ros1 #0.09
## graph arcal
grap<- read.csv("comapremetasprimaryarcal.csv",sep = ';')
dim(grap)
names(grap)
summary(grap)
grap$Group1

library(ggplot2)
library(ggsignif)
library(ggpattern)
grap$Gene1<-fct_relevel(grap$Gene,c("NRAS", "BRAF","KIT",  "BRAFV600E","NF1", "TERT",  "ROS1","PIK3CA"))
grap$Group1<-fct_relevel(grap$Group,c("Primary tumors","Metastases"))


x<-ggplot(grap, aes(Gene1, Mutation.Prevalence)) +
  geom_bar(aes(fill = Group1), stat="identity", position="dodge", width=.8) +
  geom_signif(stat="identity",
              data=data.frame(x=c(0.675, 1.675,2.675,3.675,4.675,5.675,6.675,7.675), 
                              xend=c(1.125, 2.125,3.125,4.125,5.125,6.125,7.125,8.125),
                              y=c(21.3), annotation=c("*"," ",""," ","",""," "," ")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) + theme_classic(base_size = 14)

ggsave("mutation_plot2 arcal.png", plot = x, width = 13, height = 8, dpi = 300)
