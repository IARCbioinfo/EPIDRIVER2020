#!/usr/bin/Rscript
###########################
#load data
###########################

#libraries and parameters
library(data.table)
library(reshape2)
library(Homo.sapiens)
library(stringr)
library(ggplot2)

genes426<-read.csv("/mojo/andrea/426genes.txt", sep=" ")$x

genelist_<-list(
    DNA_methylation_writers=c("DNMT1","DNMT3A","DNMT3B","DNMT3L"),
    DNA_methylation_editors=c("AICDA","TET1","TET2","TET3","IDH1","IDH2"),
    DNA_methylation_readers=c("MBD1","MBD2","MBD3","MBD4","MBD5","MECP2","UHRF1","UHRF2"),
    
    histone_acetylation_writer=c("CDYL","CDYL2","CLOCK","CREBBP","ELP3","ELP4","EP300","GTF3C4","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B","KAT7","KAT8","NCOA1","NCOA3","MSL3"),
    histone_methylation_writer=c("ASH1L","ASH2L","CARM1","DOT1L","EHMT1","EHMT2","EZH1","EZH2","MLLT10","MLLT6","NSD1","PRDM1","PRDM2","PRDM4","PRDM5","PRDM6","PRDM7","PRDM8","PRDM9","PRDM10","PRDM11","PRDM12","PRDM13","PRDM14","PRDM15","PRDM16","PRMT1","PRMT2","PRMT3","PRMT5","PRMT6","PRMT7","PRMT8","SETD1A","SETD1B","SETD2","SETD3","SETD4","SETD5","SETD6","SETD7","KMT5A","SETDB1","SETDB2","SETMAR","SMYD1","SMYD2","SMYD3","SMYD4","SMYD5","SUV39H1","SUV39H2","KMT5B","KMT5C"),
    histone_acetylation_editor=c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7"),
    histone_methylation_editor=c("JMJD1C","JMJD6","JMJD8","KDM1A","KDM1B","KDM2A","KDM2B","KDM3A","KDM3B","KDM4A","KDM4B","KDM4C","KDM4D","KDM4E","KDM5A","KDM5B","KDM5C","KDM5D","KDM6A","KDM6B","KDM7A","KDM8","UTY"),
    histone_ace_meth_pho_reader=c("PHF1","PHF3","PHF5A","PHF6","PHF7","PHF8","PHF10","PHF11","PHF12","PHF13","PHF14","PHF19","PHF2","PHF20","PHF20L1","PHF21A","PHF21B","PHF23","EP400","GATAD2A","GATAD2B","HCFC1","RAI1"),
    histone_acetylation_reader=c("ATAD2","ATAD2B","BAZ1B","BAZ1A","BAZ2A","BAZ2B","BPTF","BRD1","BRD2","BRD3","BRD4","BRD7","BRD8","BRD9","BRDT","BRPF1","BRPF3","BRWD1","BRWD3","CECR2","PBRM1","PHIP","SP100","SP110","SP140","SP140L","TAF1","TAF1L","TAF3","TRIM24","TRIM28","TRIM33","TRIM66","ZMYND11","ZMYND8"),
    
    helicase=c("ATRX","HELLS","INO80","SMARCA1","SMARCA2","SMARCA4","SMARCA5","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","CHD1","CHD1L","CHD2","CHD3","CHD4","CHD5","CHD6","CHD7","CHD8","CHD9"),
    histones=c("H2AFZ","H3F3A","HIST1H1B","HIST1H1C","HIST1H3B"),
    nonHistone=c("CBX1","CBX2","CBX3","CBX4","CBX5","CBX6","CBX7","CBX8"),
    ringFinger=c("HLTF","RING1","RNF2","RNF17","RNF20","RNF217","RNF40","PCGF1","PCGF2","PCGF5","PCGF6","PHRF1","SHPRH","MARCH5"),
    tudor=c("TDRD1","TDRD10","TDRD12","TDRD3","TDRD5","TDRD6","TDRD7","TDRD9","TDRKH","STK31","SND1","SMNDC1","AKAP1","LBR","MTF2","TP53BP1"),
    phdFinger=c("AIRE","CXXC1","DIDO1","DPF1","DPF2","DPF3","FBXL19","ING1","ING2","ING3","ING4","ING5","G2E3","JADE1","JADE2","JADE3","INTS12","PYGO1","PYGO2","TCF19"),
    pwwp=c("GLYR1","HDGF","HDGFL1","MSH6","MUM1","PSIP1","PWWP2B","ZCWPW1","ZCWPW2"),
    chromatin=c("KANSL1","ATAT1","ACTL6A","ATF7IP","SRCAP","DMAP1"),
    lysineMeth=c("AEBP2","KMT2A","KMT2B","KMT2C","KMT2D","KMT2E","NSD2","NSD3"),
    peptidylArginineDeiminases=c("PADI1","PADI2","PADI3","PADI4","PADI6"),
    ubiquitin=c("UBE2A","UBE2B","UBE2E1","UBE2I","UBR7","USP22","USP27X","USP51","BAP1")
)  
    
genelist<-list(  
    HAT=c("CLOCK","CREBBP","ELP3","EP300","GTF3C4","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B","KAT7","KAT8","ELP4","NCOA1","NCOA3","MSL3","KANSL1"),
    HDAC=c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7"),
    HA_reader=c("ATAD2","ATAD2B","BAZ1B","BAZ1A","BAZ2A","BAZ2B","BPTF","BRD1","BRD2","BRD3","BRD4","BRD7","BRD8","BRD9","BRDT","BRPF1","BRPF3","BRWD1","BRWD3","CECR2","PBRM1","PHIP","SP100","SP110","SP140","SP140L","TAF1","TAF1L","TAF3","TRIM24","TRIM28","TRIM33","TRIM66","ZMYND11","ZMYND8","DPF3"),
    HMT=c("ASH1L","ASH2L","CARM1","DOT1L","EHMT1","EHMT2","EZH1","EZH2","MLLT10","MLLT6","NSD1","PRDM1","PRDM2","PRDM4","PRDM5","PRDM6","PRDM7","PRDM8","PRDM9","PRDM10","PRDM11","PRDM12","PRDM13","PRDM14","PRDM15","PRDM16","PRMT1","PRMT2","PRMT3","PRMT5","PRMT6","PRMT7","PRMT8","SETD1A","SETD1B","SETD2","SETD3","SETD4","SETD5","SETD6","SETD7","SETD8","SETDB1","SETDB2","SETMAR","SMYD1","SMYD2","SMYD3","SMYD4","SMYD5","SUV39H1","SUV39H2","SUV420H1","SUV420H2","AEBP2","KMT2A","KMT2B","KMT2C","KMT2D","KMT2E","WHSC1","WHSC1L1"),
    HDM=c("JMJD1C","JMJD6","JMJD8","KDM1A","KDM1B","KDM2A","KDM2B","KDM3A","KDM3B","KDM4A","KDM4B","KDM4C","KDM4D","KDM4E","KDM5A","KDM5B","KDM5C","KDM5D","KDM6A","KDM6B","KDM7A","KDM8","UTY","HR","MINA"),
    HM_reader=c("CDYL","CDYL2","PHF1","PHF12","PHF19","PHF20","PHF21A","PHF21B","PHF23","CBX1","CBX3","CBX4","CBX5","CBX6","CBX7","TDRD3","SND1","MTF2","AIRE","ING1","ING2","ING3","ING4","ING5","PYGO1","PYGO2","GLYR1","MSH6","PSIP1","ZCWPW1","ZCWPW2","CCDC101","HDGFRP2","MORF4L1","MPHOSPH8","RAG2","SCML2","SFMBT1"),
    
    DNA_methylation_writers=c("DNMT1","DNMT3A","DNMT3B","DNMT3L"),
    DNA_methylation_editors=c("AICDA","TET1","TET2","TET3","IDH1","IDH2"),
    DNA_methylation_readers=c("MBD1","MBD2","MBD3","MBD4","MBD5","MECP2","UHRF1","UHRF2"),

    helicase=c("ATRX","HELLS","INO80","SMARCA1","SMARCA2","SMARCA4","SMARCA5","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","CHD1","CHD1L","CHD2","CHD3","CHD4","CHD5","CHD6","CHD7","CHD8","CHD9"),
    chromatin_remodeling_complex=c("RSF1","POLE3","CHRAC1","RBBP7","RBBP4","MTA1","MTA2","MTA3","ARID1A","ARID2","PHC3","PHC1","PHC2","BMI1","YY1","JARID2","SUZ12","EED","HNF1A","NCOR1","ARID4A","MEN1","RBBP5","WDR5","CHAF1A","CHAF1B","TDG","DAXX","SIN3A","SIN3B","ARID1B","ACTL6B"),
    other=c("PADI1","PADI2","PADI3","PADI4","PADI6","UBE2A","UBE2B","UBE2E1","UBE2I","UBR7","USP22","USP27X","USP51","BAP1")
)

genelist<-list( 
  HM_e=c("JMJD1C","JMJD6","JMJD8","KDM1A","KDM1B","KDM2A","KDM2B","KDM3A","KDM3B","KDM4A","KDM4B","KDM4C","KDM4D","KDM4E","KDM5A","KDM5B","KDM5C","KDM5D","KDM6A","KDM6B","KDM7A","KDM8","UTY","HR","MINA"),
  HM_w=c("ASH1L","ASH2L","CARM1","DOT1L","EHMT1","EHMT2","EZH1","EZH2","MLLT10","MLLT6","NSD1","PRDM1","PRDM2","PRDM4","PRDM5","PRDM6","PRDM7","PRDM8","PRDM9","PRDM10","PRDM11","PRDM12","PRDM13","PRDM14","PRDM15","PRDM16","PRMT1","PRMT2","PRMT3","PRMT5","PRMT6","PRMT7","PRMT8","SETD1A","SETD1B","SETD2","SETD3","SETD4","SETD5","SETD6","SETD7","SETD8","SETDB1","SETDB2","SETMAR","SMYD1","SMYD2","SMYD3","SMYD4","SMYD5","SUV39H1","SUV39H2","SUV420H1","SUV420H2","AEBP2","KMT2A","KMT2B","KMT2C","KMT2D","KMT2E","WHSC1","WHSC1L1"),
  HM_r=c("CDYL","CDYL2","PHF1","PHF12","PHF19","PHF20","PHF21A","PHF21B","PHF23","CBX1","CBX3","CBX4","CBX5","CBX6","CBX7","TDRD3","SND1","MTF2","AIRE","ING1","ING2","ING3","ING4","ING5","PYGO1","PYGO2","GLYR1","MSH6","PSIP1","ZCWPW1","ZCWPW2","CCDC101","HDGFRP2","MORF4L1","MPHOSPH8","RAG2","SCML2","SFMBT1"),
  HA_e=c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7"),
  HA_w=c("CLOCK","CREBBP","ELP3","EP300","GTF3C4","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B","KAT7","KAT8","ELP4","NCOA1","NCOA3","MSL3","KANSL1"),
  HA_r=c("ATAD2","ATAD2B","BAZ1B","BAZ1A","BAZ2A","BAZ2B","BPTF","BRD1","BRD2","BRD3","BRD4","BRD7","BRD8","BRD9","BRDT","BRPF1","BRPF3","BRWD1","BRWD3","CECR2","PBRM1","PHIP","SP100","SP110","SP140","SP140L","TAF1","TAF1L","TAF3","TRIM24","TRIM28","TRIM33","TRIM66","ZMYND11","ZMYND8","DPF3"),
  DM_w=c("DNMT1","DNMT3A","DNMT3B","DNMT3L"),
  DM_e=c("AICDA","TET1","TET2","TET3","IDH1","IDH2"),
  DM_r=c("MBD1","MBD2","MBD3","MBD4","MBD5","MECP2","UHRF1","UHRF2"),
  ChRC=c("RSF1","POLE3","CHRAC1","RBBP7","RBBP4","MTA1","MTA2","MTA3","ARID1A","ARID2","PHC3","PHC1","PHC2","BMI1","YY1","JARID2","SUZ12","EED","HNF1A","NCOR1","ARID4A","MEN1","RBBP5","WDR5","CHAF1A","CHAF1B","TDG","DAXX","SIN3A","SIN3B","ARID1B","ACTL6B"),
  Helicases=c("ATRX","HELLS","INO80","SMARCA1","SMARCA2","SMARCA4","SMARCA5","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","CHD1","CHD1L","CHD2","CHD3","CHD4","CHD5","CHD6","CHD7","CHD8","CHD9"),
  others=c("PADI1","PADI2","PADI3","PADI4","PADI6","UBE2A","UBE2B","UBE2E1","UBE2I","UBR7","USP22","USP27X","USP51","BAP1")
)

genes426<-read.csv("/mojo/andrea/426genes.txt", sep=" ")$x
HM_e=c("JMJD1C","JMJD6","JMJD8","KDM1A","KDM1B","KDM2A","KDM2B","KDM3A","KDM3B","KDM4A","KDM4B","KDM4C","KDM4D","KDM4E","KDM5A","KDM5B","KDM5C","KDM5D","KDM6A","KDM6B","KDM7A","KDM8","UTY")
HM_w=c("ASH1L","ASH2L","CARM1","DOT1L","EHMT1","EHMT2","EZH1","EZH2","MLLT10","MLLT6","NSD1","PRDM1","PRDM2","PRDM4","PRDM5","PRDM6","PRDM7","PRDM8","PRDM9","PRDM10","PRDM11","PRDM12","PRDM13","PRDM14","PRDM15","PRDM16","PRMT1","PRMT2","PRMT3","PRMT5","PRMT6","PRMT7","PRMT8","SETD1A","SETD1B","SETD2","SETD3","SETD4","SETD5","SETD6","SETD7","SETD8","SETDB1","SETDB2","SETMAR","SMYD1","SMYD2","SMYD3","SMYD4","SMYD5","SUV39H1","SUV39H2","SUV420H1","SUV420H2","AEBP2","KMT2A","KMT2B","KMT2C","KMT2D","KMT2E","WHSC1","WHSC1L1")
HM_r=c("CDYL","CDYL2","PHF1","PHF6","PHF19","PHF20","PHF21A","PHF21B","PHF23","CBX1","CBX3","CBX4","CBX5","CBX6","CBX7","TDRD3","SND1","MTF2","AIRE","ING1","ING2","ING3","ING4","ING5","PYGO1","PYGO2","GLYR1","MSH6","PSIP1","ZCWPW1","ZCWPW2","CCDC101","HDGFRP2","MORF4L1","MPHOSPH8","RAG2","SCML2","SFMBT1","PHF20L1", "GATAD2B", "GATAD2A")
HA_e=c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7")
HA_w=c("CLOCK","CREBBP","ELP3","EP300","EP400","GTF3C4","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B","KAT7","KAT8","ELP4","NCOA1","NCOA3","MSL3","KANSL1")
HA_r=c("ATAD2","ATAD2B","BAZ1B","BAZ1A","BAZ2A","BAZ2B","BPTF","BRD1","BRD2","BRD3","BRD4","BRD7","BRD8","BRD9","BRDT","BRPF1","BRPF3","BRWD1","BRWD3","CECR2","PBRM1","PHIP","SP100","SP110","SP140","SP140L","TAF1","TAF1L","TAF3","TRIM24","TRIM28","TRIM33","TRIM66","ZMYND11","ZMYND8","DPF3")
DM_w=c("DNMT1","DNMT3A","DNMT3B","DNMT3L")
DM_e=c("AICDA","TET1","TET2","TET3","IDH1","IDH2")
DM_r=c("MBD1","MBD2","MBD3","MBD4","MBD5","MECP2","UHRF1","UHRF2")
ChRC=c("RSF1","POLE3","CHRAC1","RBBP7","RBBP4","MTA1","MTA2","MTA3","ARID1A","ARID2","PHC3","PHC1","PHC2","BMI1","YY1","JARID2","SUZ12","EED","HNF1A","NCOR1","ARID4A","MEN1","RBBP5","WDR5","CHAF1A","CHAF1B","TDG","DAXX","SIN3A","SIN3B","ARID1B","ACTL6B","DPF1","ACTL6A", "CBX8", "CBX2", "RNF2")
Helicases=c("ATRX","HELLS","INO80","SMARCA1","SMARCA2","SMARCA4","SMARCA5","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2","SMARCD3","SMARCE1","CHD1","CHD1L","CHD2","CHD3","CHD4","CHD5","CHD6","CHD7","CHD8","CHD9")
Others=genes426[! genes426 %in% c(HM_e,HM_w,HM_r,HA_e,HA_w,HA_r,DM_w,DM_e,DM_r,ChRC,Helicases)]
genelist<-list(HM_e=HM_e,HM_w=HM_w,HM_r=HM_r,HA_e=HA_e,HA_w=HA_w,HA_r=HA_r,DM_w=DM_w,DM_e=DM_e,DM_r=DM_r,ChRC=ChRC,Helicases=Helicases,Others=Others)

load_data<-function(path="./"){
  
  print(path)
  data<-data.table()
  cancers<-list.files(paste0(path,"/cancer/"))
  
  for (cancer in cancers){
    
    print (cancer)
    print("read files")
    files<-list.files(paste0(path,"/cancer/", cancer))
    mutfile<-files[grep("mut",files)]
    cnvfile<-files[grep("CNV",files)]
    expfile<-files[grep("mRNA",files)]
    cntfile<-paste0(path,"/expression/5330539/",tolower(cancer),"-rsem-count-tcga-t.txt")
    
    print("process cnvs")
    cnv<-read.csv(paste0(path,"/cancer/", cancer, "/", cnvfile), sep="\t", header=T, comment.char = '#', row.names=2)
    cnv<-cnv[,-c(1,ncol(cnv))]
    cnv[is.na(cnv)]<-0
    cnv<-melt(as.matrix(cnv))
    cnv<-cbind(rep(cancer,nrow(cnv)), cnv)
    colnames(cnv)<-c("cancer","genes","ind","cnv")
    tmp<-data.table(cnv)
    
    print("process mutations")
    mut<-read.csv(paste0(path,"/cancer/", cancer, "/", mutfile), sep="\t", header=T, comment.char = '#', row.names=2)
    mut<-mut[,-c(1,ncol(mut))]
    mut<-melt(as.matrix(mut))
    mut<-cbind(rep(cancer,nrow(mut)), mut)
    colnames(mut)<-c("cancer","genes","ind","mut")
    tmp<-merge(tmp, data.table(mut),all.x=T, all.y=T)
    
    print("process Zscores expression")
    exp<-read.csv(paste0(path,"/cancer/", cancer, "/", expfile), sep="\t", header=T, comment.char = '#', row.names=2)
    exp<-exp[,-c(ncol(exp))]
    exp[is.na(exp)]<-0
    exp<-melt(as.matrix(exp))
    exp<-cbind(rep(cancer,nrow(exp)), exp)
    colnames(exp)<-c("cancer","genes","ind","zscore")
    tmp<-merge(tmp, data.table(exp),all.x=T, all.y=T)
    
    print("process count expression")
    if(file.exists(cntfile)){
      TCGAt<-fread(cntfile)
      TCGAt<-TCGAt[Hugo_Symbol %in% genes426,-2]
      TCGAt<-melt(TCGAt[,unique(colnames(TCGAt)), with=FALSE], id.vars="Hugo_Symbol" )
      TCGAt$variable<-substr(TCGAt$variable,1,15)
      TCGAt$variable<-gsub("-",".",TCGAt$variable)
      colnames(TCGAt)<-c("genes","ind","exp")
      TCGAt<-TCGAt[ , .(exp=mean(exp)), by=c("genes","ind") ]
      tmp<-merge(tmp, TCGAt, by=c("genes","ind"), all.x=T, all.y=T)
      
    }else if ( cancer=="COADREAD" ){
      cntfile<-paste0(path,"/expression/5330539/read-rsem-count-tcga-t.txt")
      READ<-fread(cntfile)
      READ<-READ[Hugo_Symbol %in% genes426,-2]
      READ<-melt(READ[,unique(colnames(READ)), with=FALSE], id.vars="Hugo_Symbol" )
      cntfile<-paste0(path,"/expression/5330539/coad-rsem-count-tcga-t.txt")
      COAD<-fread(cntfile)
      COAD<-COAD[Hugo_Symbol %in% genes426,-2]
      COAD<-melt(COAD[,unique(colnames(COAD)), with=FALSE], id.vars="Hugo_Symbol" )
      TCGAt<-rbind(READ, COAD, fill=T, all.y=T)
      
      TCGAt$variable<-substr(TCGAt$variable,1,15)
      TCGAt$variable<-gsub("-",".",TCGAt$variable)
      colnames(TCGAt)<-c("genes","ind","exp")
      
      TCGAt<-TCGAt[ , .(exp=mean(exp)), by=c("genes","ind") ]
      tmp<-merge(tmp, TCGAt, by=c("genes","ind"), all.x=T, all.y=T)
      
    }else { tmp$exp<-NA }
    
    print("process clinical data")
    clinfile<-fread(paste0(path,"/cancer/", cancer, "/data_bcr_clinical_data_patient.txt"), skip=4)[,c("PATIENT_ID","SEX"), ]
    clinfile$PATIENT_ID<-gsub("-",".",clinfile$PATIENT_ID)
    clinfile$PATIENT_ID<-paste0(clinfile$PATIENT_ID,".01")
    colnames(clinfile)[1]<-"ind"
    tmp<-merge(tmp,clinfile, by="ind", all.x=T, all.y=T)
    
    data<-rbind(data, tmp, fill=T)
    
  }
  #change COADREAD name
  levels(data$cancer)[6]<-"COAD/READ"
  
  #count mutations
  data$mut[data$mut=="NaN"]=NA
  data[,mutStatus:=str_count(mut, ",")+1 ]
  data$mutStatus[is.na(data$mut)]=0
  
  #mutation Status
  data$mutStatus[is.na(data$mut)]=0
  data$mutStatus[!is.na(data$mut)]=1
  
  #exon length
  k = keys(Homo.sapiens, keytype="SYMBOL")
  res <- AnnotationDbi::select(Homo.sapiens, keys = k, columns =c("TXCHROM","EXONSTART","EXONEND"), keytype="SYMBOL", multiVals='first')
  res <- data.table(res)
  res <- res[ , .( SIZE=sum( EXONEND-EXONSTART ) ), by=SYMBOL ]
  
  data<-merge(data, res, by.x="genes", by.y="SYMBOL")
  
  #RIOX1 length not database, add it manually !
  data[ genes=="RIOX1", ]$SIZE = 2462
  
  #genes exome size mean
  gsm<-mean(data[ , .( mean(SIZE) ) , by="genes" ]$V1)
  #6470.972
  
  #calculate mutation coef for size normalisation
  data$mutCoef<-data$mutStatus/data$SIZE*gsm
  
  ##########################################################################################
  #geneclass
  df<-do.call(rbind,lapply(genelist,data.frame))
  df$geneclass<-rownames(df)
  df<-data.table(df)
  df[ , geneclass:=tstrsplit(geneclass,"\\.",keep=1)]
  colnames(df)[1]<-"genes"
  data<-merge(data,df[geneclass!="drivers"], by="genes",all.x=T)
  data[is.na(geneclass),]$geneclass<-"None"
  
  save(data, file="data.rdata")
  return(data)
}


##########################################################################################

load_data_all<-function(path="./"){
  
  print(path)
  data<-data.table()
  cancers<-list.files(paste0(path,"/cancers_all_genes/"),pattern="^.{2,8}$")
  
  for (cancer in cancers){
    
    print (cancer)
    print("read files")
    files<-list.files(paste0(path,"/cancers_all_genes/", cancer))
    mutfile<-files[grep("mutations_mskcc",files)]
    cnvfile<-files[grep("data_CNA",files)]
    zscfile<-files[ grep("data_RNA_Seq.*Zscore",files, perl=T) ]
    expfile<-files[ grep("data_RNA_Seq_v2_expression_median.txt",files, perl=T) ]
    clinfile<-files[ grep("data_bcr_clinical_data_patient.txt",files, perl=T) ]
    #cntfile<-paste0(path,"/expression/5330539/",tolower(cancer),"-rsem-count-tcga-t.txt")
    
    print("process cnvs")
    cnv<-read.csv(paste0(path,"/cancers_all_genes/", cancer, "/", cnvfile), sep="\t", header=T, comment.char = '#', row.names=1)
    cnv<-cnv[,-1]
    cnv[is.na(cnv)]<-0
    cnv<-melt(as.matrix(cnv))
    cnv<-data.table(cnv)
    cnv$cancer<-cancer
    colnames(cnv)<-c("genes","ind","cnv","cancer")
    cnv$genes<-gsub("\\|.*","",cnv$genes)
    cnv<-cnv[, .( cnv=round(mean(cnv)) ) , by=c("ind","genes","cancer")]
    
    print("process mutations")
    mut<-fread(paste0(path,"/cancers_all_genes/", cancer, "/", mutfile), sep="\t", header=T)
    mut<-mut[ , .( mut=.N ) , by=c("Hugo_Symbol","Tumor_Sample_Barcode") ]
    mut$cancer=cancer
    colnames(mut)<-c("genes","ind","mut","cancer")
    mut$ind<-gsub("-",".",mut$ind)
    tmp<-merge(cnv, mut,all.x=T, all.y=T)
    
    print("process Zscores expression")
    exp<-read.csv(paste0(path,"/cancers_all_genes/", cancer, "/", zscfile), sep="\t", header=T, comment.char = '#')
    exp<-exp[ !duplicated(exp[,1]), ]
    rownames(exp)<-exp[,1]
    exp<-exp[,-c(1,2)]
    exp[is.na(exp)]<-0
    exp<-melt(as.matrix(exp))
    exp<-data.table(exp)
    exp$cancer<-cancer
    colnames(exp)<-c("genes","ind","zscore","cancer")
    tmp<-merge(tmp, exp,all.x=T, all.y=T)
    
    print("process rsemcount median expression")
    exp<-read.csv(paste0(path,"/cancers_all_genes/", cancer, "/", expfile), sep="\t", header=T, comment.char = '#')
    exp<-exp[ !duplicated(exp[,1]), ]
    exp<-exp[ !is.na(exp[,1]), ]
    rownames(exp)<-exp[,1]
    exp<-exp[,-c(1,2)]
    exp<-melt(as.matrix(exp))
    exp<-data.table(exp)
    exp$cancer<-cancer
    colnames(exp)<-c("genes","ind","rsemcount","cancer")
    tmp<-merge(tmp, exp,all.x=T, all.y=T)
    
    #print("process clinical data")
    #clinfile<-fread(paste0(path,"/cancer/", cancer, "/", clinfile), skip=4)
    #if (!is.null(clinfile$TUMOR_STATUS)){
    #  clinfile<-clinfile[,c("PATIENT_ID","TUMOR_STATUS","SEX")]
    #}else
    #{
    #  clinfile<-clinfile[,c("PATIENT_ID","SEX"), ]
    #}
    #clinfile$PATIENT_ID<-gsub("-",".",clinfile$PATIENT_ID)
    #clinfile$PATIENT_ID<-paste0(clinfile$PATIENT_ID,".01")
    #colnames(clinfile)[1]<-"ind"
    #tmp<-merge(tmp,clinfile, by="ind", all.x=T, all.y=T)
    
    #print("process count expression")
    #if(file.exists(cntfile)){
    #  TCGAt<-fread(cntfile)
    #  TCGAt<-TCGAt[Hugo_Symbol %in% genes426,-2]
    #  TCGAt<-melt(TCGAt[,unique(colnames(TCGAt)), with=FALSE], id.vars="Hugo_Symbol" )
    #  TCGAt$variable<-substr(TCGAt$variable,1,15)
    #  TCGAt$variable<-gsub("-",".",TCGAt$variable)
    #  colnames(TCGAt)<-c("genes","ind","exp")
    #  TCGAt<-TCGAt[ , .(exp=mean(exp)), by=c("genes","ind") ]
    #  tmp<-merge(tmp, TCGAt, by=c("genes","ind"), all.x=T)
      
    #}else if ( cancer=="COADREAD" ){
    #  cntfile<-paste0(path,"/expression/5330539/read-rsem-count-tcga-t.txt")
    #  READ<-fread(cntfile)
    #  READ<-READ[Hugo_Symbol %in% genes426,-2]
    #  READ<-melt(READ[,unique(colnames(READ)), with=FALSE], id.vars="Hugo_Symbol" )
    #  cntfile<-paste0(path,"/expression/5330539/coad-rsem-count-tcga-t.txt")
    #  COAD<-fread(cntfile)
    #  COAD<-COAD[Hugo_Symbol %in% genes426,-2]
    #  COAD<-melt(COAD[,unique(colnames(COAD)), with=FALSE], id.vars="Hugo_Symbol" )
    #  TCGAt<-rbind(READ, COAD, fill=T)
      
    #  TCGAt$variable<-substr(TCGAt$variable,1,15)
    #  TCGAt$variable<-gsub("-",".",TCGAt$variable)
    #  colnames(TCGAt)<-c("genes","ind","exp")
      
    # TCGAt<-TCGAt[ , .(exp=mean(exp)), by=c("genes","ind") ]
    #  tmp<-merge(tmp, TCGAt, by=c("genes","ind"), all.x=T)
      
    #}else { tmp$exp<-NA }
    
    #print("process clinical data")
    #clinfile<-fread(paste0(path,"/cancer/", cancer, "/data_bcr_clinical_data_patient.txt"), skip=4)[,c("PATIENT_ID","SEX"), ]
    #clinfile$PATIENT_ID<-gsub("-",".",clinfile$PATIENT_ID)
    #clinfile$PATIENT_ID<-paste0(clinfile$PATIENT_ID,".01")
    #colnames(clinfile)[1]<-"ind"
    #tmp<-merge(tmp,clinfile, by="ind", all.x=T)
    
    data<-rbind(data, tmp, fill=T)
    
  }
  #change COADREAD name
  data$cancer<-factor(data$cancer)
  levels(data$cancer)[6]<-"COAD/READ"
  
  #mutation Status
  data$mutStatus[is.na(data$mut)]=0
  data$mutStatus[!is.na(data$mut)]=1
  
  #exon length
  k = keys(Homo.sapiens, keytype="SYMBOL")
  res <- AnnotationDbi::select(Homo.sapiens, keys = k, columns =c("TXCHROM","EXONSTART","EXONEND"), keytype="SYMBOL", multiVals='first')
  res <- data.table(res)
  res <- res[ , .( SIZE=sum( EXONEND-EXONSTART ) ), by=SYMBOL ]
  
  data<-merge(data, res, by.x="genes", by.y="SYMBOL")
  
  #RIOX1 length not database, add it manually !
  data[ genes=="RIOX1", ]$SIZE = 2462
  
  #genes exome size mean
  gsm<-mean(data[ , .( mean(SIZE) ) , by="genes" ]$V1)
  #6470.972
  
  #calculate mutation coef for size normalisation
  data$mutCoef<-data$mutStatus/data$SIZE*gsm
  
  ##########################################################################################
  #geneclass
  df<-do.call(rbind,lapply(genelist,data.frame))
  df$geneclass<-rownames(df)
  df<-data.table(df)
  df[ , geneclass:=tstrsplit(geneclass,"\\.",keep=1)]
  colnames(df)[1]<-"genes"
  data<-merge(data,df[geneclass!="drivers"], by="genes",all.x=T)
  data[is.na(geneclass),]$geneclass<-"None"
  
  #save(data, file="data_all.rdata")
  
  return(data)
}












##########################################################################################

load_expression<-function(path="./"){
  
  dic<-list("BLCA"="Bladder",
            "BRCA"="breast",
            "CESC"="cervix",
            "COAD"="colon",
            "ESCA"="esophagus_muc",
            "KICH"="kidney",
            "KIRC"="kidney",
            "KIRP"="kidney",
            "LIHC"="Liver",
            "LUAD"="Lung",
            "LUSC"="Lung",
            "PRAD"="Prostate",
            "READ"="colon",
            "STAD"="Stomach",
            "THCA"="Thyroid",
            "UCEC"="Uterus",
            "CHOL"="cholangiocarcinoma",
            "HNSC"="headandneck"
  )
  
  expression<-data.table()
  cancers<-names(dic)
  
  for (cancer in cancers){
    
    print(cancer)
    
    expr<-data.table()
    
    if(file.exists(paste0(path,tolower(dic[cancer]),"-rsem-count-gtex.txt"))){
      GTEx<-fread(paste0(path,tolower(dic[cancer]),"-rsem-count-gtex.txt"))
      GTEx<-GTEx[Hugo_Symbol %in% genes426,which(!duplicated(colnames(GTEx))), with=FALSE]
      GTEx<-melt( GTEx[,-2] , id.vars="Hugo_Symbol" )
      GTEx$database<-"GTEx"
      GTEx$group<-"normal"
      expr<-rbind(expr,GTEx)
    }
    
    if(file.exists(paste0(path,tolower(cancer),"-rsem-count-tcga.txt"))){
      TCGA<-fread(paste0(path,tolower(cancer),"-rsem-count-tcga.txt"))
      TCGA<-TCGA[Hugo_Symbol %in% genes426,which(!duplicated(colnames(TCGA))), with=FALSE]
      TCGA<-melt( TCGA[,-2] , id.vars="Hugo_Symbol" )
      TCGA$database<-"TCGA"
      TCGA$group<-"normal"
      expr<-rbind(expr,TCGA)
    }
    
    if(file.exists(paste0(path,tolower(cancer),"-rsem-count-tcga-t.txt"))){
      TCGAt<-fread(paste0(path,tolower(cancer),"-rsem-count-tcga-t.txt"))
      TCGAt<-TCGAt[Hugo_Symbol %in% genes426,which(!duplicated(colnames(TCGAt))), with=FALSE]
      TCGAt<-melt( TCGAt[,-2] , id.vars="Hugo_Symbol" )
      TCGAt$database<-"TCGA"
      TCGAt$group<-"tumor"
      expr<-rbind(expr,TCGAt)
    }
    
    #merge data from both database
    expr$cancer<-cancer
    colnames(expr)<-c("genes","samples","exp","database","group","cancer")
    expr$samples<-gsub("-",".",expr$samples)

    rm(GTEx,TCGA,TCGAt)
    expression<-rbind(expression,expr)
  }
  
  
  
  
  save(expression,file="expression.rdata")
  return(expression)
}  
 

##########################################################################################

load_expression_all<-function(path="./"){
  
  dic<-list("BLCA"="Bladder",
            "BRCA"="breast",
            "CESC"="cervix",
            "COAD"="colon",
            "ESCA"="esophagus_muc",
            "KICH"="kidney",
            "KIRC"="kidney",
            "KIRP"="kidney",
            "LIHC"="Liver",
            "LUAD"="Lung",
            "LUSC"="Lung",
            "PRAD"="Prostate",
            "READ"="colon",
            "STAD"="Stomach",
            "THCA"="Thyroid",
            "UCEC"="Uterus",
            "CHOL"="cholangiocarcinoma",
            "HNSC"="headandneck"
  )
  
  expression<-data.table()
  cancers<-names(dic)
  
  for (cancer in cancers){
    
    print(cancer)
    
    expr<-data.table()
    
    if(file.exists(paste0(path,tolower(dic[cancer]),"-rsem-count-gtex.txt"))){
      GTEx<-fread(paste0(path,tolower(dic[cancer]),"-rsem-count-gtex.txt"))
      GTEx<-GTEx[,which(!duplicated(colnames(GTEx))), with=FALSE]
      GTEx<-melt( GTEx[,-2] , id.vars="Hugo_Symbol" )
      GTEx$database<-"GTEx"
      GTEx$group<-"normal"
      expr<-rbind(expr,GTEx)
    }
    
    if(file.exists(paste0(path,tolower(cancer),"-rsem-count-tcga.txt"))){
      TCGA<-fread(paste0(path,tolower(cancer),"-rsem-count-tcga.txt"))
      TCGA<-TCGA[,which(!duplicated(colnames(TCGA))), with=FALSE]
      TCGA<-melt( TCGA[,-2] , id.vars="Hugo_Symbol" )
      TCGA$database<-"TCGA"
      TCGA$group<-"normal"
      expr<-rbind(expr,TCGA)
    }
    
    if(file.exists(paste0(path,tolower(cancer),"-rsem-count-tcga-t.txt"))){
      TCGAt<-fread(paste0(path,tolower(cancer),"-rsem-count-tcga-t.txt"))
      TCGAt<-TCGAt[,which(!duplicated(colnames(TCGAt))), with=FALSE]
      TCGAt<-melt( TCGAt[,-2] , id.vars="Hugo_Symbol" )
      TCGAt$database<-"TCGA"
      TCGAt$group<-"tumor"
      expr<-rbind(expr,TCGAt)
    }
    
    #merge data from both database
    expr$cancer<-cancer
    colnames(expr)<-c("genes","samples","exp","database","group","cancer")
    expr$samples<-gsub("-",".",expr$samples)
    
    rm(GTEx,TCGA,TCGAt)
    expression<-rbind(expression,expr)
  }
  
  save(expression,file="expression_all.rdata")
  return(expression)
}  




load_expression_rpkm_all<-function(path="./"){
  
  dic<-list("BLCA"="Bladder",
            "BRCA"="breast",
            "CESC"="cervix",
            "COAD"="colon",
            "ESCA"="esophagus_muc",
            "KICH"="kidney",
            "KIRC"="kidney",
            "KIRP"="kidney",
            "LIHC"="Liver",
            "LUAD"="Lung",
            "LUSC"="Lung",
            "PRAD"="Prostate",
            "READ"="colon",
            "STAD"="Stomach",
            "THCA"="Thyroid",
            "UCEC"="Uterus",
            "CHOL"="cholangiocarcinoma",
            "HNSC"="headandneck"
  )
  
  expression<-data.table()
  cancers<-names(dic)
  
  for (cancer in cancers){
    
    print(cancer)
    
    expr<-data.table()
    
    if(file.exists(paste0(path,tolower(dic[cancer]),"-rsem-fpkm-gtex.txt"))){
      GTEx<-fread(paste0(path,tolower(dic[cancer]),"-rsem-fpkm-gtex.txt"))
      GTEx<-GTEx[,which(!duplicated(colnames(GTEx))), with=FALSE]
      GTEx<-melt( GTEx[,-2] , id.vars="Hugo_Symbol" )
      GTEx$database<-"GTEx"
      GTEx$group<-"normal"
      expr<-rbind(expr,GTEx)
    }
    
    if(file.exists(paste0(path,tolower(cancer),"-rsem-fpkm-tcga.txt"))){
      TCGA<-fread(paste0(path,tolower(cancer),"-rsem-fpkm-tcga.txt"))
      TCGA<-TCGA[,which(!duplicated(colnames(TCGA))), with=FALSE]
      TCGA<-melt( TCGA[,-2] , id.vars="Hugo_Symbol" )
      TCGA$database<-"TCGA"
      TCGA$group<-"normal"
      expr<-rbind(expr,TCGA)
    }
    
    if(file.exists(paste0(path,tolower(cancer),"-rsem-fpkm-tcga-t.txt"))){
      TCGAt<-fread(paste0(path,tolower(cancer),"-rsem-fpkm-tcga-t.txt"))
      TCGAt<-TCGAt[,which(!duplicated(colnames(TCGAt))), with=FALSE]
      TCGAt<-melt( TCGAt[,-2] , id.vars="Hugo_Symbol" )
      TCGAt$database<-"TCGA"
      TCGAt$group<-"tumor"
      expr<-rbind(expr,TCGAt)
    }
    
    #merge data from both database
    expr$cancer<-cancer
    colnames(expr)<-c("genes","samples","exp","database","group","cancer")
    expr$samples<-gsub("-",".",expr$samples)
    
    rm(GTEx,TCGA,TCGAt)
    expression<-rbind(expression,expr)
  }
  
  ##########################################################################################
  #geneclass
  df<-do.call(rbind,lapply(genelist,data.frame))
  df$geneclass<-rownames(df)
  df<-data.table(df)
  df[ , geneclass:=tstrsplit(geneclass,"\\.",keep=1)]
  colnames(df)[1]<-"genes"
  expression<-merge(expression,df[geneclass!="drivers"], by="genes",all.x=T)
  expression[is.na(geneclass),]$geneclass<-"None"
  
  #if expression is nul, set it to minimum to avoid infinite values
  expression[exp==0]$exp<-min(expression[exp>0 ]$exp)

  #calcul log2 exp
  expression[ , log_exp:=log2(exp), ]
  
  save(expression,file="expression_fpkm_all.rdata")
  return(expression)
}  




##########################################################################################

normalize_expression<-function(path="./"){
  
  load("expression_all.rdata")
  
  cancers<-unique(expression$cancer)
  expression<-data.table()
  for (cancer_ in cancers){
    print(cancer_)
    data<-acast(expression[cancer==cancer_,], genes~samples+database+group, value.var="exp" )
    data<-normalize.quantiles(data, copy=FALSE)
    data<-melt(data)
    data<-as.data.table(data)
    data[ , c("samples","database","group"):=tstrsplit(Var2,"_") ]
    colnames(data)<-c("genes","model","exp","samples","database","group")
    data$cancer=cancer_
    expression<-rbind(expression, data,fill=T)
  }
  expression<-tpm[ ,c("genes","samples","exp","database","group","cancer")]
  
  #calcul log2 tmp, then mean of log2 tpm
  expression[ , log_exp:=log2(exp), ]
  logex<-expression[ is.finite(log_exp) , .( mean_log_exp=mean(log_exp)), by=c("cancer","genes","group","database") ]
  logex$class<-paste0(logex$database,'_',logex$group)
  
  #density plots
  fig<-ggplot(logex, aes( x=mean_log_exp, color=class ) ) + geom_density()
  jpeg("expression_all_norm.jpeg", width=2800, height=1200, res = 200)
  fig
  dev.off()
  
  save(expression,file="expression_all_normalized.rdata")
  return(expression)
  
  #load("../data.rdata")
  #ERGs
  #epidrivers<-unique(data$genes)
  #expression<-expression[genes %in% epidrivers]
  #save(expression,file="expression_normalized.rdata")
  #return(expression)
}

normalize_edger_expression<-function(path="./"){
  
  library(reshape2)
  library(edgeR)
  expression_n<-data.table()
  
  for (cancer_ in cancers){
    
    pdata<-unique(expression[ cancer==cancer_, c("ind","TUMOR_STATUS") ])
    data<-acast(expression[cancer==cancer_,], genes~ind+TUMOR_STATUS, value.var="rsemcount" )
    
    #design
    design<-model.matrix(~0 + TUMOR_STATUS, data=pdata)
    
    #remove genes that are unexpressed or very lowly expressed in the samples
    cpm_log <- cpm(data, log = TRUE)
    median_log2_cpm <- apply(cpm_log, 1, median)
    #hist(median_log2_cpm)
    expr_cutoff <- -1
    #abline(v = expr_cutoff, col = "red", lwd = 3)
    sum(median_log2_cpm > expr_cutoff)
    data_clean <- data[median_log2_cpm > expr_cutoff, ]
    
    #log2 cpm
    cpm_log <- cpm(data_clean, log = TRUE)
    
    #differential expression analysis
    y <- DGEList(counts = data_clean, group = pdata$group)
    y <- calcNormFactors(y)
    
    y<-melt(y$count)
    y<-as.data.table(y)
    y[ , c("samples","database","group"):=tstrsplit(Var2,"_") ]
    colnames(y)<-c("genes","model","exp","samples","group")
    y$cancer=cancer_
    y<-y[, c("genes","samples","group","exp")]
    
    #save data
    expression_n<-rbind(expression, y,fill=T)
  }
  expression<-expression_n
  save(expression,file="expression_edger.rdata")

}

files<-list.files(".",pattern = "deseq2")
final=data.table()
for (file in files){
  load(file)
  final<-rbind(final, res,fill=T)
}
colnames(final)<-c("baseMean","logFC","lfcSE","stat","pvalue","FDR","genes","cancer")
save(final,file="final_deseq2_all.RData")

normalize_deseq2_expression<-function(path="./"){
  
  library(data.table)
  library(DESeq2)
  library(reshape2)
  
  load("data_all.rdata")
  data<-data[ TUMOR_STATUS %in% c("WITH TUMOR","TUMOR FREE") ]
  
  cancers<-unique(data$cancer)
  
  expression_n<-data.table()
  for (cancer_ in cancers){
  
    expression<-acast(data[cancer==cancer_,], genes~ind+TUMOR_STATUS, value.var="rsemcount" )
    #data<-round(data)
    
    pdata<-data.frame(t(matrix(unlist(strsplit(colnames(expression),"_")), nrow=2 )))
    colnames(pdata)<-c("samples","group")
    expression[is.na(expression)]=1
    expression<-round(expression)
    
    #design
    design<-model.matrix(~0 + group, data=pdata)
    dds<-DESeqDataSetFromMatrix(countData=as.matrix(expression), colData=pdata, design=~group)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq(dds)
    res <- results(dds)
    names<-rownames(res)
    res<-as.data.table(res)
    res$genes<-names
    res$cancer<-cancer_
    expression_n<-rbind(expression_n, res ,fill=T)
  }
  
  
  
  
}



normalize_deseq2_expression<-function(path="./"){
  library(sva)
  #load("expression_all.rdata")
  load("data_all.rdata")
  data<-data[ TUMOR_STATUS %in% c("WITH TUMOR","TUMOR FREE") ]
  
  for (cancer_ in cancers){
    
    data<-acast(expression[cancer==cancer_], genes~samples+database+group, value.var="exp" )
    pdata<-data.table(t(matrix(unlist(strsplit(colnames(data),"_")), nrow=3 )))
    rownames(pdata)<-colnames(data)
    colnames(pdata)<-c("ind","database","group")
    design<-model.matrix(~group,data=pdata)
    databc<-ComBat(dat=data, batch=pdata$database, mod=design, par.prior=TRUE, prior.plots=TRUE)
  
  }
  
  library(DESeq2)
  data<-acast(expression, genes~cancer+samples+database+group, value.var="exp" )
  data<-round(data)
  pdata<-data.frame(t(matrix(unlist(strsplit(colnames(data),"_")), nrow=4 )))
  rownames(pdata)<-colnames(data)
  colnames(pdata)<-c("cancer","samples","database","group")
  dds<-DESeqDataSetFromMatrix(countData=as.matrix(data), colData=pdata, design=~database)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- estimateSizeFactors(dds)
  data<-counts(dds, normalized=T)
  data<-melt(data)
  data<-data.table(data)
  data[ , c("cancer","samples","database","group"):=tstrsplit(Var2,"_") ]
  colnames(data)<-c("genes","model","exp","cancer","samples","database","group")
  expression<-rbind(expression, data,fill=T)
  
  expression$class<-paste0(expression$database,'_',expression$group)
  expression$exp[expression$exp==0]<-min( expression[exp>0]$exp  )
  expression$log_exp<-log2(expression$exp)
  logexp<-expression[ , .( log_mean_exp=mean(log_exp)), by=c("cancer","genes","class") ]
  
  fig3A<-ggplot( logexp, aes(x=cancer, y=log_mean_exp) ) + 
    geom_boxplot(aes(fill=class), outlier.shape=NA) + 
    scale_y_continuous(limits=c(6,15)) +
    scale_fill_manual(values=palette1) +
    theme(axis.text.x=element_text(angle=50,hjust=1, size=12), axis.title.x=element_blank() ) +
    theme(axis.text.y=element_text(size=12), axis.title.y=element_text(size=20)) +
    theme(legend.text=element_text(size=20), legend.title=element_text(size=20)  ) +
    ylab("mean Log2(TPM)")
  fig3A
}


limma_expression<-function(path="./"){
  
  library(limma)
  library(qqman)
  fdr=0.05
  load("expression_fpkm_all.rdata")
  cancers<-unique(expression$cancer)
  final<-data.table()
  for (cancer_ in cancers){
    
    pdata<-unique(expression[ cancer==cancer_, c("samples","database","group") ])
    data<-acast(expression[cancer==cancer_,], genes~samples+database+group, value.var="exp" )
    #data<-data[rowSums(data)>10000,]
    #design
    nblvls<<-length(levels(as.factor(pdata$database)))
    if ( nblvls > 1) { design<-model.matrix(~0 + group + database, data=pdata) }
    else { design<-model.matrix(~0 + group, data=pdata) }
    cmtx <- makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )
    
    #remove genes that are unexpressed or very lowly expressed in the samples
    
    
    #limma
    fit<-lmFit(data,design, pdata, ndups=1, method="ls")
    rownames(cmtx)<-colnames(fit)
    fitContrasts=contrasts.fit(fit,cmtx)
    eb=eBayes(fitContrasts)
    lmPvals = eb$p.value
    chisq <- qchisq(1-eb$p.value,1)
    lambda<-median(chisq)/qchisq(0.5,1)
    print(lambda)
    jpeg(paste0("qqplot_",cancer_,".jpeg"))
    qq(lmPvals,title=paste0("QQ plot: ","lambda=",lambda))
    dev.off()
    table<-topTable(eb, adjust="BH", number=Inf, p=fdr, sort.by="P")
  
    #save data
    table$cancer<-cancer_
    table$genes<-rownames(table)
    final<-rbind(final, table,fill=T)
  }
  save(final,file="final.rdata")
}


#Lambda
#[1] 23.43511
#[1] Inf
# 14.6993
#[1] 146.3878
#[1] 46.30582
#[1] 64.88445
#[1] 144.1747
#[1] 84.39973
#[1] 97.21228
#[1] Inf
#[1] Inf
#[1] Inf
#[1] 48.782
#[1] 91.73531
#[1] Inf
#[1] 66.58872
#[1] 80.66506
#[1] Inf







#run it !
#data<-load_data("/mojo/andrea/")
#data<-load_data_all("/mojo/andrea/")
#expression<-load_expression("/mojo/andrea/expression/5330539/") 
#expression<-load_expression_all("/mojo/andrea/expression/5330539/")
#expression<-normalize_expression("/mojo/andrea/") 
#expression<-load_expression_rpkm_all("/mojo/andrea/expression/2019/")

