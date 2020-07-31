#!/home/cahaisv/miniconda3/bin/Rscript
#SBATCH --job-name=expression
#SBATCH --ntasks=17
#SBATCH --cpus-per-task=1
#SBATCH --mem=80Gb

library(data.table)
library(DESeq2)
#library(reshape2)


args = commandArgs(trailingOnly=TRUE)
print(args[1])
print(args[2])

cancer_=gsub("expression_","",args[1])
cancer_=gsub(".RData","",cancer_)

load(args[1])

dds<-DESeqDataSetFromMatrix(countData=as.matrix(expression), colData=pdata, design=~group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
names<-rownames(res)
res<-as.data.table(res)
res$genes<-names
res$cancer<-cancer_

cancer_=gsub("/","",cancer_)

save(res,file=paste0("final_deseq2_",cancer_,".RData"))
