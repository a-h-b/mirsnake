log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

library(vegan)

unique.files <- snakemake@input[['uni']]
ambigue.files <- snakemake@input[['ambi']]

first <- TRUE
for(f in unlist(unique.files)){
  sam <- gsub("isomiR_SEA/","",
              gsub("/out_result_mature_unique.txt","",
                   gsub("/","__",f)))
  if(first){
    miRall <- data.frame("sample"=sam,
                         read.delim(f,stringsAsFactors = F),
                         stringsAsFactors = F)
    first <- FALSE
  }else{
    miRall <- rbind(miRall,
                    data.frame("sample"=sam,
                               read.delim(f,stringsAsFactors = F),
                               stringsAsFactors = F))
  }
}

for(f in unlist(ambigue.files)){
  sam <- gsub("/","__",
              gsub("isomiR_SEA/","",
                   gsub("/out_result_mature_ambigue.txt","",
                        f)))
  miRall <- rbind(miRall,
                  data.frame("sample"=sam,
                             read.delim(f,stringsAsFactors = F),
                             stringsAsFactors = F))
}

totaM <- tapply(miRall$X.count_tags,list(miRall$mirna_info,miRall$sample),sum) 
totaM[is.na(totaM)] <- 0
totaM <- totaM[rowSums(totaM)>0,]
rownames(totaM) <- gsub(">","",gsub(" MIMAT.+","",rownames(totaM)))
colnames(totaM) <- gsub("M","M0",colnames(totaM))
totaN <- decostand(totaM,"total",2)
totaMB <- totaM
totaMB[totaMB>0] <-1

save.image(snakemake@output[[1]])

