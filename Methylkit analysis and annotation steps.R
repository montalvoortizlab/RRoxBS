# General methylkit workflow from beta values files to differential analysis and annotation of results tables
# See https://bioconductor.riken.jp/packages/3.5/bioc/vignettes/methylKit/inst/doc/methylKit.html for most detail as this script is largely derived from that
# workflow

## Load packages

library(methylKit)
library(stringr)
library(genomation)


## oxBS raw list -- functions as your READY-TO-GO 5MC estimations

folder <-"/sc/hydra/directory_with_your_files/"
folder_names <- as.list(list.files(folder,pattern="-BS")) 


# Extracting subject ID as a string from the long file name
SAMPLES <- vector()

for (file in folder_names){
     sample.ids= str_replace(file,"") # remove unwanted strings here
     SAMPLES <- append(SAMPLES,sample.ids)
}


# Read the files to a methylRawList object: Make one list each for BS and OxBS files

BS.RawList.PTSD=methRead(folder_names,
           sample.id=SAMPLES,
           assembly="hg38",
		       treatment=treatment,
           context="CpG",
		       dbtype = "tabix",
           dbdir = "methylDB-BS")


## Load in metadata

meta <- read.csv(file="optimal.meta.csv",header=T)
treatment <- as.numeric(meta$Opioid) # set treatment variable as numeric
treatment[treatment==1] <- 0 
treatment[treatment==2] <- 1

## Adjust the BS-seq file with the OXBS to estimate hydroxymethylation

HydroxyRawList <- adjustMethylC(BS.RawList.PTSD,MethRawList.PTSD)

## Filter CpG sites with potential high PCR bias
filtered.HydroxyRawList=filterByCoverage(HydroxyRawList,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

## Normalizing coverage across subjects for each CpG
filt.and.normed <- normalizeCoverage(filtered.HydroxyRawList,method="median")
									  
## Get coverage stats for individual samples		  
Coverage <- getCoverageStats(filtered.HydroxyRawList,plot=FALSE,both.strands=FALSE)

## Unite individual samples to get data in one object ready for differential expression
united.meth=unite(filt.and.normed, destrand=FALSE)



## Differential methylation or hydroxymethylation analysis

my.diffMeth<-calculateDiffMeth(united.meth,
                               overdispersion="MN", # good to adjust for overdispersion (toward 0% methylation) in the data
                               test="Chisq",
                               covariates=meta[,c(4:6,8,9)], # covariate argument using metadata columns
                               mc.cores=4) # can set number of cores for parallel processing, also note 



## Annotation steps for genomic regions and distance to TSS

## Steps to convert an hg38 GTF to hg38 BED file suitable for this 
# run in bash >>>>>>>>>>  awk '{if($3 != "gene") print $0}' Homo_sapiens.GRCh38.84.gtf | grep -v "^#" | gtfToGenePred /dev/stdin /dev/stdout | genePredToBed stdin Homo_sapiens.GRCh38.84.bed

gene.annotation <- readTranscriptFeatures("xxx", # your BED file
                                          remove.unusual = TRUE, 
                                          up.flank = 1000,
                                          down.flank = 1000,
                                          unique.prom = TRUE)


## used to add chr to chromosome name in bed file if missing (if chr is in name of your raw files, it will be needed in the bed file)
sed 's/^/chr/' bed.in > bed.out



## Cpg islands and shores annotations

## GETTING THE ANNOTATIONS
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
   | gunzip -c \
   | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, substr($0, index($0, $7)); }' \
   | sort-bed - \
   > cpgIslandExt.hg38.bed
## remove line 11 from this bed with sed '11d' file.in > file.out

 cpg.obj=readFeatureFlank("cpgIsland.bed",
                           feature.flank.name=c("CpGi","shores"))
diffCpGann=annotateWithFeatureFlank(as(my.diffHydroxy,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")
##################################################################################3  

ANNO <- annotateWithGeneParts(as(my.diffMeth,"GRanges"),gene.annotation)

##  Extracting differential methylation table and annotations like nearest Ensembl transcript TSS
diff.table<- getData(my.diffMeth)
gene.part <- getMembers(ANNO)
cpg.annos <- getMembers(diffCpGann)
tss.dis <-  getAssociationWithTSS(ANNO)  

## Combine Diff meth table and annotations
diff.meth.results <- cbind(diff.table,gene.part,cpg.annos,tss.dis)  


## Get gene names and description info from Ensembl transcript ID 
library("biomaRt") 
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
values <- unique(diff.hydro.results[,15])
data <- getBM(attributes=c("hgnc_symbol","description","gene_biotype","ensembl_transcript_id"),
              filters = "ensembl_transcript_id", 
              values = values, 
              mart= ensembl)
final <- merge(diff.meth.results,data,by.x="feature.name",by.y="ensembl_transcript_id",all.x=T)
save(final,file="Final_data_table.RData") # will be a very large file


















