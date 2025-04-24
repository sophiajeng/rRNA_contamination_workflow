## rRNA Contamination Protocol

### Trimming

Trimming is taken care of bwa_rrna_workflow.wdl
### Determine identifiers and barcodes

This section of R code extract run identifiers from the `Stats.json` file currently provided by the sequencing core for each run.  From
this information a uniquely identifying platform unit (PU) is derived, which is a combination of flowcell ID, lane and barcode.  
For more info see this [article](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).  Each sample is associated with a UUID barcode that will be carried downstream.  For reference, these barcode
assignments are written to a `*barcodes.txt` file with the entire set of identifiers associated with the files written to a `*file_map.txt` file.
Finally, the files used for the basic STAR-Fusion processing step are written: `star_fastqs.txt` and `star_samples.txt`.

`star_fastqs.txt` is a tab-delimited file with the R1 fastq file first followed by the R2 fastq file.
`star_samples.txt` contains the corresponding sample names for the fastq files, which in this case are the barcodes.

```r

library(data.table)
library(stringr)
library(jsonlite)
library(uuid)

stats.files <- list.files(".",pattern="Stats.json", full.names=T, recursive = T)

samp.gen.info <- rbindlist(lapply(stats.files, function(cur.file){
  
  tmp.samp <- read_json(cur.file)
  
  rbindlist(lapply(tmp.samp$ConversionResults, function(x){
    
    rbindlist(lapply(x$DemuxResults, function(y){
      
      stopifnot(length(y$IndexMetrics)==1)
      
      data.table(runid=tmp.samp$RunId, flowcell=tmp.samp$Flowcell, lane=as.character(x$LaneNumber), index=y$IndexMetrics[[1]]$IndexSequence, sample=y$SampleName,  generated=y$NumberReads)
    }))
  }))
  
}))

samp.gen.info[,sample:=sub("RNA\\d+[A-Z]{2}_", "", sample)]

#this one happens to be duplicated as the stats files are the same..

samp.gen.info <- unique(samp.gen.info)

all.fastqs <- list.files(pattern=".fastq.gz", full.names = T, recursive=T)
fastq.dt <- data.table(files=all.fastqs)

fastq.dt[,c("lib","sample","sid","read","batch"):=transpose(lapply(strsplit(basename(files), "_"), function(x){
 sub(".fastq.gz","",x)
}))]

fastq.dt <- fastq.dt[grepl("Undetermined", lib)==F]

#by lane will work here..

fastq.dt.m <- merge(samp.gen.info[,.(lane=paste(sort(lane), collapse="+"), generated=sum(generated)),by=.(runid, flowcell, index, sample)], fastq.dt, by=c("sample"))

fastq.dt.m[,pu:=paste(flowcell,lane,index,sep=".")]

samp.barcodes <- fastq.dt.m[,.(sample=unique(sample))][,.(barcode=UUIDgenerate()), by=sample]

fwrite(samp.barcodes, file="RNA200428ET_barcodes.txt", sep="\t", col.names=T)

fastq.dt.m <- merge(fastq.dt.m, samp.barcodes, by="sample")

fwrite(fastq.dt.m, file="RNA200428ET_file_map.txt", sep="\t", col.names=T)

tmp.cast <- dcast(barcode+lane+pu~read, value.var="files" ,data=fastq.dt.m)
  
fwrite(tmp.cast[,.(R1, R2)], file="star_fastqs.txt", sep="\t", col.names=F, quote=F)

fwrite(tmp.cast[,.(barcode)], file="star_samples.txt", sep="\t", col.names=F, quote=F)

```

### bwa rRNA pipeline

This pipeline uses Trimmomatic to trim and uses bwa to align to rRNA sequences. It then runs samtools flagstat to determine the percentage of reads aligning to rRNA sequences.

```bash

cp star_samples.txt bwa_samples.txt
cp star_fastqs.txt bwa_fastqs.txt

sbatch -o RNA210219_bwa_rrna.out -e RNA210219_bwa_rrna.err -t 2160 -p exacloud -c 2 --mem=5G  \
--wrap "
source /home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/activate gatk4
java -jar -Djava.io.tmpdir=tmp_dir \
-Dconfig.file=/home/exacloud/gscratch/mcweeney_lab/resources/programs/exacloud.cromwell.conf /home/exacloud/gscratch/mcweeney_lab/resources/programs/cromwell-48.jar run \
wdl_sequencing/bwa_rrna_workflow.wdl \
-i wdl_sequencing/bwa_rrna_workflow_inputs.json \
-m RNA210219KW_bwa_rrna.json
"

```
### Analysing flagstat results


```r
srun --time=1440 -p light --pty bash -i
export PATH=/home/exacloud/gscratch/mcweeney_lab/resources/programs/R-3.6.1/bin/:$PATH
export LD_LIBRARY_PATH=/home/exacloud/gscratch/mcweeney_lab/resources/programs/SlurmRlibs/system_libs/:$LD_LIBRARY_PATH
export R_LIBS_USER=/home/exacloud/gscratch/mcweeney_lab/resources/programs/SlurmRlibs/x86_64-redhat-linux-gnu-library/3.6/

library(jsonlite)
library(data.table)
library(plyr)

sf.out <- read_json("RNA210219KW_bwa_rrna.json")

#make sure it actually worked
stopifnot(sf.out$status == "Succeeded")

#copy the fusion calls, note that the entire parental directory is copied
did.work.1 <- sapply(sf.out$outputs$star_fusion.out_fus, function(x) file.copy(from=dirname(x), to="/back/up/dir", recursive=T))

stopifnot(all(did.work.1))
barcode_sample<-fread("RNA210219KW_barcodes.txt")
rrna<-fread("/home/exacloud/gscratch/mcweeney_lab/resources/bwa/ncbi_rrna_9606/rrna_9606.fasta",header=F)
#keep only those that have ">"
rrna<-rrna[grep(">",V1)]
rrna$V1<-gsub(">","",rrna$V1)
rrna_annotation<-lapply(rrna$V1,function(x) {
        tmp<-unlist(strsplit(x,"\\|"))
        return(c(paste(tmp[1:4],collapse="|"),tmp[5]))
})
rrna_annotation<-do.call(rbind,rrna_annotation)
colnames(rrna_annotation)<-c("rRNA","rRNA Description")



rrna_flagstat_list<-lapply(sf.out$outputs[[1]],function(x,barcodes) {
        sample_name<-gsub(".*Sample_","",gsub("_bwa_rrna.out","",x))
        sample_name<-mapvalues(sample_name,from=barcodes$barcode,to=barcodes$sample)
        input_reads<-readLines(x)[1]
        input_reads<-gsub("\"","",gsub(" +.*","",input_reads))
        paired_reads<-readLines(x)[9]
        paired_reads<-gsub("\"","",gsub(" \\+.*","",paired_reads))
        percentage<-(as.numeric(paired_reads)/as.numeric(input_reads))*100
        flag_stat_results<-cbind(sample_name,paired_reads,input_reads,percentage)
        return(flag_stat_results)
},barcode_sample)

rrna_flagstat_dt<-do.call(rbind,rrna_flagstat_list)
unique_reads_list<-lapply(sf.out$outputs[["rrna_workflow.unique_reads_out"]],function(x,barcodes) {
        unique_reads<-fread(x,sep=" ",header=F)
        unique_reads$V1<-gsub(".*Sample_","",gsub("_bwa_rrna_sorted.bam","",unique_reads$V1))
        unique_reads$V1<-mapvalues(unique_reads$V1,from=barcodes$barcode,to=barcodes$sample)
        unique_reads$V2<-unlist(lapply(unique_reads$V2,function(y) {
                tmp<-unlist(strsplit(y,"\\|"))
                return(paste(tmp[1:4],collapse="|"))
        }))

        return(unique_reads)
},barcode_sample)

unique_reads_dt<-do.call(rbind,unique_reads_list)
unique_reads_dt_dcast<-dcast(unique_reads_dt,formula = V2 ~ V1, value.var="V3")
unique_rrna_reads_flagstat<-merge(unique_reads_dt,rrna_flagstat_dt,by.x="V1",by.y="sample_name")
unique_rrna_reads_flagstat_annotated<-merge(rrna_annotation,unique_rrna_reads_flagstat,by.x="rRNA",by.y="V2")
unique_rrna_reads_flagstat_annotated<-data.table(unique_rrna_reads_flagstat_annotated)
colnames(unique_rrna_reads_flagstat_annotated)<-c("rRNA","rRNA_Description","sample","Uniquely_Mapped_specific_rRNA","paired_reads_to_all_rrna","input_reads","all_rrna_percentage")
unique_rrna_reads_flagstat_annotated$input_reads<-as.numeric(unique_rrna_reads_flagstat_annotated$input_reads)
unique_rrna_reads_flagstat_annotated[, Primary_Alignment_rRNA_Pct := (Uniquely_Mapped_specific_rRNA/input_reads)*100]

```


