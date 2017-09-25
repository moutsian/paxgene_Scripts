#!/bin/bash

# usage for txImportData:
# /software/R-3.3.0/bin/R CMD BATCH --vanilla '--args '
# pathToDir="path/to/dir/with/files"
# pathTosampleFile="listOfSamples"
# typeOfTranscript="kallisto/salmon"
# t2g="transcript2geneFile"
#typeMatrix=c(counts, tpm, rpkm)

#LIBRARY=$1 #e.g. ISR
#RUN=$2 #e.g. 21121
R-3.3.0 CMD BATCH --vanilla '--args module="tximport" pathToDir="/lustre/scratch115/projects/paxgene/SalmonOutput/ISR_21364/" pathTosampleFile="/lustre/scratch115/projects/paxgene/SalmonOutput/allfiles.ISR.21364.txt"  typeOfTranscript="salmon" t2g="transcripts2gene.index" output="/lustre/scratch115/projects/paxgene/SalmonOutput/ISR.21364.transcript_to_gene.counts" typeMatrix="counts"' /lustre/scratch115/projects/paxgene/rnaseq.generalTools.WithArguments.R /lustre/scratch115/projects/paxgene/rnaseq.generalTools.WithArguments.Rout

