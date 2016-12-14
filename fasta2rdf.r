###fasta2rdf.r###
### Created November 2016 ###
### Author: Austin Reynolds ###
### Email: awreynolds@utexas.edu ###


#load necessary packages
if(!require(seqinr)){
  install.packages("seqinr")
}
if(!require(optparse)){
  install.packages("optparse")
}

library(seqinr)
library(optparse)

###Define options
option_list <- list(
  make_option(c("-i", "--input"), action="store",
              help="Input path/to/fasta/file"),
  make_option(c("-r", "--reference"), action="store",
              help="Reference path/to/fasta/file. must be same length as input"),
  make_option(c("-a", "--position_start"), action="store", type="numeric",
              help="What base position your FASTA sequences start on"),
  make_option(c("-z", "--position_end"), action="store", type="numeric",
              help="What base position your FASTA sequences end on"),
  make_option(c("-o", "--output"), action="store",
              help = "Output path/to/rdf/file")
)
opt<-parse_args(OptionParser(option_list=option_list))
#print(opt$input)
#print(opt$position_start)
#print(opt$position_end)
#print(opt$output)

#read in fasta file
fasta_input<-read.fasta(opt$input)
fasta_input<-data.frame(INDS=names(fasta_input), SEQS=unlist(getSequence(fasta_input, as.string=T)))

###read in reference
fasta_reference<-read.fasta(opt$reference)
fasta_reference<-data.frame(INDS=names(fasta_reference), SEQS=unlist(getSequence(fasta_reference, as.string=T)))
###split the reference
split_ref<-strsplit(as.vector(t(fasta_reference$SEQS[1])),"")
split_ref_df<-as.data.frame(split_ref)
split_ref_df<-t(split_ref_df)
rownames(split_ref_df)<-"reference"

#prepare fasta file to be split
split_seqs<-strsplit(as.vector(t(fasta_input$SEQS[1])),"")
split_seqs_df<-as.data.frame(split_seqs)
split_seqs_df<-t(split_seqs_df)

#make new dataframe
for (i in 2:nrow(fasta_input)){
  #print(fasta_input$INDS[i])
  df1<-strsplit(as.vector(t(fasta_input$SEQS[i])),"")
  df2<-as.data.frame(df1)
  df3<-t(df2)
  #print(nrow(df3))
  split_seqs_df<-rbind(split_seqs_df,df3)
  #print(nrow(split_seqs_df))
}

#assign col and row names
cnames<-seq(opt$position_start,opt$position_end)
rnames<-as.vector(fasta_input$INDS)
rownames(split_seqs_df)<-rnames
colnames(split_seqs_df)<-cnames

###combine ref with fasta
split_seqs_df<-rbind(split_ref_df,split_seqs_df)

###filter to only polymorphic sites
split_seqs_df[split_seqs_df == "?"] <- NA
unique_count<-c()
for (i in 1:ncol(split_seqs_df)){
  #print(i)
  #print(split_seqs_df[i,i])
  if(anyNA(split_seqs_df[,i])==FALSE) {newval<-length(unique(split_seqs_df[,i]))
  } else {newval<-length(unique(split_seqs_df[,i]))-1
  }
  #print(class(newval))
  unique_count<-c(unique_count, newval)
}
split_seqs_df <- data.frame(subset(split_seqs_df, select=unique_count>1),check.names = FALSE,stringsAsFactors=FALSE)
split_seqs_df <- split_seqs_df[-1,]

#Prepare first three lines to be written to rdf
opener<-"  ;1.0"
positions<-paste(as.numeric(colnames(split_seqs_df)),collapse=';')
positions<-paste(positions,";",sep="")
third<-paste(rep("10;",length(colnames(split_seqs_df))),collapse = '')

#write first three lines of rdf
write(opener,opt$output)
write(positions,opt$output,append = TRUE)
write(third,opt$output,append = TRUE)

#convert NA's back to N for printing
split_seqs_df[is.na(split_seqs_df)] <- "N"
#Prepare and write sample names of polymorphic sites to rdf
for (i in rownames(split_seqs_df)){
  #print(i)
  x<-paste(">",i,";1;;;;;;;",sep="")
  #print(x)
  write(x,opt$output,append = TRUE)
  y<-toupper(paste(as.vector(as.matrix(split_seqs_df[i,])),collapse=""))
  #print(y)
  write(y,opt$output,append = TRUE)
}
