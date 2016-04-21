# Load ExomeDepth library (without warnings)
suppressMessages(library(cn.mops))

# Import parameters from xml wrapper (args_file)
args  <- commandArgs(trailingOnly=TRUE)
param <- read.table(args[1],sep="=", as.is=TRUE)

# Set common parameters
target       <- read.table(param[match("target",param[,1]),2], sep="\t", as.is=TRUE)
padding      <- as.numeric(param[match("padding",param[,1]),2])
mapping_mode <- param[match("mapping_mode",param[,1]),2]
output       <- param[match("output",param[,1]),2]

# Set advanced parameters
advanced_mode <- ifelse("advanced_mode" %in% param[,1], TRUE, FALSE)
if (advanced_mode){
	prior_impact    <- as.integer(param[match("prior_impact",param[,1]),2])
	cyc             <- as.integer(param[match("cyc",param[,1]),2])
	norm_type       <- param[match("norm_type",param[,1]),2]
	norm            <- as.logical(param[match("norm",param[,1]),2])
	upper_threshold <- as.numeric(param[match("upper_threshold",param[,1]),2])
	lower_threshold <- as.numeric(param[match("lower_threshold",param[,1]),2])
	min_width       <- as.integer(param[match("min_width",param[,1]),2])
	seq_alg         <- param[match("seq_alg",param[,1]),2]
	min_read_count  <- as.integer(param[match("min_read_count",param[,1]),2])
	use_median      <- as.logical(param[match("use_median",param[,1]),2])
}

# Create symbolic links for multiple bam and bai 
bam         <- param[param[,1]=="bam",2]
bam_bai     <- param[param[,1]=="bam_bai",2]
bam_label   <- param[param[,1]=="bam_label",2]
bam_label   <- gsub(" ", "_", bam_label)

for(i in 1:length(bam)){
  stopifnot(file.symlink(bam[i], paste(bam_label[i], "bam", sep=".")))
  stopifnot(file.symlink(bam_bai[i], paste(bam_label[i], "bam.bai", sep=".")))
}

# Create genomic ranges
gr  <- GRanges(target[,1],IRanges(target[,2]-padding,target[,3]+padding))
# Merge overlapping segments (make sense if padding != 0)
gr <- reduce(gr)


# Get read counts
BAMFiles <- paste(bam_label, "bam", sep=".")
X <- suppressMessages(getSegmentReadCountsFromBAM(BAMFiles,GR=gr,
												  mode=mapping_mode,
												  sampleNames=bam_label))

if (advanced_mode){
	resCNMOPS <- suppressMessages(exomecn.mops(X,
									priorImpact=prior_impact,
									cyc=cyc,
									normType=norm_type,
									norm=norm,
									upperThreshold=upper_threshold,
									lowerThreshold=lower_threshold,
									minWidth=min_width,
									seqAlgorithm=seq_alg,
									minReadCount=min_read_count,
									useMedian=use_median))
}else{
	resCNMOPS <- suppressMessages(exomecn.mops(X))
}

resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

# Extract individual CNV calls
# Legend for CN values is as follows:
# The copy number classes default to CN0, CN1, CN2, CN3, .., CN8.
# CN2 is the normal copy number for diploid samples.
# CN1 is a heterozygous deletion and 
# CN0 is a homozygous deletion.
# CN3 thru CN8 are amplifications.
# For non-tumor samples the highest copy number class is 8 - higher copy numbers have not been reported.
# CN8 is expected to have 4 times as many reads (for times as high coverage) as CN2.
# For tumor samples very high copy numbers have been observed (e.g. CN64),
# therefore the parameters of cn.mops have to be adjusted to allow for high copy numbers.
# A way to set the parameters is given in https://gist.github.com/gklambauer/8955203
res <- cnvs(resCNMOPS)

# Convert results to data.frame
df <- data.frame(chr=as.character(seqnames(res)),
  				 starts=start(res),
  				 #ranges(res),
  				 ends=end(res),
  				 as.data.frame(elementMetadata(res)))

# Remove CN2 (= normal copy number), if any
df <- df[df$CN != "CN2",]

# Write results
write.table(df, sep='\t', quote=FALSE, file = output,
            row.names = FALSE, dec=".")
