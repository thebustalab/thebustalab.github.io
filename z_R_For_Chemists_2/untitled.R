******* YOU NEED TO ALSO BE BLASTING WITH THE QUERYS COMPLIMENT, REVERSE, AND REVERSE COMPLIMENT,
THEN ROTATE ALL OF THOSE ACCORDING TO HOW YOU WOULD EXPECT THEM TO BE ROTATED SO THAT THEY MATCH THE
ORIGINAL QUERY. TRYING TO DO THAT ROTATION COMPUTAITONALLY IS TOO INTENSE.
THEN MERGE ALL THOSE READS INTO ONE FILE 
DO NOT USE A KMER ASEMBLY METHOD - THAT WILL NOT REPOSND TO SNPS
INSTEAD, ALIGN WITH MSA. EVERYTIHNG WILL ALLONG - IT IS A GENE FAMILY!
THEN, PER ALIGNMENT DEISTANCE MATRIX, MAKE CLUSTERS AND GET CONSENSUS TO RECONSTRUCT FAMILY MEMBERS



fastqToFasta <- function(file_in_path, file_out_path) {
	data <- readLines(file_in_path)
	fasta <- DNAStringSet(data[seq(1,length(data),4)+1])
	names(fasta) <- data[seq(1,length(data),4)]
	writeXStringSet(fasta, file_out_path)
}

fastqToFasta(
	file_in_path = "/Users/bust0037/Desktop/test/supplements/OSC_fqs/OSC_reads.fq",
	file_out_path = "/Users/bust0037/Desktop/test/supplements/OSC_fqs/OSC_reads.fa"
)

readFastq <- function(file_in_path) {

	data <- readLines(file_in_path)

	out <- list()
	for (i in 1:(length(data)/4)) {
		base <- paste(data[-2+(i*4)])
		score <- data[0+(i*4)]
		this_read <- data.frame(
			base = as.vector(str_split_fixed(base, pattern = "", n = nchar(base))),
			score = as.vector(str_split_fixed(score, pattern = "", n = nchar(score)))
		)
		this_read$header <- data[-3+(i*4)]
		this_read$score <- phred33_lookup$X2[match(this_read$score, phred33_lookup$X1)]
		out[[i]] <- this_read
	}
	out <- do.call(rbind, out)
	return(out)
}
	
findKmers <- function(fastq) {

	kmer_size <- 19
	kmers <- list()

	pb <- progress::progress_bar$new(total = length(unique(fastq$header)))

	for (i in 1:length(unique(fastq$header))) {
	
		data <- filter(fastq, header == unique(fastq$header)[i])
		
		for (j in 1:(length(data$base)-kmer_size)) {
			kmers <- c(kmers, paste(data$base[j:(j+kmer_size)], collapse = ""))
		}

		pb$tick()

	}
	kmers <- unique(unlist(kmers))
	return(kmers)
}

kmerCompare <- function(fasta, kmers) {

	pb <- progress::progress_bar$new(total = length(kmers))
	out <- list()	
	for (i in 1:length(kmers)) {
		this_kmer <- data.frame(
			presence = rep(0,length(names(fasta)))
		)
		this_kmer$presence[grep(kmers[i], fasta)] <- 1
		out[[i]] <- this_kmer
		pb$tick()
	}
	out <- as.data.frame(do.call(cbind, out))
	if(length(which(colsums(as.matrix(out)) > 1)) == 0 ) {stop("No shared k-mers.")}
	# str(out, strict.width = "cut")
	colnames(out) <- kmers
	out2 <- cbind(
		# data.frame(accession = names(fasta)),
		data.frame(accession = as.character(seq(1,length(names(fasta))))),
		out[,which(colsums(as.matrix(out)) > 1)]
	)
	return(out2)
}

out2 <- kmerCompare(
	fasta = readDNAStringSet("/Users/bust0037/Desktop/test/supplements/OSC_fqs/OSC_reads.fa"),
	kmers = findKmers(readFastq(file_in_path = "/Users/bust0037/Desktop/test/supplements/OSC_fqs/OSC_reads.fq"))
)

clusters <- runMatrixAnalysis(
    data = out2,
    analysis = c("hclust"),
    column_w_names_of_multiple_analytes = NULL,
    column_w_values_for_multiple_analytes = NULL,
    columns_w_values_for_single_analyte = colnames(out2)[2:dim(out2)[2]],
    columns_w_additional_analyte_info = NULL,
    columns_w_sample_ID_info = "accession",
    kmeans = "auto",
    na_replacement = "drop",
    output_format = "wide"
)

ggtree(clusters) + geom_tiplab(offset = 1) +
	theme_classic() +
	geom_tippoint(aes(fill = kmeans_cluster), shape = 21, size = 5)