# in shell:

grep "^#" indels_filtered_mac10.recode.vcf > informative_big_deletions.vcf

grep "^#" indels_filtered_mac10.recode.vcf > informative_big_insertions.vcf

grep "^#" indels_filtered_mac10.recode.vcf > informative_small_deletions.vcf

grep "^#" indels_filtered_mac10.recode.vcf > informative_small_insertions.vcf

grep "^#" snps_filtered_mac10.recode.vcf > informative_snps.vcf



# in R:

options(scipen=999)

# read in a small vcf (don't use for large vcf files)
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}
#########################
######### check indels
#########################
x <- read_vcf("indels_filtered_mac10.recode.vcf")

genotypes <- x[,10:ncol(x)]

allo_virens <- genotypes[ ,sapply(strsplit(colnames(genotypes), "__"), "[[", 1) == "allo" & sapply(strsplit(colnames(genotypes), "__"), "[[", 2) == "virens"]

allo_sordid <- genotypes[ ,sapply(strsplit(colnames(genotypes), "__"), "[[", 1) == "allo" & sapply(strsplit(colnames(genotypes), "__"), "[[", 2) == "sordidulus"]

# loop for each variant
keep <- list()
for(a in 1:nrow(allo_sordid)) {
	a_rep <- substr(allo_sordid[a,], 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- c(substr(a_rep, 1, 1), substr(a_rep, 3, 3))
	if(length(unique(a_rep)) == 1) {
		a_rep2 <- substr(allo_virens[a,], 1, 3)
		a_rep2 <- a_rep2[a_rep2 != "./."]
		a_rep2 <- c(substr(a_rep2, 1, 1), substr(a_rep2, 3, 3))
		if(length(unique(a_rep2)) == 1 & a_rep[1] != a_rep2[1]) {
			keep[[a]] <- a
		}
	}
}
keep <- unlist(keep)

x2 <- x[keep,]

del <- x2[nchar(x2$REF) > 1 & nchar(x2$ALT) == 1,]
ins <- x2[nchar(x2$REF) == 1 & nchar(x2$ALT) > 1,]

big_del <- del[nchar(del$REF) >= 50, ]
small_del <- del[nchar(del$REF) < 50, ]
big_ins <- ins[nchar(ins$ALT) >= 50, ]
small_ins <- ins[nchar(ins$ALT) < 50, ]

dim(big_del)
dim(small_del)
dim(big_ins)
dim(small_ins)

write.table(big_del, file="informative_big_deletions.vcf", quote=F, row.names=F, col.names=F, sep="\t", append=T)

write.table(big_ins, file="informative_big_insertions.vcf", quote=F, row.names=F, col.names=F, sep="\t", append=T)

write.table(small_del, file="informative_small_deletions.vcf", quote=F, row.names=F, col.names=F, sep="\t", append=T)

write.table(small_ins, file="informative_small_insertions.vcf", quote=F, row.names=F, col.names=F, sep="\t", append=T)


# get fasta files for the big insertions and deletions to blast the variants

big_del_names <- paste0(">", sapply(strsplit(big_del[,8], ";"), "[[", 3))
big_del_seqs <- big_del[,4]
for(a in 1:length(big_del_names)) {
	if(a == 1) {
		write(big_del_names[a], file="informative_big_deletions.fasta")
	} else {
		write(big_del_names[a], file="informative_big_deletions.fasta", append=T)
	}
	write(big_del_seqs[a], file="informative_big_deletions.fasta", append=T)
}

big_ins_names <- paste0(">", sapply(strsplit(big_ins[,8], ";"), "[[", 3))
big_ins_seqs <- big_ins[,4]
for(a in 1:length(big_ins_names)) {
	if(a == 1) {
		write(big_ins_names[a], file="informative_big_insertions.fasta")
	} else {
		write(big_ins_names[a], file="informative_big_insertions.fasta", append=T)
	}
	write(big_ins_seqs[a], file="informative_big_insertions.fasta", append=T)
}


# which allele is ancestral (using outgroup)

temp_matrix <- big_del
temp_matrix <- temp_matrix[,10:ncol(temp_matrix)]
ancestral <- list()
for(a in 1:nrow(temp_matrix)) {
	a_rep <- temp_matrix[a ,sapply(strsplit(colnames(temp_matrix), "__"), "[[", 1) == "outgroup" ]
	a_rep <- substr(a_rep, 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- unique(c(substr(a_rep, 1, 1), substr(a_rep, 3, 3)))
	if(length(a_rep) > 0) {
		if(length(a_rep) == 1) {
			if(a_rep[1] == 0) {
				ancestral[[a]] <- "ref"
			} else if(a_rep[1] == 1) {
				ancestral[[a]] <- "alt"
			}
		} else {
			ancestral[[a]] <- "het"
		}
	} else {
		ancestral[[a]] <- "missing"
	}
}
ancestral <- unlist(ancestral)
write(ancestral, file="informative_big_deletions.ancestral", ncolumns=1)



temp_matrix <- big_ins
temp_matrix <- temp_matrix[,10:ncol(temp_matrix)]
ancestral <- list()
for(a in 1:nrow(temp_matrix)) {
	a_rep <- temp_matrix[a ,sapply(strsplit(colnames(temp_matrix), "__"), "[[", 1) == "outgroup" ]
	a_rep <- substr(a_rep, 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- unique(c(substr(a_rep, 1, 1), substr(a_rep, 3, 3)))
	if(length(a_rep) > 0) {
		if(length(a_rep) == 1) {
			if(a_rep[1] == 0) {
				ancestral[[a]] <- "ref"
			} else if(a_rep[1] == 1) {
				ancestral[[a]] <- "alt"
			}
		} else {
			ancestral[[a]] <- "het"
		}
	} else {
		ancestral[[a]] <- "missing"
	}
}
ancestral <- unlist(ancestral)
write(ancestral, file="informative_big_insertions.ancestral", ncolumns=1)


temp_matrix <- small_del
temp_matrix <- temp_matrix[,10:ncol(temp_matrix)]
ancestral <- list()
for(a in 1:nrow(temp_matrix)) {
	a_rep <- temp_matrix[a ,sapply(strsplit(colnames(temp_matrix), "__"), "[[", 1) == "outgroup" ]
	a_rep <- substr(a_rep, 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- unique(c(substr(a_rep, 1, 1), substr(a_rep, 3, 3)))
	if(length(a_rep) > 0) {
		if(length(a_rep) == 1) {
			if(a_rep[1] == 0) {
				ancestral[[a]] <- "ref"
			} else if(a_rep[1] == 1) {
				ancestral[[a]] <- "alt"
			}
		} else {
			ancestral[[a]] <- "het"
		}
	} else {
		ancestral[[a]] <- "missing"
	}
}
ancestral <- unlist(ancestral)
write(ancestral, file="informative_small_deletions.ancestral", ncolumns=1)


temp_matrix <- small_ins
temp_matrix <- temp_matrix[,10:ncol(temp_matrix)]
ancestral <- list()
for(a in 1:nrow(temp_matrix)) {
	a_rep <- temp_matrix[a ,sapply(strsplit(colnames(temp_matrix), "__"), "[[", 1) == "outgroup" ]
	a_rep <- substr(a_rep, 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- unique(c(substr(a_rep, 1, 1), substr(a_rep, 3, 3)))
	if(length(a_rep) > 0) {
		if(length(a_rep) == 1) {
			if(a_rep[1] == 0) {
				ancestral[[a]] <- "ref"
			} else if(a_rep[1] == 1) {
				ancestral[[a]] <- "alt"
			}
		} else {
			ancestral[[a]] <- "het"
		}
	} else {
		ancestral[[a]] <- "missing"
	}
}
ancestral <- unlist(ancestral)
write(ancestral, file="informative_small_insertions.ancestral", ncolumns=1)



#########################
######### check snps
#########################

x <- read_vcf("snps_filtered_mac10.recode.vcf")

genotypes <- x[,10:ncol(x)]

allo_virens <- genotypes[ ,sapply(strsplit(colnames(genotypes), "__"), "[[", 1) == "allo" & sapply(strsplit(colnames(genotypes), "__"), "[[", 2) == "virens"]

allo_sordid <- genotypes[ ,sapply(strsplit(colnames(genotypes), "__"), "[[", 1) == "allo" & sapply(strsplit(colnames(genotypes), "__"), "[[", 2) == "sordidulus"]

# loop for each variant
keep <- list()
for(a in 1:nrow(allo_sordid)) {
	a_rep <- substr(allo_sordid[a,], 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- c(substr(a_rep, 1, 1), substr(a_rep, 3, 3))
	if(length(unique(a_rep)) == 1) {
		a_rep2 <- substr(allo_virens[a,], 1, 3)
		a_rep2 <- a_rep2[a_rep2 != "./."]
		a_rep2 <- c(substr(a_rep2, 1, 1), substr(a_rep2, 3, 3))
		if(length(unique(a_rep2)) == 1 & a_rep[1] != a_rep2[1]) {
			keep[[a]] <- a
		}
	}
}
keep <- unlist(keep)

x2 <- x[keep,]

write.table(x2, file="informative_snps.vcf", quote=F, row.names=F, col.names=F, sep="\t", append=T)


temp_matrix <- x2
temp_matrix <- temp_matrix[,10:ncol(temp_matrix)]
ancestral <- list()
for(a in 1:nrow(temp_matrix)) {
	a_rep <- temp_matrix[a ,sapply(strsplit(colnames(temp_matrix), "__"), "[[", 1) == "outgroup" ]
	a_rep <- substr(a_rep, 1, 3)
	a_rep <- a_rep[a_rep != "./."]
	a_rep <- unique(c(substr(a_rep, 1, 1), substr(a_rep, 3, 3)))
	if(length(a_rep) > 0) {
		if(length(a_rep) == 1) {
			if(a_rep[1] == 0) {
				ancestral[[a]] <- "ref"
			} else if(a_rep[1] == 1) {
				ancestral[[a]] <- "alt"
			}
		} else {
			ancestral[[a]] <- "het"
		}
	} else {
		ancestral[[a]] <- "missing"
	}
}
ancestral <- unlist(ancestral)
write(ancestral, file="informative_snps.ancestral", ncolumns=1)








