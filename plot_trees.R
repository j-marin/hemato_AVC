library(ape)
library(ggtree)
library(gclus)
library(RColorBrewer)
require(pheatmap)
library(magrittr)
library(Pigengene)

## functions
list_snp <- function(samples) {
	snp.d <- integer(0)
		for (i in 1:length(samples)) {
			snps <- snp$SNP_ID[which(snp$Sample == samples[i])]
			snp.d <- c(snp.d, list(snps))
			}	
	return(snp.d)
	}

dist_mat <- function(samples, snp.list) {
	cb <- combn(length(samples), 2)
	mat <- matrix(NA, nrow=length(samples), ncol=length(samples))

	for (i in 1:dim(cb)[2]) {
		d1 <- length(setdiff(snp.list[[cb[1,i]]], snp.list[[cb[2,i]]]))
		d2 <- length(setdiff(snp.list[[cb[2,i]]], snp.list[[cb[1,i]]]))
		mat[cb[2,i], cb[1,i]] <- sum(d1,d2)
		}
	mat <- as.matrix(Matrix::forceSymmetric(mat, uplo="L"))
	mat[is.na(mat)] <- 0
	rownames(mat) <- samples
	colnames(mat) <- samples
	
	return(mat)
	}

## open tables
recap <- read.table("../recap.txt", header=TRUE)
head(recap)
recap$Prélèvement[recap$Prélèvement == "ECBU"] <- "urine"
recap$Prélèvement[recap$Prélèvement == "hemoc"] <- "blood"
recap$Prélèvement[recap$Prélèvement == "selles"] <- "stool"

strains <- read.table("../recapSeqFinal.txt", header=TRUE)

snp <- read.table("../snippy/v2/recap_snp_patients_tab.txt", header=TRUE)
head(snp)
snp$Sample <- tolower(snp$Sample)

## patient 2 (A)
labs <- strains$Sample[strains$Patient == 2]

recap2 <- recap[match(tolower(labs), recap$Nom),]
samples2 <- recap2$Nom

snp2 <- list_snp(tolower(labs))
mat2 <- dist_mat(samples2, snp2)

t2 <- nj (as.dist(mat2)) 

tip_data <- recap[match(tolower(labs), recap$Nom), c(4,3)]
colnames(tip_data) <- c("Strain", "Site")
t2$tip.label <- tolower(t2$tip.label)

p <- ggtree(t2, layout="equal_angle")
p <- p %<+% tip_data + geom_tippoint(aes(color= Site), size=3) + scale_color_manual(values=c("firebrick1", "chocolate3","gold1")) + geom_treescale(offset=0.15, width=10)

# save figure
pdf(file = "../trees/tree_patientA.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

plot(p)

dev.off()


## patient 16 (C)
labs <- strains$Sample[strains$Patient == 16]

recap2 <- recap[match(tolower(labs), tolower(recap$Nom)),]
samples2 <- recap2$Nom

snp2 <- list_snp(tolower(labs))
mat2 <- dist_mat(samples2, snp2)

t2 <- nj (as.dist(mat2)) 

tip_data <- recap[match(tolower(labs), tolower(recap$Nom)), c(4,3)]
colnames(tip_data) <- c("Strain", "Site")

p <- ggtree(t2, layout="equal_angle")
p <- p %<+% tip_data + geom_tippoint(aes(color= Site), size=3) + scale_color_manual(values=c("firebrick1", "chocolate3","gold1")) + geom_treescale(offset=0.15, width=10)

# save figure
pdf(file = "../trees/tree_patientC.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

plot(p)

dev.off()


## patient 92 (B)
labs <- strains$Sample[strains$Patient == "6"]

recap2 <- recap[match(tolower(labs), tolower(recap$Nom)),]
samples2 <- recap2$Nom

snp2 <- list_snp(tolower(labs))
mat2 <- dist_mat(samples2, snp2)

t2 <- nj (as.dist(mat2)) 

tip_data <- recap[match(tolower(labs), tolower(recap$Nom)), c(4,3)]
colnames(tip_data) <- c("Strain", "Site")

p <- ggtree(t2, layout="equal_angle")
p <- p %<+% tip_data + geom_tippoint(aes(color= Site), size=3) + scale_color_manual(values=c("firebrick1", "chocolate3")) + geom_treescale(offset=0.15, width=10)

# save figure
pdf(file = "../trees/tree_patientB.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches

plot(p)

dev.off()

# sub-group
plot(t2)
nodelabels()

t3 <- root(t2, node=27)
plot(t3)
nodelabels()

t2.zoom <- extract.clade(t3, 25)
tip_data.zoom <- tip_data[match(t2.zoom$tip.label, tip_data$Strain),]

p <- ggtree(t2.zoom, layout="equal_angle")
p <- p %<+% tip_data.zoom + geom_tippoint(aes(color= Site), size=3) + scale_color_manual(values=c("firebrick1", "chocolate3")) + geom_treescale(offset=0.15, width=10)

# save figure
pdf(file = "../trees/tree_patientB_zoom.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches

plot(p)

dev.off()


