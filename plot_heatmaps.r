library(ape)
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

## order by time and site
source("../../scripts/pheatmap.type.R")

## open tables
recap <- read.table("../../recap.txt", header=TRUE)
head(recap) # to complete with petanc results
recap$Prélèvement[recap$Prélèvement=="ECBU"] <- "urine"
recap$Prélèvement[recap$Prélèvement=="selles"] <- "stool"
recap$Prélèvement[recap$Prélèvement=="hemoc"] <- "blood"

snp <- read.table("recap_snp_patients_tab.txt", header=TRUE)
head(snp)
snp$Sample <- tolower(snp$Sample)


## patient 9 clade 2 (B)
aln <- read.dna("patient9_cl2/clean.full.aln", format="fasta")
labs <- gsub("[[:blank:]]", "", labels(aln))

recap9 <- recap[match(tolower(labs[1:17]), recap$Nom),] # 9B4aeC1 ??
samples9 <- recap9$Nom

snp9 <- list_snp(samples9)
mat9 <- dist_mat(samples9, snp9)
# 42 snp ref / ref

sample_ann <- data.frame(Site = recap$Prélèvement[match(samples9, recap$Nom)], Date = recap$Date[match(samples9, recap$Nom)], Phylogroup = recap$Phylogroup[match(samples9, recap$Nom)], ST = as.character(recap$ST[match(samples9, recap$Nom)]), Serotype = recap$Serotype[match(samples9, recap$Nom)], fimH = recap$fimH[match(samples9, recap$Nom)])
row.names(sample_ann) <- samples9

# choose colors
Site <- c("firebrick1", "blue")
names(Site) <- levels(factor(sample_ann$Site))
Date <- brewer.pal(n = 3, name = "Greens") [1]
names(Date) <- levels(factor(sample_ann$Date))
Phylogroup <- "forestgreen"
names(Phylogroup) <- levels(factor(sample_ann$Phylogroup))
ST <- brewer.pal(n = 3, name = "Set2") [1]
names(ST) <- levels(factor(sample_ann$ST))
Serotype <- brewer.pal(n = 7, name = "Set2") [-c(1:3)] [1]
names(Serotype) <- levels(factor(sample_ann$Serotype))
fimH <- brewer.pal(n = 3, name = "Set3") [2]
names(fimH) <- levels(factor(sample_ann$fimH))

my_colour <- list(Site = Site, Date = Date, Phylogroup = Phylogroup, ST = ST, Serotype = Serotype, fimH = fimH)

pheatmap(mat9, show_colnames=T, labels_col=rep("", length(samples9)), annotation_row=sample_ann, clustering_method="complete", color=brewer.pal(n = 9, name = "Blues"), annotation_colors=my_colour, angle_col=45)

# get dendogram order and annotation table ordered
my_heatmap <- pheatmap(mat9, clustering_method="complete", silent = TRUE)
s <- my_heatmap$tree_row$labels[my_heatmap$tree_row$order]
sample_ann_o <- sample_ann[match(s,rownames(sample_ann)),]

write.table(sample_ann_o, "annotation_p9_cl2.txt")

## order by time and site
source("../../scripts/pheatmap.type.R")

sample_ann <- sample_ann[order(sample_ann[, 1]),]

snp2 <- list_snp(tolower(rownames(sample_ann)))
mat2 <- dist_mat(rownames(sample_ann), snp2)

makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

cutoff.distance <- 14 
cols <- makeColorRampPalette(c("white", "dodgerblue4",    # distances 0 to 3 colored from white to red
                               "gray", "black"), # distances 3 to max(distmat) colored from green to black
                             cutoff.distance / max(mat2),
                             1000)
                             
# save figure
pdf(file = "heatmap_patientB_final.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

pheatmap.type(mat2, annRow= sample_ann, clustering_method="complete", type=colnames(sample_ann)[c(1)], color=cols, doTranspose=FALSE,annotation_colors=my_colour, annotation_col=sample_ann[1], main="Patient B") 

dev.off()


## patient 2 (A)
aln <- read.dna("patient2/clean.full.aln", format="fasta")
labs <- gsub("[[:blank:]]", "", labels(aln))
labs <- labs[-which(labs == "Reference")] # remove Ref
labs <- labs[-which(labs == "2A3aec1")] # remove Ref
labs <- labs[-which(labs == "2b2c2")] # remove Ref
labs <- labs[-which(labs == "3A1Aec1")] # remove Ref

recap2 <- recap[match(tolower(labs), recap$Nom),]
samples2 <- recap2$Nom

snp2 <- list_snp(tolower(labs))
mat2 <- dist_mat(samples2, snp2)

sample_ann <- data.frame(Site = recap$Prélèvement[match(samples2, recap$Nom)], Date = recap$Date[match(samples2, recap$Nom)], Phylogroup = recap$Phylogroup[match(samples2, recap$Nom)], ST = as.character(recap$ST[match(samples2, recap$Nom)]), Serotype = recap$Serotype[match(samples2, recap$Nom)], fimH = recap$fimH[match(samples2, recap$Nom)])
row.names(sample_ann) <- samples2
sample_ann$Date <- factor(sample_ann$Date, levels = c("10/03/2015","11/03/2015","13/03/2015","27/03/2015","02/04/2015","03/04/2015","22/04/2015","04/05/2015","05/05/2015"))

# choose colors
Site <- c("gold1", "firebrick1", "chocolate3")
names(Site) <- levels(factor(sample_ann$Site)) [c(3,1,2)]
Date <- brewer.pal(n = 9, name = "Greens")
names(Date) <- levels(sample_ann$Date)
Phylogroup <- "dodgerblue"
names(Phylogroup) <- levels(factor(sample_ann$Phylogroup))
ST <- brewer.pal(n = 3, name = "Set2") [c(1)]
names(ST) <- levels(factor(sample_ann$ST))
Serotype <- brewer.pal(n = 7, name = "Set2") [-c(1:2)] [c(1)]
names(Serotype) <- levels(factor(sample_ann$Serotype))
fimH <- brewer.pal(n = 3, name = "Set3") [2]
names(fimH) <- levels(factor(sample_ann$fimH))

my_colour <- list(Site = Site, Date = Date, Phylogroup = Phylogroup, ST = ST, Serotype = Serotype, fimH = fimH)


pheatmap(mat2, show_colnames=T, labels_col=rep("", length(samples2)), annotation_row=sample_ann, clustering_method="complete", color=brewer.pal(n = 9, name = "Blues"), annotation_colors=my_colour, angle_col=45)

# get dendogram order and annotation table ordered
my_heatmap <- pheatmap(mat2, clustering_method="complete", silent = TRUE)
s <- my_heatmap$tree_row$labels[my_heatmap$tree_row$order]
sample_ann_o <- sample_ann[match(s,rownames(sample_ann)),]

write.table(sample_ann_o, "annotation_p2.txt")


## order by time and site
source("../../scripts/pheatmap.type.R")

sample_ann <- sample_ann[order(sample_ann[, 1]),]

snp2 <- list_snp(tolower(rownames(sample_ann)))
mat2 <- dist_mat(rownames(sample_ann), snp2)

cutoff.distance <- 14
cols <- makeColorRampPalette(c("white", "dodgerblue4",    # distances 0 to 3 colored from white to red
                               "gray", "black"), # distances 3 to max(distmat) colored from green to black
                             cutoff.distance / 14,
                             1000)
# save figure
pdf(file = "heatmap_patientA_final.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

pheatmap.type(mat2, annRow= sample_ann, clustering_method="complete", type=colnames(sample_ann)[c(1)], color=cols, doTranspose=FALSE,annotation_colors=my_colour, annotation_col=sample_ann[1], main="Patient A") 

dev.off()


## patient 16 (C)
aln <- read.dna("patient16/clean.full.aln", format="fasta")
labs <- gsub("[[:blank:]]", "", labels(aln))
labs <- labs[!labs %in% c("Reference", "EC30FP", "EC30tris")]

# replace patient16 names
r16 <- recap[recap$Patient == "16",]
for (i in 1:dim(r16)[1]) {
	labs[which(labs== r16$correspondance[i])] <- tolower(r16$Nom[i])
	}

recap2 <- recap[match(tolower(labs), tolower(recap$Nom)),]
samples2 <- recap2$Nom

snp2 <- list_snp(tolower(labs))
mat2 <- dist_mat(samples2, snp2)

sample_ann <- data.frame(Site = recap$Prélèvement[match(samples2, recap$Nom)], Date = recap$Date[match(samples2, recap$Nom)], Phylogroup = recap$Phylogroup[match(samples2, recap$Nom)], ST = as.character(recap$ST[match(samples2, recap$Nom)]), Serotype = recap$Serotype[match(samples2, recap$Nom)], fimH = recap$fimH[match(samples2, recap$Nom)])
row.names(sample_ann) <- samples2
sample_ann$Date <- factor(sample_ann$Date, levels = c("03/04/2017","04/04/2017"))

# choose colors
Site <- c("gold1", "firebrick1", "chocolate3")
names(Site) <- levels(factor(sample_ann$Site)) [c(3,1,2)]
names(Site) <- levels(factor(sample_ann$Site)) [c(2,3,1)]
Date <- brewer.pal(n = 5, name = "Greens") [c(1,3)]
names(Date) <- levels(sample_ann$Date)
Phylogroup <- "dodgerblue"
names(Phylogroup) <- levels(factor(sample_ann$Phylogroup))
ST <- brewer.pal(n = 3, name = "Set2") [c(1)]
names(ST) <- levels(factor(sample_ann$ST))
Serotype <- brewer.pal(n = 7, name = "Set2") [-c(1:2)] [c(1)]
names(Serotype) <- levels(factor(sample_ann$Serotype))
fimH <- brewer.pal(n = 3, name = "Set3") [2:3][1]
names(fimH) <- levels(factor(sample_ann$fimH))

my_colour <- list(Site = Site, Date = Date, Phylogroup = Phylogroup, ST = ST, Serotype = Serotype, fimH = fimH)


pheatmap(mat2, show_colnames=T, labels_col=rep("", length(samples2)), annotation_row=sample_ann, clustering_method="complete", color=brewer.pal(n = 9, name = "Blues"), annotation_colors=my_colour, angle_col=45)

# get dendogram order and annotation table ordered
my_heatmap <- pheatmap(mat2, clustering_method="complete", silent = TRUE)
s <- my_heatmap$tree_row$labels[my_heatmap$tree_row$order]
sample_ann_o <- sample_ann[match(s,rownames(sample_ann)),]

write.table(sample_ann_o, "annotation_p16.txt")


## order by time and site
source("../../scripts/pheatmap.type.R")

sample_ann <- sample_ann[order(sample_ann[, 1]),]

snp2 <- list_snp(tolower(rownames(sample_ann)))
mat2 <- dist_mat(rownames(sample_ann), snp2)


# save figure
pdf(file = "heatmap_patientC_final.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

pheatmap.type(mat2, annRow= sample_ann, clustering_method="complete", type=colnames(sample_ann)[c(1)], color=cols, doTranspose=FALSE,annotation_colors=my_colour, annotation_col=sample_ann[1], main="Patient C") 

dev.off()


