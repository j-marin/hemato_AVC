library(stringr)
library(ggplot2)
library(ape)
library(seqinr)

a <- read.table("../recapSeqFinal.txt", header=TRUE)
head(a)
r <- read.table("../recap.txt", header=TRUE)
head(r)

### Resistance
# open data
seqR <- read.fasta("RefSeq_resistance.fasta", as.string=FALSE)
g <- names(seqR)
l <- lapply(seqR, length)

files <- list.files("resistance_mapping")

f <- unlist(lapply(strsplit(files, "-"), `[`, 1))
f <- unlist(lapply(strsplit(f, "_"), `[`, 1))
setdiff(a$Sample, f)

matR <- matrix(0, ncol=length(g), nrow=length(f))
for (i in 1:length(files)) {
	d <- read.table(paste("resistance_mapping/", files[i], sep=""))
	for (j in 1:length(g)) {
		p <- length(which(d[which(d[,1] == g[j]),3] >= 1)) / unlist(l[j])
		if (p>0.8) {
			matR[i,j] <- 1
		}
	}
print(i)}

colnames(matR) <- g
rownames(matR) <- f
write.table(matR, quote=FALSE, "resistance_pres_abs_mapping_d1.txt")

## make Presence / Absence plots
#prep data
Sample <- rep(f, length(g))
Site <- a$Site[match(Sample, a$Sample)]
Patient <- a$Patient[match(Sample, a$Sample)]
Gene <- rep(g, each = length(f))
presence <- as.vector(matR)

tab <- data.frame(Patient = Patient, Sample = Sample, Site = Site, Gene = Gene, Presence = presence)
tab <- tab[-which(presence == 0),]

#plot
tab$Patient[tab$Patient == 2] <- "A"
tab$Patient[tab$Patient == 6] <- "B"
tab$Patient[tab$Patient == 16] <- "C"
tab$Patient[tab$Patient == 19] <- "D"
tab$Patient[tab$Patient == 22] <- "D"

tab$Site <- as.factor(tab$Site)
tab$Site <- factor(tab$Site, levels = c("stool", "blood", "urine"))

tab <- tab[order(tab$Gene),] # order by Gene
tab <- tab[order(tab$Sample),] # order by Sample
tab <- tab[order(tab$Site),] # order by Site
tab <- tab[order(tab$Patient),] # order by Patient

tab$Gene <- as.factor(tab$Gene)
tab$Gene <- factor(tab$Gene, levels = unique(tab$Gene))

tab$Sample <- as.factor(tab$Sample)
tab$Sample <- factor(tab$Sample, levels = unique(tab$Sample))

ggplot(tab, aes(x=Gene, y=Sample, fill=Site)) + 
      geom_tile() + 
      facet_grid(rows=vars(Patient), scales="free", space="free") + 
      theme(panel.background=element_rect(colour="white", fill="white"), panel.border=element_rect(colour="gray", fill=NA), axis.text.x=element_text(angle=45, vjust=0.9, hjust=0.9), axis.text.y=element_text(size=5)) +
      scale_fill_manual(values=c("chocolate3", "firebrick1", "gold1")) +
      ggtitle("Resistance genes")


### Virulence
# open data
seqV <- read.fasta("RefSeq_virulence.fasta", as.string=FALSE)
g <- names(seqV)
l <- lapply(seqV, length)

files <- list.files("virulence_mapping")

f <- unlist(lapply(strsplit(files, "-"), `[`, 1))
f <- unlist(lapply(strsplit(f, "_"), `[`, 1))
setdiff(a$Sample, f)

matV <- matrix(0, ncol=length(g), nrow=length(f))
for (i in 1:length(files)) {
	d <- read.table(paste("virulence_mapping/", files[i], sep=""))
	for (j in 1:length(g)) {
		p <- length(which(d[which(d[,1] == g[j]),3] >= 1)) / unlist(l[j])
		if (p>0.8) {
			matV[i,j] <- 1
		}
	}
print(i)}

colnames(matV) <- g
rownames(matV) <- f
write.table(matV, quote=FALSE, "virulence_pres_abs_mapping_d1.txt")

## make Presence / Absence plots
#prep data
Sample <- rep(f, length(g))
Site <- a$Site[match(Sample, a$Sample)]
Patient <- a$Patient[match(Sample, a$Sample)]
Gene <- rep(g, each = length(f))
presence <- as.vector(matV)

tab <- data.frame(Patient = Patient, Sample = Sample, Site = Site, Gene = Gene, Presence = presence)
tab <- tab[-which(presence == 0),]

#plot
tab$Patient[tab$Patient == 2] <- "A"
tab$Patient[tab$Patient == 6] <- "B"
tab$Patient[tab$Patient == 16] <- "C"
tab$Patient[tab$Patient == 19] <- "D"
tab$Patient[tab$Patient == 22] <- "D"

tab$Site <- as.factor(tab$Site)
tab$Site <- factor(tab$Site, levels = c("stool", "blood", "urine"))

tab <- tab[order(tab$Gene),] # order by Gene
tab <- tab[order(tab$Sample),] # order by Sample
tab <- tab[order(tab$Site),] # order by Site
tab <- tab[order(tab$Patient),] # order by Patient

tab$Gene <- as.factor(tab$Gene)
tab$Gene <- factor(tab$Gene, levels = unique(tab$Gene))

tab$Sample <- as.factor(tab$Sample)
tab$Sample <- factor(tab$Sample, levels = unique(tab$Sample))

ggplot(tab, aes(x=Gene, y=Sample, fill=Site)) + 
      geom_tile() + 
      facet_grid(rows=vars(Patient), scales="free", space="free") + 
      theme(panel.background=element_rect(colour="white", fill="white"), panel.border=element_rect(colour="gray", fill=NA), axis.text.x=element_text(angle=45, vjust=0.9, hjust=0.9), axis.text.y=element_text(size=5)) +
      scale_fill_manual(values=c("chocolate3", "firebrick1", "gold1")) +
      ggtitle("Virulence genes")



