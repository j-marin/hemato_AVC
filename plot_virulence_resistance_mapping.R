library(stringr)
library(ggplot2)

a <- read.table("../recapSeqFinal.txt", header=TRUE)
head(a)
r <- read.table("../recap.txt", header=TRUE)
head(r)

### Resistance
# open data
res1 <- read.table("v1/Resistance_genes.txt", header=TRUE)
head(res1)
res2 <- read.table("v2/Resistance_genes.txt", header=TRUE)
head(res2)
res3 <- read.table("v2.2/Resistance_genes.txt", header=TRUE)
head(res3)

# harmonize names
res2$FILE <- sapply(str_split(res2$FILE, "-"), `[`, 1)
res2$FILE <- r$Nom[match(tolower(res2$FILE), tolower(r$ngs))]

res3$FILE <- sapply(str_split(res3$FILE, "_"), `[`, 1)
res3$FILE <- r$Nom[match(tolower(res3$FILE), tolower(r$correspondance))]

# combine res1, res2 and res3
res1 <- res1 [!startsWith(res1$FILE, '22'), ]
res1 <- res1 [!startsWith(res1$FILE, '9'), ]

res <- rbind(res1,res2,res3)
res$FILE <- tolower(res$FILE)

res <- res[res$FILE %in% a$Sample,]
s <- unique(res$FILE) # 94

# make presence/absence table
g <- unique(c(res$PRODUCT)) # 21

res_tab <- matrix(0, ncol=length(g), nrow=length(s))
for (i in 1:length(s)) {
		res_tab[i,!is.na(match(g, res$PRODUCT[which(res$FILE == s[i])]))] <- 1	
	}
colnames(res_tab) <- g
rownames(res_tab) <- s

#write.table(res_tab, "ResistanceGenesPresAbs.txt")

#prep data
Site <- a$Site[match(res$FILE, a$Sample)]
Patient <- a$Patient[match(res$FILE, a$Sample)]
tab <- data.frame(Patient = Patient, Sample = res$FILE, Site = Site, Gene=res$PRODUCT, Contig=res$SEQUENCE, Pos=paste(res$START, res$END, sep="-"))
dim(tab)

# add localisation
loc <- integer(0)
for (i in 1:dim(tab)[1]) {
	loc[i] <- unique(res$PLASCOPES[res$FILE == tab$Sample[i] & res$PRODUCT == tab$Gene[i]])
	}
tab["loc"] <- loc
write.table(tab, "petanc_resistance_table.txt", quote=FALSE, row.names=FALSE, sep="	")

# Now mapped genes (coverage > 80%, depth > 10)
recov <- read.table("resistance_pres_abs_mapping_d1.txt", header=TRUE)
colnames(recov) = c("dfrA5","mdf(A)","aph(4)-Ia","aac(3)-IVa","tet(B)","sul2","aph(3'')-Ib","aph(6)-Id","aph(3')-Ia","dfrA1","tet(A)","blaTEM-35","blaTEM-1B","dfrA14","blaTEM-1C","qnrS1","blaCTX-M-15","dfrA17","aadA5","sul1","mph(A)")

tab2 <- data.frame(Patient=NA, Sample=NA, Site=NA, Gene=NA)
p <- unique(Patient)
for (i in 1:length(p)) {
	s <- unique(tab$Sample[which(Patient == p[i])])
	g <- unique(tab$Gene[which(Patient == p[i])])
	for (k in 1:length(s)) {
		v <- recov[which(rownames(recov) == s[k]), match(g, colnames(recov))]
		v <- v[which(v > 0)]
		tabx <- data.frame(Patient=p[i], Sample=s[k], Site=a$Site[match(s[k], a$Sample)], Gene=names(v))
		tab2 <- rbind(tab2, tabx)		
		} 	
	}
tab2 <- tab2[-1,]
head(tab2)

# make plot
tab2$Patient[tab2$Patient == 2] <- "A"
tab2$Patient[tab2$Patient == 6] <- "B"
tab2$Patient[tab2$Patient == 16] <- "C"
tab2$Patient[tab2$Patient == 19] <- "D"
tab2$Patient[tab2$Patient == 22] <- "D"

# make data to plot
g <- unique(tab2$Gene)
sample_order <- tolower(c("2c1c5","2c1c4","2c1c3","2c1c2","2b1c2","2c1c1","2b1c1","2b1c4","2ana3","2a1ae1","2ae3","2b3c5","2b3c4","2b3c3","2b3c2","2b3c1","2b2c5","2b2c4","2b2c3","2b2c1","6a3c1","3a2anac1","3a2aec1","2a4aec1","3a1anac1","6a2c1","2b1c5","4a1aec1","6a1c1","6a4c1","9a5aec1","9a5aec2","9b4c10","9b4c2","9b4c5","9b4c4","9b4c6","9b4c3","9b4c9","9b4c8","9a4anac1","9b4c7","9a4anac2","9a4anac3","9a4aec2","9a4aec3","9b4c1","19anaeroC5","19C1C10","19anaeroC4","19aero1C2","19anaeroC1","19anaeroC2","19C1C5","19aero1C1","19anaeroC3","19B1C10","19C1C4","19C1C9","19B1C8","19B1C4","19B1C6","19B1C9","19C1C1","19C1C6","19C1C7","19C1C8","19aero1C5","19B1C3","19C1C3","19B1C2","19B1C1","19C1C2","19aero1C3","19aero1C4","19B1C5","19B1C7","22a1anac10","22b1c6","22a1anac3","22a1anac5","22b1c3","22b1c5","22b1c1","22b1c4","22b1c110","22b1c2","22b1c8","22a1anac4","22a1anac6","22a1anac2","22a1anac8","22a1anac7","22a1anac9"))
gene <- rep(g, each=length(sample_order))
sample_order <- rep(sample_order, length(g))
patient <- tab2$Patient[match(sample_order, tab2$Sample)]

plot_data <- data.frame(patient=patient, sample=sample_order, gene=gene)

values <- sapply(1:dim(plot_data)[1], function (x) dim(tab2[tab2$Sample == plot_data$sample[x] & tab2$Gene == plot_data$gene[x],])[1])
plot_data["values"] <- values

# add sites
sites <- plot_data[1:94,]
sites$gene <- "Site"
sites$values <- tab2$Site[match(sites$sample, tab2$Sample)]

# modify levels
plot_data <- rbind(plot_data, sites)
plot_data$values <- as.factor(plot_data$values)
levels(plot_data$values) <- c(levels(plot_data$values ), "") # add artificial levels to split the legend into 2 columns
levels(plot_data$values) <- c(levels(plot_data$values ), " ") 
plot_data$values <- factor(plot_data$values, levels(plot_data$values)[c(1,2,6,7,3:5)])

plot_data$gene <- as.factor(plot_data$gene)
levels_gene <- c("Site","aph(4)-Ia","aac(3)-IVa","qnrS1","blaCTX-M-15","blaTEM-1B","aph(6)-Id","aph(3'')-Ib","sul2","mph(A)","sul1","aadA5","dfrA17","tet(A)","dfrA5","mdf(A)","tet(B)","dfrA1","aph(3')-Ia","blaTEM-35","dfrA14","blaTEM-1C")
plot_data$gene <- factor(plot_data$gene, levels_gene)

levels_sample <-tolower(c("2c1c5","2c1c4","2c1c3","2c1c2","2b1c2","2c1c1","2b1c1","2b1c4","2ana3","2a1ae1","2ae3","2b3c5","2b3c4","2b3c3","2b3c2","2b3c1","2b2c5","2b2c4","2b2c3","2b2c1","6a3c1","3a2anac1","3a2aec1","2a4aec1","3a1anac1","6a2c1","2b1c5","4a1aec1","6a1c1","6a4c1","9a5aec1","9a5aec2","9b4c10","9b4c2","9b4c5","9b4c4","9b4c6","9b4c3","9b4c9","9b4c8","9a4anac1","9b4c7","9a4anac2","9a4anac3","9a4aec2","9a4aec3","9b4c1","19anaeroC5","19C1C10","19anaeroC4","19aero1C2","19anaeroC1","19anaeroC2","19C1C5","19aero1C1","19anaeroC3","19B1C10","19C1C4","19C1C9","19B1C8","19B1C4","19B1C6","19B1C9","19C1C1","19C1C6","19C1C7","19C1C8","19aero1C5","19B1C3","19C1C3","19B1C2","19B1C1","19C1C2","19aero1C3","19aero1C4","19B1C5","19B1C7","22a1anac10","22b1c6","22a1anac3","22a1anac5","22b1c3","22b1c5","22b1c1","22b1c4","22b1c110","22b1c2","22b1c8","22a1anac4","22a1anac6","22a1anac2","22a1anac8","22a1anac7","22a1anac9"))
plot_data$sample <- factor(plot_data$sample, levels_sample)

# plot 
colors <- c("firebrick1", "blue", "gold1")
sizes <- plot_data$loc
sizes <- gsub("chromosome", 0.5, sizes)
sizes <- gsub("plasmid", 0.5, sizes)
sizes <- gsub("unclassified", 0, sizes)
sizes <- as.numeric(sizes)

ggplot(plot_data[plot_data$patient != "D",], aes(gene, sample)) + 
      geom_tile(aes(fill=values), show.legend=TRUE) + 
      facet_grid(rows=vars(patient), scales="free", space="free") + 
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=0.9, size=5), axis.text.y=element_text(size=5), legend.key = element_blank()) +
      scale_fill_manual(values = c("gray90", "gray42","white","white", colors), "Gene Site", drop=FALSE) + 
      #scale_color_manual(values = c("black", "white","pink","pink"), "Localisation") + 
      ggtitle("Resistance genes") +
      guides (fill = guide_legend(ncol=2))


######## ############ ####### #######
### Virulence
# open data
vir1 <- read.table("v1/Virulence_genes.txt", header=TRUE)
head(vir1)
vir2 <- read.table("v2/Virulence_genes.txt", header=TRUE)
head(vir2)
vir3 <- read.table("v2.2/Virulence_genes.txt", header=TRUE)
head(vir3)

# harmonize names
vir2$FILE <- sapply(str_split(vir2$FILE, "-"), `[`, 1)
vir2$FILE <- r$Nom[match(tolower(vir2$FILE), tolower(r$ngs))]

vir3$FILE <- sapply(str_split(vir3$FILE, "_"), `[`, 1)
vir3$FILE <- r$Nom[match(tolower(vir3$FILE), tolower(r$correspondance))]

# combine vir1, vir2 and vir3
vir1 <- vir1 [!startsWith(vir1$FILE, '22'), ]
vir1 <- vir1 [!startsWith(vir1$FILE, '9'), ]

vir <- rbind(vir1,vir2,vir3)
vir$FILE <- tolower(vir$FILE)

vir <- vir[vir$FILE %in% a$Sample,]
s <- unique(vir$FILE) # 94

# make presence/absence table
g <- unique(c(vir$GENE)) # 146

vir_tab <- matrix(0, ncol=length(g), nrow=length(s))
for (i in 1:length(s)) {
		vir_tab[i,!is.na(match(g, vir$GENE[which(vir$FILE == s[i])]))] <- 1	
	}
colnames(vir_tab) <- g
rownames(vir_tab) <- s

#write.table(vir_tab, "VirulenceGenesPresAbs.txt")

#prep data
Site <- a$Site[match(vir$FILE, a$Sample)]
Patient <- a$Patient[match(vir$FILE, a$Sample)]
tab <- data.frame(Patient = Patient, Sample = vir$FILE, Site = Site, Gene=vir$GENE, Contig=vir$SEQUENCE, Pos=paste(vir$START, vir$END, sep="-"))
dim(tab)

# add localisation
loc <- integer(0)
for (i in 1:dim(tab)[1]) {
	loc[i] <- unique(vir$PLASCOPES[vir$FILE == tab$Sample[i] & vir$GENE == tab$Gene[i]])
	}
tab["loc"] <- loc
write.table(tab, "petanc_virulence_table.txt", quote=FALSE, row.names=FALSE, sep="	")


# add mapped genes (coverage > 80%, depth > 10)
recov <- read.table("virulence_pres_abs_mapping.txt", header=TRUE)
colnames(recov) = c("hra","espX5","espX4","espL4","aslA","entA","entB","entE","entC","fepB","entS","fepD","fepG","fepC","entF","fes","fepA","entD","malX","espL1","espR4","yagV/ecpE","yagW/ecpD","yagX/ecpC","yagY/ecpB","yagZ/ecpA","ykgK/ecpR","fdeC","fimH","fimG","fimF","fimD","fimC","fimI","fimA","fimE","fimB","espY1","espX1","cvaC","mchF_14","iroN_6","iroE","iroD","iroC","iroB","iss_12","ompA","Episomal_ompT","hlyF","cba_6","cma_16","iutA","iucD","iucC","iucB","iucA","traT","gad_17","gad_27","etsC","sepA_2","espT","nleC_5","escO","escN","escV","eae_45","nleB_3","escJ","escI","tir_34","espA_22","nleA_3","cesD","cesD2","espB_13","map","espG","espH","escS","paa","espF_2","cesF","sepQ/escQ","escP","cesL","cesAB","escE","nleH2","nleF","espD_E2348/69","escF","escG","sepD","escC","escT","etgA","escL","escD","cesT","escU","escR","espJ_1","sepL","gad_37","gad_15","lpfA_8","gad_42","capU_5","iss_13","gad_13","gad_11","gad_20","gad_31","daaF","afaA","nfaE_4","afaC-I","afaD","draP","afaE-I","Chromosomal_ompT","kpsM","fyuA","irp2","chuS","chuA","chuT","chuW","chuX","chuY","chuU","chuV","gad_57","iha_11","iss_11","kpsD","kpsE","senB_3","papI","papB","sat_3","gad_64","papX","usp")

tab2 <- data.frame(Patient=NA, Sample=NA, Site=NA, Gene=NA)
p <- unique(Patient)
for (i in 1:length(p)) {
	s <- unique(tab$Sample[which(Patient == p[i])])
	g <- unique(tab$Gene[which(Patient == p[i])])
	for (k in 1:length(s)) {
		v <- recov[which(rownames(recov) == s[k]), match(g, colnames(recov))]
		v <- v[which(v > 0)]
		tabx <- data.frame(Patient=p[i], Sample=s[k], Site=a$Site[match(s[k], a$Sample)], Gene=names(v))
		tab2 <- rbind(tab2, tabx)		
		} 	
	}
tab2 <- tab2[-1,]
head(tab2)

tab2$Patient[tab2$Patient == 2] <- "A"
tab2$Patient[tab2$Patient == 6] <- "B"
tab2$Patient[tab2$Patient == 16] <- "C"
tab2$Patient[tab2$Patient == 19] <- "D"
tab2$Patient[tab2$Patient == 22] <- "D"

# make data to plot
g <- unique(tab2$Gene)
sample_order <- tolower(c("2c1c5","2c1c4","2c1c3","2c1c2","2b1c2","2c1c1","2b1c1","2b1c4","2ana3","2a1ae1","2ae3","2b3c5","2b3c4","2b3c3","2b3c2","2b3c1","2b2c5","2b2c4","2b2c3","2b2c1","6a3c1","3a2anac1","3a2aec1","2a4aec1","3a1anac1","6a2c1","2b1c5","4a1aec1","6a1c1","6a4c1","9a5aec1","9a5aec2","9b4c10","9b4c2","9b4c5","9b4c4","9b4c6","9b4c3","9b4c9","9b4c8","9a4anac1","9b4c7","9a4anac2","9a4anac3","9a4aec2","9a4aec3","9b4c1","19anaeroC5","19C1C10","19anaeroC4","19aero1C2","19anaeroC1","19anaeroC2","19C1C5","19aero1C1","19anaeroC3","19B1C10","19C1C4","19C1C9","19B1C8","19B1C4","19B1C6","19B1C9","19C1C1","19C1C6","19C1C7","19C1C8","19aero1C5","19B1C3","19C1C3","19B1C2","19B1C1","19C1C2","19aero1C3","19aero1C4","19B1C5","19B1C7","22a1anac10","22b1c6","22a1anac3","22a1anac5","22b1c3","22b1c5","22b1c1","22b1c4","22b1c110","22b1c2","22b1c8","22a1anac4","22a1anac6","22a1anac2","22a1anac8","22a1anac7","22a1anac9"))
gene <- rep(g, each=length(sample_order))
sample_order <- rep(sample_order, length(g))
patient <- tab2$Patient[match(sample_order, tab2$Sample)]

plot_data <- data.frame(patient=patient, sample=sample_order, gene=gene)

values <- sapply(1:dim(plot_data)[1], function (x) dim(tab2[tab2$Sample == plot_data$sample[x] & tab2$Gene == plot_data$gene[x],])[1])
plot_data["values"] <- values

# add sites
sites <- plot_data[1:94,]
sites$gene <- "Site"
sites$values <- tab2$Site[match(sites$sample, tab2$Sample)]

# modify levels
plot_data <- rbind(plot_data, sites)
plot_data$values <- as.factor(plot_data$values)
levels(plot_data$values) <- c(levels(plot_data$values ), "") # add artificial levels to split the legend into 2 columns
levels(plot_data$values) <- c(levels(plot_data$values ), " ") 
plot_data$values <- factor(plot_data$values, levels(plot_data$values)[c(1,2,6,7,3:5)])

levels_sample <-tolower(c("2c1c5","2c1c4","2c1c3","2c1c2","2b1c2","2c1c1","2b1c1","2b1c4","2ana3","2a1ae1","2ae3","2b3c5","2b3c4","2b3c3","2b3c2","2b3c1","2b2c5","2b2c4","2b2c3","2b2c1","6a3c1","3a2anac1","3a2aec1","2a4aec1","3a1anac1","6a2c1","2b1c5","4a1aec1","6a1c1","6a4c1","9a5aec1","9a5aec2","9b4c10","9b4c2","9b4c5","9b4c4","9b4c6","9b4c3","9b4c9","9b4c8","9a4anac1","9b4c7","9a4anac2","9a4anac3","9a4aec2","9a4aec3","9b4c1","19anaeroC5","19C1C10","19anaeroC4","19aero1C2","19anaeroC1","19anaeroC2","19C1C5","19aero1C1","19anaeroC3","19B1C10","19C1C4","19C1C9","19B1C8","19B1C4","19B1C6","19B1C9","19C1C1","19C1C6","19C1C7","19C1C8","19aero1C5","19B1C3","19C1C3","19B1C2","19B1C1","19C1C2","19aero1C3","19aero1C4","19B1C5","19B1C7","22a1anac10","22b1c6","22a1anac3","22a1anac5","22b1c3","22b1c5","22b1c1","22b1c4","22b1c110","22b1c2","22b1c8","22a1anac4","22a1anac6","22a1anac2","22a1anac8","22a1anac7","22a1anac9"))
plot_data$sample <- factor(plot_data$sample, levels_sample)

levels_gene <-c("Site","espX5","espX4","espL4","entA","entB","entE","entC","fepB","entS","fepD","fepG","fepC","entF","fes","fepA","entD","gad_27","malX","espL1","espR4","yagV/ecpE","yagW/ecpD","yagX/ecpC","yagY/ecpB","yagZ/ecpA","ykgK/ecpR","fdeC","fimH","fimG","fimF","fimD","fimC","fimI","fimA","fimE","fimB","espY1","espX1","chuS","chuA","chuT","chuW","chuX","chuY","chuU","chuV","gad_57","kpsM","kpsD","kpsE","daaF","afaA","nfaE_4","afaC-I","afaD","draP","afaE-I","sat_3","papI","papB","irp2","fyuA","hra","aslA","ompA","gad_17","lpfA_8","gad_42","capU_5","iss_13","gad_13","gad_20","Chromosomal_ompT","iha_11","iss_11","gad_64","papX","usp","hlyF","Episomal_ompT","etsC","iss_12","iroB","iroC","iroD","iroE","iroN_6","mchF_14","cvaC","cba_6","cma_16","iutA","iucD","iucC","iucB","iucA", "traT", "senB_3")

plot_data$gene <- factor(plot_data$gene, levels_gene)

# plot 
colors <- c("firebrick1", "blue", "gold1")

ggplot(plot_data[plot_data$patient != "D",], aes(gene, sample, fill=values)) + 
      geom_tile() + 
      geom_tile(show.legend=FALSE) +
      facet_grid(rows=vars(patient), scales="free", space="free") + 
      theme(panel.background=element_rect(colour="white", fill="white"), panel.border=element_rect(colour="gray", fill=NA), axis.text.x=element_text(angle=45, vjust=1, hjust=0.9, size=5), axis.text.y=element_text(size=5), legend.key = element_blank(),  strip.background = element_rect(fill = "white", colour = "white", linewidth = rel(2))) +
      scale_fill_manual(values = c("gray90", "gray42","white","white", colors), "Gene Site", drop=FALSE) + 
      ggtitle("Virulence genes") +
      guides (fill = guide_legend(ncol=2))
        
      