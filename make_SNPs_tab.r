# make result table
tab <- rbind(read.table("patient2_complete.tab", header=TRUE), 
read.table("patient9_cl2_complete.tab", header=TRUE),
read.table("patient16_complete.tab", header=TRUE))
recap <- read.table("../../recap.txt", header=TRUE)
head(recap)

recap$Patient[c(which(recap$Patient == 3), which(recap$Patient == 4), which(recap$Patient == 6))] <- 2 # patient 2 = 2, 3, 4, 6

tab <- tab[order(tab[,2]),]
tab <- tab[order(tab[,4]),] #REF
tab <- tab[order(tab[,5]),] #ALT
tab <- tab[order(tab[,1]),] #CHROM
tab <- tab[order(tab[,2]),] # POS
head(tab)

# check sample names
strains <- read.table("../../recapSeqFinal.txt", head=TRUE)$Sample 
tab <- tab[tolower(tab$NOM) %in% strains,]
dim(tab)

# make SNP ID
tabx <- tab[,c(1,2,4,5)]
head(tabx)
d <- duplicated(tabx)
rg <- c(which(d == "FALSE"), length(d)+1)
id <- paste("snp_", 1:length(rg), sep="")

for (i in 1:length(rg)) {
	d[rg[i]:(rg[i+1]-1)] <- id[i]
	}

tab2 <- data.frame(tab["PATIENT"])
tab2["Sample"] <- tab$NOM
tab2["Site"] <- recap$Prélèvement[match(tolower(tab$NOM), recap$Nom)]
tab2["SNP_ID"] <- d
tab2["ref"] <- tab$REF
tab2["alt"] <- tab$ALT
tab2["count"] <- tab$EVIDENCE
tab2["type"] <- tab$TYPE
tab2["contig"] <- tab$CHROM
tab2["Pos"] <- tab$POS
tab2["Gene"] <- tab$GENE
tab2["Type"] <- sapply(strsplit(tab$EFFECT,"/"), `[`, 1)
tab2["AA"] <- sapply(strsplit(tab$EFFECT,"/"), `[`, 3)
tab2["Product"] <- tab$PRODUCT
tab2["Alt_count"] <- as.numeric(sapply(strsplit(sapply(strsplit(tab2$count, ":"), `[`, 2), "_"),`[`, 1))
tab2["Ref_count"] <- as.numeric(sapply(strsplit(tab2$count, ":"), `[`, 3))
head(tab2)

write.table(tab2, "recap_snp_patients_tab_final.txt", row.name=FALSE, sep="\t")
