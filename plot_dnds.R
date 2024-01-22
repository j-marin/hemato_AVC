library(ggplot2)
library(treemapify)
library(wesanderson)
library(dichromat)
library(recolorize)

rm(list=(ls()))

a <- read.table("../SelectionRes.txt", header=TRUE)
#pal <- wes_palette("Zissou1", 5)

pal <- RColorBrewer::brewer.pal(11,"RdYlGn")
recolorize::plotColorPalette(pal)


p2 <- ggplot(a, aes(area=Nb, fill=kaks, label=paste(site, paste("(n=", Nb, ")", sep=""), sep="\n"), subgroup=patient)) + 
    geom_treemap() + 
    geom_treemap_text(colour="black", place="centre", size=15) + 
    geom_treemap_subgroup_border(colour="white", size=3) + 
    geom_treemap_subgroup_text(place="centre", colour="black") +
    ggtitle("dN/dS") + 
    #scale_fill_gradient2(low = "darkgreen", mid = pal[6], high= pal[1], midpoint=1, na.value = "grey50", space="Lab") +
    scale_fill_gradientn(colors=c(pal[9], pal[6], pal[6], pal[1]), values=c(0,0.08,0.30,1)) +
	geom_blank()
p2 



