# Set wd

setwd("C:/Users/mrwil/OneDrive/Documents/NCSU/CCRP/Hackathon_2023/STID12/ProcessedData")

# Compare the distance between the Lotus ecotypes making a kinship matrix

# Packages
require("AGHmatrix")
require("vegan")
require("ggplot2")
require("stringr")
require("tidyverse")
require(phyloseq)

#  Microbial association network construction packages
library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
#install_github("hallucigenia-sparsa/seqtime")
library(seqtime)


# Load the SNP data from a file
snp_data <- read.csv(file = "C:/Users/au672312/Documents/Hackhaton/STID12/ProcessedData/lotus_snps.csv", header = TRUE)

snp_only=snp_data[,-c(1,2)]
G_VrandenLotus=Gmatrix(SNPmatrix = t(snp_only), missingValue=-9, maf=0.05, method="VanRaden")

hm<- heatmap(x = G_VrandenLotus, scale="none")

# Generate 
clust=hclust(d = as.dist(G_VrandenLotus), method = "complete" )
plot(lala)
grp <- cutree(clust, k=4)
table(grp)

#  Merging metadata with groups found attempt but not used in the end

# metadata <- read.csv(file = "C:/Users/au672312/Documents/Hackhaton/STID12/ProcessedData/metadata.csv", header = TRUE, row.names = 1)
# 
# aaa = strsplit(metadata$Genotype, split = 'G')[1:495]
# aaa = sapply(aaa, function(x) x[[2]])
# bbb = stringr::str_pad(string = aaa[4:495], width = 3, side = 'left', pad = '0')
# bbb = paste0('MG', bbb)
# bbb = c(rep('Gifu',3), bbb)


metadata <- read.csv(file = "metadata_grpG.csv", header = TRUE)
rownames(metadata)=metadata$Sample

## Beta diversity plots

amplicon_table_nofilter <- read.table(file = "dada2_table.txt", header = TRUE, row.names = 1, sep="\t")

amplicon_table_filtered <- read.table(file = "data2_table_isolates_only.txt", header = TRUE, row.names = 1, sep="\t")




transpose=t(amplicon_table_nofilter)

Bray_curtis_df <- as.matrix(vegdist(transpose, binary=TRUE, method = "bray", diag = T))



#Make PCoA plot for Bray Curtis Distance matrix
pcoa = cmdscale(Bray_curtis_df, k=3, eig=T)
points = as.data.frame(pcoa$points)
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig

pcoa_df = merge(points,metadata, by.x = 0, by.y = "Sample")
pcoa_df_nona = na.omit(pcoa_df)



#  Plot Beta div replicate per replicate : 
p = ggplot(pcoa_df, aes(x=x, y=y, color=Replicate, stroke=1))+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  geom_point(alpha=1, size=2)+stat_ellipse(aes(group = Replicate))+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  theme(axis.title = element_text(size = 14, hjust = 1), axis.text = element_text(size = 12, hjust = 1))+
  guides(fill = guide_legend(override.aes=list(shape=21, size=5)))
p

#  Plot Beta div and color by genetic group obtained with kinship matrix and hcl

pcoa_df_nona$grpVR=as.factor(pcoa_df_nona$grpVR)

q = ggplot(pcoa_df_nona, aes(x=x, y=y, color=grpVR, stroke=1))+
  theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line.x=element_line(size=.5, colour="black"),axis.line.y=element_line(size=.5, colour="black"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black", size=7),legend.position="right",legend.background=element_blank(),legend.key=element_blank(),legend.text= element_text(size=10),text=element_text(family="sans", size=10))+
  geom_point(alpha=1, size=2)+stat_ellipse(aes(group = grpVR))+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  theme(axis.title = element_text(size = 14, hjust = 1), axis.text = element_text(size = 12, hjust = 1))+
  guides(fill = guide_legend(override.aes=list(shape=21, size=5)))+facet_wrap(facets = "Replicate")

q

# Using phyloseq prior to construct the association network
# tax table
tax_df <- read.table("taxonomy_FL.txt", header=TRUE, sep=";")
rownames(tax_df) = tax_df$Isolate

tax_df_filt <- subset(tax_df, rownames(tax_df) %in% rownames(amplicon_table_filtered))



OTU=phyloseq::otu_table(object = amplicon_table_filtered, taxa_are_rows = T)
TAX = tax_table(as.matrix(tax_df_filt))
samples = sample_data(metadata)
phylo = phyloseq(OTU,TAX, samples)

#  Microbial association network construction
library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)
# install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)


filterobj=filterTaxonMatrix(OTU,minocc=20,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f=filterobj$mat
taxa.f=TAX[setdiff(1:nrow(TAX),filterobj$filtered.indices),]



spiec.out=spiec.easi(phylo, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit$stars, vertex.attr=list(name=taxa_names(phylo)))
plot_network(spiec.graph, phylo, type='taxa', color="Family", label="Isolate")
degree(spiec.graph)

#Rep1
phylo1 = subset_samples(phylo, Replicate=="R1")
filterobj=filterTaxonMatrix(otu_table(phylo1),minocc=20,keepSum = TRUE, return.filtered.indices = TRUE)
otus.f=filterobj$mat
taxa = tax_table(phylo1)
taxa.f=taxa[setdiff(1:nrow(taxa),filterobj$filtered.indices),]

spiec.out=spiec.easi(phylo1, method="mb",icov.select.params=list(rep.num=20))
spiec.out1 = spiec.out
spiec.graph1=adj2igraph(spiec.out1$refit$stars, vertex.attr=list(name=taxa_names(phylo1)))
plot_network(spiec.graph1, phylo1, type='taxa', color="Family", label="Isolate")
degree(spiec.graph)

phylo2 = subset_samples(phylo, Replicate=="R2")
spiec.out2=spiec.easi(phylo2, method="mb",icov.select.params=list(rep.num=20))
spiec.graph2=adj2igraph(spiec.out2$refit$stars, vertex.attr=list(name=taxa_names(phylo1)))
plot_network(spiec.graph2, phylo2, type='taxa', color="Family", label="Isolate")
degree(spiec.graph)

ps = subset_samples(phylo, Genotype!="Stick")
ps = subset_samples(ps, Genotype!="Input")

spiec.out=spiec.easi(ps, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit$stars, vertex.attr=list(name=taxa_names(ps)))
plot_network(spiec.graph, ps, type='taxa', color="Family", label="Isolate")
degree(spiec.graph)

ps1 = subset_samples(ps, Replicate=="R1")
ps2 = subset_samples(ps, Replicate=="R2")
ps3 = subset_samples(ps, Replicate=="R3")

spiec.out1=spiec.easi(ps1, method="mb",icov.select.params=list(rep.num=20))
spiec.graph1=adj2igraph(spiec.out1$refit$stars, vertex.attr=list(name=taxa_names(ps1)))
plot_network(spiec.graph1, ps1, type='taxa', color="Family", label="Isolate")
degree(spiec.graph1)

spiec.out2=spiec.easi(ps2, method="mb",icov.select.params=list(rep.num=20))
spiec.graph2=adj2igraph(spiec.out2$refit$stars, vertex.attr=list(name=taxa_names(ps2)))
plot_network(spiec.graph2, ps2, type='taxa', color="Family", label="Isolate")
degree(spiec.graph2)

spiec.out3=spiec.easi(ps3, method="mb",icov.select.params=list(rep.num=20))
spiec.graph3=adj2igraph(spiec.out3$refit$stars, vertex.attr=list(name=taxa_names(ps3)))
plot_network(spiec.graph3, ps3, type='taxa', color="Family", label="Isolate")
degree(spiec.graph3)

rbind(degree(spiec.graph1), degree(spiec.graph2), degree(spiec.graph3))

cor(degree(spiec.graph1), degree(spiec.graph2), method='spearman')
cor(degree(spiec.graph1), degree(spiec.graph3), method='spearman')
cor(degree(spiec.graph2), degree(spiec.graph3), method='spearman')

library(ggplot2)
p1 = qplot(x=degree(spiec.graph1), y = degree(spiec.graph2)) + geom_smooth(method = 'lm') +
  xlab("Rep 1 degree centrality") + ylab("Rep 2 degree centrality") + 
  labs(title="Spearman correlation coefficient = 0.402")
p2 = qplot(x=degree(spiec.graph1), y = degree(spiec.graph3)) + geom_smooth(method = 'lm') +
  xlab("Rep 1 degree centrality") + ylab("Rep 3 degree centrality") +
  labs(title="Spearman correlation coefficient = 0.379")
p3 = qplot(x=degree(spiec.graph2), y = degree(spiec.graph3)) + geom_smooth(method = 'lm') +
  xlab("Rep 1 degree centrality") + ylab("Rep 2 degree centrality") +
  labs(title="Spearman correlation coefficient = 0.643")

ggsave("cor12.png", p1)
ggsave("cor13.png", p2)
ggsave("cor23.png", p3)



identical(names(degree(spiec.graph1)), names(degree(spiec.graph2)))
identical(names(degree(spiec.graph1)), names(degree(spiec.graph3)))
identical(names(degree(spiec.graph2)), names(degree(spiec.graph3)))

int123 = intersection(spiec.graph1, spiec.graph2, spiec.graph3) %>% as_long_data_frame
int12 = intersection(spiec.graph1, spiec.graph2)
int13 = intersection(spiec.graph1, spiec.graph3)
int23 = intersection(spiec.graph2, spiec.graph3)

write.csv(int123, "3-way_intersection.csv")

sp1 = spiec.graph1 %>% as_long_data_frame %>% as.data.frame
colnames(sp1)[4:5] <- c("v1","v2")
sp2 = spiec.graph2 %>% as_long_data_frame
sp3 = spiec.graph3 %>% as_long_data_frame
colnames(sp2)[4:5] <- c("v1","v2")
colnames(sp3)[4:5] <- c("v1","v2")


edges1 = paste(sp1$v1, sp1$v2, sep="_")
edges2 = paste(sp2$v1, sp2$v2, sep="_")
edges3 = paste(sp3$v1, sp3$v2, sep="_")

intersection(edges1, edges2)

require(VennDiagram)
venn.diagram(list(edges1, edges2, edges3), filenam="venn.png")

library(RVenn)
ggvenn(Venn(list(edges1, edges2, edges3)))



difference(spiec.graph1, spiec.graph2)
difference(spiec.graph2, spiec.graph1)


set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
