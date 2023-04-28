library(cluster)
library(compositions)
library(permute)
library(phyloseq)
library(MiHC)
library(tidyverse)
library(qvalue)
library(ggrepel)

#https://github.com/hk1785/MiHC
#Tried MiHC but only for binary y data (so not continuous phenotype)
setwd("D:/Google_Drive/Post-doc/Project_MATRIX/Hackathon_CCRP/ProcessedDataSTID12/ProcessedData/")

#data(phy)
#otu.tab <- otu_table(phy)
#tree <- phy_tree(phy)
#y <- sample_data(phy)$y
#covs <- data.frame(matrix(NA, length(y), 2))
#covs[,1] <- as.numeric(sample_data(phy)$x1)
#covs[,2] <- as.factor(sample_data(phy)$x2)

#set.seed(123)
#out <- MiHC(y, covs=covs, otu.tab=otu.tab, tree=tree, model="binomial")
#out
#MiHC.plot(out)

#https://www.rdocumentation.org/packages/MiRKAT/versions/1.2.3/topics/MMiRKAT
library(MiRKAT)
library(vegan)

#Did filtering using sum, but should be improved upon
asv = read.csv("isolatesXPlantGeno.csv") %>% pivot_longer(-"X") %>%
  separate(name, into = c("genotype", "replicate")) %>% group_by(X, genotype) %>%
  mutate(sum = sum(value)) %>% filter(sum > 100)

#replaced zeros with very small number because distance matrix doesn't take in zeros, however should be double checked
asv_matrix = asv %>% select(X, genotype, sum) %>% unique() %>% 
  pivot_wider(names_from = X, values_from = sum, values_fill = 0) %>% 
  filter(genotype != "SC", genotype != "Stick") %>%
  mutate(across(everything(), ~replace(., . == 0, 1e-6)))

#removed few rows with lots of NA and replaced few with zero, should be changed -- possiblity to take average
pheno = read.csv("PhenotypicMeasurements.csv") %>% select(-Batch) %>%
  t() %>% data.frame() %>% select(-6, -7, -13, -14, -20, -21) %>% 
  mutate_all(~if_else(is.na(.), 0, .)) %>% as.matrix()

#tried doing without for loop using map() but got error
#p = asv_matrix %>% pivot_longer(-genotype) %>% group_by(name) %>% nest() %>% ungroup() %>%
#  mutate(D = map(data, ~ .x %>% column_to_rownames("genotype") %>% vegdist(method="bray") %>% as.matrix())) %>%
#  mutate(K.BC = map(D, ~ .x %>% D2K()))

names = asv %>% ungroup() %>% select(X) %>% unique() %>%
   unlist()

p_all = list()

#MMiRKAT also takes in a list of kernels, so this could be improved upon (for loop not needed anymore)
for (i in names) {
  test = asv_matrix %>% pivot_longer(-genotype) %>%
    filter(name == i) %>% select(-name) %>% column_to_rownames("genotype")
  test_D= as.matrix(vegdist(test, method="bray"))
  test_K.BC = D2K(test_D)
  test_p = MMiRKAT(pheno, X=NULL, test_K.BC)
  p_all = c(p_all, test_p)
}

p_all = p_all %>% as.data.frame() %>% t() %>% as.data.frame()

plot_data = names %>% as.data.frame() %>% bind_cols(p_all) %>%
  rename(name = 1, p_value = 2) %>% mutate(qvalue = qvalue(.$p_value, lfdr.out = TRUE)[["qvalues"]]) %>%
  mutate(log = -log10(qvalue))
label_data = plot_data %>% filter(qvalue < 0.05) %>% filter(str_detect(name, "^L"))

p = ggplot(plot_data, aes(name, log, color = log)) +
  geom_point(size = 1) +
  scale_color_continuous() +
  labs(x = "isloates", y = "-log10(q value)") +
  geom_hline(yintercept = c(-log10(0.05)), color = "firebrick") +
  theme_light() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank()) +
  geom_text_repel(data = label_data, aes(name, log, label = name))
p

ggsave("manhattan_plot.png", p, device = "png", dpi = 600)
