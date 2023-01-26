#!/usr/bin/Rscript

################################################################################
## Matthew Wersebe
## University of Oklahoma
## Phylogenomics and Evolution of the Daphnia pulex Complex
## Jan 21, 2023
################################################################################
setwd("/home/giovanni/Phylogeny/MaxLhood")

## load requires libraries:

library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)

## load tree data:

All_SNPs <- ape::read.tree("PulexGroup_Phylogeny.min4.phy.varsites.phy.nwk")

Island_SNPs <- ape::read.tree("Fst_Islands.min4.phy.varsites.phy.nwk")

phangorn::RF.dist(All_SNPs, Island_SNPs, check.labels = T, normalize = T)

## Quickly plot a tangle gram

association <- cbind(All_SNPs$tip.label, Island_SNPs$tip.label)

cophyloplot(All_SNPs, Island_SNPs, assoc = association,
            length.line = 4, space = 500, gap = 2)

## Plot with gg tree:
All_SNPs <- treeio::read.newick("PulexGroup_Phylogeny.min4.phy.varsites.phy.nwk", node.label = "support")

Island_SNPs <- treeio::read.newick("Fst_Islands.min4.phy.varsites.phy.nwk", node.label = "support")

p1 <- ggtree(All_SNPs, aes(color = support)) +
  geom_tiplab(offset = 0.01, size = 2.0)+
  geom_treescale(y = 35, x = 0)+
  xlim(0, 1.0)
p1

p2 <- ggtree(Island_SNPs, aes(color = support)) +
  geom_tiplab(offset = 0.01, size = 2.0)+
  ggtitle("Fst Island SNPs Phylogeny")+
  geom_treescale(y = 30, x = 0)+
  xlim(0, 1.0)
p2

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2) +
  geom_tiplab(data = d2, offset = -0.22, size = 2.0) +
  xlim(0, 2.0)
pp

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

pp + geom_line(aes(x, y, group=label), data=dd, color='grey')+
  ggtitle("All SNPs vs. Fst Island SNPs Tangle Tree")+
  scale_color_continuous(name = "Bootstrap\nSupport")
  