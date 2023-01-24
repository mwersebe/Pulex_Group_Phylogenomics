#!/usr/bin/Rscript

################################################################################
## Matthew Wersebe
## University of Oklahoma
## Phylogenomics and Evolution of the Daphnia pulex Complex
## Jan 21, 2023
################################################################################
setwd("/home/giovanni/Phylogeny/MitoPhylo")
## load requires libraries:

library(ape)
library(treeio)
library(ggtree)
library(ggplot2)

## Read in Iqtree Phylogeny:

tree <- treeio::read.newick("Mito_Genes.nwk", node.label = "support")

## Get Clades for the

plot(tree@phylo)

getMRCA(tree@phylo, tip = c("SRR2148391", "SRR2148403"))

getMRCA(tree@phylo, tip = c("SRR2148408", "B119_12-1_S32"))

getMRCA(tree@phylo, tip = c("MP_4-4_S116", "C31_12-1_S40"))

getMRCA(tree@phylo, tip = c("SRR10160692", "Sedlec_02_S60"))

getMRCA(tree@phylo, tip = c("SRR10160660", "Daphnia_mitsukuri"))

getMRCA(tree@phylo, tip = c("SRR10160568", "Daphnia_obtusa"))

tree@phylo <- root.phylo(tree@phylo, node = 86)

## Plot Tree

ggtree(tree,aes(color = support)) + geom_tiplab(offset = 0.032, size = 2.5)+
  scale_color_continuous(low='dark blue', high='orange', name = "Bootstrap\nSupport") +
  theme(legend.position="right")+
  ggtitle("Mitochondrial Genome Phylogeny")+
  geom_cladelab(node = 92, label = "Pan-Arctic\nPulex Group", offset = 0.21)+
  geom_cladelab(node = 142, label = "Pulicaria\nGroup", offset = 0.2)+
  geom_cladelab(node = 155, label = "Tenebrosa\nGroup", offset = 0.27)+
  geom_cladelab(node = 162, label = "European Pulex", offset = 0.2)+
  geom_cladelab(node = 166, label = "Mitsukuri Group", offset = 0.22)+
  geom_treescale()+
  xlim(0, 1.2)
  



