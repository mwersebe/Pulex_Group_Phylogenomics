#!/usr/bin/Rscript

################################################################################
## Matthew Wersebe
## University of Oklahoma
## Phylogenomics and Evolution of the Daphnia pulex Complex
## Jan 21, 2023
################################################################################

## load requires libraries:

suppressMessages(library(optparse))
suppressMessages(library(Biostrings))
suppressMessages(library(DECIPHER))

## Parse Options list:

option_list = list(
  make_option(c("-i", "--input"), type = "character", default=NA, action = "store",
              help="Input fasta file"),
  make_option(c("-o", "--output"), type = "character", default=NA, action = "store",
              help="Output fasta file"),
  make_option(c("-g", "--outgroup"), type = "character", default=NA, action = "store",
              help="outgroup sequence name"),
  make_option(c("-c", "--genetic_code"), type = "character", default=NA, action = "store",
               help="genetic code to use"),
  make_option(c("-t", "--threads"), type = "numeric", default=8, action = "store",
              help="threads to use"))

opt = parse_args(OptionParser(option_list=option_list))

## Read in unaligned DNA sequences, orient everything correctly:

unaligned_seqs <- Biostrings::readDNAStringSet(opt$input)

unaligned_seqs <- OrientNucleotides(unaligned_seqs, processors = 8)

## Check for Frameshifts and fix:

reference <- unaligned_seqs[opt$outgroup]

uncorrected <- CorrectFrameshifts(unaligned_seqs, translate(reference, genetic.code = getGeneticCode(opt$genetic_code), if.fuzzy.codon = "solve"), processors = opt$threads)


print("Table of Instertions")
print(table(sapply(uncorrected, function(uncorrected) length(uncorrected[["deletions"]]))))

print("Table of Deletions")
print(table(sapply(uncorrected, function(uncorrected) length(uncorrected[["insertions"]]))))

corrected <- CorrectFrameshifts(unaligned_seqs, translate(reference, genetic.code = getGeneticCode(opt$genetic_code), if.fuzzy.codon = "solve"), processors = opt$threads, type = "sequences")


## Perform Codon Aware Alignment:

Aligned_seqs <- AlignTranslation(corrected, geneticCode = getGeneticCode(opt$genetic_code), processors = opt$threads)

## Write the Alignment

writeXStringSet(Aligned_seqs, filepath = opt$output, format = "fasta")

## DONE

print("Done with codon aware alignment")

