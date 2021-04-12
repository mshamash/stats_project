# BINF 531 - Statistical Bioinformatics
# McGill University
# Final Project - Source Code
# Michael Shamash, Wen Da Lu, Garrie Peng

library(plyr)
library(tidyverse)
library(phyloseq)
library(vegan)

# Diversity Analysis Section - Michael ----------------------------------------
source("HelperFunctions.R")

# Import feature table, taxonomy table, and metadata table
features <- read.csv("features.tsv", header = T, row.names = 1, as.is = T, sep = "\t", skip = 1)
taxonomy <- read.csv("taxonomy.tsv", header = T, row.names = 1, as.is = T, sep = "\t")
metadata <- read.csv("metadata.csv", header = T, row.names = 1, as.is = T)

# Process taxonomy data into correct format
taxaLevels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy.mat <- matrix(nrow = nrow(taxonomy), ncol = length(taxaLevels))
rownames(taxonomy.mat) <- rownames(taxonomy)
colnames(taxonomy.mat) <- taxaLevels
for (i in 1:nrow(taxonomy)) for (j in 1:7) taxonomy.mat[i, j] <- substring(unlist(strsplit(taxonomy[i, 1], "; ")), 4)[j]

# Create phyloseq object with above data
OTU <- otu_table(features, taxa_are_rows = T)
TAX <- tax_table(taxonomy.mat)
DATA <- sample_data(metadata)
ps <- merge_phyloseq(OTU, TAX, DATA)#, TREE)

# Filter out contaminating mitochondria and chloroplasts, remove negative controls
ps.filt <- ps %>% subset_taxa(
  Domain == "Bacteria" &
  Family != "mitochondria" &
  Class != "Chloroplast"
) %>% subset_samples(
  Genotype != "negative"
)

# Plot basic rarefaction curve
rarecurve(t(otu_table(ps.filt)), step = 50, ylab = "Richness (# of OTUs)", xlab = "# of reads", label = F)

# Clean up relative abundance data to merge low-abundance taxa into "Other"
phy <- transform_sample_counts(ps.filt, function(x) x/sum(x)) # Get relative abundances
glom <- tax_glom(phy, taxrank = "Genus") # Agglomerate taxa
dat <- psmelt(glom) # Create dataframe from phyloseq object
dat$Genus <- as.character(dat$Genus) # Convert genus to a character vector from a factor
medians <- ddply(dat, ~ Genus, function(x) c(median = median(x$Abundance))) # Group dataframe by genus, calculate median rel. abund.
other <- medians[medians$median <= 0.01, ]$Genus # Find genera whose rel. abund. is less than 1%
dat[dat$Genus %in% other, ]$Genus <- "Other" # Change their name to "Other"

# Convert Age and Genotype to factors
dat$Age <- factor(dat$Age, levels = c("week6", "week12"))
dat$Genotype <- factor(dat$Genotype, levels = c("WT", "KO"))

# Boxplot of genus-level relative abundances
genera <- c("Acetobacter", "Allobaculum", "Bacteroides", "Campylobacter", "Propionibacterium", "Ruminococcus", "Staphylococcus", "Streptococcus", "Other")
#phyla <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Other")
ggplot(dat, aes(fill = factor(Genus, levels = genera), x = Sample, y = Abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set3", name = "Genus") +
  theme_bw()+
  theme(axis.text.x = element_blank()) +
  facet_grid(. ~ Age + Genotype, scales = "free", space="free")

# Calculate alpha diversity (Shannon diversity and observed richness) for each sample
samp.div <- estimate_richness(ps.filt, measures = c("Shannon", "Observed"))
samp.div$Genotype <- rep("", 30)
samp.div$Age <-  rep("", 30)
for(i in 1:nrow(samp.div)){
  tmp <- sample_data(ps.filt)[match(rownames(samp.div)[i], rownames(sample_data(ps.filt))), ]
  samp.div[i, 3] <- tmp$Genotype
   samp.div[i, 4] <- tmp$Age
}
samp.div$Genotype <- factor(samp.div$Genotype, levels = c("WT", "KO"))
samp.div$Age <- factor(samp.div$Age, levels = c("week6", "week12"))

# Plot alpha diversity boxplots by age and genotype
col <- c(WT="#377EB8", KO="#E41A1C")
ggplot(samp.div, aes(fill = Genotype, x = Genotype, y = Observed)) + geom_boxplot() + facet_grid(. ~ Age) + theme_bw() + scale_fill_manual(values = col)
ggplot(samp.div, aes(fill = Genotype, x = Genotype, y = Shannon)) + geom_boxplot() + facet_grid(. ~ Age) + theme_bw() + scale_fill_manual(values = col)

# Statistical testing of alpha diversity metrics between groups
# Richness stats
t.test(subset(samp.div, Age == "week6" & Genotype == "WT")$Observed, subset(samp.div, Age == "week6" & Genotype == "KO")$Observed)
t.test(subset(samp.div, Age == "week12" & Genotype == "WT")$Observed, subset(samp.div, Age == "week12" & Genotype == "KO")$Observed)

# Shannon stats
Hutcheson_t_test(subset(samp.div, Age == "week6" & Genotype == "WT")$Shannon, subset(samp.div, Age == "week6" & Genotype == "KO")$Shannon)
Hutcheson_t_test(subset(samp.div, Age == "week12" & Genotype == "WT")$Shannon, subset(samp.div, Age == "week12" & Genotype == "KO")$Shannon)

# Calculate weighted UniFrac distances between samples and plot on PCoA
ps.bdiv <- ordinate(ps.filt, "PCoA", "bray")
plot_ordination(ps.filt, ps.bdiv, type = "samples", color = "Genotype", shape = "Age") +
  geom_point(size = 4) +
  scale_colour_brewer(palette = "Set1") +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15, aes(fill = Genotype)) +
  theme_bw() +
  ggtitle("PCoA on Bray Curtis distances of all samples")

# Statistical testing (PERMANOVA) of beta diversity between groups
set.seed(123456)
adonis(distance(ps.filt, method = "bray") ~ Genotype, data = data.frame(sample_data(ps.filt)))
adonis(distance(ps.filt, method = "bray") ~ Age, data = data.frame(sample_data(ps.filt)))
adonis(distance(ps.filt, method = "bray") ~ Cage, data = data.frame(sample_data(ps.filt)))


# Correlation & Permutation Test Section - Wen Da -----------------------------
library(reshape2)

# Generate Pearson correlation heatmap
corr.plot <- function(features) {
  cormat <- abs(round(cor(features), 2))
  
  get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  upper_tri <- get_lower_tri(cormat)
  
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  return(ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
           geom_tile(color = "white") +
           scale_fill_gradient2(low = "blue", high = "red", mid = "green", 
                                midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                                name ="Pearson\nCorrelation") +
           theme_minimal() +
           theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                            size = 10, hjust = 1)) +
           coord_fixed())
}
corr.plot(features)

# Perform permutation tests
data1 <- features[, 1:5]
data2 <- features[, 6:15]
data3 <- features[, 16:20]
data4 <- features[, 21:30]
data5 <- features[, 1:15]
data6 <- features[, 16:30]

perm <- function(x,y) {
  t.res <- rep(0, length(x[, 1]))
  
  P <- length(x[, 1])
  
  for (i in 1:P) {
    test.stat <- abs(sum(x[i, ])/length(x) - sum(y[i, ])/length(y)) # Get mean difference
    
    permut <- 1000 # Permute 1000 times
    vari <- cbind(x[i, ], y[i, ]) # All total variables
    
    PermSamples <- matrix(0, nrow = length(x) + length(y), ncol = permut) #sample storing
    Perm.test.stat <- rep(0, permut) #mean storing
    
    for(j in 1:permut) {
      temp <- sample(vari, 
                     size = length(x) + length(y), 
                     replace = FALSE)
      
      names(temp) <- NULL
      
      PermSamples[, j] <- t(temp) 
      
      Perm.test.stat[j] <- abs(sum(PermSamples[1:length(x), j])/length(x) 
                               - sum(PermSamples[(length(x)+1):(length(x) + length(y)), j])/length(y))
    }
    
    p.val <- mean(Perm.test.stat >= test.stat)
    
    t.res[i] <- p.val
  }
  return(t.res)
}

comp.data1.data2 <- perm(data1, data2)
comp.data1.data3 <- perm(data1, data3)
comp.data2.data4 <- perm(data2, data4)
comp.data3.data4 <- perm(data3, data4)
comp.data5.data6 <- perm(data5, data6)

a <- comp.data1.data2 <= 0.05
b <- comp.data1.data3 <= 0.05
c <- comp.data2.data4 <= 0.05
d <- comp.data3.data4 <= 0.05

sig <- rowSums(cbind(a, b, c, d))

e <- as.data.frame(comp.data5.data6 <= 0.05)
colnames(e) <- "sig"

temp <- cbind(features, e)
temp <- subset(temp, sig == TRUE) 
temp <- subset(temp, select = -sig)

barplot(sig, ylab = "Number of Significant Match", xlab = "OTU Number",
        main = "Number of Significant Match between the Four Groups Depending on the OTU")

corr.plot(temp) # Correlation heatmap of the new matrix

rownames(temp) # Obtain the OTU names


# Random Forests Section - Garrie ---------------------------------------------
library(randomForest)
