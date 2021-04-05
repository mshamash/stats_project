# BINF 531 - Statistical Bioinformatics
# McGill University
# Final Project - Source Code
# Michael Shamash

library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(plyr)

set.seed(1)

#Import feature table, taxonomy table, and metadata table
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
ps <- merge_phyloseq(OTU, TAX, DATA)

# Filter out contaminating mitochondria and chloroplasts, remove negative controls
ps.filt <- ps %>% subset_taxa(
  Domain == "Bacteria" &
  Family != "mitochondria" &
  Class != "Chloroplast"
) %>% subset_samples(
  Genotype != "negative"
)

# Create phylogenetic tree and attach to phyloseq object
TREE <- rtree(ntaxa(ps.filt), rooted = T, tip.label = taxa_names(ps.filt))
ps.filt <- merge_phyloseq(ps.filt, TREE)

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

# Boxplot of phyla-level relative abundances, need to change code above to Phyla instead of Genus to use this
#phyla <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Other")
#ggplot(dat, aes(fill = factor(Phylum, levels = phyla), x = Sample, y = Abundance)) +
#  geom_bar(stat = "identity", position = "fill") +
#  scale_fill_brewer(palette = "Set3", name = "Phylum") +
#  theme(axis.text.x = element_blank()) +
#  facet_grid(. ~ Age + Genotype, scales = "free", space="free")

# Boxplot of genus-level relative abundances
genera <- c("Acetobacter", "Allobaculum", "Bacteroides", "Campylobacter", "Propionibacterium", "Ruminococcus", "Staphylococcus", "Streptococcus", "Other")
ggplot(dat, aes(fill = factor(taxon, levels = genera), x = Sample, y = Abundance)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set3", name = taxon) +
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
t.test(subset(samp.div, Age == "week12" & Genotype == "WT")$Observed, subset(samp.div, Age == "week12" & Genotype == "KO")$Observed)
t.test(subset(samp.div, Age == "week6" & Genotype == "WT")$Observed, subset(samp.div, Age == "week6" & Genotype == "KO")$Observed)

# Shannon stats
Hutcheson_t_test(subset(samp.div, Age == "week6" & Genotype == "WT")$Shannon, subset(samp.div, Age == "week6" & Genotype == "KO")$Shannon)
Hutcheson_t_test(subset(samp.div, Age == "week12" & Genotype == "WT")$Shannon, subset(samp.div, Age == "week12" & Genotype == "KO")$Shannon)

# Calculate weighted UniFrac distances between samples and plot on PCoA
ps.bdiv <- ordinate(ps.filt, "PCoA", "wunifrac")
plot_ordination(ps.filt, ps.bdiv, type = "samples", color = "Genotype", shape = "Age") +
  geom_point(size = 4) +
  scale_colour_brewer(palette = "Set1") +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15, aes(fill = Genotype)) +
  theme_bw() +
  ggtitle("PCoA on weighted UniFrac distances of all samples")

# Statistical testing (PERMANOVA) of beta diversity between groups
adonis(distance(ps.filt, method = "wunifrac") ~ Age, data = data.frame(sample_data(ps.filt)))
adonis(distance(ps.filt, method = "wunifrac") ~ Genotype, data = data.frame(sample_data(ps.filt)))
adonis(distance(ps.filt, method = "wunifrac") ~ Cage, data = data.frame(sample_data(ps.filt)))


# OTHER CODE GOES HERE


# Custom function for Hutcheson t-test, derived from ecolTest package
Hutcheson_t_test <- function(x, y, shannon.base = exp(1),
                             alternative = "two.sided",
                             difference = 0) {
  dname <-paste ((deparse(substitute(x))), ", ", (deparse(substitute(y))))
  x <- drop(as.matrix(x))
  y <- drop(as.matrix(y))
  if (!is.numeric(x) | !is.numeric(y)) {
    stop("x and y must be numeric")
  }
  if (any(c(x,y) < 0, na.rm = TRUE)) {
    stop("x and y must be non-negative")
  }
  if (any(c(length(x) < 2, length(y) < 2))) {
    stop("x and y must contain at least two elements")
  }
  if (any(c(sum(x, na.rm = TRUE) < 3, sum(y, na.rm = TRUE) < 3))) {
    stop("x and y abundance must be at least two")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop('Package "stats" is needed')
  }
  
  if (any(is.na(c(x, y)))) {
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
    warning("missing values in x and y replaced with zeroes")
  }
  
  alternative <- char.expand(alternative, c("two.sided",
                                            "less", "greater", "auto"))
  if (length(alternative) > 1L || is.na(alternative)) {
    stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
  }
  length_diff <- length(x)-length(y)
  if (length_diff > 0){
    y <- c(y, rep(0, length_diff))
  }
  else if(length_diff < 0){
    x <- c(x, rep(0, abs(length_diff)))
  }
  xy <- matrix(c(x, y), ncol=2)
  N <- colSums(xy)
  H <- (N*log(N, shannon.base)-colSums(xy*log(xy, shannon.base),
                                       na.rm = TRUE))/N
  S<-(colSums(xy*log(xy, shannon.base)**2, na.rm = TRUE) -
        ((colSums(xy*log(xy, shannon.base), na.rm = TRUE)**2)/N))/(N**2)
  Hutchesontstat <- (diff(H[c(2, 1)])-difference)/sqrt(sum(S))
  df <- (sum(S)**2)/(sum(S**2/N))
  estimate_dif <- diff(H[c(2,1)])
  if (alternative == "auto") {
    alternative <- if (estimate_dif < 0) {
      "less"
    }else{
      "greater"
    }
  }
  
  if (alternative == "less") {
    pval <- pt(Hutchesontstat, df)
  }
  else if (alternative == "greater") {
    pval <- pt(Hutchesontstat, df, lower.tail = FALSE)
  }
  else {
    pval <- 2 * pt(-abs(Hutchesontstat), df)
  }
  names(Hutchesontstat) <- "Hutcheson t-statistic"
  names(df) <- "df"
  names(H) <- c("x", "y")
  mu <- difference
  names(mu) <- "difference in H'"
  rval <- list(statistic = Hutchesontstat, parameter = df, p.value = pval,
               estimate = H, null.value = mu,
               method="Hutcheson t-test for two communities",
               alternative = alternative, data.name=dname)
  class(rval) <- "htest"
  return(rval)
}