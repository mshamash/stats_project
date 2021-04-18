# BINF 531 - Statistical Bioinformatics
# McGill University
# Final Project - Analysis Code
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

# Clean up OTU names
taxa_names(ps.filt) <- paste("OTU", 1:664, sep="_")

# Clean up taxonomy and replace blanks with previous known taxon
tax.clean <- data.frame(tax_table(ps.filt))
for (i in 1:7) tax.clean[, i] <- as.character(tax.clean[, i])
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i, 2] == "") {
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i, 3] == "") {
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i, 4] == "") {
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i, 5] == "") {
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i, 6] == "") {
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i, 7] == "") {
    tax.clean$Species[i] <- paste("Genus", tax.clean$Genus[i], sep = "_")
  }
}

tax_table(ps.filt) <- as.matrix(tax.clean)

# Filter 6week and 12week samples
ps.filt.6wk <- ps.filt %>% subset_samples(Age == "week6")
ps.filt.12wk <- ps.filt %>% subset_samples(Age == "week12")

# Differential abundance analysis
library(DESeq2)
library(ggpubr)

# Calculate geometric means
gm_mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
alpha = 0.001

dds.age <- phyloseq_to_deseq2(ps.filt, ~ Cage + Genotype + Age)
dds.age$Age <- relevel(dds.age$Age, ref = "week6")
geoMeans = apply(counts(dds.age), 1, gm_mean)
dds.age = estimateSizeFactors(dds.age, geoMeans = geoMeans)
dds.age = DESeq(dds.age, fitType="local")
res.age = results(dds.age)
res.age = res.wt.age[order(res.age$padj, na.last=NA), ]
sigtab.age = res.wt.age[(res.age$padj < alpha), ]
sigtab.age = cbind(as(sigtab.age, "data.frame"), as(tax_table(ps.filt)[rownames(sigtab.age), ], "matrix"))
sigtab.age

x = tapply(sigtab.age$log2FoldChange, sigtab.age$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.age$Genus = factor(as.character(sigtab.age$Genus), levels=names(x))

sigtab.age$OTU <- rownames(sigtab.age)
ggbarplot(sigtab.age, x = "OTU", y = "log2FoldChange",
          fill = "Genus",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "Set3",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          ylab = "log2FoldChange",
          legend.title = "Genus",
          title = "Differentially abundant OTUs in 12-week-old mice compared to 6-week-old mice",
          rotate = TRUE,
          ggtheme = theme_classic()
)

dds.geno <- phyloseq_to_deseq2(ps.filt, ~ Cage + Age + Genotype)
dds.geno$Genotype <- relevel(dds.geno$Genotype, ref = "WT")
geoMeans = apply(counts(dds.geno), 1, gm_mean)
dds.geno = estimateSizeFactors(dds.geno, geoMeans = geoMeans)
dds.geno = DESeq(dds.geno, fitType="local")
res = results(dds.geno)
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.filt)[rownames(sigtab), ], "matrix"))
head(sigtab)

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

sigtab$OTU <- rownames(sigtab)
ggbarplot(sigtab, x = "OTU", y = "log2FoldChange",
          fill = "Genus",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "Set3",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          ylab = "log2FoldChange",
          legend.title = "Genus",
          title = "Differentially abundant OTUs in CD109-KO mice compared to WT mice",
          rotate = TRUE,
          ggtheme = theme_classic()
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

# Calculate Bray Curtis distances between samples and plot on PCoA
ps.bdiv <- ordinate(ps.filt, "PCoA", "bray")
ps.bdiv.6wk <- ordinate(ps.filt.6wk, "PCoA", "bray")
ps.bdiv.12wk <- ordinate(ps.filt.12wk, "PCoA", "bray")

plot_ordination(ps.filt, ps.bdiv, type = "samples", color = "Genotype", shape = "Age") +
  geom_point(size = 4) +
  scale_colour_brewer(palette = "Set1") +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15, aes(fill = Genotype)) +
  theme_bw() +
  ggtitle("PCoA on Bray Curtis distances of all samples")

plot_ordination(ps.filt.6wk, ps.bdiv, type = "samples", color = "Genotype") +
  geom_point(size = 4) +
  scale_colour_brewer(palette = "Set1") +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15, aes(fill = Genotype)) +
  theme_bw() +
  ggtitle("PCoA on Bray Curtis distances of 6 week samples")

plot_ordination(ps.filt.12wk, ps.bdiv, type = "samples", color = "Genotype") +
  geom_point(size = 4) +
  scale_colour_brewer(palette = "Set1") +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15, aes(fill = Genotype)) +
  theme_bw() +
  ggtitle("PCoA on Bray Curtis distances of 12 week samples")

# Statistical testing (PERMANOVA) of beta diversity between groups
set.seed(123456)
adonis(distance(ps.filt, method = "bray") ~ Genotype, data = data.frame(sample_data(ps.filt)))
adonis(distance(ps.filt, method = "bray") ~ Age, data = data.frame(sample_data(ps.filt)))
adonis(distance(ps.filt, method = "bray") ~ Cage, data = data.frame(sample_data(ps.filt)))

adonis(distance(ps.filt.6wk, method = "bray") ~ Genotype, data = data.frame(sample_data(ps.filt.6wk)))
adonis(distance(ps.filt.6wk, method = "bray") ~ Cage, data = data.frame(sample_data(ps.filt.6wk)))

adonis(distance(ps.filt.12wk, method = "bray") ~ Genotype, data = data.frame(sample_data(ps.filt.12wk)))
adonis(distance(ps.filt.12wk, method = "bray") ~ Cage, data = data.frame(sample_data(ps.filt.12wk)))

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
library(knitr)
library(pROC)
library(kableExtra)

# Random forest classifier for genotype
response.geno <- as.factor(sample_data(ps.filt)$Genotype)
predictors <- t(otu_table(ps.filt))

rf.data.geno <- data.frame(response.geno, predictors)
set.seed(2)
ps.classify.geno <- randomForest(response.geno ~ ., data = rf.data.geno, ntree = 100)
print(ps.classify.geno)
plot(ps.classify.geno)

roc.geno <- roc(rf.data.geno$response.geno, ps.classify.geno$votes[, 2], levels = c("WT", "KO"))
plot(roc.geno)
auc(roc.geno)

imp.geno <- importance(ps.classify.geno)
imp.geno <- data.frame(predictors = rownames(imp.geno), imp.geno)

imp.sort.geno <- arrange(imp.geno, desc(MeanDecreaseGini))
imp.sort.geno$predictors <- factor(imp.sort.geno$predictors, levels = imp.sort.geno$predictors)
imp.20.geno <- imp.sort.geno[1:20, ]

ggplot(imp.20.geno, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n based on Genotype") +
  xlab("OTU ID") +  
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(x=10, y=0.4, label="56.67%")

otunames.geno <- imp.20.geno$predictors
r.geno <- rownames(tax_table(ps.filt)) %in% otunames.geno
imp.otu.g <- tax_table(ps.filt)[r.geno, ] %>% 
  kbl(caption = "") %>% 
  kable_styling()
save_kable(
  imp.otu.g,
  "kable.geno.png",
  bs_theme = "simplex",
  self_contained = TRUE,
  extra_dependencies = NULL,
  latex_header_includes = NULL,
  keep_tex = FALSE,
  density = 300
)
dev.off()

# Random forest classifier for age
response.age <- as.factor(sample_data(ps.filt)$Age)
predictors <- t(otu_table(ps.filt))

rf.data.age <- data.frame(response.age, predictors)
set.seed(2)
ps.classify.age <- randomForest(response.age ~ ., data = rf.data.age, ntree = 100)
print(ps.classify.age)

roc.age <- roc(rf.data.age$response.age, ps.classify.age$votes[, 2])
plot(roc.age)
auc(roc.age)

imp.age <- importance(ps.classify.age)
imp.age <- data.frame(predictors = rownames(imp.age), imp.age)

imp.sort.age <- arrange(imp.age, desc(MeanDecreaseGini))
imp.sort.age$predictors <- factor(imp.sort.age$predictors, levels = imp.sort.age$predictors)
imp.20.age <- imp.sort.age[1:20, ]

ggplot(imp.20.age, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n based on Age") +
  xlab("OTU ID") +  
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text(x=10, y=0.4, label="23.33%")

otunames.age <- imp.20.age$predictors
r.age <- rownames(tax_table(ps.filt)) %in% otunames.age
imp.otu.age <- tax_table(ps.filt)[r.age, ] %>% 
  kbl(caption = "") %>% 
  kable_styling()
save_kable(
  imp.otu.age,
  "kable.age.png",
  bs_theme = "simplex",
  self_contained = TRUE,
  extra_dependencies = NULL,
  latex_header_includes = NULL,
  keep_tex = FALSE,
  density = 300
)
dev.off()

# Generate ROC for RF classifiers of genotype and age, calculate AUC-ROC
png("roc_curve.png", width = 1000,
    height = 1000,
    units = "px",
    pointsize = 25, bg = "white", res = NA)
par(mar = c(5.1, 5.1, 4.1, 3.1))
plot(roc.geno, col = "red", 
     xlab = "Specificity (False positive rate)",
     ylab = "Sensitivity (True positive rate)", 
     main = "ROC Curves of Random Forest Models Classifying Samples\n Based on Genotype and Age")
lines(roc.age, col = "blue")
legend("bottomright", legend = c("Genotype", "Age"),
       col = c("red", "blue"), pch = 1)
text(0.3, 0.77, paste("AUC:", formatC(auc(roc.geno))), col = "red", cex = 1)
text(0.3, 0.87, paste("AUC:", formatC(auc(roc.age))), col = "blue", cex = 1)
dev.off()
