
library("phyloseq")
library("ggplot2")
library("vegan")
library("agricolae")
library("cowplot")
library("pairwiseAdonis")
library("tidyr")
library("emmeans")

## reading in ASV/OTU table RDS file generated from dada2 
# otu_file <- read.csv(file="5E_OTU_table.csv", row.names=1)
# saveRDS(otu_file, "5E_OTU_table.RDS")
otu_file <- readRDS(file="./data/5E_OTU_table.RDS") #READ IN RDS B/C FASTER
## convert to phyloseq object
OTU<-otu_table(otu_file, taxa_are_rows = FALSE)

## reading in Taxonomy RDS file generated from dada2 
# tax_file <-read.csv(file="5E_TAX_table.csv", row.names=1)
# tax_file <- as.matrix(tax_file)
# saveRDS(tax_file, "5E_TAX_table.RDS")
tax_file <- readRDS(file="./data/5E_tax_table.RDS") #READ IN RDS B/C FASTER
TAX<-tax_table(tax_file)

## reading in META DATA
meta_data <- read.table(file="./data/5E_METADATA.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)
meta_data$trt <- factor(meta_data$trt, levels = c("No inoculation","Disrupted","Live"))
meta_data$geno_trt <- factor(meta_data$geno_trt, levels = c("FPsc_no.inoculation","FPsc_disrupted","FPsc_live","R500_no.inoculation","R500_disrupted","R500_live"))
META_DATA <- sample_data(meta_data)

#DOUBLE CHECKING ORDER OF DATA AND NAME MATCHES
sample_names(META_DATA) == sample_names(OTU)

## MERGING INTO PHYLOSEQ OBJECT
data <- merge_phyloseq(OTU,TAX,META_DATA)

## DATA OVERVIEW
ntaxa(data)[1]
nsamples(data)
sample_names(data)
sample_variables(data)

## LOOKING AT READ NUMBERS BY SAMPLE
histogram(sort(sample_sums(OTU)))
head(sort(sample_sums(data)))  ##46406 is lowest number of reads
summary((sort(sample_sums(data))))

## LOOKING AT LEVEL OF CHL AND MT CONTAMINATION
rank_names(data)  ##double check naming hierarchy in Taxonomy
sort(get_taxa_unique(data,"Family"))  ## looking for contaminants

## subset out contaminants to see how prevalent they are
mitochondria_contaminants <- subset_taxa(data, Family == "Mitochondria")
chloroplast_contaminants <- subset_taxa(data, Class == "Chloroplast")

archaea_contaminants <- subset_taxa(data, Kingdom == "Archaea")
eukaryota_contaminants <- subset_taxa(data, Kingdom == "Eukaryota")
# length(tax_table(archaea_contaminants))


##REMOVING UNWANTED / CONTAMINANTS
##Note I DECIDED TO LEAVE ARCHAEA IN THE ANALYSIS
data <- subset_taxa(data, (Class!="Chloroplast") | is.na(Class))  ##is.na(Class) retains NAs
data <- subset_taxa(data, (Family!= "Mitochondria") | is.na(Family))
data <- subset_taxa(data, (Kingdom == "Bacteria") | Kingdom == "Archaea")


##RAREFACTION;  There are fancier ways to analyze these data w/o throwing reads away
##  but for a 30000 ft view of the data, testing if treatments are different
##  from ea. other, this should do just fine.
head(sort(sample_sums(data))) ## lowest number of reads is 45956
## mtb set seed for reproducability and is sampling W/O Replacement as suggested in 
## WEISS et al. 2017
data_rare <- rarefy_even_depth(data, sample.size= 45956, rngseed=(12), replace=F)
# sort(sample_sums(data_rare)) ## double check

## SIMPLE VERSION OF CORE COMMUNITY BASED ON LUNDBERG et al. 2012
data_rare_25x5 <- filter_taxa(data_rare,function(x) sum(x>25) > 5,prune=T)


### ALPHA DIVERSITY METRICS 
## ALPHA DIVERSITY STATS
a_richness_est <- estimate_richness(data_rare, measures=c("Observed","Shannon","Simpson"))
data_META <- data.frame(sample_data(data_rare))
row.names(a_richness_est) == row.names(data_META)

## CHANGE ALPHA METRIC AMONG: "Observed","Shannon","Simpson"
richness.lm <- lm(a_richness_est$Simpson ~ blk + trt + genotype + genotype:trt, data=data_META)
summary(aov(richness.lm))

pgp1_lf_ln_emm <- emmeans(richness.lm, "trt")
pgp1_lf_ln_emm
pairs(pgp1_lf_ln_emm)





## ALPHA DIVERSITY PLOTS
theme_set(theme_bw())
richness.plot <- plot_richness(data_rare, x = "trt", measures=c("Observed", "Shannon"), color='trt', shape='genotype')
richness.plot <- richness.plot + geom_point(size=3) + theme(strip.text.x=element_text(size=11))
richness.plot <- richness.plot + labs(color="Soil Treatment", shape ="Genotype")
richness.plot <- richness.plot + scale_color_manual(values = c('#0d2c54','#0d7f96','#ADD12A','#0d2c54','#0d7f96','#ADD12A'))
richness.plot <- richness.plot + theme(axis.title.x=element_blank())
richness.plot


# mtb_richness <- estimate_richness(data_rare, measures=c("Observed", "Shannon", "Simpson"))
# row.names(mtb_richness) == row.names(meta_data)
# mtb_richness <- cbind(mtb_richness,meta_data)
# mtb.obs <- ggplot(mtb_richness, aes(x=trt,y=Observed,color=trt, shape=genotype)) + geom_point()
# mtb.obs
# mtb.shan <- ggplot(mtb_richness, aes(x=trt,y=Shannon,color=trt, shape=genotype)) + geom_point()
# mtb.shan
# mtb.simp <- ggplot(mtb_richness, aes(x=trt,y=Simpson,color=trt, shape=genotype)) + geom_point()
# mtb.simp
# gridExtra::grid.arrange(mtb.obs,mtb.shan,mtb.simp,ncol=3)


# mtb_richness_long <- gather(mtb_richness,Measurement,Value,Observed,Shannon,Simpson)
# richness.plot <- ggplot(mtb_richness_long, aes(x=trt, y=Value, color=trt)) + geom_point(na.rm=T)
# richness.plot <- richness.plot + facet_wrap(~Measurement,scales="free_y")
# richness.plot
# dat_text <- data.frame(
#  label=c("A","A","B"),
#  Measurment=c("Observed","Observed","Observed"),
#  trt=c("No inoculation","Disrupted","Live"),
#  y=c(500,600,1251)
# )
# richness.plot <- richness.plot + geom_text(data=dat_text,mapping=aes(x=trt,y=y,label=label))
# richness.plot

## INDIVIDUAL PLOT OF ALPHA DIVERSITY
## NOW A PANEL IN FIG_SXXX
# pdf(file="./figures/plot_5E_alpha_new.pdf", width=6, height=4, paper="letter", useDingbats=FALSE)
# richness.plot
# dev.off()

 
richness.plot.geno <- plot_richness(data_rare, x = "genotype", measures=c("Observed", "Shannon", "Simpson"), color='genotype')
richness.plot.geno <- richness.plot.geno + geom_point(size=3)
richness.plot.geno <- richness.plot.geno + scale_color_manual(values = c('#d9601a','#0396a6'))

richness.plot.geno


### BETA DIVERSITY METRICS
## BETA DIVERSITY STATS
data_bc <-phyloseq::distance(data_rare, method = "bray")
data_META <- data.frame(sample_data(data_rare))
adonis(data_bc ~ blk + trt + genotype + genotype:trt, data=data_META)

pairwise.adonis(data_bc, factors= data_META$trt, p.adjust="bonferroni")


data_bj <-phyloseq::distance(data_rare, method = "jaccard", binary=T)
data_META <- data.frame(sample_data(data_rare))
adonis(data_bj ~ blk + trt + genotype + genotype:trt, data=data_META) 

head(sample_data(data_rare))
colnames(sample_data(data_rare)[,2]) <- "Treatment"


## BETA DIVERSITY PLOTS
## BRAY CURTIS DISTANCE MATRIX
bc.ord <- ordinate(data_rare, method="PCoA", distance="bray")
bc.plot <- plot_ordination(data_rare, bc.ord, color="trt", shape="genotype")
bc.plot <- bc.plot + labs(color="Soil Treatment", shape="Genotype")
bc.plot <- bc.plot + scale_colour_manual(values = c('#0d2c54','#0d7f96','#ADD12A'))
bc.plot <- bc.plot + geom_point(size=3)
bc.plot <- bc.plot + stat_ellipse(type = "norm", aes(group=trt))
bc.plot <- bc.plot + theme_bw()
bc.plot

# pdf(file="./figures/plot_5E_beta_bc.pdf", width=7, height=5)
# bc.plot
# dev.off()

## BINARY JACCARD DISTANCE MATRIX
bj.ord <- ordinate(data_rare, method="PCoA", distance="jaccard", binary=T)
bj.plot <- plot_ordination(data_rare, bj.ord, color="trt", shape="genotype")
bj.plot <- bj.plot + scale_colour_manual(values = c('#0d2c54','#0d7f96','#ADD12A'))
bj.plot <- bj.plot + geom_point(size=3)
bj.plot <- bj.plot + stat_ellipse(type = "norm", aes(group=trt))
bj.plot


###STACKED BAR PLOT AT PHYLUM LEVEL
data_rel_abund <- read.csv(file="./data/5E_phylum_rel_abund_stacked.csv")
data_rel_abund$Phylum <- factor(data_rel_abund$Phylum, levels = c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Cyanobacteria","Firmicutes","Planctomycetes","Proteobacteria","Verrucomicrobia","< 1% abundance"))
data_rel_abund$trt <- factor(data_rel_abund$trt, levels = c("No Inoculation", "Disrupted", "Live"))

mtb_palette <- c("#ADD12A","#332288","#cc6677","aquamarine4",
                 "#661100","#ddcc77","#aa4499","#88ccee","#44aa99","#882255")
#mtb_palette <- mtb_palette[sample(1:10)]
#mtb_palette
stacked.plot <- ggplot(data=data_rel_abund, aes(x=sample_3, y=Abundance, fill=Phylum)) +
  facet_grid(~trt, scales = "free") + theme(strip.text.x=element_text(size=11))
stacked.plot <- stacked.plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mtb_palette) +
  theme(legend.position="right", axis.text.x = element_text(angle = -90), strip.text=element_text(size=6)) + guides(fill=guide_legend(ncol=1)) +
  theme(axis.title.x=element_blank())
stacked.plot

# pdf(file="./figures/plot_5E_stacked_bar.pdf", width=7, height=5)
# stacked.plot
# dev.off()


## COMBINE PANELS FOR FIGURE_SXXX
fig_sXXX <- plot_grid(stacked.plot,richness.plot, ncol=2, labels="AUTO")

# pdf(file="./figures/fig_sXXX.pdf", height=4, width=12)
# fig_sXXX
# dev.off()

fig_sXXX_vert <- plot_grid(stacked.plot,richness.plot, ncol=1, labels="AUTO")

# pdf(file="./figures/Fig_S3_microbe_diversity.pdf", height=8, width=6)
# fig_sXXX_vert
# dev.off()

# png(file="./figures/Fig_S3_microbe_diversity.png", height=4800, width=3600, res=600)
# fig_sXXX_vert
# dev.off()

