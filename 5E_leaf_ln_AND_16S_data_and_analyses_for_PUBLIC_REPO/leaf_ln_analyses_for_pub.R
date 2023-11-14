setwd("../5E_leaf_ln_AND_16S_data_for_public_repo")
library(emmeans)
library(ggplot2)
library(cowplot)

#### READ IN DATA FROM EXP 1 (AKA PGP1_REDUX)
data_exp1 <- read.csv(file="./data/pot_level_leaf_ln_PGP1_redux_to_JM.csv")
head(data_exp1)

## CHANGING SEVERAL NUMERIC VARIABLES TO FACTORS AND GETTING THE ORDER OF SOIL TREATMENTS RIGHT FOR FIGURES
data_exp1$bench <- as.factor(data_exp1$bench)
data_exp1$block <- as.factor(data_exp1$block)
data_exp1$X5E_soil_trts <- factor(data_exp1$X5E_soil_trts, levels=c("No inoculation", "Disrupted","Live"))

## ANOVA MODEL FOR AVERAGE LEAF LENGTH 
pgp1_lf_ln <- lm(avg_lf_03.24_no_germ ~bench + bench/block + genotype + X5E_soil_trts + genotype:X5E_soil_trts, data=data_exp1)
anova(pgp1_lf_ln)

##  LSMEANS (NOW EMMEANS) OF THREE SOIL TREATMENTS IN A POST HOC TUKEY'S TEST
pgp1_lf_ln_emm <- emmeans(pgp1_lf_ln, "X5E_soil_trts")
pgp1_lf_ln_emm
pairs(pgp1_lf_ln_emm)

##  LSMEANS (NOW EMMEANS) OF GENOTYPES
pgp1_lf_ln_emm.geno <- emmeans(pgp1_lf_ln, "genotype")
pgp1_lf_ln_emm.geno

## GET RESIDUALS OF LEAF_LN AFTER REMOVING SPATIAL BLOCK EFFECTS FOR ALTERNATIVE FIG1A
pgp1_lf_ln_noblk_nogeno <- lm(avg_lf_03.24_no_germ ~bench + bench/block + genotype, data=data_exp1)
anova(pgp1_lf_ln_noblk_nogeno)

pgp1_resid <- residuals(lm(avg_lf_03.24_no_germ ~bench + bench/block + genotype, data=data_exp1))
length(pgp1_resid)
pgp1_mean <- mean(data_exp1$avg_lf_03.24_no_germ, na.rm=T)
pgp1_resid <- pgp1_resid+pgp1_mean ## add residuals to mean leaf_ln of all plants to get back on scale
data_exp1 <- cbind(data_exp1, pgp1_resid)


#### READ IN DATA FROM FOLLOW.UP EXPERIMENT (AKA DENSITY EXP)
data_validation <- read.csv(file="./data/pot_level_leaf_ln_followup_exp.csv")
head(data_validation)


## CHANGING SEVERAL NUMERIC VARIABLES TO FACTORS AND GETTING THE ORDER OF SOIL TREATMENTS RIGHT FOR FIGURES
data_validation$bench <- as.factor(data_validation$bench)
data_validation$block <- as.factor(data_validation$block)
data_validation$X5E_soil_trts <- factor(data_validation $X5E_soil_trts, levels=c("Disrupted","Live"))

## ANOVA MODEL FOR AVERAGE LEAF LENGTH 
valid_lf_ln <- lm(raw_avg_lf_11.07.17 ~ bench + bench/block + X5E_soil_trts, data= data_validation)
anova(valid_lf_ln)

##  LSMEANS (NOW EMMEANS) OF THREE SOIL TREATMENTS IN A POST HOC TUKEY'S TEST
valid_lf_ln_emm <- emmeans(valid_lf_ln, "X5E_soil_trts")
valid_lf_ln_emm
pairs(valid_lf_ln_emm)

## GET RESIDUALS OF LEAF_LN AFTER REMOVING SPATIAL BLOCK EFFECTS FOR ALTERNATIVE FIG1A
valid_lf_ln_no_blk <- lm(raw_avg_lf_11.07.17 ~ bench + bench/block, data= data_validation, na.action=na.exclude)
anova(valid_lf_ln_no_blk)
summary(valid_lf_ln_no_blk)

valid_lf_ln_mean <- mean(data_validation$raw_avg_lf_11.07.17, na.rm=T)

valid_lf_ln_resid <- residuals(lm(raw_avg_lf_11.07.17 ~ bench + bench/block, data=data_validation, na.action=na.exclude))
valid_lf_ln_resid <- valid_lf_ln_resid + valid_lf_ln_mean
data_validation  <- cbind(data_validation, valid_lf_ln_resid)


## COMBINING EMM DATA FOR FIG1
combined_data <- summary(pgp1_lf_ln_emm)[,1:3]
combined_data <- rbind(combined_data, summary(valid_lf_ln_emm)[1:2,1:3])
combined_data
combined_data$exp_group <- c("PGP Exp.","PGP Exp.","PGP Exp.","Validation Exp.", "Validation Exp.")
combined_data

lf_ln_plot <- ggplot(combined_data, aes(x=exp_group, y=emmean, fill=X5E_soil_trts)) 
lf_ln_plot <- lf_ln_plot + geom_bar(stat="identity", color="black", position=position_dodge()) + theme_bw()
lf_ln_plot <- lf_ln_plot + ylim(0,5.1)
lf_ln_plot <- lf_ln_plot + geom_errorbar(aes (ymin=emmean-SE, ymax=emmean+SE), width=0.2, position=position_dodge(0.9))
lf_ln_plot <- lf_ln_plot + scale_fill_manual(values=c('#0d2c54','#0d7f96','#ADD12A'))
lf_ln_plot <- lf_ln_plot + ylab("Leaf Length (cm)")
lf_ln_plot <- lf_ln_plot + guides(fill=guide_legend("Soil Treatment"))
lf_ln_plot <- lf_ln_plot + theme(axis.title.x=element_blank(), legend.title=)
lf_ln_plot <- lf_ln_plot + geom_segment(aes(x=0.6, y=4, xend=1.1,yend=4))
lf_ln_plot <- lf_ln_plot + geom_segment(aes(x=1.2, y=4, xend=1.4, yend=4))
lf_ln_plot <- lf_ln_plot + annotate("text", x=1.15, y=4.1, label="***") 
lf_ln_plot <- lf_ln_plot + geom_segment(aes(x=1.6, y=5, xend=2.1,yend=5))
lf_ln_plot <- lf_ln_plot + geom_segment(aes(x=2.2, y=5, xend=2.4, yend=5))
lf_ln_plot <- lf_ln_plot+ annotate("text", x=2.15, y=5.1, label="***")
lf_ln_plot


head(data_validation)
dim(data_validation)
## COMBINING RESIDUAL DATA FOR FIG1
residuals_pgp1 <- data_exp1[,c(10,18)]
residuals_pgp1$exp_group <- "PGP Exp."
colnames(residuals_pgp1)[2] <- "leaf_ln_resid"
 
residuals_valid <- data_validation[,c(9,12)]
residuals_valid$exp_group <- "Validation Exp."
colnames(residuals_valid)[2] <- "leaf_ln_resid"

combined_resid <- rbind(residuals_pgp1, residuals_valid)
combined_resid_noNAs <- combined_resid[!is.na(combined_resid$leaf_ln_resid),]
tail(combined_resid_noNAs)
combined_resid_noNAs <- combined_resid_noNAs[!(combined_resid_noNAs$X5E_soil_trts=="No inoculation" & 
                                               combined_resid_noNAs$exp_group=="Validation Exp."),]

## DIFF TYPE OF PLOT TO MATCH OTHER IN MANUSCRIPT
lf_ln_plot_vs2 <- ggplot(combined_resid_noNAs, aes(x=exp_group, y=leaf_ln_resid, fill=X5E_soil_trts)) 
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + geom_boxplot(aes(middle = mean(leaf_ln_resid)), color="grey50") + theme_bw()
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + scale_fill_manual(values=c('#0d2c54','#0d7f96','#ADD12A'))
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + ylab("Leaf Length (cm)")
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + guides(fill=guide_legend("Soil Treatment"))
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + theme(axis.title.x=element_blank(), legend.title=)
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + geom_segment(aes(x=0.65, y=4.75, xend=1.05,yend=4.75))
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + geom_segment(aes(x=1.2, y=4.75, xend=1.365, yend=4.75))
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + annotate("text", x=1.125, y=4.8, label="***") 
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + geom_segment(aes(x=1.65, y=5.5, xend=2.05,yend=5.5))
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + geom_segment(aes(x=2.2, y=5.5, xend=2.3725, yend=5.5))
lf_ln_plot_vs2 <- lf_ln_plot_vs2 + annotate("text", x=2.125, y=5.55, label="***")
lf_ln_plot_vs2



## JOINING LEAF_LN_PLOT WITH BRAY_CURTIS_PLOT OF 16S
Fig_1 <- plot_grid(lf_ln_plot, bc.plot, ncol=2, labels="AUTO")
Fig_1

# pdf(file="./figures/Fig_1_new.pdf", height=4, width=10)
# Fig_1
# dev.off()

# png(file="./figures/Fig_1.png", height=2400, width=6000, res=600)
# Fig_1
# dev.off()


## JOINING LEAF_LN_PLOT WITH BRAY_CURTIS_PLOT OF 16S
Fig_1_vs2 <- plot_grid(lf_ln_plot_vs2, bc.plot, ncol=2, labels="AUTO")
Fig_1_vs2

# pdf(file="./figures/Fig_1_vs2_new.pdf", height=4, width=10)
# Fig_1_vs2
# dev.off()

# png(file="./figures/Fig_1_vs2.png", height=2400, width=6000, res=600)
# Fig_1_vs2
# dev.off()



