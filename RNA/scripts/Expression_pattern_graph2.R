# required libraries
library(edgeR)
library(tidyverse)
library(readr)
library(readxl)
library(cowplot) # for plotting both genotypes or density

# load reads mapped to Brassica genome (either v1.5 annotation or v3.0 annotation)
getwd()
## for exp1 v3.0annotation
# counts.exp1.v3.0 <- readr::read_csv(file.path("RNA","input","wyo001_V3.0_raw_counts.csv.gz"),col_names=TRUE)
### cpm
#cpm.exp1.leaf.v3.0 <- readr::read_csv(file.path("RNA","v3.0annotation","20170617-samples","output","cpm_wide_20170617_leaf_samples.csv.gz"),col_names=TRUE)
#cpm.exp1.root.v3.0 <- readr::read_csv(file.path("..","v3.0annotation","20170617-samples","output","cpm_wide_20170617_root_samples.csv.gz"),col_names=TRUE)

### voom
# plus R500 in residual
root.voom5 <- readr::read_tsv(file.path("RNA","output","voom_expression.e1and3.resid.exp_gt.root.plusR500.txt.gz"))
leaf.voom5 <- readr::read_tsv(file.path("RNA","output","voom_expression.e1and3.resid.exp_gt.leaf.plusR500.txt.gz"))

#counts.exp1.v3.0 # make sure this is v3.0 annotation (look target_id column)
## for exp3 v3.0annotation
#counts.exp3.v3.0 <- readr::read_csv(file.path("..","v3.0annotation","20180202-samples","input","20180202_V3.0_raw_counts.csv.gz"),col_names=TRUE)
#counts.exp3.v3.0
### cpm
#cpm.exp3.leaf.v3.0<-readr::read_csv(file.path("..","v3.0annotation","20180202-samples","output","cpm_wide_20180202_leaf_samples.csv.gz"),col_names=TRUE)
#cpm.exp3.root.v3.0<-readr::read_csv(file.path("..","v3.0annotation","20180202-samples","output","cpm_wide_20180202_root_samples.csv.gz"),col_names=TRUE)

# sample files
# exp1 (20170617-samples)
sample.description.e1<-readr::read_csv(file.path("plant","output","Br.mbio.e1.sample.description.csv")) # does not work
# exp3 (20180202-samples)
sample.description.e3<-readr::read_csv(file.path("plant","output","Br.mbio.e3.sample.description.csv"))

### copy and paste from 13_WGCNA_merged_moduleColors_signedhybrid.Rd below
sample.description.e1 %>% summarize(n_distinct(group))
sample.description.e1 %>% group_by(trt) %>% summarize(n())
##
root.voom5.e1 <- tibble(root.voom5 = colnames(root.voom5)[-1]) %>% inner_join(sample.description.e1[,c("sample","genotype","tissue","trt")],by=c("root.voom5"="sample"))
root.voom5.e3 <- tibble(root.voom5 = colnames(root.voom5)[-1]) %>% mutate(genotype="FPsc") %>% inner_join(sample.description.e3[,c("sample","tissue","trt")],by=c("root.voom5"="sample"))  
root.voom5.sample <- bind_rows(root.voom5.e1,root.voom5.e3) %>% mutate(trt=str_remove(trt,"5E_"))
# for Whitney
# root.voom5.sample <- bind_rows(root.voom5.e1,root.voom5.e3)
# root.voom5.sample$trt <- gsub("5E_","",root.voom5.sample$trt)
root.voom5.sample %>% group_by(genotype) %>% summarise(num=n())


# leaf.voom5.sample
leaf.voom5.e1 <- tibble(leaf.voom5 = colnames(leaf.voom5)[-1]) %>% inner_join(sample.description.e1[,c("sample","genotype","tissue","trt")],by=c("leaf.voom5"="sample"))
leaf.voom5.e3 <- tibble(leaf.voom5 = colnames(leaf.voom5)[-1]) %>% mutate(genotype="FPsc") %>% inner_join(sample.description.e3[,c("sample","tissue","trt")],by=c("leaf.voom5"="sample"))  
leaf.voom5.sample <- bind_rows(leaf.voom5.e1,leaf.voom5.e3) %>% mutate(trt=str_remove(trt,"5E_")) # str_remove did not work in Whitney (041721)
# for Whitney
# leaf.voom5.sample <- bind_rows(leaf.voom5.e1,leaf.voom5.e3)
# leaf.voom5.sample$trt <- gsub("5E_","",root.voom5.sample$trt)
leaf.voom5.sample %>% group_by(genotype) %>% summarise(num=n())

### copied from 13_WGCNA_merged_moduleColors_signedhybrid.Rd


# functions for drawing expression pattern

# using voom data (tissue specific)
#### under construction ####
expression.pattern.Br.graph.exp1and3.v3.0annotation.voom <- 
  function(target.genes="",
           title.e1=" ",
           title.e3=" ",
           tissue.type="leaf",
           FDR.e1, FDR.e3,# "N/A" if not available.
           geno=c("R500","FPsc","both"),
           exp=c("exp1and3","exp1","exp3")){
  if (tissue.type=="leaf") {
      data <- leaf.voom5
    } else {data <- root.voom5}
  print(paste("data is",data[1:10,]))
  print(paste("tissue.type is",tissue.type))
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  # data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
  #   inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type)
  # 
  # using FDR
  # if target_id exist (under construction)
  # if(data.temp$target_id) 
  data.temp.e1 <- data  %>%
    filter(genes %in% target.genes) %>% 
    pivot_longer(cols=-genes,names_to="sample",values_to="value") %>%
    #gather(sample,value,-genes) %>%
    inner_join(sample.description.e1, by="sample") %>% #View()
    filter(tissue==tissue.type) 
    # add FDR info (if any)
  if(FDR.e1=="N/A") {data.temp.e1 <- data.temp.e1} else {
    print("adding FDR.e1 to data.temp.e1.")
    data.temp.e1 <- data.temp.e1 %>%
    inner_join(FDR.e1, by="genes") %>% mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(genes.FDR,genes,FDR,sep="\n FDR ")
    data.temp.e1
  }
  data.temp.e3 <- data  %>%
    filter(genes %in% target.genes) %>% 
    pivot_longer(cols=-genes,names_to="sample",values_to="value") %>%
    #gather(sample,value,-genes) %>%
    inner_join(sample.description.e3, by="sample") %>% 
    mutate(genotype="FPsc") %>%
    filter(tissue==tissue.type)  
  if(FDR.e3=="N/A") {data.temp.e3 <- data.temp.e3} else {
    print("adding FDR.e3 to data.temp.e3.")
    data.temp.e3 <- data.temp.e3 %>%
      inner_join(FDR.e3, by="genes") %>% mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
      unite(genes.FDR,genes,FDR,sep="\n FDR ")
  }
  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    if(FDR.e1=="N/A") {
      p.e1 <- data.temp.e1 %>% ggplot(aes(x=trt,y=value))  + 
        geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2,size=3)  + 
        theme_bw() +
        facet_grid(genes~genotype,scales="free") + 
        theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
        theme(legend.position="bottom") + labs(title=title.e1)
    } else {
    p.e1 <- data.temp.e1 %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(genes.FDR~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=title.e1)
    }
    if(FDR.e3=="N/A") {
    p.e3 <- data.temp.e3 %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(genes~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=title.e3)
    } else {
    p.e3 <- data.temp.e3 %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt,shape=as.character(block)),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(genes.FDR~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=title.e3)
    }
    # plot
    if(exp=="exp1and3") {
      # combine e1.p and e3.p by cowplot
        p <- plot_grid(p.e1, p.e3) 
      return(p)
      } else if(exp=="exp1") {
        p <- p.e1
        return(p)
    } else if(exp=="exp3") {
       p <- p.e3
       return(p)
    } else {print("Which experiment?");stop}
  } else if(geno=="FPsc"|geno=="R500") {
    p.e1 <- data.temp.e1 %>% filter(genotype==geno) %>% 
      ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt))  
    + theme_bw() + facet_grid(genes~.,scales="free") 
    + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) 
    + theme(legend.position="bottom") 
    + labs(title=title.e1)
    p.e3 <- data.temp.e3 %>% filter(genotype==geno) %>% 
      ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt))  
    + theme_bw() + facet_grid(genes~.,scales="free") 
    + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) 
    + theme(legend.position="bottom") 
    + labs(title=title.e3)    
    # plot
    if(exp=="exp1and3") {
      # combine e1.p and e3.p by cowplot
      p <- plot_grid(p.e1, p.e3) 
      p
      return(p)
    } else if(exp=="exp1") {
      p <- p.e1
      return(p)
    } else if(exp=="exp3") {
      p <- p.e3
      return(p)
    } else {print("Which experiment?");stop}
  }  else {print("Specify genotype.");stop}
}

# 
# using voom data (tissue specific), exp1 and 3
expression.pattern.Br.graph.exp1and3.v3.0annotation.voom2 <- function(target.genes="",FDR="",sample.description="",title="",tissue.type="leaf",geno=c("R500","FPsc","both")){
  if (tissue.type=="leaf") {
    data <- leaf.voom5
    sample.description <- leaf.voom5.sample
    data.temp <- data  %>%
      filter(genes %in% target.genes) %>% 
      pivot_longer(cols=-genes,names_to="sample",values_to="value") %>%
      #gather(sample,value,-genes) %>%
      inner_join(sample.description, by=c("sample"="leaf.voom5"))
  } else {
    data <- root.voom5
    sample.description <- root.voom5.sample
    data.temp <- data  %>%
      filter(genes %in% target.genes) %>% 
      pivot_longer(cols=-genes,names_to="sample",values_to="value") %>%
      #gather(sample,value,-genes) %>%
      inner_join(sample.description, by=c("sample"="root.voom5"))
  }
  print(paste("data is",data[1:10,]))
  print(paste("tissue.type is",tissue.type))
  
  data[is.na(data)] <- 0 #
  # select genes and add sample info
  # data.temp<-data %>% filter(target_id %in% target.genes) %>% gather(sample,value,-target_id) %>%
  #   inner_join(sample.description, by="sample") %>% filter(tissue==tissue.type)
  # 
  # using FDR
  data.temp <- data.temp %>% 
    filter(tissue==tissue.type) %>%
    right_join(FDR, by="genes") %>% # add FDR info
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(genes.FDR,genes,FDR,sep="\n FDR ")
  
  # 
  if(geno=="both") { # needs to impove this
    # ggplot(data.temp, aes(x=genotype,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt,shape=tissue) )  + theme_bw() + facet_grid(target_id~tissue,scales="free") + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) + theme(legend.position="bottom") + labs(title=title)
    p<-data.temp %>% ggplot(aes(x=trt,y=value))  + 
      geom_jitter(alpha = 0.5,aes(colour=trt),width=0.2,size=3)  + 
      theme_bw() +
      facet_grid(genes.FDR~genotype,scales="free") + 
      theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
      theme(legend.position="bottom") + labs(title=title)
    p
  } else if(geno=="FPsc"|geno=="R500") {
    data.temp %>% filter(genotype==geno) %>% 
      ggplot(aes(x=trt,y=value))  + geom_jitter(alpha = 0.5,aes(colour=trt))  
    + theme_bw() + facet_grid(genes~.,scales="free") 
    + theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) 
    + theme(legend.position="bottom") 
    + labs(title=title)
  }  else {print("Specify genotype.");stop}
}

# Arabidopsis RNAseq version
expression.pattern.AT.graph.cpm<-function(data="",sample.description="",target.genes.FDR="",title=""){
  data[is.na(data)] <- 0 #
  # using FDR
  data.temp<-data  %>% 
    dplyr::rename(target_id=AGI) %>%
    filter(target_id %in% target.genes.FDR$AGI) %>% dplyr::select(-transcript_ID,-variant) %>%
    gather(sample,value,-target_id) %>%
    inner_join(sample.description, by="sample")  %>%
    right_join(target.genes.FDR[,c("AGI","FDR")], by=c("target_id"="AGI")) %>%
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(target_id.FDR,target_id,FDR,sep="\n FDR ") 
  # print(data.temp)
  p<-data.temp %>% ggplot(aes(x=group,y=value))  + 
    geom_jitter(alpha = 0.5,aes(color=group))  + 
    theme_bw() + facet_grid(target_id.FDR~.,scales="free") +
    theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
    theme(legend.position="none") + labs(title=title)
  p
}


# Arabidopsis ATH1 version (under construction)
expression.pattern.AT.graph.ATH1<-function(data="",sample.description="",target.genes.FDR="",title=""){
  data[is.na(data)] <- 0 #
  # using FDR
  data.temp<-data  %>% 
#    dplyr::rename(target_id=AGI) %>%
    filter(AGI %in% target.genes.FDR$AGI) %>%
    # under construction #####
    # gather(sample,value,-target_id) %>%
    # inner_join(sample.description, by="sample")  %>%
    right_join(target.genes.FDR[,c("AGI","FDR")], by="AGI") %>%
    mutate(FDR=format(FDR,digits=2,scientific=TRUE)) %>%
    unite(AGI.FDR,AGI,FDR,sep="\n FDR ") %>%
    mutate(trt=fct_relevel(trt))
  # print(data.temp)
  p<-data.temp %>% ggplot(aes(x=trt,y=value))  + 
    geom_jitter(alpha = 0.5,aes(color=trt),width=0.25)  + 
    theme_bw() + facet_grid(AGI.FDR~.,scales="free") +
    theme(strip.text.y=element_text(angle=0),axis.text.x=element_text(angle=90)) +
    theme(legend.position="none") + labs(title=title)
  p
}

