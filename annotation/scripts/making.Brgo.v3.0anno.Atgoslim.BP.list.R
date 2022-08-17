# common tools for Brassica microbiome projects
## bug fixed. GO terms in each gene are sometimes redundant
## eg. BraA01g000010.3C has three GO:0006355.
## so remove the redundancy (051420)
## modify thi sscript for Br_Mirobe_Paper_2021 (081622)
## home dir is this file location
  library(tidyverse)
  ## making Brgo.v3.0anno.Atgoslim.BP.list (run once)
  ### prerequisit (bug fixed with max_target_seq option problem in blastn)
  # Br.v3.0anno.At.BLAST<-readr::read_csv(file.path("Annotation","input","v3.0annotation","Brapa3.0_vs_At_dc-megablast_out.csv"),col_names=FALSE) # BLASTed by Julin
  # colnames(Br.v3.0anno.At.BLAST) <- c("query","subject","perc_ID","aln_length","mismatch","gap_open","qstart","qend","sstart","send","eval","score")
  # Br.v3.0anno.At.BLAST <- Br.v3.0anno.At.BLAST %>% tidyr::separate(subject, c("AGI","splicing"),sep="\\.") #%>% select(AGI)
  # head(Br.v3.0anno.At.BLAST)
  # summary(Br.v3.0anno.At.BLAST)
  # 
  # read annotated version
  Br.v3.0anno.At.BLAST<-readr::read_csv(file.path("..","output","Brapa_V3.0_annotated.csv"),col_names=TRUE) # BLASTed by Julin
  # there are multiple Br genes found
  Br.v3.0anno.At.BLAST %>% group_by(name) %>% summarise(n=n()) %>% arrange(desc(n))
  # pick only one with highest score within the same Br gene name
  Br.v3.0anno.At.BLAST.highscore <- Br.v3.0anno.At.BLAST %>% group_by(name) %>% arrange(desc(score)) %>% slice(1)
  Br.v3.0anno.At.BLAST.highscore %>% group_by(name) %>% summarise(n=n()) %>% arrange(desc(n))
  
  # read At GO list
  load(file.path("..","input","Atgoslim.TAIR.BP.list.Rdata"))
  # asign At GO into corresponding Br genes
  Brgo.v3.0anno.Atgoslim.BP.list<-list()

  for(i in 1:length(Br.v3.0anno.At.BLAST.highscore$name)) {
    if(is.null(Atgoslim.TAIR.BP.list[[as_vector(Br.v3.0anno.At.BLAST.highscore[i,"AGI"])]])) next else {
      Brgo.v3.0anno.Atgoslim.BP.list[[i]]<-Atgoslim.TAIR.BP.list[[as_vector(Br.v3.0anno.At.BLAST.highscore[i,"AGI"])]]
      names(Brgo.v3.0anno.Atgoslim.BP.list)[[i]]<-as_vector(Br.v3.0anno.At.BLAST.highscore[i,"name"])
    }
  }
  table(sapply(Brgo.v3.0anno.Atgoslim.BP.list,is.null))
  Brgo.v3.0anno.Atgoslim.BP.list<-Brgo.v3.0anno.Atgoslim.BP.list[!sapply(Brgo.v3.0anno.Atgoslim.BP.list,is.null)]
  table(sapply(Brgo.v3.0anno.Atgoslim.BP.list,is.null))
  # bug fix (051420)
  Brgo.v3.0anno.Atgoslim.BP.list[[1]] %>% duplicated() %>% table() # False 3, TRUE 2
  Brgo.v3.0anno.Atgoslim.BP.list[[1]] %>% unique() %>% duplicated() %>% table() # False 3
  # so
  Brgo.v3.0anno.Atgoslim.BP.list.unique <- Brgo.v3.0anno.Atgoslim.BP.list %>% map(.,unique)
  Brgo.v3.0anno.Atgoslim.BP.list.unique[[1]] %>% duplicated() %>% table() # False 3
  
  save(Brgo.v3.0anno.Atgoslim.BP.list.unique,file=file.path("..","output","Brgo.v3.0anno.Atgoslim.BP.list.unique.Rdata"))



  