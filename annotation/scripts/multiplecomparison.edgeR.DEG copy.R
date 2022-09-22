# simple DEG function for RNAseq
multiplecomparison.edgeR.DEG<-function(counts="",
                                    group="",
                                    comparison.num=13,
                                    DEGfilename="",
                                    dge.object.filepath.name="",
                                    expression.plot="") {
  # format counts into counts2
  counts2 <- counts %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "target_id") %>%
    as.matrix() %>%
    round(0)
  # dge 
  dge <- DGEList(counts=counts2, 
               samples=colnames(counts2), 
               group=group)
  # normalize
  dge <- calcNormFactors(dge)
  # save
  save(dge,file=dge.object.filepath.name)
  # DEG
    design <- model.matrix(~group) 
    y <- estimateDisp(dge,design)
    fit <- glmFit(y,design)
    
    #### under construction #####
    # glmLRT for multiple comparison 
    lrt.map <- map(2:comparison.num,function(x) glmLRT(fit,coef=x))
    #lrt.list <- lapply(2:comparison.num,function(x) glmLRT(fit,coef=x))
    names(lrt.map) <- levels(group)[-1]
    ###
    DEGs.map <- map(lrt.map,function(x) topTags(x,n = Inf,p.value = 1)$table)
  # annotation
    DEGfilename.lists <- paste("minus",levels(group)[-1],"_",DEGfilename,sep="") %>% as.list()
   # pmap(DEGs.map ,DEGfilename.lists,function(x,y) addAGIAnno(DEG=x,DEGfilename2=y)) # does not work
  for(i in 1:length(DEGfilename.lists)) {
    addAGIAnno(DEG=DEGs.map[[i]],DEGfilename2=DEGfilename.lists[[i]]) # default file path is file.path("custom_categories")
  }
# Expression pattern
# load saved DEG files
# csv (full name)
    DGE.file <- list.files(path=file.path("custom_categories"),full.names=TRUE,recursive=TRUE, pattern=DEGfilename)
    DGE.file
    # read csv file
    DGEs<-lapply(DGE.file, function(x) read_csv(file=file.path(x)))
    # name
    DGE.file2<- list.files(path=file.path("custom_categories"),full.names=FALSE,recursive=TRUE, pattern=DEGfilename)
    names(DGEs) <- DGE.file2 %>% str_remove_all("_root.DEGs.csv.gz")
  # load dge object for calculating cpm
    load(dge.object.filepath.name)
    # # Expression pattern coef=2
  #   DEGs1.all.anno <- read_csv(file.path("custom_categories",DEG1filename))
  expression.cpm <- cpm(dge) %>% as_tibble() %>% bind_cols(data.frame(transcript_ID=rownames(dge$counts)),.) %>% separate(transcript_ID,into=c("AGI","variant"),remove=FALSE)
     sample.description<-dge$samples %>% as_tibble() %>% dplyr::rename(sample=samples)
  # for loop method
  for(i in 1:length(DGE.file2)) {
    DEG.title = names(DGEs)[i]
  # upregulated genes
     gene.of.interest.FDR.up <- DGEs[[i]] %>% filter(FDR< 0.05,logFC>0) %>% dplyr::select(AGI,logFC,FDR) %>% arrange(FDR)
       print("logFC>0")
       if(dim(gene.of.interest.FDR.up)[1]==0) {print(paste("Skip ",DEG.title,".up",sep=""));next} else {
       
     p.up<-expression.pattern.AT.graph.cpm(data=expression.cpm,target.genes.FDR=gene.of.interest.FDR.up[1:5,],sample.description=sample.description) + labs(title=paste(DEG.title,".up",sep=""))
     print(p.up)   
   ggsave(p.up,file=file.path("custom_categories",paste(names(DGEs)[i],".up.plot.png",sep="")))
  }
  # # down-regulated genes
      gene.of.interest.FDR.down <- DGEs[[i]] %>% filter(FDR< 0.05,logFC<0) %>% dplyr::select(AGI,logFC,FDR) %>% arrange(FDR)
      print("logFC<0")
      if(dim(gene.of.interest.FDR.up)[1]==0) {print(paste("Skip ",DEG.title,".down",sep=""));next} else {
      p.down<-expression.pattern.AT.graph.cpm(data=expression.cpm,target.genes.FDR=gene.of.interest.FDR.down[1:5,],sample.description=sample.description) + labs(title=paste(DEG.title,".down",sep=""))
      print(p.down)   
      ggsave(p.down,file=file.path("custom_categories",paste(names(DGEs)[i],".down.plot.png",sep="")))
      }
  }
# under construction #####
 
# tidyverse method (map? nest?)
      
}


