library(rhdf5)
library(gplots)
library(reshape2)
library(tidyr)
setwd('~/projects/Auto_Bionf/scLVM2/R')


dataCounts = h5dump('../data/normCountsBuettnerEtAl.h5f')
idx_het = which(dataCounts$genes_heterogen==1)

#log transformed and normalized expression matirx
Y = dataCounts$LogNcountsHum
countsMmusFilter= Y[idx_het,]
rownames(countsMmusFilter) = dataCounts$sym_names_het

nEspressed = colSums(dataCounts$countsMmusAll>0)
modelMat = model.matrix(~as.factor(dataCounts$labels))[,2:3]
ERCCPCs = prcomp(log2(t(dataCounts$nCountsERCC)+1))$x
Known = cbind(modelMat,nEspressed,ERCCPCs[,1:2])

colnames(Known) = c(paste("ccPhase",c("S", "G2M"),sep=""),'nExpressed', paste("PC",seq(1,2),sep=""))


#if you have custom gene sets add them here
#setsAdd = list()
#setsAdd[[1]]=custom_genes
#setsAdd_names = 'XXX'

write_scLVM2(countsMmusFilter, '../data/Buettneretal_.hdf5', Known=Known)



write_scLVM2 <- function(countsMmusFilter, filename, minGenes=15, Known=NULL, is_counts = FALSE, setsAdd=NULL, setsAdd_names=NULL){
  require(GSEABase)
  require(limma)
  fname = '../data/h.all.v5.0.symbols.gmt.txt'
  c2_set <- getGmt("../data/c2.cp.reactome.v4.0.symbols.gmt")
  files = Sys.glob('../data/wikipathways/Mm*txt')
  
  
  gene_lists = list()
  sets=c()
  i=1
  
  require(org.Mm.eg.db)

  x <- org.Mm.egSYMBOL
  xx <- as.list(x[mappedkeys(x)])
  
  for(file in files){
    pway = read.table(file, skip=1, sep='\t')
    gtype = as.character(pway[,2])
    
    if(sum(gtype=="Entrez Gene")/length(gtype)>.5 & sum(gtype=="Entrez Gene")>5){
      EG_names = as.character(pway[gtype=="Entrez Gene",1])
      
      gene_lists[[i]] = unlist(xx[EG_names])
      sets_ = unlist(strsplit(file,'_'))
      sets[i] = paste(sets_[3:(length(sets_)-2)],collapse = '_')
      i=i+1
    }
    
  }
  names(gene_lists) = sets
  
  sets_indicesWiki <- ids2indices(gene_lists, rownames(countsMmusFilter))
  sets_indicesWiki15 = sets_indicesWiki[unlist(lapply(sets_indicesWiki,length))>minGenes]
  
  gene_ids <- geneIds(c2_set)
  sets_indices <- ids2indices(gene_ids, toupper(rownames(countsMmusFilter)))

  sets_indices20 = sets_indices[unlist(lapply(sets_indices,length))>20]
  
  sym_names = rownames(countsMmusFilter)

  Yhet = t(countsMmusFilter)
  gene_list = list()
  idx_list = list()
  terms = c()
  for(i in seq(1,50)){
    line_i = scan(fname, sep='\t', flush=F,what = "list", ,multi.line = F, fill=T,nlines = 1, skip=i-1)  
    endl = length(line_i)
    gene_list[[i]] = unlist(line_i)[3:endl]
    idx_list[[i]] = na.omit(match(tolower(gene_list[[i]]), tolower(sym_names)))
    terms = c(terms, paste(strsplit(line_i[1],'_')[[1]][-1],collapse='_'))
  }
  
  if(!is.null(sets)){
    idx_list = c(idx_list,ids2indices(setsAdd, rownames(countsMmusFilter)))
    terms = c(terms, setsAdd_names)    
  }
  
  Pi = matrix(0.001,dim(Yhet)[2],length(idx_list))
  for(i in 1:length(idx_list)){
    Pi[idx_list[[i]],i]=0.99 
  }
  
  PiR20 = matrix(0.001,dim(Yhet)[2],length(sets_indices20))
  for(i in 1:length(sets_indices20)){
    PiR20[sets_indices20[[i]],i]=0.99 
  }
  termsR = names(sets_indices20)
  
  
  PiW15 = matrix(0.001,dim(Yhet)[2],length(sets_indicesWiki15))
  for(i in 1:length(sets_indicesWiki15)){
    PiW15[sets_indicesWiki15[[i]],i]=0.99 
  }
  termsW = names(sets_indicesWiki15)
  
  if(!is.null(dim(Known))){
    known_names = colnames(Known)}else{
    known_names = "nExpressed"
  }
  
  h5save(Yhet, sym_names, terms, termsR, termsW, Pi, PiR20,PiW15,Known,known_names,file = paste('../data/', filename,sep=""))
  H5close()
    
}
