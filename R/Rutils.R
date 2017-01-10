#requires c2.cp.reactome.v4.0.symbols.gmt file and h.all.v5.0.symbols.gmt.txt file
write_scLVM2 <- function(countsMmusFilter, filename, minGenes=15, Known=NULL, is_counts = FALSE, 
                         setsAdd=NULL, setsAdd_names=NULL, data_dir = '../data/', species = 'Hs'){
  require(GSEABase)
  require(limma)
  fname = paste0(data_dir,'h.all.v5.0.symbols.gmt.txt')
  c2_set <- getGmt(paste0(data_dir,"c2.cp.reactome.v4.0.symbols.gmt"))
  
  if(species=='Hs'){
  require(org.Hs.eg.db)  
  x <- org.Hs.egSYMBOL}else{
    require(org.Mm.eg.db)  
    x <- org.Mm.egSYMBOL    
  }
  
  xx <- as.list(x[mappedkeys(x)])
    
  gene_ids <- geneIds(c2_set)
  sets_indices <- ids2indices(gene_ids, toupper(rownames(countsMmusFilter)))
  
  sets_indices20 = sets_indices[unlist(lapply(sets_indices,length))>20]
  termsR = names(sets_indices20)
  
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
  if(!is.null(setsAdd)){  
    if(length(setsAdd)==1){
      idx_list = c(idx_list,ids2indices(setsAdd, rownames(countsMmusFilter)))
      terms = c(terms, setsAdd_names)    
    }else{
      for(add_i in seq(1,length(setsAdd))){
        idx_list = c(idx_list,ids2indices(setsAdd[[add_i]], rownames(countsMmusFilter)))
        terms = c(terms, setsAdd_names[add_i]) 
        
        sets_indices20 = c(sets_indices20,ids2indices(setsAdd[[add_i]], rownames(countsMmusFilter))) 
        termsR = c(termsR, setsAdd_names[add_i]) 
      }
    }
  }
  
  IMSigDB = matrix(0.0,dim(Yhet)[2],length(idx_list))
  for(i in 1:length(idx_list)){
    IMSigDB[idx_list[[i]],i]=1.0
  }
  
  IR20 = matrix(0.0,dim(Yhet)[2],length(sets_indices20))
  for(i in 1:length(sets_indices20)){
    IR20[sets_indices20[[i]],i]=1.0 
  }
    
  
  if(!is.null(dim(Known))){
    known_names = colnames(Known)}else{
      known_names = "nExpressed"
    }
  IMSigDB = t(IMSigDB)
  IR20 = t(IR20)
  h5save(Yhet, sym_names, terms, termsR, IMSigDB, IR20,Known,known_names,file = paste(data_dir, filename,sep=""))
  H5close()
  
}
