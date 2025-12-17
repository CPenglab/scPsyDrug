# scPsyDrug
we developed a cell-type specific network based computational pipeline scPsyDrug to prioritize the drugs for psychiatric diseases by integrating the single-cell transcriptomics, protein-protein interactions, drug-target interactions and psychiatric risk genes. Single cell transcriptomics were integrated with reference protein-protein interactions to construct cell-type specific networks. The cell-type specific differentially expressed genes (DEGs) and psychiatric disease risk (PSD) genes were used to define seed genes in cell-type specific network. Network proximity between drug targets and the seed genes were calculated in each cell-type specific module, which was used to rank the drugs.

**Installation**
1.	Install dependencies (In Terminal):
```
conda install -c conda-forge r-tidyr r-stringr r-ggplot2 r-pbapply r-igraph r-tidygraph r-ggraph r-ggrastr r-ggalluvial r-ggrepel r-magrittr r-matrix
conda install -c conda-forge hdf5 suitesparse openssl libcurl libxml2
conda install r::r-devtools
```

2.	Install ACTIONet & SCINET (in R)
```
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
devtools::install_github("shmohammadi86/SCINET") # More details see https://github.com/shmohammadi86/SCINET
```

3. Install scPsyDrug

R Version >=4.3.3
```
install.packages("scPsyDrug_0.0.0.9000.tar.gz", repos = NULL, type = "source") 
```

**Quick Start**

scPsyDrug performs drug repurposing based on single-cell transcriptomics. The required input data are as follows:

(1) A list of differentially expressed genes (DEGs) between the disease (abnormal) condition and the normal control. It is recommended to name the columns using the following rules:

NewType: cell type

avg_log2FC: log2 fold change (log2FC)

gene: gene symbol

(2) A gene expression matrix for the disease (abnormal) condition, formatted as gene × cell.

(3) Metadata for the single-cell dataset, including cell names and cell types.

For optimal performance, run this software on a high-memory server.

Using the demo as an example.
```
#load data
load('demo.rda')

#run pepline for every celltype
mysigDEG <- DEGfilt
log2fc <- 0.1
SamType= 'Suicide'
myexp <- myexp[,subset(NewType,Type== SamType)[,'Cell']]
for(celltype in c('Astros','Excitatory_Neurons_L4')){
  print(celltype)
  if(grepl('\\/',celltype)){
    filename <- gsub('\\/','',celltype)
  } else{
    filename <- celltype
  }

  #Extract cells from one celltype
  ctcell <- subset(NewType,NewType==celltype & Type==SamType)[,'cell']

  #Network Construction
  ctnets <- ConstructCTnet(species='human',myexp=myexp,ctcell=ctcell,log2fc=log2fc,celltype=celltype,imput=TRUE,mysigDEG=mysigDEG)
  save(ctnets,file=paste0(filename,'nets.rda'))

  #Network visulization
  netVisual(ctnet=ctnets,filename=filename)

  #Drug Repurposing Score computation
  #cutoff=1.5
  drugEffect <- DrugEft(species='human',ctnet=ctnets,mysigDEG=mysigDEG,celltype=celltype,zcutoff=-1.5)
  write.table(drugEffect,file=paste0(filename,'.drug_rank_1.5.txt'),sep='\t',quote=F,row.names=F)
}
```
When the parameter imput = TRUE, the ConstructCTnet function may take a long time to run and will generate a graphdata.rda file. Therefore, if you have already run ConstructCTnet once and obtained graphdata.rda, but want to modify other parameters (e.g., log2fc) and run ConstructCTnet again, you can provide the absolute path to your existing graphdata.rda file and set imput = FALSE. This can significantly reduce the computation time.
```
graphdata=paste0(filename,'.graphdata.rda')
ctnets <- ConstructCTnet(species='human',myexp=myexp,ctcell=ctcell,log2fc=log2fc,celltype=celltype,imput=FALSE, graphdata=graphdata,mysigDEG=mysigDEG)
```

For the DrugEft function, the default setting is prox = TRUE, which generates a proximity_raw.rda file. This means that if you want to rerun DrugEft using a different cutoff value, you can set prox = FALSE and provide the absolute path to your proximity_raw.rda file via myprox. For example, myprox = "Astros.proximity_raw.rda". This can speed up the execution of DrugEft.
```
for(celltype in c('Astros','Excitatory_Neurons_L4'){
  print(celltype)
  if(grepl('\\/',celltype)){
    filename <- gsub('\\/','',celltype)
  } else{
    filename <- celltype
  }

  load(paste0(filename,'nets.rda'))
  drugEffect <- DrugEft(species='human',ctnet=ctnets,mysigDEG=mysigDEG,celltype=celltype,zcutoff=-1,prox=FALSE,myprox=paste0(filename,'.proximity_raw.rda'))
  write.table(drugEffect,file=paste0(filename,'.drug_rank_1.txt'),sep='\t',quote=F,row.names=F)

}
```
