# scPsyDrug
For optimal performance, run this software on a high-memory server.

Installation
1.	Install dependencies (In Terminal):
```
conda install -c conda-forge r-tidyr r-stringr r-ggplot2 r-pbapply r-igraph r-tidygraph r-ggraph r-ggrastr r-ggalluvial r-ggrepel r-magrittr r-matrix
conda install -c conda-forge hdf5 suitesparse openssl libcurl libxml2
```

3.	Install ACTIONet & SCINET (in R)
```
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
devtools::install_github("shmohammadi86/SCINET")

4. Install scPsyDrug
install.packages("scPsyDrug_0.0.1.tar.gz", repos = NULL, type = "source") 
```
More details see https://github.com/cure-lab/SCINet

Using the MDD dataset as an example
```
#load data
load('demo.rda')

#run pepline for every celltype
mysigDEG <- DEGfilt
log2fc <- 0.1
SamType= 'Suicide'
myexp <- demo@assays$RNA$data
celltypes <- unique(NewType$NewType)
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
  drugEffect <- DrugEft(species='human',ctnet=ctnets,mysigDEG=mysigDEG,celltype=celltype,zcutoff=-1.5,prox=FALSE,myprox=paste0(filename,'.proximity_raw.rda'))
  write.table(drugEffect,file=paste0(filename,'.drug_rank_1.5.txt'),sep='\t',quote=F,row.names=F)
}
```
When the parameter imput = TRUE, the ConstructCTnet function may take a long time to run and will generate a graphdata.rda file. Therefore, if you have already run ConstructCTnet once and obtained graphdata.rda, but want to modify other parameters (e.g., log2fc) and run ConstructCTnet again, you can provide the absolute path to your existing graphdata.rda file and set imput = FALSE. This can significantly reduce the computation time.
```
graphdata=paste0(filename,'.graphdata.rda')
ctnets <- ConstructCTnet(species='human',myexp=myexp,ctcell=ctcell,log2fc=log2fc,celltype=celltype,imput=FALSE, graphdata=graphdata,intersect=0.01,mysigDEG=mysigDEG)
```
