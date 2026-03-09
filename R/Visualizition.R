#' plot your sanky to visual drug-target-module
#' @importFrom magrittr %>%
#' @importFrom stringr str_split_fixed
#' @importFrom ggalluvial to_lodes_form
#' @importFrom ggplot2 ggplot scale_x_discrete geom_point geom_text geom_density geom_bar scale_fill_manual scale_colour_manual scale_size_continuous labs xlab ylab theme theme_bw element_text element_blank ggtitle guides ggsave
#' @importFrom ggraph geom_node_point
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggalluvial geom_flow geom_stratum
#' @param dfdrug dataframe, dataframe yielding by DrugEft containing drugs you want to plot
#' @param ctnets ctnet CTNet Class, output by ConstructCTnet
#' @param filename character, define your filename by your celltype
#' @param pattle character vector, default is NULL with the colors we provide. You can also give your colors pattle.
#' @param width output pdf width
#' @param height output pdf height
#' @return none
#' @export
plotsanky <- function(dfdrug,ctnets,filename='Celltype',pattle=NULL,width=6,height=10){
    mdsdf <- ctnets@ctmodule$modules
    df1 <- do.call(rbind,lapply(1:nrow(dfdrug),function(i){
        tmp <- dfdrug[i,'target']
        tmp <- stringr::str_split_fixed(string=tmp,pattern=';',n=Inf) %>% unlist(.) %>% as.vector(.)
        data.frame(Drug=dfdrug[i,'Drug'],target=tmp)
    }))
    df2 <- subset(mdsdf,node %in% df1$target)
    df <- merge(df1,df2,by.x='target',by.y='node')
    df$module <- paste0('Module',df$module)
    df <- df[,c('Drug','target','module')]
    df$Drug <- factor(df$Drug,levels=dfdrug$Drug)
    UCB_lodes <- ggalluvial::to_lodes_form(df,axes = 1:3,id = "Cohort")
    ypos <- 1/table(UCB_lodes$stratum)
    UCB_lodes$weight <- ypos[UCB_lodes$stratum]
    UCB_lodes$weight <- ifelse(UCB_lodes$x %in% c('module','target'),nrow(dfdrug)/nrow(df),UCB_lodes$weight)
    if(is.null(pattle)){
        pattle <- c('#85AEC9','#A1E5D7','#FACBE5','#FCF6B8','#FDB46F','#FF887E','#C2BEE3','#C7E873')
        colourcount <- length(unique(UCB_lodes$stratum))
        getPalette = colorRampPalette(pattle)
        sankycols <- getPalette(colourcount)
    } else{
        colourcount <- length(unique(UCB_lodes$stratum))
        getPalette = colorRampPalette(pattle)
        sankycols <- getPalette(colourcount)
    }
    p <- ggplot2::ggplot(UCB_lodes,ggplot2::aes(x = x, y=weight,stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
    ggplot2::scale_x_discrete(expand = c(0, 0.1)) + 
    ggalluvial::geom_flow(width = 1/3) + #线跟方块间空隙的宽窄
    ggalluvial::geom_stratum(alpha = 0.7,width = 0.85) + #方块的透明度、宽度
    ggplot2::geom_text(ggplot2::aes(label = stratum),stat = ggalluvial::StatStratum, size = 6,color="black") + #文字大小、颜色
    ggplot2::scale_fill_manual(values=sankycols)+
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_bw() + #去除背景色
    ggplot2::theme(panel.grid =element_blank()) + #去除网格线
    ggplot2::theme(panel.border = element_blank()) + #去除外层边框
    ggplot2::theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
    ggplot2::ggtitle("")+
    ggplot2::guides(fill = 'none') 
    ggplot2::ggsave(p,file=paste0(filename,'.Drugsanky.pdf'),width = width,height = height)
}

#' plot NetworkModule PPI
#' @importFrom magrittr %>%
#' @importFrom stringr str_split_fixed
#' @importFrom ggplot2 ggplot geom_point geom_text geom_density geom_bar scale_fill_manual scale_colour_manual scale_size_continuous labs xlab ylab theme theme_bw element_rect element_text element_blank ggtitle guides ggsave
#' @importFrom ggraph geom_node_point geom_edge_link create_layout ggraph
#' @importFrom ggrepel geom_label_repel
#' @importFrom igraph induced_subgraph layout_with_kk V   
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom ggrastr rasterise
#' @param ctnet CTNet Class, output by ConstructCTnet
#' @param filename character, define your filename by your celltype
#' @param SamType character, define your filename by your sample phenotype
#' @param labelgene character, default is NULL with the seed/target genes ploted. You can provide the label genes you want to plot
#' @param cols character, default is NULL with the colors we provide. You can provide colors as you like
#' @param width output pdf width
#' @param height output pdf height
#' @return none
#' @export
netVisual <- function(ctnet=NULL,filename=NULL,SamType=NULL,labelgene=NULL,cols=NULL,width=5,height=5){
    mdsdf <- ctnet@ctmodule$modules
    phyerMod <- ctnet@ctmodule$phyerMod
    phyerMod <- subset(phyerMod,pvalue<0.05)[,'module']
    #print module network
    for(md in phyerMod){
        submdsdf <- subset(mdsdf,module==md)
        print(paste0('Plot Module',md,' network'))
        #define visual hub by seed and target
        g <- igraph::induced_subgraph(ctnet@ctnets$graph, vids = submdsdf$node)
        netfr_layout <- igraph::layout_with_kk(g)
        colnames(netfr_layout)[c(1,2)] <- c('x','y')
        netfr_layout <- data.frame(netfr_layout,check.names=F)
        netfr_layout$name <- igraph::V(g)$name
        rownames(netfr_layout) <- netfr_layout$name
        stgene <- intersect(ctnet@features$seedgene,ctnet@features$target) %>% intersect(.,netfr_layout$name)
        #seed/target;seed;target;other
        seedonly <- intersect(ctnet@features$seedgene,netfr_layout$name) %>% setdiff(.,stgene)
        targetonly <- intersect(ctnet@features$target,netfr_layout$name) %>% setdiff(.,stgene)
        other <- netfr_layout$name[!netfr_layout$name %in% c(stgene,seedonly,targetonly)]
        if(is.null(labelgene)){
            mylabelgene <- stgene
        } else{
            mylabelgene <- labelgene
        }
        labeldf <- c(rep('seed/target',length(stgene)),rep('seedonly',length(seedonly)),
                            rep('targetonly',length(targetonly)),rep('other',length(other)))
        names(labeldf) <- c(stgene,seedonly,targetonly,other)
        netfr_layout$label <- labeldf[netfr_layout$name]
        #add degree to umap layout
        netfr_layout$size <- ifelse(netfr_layout$label=='seed/target',1,
                                    ifelse(netfr_layout$label %in% c('seedonly','targetonly'),0.5,0.1))
        graph <- g %>% tidygraph::as_tbl_graph(directed=FALSE) %>% 
                tidygraph::activate(nodes)
        igraph::V(graph)$anno <- ifelse(igraph::V(graph)$name %in% mylabelgene, igraph::V(graph)$name, "")
        netfr_layout <- netfr_layout[igraph::V(graph)$name,]
        #create layout for ggraph
        lay <- ggraph::create_layout(graph, netfr_layout)
        lay$anno <- igraph::V(graph)$anno
        nodecols <- c('#E25508','#DDF9DC','#7BCE9F','#56519B')
        names(nodecols) <- c('seed/target','other','seedonly','targetonly')
        p<- ggraph::ggraph(lay) + 
        ggrastr::rasterise(
        ggplot2::geom_point(inherit.aes=FALSE, data=netfr_layout, ggplot2::aes(x=x, y=y), alpha=0.9, size=0.1),dpi=300) +
        ggraph::geom_edge_link(alpha=0.5,width=0.5,colour='Bisque') +
        ggraph::geom_node_point(data=subset(lay, label != 'seed/target'), ggplot2::aes(size=size,fill=label,colour=label), shape=21) + 
        ggraph::geom_node_point(data=subset(lay, label == 'seed/target'), ggplot2::aes(size=size,fill=label,colour=label), shape=23) +
        ggrepel::geom_label_repel(data = lay, ggplot2::aes(x = x, y = y, label = anno), max.overlaps=Inf,
                                box.padding = 0.8, point.padding = NA,segment.color = "black")+
        ggplot2::scale_fill_manual(values=nodecols)+
        ggplot2::scale_colour_manual(values=nodecols)+
        ggplot2::scale_size_continuous(range = c(1,3))+
        ggplot2::theme(panel.background = element_rect(fill='white'),plot.title = element_text(hjust=0.5)) + 
        ggplot2::ggtitle(paste0('Predicted PPI in ',filename,' of Module',md))
        if(grepl('\\/',filename)){
            filename <- gsub('\\/','',filename)
        } else{
            filename <- filename
        }
        ggplot2::ggsave(p,file=paste0(filename,'.',SamType,'.Module',md,'.PPI.pdf'),width=width,height=height,bg = "transparent")   
    }

    #print all network,colored by module
    p <- ggplot2::ggplot(ctnet@ctnets$nets,ggplot2::aes(x = weight)) +
        ggplot2::geom_density(fill = "#FFEFFF", color='#7030A0',alpha = 0.5,linewidth = 1.5) +  # 绘制带填充色的密度曲线
        ggplot2::labs(title = "Weight", x = "Value", y = "Density") +
        ggplot2::theme_bw()
    ggplot2::ggsave(p,file=paste0(filename,'.',SamType,'.PPI.density.pdf'),width=3,height=3)

    modulenum <- table(mdsdf$module) %>% data.frame(.,check.names=F)
    colnames(modulenum) <- c('Module','Nodes')
    modulenum$Nodes <- log2(modulenum$Nodes)
    modulenum$Module <- paste0('Module',modulenum$Module)
    if(is.null(cols)){
        colourCount <- length(unique(mdsdf$module))
        if(length(unique(mdsdf$module))>=8){
            getPalette = colorRampPalette(c('#C7E873','#A1E5D7','#C2BEE3','#FF887E','#85AEC9','#FACBE5','#FDB46F','#FCF6B8'))
            cols <- getPalette(colourCount)
        }else{
            cols <- c('#C7E873','#A1E5D7','#C2BEE3','#FF887E','#85AEC9','#FACBE5','#FDB46F','#FCF6B8')
        }
    } else{
        cols <- cols
    }
    names(cols) <- paste0('Module',ctnet@ctmodule$phyerMod$module)
    p <- ggplot2::ggplot(modulenum, ggplot2::aes(x=Module, y=Nodes,fill=Module)) + 
    ggplot2::geom_bar(stat='identity',position = 'dodge',width=0.8)+theme_bw()+
    ggplot2::scale_fill_manual(values=cols)+
    ggplot2::theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5))+
    ggplot2::labs(y='Log2(Nodes) of Modules')
    ggplot2::ggsave(p,filename=paste0(filename,'.',SamType,'.PPI.ModuleNodes.pdf'),width=3,height=3)
}


#' Targetppi, not export
#' @importFrom magrittr %>%
#' @importFrom igraph distances induced_subgraph as_data_frame
#' @importFrom tidyr pivot_longer
#' @param module character, module
#' @param targets character vector, target genes
#' @param hop numeric, hop step
#' @param ctnet CTNet Class
Targetppi <- function(module=NULL,targets=NULL,hop=hop,ctnet=ctnets){
    dismat <- igraph::distances(ctnet@ctnets$graph,weight=NA)
    dismat <- data.frame(dismat,check.names=F)
    dismat$gene1 <- rownames(dismat)
    dismat <- dismat %>% tidyr::pivot_longer(cols=1:(ncol(dismat)-1),
                names_to = 'gene2',values_to = 'distances') %>% data.frame(.,check.names=F)
    edis <- subset(dismat,distances!=0)
    mdfinals <- do.call(rbind,lapply(targets,function(target){
        edisfilt <- subset(edis,gene1 == target & distances<=hop)
        mdsdf <- ctnet@ctmodule$modules
        mdnodes <- subset(mdsdf,module==module)[,'node']
        mdseed <- intersect(mdnodes,ctnet@features$seedgene)
        mdfinal <- intersect(mdseed,c(target,edisfilt$gene2))
        if(length(mdfinal)>0){
            mdfinal <- data.frame(Module=module,target=mdfinal)
        } else{
            mdfinal <- data.frame()
        }
        return(mdfinal)
    }))
    
    ppis <- do.call(rbind,lapply(targets,function(target){
        if(hop>1){
            addgene <- subset(edis,gene1 == target & distances<=(hop-1))[,'gene2']
        } else{
            addgene <- c()
        }
        mdfinals <- c(addgene,unique(mdfinals$target)) %>% unique(.)
        ppi <- igraph::induced_subgraph(ctnet@ctnets$graph, vids = unique(c(target,mdfinals)))
        ppi <- igraph::as_data_frame(ppi, what = "edges")
        return(ppi)
    }))
    reslist <- list(mdfinals,ppis)
    return(reslist)
}

#' output for cytoscape
#' @importFrom magrittr %>%
#' @param drugs character vector, give the drugs you want to output
#' @param proximity dataframe, the proximity file yielding by DrugEft function
#' @param ctnets ctnet CTNet Class, yielding by ConstructCTnet
#' @param hop numeric, default is 3. Maximum steps considered for computing network proximity.
#' @export
DB2Cyto <- function(drugs=NULL,proximity=NULL,ctnets=NULL,hop=3){
    myproximity <- subset(proximity,Molecular %in% drugs)
    df <- myproximity[,c('target','module','Molecular')]
    df <- unique(df)
    df$weight <- 1
    df$Type <- 'Drug-Target'
    
    mds <- unique(df$module)
    enrichgenes <- c()
    cytoppi <- data.frame()
    for(md in unique(df$module)){
        df1 <- subset(df,module==md)
        mdnum <- gsub('Module','',md)
        myreslist <- Targetppi(module=mdnum,targets=unique(df1$target),hop=hop,ctnet=ctnets)
        genetmp <- unique(myreslist[[1]]$target)
        enrichgenes <- c(enrichgenes,genetmp)
        ppitmp <- myreslist[[2]]
        cytoppi <- rbind(cytoppi,ppitmp)
    }
    finalgene <- unique(enrichgenes)

    finalppi <- unique(cytoppi)
    df2 <- unique(df[,c('target','Molecular','weight','Type')])
    colnames(df2)[1:2] <- c('from','to')
    finalppi$Type <- 'PPI'
    finalppi <- rbind(df2,finalppi)
    allnode <- c(finalppi$from,finalppi$to) %>% unique(.)
    ppiky <- data.frame(gene=allnode,keytype=ifelse(allnode %in% unique(df$target),'Targets',
                                                    ifelse(allnode %in% ctnets@features$seedgene,'Seed',
                                                        ifelse(allnode %in% drugs,'Drugs','Others'))))
    cytores <- list(finalgene,finalppi,ppiky)
    return(cytores)
}

#' HERB output for cytoscape
#' @importFrom magrittr %>%
#' @param drugs character vector, give the drugs you want to output
#' @param proximity dataframe, the proximity file yielding by DrugEft function
#' @param ctnets ctnet CTNet Class, yielding by ConstructCTnet
#' @param hop numeric, default is 3. Maximum steps considered for computing network proximity.
#' @export
HERB2Cyto <- function(drugs=NULL,proximity=NULL,ctnets=NULL,hop=3){
    myproximity <- subset(proximity,Molecular %in% drugs)
    df <- myproximity[,c('target','module','Molecular','RelatedDisease')]
    df <- unique(df)
    df$weight <- 1
    df$Type <- 'Drug-Target'

    mds <- unique(df$module)
    enrichgenes <- c()
    cytoppi <- data.frame()
    for(md in unique(df$module)){
        df1 <- subset(df,module==md)
        mdnum <- gsub('Module','',md)
        myreslist <- Targetppi(module=mdnum,targets=unique(df1$target),hop=hop,ctnet=ctnets)
        genetmp <- unique(myreslist[[1]]$target)
        enrichgenes <- c(enrichgenes,genetmp)
        ppitmp <- myreslist[[2]]
        cytoppi <- rbind(cytoppi,ppitmp)
    }
    finalgene <- unique(enrichgenes)

    finalppi <- unique(cytoppi)
    df2 <- unique(df[,c('target','Molecular','weight','Type')])
    colnames(df2)[1:2] <- c('from','to')
    finalppi$Type <- 'PPI'
    finalppi <- rbind(df2,finalppi)

    #herbclass
    mydieasecla <- herblib
    mydieasecla <- mydieasecla$RelatedDisease %>% str_split_fixed(.,';',n=Inf) %>% as.vector(.) %>% unique(.)
    mydieasecla <- data.frame(Disease=mydieasecla) %>% merge(.,DiseaseClass,by='Disease')

    #add target-disease to finalppi
    addppi <- do.call(rbind,lapply(drugs,function(d){
        addtmp <- subset(herblib,Molecular == d)[,c('Molecular','RelatedDisease')] %>% unique(.)
        addtmp <- addtmp$RelatedDisease %>% str_split_fixed(.,';',n=Inf) %>% as.vector(.) %>% unique(.)
        addtmp <- data.frame(from=d,to=addtmp,weight=1,Type='Drug-Diseases')
        return(addtmp)
    }))
    finalppi <- rbind(finalppi,addppi)
    #add diseaseclass to ppiky
    addky <- merge(data.frame(Disease=addppi$to),mydieasecla,by='Disease')
    allnode <- c(finalppi$from,finalppi$to) %>% unique(.)
    ppiky <- data.frame(gene=allnode,keytype=ifelse(allnode %in% unique(df$target),'Targets',
                                                    ifelse(allnode %in% ctnets@features$seedgene,'Seed',
                                                        ifelse(allnode %in% drugs,'Drugs',
                                                            ifelse(allnode %in% addky$Disease,'Disease','Other')))))
    cytores <- list(finalgene,finalppi,ppiky)
    return(cytores)
}

#' extract drugs
#' @importFrom magrittr %>%
#' @importFrom utils read.delim2
#' @param key character, the phenotype you want to grep
#' @param keytype character, grep key of column: Description or Details
#' @param celltype character
#' @param cutoff numeric
#' @export

keygrep <- function(key=NULL,keytype='Description',celltype=NULL,cutoff){
  if(grepl('\\/',celltype)){
    filefix <- gsub('\\/','',celltype)
  } else{
    filefix <- celltype
  }
  filename <- paste0(filefix,'.drug_rank_',cutoff,'.txt')
  drug_sta <- read.delim2(filename,sep='\t',check.names=F,header=T)
  drug_sta$DRScore <- as.numeric(drug_sta$DRScore)
  drug_sta$ranknum <- rank(drug_sta$DRScore, ties.method = "min")
  myindex <- grepl(key,drug_sta[,keytype],ignore.case=TRUE)
  drugsearch <- drug_sta[myindex,]
  if(nrow(drugsearch)>0){
    drugsearch$Celltype <- celltype
    drugsearch$cutoff <- cutoff
  } else{
    drugsearch <- data.frame()
  }
  return(drugsearch)
}