#' proximity process
#' @description
#' Specify your cutoff of proximity value computed by drugbank, not export
#' @importFrom stringr str_split
#' @param proximity dataframe
#' @param zcutoff numeric
#' @param ctnet CTNet
#' @param filename character
cutproximity_drugdb <- function(proximity=NULL,zcutoff=NULL,ctnet=NULL,filename=NULL){
        proximity <- subset(proximity,Zproximity<zcutoff)
        if(nrow(proximity)<1){
            stop('Error! Please increase your zcutoff to left more drugs!')
        }
        save(proximity,file=paste0(filename,'.proximity_',abs(zcutoff),'.rda'))
        write.table(proximity,file=paste0(filename,'.proximity_',abs(zcutoff),'.txt'),sep='\t',quote=F,row.names=F)
        
        #读取rwr
        myrwr <- ctnet@RWR$rwr
        rownames(myrwr) <- myrwr$node

        drug_sta <- aggregate(proximity$Zproximity,by=list(proximity$Molecular,proximity$module),mean)#其实就是直接求该模块对应靶点的平均
        colnames(drug_sta) <- c('Drug','Module','value')
        mdsta <- aggregate(drug_sta$value,by=list(drug_sta$Drug),length)
        colnames(mdsta) <- c('Drug','ModuleNumber')
        mdsta2 <- aggregate(drug_sta$value,by=list(drug_sta$Drug),mean)
        colnames(mdsta2) <- c('Drug','DRScore')
        mdsta <- merge(mdsta,mdsta2,by='Drug')
        drug_sta <- merge(drug_sta,mdsta,by='Drug')
        drug_sta <- drug_sta[,c('Drug','DRScore')]
        drug_sta <- unique(drug_sta)

        drug_sta2 <- do.call(rbind,lapply(unique(proximity$Molecular),function(drug){
            tmp <- subset(proximity,Molecular==drug)[,c('Molecular','target','direct')]
            tmp <- unique(tmp)
            drug_sta2tmp <- data.frame(Drug=drug,direct=paste0(round(tmp$direct,2),collapse=';'))
        }))
        drug_sta <- merge(drug_sta,drug_sta2,by='Drug')
        prodirect <- unique(proximity[,c('Molecular','Group','target','direct','Regulation','DrugBankID','KEGGDrugID','DrugClass','Description','Details')])
        drug_sta$target <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            target <- paste0(subset(prodirect,Molecular==drug)[,'target'],collapse=';')
            return(target)
        })  %>% unlist(.) %>% as.vector(.)

        drug_sta$Group <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            group <- unique(prodirect[prodirect$Molecular== drug,'Group'])
            return(group)
        })  %>% unlist(.) %>% as.vector(.) 

        drug_sta$RWR <- lapply(1:nrow(drug_sta),function(i){
            target <- drug_sta[i,'target']
            irwr <- myrwr[c(stringr::str_split(target,';') %>% unlist(.)),'rwr'] %>% round(.,2)
            irwr <- paste0(irwr,collapse=';')
            return(irwr)
        })  %>% unlist(.) %>% as.vector(.)      

        drug_sta$direct <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            direct <- paste0(subset(prodirect,Molecular==drug)[,'direct'] %>% round(.,digits= 2),collapse=';')
            return(direct)
        })  %>% unlist(.) %>% as.vector(.)

        drug_sta$Regulation <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            target <- drug_sta[i,'target']
            reg <- prodirect[prodirect$Molecular== drug & prodirect$target %in% c(stringr::str_split(target,';') %>% unlist(.)),'Regulation']
            reg <- paste0(reg,collapse=';')
            return(reg)
        })  %>% unlist(.) %>% as.vector(.)

        drug_sta$DrugBankID <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            drugid <- unique(prodirect[prodirect$Molecular== drug,'DrugBankID'])
            return(drugid)
        })  %>% unlist(.) %>% as.vector(.)        

        drug_sta$Description <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            info <- unique(prodirect[prodirect$Molecular== drug,'Description'])
            return(info)
        })  %>% unlist(.) %>% as.vector(.)

        drug_sta$KEGGDrugID <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            kegg <- unique(prodirect[prodirect$Molecular== drug,'KEGGDrugID'])
            return(kegg)
        })  %>% unlist(.) %>% as.vector(.)    

        drug_sta$class <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            class <- unique(prodirect[prodirect$Molecular== drug,'DrugClass'])
            return(class)
        })  %>% unlist(.) %>% as.vector(.)

        drug_sta$Details <- lapply(1:nrow(drug_sta),function(i){
            drug <- drug_sta[i,'Drug']
            Details <- unique(prodirect[prodirect$Molecular== drug,'Details'])
            return(Details)
        })  %>% unlist(.) %>% as.vector(.)
        return(drug_sta)
}

#' proximity process
#' @description
#' Specify your cutoff of proximity value computed by herb, not export
#' @importFrom stringr str_split
#' @param proximity dataframe
#' @param zcutoff numeric
#' @param ctnet CTNet
#' @param filename character
cutproximity_herb <- function(proximity=NULL,zcutoff=NULL,ctnet=NULL,filename=NULL){
    proximity <- subset(proximity,Zproximity<zcutoff)
    if(nrow(proximity)<1){
        stop('Error! Please increase your zcutoff to left more drugs!')
    }
    save(proximity,file=paste0(filename,'.proximity_',abs(zcutoff),'.rda'))
    write.table(proximity,file=paste0(filename,'.proximity_',abs(zcutoff),'.txt'),sep='\t',quote=F,row.names=F)
    
    myrwr <- ctnet@RWR$rwr
    rownames(myrwr) <- myrwr$node
    
    drug_sta <- aggregate(proximity$Zproximity,by=list(proximity$Molecular,proximity$module),mean)
    colnames(drug_sta) <- c('Drug','Module','value')
    mdsta <- aggregate(drug_sta$value,by=list(drug_sta$Drug),length)
    colnames(mdsta) <- c('Drug','ModuleNumber')
    mdsta2 <- aggregate(drug_sta$value,by=list(drug_sta$Drug),mean)
    colnames(mdsta2) <- c('Drug','DRScore')
    mdsta <- merge(mdsta,mdsta2,by='Drug')
    drug_sta <- merge(drug_sta,mdsta,by='Drug')
    drug_sta <- drug_sta[,c('Drug','DRScore')]
    drug_sta <- unique(drug_sta)
    drug_sta2 <- do.call(rbind,lapply(unique(proximity$Molecular),function(drug){
        tmp <- subset(proximity,Molecular==drug)[,c('Molecular','target','direct')]
        tmp <- unique(tmp)
        drug_sta2tmp <- data.frame(Drug=drug,direct=paste0(round(tmp$direct,2),collapse=';'))
    }))
    drug_sta <- merge(drug_sta,drug_sta2,by='Drug')
    prodirect <- unique(proximity[,c('Molecular','target','direct','RelatedDisease')])
    drug_sta$target <- lapply(1:nrow(drug_sta),function(i){
        drug <- drug_sta[i,'Drug']
        target <- paste0(subset(prodirect,Molecular==drug)[,'target'],collapse=';')
        return(target)
    })  %>% unlist(.) %>% as.vector(.)

    drug_sta$RWR <- lapply(1:nrow(drug_sta),function(i){
        target <- drug_sta[i,'target']
        irwr <- myrwr[c(stringr::str_split(target,';') %>% unlist(.)),'rwr'] %>% round(.,2)
        irwr <- paste0(irwr,collapse=';')
        return(irwr)
    })  %>% unlist(.) %>% as.vector(.)   
        
    drug_sta$direct <- lapply(1:nrow(drug_sta),function(i){
        drug <- drug_sta[i,'Drug']
        direct <- paste0(subset(prodirect,Molecular==drug)[,'direct'] %>% round(.,digits= 2),collapse=';')
        return(direct)
    })  %>% unlist(.) %>% as.vector(.)

    drug_sta$RelatedDisease <- lapply(1:nrow(drug_sta),function(i){
        drug <- drug_sta[i,'Drug']
        target <- drug_sta[i,'target']
        reg <- prodirect[prodirect$Molecular== drug & prodirect$target %in% c(stringr::str_split(target,';') %>% unlist(.)),'RelatedDisease']
        reg <- stringr::str_split(reg,pattern=';') %>% unlist(.) %>% unique(.)
        reg <- paste0(reg,collapse=';')
        return(reg)
    })  %>% unlist(.) %>% as.vector(.)
    return(drug_sta)
}

#' Network Proximity computation
#' @importFrom igraph distances
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom stringr str_split_fixed str_split
#' @importFrom ggalluvial to_lodes_form geom_flow geom_stratum
#' @importFrom ggplot2 ggplot scale_x_discrete geom_text xlab ylab theme theme_bw element_blank ggtitle guides ggsave
#' @importFrom pbapply pbapply
#' @param species character,define the species, should be 'mmu' or 'human'
#' @param s numeric, default is 0.2. Shrinkage coefficient, ranging from 0 to 1. The smaller the value, the smaller the difference in weights between different scoring levels.
#' @param log2fc numeric, default is 0.1. Genes whose |log2fc| greater than this threshold are considered as differentially expressed genes
#' @param lib character, default is NULL. The lib must be NULL (Drugbank database) or herb (HERB database)
#' @param ctnet CTNet Class, output by ConstructCTnet
#' @param celltype character, define celltype, must be consistent with the celltype you provided to mysigDEG
#' @param prox logical, whether computing proximity, once you have run this function, you will get the raw proximity value without cutoff (named yourcelltype.proximity_raw.rda). If you want to change the cutoff but do not compute proximity value again for reducing computational time, you can just define prox=FALSE and give the proximity_raw.rda filepath to myprox.
#' @param zcutoff numeric, default is 0.5. Cutoff of proximity value to compute Drug repurposing score.
#' @param myprox character, yourcelltype.proximity_raw.rda filepath
#' @param mysigDEG dataframe, a dataframe given your differentially expressed genes which must has column: gene,NewType,avg_log2FC in which NewType column defining your celltypes
#' @param hop numeric, default is 3. Maximum steps considered for computing network proximity.
#' @param sankyn numeric, drug number for sankyplot
#' @param sankycols character vector, default is NULL for using the colors we provide. You can provide your color codes to sankycols
#' @return data.frame
#' @export
DrugEft <- function(species='mmu',s=0.2,log2fc=0.1,lib=NULL,ctnet=NULL,celltype=NULL,prox=TRUE,myprox=NULL,
                    mysigDEG=mysigDEG,hop=3,zcutoff=-0.5,sankyn=20,sankycols=NULL){
    #Compute DrugEffect
    if(grepl('\\/',celltype)){
        filename <- gsub('\\/','',celltype)
    } else{
        filename <- celltype
    }
    mdsdf <- ctnet@ctmodule$modules
    if(prox){
        if(is.null(lib)){
            druglib <- druglib
            if(species=='mmu'){
                druglib$gene_symbol <- str_to_title(druglib$gene_symbol)
            }
        }else{
            druglib <- herblib
            if(species=='mmu'){
                druglib$gene_symbol <- str_to_title(druglib$gene_symbol)
            }
        }
        if(species=='mmu'){
            NDGs <- NDDgene$gene_symbol %>% str_to_title(.) %>% unique(.)
        }
        if(species=='human'){
            NDGs <- NDDgene$gene_symbol %>% unique(.)
        }

        print('Compute Shortest Distance')
        dismat <- igraph::distances(ctnet@ctnets$graph,weight=NA)
        dismat <- data.frame(dismat,check.names=F)
        dismat$gene1 <- rownames(dismat)
        dismat <- dismat %>% tidyr::pivot_longer(cols=1:(ncol(dismat)-1),
                    names_to = 'gene2',values_to = 'distances') %>% data.frame(.,check.names=F)
        edis <- subset(dismat,distances!=0)
        edisfilt <- subset(edis,distances<=hop)
    
        #comput drug effect
        edis_seed <- subset(edisfilt,gene1 %in% ctnet@features$target) %>% subset(.,gene2 %in% ctnet@features$seedgene)
        sigmodule <- subset(ctnet@ctmodule$phyerMod,pvalue<0.05)[,'module']
        
        if(length(sigmodule)>0){
            proximity2 <- do.call(rbind,lapply(sigmodule,function(md){
                print(paste0('Comput Proximity between all targets and Module',md))
                mdnodes <- subset(mdsdf,module==md)[,'node']
                mdseed <- intersect(mdnodes,ctnet@features$seedgene)
                if(length(mdseed)>0){
                    proximity2tmp <- do.call(rbind,pbapply::pblapply(ctnet@features$target,function(tt){
                        editmp2 <- subset(edisfilt,gene1 == tt) %>% subset(.,gene2 %in% mdseed)
                        editmp2 <- suppressWarnings(subset(editmp2,distances==min(editmp2$distances)))
                        
                        if(nrow(editmp2)>1){
                            editmp2 <- editmp2[1,]
                        } else{
                            editmp2 <- editmp2
                        }
                        return(editmp2)        
                    }))
                    if(nrow(proximity2tmp)>0){
                        values2 <- sum(proximity2tmp$distances)
                    } else{
                        values2 <- NA
                    }
                    proximity2tmp <- data.frame(module=paste0('Module',md),target='all',key='Proximity2',value=values2)
                } else{
                    proximity2tmp <- data.frame(module=paste0('Module',md),target='all',key='Proximity2',value=NA)
                }
                return(proximity2tmp)
            }))
        } else{
            stop('There is no module enriched by seed genes, please adjust valuek or valuen according to your network in ConstructCTnet and try again')
        }
        proximity2 <- na.omit(proximity2)

        #合并第一部分+第二部分####
        if(is.null(celltype)){
            ctDEG <- mysigDEG
            celltype <- 'RNA'
        } else{
            ctDEG <- subset(mysigDEG,NewType==celltype & abs(avg_log2FC)>log2fc)
            cellype <- celltype
        }
        proximity <- do.call(rbind,pbapply::pblapply(ctnet@features$target,function(t){
            proximitymd <- data.frame()
            usemd <- gsub('Module',"",unique(proximity2$module))
            for(md in usemd){
                #Comput Proximity between one target and Module
                mdnodes <- subset(mdsdf,module==md)[,'node']
                mdseed <- intersect(mdnodes,ctnet@features$seedgene)
                edistmp <- subset(edis_seed,gene1 == t) %>% subset(.,gene2 %in% mdseed)
                if(nrow(edistmp)>0){
                    values <- sum(edistmp$distances)
                    value2 <- subset(proximity2,module==paste0('Module',md))[,'value']
                    value <- (values+value2)/(length(mdseed)+length(ctnet@features$target))
                    w1 <- ifelse(t %in% ctnet@features$seedgene,1-s*1,
                                        ifelse(t %in% NDGs,1-s*0.5,1-s*0))
                    proximitytmp <- data.frame(module=paste0('Module',md),target=t,
                                                proximity1=values,proximity2=value2,
                                                proximity3=value,w1=w1)
                } else {
                    proximitytmp <- data.frame() 
                }

                if(t %in% ctDEG$gene){
                    values <- subset(ctDEG,gene==t)[,'avg_log2FC'] * -log10(subset(ctDEG,gene==t)[,'p_val_adj'])
                    direct <- data.frame(module=paste0('Module',md),target=t,key='Direct',value=values)
                } else{
                    direct <- data.frame(module=paste0('Module',md),target=t,key='Direct',value=NA)
                }

                if(nrow(proximitytmp)>0){
                    proximitym <- cbind(proximitytmp,direct=direct$value)
                    proximitym$NewType <- celltype
                } else {
                    proximitym <- data.frame()
                }
                proximitymd <- rbind(proximitymd,proximitym)
            }
            return(proximitymd)
        }))
        if(nrow(proximity)<1){
            stop('Error！No seed genes cause of too less DEGs!')
        }

        #Drugbank drug
        if(is.null(lib)){
            druglib <- druglib
            if(species=='mmu'){
                druglib$gene_symbol <- str_to_title(druglib$gene_symbol)
            }
            proximity$w2 <- proximity$w1/mean(proximity$w1)
            proximity$proximity <- proximity$proximity3 * proximity$w2
            print('Compute Drug Effect')
            proximity <- merge(proximity,druglib[,c('gene_symbol','Molecular','Group','Regulation','DrugBankID','KEGGDrugID','DrugClass','Description','Details')],by.x='target',by.y='gene_symbol')
            proximity <- unique(proximity)
            proximity$Zproximity <- (proximity$proximity-mean(proximity$proximity))/sd(proximity$proximity)
            if(grepl('\\/',celltype)){
                filename <- gsub('\\/','',celltype)
            } else{
                filename <- celltype
            }
            save(proximity,file=paste0(filename,'.proximity_raw.rda'))
            drug_sta <- cutproximity_drugdb(proximity=proximity,zcutoff=zcutoff,ctnet=ctnet,filename=filename)

        } else{
            if(lib =='herb'){
                druglib <- herblib
            } else{
                stop('Error! The lib must be NULL or herb.')
            }        
            if(species=='mmu'){
                druglib$gene_symbol <- str_to_title(druglib$gene_symbol)
            }
            proximity$w2 <- proximity$w1/mean(proximity$w1)
            proximity$proximity <- proximity$proximity3 * proximity$w2
            print('Compute Drug Effect')
            proximity <- merge(proximity,druglib[,c('gene_symbol','Molecular','RelatedDisease')],by.x='target',by.y='gene_symbol')
            proximity <- unique(proximity)
            if(grepl('\\/',celltype)){
                filename <- gsub('\\/','',celltype)
            } else{
                filename <- celltype
            }
            proximity$Zproximity <- (proximity$proximity-mean(proximity$proximity))/sd(proximity$proximity)
            save(proximity,file=paste0(filename,'.proximity_raw.rda'))
            drug_sta <- cutproximity_herb(proximity=proximity,zcutoff=zcutoff,ctnet=ctnet,filename=filename)
        }
    } else{
        if(is.null(myprox)){
            stop('Error, if prox=FALSE, you must assign your filename of proximity_raw.rda to myprox')
        }
        if(is.null(lib)){
            load(myprox)
            drug_sta <- cutproximity_drugdb(proximity=proximity,zcutoff=zcutoff,ctnet=ctnet,filename=filename)
        } else{
            load(myprox)
            drug_sta <- cutproximity_herb(proximity=proximity,zcutoff=zcutoff,ctnet=ctnet,filename=filename)
        }
    }

    #排序
    drug_sta <- drug_sta[order(drug_sta$DRScore,decreasing=F),]
    dfdrug <- drug_sta[1:sankyn,]
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
    if(is.null(sankycols)){
        pattle <- c('#85AEC9','#FACBE5','#FDB46F','#FF887E','#C2BEE3','#A1E5D7','#FCF6B8','#C7E873')
        colourcount <- length(unique(UCB_lodes$stratum))
        getPalette = colorRampPalette(pattle)
        sankycols <- getPalette(colourcount)
    } else{
        colourcount <- length(unique(UCB_lodes$stratum))
        getPalette = colorRampPalette(sankycols)
        sankycols <- getPalette(colourcount)
    }
    p <- ggplot2::ggplot(UCB_lodes,ggplot2::aes(x = x, y=weight,stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) + 
    ggplot2::scale_x_discrete(expand = c(0, 0.1)) + 
    ggalluvial::geom_flow(width = 1/3) + 
    ggalluvial::geom_stratum(alpha = 0.7,width = 0.85) + 
    ggplot2::geom_text(ggplot2::aes(label = stratum),stat = ggalluvial::StatStratum, size = 6,color="black") + 
    ggplot2::scale_fill_manual(values=sankycols)+
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid =element_blank()) + 
    ggplot2::theme(panel.border = element_blank()) + 
    ggplot2::theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + 
    ggplot2::ggtitle("")+
    ggplot2::guides(fill = 'none') 
    ggplot2::ggsave(p,file=paste0(filename,'.Drugsanky_',abs(zcutoff),'.pdf'),width = 6,height = 10)   
    return(drug_sta)
}