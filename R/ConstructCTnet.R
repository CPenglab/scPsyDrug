#' CTNet Class
#' @description
#' S4 class for storing cell-type network.
#' @slot features A list of features.
#' @slot ctnets A list of cell-type-specific networks.
#' @slot ctmodule A list of modules.
#' @slot RWR A list of random walk results.
#' @exportClass CTNet
setClass(
  "CTNet",
  slots = list(
    features = "list",
    ctnets = "list",
    ctmodule = "list",
    RWR = "list"
  )
)

#' ConstructCTnet for construct cell-specific PPI
#' @importFrom stringr str_to_title
#' @importFrom pbapply pbapply
#' @importFrom igraph graph_from_data_frame as_edgelist as_adjacency_matrix V cluster_louvain membership
#' @importFrom ACTIONet decomp.ACTIONMR
#' @importFrom SCINET compute_gene_activities pair.datasets construct_cell_networks
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot scale_fill_manual theme theme_bw element_text element_blank geom_bar geom_line geom_vline geom_hline annotate labs ggsave 
#' @importFrom magrittr %>%
#' @import Matrix
#' @param species character,define the species, should be 'mmu' or 'human' 
#' @param min.pct numeric, default is 0.1 (10%). Percentage of gene expression: genes with an expression percentage below this value will be filtered out
#' @param intersect numeric, default is 0.05 (5%). Percentage of gene co-expression: gene pairs with a co-expression percentage below this value will be filtered out 
#' @param lib character, default is NULL. The lib must be NULL (Drugbank database) or herb (HERB database)
#' @param myexp gene expression matrix (dgMatrix) extracted from SeuratObject
#' @param ctcell character vector, cells within one cell-type
#' @param log2fc numeric, default is 0.1. Genes whose |log2fc| greater than this threshold are considered as differentially expressed genes 
#' @param imput logical, whether compute imputation. If imputation has been completed and you want to change log2fc or other parameters and rerun this function, define the imput=FALSE and give your graphdata.rda filepath to graphdata.
#' @param graphdata character vector, your graphdata.rda filepath. If imput=FALSE, you must give your graphdata.rda filepath to graphdata
#' @param seedgene character vector, if seedgene is NULL, we use differentially expressed psychiatry diseases genes as seedgene, you can also provide a geneset as seedgene 
#' @param mysigDEG dataframe, a dataframe given your differentially expressed genes which must has column: gene,NewType,avg_log2FC in which NewType column defining your celltypes
#' @param celltype character, define celltype, must be consistent with the celltype you provided to mysigDEG
#' @param k numeric, default is 15.  Number of nearest neighbors used in the ACTIONet decomposition. for improving computational speed we use one k without performing a search over a range of k. the larger the value of k, the slower the computation
#' @param thread_no numeric, default is 2. Number of CPU threads to use for parallel computation
#' @param maxiter numeric, default is 1000. Maximum number of iterations allowed for the random walk with restart (RWR) algorithm. The iteration stops earlier if convergence is reached, as determined by the epsilon threshold.
#' @param reso1 numeric, default is 0.6. Resolution of louvain algorithm identifying modules 
#' @param valuek numeric, default is 50. The expected number of genes drawn as seed genes each time
#' @param valuen numeric, default is 100. Total number of genes drawn each time
#' @return CTNet Class
#' @export
ConstructCTnet <- function(species='mmu',min.pct=0.1,intersect=0.05,lib=NULL,myexp=NULL,ctcell=NULL,log2fc=0.1,
                    imput=TRUE,graphdata=NULL,seedgene=NULL,mysigDEG=NULL,celltype=NULL,k=15,thread_no=2,maxiter=1000,reso1 = 0.6,valuek=50,valuen=100){
    
    ######SCINET to contruct cell-specific PPI###########
    #for high speed, the input expression matrix of decomp.ACTIONMR should be prexp
    ## Imputation
    if(imput){
        print('construct cell-specific network yourself')
        prexp <- myexp[,ctcell]
        rm(myexp)
        prexp <- prexp[rowSums(prexp>0) >= length(ctcell)*min.pct,]
        netgene <- rownames(prexp)
        print(paste0('select ',length(netgene),' network genes'))

        ## Load 2019Cell System cell-specific network as the reference interactome 
        if(species=='mmu'){
            csnet$geneA <- stringr::str_to_title(csnet$geneA)
            csnet$geneB <- stringr::str_to_title(csnet$geneB)
        }
        csnet <- subset(csnet,geneA %in% netgene) %>% subset(.,geneB %in% netgene)

        print('test gene coexpression')
        thresholds <- ncol(prexp)*intersect
        netcells <- colnames(prexp)
        csnet <- do.call(rbind,pbapply::pblapply(1:nrow(csnet),function(i){
            tmpnet <- csnet[i,]
            netexp <- prexp[c(tmpnet$geneA,tmpnet$geneB),]
            if(length(netcells[netexp[1,]>0 & netexp[2,]>0])>=thresholds){
                return(tmpnet)
            }
        }))
        prexp <- prexp[unique(c(csnet$geneA,csnet$geneB)),]
        netgene <- rownames(prexp)
        print(paste0('select ',length(netgene),' network genes'))
        G = igraph::graph_from_data_frame(csnet[,c('geneA','geneB')],directed=FALSE)
        print('imputation')
        imputout <- ACTIONet::decomp.ACTIONMR(as.matrix(prexp),k,k,thread_no=thread_no)
        if(grepl('\\/',celltype)){
            filename <- gsub('\\/','',celltype)
        } else{
            filename <- celltype
        }
        
        A <- imputout$W %*% imputout$H
        ## Estimate gene activity score within sub-population
        # It can be any imputed profile
        activity.scores = SCINET::compute_gene_activities(A = A, samples = 1:ncol(prexp), thread_no = thread_no)
        rownames(activity.scores) <- rownames(A)
        paired.datasets = SCINET::pair.datasets(G, activity.scores)
        EL = igraph::as_edgelist(G, names = FALSE)
        G.adj = as(igraph::as_adjacency_matrix(paired.datasets$net), 'TsparseMatrix')
        edge.idx = (G.adj@i)*nrow(G.adj) + (G.adj@j+1)
        ## Compute topological-specificity of genes in each network
        cellstate.nets = SCINET::construct_cell_networks(net = G.adj, gene_activities = paired.datasets$activity.scores, thread_no = thread_no)
        print('compute mean matrix')
        rm(prexp,csnet,A,activity.scores,paired.datasets,EL,edge.idx)
        cellstate.nets.list <- lapply(1:nrow(cellstate.nets),function(i){
            inet <- cellstate.nets[i,1][[1]]
            rownames(inet) <- rownames(G.adj)
            colnames(inet) <- rownames(G.adj)
            return(inet)
        })
        sum_cellstate.nets <- Reduce(`+`, cellstate.nets.list)
        mean_cellstate.nets <- sum_cellstate.nets / length(cellstate.nets.list)
        aggmatrix <- mean_cellstate.nets
        head(aggmatrix)
        aggmatrix[!upper.tri(aggmatrix,diag=FALSE)] <- NA
        aggnet <- data.frame(as.matrix(aggmatrix),check.names=F)
        aggnet$geneA <- rownames(aggnet)
        aggnet <- aggnet %>% tidyr::pivot_longer(cols=1:(ncol(aggnet)-1),
            names_to = 'geneB',values_to = 'weight') %>% data.frame(.,check.names=F)
        aggnet <- subset(aggnet,weight>0)
        head(aggnet)
        #排序
        aggnet <- aggnet[order(aggnet$geneA),]
        aggnet <- aggnet[order(aggnet$geneB),]
        mygraph<- igraph::graph_from_data_frame(aggnet[,c('geneA','geneB','weight')],directed=FALSE)
        save(aggnet,mygraph,file=paste0(filename,'.graphdata.rda'))

    } else{
        print('imput=FALSE, use your imputation output')
        if(is.null(graphdata)){
            stop('Error! You have to use your graphdata when imput=FALSE, please give your filename of graphdata.rda to graphdata')
        }
        load(graphdata)
    }

    #herb drug
    if(is.null(lib)){
        druglib <- druglib
    }  else if(lib=='herb'){
        druglib <- herblib
    } else{
        stop('Error! The lib must be NULL or herb.')
    }
    
    if(species=='mmu'){
        NDGs <- NDDgene$gene_symbol %>% stringr::str_to_title(.) %>% unique(.)
        druglib$gene_symbol <- stringr::str_to_title(druglib$gene_symbol)
        Targets <- unique(druglib$gene_symbol)
    }
    if(species=='human'){
        NDGs <- NDDgene$gene_symbol %>% unique(.)
        Targets <- unique(druglib$gene_symbol)
    }
    if(! 'gene' %in% colnames(mysigDEG)|!'NewType' %in% colnames(mysigDEG)|!'avg_log2FC' %in% colnames(mysigDEG)){
        stop('Error! mysigDEG must have coloumns: gene, NewType and avg_log2FC!')
    }
    if(is.null(seedgene)){
        seedgene <- unique(intersect(subset(mysigDEG,abs(avg_log2FC)>log2fc & NewType==celltype)[,'gene'],NDGs))
        seedgene <- seedgene[seedgene %in% names(igraph::V(mygraph))]
    } else{
        seedgene <- seedgene
        seedgene <- seedgene[seedgene %in% names(igraph::V(mygraph))]
    }
    targets <- Targets[Targets %in% names(igraph::V(mygraph))]
    print(paste0('network include ',length(seedgene),' seed genes'))
    print(paste0('network include ',length(targets),' target genes'))
 
    #louvain,fast
    set.seed(666)
    mds <- igraph::cluster_louvain(mygraph, weights = NULL, resolution = reso1)
    mdsdf <- data.frame(node=names(igraph::membership(mds)),module=as.vector(igraph::membership(mds)))
    mdsdf <- mdsdf[order(mdsdf$module),]
    hyp <- data.frame()
    
    for(md in unique(mdsdf$module)){
        valuek <- valuek
        mdnode <- subset(mdsdf,module==md)[,'node']
        M <- length(intersect(seedgene,mdnode))
        N <- length(mdnode)
        n <- valuen
        if(n >= valuek & N > M){
            prob <- suppressWarnings(phyper(q=valuek-1, m=M, n=N-M, k=n,lower.tail=FALSE))
            prob <- data.frame(module=md,pvalue=prob)
            hyp <- rbind(hyp,prob)
        } else{
            print(paste0('Warning, the number of nodes in module',md,' is less than valueq, ignore module',md,' ;number of module genes is too low to conduct phyper, if you want less modules and more module genes, try to reduce your reso1'))
        }
        #针对超几何分布画图
        if(grepl('\\/',celltype)){
            filename <- gsub('\\/','',celltype)
        } else{
            filename <- celltype
        }
        x <- 0:valuen
        cdf <- suppressWarnings(phyper(x, m=M, n=N-M, k=valuen,lower.tail=FALSE))
        df <- data.frame(seed_genes = x, Cumulative_Prob = cdf)
        x_obs <- valuek
        df$Right_Tail <- 1 - cdf
        p_value <- suppressWarnings(phyper(x_obs - 1, m=M, n=N-M, k=valuen, lower.tail = FALSE))
        tryCatch({
            if(p_value<0.05){
                p <- ggplot2::ggplot(df, ggplot2::aes(x = seed_genes, y = Right_Tail)) +
                ggplot2::theme_bw()+
                ggplot2::geom_line(color = "red", size = 1) +
                ggplot2::geom_vline(xintercept = x_obs, linetype = "dashed", color = "blue") +
                ggplot2::annotate("text", x = x_obs, y = 0.05, label = paste0("P(X>=",valuek,')'), hjust = -0.1, color = "blue") +
                ggplot2::annotate("text", x = x_obs, y = 0.1, label = paste0('Pvalue=',formatC(p_value, format = "e", digits = 2)), hjust = -0.1, color = "blue") +
                ggplot2::labs(title = paste0("Hypergeometric Right-Tail Probability: P(X>=",valuek,')'),
                    x = paste0("Number of Seed Genes in Sample (",valuek,') of Module',md),
                    y = "Right Tail Probability") 
                ggplot2::ggsave(p,file=paste0(filename,'.Module',md,'.Hyp.pdf'),width=5,height=5)
        }
        },error=function(e){})

        #用hyp画一个条形图
        if(nrow(hyp)>=8){
            colourCount <- nrow(hyp)
            getPalette = colorRampPalette(c('#C7E873','#A1E5D7','#C2BEE3','#FF887E','#85AEC9','#FACBE5','#FDB46F','#FCF6B8'))
            mycol <- getPalette(colourCount)
        } else{
            mycol <- c('#C7E873','#A1E5D7','#C2BEE3','#FF887E','#85AEC9','#FACBE5','#FDB46F','#FCF6B8')
        }
        names(mycol) <- hyp$module
        data <- na.omit(hyp)
        data$y <- -log10(data$pvalue)
        data <- data[order(data$y,decreasing=T),]
        data$module <- factor(data$module,levels=data$module)
        p <- ggplot2::ggplot(data, ggplot2::aes(x = module, y = y, fill = module)) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.6) +
        ggplot2::scale_fill_manual(values=mycol) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")+
        ggplot2::labs(title = "Modules enriched with seed genes",
            x = "Module",y = "-log10(pvalue)",fill = "module") +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(
            plot.title = element_text(hjust = 0.5,size=18),
            axis.text = element_text(color = "black"),
            panel.grid.major.x = element_blank(),
            legend.position='none'
        )
        ggplot2::ggsave(p,filename=paste0(filename,'.sigModule.pdf'),width=5,height=5)
    }
    
    # ----- RWR -----
    rwr <- function(W, p0, gamma = 0.5, epsilon = 1e-10, max_iter = maxiter) {
    p_prev <- p0
    for (i in 1:max_iter) {
        p_next <- (1 - gamma) * W %*% p_prev + gamma * p0
        if (sum(abs(p_next - p_prev)) < epsilon) {
        break
        }
        p_prev <- p_next
    }
    return(as.vector(p_next))
    }
    adj_matrix <- igraph::as_adjacency_matrix(mygraph)
    normalize_columns <- function(mat) {
    col_sums <- colSums(mat)
    mat / matrix(rep(col_sums, each = nrow(mat)), nrow = nrow(mat))
    }
    W <- normalize_columns(adj_matrix)
    p0 <- ifelse(colnames(W) %in% seedgene,1,0)
    names(p0) <- colnames(W)
    # ----- rn RWR -----
    set.seed(666)    
    p_final <- rwr(W, p0, gamma = 0.5)
    names(p_final) <- names(p0)
    p_final <- data.frame(node=names(p_final),rwr=as.vector(p_final))
    p_final <- p_final[order(p_final$rwr,decreasing=T),]
    sigmd <- subset(hyp,pvalue<0.05)[,'module']
    mdtarget <- subset(mdsdf,module %in% sigmd)[,'node'] %>% intersect(.,targets)
    ctnetclass <- new('CTNet',features=list(seedgene=seedgene,target=mdtarget),
                                ctnets=list(nets=aggnet,graph=mygraph),
                                ctmodule=list(modules=mdsdf,phyerMod=hyp),
                                RWR=list(rwr=p_final))
    return(ctnetclass)
}

#' @describeIn CTNet show method
#' @param object CTNet class
#' @export
setMethod(
  "show",
  "CTNet",
  function(object) {
    cat("An object of class 'CTNet'\n")
    cat(" Slots:\n")
    cat("  features:", length(object@features), "items\n")
    cat("  ctnets:", length(object@ctnets), "items\n")
    cat("  ctmodule:", length(object@ctmodule), "items\n")
    cat("  RWR:", length(object@RWR), "items\n")
  }
)