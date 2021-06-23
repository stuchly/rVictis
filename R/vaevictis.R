#' \code{vaevictis}
#' @param x numeric data matrix.
#' @param dim integer (default 2); dimension of reduced data
#' @param vsplit double (default 0.1); percentage of data used as validation step in "vaevictis".
#' By default split is performed before shuffling - see \code{shuffle}
#' @param enc_shape integer vector (default \code{c(128,128,128)}); shape (depth and wisth) of the encoder.
#' @param dec_shape integer vector (default \code{c(128,128,128)}); shape (depth and wisth) of decoder.
#' @param perplexity double (default 10.); perplexity for tsne regularisation see https://www.nature.com/articles/s41467-018-04368-5.
#' @param batch_size integer (default 512); batch size for "vaevictis" training.
#' @param epochs integer; maximum number of epochs for "vaevictis" training.
#' @param patience integer; maxim patience for for "vaevictis" training (early stopping).
#' @param ww vector double; weights for vaevictis in this order - tsne regularization, cense regularization, umap regularization, ivis pn loss, reconstruction error, KL divergence
#' @param margin double; ivis pn loss margin
#' @param shuffle logical; shuffle data before validation split; involves recomputation of KNN matrix
#' @param load_model character vector of 2 components; paths to files created by by vae$save(file1,file2) - model is loaded and applied
#' @param upsample named list \code{list(labels=,N=,takeall=, cluster=, sampsize=, samples=,method=,pheno_K=,trees=)} or \code{NULL};  sample events by labels;
#' if NULL nothing happens, \code{N} events per label, \code{labels} vector of labels, takes all events from labels in character vector \code{takeall},
#' \code{cluster} number of clusters to create vector of labels. If \code{upsample} is not \code{NULL}, \code{N} and unless \code{method}=="phenograph" either \code{labels} or \code{cluster} must be
#' provided. Clustering is performed if \code{cluster} is not \code{NULL}. \code{samples} and \code{sampsize} are parameters for cluster::clara,
#' if not provided default values are used.
#' if \code{method}=="phenograph", parameter \code{pheno_K} (number of NN for phenograph) could be provided (by default \code{pheno_K}==\code{k})
#' @param k integer - number of NN for ivis regularisation
#'
#' @param map_only logical; return only the mapping not layout
#'
#' @details see https://github.com/stuchly/vaevictis
#'
#' @return  returns named list list(layout,vae,map) - reduced data, model from python module vaevictis and a pointer to python function to perform the reduction
#'
#' @examples
#' library(rVictis)
#'
#' res <- vaevictis(iris[,1:4],upsample=list(N=100,cluster=10)) ## 10 clusters by clara, 100 examples from each cluster
#'
#' plot(res[[1]],col=as.factor(iris[,5])) ## plot reduced data
#'
#' rediris <- res[[3]](as.matrix(iris[,1:4])) ##apply train reduction on "new" data, now we pass the data to python function - must be matrix
#'
#' plot(rediris,col=as.factor(iris[,5]))
#'
#' @export
vaevictis <-function(x,
                     dim = 2,
                     vsplit = 0.1,
                     enc_shape = c(128, 128, 128),
                     dec_shape = c(128, 128, 128),
                     perplexity = 10.,
                     batch_size = 512L,
                     epochs = 100L,
                     patience = 0L,
                     ivis_pretrain=0,
                     ww=c(10.,10.,1.,1.,0.,1.),
                     margin=1.,
                     shuffle=FALSE,
                     load_model=NULL,
                     upsample=NULL,
                     k=30,
                     map_only=FALSE)
{

    x<-as.matrix(x)
    vv = reticulate::import("vaevictis")
    if (!is.null(load_model)){
      model <- vv$loadModel(config_file = load_model[1],weights_file = load_model[2])
      layout <- model[[2]](x)
      vae <- model
      map<-function(x) model[[2]](x)
    } else {
        if (!is.null(upsample)){
            if (is.null(upsample$method)) upsample$method<-"clara"
            if (is.null(upsample$cluster) & upsample$method!="phenograph") if (!is.null(upsample$labels)) labl<-as.factor(upsample$labels) else  stop("labels not provided!")
            else {
                if (upsample$method=="clara"){
                    message("running clara clustering...")
                    sampsize<-ifelse(is.null(upsample$sampsize),min(nrow(x),40+2*upsample$cluster),upsample$sampsize)
                    samples<-ifelse(is.null(upsample$samples),5,upsample$samples)
                    labl<-as.factor(cluster::clara(x,k = upsample$cluster,sampsize = sampsize,samples = samples)$clustering)
                    message("~done\n")
                } else if (upsample$method=="phenograph"){
                    if (is.null(upsample$pheno_K)) upsample$pheno_K<-k
                    if(is.null(upsample$trees)) upsample$trees<-150
                    knn<-knn.annoy(x,K=upsample$pheno_K,trees=upsample$trees)[,-1]+1
                    labl<-as.factor(Rphenoannoy(k,knn)$membership)
    
                } else if (upsample$method=="kmeans"){
                    labl<-as.factor(kmeans(x,centers=upsample$cluster)$cluster)
                } else stop("not implemented")
    
            }
            ss<-.upsample.labels(labl,N=upsample$N,takeall = upsample$takeall)
                 knn<-knn.annoy(x[ss,],K=k)
            layout = vv$dimred(
              x[ss,],
              as.integer(dim),
              vsplit,
              enc_shape,
              dec_shape,
              perplexity,
              as.integer(batch_size),
              as.integer(epochs),
              as.integer(patience),
              as.integer(ivis_pretrain),
              ww,
              "euclidean",
              margin,
              NULL,
              NULL,
              knn
            )
        } else {
            if (shuffle) sshuf<-sample(nrow(x)) else sshuf<-1:nrow(x)
            knn<-knn.annoy(x[sshuf,],K=k)
            layout = vv$dimred(
              x[sshuf,],
              as.integer(dim),
              vsplit,
              enc_shape,
              dec_shape,
              perplexity,
              as.integer(batch_size),
              as.integer(epochs),
              as.integer(patience),
              as.integer(ivis_pretrain),
              ww,
              "euclidean",
              margin,
              NULL,
              NULL,
              knn
            )
        }

        map <- layout[[2]]
        vae <- layout[[3]]
        if (!map_only) layout <- layout[[2]](x) else layout<-NULL

    }
    return(list(layout=layout,vae=vae,map=map))
}



Rphenoannoy <-  function(k = 30, neighborMatrix,  threshold=0L) {


    cat( " Compute jaccard coefficient between nearest-neighbor sets...\n")
    message("Presorting knn...\n")
    nbh <- neighborMatrix[, 1:min(k, ncol(neighborMatrix))]
    t21 <- system.time(nbh <- t(apply(nbh, 1, sort)))
    cat("presorting DONE ~", t21[3], "s\n", " Start jaccard\n")
    t2 <- system.time(links <- jaccard_coeff_true_parallel(nbh,threshold))

    cat("DONE ~",
        t2[3],
        "s\n",
        " Build undirected graph from the weighted links...")
    links <- links[links[, 1] > 0,]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
                                        # t3 <-
                                        #   system.time(g <- graph.data.frame(relations, directed = FALSE))
    t3<-system.time({
        N<-nrow(nbh)
        A<-Matrix::sparseMatrix(i=relations$from, j=relations$to, x=relations$weight, dims=c(N,N))
        rm(relations)
        g<-igraph::graph_from_adjacency_matrix(A,weighted = TRUE,mode = "max")
        rm(A)
    })
                                        # Other community detection algorithms:
                                        #    cluster_walktrap, cluster_spinglass,
                                        #    cluster_leading_eigen, cluster_edge_betweenness,
                                        #    cluster_fast_greedy, cluster_label_prop
    cat("DONE ~", t3[3], "s\n", " Run louvain clustering on the graph ...")
    t4 <- system.time(community <- igraph::cluster_louvain(g))
    cat("DONE ~", t4[3], "s\n")

    message("Run Rphenograph DONE, totally takes ", sum(c(t2[3], t3[3], t4[3])), "s.")
    cat("  Return a community class\n  -Modularity value:",
        igraph::modularity(community),
        "\n")
    cat("  -Number of clusters:", length(unique(igraph::membership(community))),"\n")

    return(community)
}
