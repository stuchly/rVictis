#' K Nearest Neighbour Search
#'
#' Uses a annoylib to find K nearest neighbors in dataset
#' see https://github.com/spotify/annoy/
#'
#' @param data matrix; input data matrix
#' @param K integer; number of nearest neighbours
#' @param trees integer; number of trees in annoylib; more trees - more precise result
#' @param return_dist bool; if TRUE a list of indices and distances is returned
#'
#' @details If used as precomputed matrix for Rphenoannoy, discard first column - 1st nearest 
#' neigbor is the point itself
#' 
#' @return a n-by-k matrix of neighbors indices or list(IND,DIST) see parameter return_dist 
#'
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' neighbors <- knn.annoy(data, k=10)
#'
#' @export
knn.annoy<-function(data,K=30,trees=150){
  res<-knn_annoy(data,K,trees)
  return(matrix(as.integer(do.call("rbind",res[[1]])),ncol=K+1))
}