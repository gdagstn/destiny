#' kNN search
#' 
#' Approximate k nearest neighbor search with flexible distance function.
#' 
#' @param data      Data matrix
#' @param query     Query matrix. Leave it out to use \code{data} as query
#' @param k         Number of nearest neighbors
#' @param distance  Distance metric to use. Allowed measures: Euclidean distance (default), cosine distance (\eqn{1-corr(c_1, c_2)}) or rank correlation distance (\eqn{1-corr(rank(c_1), rank(c_2))})
#' @param method    Method to use. \code{'kmknn'} and \code{'covertree'} are exact methods; \code{'hnsw'} is approximate. (default: kmknn)
#' @param sym       Return a symmetric matrix (as long as query is NULL)?
#' @param verbose   Show a progressbar? (default: FALSE)
#' @param BNPARAM	A \code{\link[BiocNeighbors]{BiocNeighborParam}} object specifying the algorithm to use. This can be left empty if \code{'method'} is specified.
#' 
#' @return A \code{\link{list}} with the entries:
#' \describe{
#'   \item{\code{index}}{A \eqn{nrow(data) \times k} \link{integer} \link{matrix} containing the indices of the k nearest neighbors for each cell.}
#'   \item{\code{dist}}{A \eqn{nrow(data) \times k} \link{double} \link{matrix} containing the distances to the k nearest neighbors for each cell.}
#'   \item{\code{dist_mat}}{
#'     A \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} if \code{sym == TRUE},
#'     else a \code{\link[Matrix:dsCMatrix-class]{dsCMatrix}} (\eqn{nrow(query) \times nrow(data)}).
#'     Any zero in the matrix (except for the diagonal) indicates that the cells in the corresponding pair are close neighbors.
#'   }
#' }
#' 
#' @rdname knn
#' @importFrom BiocNeighbors findKNN queryKNN HnswParam KmknnParam
#' @importFrom matrixStats colRanks
#' @export
find_knn <- function(
	data, k,
	query = NULL,
	distance = c('euclidean', 'cosine', 'rankcor', 'l2'),
	method = c('kmknn', 'covertree', 'hnsw'),
	sym = TRUE,
	verbose = FALSE,
	BNPARAM = NULL
) {
	#p <- utils::modifyList(formals(RcppHNSW::hnsw_knn), list(...))
	if(!is.null(method)) {
		method <- match.arg(method)
	} else if(is.null(method) & is.null(BNPARAM)) {
		stop('You must either specify a method or a BNPARAM')
	} 

	distance <- match.arg(distance)

	if(!is.null(method)) {
		if (!is.double(data) & method == 'covertree') {	
			warning('find_knn does not yet support sparse matrices, converting data to a dense matrix.')
			data <- as.matrix(data)
		}
		if (method == 'covertree') {
			message("Using method \'covertree\'")
			return(knn.covertree::find_knn(data, k, query = query, distance = distance, sym = sym))
		}
	}
	
	if (distance == 'rankcor') {
		distance <- 'cosine'
		data <- t(colRanks(data, ties.method = "random", preserveShape = FALSE, decreasing = TRUE))
		if (!is.null(query)) query <- colRanks(query)
	}
	
    dchar = unlist(strsplit(distance, ''))
    distance <- paste0(toupper(dchar[1L]), paste0(dchar[-1L], collapse = ''))

    if(is.null(BNPARAM)) {
		if(method == 'kmknn') {
            knnparam = KmknnParam(distance = distance)
        } else if(method == 'hnsw') {
            knnparam = HnswParam(distance = distance)
        }
	} else {
		knnparam = BNPARAM
	}

	if (is.null(query)) {
        knn <- findKNN(X = data, k = k, BNPARAM = knnparam)
	} else {
        knn <- queryKNN(X = data, query = query, k = k, BNPARAM = knnparam)
	}
  
	# R matrices are column-major, so as.vector(m) == c(m[, 1], m[, 2], ...)
	knn$dist_mat <- sparseMatrix(
		rep(seq_len(nrow(knn$index)), k),
		as.vector(knn$index),
		x = as.vector(knn$distance),
		dims = c(nrow(if (is.null(query)) data else query), nrow(data))
	)
	if (is.null(query)) {
		if (sym) knn$dist_mat <- symmetricise(knn$dist_mat)
		nms <- rownames(data)
	} else {
		nms <- rownames(query)
	}
	rownames(knn$dist_mat) <- rownames(knn$index) <- rownames(knn$distance) <- nms
	colnames(knn$dist_mat) <- rownames(data)
	names(knn)[2] <- 'dist'
	knn
}


# (double generic columnsparse to ... symmetric ...: dgCMatrix -> dsCMatrix)
# retain all differences fully. symmpart halves them in the case of trans_p[i,j] == 0 && trans_p[j,i] > 0
# TODO: could be more efficient
symmetricise <- function(dist_asym) {
	dist_sym <- symmpart(dist_asym) + abs(forceSymmetric(skewpart(dist_asym), 'U'))
	as(dist_sym, 'symmetricMatrix')
}
