% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn.r
\name{find_knn}
\alias{find_knn}
\title{kNN search}
\usage{
find_knn(
  data,
  k,
  query = NULL,
  distance = c("euclidean", "cosine", "rankcor", "l2"),
  method = c("kmknn", "covertree", "hnsw"),
  sym = TRUE,
  verbose = FALSE,
  BNPARAM = NULL
)
}
\arguments{
\item{data}{Data matrix}

\item{k}{Number of nearest neighbors}

\item{query}{Query matrix. Leave it out to use \code{data} as query}

\item{distance}{Distance metric to use. Allowed measures: Euclidean distance (default), cosine distance (\eqn{1-corr(c_1, c_2)}) or rank correlation distance (\eqn{1-corr(rank(c_1), rank(c_2))})}

\item{method}{Method to use. \code{'kmknn'} and \code{'covertree'} are exact methods; \code{'hnsw'} is approximate. (default: kmknn)}

\item{sym}{Return a symmetric matrix (as long as query is NULL)?}

\item{verbose}{Show a progressbar? (default: FALSE)}

\item{BNPARAM}{A \code{\link[BiocNeighbors]{BiocNeighborParam}} object specifying the algorithm to use. This can be left empty if \code{'method'} is specified.}
}
\value{
A \code{\link{list}} with the entries:
\describe{
  \item{\code{index}}{A \eqn{nrow(data) \times k} \link{integer} \link{matrix} containing the indices of the k nearest neighbors for each cell.}
  \item{\code{dist}}{A \eqn{nrow(data) \times k} \link{double} \link{matrix} containing the distances to the k nearest neighbors for each cell.}
  \item{\code{dist_mat}}{
    A \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} if \code{sym == TRUE},
    else a \code{\link[Matrix:dsCMatrix-class]{dsCMatrix}} (\eqn{nrow(query) \times nrow(data)}).
    Any zero in the matrix (except for the diagonal) indicates that the cells in the corresponding pair are close neighbors.
  }
}
}
\description{
Approximate k nearest neighbor search with flexible distance function.
}
