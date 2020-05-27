# data.R - DESC
# mydas/R/data.R
#' mydas datasets
#'
#' Example datasets for the classes defined in FLCore.
#'
#' \item{om, \code{\link{FLStock}}}{Operation Model conditioned on turbot life history parameters'  \code{lh}}
#' \item{eq, \code{\link{FLBRP}}}{Expected dynamics'  \code{lh}}
#' 
#' See [mydas_turbot](https://3o2y9wugzp1kfxr5hvzgzq-on.drv.tw/MyDas/om/mydas_turbot.html) vignette for more infomation
#' 
#' Datasets can be loaded by issuing the \code{data} command, as in
#' \code{data(om)}
#' \code{data(eq)}
#' 
#' @name datasets
#' @aliases 
#' @seealso \linkS4class{FLStock}, \linkS4class{FLBRP}, \linkS4class{FLPar}
#' @keywords datasets
#' @examples
#' \dontrun{
#' data(om)
#' summary(om)
#' }
NULL
