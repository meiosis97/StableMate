# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Built-in \code{st2e} proposal selection functions.
#'
#' Given objective scores of multiple proposals calculated by objective functions
#'  \code{\link[StableMate]{obj_fun}} provided to \code{\link[StableMate]{st2e}},
#' defines how to select the best proposal to make a step. The customized selection function should have same input as described below,
#' and output the index of the best proposal.
#'
#' @param scores A numerical vector of objective scores.
#'
#' @return The index of the best proposal.
#'
#' @export
#'
#' @seealso \code{\link[StableMate]{obj_fun}}
#' @name sel_fun
#' @rdname sel_fun
sel_min <- function(scores){
  which.min(scores)
}
