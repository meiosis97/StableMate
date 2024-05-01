# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' print function for \code{st2e} object.
#'
#' @param st2e_obj A \code{st2e} object returned by \code{\link[StableMate]{st2e}}.
#' @param imp_type Either 'conditional' or 'joint' to make selections based on conditional importance scores
#' or joint importance scores.
#' @param sigthresh Numerical; Significance threshold.
#'
#' @export
#'
#' @method print st2e
#' @rdname print.st2e
print.st2e <- function(st2e_obj, imp_type = c('conditional','joint'), sigthresh = 1-exp(-5)){
  cat('----------------------------------------------------\n')
  cat("Summary of the objective scores of selections in the ensemble:\n")
  print(summary(st2e_obj$obj_scores))
  cat('----------------------------------------------------\n')

  # Check the type of importance scores (and their significance) used to make statistically significant selections
  imp_type <- match.arg(imp_type)

  # Check whether importance scores were calculated
  if(is.null(st2e_obj$imp_scores)){
    cat("Importance scores were not found in the st2e object. Run calc_imp before to make statistically significant selections.")
    selected <- NULL

  }else{
    # Extract importance scores
    imp_obj <- st2e_obj$imp_scores[[imp_type]]

    imp <- imp_obj$importance[-1]
    sig <- imp_obj$significance[-1]

    cat('Selections are made based on', imp_type, 'importance scores \n')
    if(sum(sig>sigthresh)>0){
      cat('Predictors selected significantly more often than the pseudo-predictor are: \n')
      selected <- names(sig)[sig>sigthresh]
      cat(selected,'\n\n')
    }else{
      cat('Predictors with top 5 importancescores are: \n')
      selected <- names(imp)[rank(-imp, ties.method = 'min')<=5]
      cat(selected,'\n\n')
    }

  }

  invisible(list(selected = selected))

}


#' plot function for \code{st2e} object.
#'
#' @param st2e_obj A \code{st2e} object returned by \code{\link[StableMate]{st2e}}.
#' @param imp_type Either 'conditional' or 'joint' to show conditional importance scores
#' or joint importance scores.
#' @param fill_by User provided filling, should be a data frame of the same row length as the column size of the ensemble (exclude the pseudo-predictor),
#' Rows should be labelled by predictor names. The column name will be used as the legend title.
#' @param colors A vector of color values passed to \code{\link[ggplot]{scale_fill_manual}} and \code{\link[ggplot]{scale_color_manual}}.
#' @param sigthresh Numerical; Significance threshold.
#' @param plot_density Logical; If TRUE, bootstrap density of importance scores will be shown.
#' @param bt_prop Numerical; Proportion of bootstrap importance scores to show. Should be a value between 0 and 1.
#' @param base_size base_size passed to \code{\link[ggplot]{theme_classic}} .
#' @param show_labels Logical; If TRUE, show labels of predictors on x-axis.
#' @param labels User provided x label, It should be a vector named by the predictor names.
#' @param label_size label_size will be passed to \code{\link[ggplot]{element_text}} to control x axis label size.
#' @param box_size \code{\link[ggplot]{geom_boxplot}} size aesthetic.
#'
#' @export
#'
#' @method plot st2e
#' @rdname plot.st2e
plot.st2e <- function(st2e_obj, imp_type = c('conditional','joint'),
                      fill_by = NULL, colors = NULL, sigthresh = 1-exp(-5),
                      plot_density = FALSE, bt_prop = 0.1, base_size = 15,
                      show_labels = TRUE, labels = NULL, label_size = 10, box_size = 0.1){

  # Check the type of importance scores (and their significance) used to make statistically significant selections
  imp_type <- match.arg(imp_type)

  # Check whether importance scores were calculated
  if(is.null(st2e_obj$imp_scores)){
    cat("Importance scores were not found in the st2e object. Run calc_imp before to make statistically significant selections.")
  }

  # Extract importance scores
  imp_obj <- st2e_obj$imp_scores[[imp_type]]

  # Name of predictors
  pred_names <- names(imp_obj$importance)[-1]

  # The number of predictors
  P <- length(pred_names)

  # The number of bootstraps
  B <- nrow(imp_obj$bt_imp)


  ############### Code for label tuning ###############
  # Whether to show column labels (x-ticks)
  if(show_labels){
    # User provided x label, It should be a vector named by the predictor names
    # to be shown
    if(!is.null(labels)){
      if(is.null(names(labels))) names(labels) <- labels
      tmp <- rep('',P)
      names(tmp) <- pred_names
      tmp[names(labels)] <- labels
      labels <- tmp

      # By default, show the labels of selected predictors
    }else{
      tmp <- rep('',P)
      names(tmp) <- pred_names
      tmp[imp_obj$significance[-1]>sigthresh] <- pred_names[imp_obj$significance[-1]>sigthresh]
      labels <- tmp

    }

  }


  ############### Code for the boxplot ###############
  # Plot the bootstrapped distribution of the differences in importance scores between the pseudo predictor and each true predictor
  if(plot_density){
    # Randomly choose p percent of boostrap samples to show
    if(!is.null(bt_prop)) bt_id <- sample(1:B, bt_prop * B) else bt_id <- 1:B

    bt_imp_dff <- (imp_obj$bt_imp[bt_id,-1] - imp_obj$bt_imp[bt_id,1]) %>% as.matrix() %>% as.data.frame() %>%
      reshape2::melt(value.name = 'imp_diff')

    # Make a boxplot of the differences in importance scores between the pseudo predictor and each true predictor
    g <- ggplot2::ggplot(data = bt_imp_dff) + ggplot2::ylim(-1.001,1.001)  +
      ggplot2::theme_grey(base_size = base_size) +
      ggplot2::geom_boxplot(ggplot2::aes(variable, imp_diff, color = fill_by, fill = fill_by), size = box_size) +
      ggplot2::xlab('Predictors') + ggplot2::ylab(paste('Differences in', imp_type,'importance scores (Pseudo vs True)')) +
      ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')


    ############### Code for color tuning ###############
    # Create default filling according to the final selection result.
    if(is.null(fill_by)){
      fill_by <- rep('Not Selected', P)
      fill_by[imp_obj$significance[-1]>sigthresh] <- 'Selected'
      names(fill_by) <- pred_names
      fill_by <- fill_by[as.character(bt_imp_dff$variable)]

      g <- g + ggplot2::scale_fill_manual('Selection',
                                          values = c('lightgrey','red'))+
        ggplot2::scale_color_manual('Selection', values = c('lightgrey','red'))

    }else{
      # User provided filling, should be a data frame of the same row length as
      # the column size of the ensemble (exclude the pseudo-predictor),
      # Rows should be labelled by predictor names. The column name will be
      # used as the legend title
      legend_title <- colnames(fill_by)
      fill_by <- fill_by[as.character(bt_imp_dff$variable),1]

      # Color schemes
      if(is.null(colors)){
        g <- g + ggplot2::scale_fill_discrete(legend_title)+
          ggplot2::scale_color_discrete(legend_title)

      }else{
        g <- g + ggplot2::scale_fill_manual(legend_title,
                                            values = colors)+
          ggplot2::scale_color_manual(legend_title, values = colors)

      }

    }


    ############### Code for the barplot ###############
  }else{
    g <- ggplot2::ggplot() + ggplot2::ylim(0,1) + ggplot2::theme_classic(base_size = base_size) +
      ggplot2::geom_col(ggplot2::aes(x = pred_names, y = imp_obj$importance[-1], fill = fill_by), size = box_size)+
      ggplot2::geom_hline(yintercept = imp_obj$importance[1], linetype = 'dashed') +
      ggplot2::xlab('Predictors') +
      if(imp_type == 'conditional') ggplot2::ylab('Conditional mportance scores') else ggplot2::ylab('Joint importance scores')


    ############### Code for color tuning ###############
    # Create default filling according to the final selection result.
    if(is.null(fill_by)){
      fill_by <- rep('Not Selected', P)
      fill_by[imp_obj$significance[-1]>sigthresh] <- 'Selected'

      g <- g + ggplot2::scale_fill_manual('Selection',
                                          values = c('lightgrey','red'))+
        ggplot2::scale_color_manual('Selection', values = c('lightgrey','red'))

    }else{
      # User provided filling, should be a data frame of the same row length as
      # the column size of the ensemble (exclude the pseudo-predictor),
      # Rows should be labelled by predictor names. The column name will be
      # used as the legend title
      legend_title <- colnames(fill_by)
      fill_by <- fill_by[pred_names,1]

      # Color schemes
      if(is.null(colors)){
        g <- g + ggplot2::scale_fill_discrete(legend_title)+
          ggplot2::scale_color_discrete(legend_title)

      }else{
        g <- g + ggplot2::scale_fill_manual(legend_title,
                                            values = colors)+
          ggplot2::scale_color_manual(legend_title, values = colors)

      }

    }

  }

  g + ggplot2::theme(axis.ticks.x =  ggplot2::element_blank(),
                     axis.text.x= ggplot2::element_text(size = label_size)) +
    ggplot2::scale_x_discrete(labels = labels)

}


#' predict function for \code{st2e} object.
#'
#' @param st2e_obj A \code{st2e} object returned by \code{\link[StableMate]{st2e}}.
#' @param X New predictor data that can be recognized by \code{st2e} regression models fitted by the \code{reg_fun} argument passed to \code{\link[StableMate]{st2e}}.
#' @param prune Logical; If TRUE, weight selections by ranking quantiles of their objective scores when aggregating regression models.
#' @param prefilt_scores A numerical vector that contains some weighting measurements of predictor pools. Predictor pools will be
#' weighted by ranking quantiles of \code{prefilt_scores} when aggregating regression models.
#' @param pred_fun The function that makes prediction on the response given new predictor data \code{X} and regression models
#' in \code{st2e}. By default, expect regression models fitted by code{\link[StableMate]{reg_ols}}.
#' See \code{\link[StableMate]{reg_fun}} for more options and how to customize this function.
#'
#' @export
#'
#' @method predict st2e
#' @rdname predict.st2e
predict.st2e <- function(st2e_obj, X, prune = TRUE, prefilt_scores = NULL, pred_fun = pred_ols){
  K <- nrow(st2e_obj$ensemble)
  # Calculate the weights of st2e_obj selections if prune the ensemble
  if(prune) q_sel <- rank(-st2e_obj$obj_scores, ties.method = 'average')/K else q_sel <- rep(1,K)
  # Calculate the weights of prior selections if prune the ensemble
  if(prune & !is.null(prefilt_scores)) q_prefilt <- rank(-prefilt_scores, ties.method = 'average')/K else q_prefilt <-  rep(1,K)
  q_sel <- q_sel * q_prefilt

  # Build prediction ensemble
  prediction_ensemble <- matrix(0, ncol = nrow(X), nrow = K)
  pred_selected <- names(which(colSums(st2e_obj$ensemble[,-1])>0)) # Check whether all predictors that were used to build models are also in the new data
  if( !all( pred_selected %in% colnames(X) )) stop('At least one predictor in the ensemble was not found in the data')

  #Each model makes prediction
  for(i in 1:K){
    pred_selected <- names(which(st2e_obj$ensemble[i,-1]==1))
    prediction_ensemble[i,] <- pred_fun(st2e_obj$model[[i]], X[, pred_selected, drop =F])
  }

  aggregated_prediction <- colSums(prediction_ensemble*q_sel)/sum(q_sel)

  attr(aggregated_prediction, which = 'ensemble') <- prediction_ensemble

  aggregated_prediction
}
