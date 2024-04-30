print.stablemate <- function(stbm_obj, pred_imp_type = c('joint','conditional'),
                             stab_imp_type = c('conditional','joint'),  sigthresh = 1-exp(-5)){

  # Print out predictors selected as predictive
  cat('----------------------------------------------------\n')
  cat("Summary of the objective scores of the selections in the prediction ensemble:\n")
  print(summary(stbm_obj$prediction_ensemble$obj_scores))
  cat('----------------------------------------------------\n')
  pred_selected <- NULL

  # Check whether importance scores of predictivity selections were calculated
  if(is.null(stbm_obj$prediction_ensemble$imp_scores)){
    cat("Importance scores were not calculated on the prediction ensemble. Run calc_imp before to make statistically significant selections.")

  }else{
    # Check the type of importance scores (and their significance) used to make statistically significant selections
    pred_imp_type <- match.arg(pred_imp_type)

    # Extract the importance scores (and their significance) of predictive predictors
    pred_imp_obj <- stbm_obj$prediction_ensemble$imp_scores[[pred_imp_type]]
    pred_sel_imp <- pred_imp_obj$importance[-1]
    pred_sig <- pred_imp_obj$significance[-1]

    # Print out predictors selected as predictive
    cat('Selections of predictive predictors are made based on', pred_imp_type, 'importance scores \n')
    if(sum(pred_sig>sigthresh)>0){
      cat('Predictors selected significantly more often than the pseudo-predictor in the predictivity selection are: \n')
      pred_selected <- names(pred_sig)[pred_sig>sigthresh]
      cat(pred_selected,'\n\n')
    }else{
      cat('No predictor was found selected significantly more often than the pseudo-predictor in the predictivity selection. No selection will be made on stable or unstable predictors. \n\n')
    }

  }

  # Print out predictors selected as stable and predictive
  cat('----------------------------------------------------\n')
  cat("Summary of the objective scores of the selections in the stable ensemble:\n")
  print(summary(stbm_obj$stable_ensemble$obj_scores))
  cat('----------------------------------------------------\n')
  stab_selected = NULL;  nstab_selected = NULL

  # Check whether importance scores of stability selections were calculated
  if(is.null(stbm_obj$stable_ensemble$imp_scores)){
    cat("Importance scores were not calculated on the stable ensemble. Run calc_imp before to make statistically significant selections.")

  }else if(is.null(stbm_obj$prediction_ensemble$imp_scores)){
    cat("Selections on predictive and stable predictors cannot be made without calculating importance scores on the prediction ensemble.")

  }else{
    # Check the type of importance scores (and their significance) used to make statistically significant selections
    stab_imp_type <- match.arg(stab_imp_type)

    # Extract the importance scores (and their significance) of stable predictors
    stab_imp_obj <- stbm_obj$stable_ensemble$imp_scores[[stab_imp_type]]
    stab_sel_imp <- stab_imp_obj$importance[-1]
    stab_sig <- stab_imp_obj$significance[-1]
    nstab_sig <- stab_imp_obj$lt_significance[-1]

    # Print out predictors selected as predictive and stable
    if(!is.null(pred_selected)){
      cat('Selections of stale predictors are made based on', stab_imp_type, 'importance scores \n')
      if(sum(stab_sig>sigthresh)>0){
        stab_selected <- names(stab_sig)[stab_sig>sigthresh]
        stab_selected <- intersect(pred_selected, stab_selected)
        if(length(stab_selected) > 0){
          cat('Predictors selected as predictive and selected significantly more often than the pseudo-predictor in the stability selection are: \n')
          cat(stab_selected,'\n\n')
        }else{
          cat('No predictive predictor was found selected significantly more often than the pseudo-predictor in the stability selection. \n\n')
        }
      }else{
        cat('No predictor was found selected significantly more often than the pseudo-predictor in the stability selection. \n\n')
      }

      # Print out predictors selected as predictive but unstable
      if(sum(nstab_sig>sigthresh)>0){
        nstab_selected <- names(nstab_sig)[nstab_sig>sigthresh]
        nstab_selected <- intersect(pred_selected, nstab_selected)
        if(length(nstab_selected) > 0){
          cat('Predictors selected as predictive and selected significantly less often than the pseudo-predictor in the stability selection are: \n')
          cat(nstab_selected,'\n\n')
        }else{
          cat('No predictive predictor was found selected significantly less often than the pseudo-predictor in the stability selection. \n\n')
        }
      }else{
        cat('No predictor was found selected significantly less often than the pseudo-predictor in the stability selection. \n\n')
      }

    }

  }

  invisible(list(pred_selected = pred_selected, stab_selected = stab_selected,  nstab_selected = nstab_selected))

}


plot.stablemate <- function(stbm_obj, pred_imp_type = c('joint','conditional'),
                            stab_imp_type = c('conditional','joint'),
                            sigthresh = 1-exp(-5),
                            color_by = NULL, colors = NULL,
                            base_size = 15, point_size  = 2,
                            show_labels = TRUE, labels = NULL, label_size = 3, parse = FALSE,
                            pred_cutoff = NULL, show_density_of = NULL){


  ############### Initialization ###############

  # Check the type of importance scores (and their significance) used to make statistically significant selections
  pred_imp_type <- match.arg(pred_imp_type)
  stab_imp_type <- match.arg(stab_imp_type)

  # Check whether importance scores of predictivity selections were calculated
  if(is.null(stbm_obj$prediction_ensemble$imp_scores)){
    stop("Importance scores were not calculated on the prediction ensemble. Run calc_imp before to make statistically significant selections.")
  }

  # Check whether importance scores of stability selections were calculated
  if(is.null(stbm_obj$stable_ensemble$imp_scores)){
    stop("Importance scores were not calculated on the stable ensemble. Run calc_imp before to make statistically significant selections.")
  }


  ############### Calculate plotting elements ###############
  # Extract the importance scores (and their significance) of predictive predictors
  pred_imp_obj <- stbm_obj$prediction_ensemble$imp_scores[[pred_imp_type]]
  pred_imp <- pred_imp_obj$importance[-1]
  pred_sig <- pred_imp_obj$significance[-1]

  # Extract the importance scores (and their significance) of stable predictors
  stab_imp_obj <- stbm_obj$stable_ensemble$imp_scores[[stab_imp_type]]
  stab_imp <- stab_imp_obj$importance[-1]
  stab_sig <- stbm_obj$stable_ensemble$imp_scores[['conditional']]$significance[-1]
  nstab_sig <- stbm_obj$stable_ensemble$imp_scores[['conditional']]$lt_significance[-1]

  # Log transform the significance of the importance scores in stability selection
  # which can be used to color predictors in the variable selection plot.
  # Its a measure of how significant, either stable or unstable, each predictor is
  log_stab_sig <- 1 - (pmax(stab_sig, nstab_sig) - pmin(stab_sig, nstab_sig))
  log_stab_sig <- -log10(log_stab_sig)
  # Unstable predictors receive negative log-significance, whereas stable predictors receive
  # positive log-significance. The larger the absolute value is, the more stable (or unstable)
  # a predictors
  log_stab_sig[nstab_sig > stab_sig] <- -log_stab_sig[nstab_sig > stab_sig]
  log_stab_sig[nstab_sig == stab_sig] <- 0

  # Make the final selection based on the significance threshold provided
  is_pred <- pred_sig > sigthresh # Whether a predictor is predictive
  is_stab <- rep('Non-significant',length(stab_sig))
  is_stab[stab_sig > sigthresh] <- 'Stable'
  is_stab[nstab_sig > sigthresh] <- 'Unstable'
  is_stab <- factor(is_stab, levels = c('Stable','Non-significant','Unstable'))

  # Set the upper bound of log-pvalue
  B <- nrow(stab_imp_obj$bt_imp) # The size of the ensemble
  p_upbd <- log10(B+1)
  log_stab_sig[log_stab_sig == Inf] <- p_upbd
  log_stab_sig[log_stab_sig == -Inf] <- -p_upbd

  # Set up a color range, where the predictors with significance higher than the
  # threshold will be marked with specific color coding for stability
  q <- (p_upbd - min(-log10(sigthresh), p_upbd))/(2*p_upbd)

  # Set up predictors to color
  pred_names <- names(pred_imp)
  alpha <- rep(0, length(pred_names))
  names(alpha) <- pred_names
  if(!is.null(pred_cutoff)){
    alpha[pred_imp>pred_cutoff] <- 1
  }else{
    alpha[pred_sig>sigthresh] <- 1

  }

  # Set up labels
  if(show_labels){
    # User provided x label, It should be a vector named by the predictor names
    # to be shown
    if(!is.null(labels)){
      label_which <- ifelse(is.null(names(labels)), labels, names(labels))

      # By default, show the labels of predictive predictors
    }else if(is.null(pred_cutoff)){
      label_which <- labels <- pred_names[pred_sig>sigthresh]

    }else{
      label_which <- labels <- pred_names[pred_imp>pred_cutoff]

    }

  }


  ############### Make plot ###############
  g <- ggplot2::ggplot()

  # Draw bootstrap quantiles of the pseudo-predictor
  if(stab_imp_type == 'conditional'){
    xintercept <- quantile(pred_imp_obj$bt_imp[,1], c(1-sigthresh, sigthresh))
    yintercept <- quantile(stab_imp_obj$bt_imp[,1], c(1-sigthresh, sigthresh))

    g <- g + ggplot2::geom_vline(xintercept = xintercept, linetype = 'dashed', col = 'blue', alpha = 0.5) +
      ggplot2::geom_hline(yintercept = yintercept, linetype = 'dashed', col = 'blue', alpha = 0.5)+
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= 0, ymax = 1), fill = 'green', alpha = 0.2) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= yintercept[2], ymax = 1), fill = 'green', alpha = 0.2) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= 0, ymax = yintercept[1]), fill = 'blue', alpha = 0.1)

  }else{
    xintercept <- quantile(stbm_obj$prediction_ensemble$imp_scores$conditional$bt_imp[,1], c(1-sigthresh, sigthresh))
    slope <- quantile(stbm_obj$stable_ensemble$imp_scores$conditional$bt_imp[,1], c(1-sigthresh, sigthresh))
    yintercept <- c(0,0)

    g <- g + ggplot2::geom_vline(xintercept = xintercept, linetype = 'dashed', col = 'blue', alpha = 0.5) +
      ggplot2::geom_abline(slope = slope, intercept = yintercept, linetype = 'dashed', col = 'blue', alpha = 0.5)+
      ggplot2::geom_rect(ggplot2::aes(xmin = xintercept[2], xmax = 1, ymin= 0, ymax = 1), fill = 'green', alpha = 0.2) +
      ggplot2::geom_ribbon(ggplot2::aes(x = seq(xintercept[2], 1, length.out = 1000),
                                        ymin = seq(xintercept[2], 1, length.out = 1000) * slope[2] + yintercept[2],
                                        ymax = 1), fill = 'green', alpha = 0.2)+
      ggplot2::geom_ribbon(ggplot2::aes(x = seq(xintercept[2], 1, length.out = 1000),
                                        ymax = seq(xintercept[2], 1, length.out = 1000) * slope[1] + yintercept[1],
                                        ymin = 0), fill = 'blue', alpha = 0.1)
  }


  ############### Density plot ###############
  # Plot the 2d density of the bootstrapped importance scores of user selected predictors and the pseudo-predictor
  if(!is.null(show_density_of)){
    # Extract the bootstrapped importance scores of user selected predictors
    message(cat('Chosen to plot the bootstrapped selection probabilities of', show_density_of))
    bt_pred_imp <- pred_imp_obj$bt_imp[,show_density_of,drop=FALSE]
    bt_stab_imp <- stab_imp_obj$bt_imp[,show_density_of,drop=FALSE]

    getLevel <- function(x,y,prob= 2*sigthresh-1) {
      kk <- MASS::kde2d(x,y)
      dx <- diff(kk$x[1:2])
      dy <- diff(kk$y[1:2])
      sz <- sort(kk$z)
      c1 <- cumsum(sz) * dx * dy
      approx(c1, sz, xout = 1 - prob)$y
    }

    # Draw the density line of the pseudo-predictor
    L95 <- getLevel(pred_imp_obj$bt_imp[,1], stab_imp_obj$bt_imp[,1])
    g <- g +  ggplot2::geom_density_2d(ggplot2::aes(pred_imp_obj$bt_imp[,1], stab_imp_obj$bt_imp[,1]), color = 'black', breaks = L95)

    for(i in show_density_of){
      # Draw the density line of the ith user selected predictor
      L95 <- getLevel(bt_pred_imp[,i],bt_stab_imp[,i])
      g <- g +  ggplot2::geom_density_2d(ggplot2::aes_string(bt_pred_imp[,i], bt_stab_imp[,i]), color = 'purple', breaks = L95)

      # If significance was evaluated at per-predictor level (i.e, each true predictor is
      # benchmarked against the bootstrap distribution of the pseudo-predictor selected
      # under the pool of the true predictor)
      if(!is.null(pred_imp_obj$benchmark) | !is.null(stab_imp_obj$benchmark)){
        if(is.null(pred_imp_obj$benchmark)){
          benchmark_x <- pred_imp_obj$bt_imp[,1]

        }else{
          benchmark_x <- pred_imp_obj$benchmark[,i]

        }

        if(is.null(stab_prob_obj$benchmark)){
          benchmark_y <- stab_imp_obj$bt_imp[,1]

        }else{
          benchmark_y <- stab_imp_obj$benchmark[,i]

        }

        # Draw the density line of the pseudo-predictor selected under the pool of the true predictor
        L95 <- getLevel(benchmark_x, benchmark_y)
        g <- g +  ggplot2::geom_density_2d(ggplot2::aes_string(benchmark_x, benchmark_y), color = 'black', breaks = L95)
        g <- g +  ggplot2::annotate('text', x = mean(benchmark_x), y = mean(benchmark_y), label = i,
                                    color = 'black',size = label_size)

      }

    }

  }


  ############### The main plot ###############
  g <- g + ggplot2::scale_x_continuous(limits = c(-0.001,1.001), expand = c(0, 0.01)) +
    ggplot2::scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0.01)) + ggplot2::theme_grey(base_size = base_size) +
    ggplot2::geom_point(ggplot2::aes(pred_imp, stab_imp, shape = is_stab),
                        color = 'lightgrey',fill = 'lightgrey', alpha = 1-alpha, size = point_size)+
    ggrepel::geom_text_repel(ggplot2::aes(x = pred_imp[label_which], y = stab_imp[label_which], label = labels),
                             segment.color = 'grey', force = 1,max.overlaps = 20, size = label_size, segment.linetype = 'dashed',
                             parse = parse) +
    ggplot2::scale_shape_manual('Selection', values = c('Stable' = 24, 'Non-significant' = 21, 'Unstable' = 25))

  g <- g + ggplot2::geom_point(ggplot2::aes(pred_imp, stab_imp,
                                            fill = color_by, shape = is_stab, col = color_by), alpha = alpha, size = point_size)


  if(is.null(color_by)){
    color_by <- is_stab
    g <- g + ggplot2::scale_fill_manual('Selection', values = c('Stable' = "blue", 'Non-significant' = "orange", 'Unstable' = "red")) +
      ggplot2::scale_color_manual('Selection', values = c('Stable' = "blue", 'Non-significant' = "orange", 'Unstable' = "red"))

  }else if(length(color_by) == 1 & color_by[1] == 'p-val'){
    color_by <- log_stab_sig
    g <- g + ggplot2::scale_fill_gradientn('Signed Log10 P-val',colours =  c("red", "orange","orange","blue"),
                                           values = c(0, q, 1-q, 1), limits = c(-p_upbd,p_upbd))+
      ggplot2::scale_color_gradientn('Signed Log10 P-val',colours =  c("red", "orange", "orange","blue"),
                                     values = c(0, q, 1-q, 1), limits = c(-p_upbd,p_upbd))

  }else{
    # User provided coloring, should be a data frame of the same row length as
    # the column size of the ensemble (exclude the pseudo-predictor),
    # Rows should be labelled by predictor names. The column name will be
    # used as the legend title
    legend_title <- colnames(color_by)
    color_by <- color_by[as.character(pred_names),1]

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

  x <- "Importance: Predictive"
  if(stab_imp_type == 'conditional'){
    y <- "Importance: Stable"
  }else{
    y <- "Importance: Stable and Predictive"
  }

  g + ggplot2::xlab(x) + ggplot2::ylab(y)

}


predict.stablemate <- function(stbm_obj, X, assay = c('stable_ensemble','prediction_ensemble'),
                               prune = TRUE, prefilt_scores = NULL, pred_fun = pred_ols){
  assay <- match.arg(assay)

  # Check whether the stability selection is based on the predictivity selection
  if(stbm_obj$matched){
    K <- nrow(stbm_obj$stable_ensemble$ensemble)

    if(prune){
      # Calculate the weights of stability selections if prune the ensemble (In stablemate, its the ranking quantiles of pss scores)
      q_sel <- rank(-stbm_obj$stable_ensemble$obj_scores, ties.method = 'average')/K
      # Calculate the weights of st2e_obj2 selections if prune the ensemble (In stablemate, its the ranking quantiles of bic scores)
      q_sel2 <- rank(-stbm_obj$prediction_ensemble$obj_scores, ties.method = 'average')/K
      # Calculate the weights of prior selections if prune the ensemble (In stablemate, its the ranking quantiles of the objective scores of the pre-filtering method applied)
      if(!is.null(prefilt_scores)) q_prefilt <- rank(-prefilt_scores, ties.method = 'average')/K else q_prefilt <-  rep(1,K)

    }else{
      q_sel <- q_sel2 <- q_prefilt <- rep(1,K)

    }

    # Weight selections by whether they are both important for the pre-filtering and predictivity selection
    q_sel2 <- q_sel2*q_prefilt
    # Weight selections by whether they are both important for the pre-filtering and predictivity selection and stability selection
    q_sel <- q_sel*q_sel2

    # Build prediction ensemble
    ensemble <- matrix(0, ncol = nrow(X), nrow = K)
    pred_selected <- names(which(colSums(stbm_obj[[assay]]$ensemble[,-1])>0)) # Check whether all predictors that were used to build models are also in the new data
    if( !all( pred_selected %in% colnames(X) )) stop('At least one predictor in the ensemble was not found in the data')

    #Each model makes prediction
    for(i in 1:K){
      pred_selected <- names(which(stbm_obj[[assay]]$ensemble[i,-1]==1))
      ensemble[i,] <- pred_fun(stbm_obj[[assay]]$model[[i]], X[, pred_selected, drop =F])
    }

    aggregated_prediction <- if(assay == 'stable_ensemble') colSums(ensemble*q_sel)/sum(q_sel) else
      colSums(ensemble*q_sel2)/sum(q_sel2)

    attr(aggregated_prediction, which = 'ensemble') <- ensemble
    aggregated_prediction

  }else{
    predict.st2e(stbm_obj[[assay]], X, prune = prune, prefilt_scores = prefilt_scores, pred_fun = pred_fun)

  }

}
