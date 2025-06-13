#' COnvex subspace Shrinkage via Exponential Tilting
#'
#' COSET() will perform importance sampling
#' using an exponential tilted prior that shrinks the
#' parameters towards a convex subspace.
#'
#'
#'
#'
#' @param draws0 named list with names 'posterior', 'prior', and 'logLik' (if using WAIC).
#' Each ought to be a matrix of posterior, prior, or log likelihood draws (each row is a draw),
#' all of which should correspond to the base prior.
#' @param Con_Matrix Definition of convex subspace which describes prior relationships
#' between parameters.
#' @param lb numeric. Lower bound of convex subspace. Can be \eqn{-\infty}.
#' @param ub numeric. Upper bound of convex subspace. Can be \eqn{\infty}.
#' @param CI_level numeric in (0,1). Level of the
#' posterior credible interval.
#' @param model_selection character. One of "bayes_factor" or "waic". Method of picking level of
#' shrinkage via nu.
#' @param nu shrinkage parameter.  If missing, the optimal nu via Bayes factors will be used
#' @param nu_max maximum value of nu to be considered.
#' @param min_ESS integer.  Minimum number of effective sample size in considering values of nu.
#' @param verbose logical. Should any output be provided?
#' @param cl optional object of class c("SOCKcluster", "cluster") generated from parallel::makeCluster.
#'
#'@return Object of class "coset", with the following named elements:
#' \itemize{
#' \item summary, a data frame with the following structure:
#'    \itemize{
#'      \item variable - Variable names extracted from colnames(Sigma0)
#'      \item mean - Posterior mean
#'      \item [(1-CI_level)/2]% - lower CI bound
#'      \item [1 - (1-CI_level)/2]% - upper CI bound
#'    }
#' \item nu the SUBSET shrinkage parameter either fixed by the user or selected by Bayes factor
#' \item ESS the effective sample size
#' \item BF_favoring_nu The Bayes factor in favor of shrinkage vs. no shrinkage (if Bayes Factor was chosen as model selection method)
#' \item WAIC A vector with the following named elements (if WAIC was chosen as model selection method):
#'    \itemize{
#'    \item WAIC_COSET The WAIC of the model using the COSET prior
#'    \item WAIC_base The WAIC of the model using the base prior
#'    }
#' \item proposal_draws Draws input by user (posterior draws under the base prior).
#' \item is_weights Importance sampling weights corresponding to the proposal draws under the
#' SUBSET prior
#' \item base_prior_draws Draws input by the user (prior draws under the base prior).
#' }
#'
#'
#'
#' @import parallel
#' @import graphics
#' @import Matrix
#' @import stats
#' @import utils
#' @import osqp
#' @import methods
#' @export
#' @exportClass coset



COSET= function(draws0,
                Con_Matrix = diag(ncol(draws0[[1]])),
                lb = -rep(Inf,ncol(draws0[[1]])),
                ub = rep(Inf,ncol(draws0[[1]])),
                CI_level = 0.95,
                model_selection = c("bayes_factor", "waic")[1],
                nu,
                nu_max,
                min_ESS = 1/10 * nrow(draws0$posterior),
                verbose = TRUE,
                cl){
  ndraws = sapply(draws0,nrow)
  CI_level = 1 - CI_level
  p = ncol(draws0$posterior)

  model_selection = tolower(model_selection)

  if(verbose){
    if(!model_selection %in% c("waic", "bayes_factor", "bayes factor")){
      if("logLik" %in% names(draws0)){
        model_selection = "waic"
      }else{
        model_selection = "bayes_factor"
        cat("\nNo model selection method chosen. Picking level of shrinkage via Bayes Factor.")
      }
    }
  }

  if(model_selection == "bayes factor") model_selection = "bayes_factor"

  if(verbose)if(missing(nu))cat("\nNo nu provided.  Picking level of shrinkage via Bayes Factor.\n")

  if(verbose){
    if(model_selection == "bayes_factor" && !"logLik" %in% names(draws0)) {
      bayes_factor = TRUE
      cat("\nNo logLikelihood draws provided.  Picking level of shrinkage via Bayes Factor.")
    }
  }

  if(!missing(cl)) if(is.numeric(cl)) cl = makeCluster(min(cl,detectCores() - 1))


  find_dist = function(theta){
    settings <- osqpSettings(verbose = FALSE)
    model <- osqp(P = diag(p), q = -theta, A =  Con_Matrix, l =  lb, u = ub, settings)
    res = model$Solve()
    dist(round(rbind(res$x, theta),5))
  }

  # Find optimal nu
  ## Compute w_k(1)
  if(missing(cl)){
    dist_w =
      apply(draws0$posterior, 1, function(x) {ifelse(sum(Con_Matrix %*% x < lb) > 0, find_dist(x), 0)})
    w_1 =
      exp(-.5 * dist_w)
    dist_tilde_w =
      apply(draws0$prior, 1, function(x) {ifelse(sum(Con_Matrix %*% x < lb) > 0, find_dist(x), 0)})
    wtilde_1 =
      exp(-.5 * dist_tilde_w)
  }else{
    clusterExport(cl,c("draws0","I_m_P"),envir = environment())
    dist_w =
      parApply(cl, draws0$posterior, 1, function(x) {ifelse(sum(Con_Matrix %*% x < lb) > 0, find_dist(x), 0)})
    w_1 =
      exp(-.5 * dist_w)
    dist_tilde_w =
      parApply(cl, draws0$prior, 1, function(x) {ifelse(sum(Con_Matrix %*% x < lb) > 0, find_dist(x), 0)})
    wtilde_1 =
      exp(-.5 * dist_tilde_w)
  }

  if(sum(dist_w) == 0){
    cat("\n No posterior samples outside of subspace. No shrinkage required.")
    return(NULL)
  }

  ## Create function that creates w_k(\nu)
  if(missing(cl)){
    get_w_k = function(nu){
      w_1^nu
    }
    get_wtilde_k = function(nu){
      wtilde_1^nu
    }
  }else{
    clusterExport(cl,c("w_1","wtilde_1"),envir = environment())
    get_w_k = function(nu){
      parSapply(cl,w_1,function(x)x^nu)
    }
    get_wtilde_k = function(nu){
      parSapply(cl,wtilde_1,function(x)x^nu)
    }
  }

  if(missing(nu)){
    ## Find maximum \nu allowed by ESS
    find_max_nu = function(x){
      w_k = get_w_k(x)
      w_k_sum = sum(w_k)
      w_k = w_k / ifelse(w_k_sum == 0,1,w_k_sum)

      ESS = 1 / sum(w_k^2)

      (ESS - min_ESS)^2
    }
    if(missing(nu_max)){
      nu_lower_bound = 0
      nu_max = nu_upper_bound = 5
      safety = 0
      if(verbose) cat("\n---Finding maximum nu to satisfy ESS constraints\n")
      while( ( abs(nu_max - nu_upper_bound) / nu_upper_bound < 1e-3) &
             (safety < 25) ){
        nu_upper_bound = 2 * nu_upper_bound
        nu_max =
          optimize(find_max_nu,
                   interval = c(nu_lower_bound,
                                nu_upper_bound))$min
        nu_lower_bound = nu_upper_bound
        safety = safety + 1
      }
    }

    if(model_selection == "bayes_factor"){
      ## Get optimal nu according to bayes factor
      best_nu_bf = function(nu){
        # Give log( BF(0,nu) )
        log(mean(get_wtilde_k(nu))) -
          log(mean(get_w_k(nu)))
      }
      nu = optimize(best_nu_bf,
                    interval = c(0,nu_max))$min
      nu = ifelse(best_nu_bf(0) <= best_nu_bf(nu),0,nu)
      if(nu == 0) cat("\n 0 Chosen as Best Fitting nu \n")
    }
    else{
      ## Get IS weights
      get_w_k_is = function(nu){
        w_k = get_w_k(nu)
        w_k_sum = sum(w_k)
        w_k = w_k / ifelse(w_k_sum == 0,1,w_k_sum)
        w_k
      }
      ## Get optimal nu according to WAIC
      get_lppd_nu = function(ppd_nu){
        sum(log(colSums(ppd_nu)))
      }

      get_lppd_var_nu = function(ppd_nu){
        sum(colSums((ppd_nu)^2) - (colSums(ppd_nu))^2)
      }

      best_nu_waic = function(nu){
        w_k = get_w_k_is(nu)
        ppd = exp(draws0$logLik)*w_k
        lppd = get_lppd_nu(ppd)
        lppd_var = get_lppd_var_nu(ppd)
        waic = -2*lppd + 2*lppd_var
        waic
      }
      nu = optimize(best_nu_waic,
                    interval = c(0, nu_max))$min
      nu = ifelse(best_nu_waic(0) <= best_nu_waic(nu),0,nu)
      if(nu == 0) cat("\n 0 Chosen as Best Fitting nu \n")
    }
  }



  if(!missing(cl)) clusterExport(cl,"nu",envir=environment())

  if(model_selection == "bayes_factor"){
    # Get final IS weights
    w_k = get_w_k(nu)
    BF_nu_0 =
      exp(
        log(mean(get_w_k(nu))) -
          log(mean(get_wtilde_k(nu)))
      )
    w_k = w_k / sum(w_k)
  }else{
    w_k = get_w_k_is(nu)
    waic_nu = best_nu_waic(nu)
    waic_null = best_nu_waic(0)
  }



  # Get ESS
  ESS = 1 / sum(w_k^2)

  # Get point estimates
  theta_hat =
    apply(draws0$posterior,2,weighted.mean,w = w_k)

  # Get CI
  get_CI_from_IS = function(x,w){
    ord = order(x)
    w = w[ord]
    x = x[ord]
    w_cum = cumsum(w)

    c(lower = x[min(which(w_cum >= CI_level/2))],
      upper = x[min(which(w_cum >= 1 - CI_level/2))]
    )
  }

  # Put results together
  if(is.null(colnames(draws0$posterior))){
    if(!is.null(colnames(draws0$prior))){
      colnames(draws0$posterior) =
        colnames(draws0$prior)
    }else{
      colnames(draws0$posterior) =
        colnames(draws0$prior) =
        paste("theta",1:p,sep="_")
    }
  }
  results = list()
  results$summary =
    data.frame(variable = colnames(draws0$posterior),
               posterior_mean = theta_hat,
               lower = numeric(p),
               upper = numeric(p))
  for(j in 1:p){
    temp =
      get_CI_from_IS(draws0$posterior[,j],w_k)
    results$summary$lower[j] = temp["lower"]
    results$summary$upper[j] = temp["upper"]
  }
  colnames(results$summary)[3:4] =
    c(paste0(CI_level/2 * 100,"%"),
      paste0((1-CI_level/2) * 100,"%"))

  results$nu = nu
  results$ESS = ESS
  if(model_selection == "bayes_factor"){
    results$BF_favoring_nu = BF_nu_0
  }else{
    results$WAIC= c("WAIC_COSET" = waic_nu,
                    "WAIC_base" = waic_null)
  }

  results$proposal_draws = draws0$posterior
  results$is_weights = w_k
  results$base_prior_draws = draws0$prior
  if("logLik" %in% names(draws0)) results$logLik = draws0$logLik

  class(results) = "coset"

  rownames(results$summary) = NULL
  return(results)
}

