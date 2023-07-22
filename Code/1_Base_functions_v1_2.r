## Base functions for ML-SAC manuscript
## Written by Insang Song (sigmafelix@hotmail.com)
## Last update: 06/16/2023: data splitting


if (!require(pacman)) { 
    install.packages("pacman")
    library(pacman) }
p_load(tidyverse, spdep, vegan, sf, sp, stringr, readxl, kernlab, ranger, foreach, vip, doParallel, CAST, mlr, parallelMap, doSNOW, corrplot, devtools)

# landmap package for RFSp is in a GitHub repository
devtools::install_github("envirometrix/landmap")


## Utility functions
replace_nan <- function(x) ifelse(is.nan(x), 0, x)
replace_na2 <- function(x) ifelse(is.na(x), 0, x)
replace_inf <- function(x) ifelse(is.infinite(x), 0, x)
rmse_direct = function(y, yhat) sqrt(mean((y-yhat)^2))

gen_form <- function(form.origin, vectors){
    formo <- as.character(form.origin)
    vcn <- colnames(vectors)
    form <- as.formula(paste(formo[2], formo[1], paste(formo[3]), '+', paste(vcn, collapse = '+')))
    return(form)
}

gen_dat <- function(dat.origin, me){
    dat <- bind_cols(dat.origin, me$vectors %>% data.frame)
    return(dat)
}
gen_dat.sf <- function(dat.origin, sf){
  dat <- bind_cols(dat.origin, fitted(sf) %>% data.frame)
  return(dat)
}



### varpart
ME.varpart <- function(dat, xvars = imix.xvec){
  # dat: x by y data.frame. which contains ME vectors and X variables
  library(vegan)
  lc <- length(colnames(dat))
  datcols <- colnames(dat)
  form.xs <- as.formula(paste('~', paste(xvars[-1], collapse = '+')))
  form.me <- as.formula(paste('~', paste(datcols[grep('^vec.', datcols)], collapse = '+')))
  #full_df <- ;
  MEV <- varpart(dat[,lc], form.xs, form.me, data = dat)
  
  return(MEV)
}

treat_analysis <- function(treat,
                           yvec,
                           xvec,
                           mode = 'area'){
  treat.st <- treat %>% as('Spatial')
  if (mode == 'area'){
    treat.nb <- treat %>% as('Spatial') %>% poly2nb
  }
  if (mode == 'point'){
    treat.nb <- treat %>% as('Spatial') %>% knearneigh(k = floor(nrow(treat) * 0.1)) %>% knn2nb
  }
  treat.lw <- treat.nb %>% nb2listw(zero.policy = TRUE)
  
  
  treat.yvec <- sort(yvec)
  treat.xvec <- xvec
  treat.mi <- treat.yvec %>% split(.,.) %>% 
    lapply(function(x) moran.test(unlist(treat.st@data[,x]), listw = treat.lw, zero.policy = TRUE)$estimate[1]
    )
  treat.mi <- treat.mi %>% do.call(c, .)
  
  treat.forms <- treat.yvec %>% split(.,.) %>% lapply(function(x) as.formula(paste(x, '~', paste(treat.xvec, collapse = '+'))))
  treat.lms <- treat.forms %>% lapply(function(x) lm(formula = x, data = treat.st@data))
  
  treat.sf <- treat.forms %>% 
    lapply(function(x) SpatialFiltering(formula = x, nb = treat.nb, ExactEV = FALSE, data = treat.st@data, zero.policy = TRUE))
  # augmented data
  treat.sfa <- lapply(treat.sf, function(x) return(gen_dat.sf(treat.st@data, x)))
  treat.form1 <- mapply(function(x, y) gen_form(x, y$dataset), treat.forms, treat.sf, SIMPLIFY = FALSE)
  treat.sfl <- mapply(function(x, y) lm(x, data = y), treat.form1, treat.sfa, SIMPLIFY = FALSE)
  
  # rsq1
  treat.sfl.rsq <- treat.sfl %>% lapply(function(x) summary(x)$adj.r.squared)
  treat.lms.rsq <- treat.lms %>% lapply(function(x) summary(x)$adj.r.squared)
  # ressac
  treat.sf.rs <- treat.sfl %>% 
    lapply(residuals) %>% 
    lapply(function(x) moran.test(x = x, listw = treat.lw, zero.policy = TRUE)$estimate[1])
  treat.lm.rs <- treat.lms %>% 
    lapply(residuals) %>% 
    lapply(function(x) moran.test(x = x, listw = treat.lw, zero.policy = TRUE)$estimate[1])
  
  
  treat.result <- data.frame(
    MoranI = treat.mi,
    LMresI = treat.lm.rs %>% do.call(c, .),
    SFresI = treat.sf.rs %>% do.call(c, .),
    LMrsq = treat.lms.rsq %>% do.call(c, .),
    SFrsq = treat.sfl.rsq %>% do.call(c, .)
  )
  return(treat.result)
}


## make listw object (compatible with spdep::listw)
get_nblw = function(spdata, mode = 'area', pp = 0.1, maxval = 30, gdbpath = NULL) {
  if (mode == 'area'){
    treat.nb <- spdata %>% poly2nb
  } else if (mode == 'point'){
    treat.nb <- spdata %>%  
      knearneigh(k = floor(pp*nrow(.))) %>% knn2nb
  } else if (mode == 'network'){
    if (is.null(gdbpath)) 
    {  stop('No input data. Please check your network or sites.')  }
    if (any(class(gdbpath) == 'matrix'))
    {
      treat.nb <- gdbpath
      treat.nb <- treat.nb %>% 
        st_distmat_to_swm(., threshold) %>% 
        spdep::mat2listw(.) %>% 
        .$neighbours
    } else {
      treat.nb <- gdbpath %>% 
        spdep::mat2listw(.) %>% 
        .$neighbours
    }
  } else {
    stop("mode should be one of 'area', 'point', and 'network.'")
  }
  return(treat.nb)
}



# modification: 200623 (finished)
# 200919: sorting problem in split(.,.) part: fixed - may affect the whole result
# 220209: vip::vi permutation based variable importance for both models
# Spatial filtering analysis
sf_analysis <- function(treat,
                        yvec,
                        xvec,
                        mode = 'area',
                        std = TRUE,
                        return.residuals = FALSE,
                        extract.coef = FALSE,
                        exact = FALSE,
                        maxval = 30, # obsolete afterwards
                        pp = 0.1,
                        gdbpath = NULL,
                        threshold = 1e4){
  treat.st <- treat #%>% as('Spatial')
  if (std){
    treat.st <- treat.st %>% 
      mutate_if(is.numeric, list(~as.vector(scale(.))))
  }

  treat.nb = get_nblw(treat.st, mode = mode, maxval = maxval, pp = pp)
  treat.lw <- treat.nb %>% nb2listw(zero.policy = TRUE)
  treat.st = st_drop_geometry(treat.st)
  
  treat.yvec <- sort(yvec)
  treat.xvec <- xvec
  treat.mi <- treat.yvec %>% split(.,.) %>% 
    lapply(function(x) moran.test(unlist(treat.st[,x]), listw = treat.lw, zero.policy = TRUE)$estimate[1]
    )
  treat.mi <- treat.mi %>% do.call(c, .)
  
  treat.forms <- treat.yvec %>% split(.,.) %>% lapply(function(x) as.formula(paste(x, '~', paste(treat.xvec, collapse = '+'))))
  treat.lms <- treat.forms %>% lapply(function(x) lm(formula = x, data = treat.st, y = TRUE))
  
  treat.sf <- treat.forms %>% 
    lapply(function(x) SpatialFiltering(formula = x, nb = treat.nb, ExactEV = exact, data = treat.st, zero.policy = TRUE))
  # augmented data
  treat.sfa <- lapply(treat.sf, function(x) return(gen_dat.sf(treat.st, x)))
  treat.form1 <- mapply(function(x, y) gen_form(x, y$dataset), treat.forms, treat.sf, SIMPLIFY = FALSE)
  treat.sfl <- mapply(function(x, y) lm(x, data = y, y = TRUE), treat.form1, treat.sfa, SIMPLIFY = FALSE)
  
  # rsq1
  treat.sfl.rsq <- treat.sfl %>% lapply(function(x) summary(x)$adj.r.squared)
  treat.lms.rsq <- treat.lms %>% lapply(function(x) summary(x)$adj.r.squared)
  treat.sfl.prsq <- treat.sfl %>% lapply(function(x) summary(x)$r.squared)
  treat.lms.prsq <- treat.lms %>% lapply(function(x) summary(x)$r.squared)
  # ressac
  treat.sf.rs <- treat.sfl %>% 
    lapply(residuals)
  treat.sf.rsi <- treat.sf.rs %>% 
    lapply(function(x) moran.test(x = x, listw = treat.lw, zero.policy = TRUE)$estimate[1])
  treat.lm.rs <- treat.lms %>% 
    lapply(residuals)
  treat.lm.rsi <- treat.lm.rs %>% 
    lapply(function(x) moran.test(x = x, listw = treat.lw, zero.policy = TRUE)$estimate[1])
  
  # accuracy measures (200623)
  treat.sf.am <- treat.sfl %>% 
    lapply(function(x) 
      data.frame(RMSE = mlr::measureRMSE(x$y, x$y + residuals(x)),
        NRMSE = mlr::measureRMSE(x$y, x$y + residuals(x))/(max(x$y)-min(x$y)),
        MAE = mlr::measureMAE(x$y, x$y + residuals(x)),
        MAPE = mlr::measureMAPE(x$y, x$y + residuals(x)))) 
  treat.lm.am <- treat.lms %>% 
    lapply(function(x) 
      data.frame(RMSE = mlr::measureRMSE(x$y, x$y + residuals(x)),
                 MAE = mlr::measureMAE(x$y, x$y + residuals(x)),
                 MAPE = mlr::measureMAPE(x$y, x$y + residuals(x))))

  
  # varimp (relative rank of regression coefficients)  
  treat.sfl.vi <- treat.sfl %>% lapply(function(x) x[1]$coefficients[2:(length(treat.xvec)+1)] %>% t %>% as.data.frame)
  treat.lms.vi <- treat.lms %>% lapply(function(x) x[1]$coefficients[2:(length(treat.xvec)+1)] %>% t %>% as.data.frame)
  varimp.result <- mapply(FUN = function(x, y){
    vis <- bind_rows(x, y) %>% 
      mutate(Model = c('SF', 'OLS')) %>% 
      pivot_longer(names_to = 'independent', values_to = 'coef', cols = 1:length(treat.xvec))
    return(vis)},
    treat.sfl.vi, treat.lms.vi, SIMPLIFY = FALSE) %>% 
    mapply(function(x, y){
      x1 <- x %>% 
        mutate(dependent = y)
      return(x1)
    }, ., treat.yvec %>% split(.,.), SIMPLIFY = FALSE) %>% 
    do.call(bind_rows, .) %>% 
    group_by(Model, dependent) %>% 
    mutate(Rank = n() + 1 - rank(abs(coef))) %>% 
    ungroup
  
  # just in case: using permutation based vi
  varimp.result = 
    mapply(function(res_lm, res_sf, y) {
        extract_vi = function(mod, yvar, metric = 'rmse') {
            vip::vi(mod, train = mod$model,
                    nsim = prod(length(xvec), length(xvec) -1 ),
                    target = yvar, pred_wrapper = predict, 
                    method = 'permute', metric = c('rmse'))
        }
        vi_lm = extract_vi(res_lm, y) %>%
            mutate(Model = 'OLS')
        vi_sf = extract_vi(res_sf, y) %>%
            filter(Variable %in% xvec) %>%
            mutate(Model = 'SF')
        vi_lmsf = bind_rows(vi_lm, vi_sf) %>%
            rename(coef = Importance,
                   independent = Variable) %>%
            mutate(dependent = y) %>%
            group_by(Model, dependent) %>%
            mutate(Rank = rank(-abs(coef)))
        return(vi_lmsf)
    }, treat.lms, treat.sfl, treat.yvec %>% split(.,.), SIMPLIFY = FALSE) %>%
    do.call(bind_rows, .)

  treat.result1 <- data.frame(
    dependent = treat.yvec,
    MoranI_Y = treat.mi,
    Model = 'OLS',
    MoranI_resid = treat.lm.rsi %>% do.call(c, .),
    Rsquared_adj = treat.lms.rsq %>% do.call(c, .),
    Rsquared_plain = treat.lms.prsq %>% do.call(c,.)
  ) %>% 
    bind_cols(do.call(rbind, treat.lm.am))
  treat.result2 <- data.frame(
    dependent = treat.yvec,
    MoranI_Y = treat.mi,
    Model = 'SF',
    MoranI_resid = treat.sf.rsi %>% do.call(c, .),
    Rsquared_adj = treat.sfl.rsq %>% do.call(c, .),
    Rsquared_plain = treat.sfl.prsq %>% do.call(c,.)
  ) %>% 
    bind_cols(do.call(rbind, treat.sf.am))
  treat.result <-
    bind_rows(treat.result1, treat.result2)
  if (extract.coef){
    return(varimp.result)
  } else {
    if (return.residuals){
      return(list(lmsfres = treat.result, lmres = treat.lm.rs, sfres = treat.sf.rs))
    } else {
      return(treat.result)
    }
  }
}







X.sac <- function(treat, xvars, mode){
  #treat.st <- treat %>% as('Spatial')
  if (mode == 'area'){
    treat.nb <- treat %>% as('Spatial') %>% poly2nb
  }
  if (mode == 'point'){
    treat.nb <- treat %>% as('Spatial') %>% knearneigh(k = floor(nrow(treat) * 0.1)) %>% knn2nb
  }
  treat.lw <- treat.nb %>% nb2listw(zero.policy = TRUE)
  
  xl <- xvars %>% split(.,.)
  xls.mt <- xl %>% 
    lapply(function(x) moran.test(x = treat[,x] %>% st_set_geometry(NULL) %>% unlist, 
                                  listw = treat.lw, randomisation = F,
                                  zero.policy = TRUE)$estimate[1])
  xls.mt <- xls.mt %>% do.call(cbind, .) %>% data.frame
  colnames(xls.mt) <- xvars
  return(xls.mt)
    
}

### residual SAC
res.sac <- function(dat, xvars = imix.xvec, nbs = imix.nb[[1]]){
  # dat: x by y data.frame. which contains ME vectors and X variables
  lc <- length(colnames(dat))
  datcols <- colnames(dat)
  yvec <- datcols[lc]
  form.full <- as.formula(paste(yvec, '~', paste(datcols[-c(1, lc)], collapse = '+')))
  form.o <- as.formula(paste(yvec, '~', paste(datcols[-c(grep('^(vec).', datcols), lc)], collapse = '+')))
  mod.full <- lm(formula = form.full, data = dat)
  mod.o <- lm(formula = form.o, data = dat)
  aic.full <- AIC(mod.full)
  aic.o <- AIC(mod.o)
  
  res.full <- residuals(mod.full)
  res.o <- residuals(mod.o)
  mor.full <- moran.test(res.full, listw = nbs)$estimate[1]
  mor.o <- moran.test(res.o, listw = nbs)$estimate[1]
  
  mor.diff <- list(mor.full, mor.o)
  names(mor.diff) <- c('MoranI.full', 'MoranI.O')
  lap <- mor.diff
  return(lap)
}

### r.squared
rsqs <- function(dat, xvars = imix.xvec, nbs = imix.nb[[1]]){
  # dat: x by y data.frame. which contains ME vectors and X variables
  lc <- length(colnames(dat))
  datcols <- colnames(dat)
  yvec <- datcols[lc]
  form.full <- as.formula(paste(yvec, '~', paste(datcols[-c(1, lc)], collapse = '+')))
  form.o <- as.formula(paste(yvec, '~', paste(datcols[-c(grep('^(vec).', datcols), lc)], collapse = '+')))

  mod.full <- lm(formula = form.full, data = dat)
  mod.o <- lm(formula = form.o, data = dat)
  rsq.full <- summary(mod.full)$adj.r.squared
  rsq.o <- summary(mod.o)$adj.r.squared
  
  mor.diff <- list(rsq.full, rsq.o)
  names(mor.diff) <- c('rsq_full', 'rsq_O')
  lap <- mor.diff
  return(lap)
}



### AIC computation
aic.diff <- function(dat, xvars = imix.xvec){
  # dat: x by y data.frame. which contains ME vectors and X variables
  lc <- length(colnames(dat))
  datcols <- colnames(dat)
  yvec <- datcols[lc]
  form.full <- as.formula(paste(yvec, '~', paste(datcols[-c(1, lc)], collapse = '+')))
  form.o <- as.formula(paste(yvec, '~', paste(datcols[-c(grep('^(vec).', datcols), lc)], collapse = '+')))
  
  mod.full <- lm(formula = form.full, data = dat)
  mod.o <- lm(formula = form.o, data = dat)
  
  aic.full <- AIC(mod.full)
  aic.o <- AIC(mod.o)
  
  aic.diff <- list(aic.full, aic.o)#aic.o - aic.full
  lap <- aic.diff
  return(lap)
}



### MLR_LEARN ####
# 200708: added deep learning (by h2o)
# 210527: fixed set.seed (mode="L'Ecuyer")
# 210529: ranger added; resampling strategy from cv to nested cv
# 220208: ksvm added; 
mlr_learn <- function(input_data, algs = list(), yvec, xvec, 
                      cvmethod = 'CV', ncore = floor(parallel::detectCores()/2.5), 
                      std = FALSE, ncv = 10, seed.num = 202006){
  #library(mlrCPO)
  if (length(yvec) != 1){
    cat('The task cannot run with yvec longer than 1.')
    break
  }
  set.seed(seed.num, "L'Ecuyer")
  cat(sprintf('\nModeling %s against [%s] with %d-fold cross-validation......\n\n', yvec, str_c(xvec, collapse = ', '), ncv))
  input_std <- input_data %>% 
    st_set_geometry(NULL) %>% 
    dplyr::select(yvec, xvec) #%>% 
  if (std){
    input_std <- input_std %>% 
      mutate_if(is.numeric, list(~scale(.) %>% as.vector))
  } 
    #
  mlr_learn <- algs %>% split(.,.) %>% 
    lapply(function(x){
      library(parallelMap)
      parallelStartSocket(ncore)
      cat(sprintf("Running %s now...\n", x))
      # Make learner and define the regression task
      lx <- makeLearner(cl = x, predict.type = 'response')
      lre <- makeRegrTask(data = input_std, 
                          target = yvec,
                          coordinates = st_centroid(input_data, of_largest_polygon = TRUE) %>% st_coordinates %>% data.frame) 
      
      ## 06172023: split (2:1 or 3:1)
      ## lre_resamp -> lre_train, lre_task
      lre_resamp = makeResampleInstance("Holdout", lre)
      lre_train = subsetTask(lre, lre_resamp$train.inds[[1]])
      lre_test = subsetTask(lre, lre_resamp$test.inds[[1]])


      # parset
      if (grepl('*.(lm)$', x)){
        pset <- makeParamSet(
          makeIntegerParam('tol', lower = -10, upper = -5, trafo = function(x) 10 ^ x) 
        )
        ctrl <- makeTuneControlGrid()
        
      }
      if (grepl('*.(xgboost)$', x)) {
        pset = makeParamSet(
          makeDiscreteParam('subsample', values = seq(0.1, 0.5, 0.1)),
          makeDiscreteParam('normalize_type', values = c('tree', 'forest')),
          makeDiscreteParam('grow_policy', values = c('depthwise', 'lossguide'))
        )
        ctrl <- makeTuneControlGrid()
      }
      if (grepl('*.(deeplearning)$', x)){
        pset <- makeParamSet(
          makeDiscreteParam('activation', values = c('RectifierWithDropout')),
          makeDiscreteParam('epochs', values = c(50)),
          makeDiscreteParam('hidden', values = list(set1 = c(50, 50), set2 = c(100, 100))),
          makeDiscreteParam(id = 'hidden_dropout_ratios', values = list(vec1 = c(0.25, 0.25), vec2 = c(0.5, 0.5), vec3 = c(0.75, 0.75)))
        )
        ctrl <- makeTuneControlGrid()
      }
      if (grepl('*.(svm)$', x)){
        pset <- makeParamSet(
          makeIntegerParam("nu", lower = -6, upper = -1, trafo = function(x) 2 ** x),
          makeIntegerParam("cost", lower = -4, upper = 3, trafo = function(x) 10 ^ x),
          makeDiscreteParam("kernel", values = c("polynomial", "radial basis", "linear")),
          makeIntegerParam("gamma", lower = -5, upper = 0, trafo = function(x) 2^x,
                           requires = quote(kernel == "radial basis")),
          makeIntegerParam("degree", lower = 2L, upper = 4L,
                           requires = quote(kernel == "polynomial"))
        )
        
        ctrl <- makeTuneControlGrid()
        
      } 
      if (grepl('*.(ksvm)$', x)) {
        pset <- makeParamSet(
          makeDiscreteParam("type", values = c("nu-svr")),
          makeDiscreteParam("kernel", values = c("polydot", "rbfdot", "vanilladot")),
          #makeIntegerParam("nu", lower = -6, upper = -1, trafo = function(x) 2 ** x),
          #makeIntegerParam("C", lower = -4, upper = 3, trafo = function(x) 10 ^ x),
          makeIntegerParam("sigma", lower = -5, upper = 1, trafo = function(x) 2^x,
                            requires = quote(kernel == "rbfdot")),
          makeIntegerParam("degree", lower = 2L, upper = 4L,
                            requires = quote(kernel == "polydot"))
        )
        ctrl <- makeTuneControlGrid()
        
      } else if (grepl('*.(ranger)$', x)){
        pset <- makeParamSet(
          makeIntegerParam("mtry", lower = floor(length(xvec)/3), upper = floor(length(xvec) - 1)),
          makeDiscreteParam("num.trees", values = c(300, 500, 1000))#,
          #makeDiscreteParam('sample.fraction', values = c(0.5, 0.6, 0.7, 0.8))
        )
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(nnet)$', x)) {
        pset <- makeParamSet(
          makeIntegerParam("size", lower = 2, 
                           upper = ifelse(floor(nrow(input_std)/30) <= 3, 10, 
                                          min(floor(nrow(input_std)/30), 16)))
        )
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(SVR)$', x)){
        pset <- makeParamSet(
          makeIntegerParam('svr_eps', lower = -4, upper = 4, trafo = function(x) 10^x),
          makeIntegerParam('type', lower = 11, upper = 12)
        )
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(brnn)$', x)) {
        pset <- makeParamSet(
          makeIntegerParam('neurons', lower = 2, upper = max(7, floor(nrow(input_std)^(1/3))))
          )
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(RRF)$', x)){
        pset <- makeParamSet(
          makeIntegerParam('mtry', lower = floor(length(xvec)/2), upper = length(xvec)),
          makeIntegerParam("maxnodes", lower = floor(nrow(input_std)/5), upper = min(floor(nrow(input_std)*0.75), 1000))
          )
        ctrl <- makeTuneControlGrid()
      }

      if (grepl('^Sp', cvmethod)){
          cvinner = 'SpRepCV'
          cvouter = 'SpCV'
          lrs_inner = makeResampleDesc(method = cvinner, 
                              iters = 2, 
                              predict = 'both')

          lrs_outer = makeResampleDesc(method = cvouter,
                                  iters = 5,
                                  predict = 'both') 
          # Main tuning process
          innerwrapper = makeTuneWrapper(learner = lx,
                                        #task = lre,
                                        par.set = pset,
                                        resampling = lrs_inner,
                                        control = ctrl,
                                        show.info = TRUE,
                                        measures = list(mlr::rmse, mlr::mae, mlr::mape, mlr::spearmanrho))
                                        #expvar, spearmanrho
          outerwrapper = tuneParams(innerwrapper, task = lre,
                                  resampling = lrs_outer,
                                  control = ctrl,
                                  par.set = pset,
                                  show.info = TRUE,
                                  measures = list(mlr::rmse, mlr::mae, mlr::mape, mlr::spearmanrho))
      } else {
          cvinner = 'Subsample'
          cvouter = 'CV'
          lrs_outer = makeResampleDesc(method = cvouter,
                                  iters = ncv,
                                  predict = 'both') 
        # task=lre to lre_train (06172023)
          outerwrapper = tuneParams(learner = lx, task = lre_train,
                                  resampling = lrs_outer,
                                  control = ctrl,
                                  par.set = pset,
                                  show.info = TRUE,
                                  measures = list(mlr::rmse, mlr::mae, mlr::mape, mlr::spearmanrho))
      }

      #res <- tuneParams(learner = lx, task = lre, 
      #                  resampling= lrs, 
      #                  control = ctrl,
      #                  par.set = pset,
      #                  measures = list(rmse, mae, mape, expvar, spearmanrho))

      hyper_optim = setHyperPars(learner = lx,
                                   par.vals = outerwrapper$x)
      # task=lre to lre_train (06172023)
      trained_optim = mlr::train(learner = hyper_optim, task = lre_train)
      
      # entire data prediction
      lres = predict(trained_optim, task = lre)
      lres = lres$data$truth - lres$data$response
      #getTaskData(lres)[,yvec] %>% unlist %>% as.vector

      if (grepl('*.(deeplearning)$', x)){
        vimp <- NA
      } else {
        vimp <- generateFeatureImportanceData(lre, learner = lx, nmc=100L)
        #vimp = getFeatureImportance(trained_optim)
      }
      # Get errors (residuals) from the optimal model

      #lres = lre %>>% cpoRegrResiduals(learner = x, 
      #                                 crr.train.residuals = 'resample')
      parallelStop()
      # res to outerwrapper
      res <- list(tuneResult = outerwrapper, residuals = lres, varimp = vimp, y = unlist(input_std[,yvec]))
      return(res)
    }
    )
  return(mlr_learn)
  }

mlr_learn_tr <- function(input_data, algs = list(), yvec, xvec, 
                      ncore = parallel::detectCores()-2, std = FALSE, ncv = 10,
                      seed.num = 202006){
  library(mlrCPO)
  if (length(yvec) != 1){
    cat('The task cannot run with yvec longer than 1.')
    break
  }
  
  set.seed(seed.num, "L'Ecuyer")
  cat(sprintf('Modeling %s against [%s] with the entire data...\n', yvec, str_c(xvec, collapse = ', ')))
  input_std <- input_data %>% 
    st_set_geometry(NULL) %>% 
    dplyr::select(yvec, xvec) #%>% 
  if (std){
    input_std <- input_std %>% 
      mutate_if(is.numeric, list(~scale(.) %>% as.vector))
  } 
  #
  mlr_learn <- algs %>% split(.,.) %>% 
    lapply(function(x){
      library(parallelMap)
      cl <- parallelStartSocket(ncore)
      registerDoSNOW(cl)
      lx <- makeLearner(cl = x, predict.type = 'response')
      lre <- makeRegrTask(data = input_std, 
                          target = yvec,
                          coordinates = st_centroid(input_data, of_largest_polygon = TRUE) %>% st_coordinates %>% data.frame) 

      if (grepl('*.(lm)$', x)){
        pset <- makeParamSet(
          makeIntegerParam('tol', lower = -10, upper = -5, trafo = function(x) 10 ^ x) 
        )
        ctrl <- makeTuneControlGrid()
        
      }
      if (grepl('*.(svm)$', x)){
        pset <- makeParamSet(
          makeIntegerParam("nu", lower = -4, upper = 4, trafo = function(x) 10 ** x),
          makeDiscreteParam("kernel", values = c("polynomial", "radial basis", "linear")),
          makeIntegerParam("gamma", lower = -6, upper = 6, trafo = function(x) 2^x,
                           requires = quote(kernel == "radial basis")),
          makeIntegerParam("degree", lower = 2L, upper = 4L,
                           requires = quote(kernel == "polynomial"))
        )
        ctrl <- makeTuneControlGrid()#makeTuneControlRandom(maxit = 100L)
        
      } else if (grepl('*.(ranger)$', x)){
        pset <- makeParamSet(
          makeIntegerParam("mtry", lower = floor(length(xvec)/2), upper = floor(length(xvec))),
          makeIntegerVectorParam("num.trees", lower = 300, upper = 1000, len = 8)#,
          )#,
          #makeDiscreteParam('sample.fraction', values = c(0.5, 0.6, 0.7, 0.8))
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(nnet)$', x)) {
        pset <- makeParamSet(
          makeIntegerParam("size", lower = 2, 
                           upper = ifelse(floor(nrow(input_std)/30) <= 3, 10, 
                                          min(floor(nrow(input_std)/30), 16)))
        )
        ctrl <- makeTuneControlGrid()#makeTuneControlRandom(maxit = 30L)
        #ctrl <- makeTuneControlGrid(resolution = 6L)
      } else if (grepl('*.(SVR)$', x)){
        pset <- makeParamSet(
          makeIntegerParam('svr_eps', lower = -4, upper = 4, trafo = function(x) 10^x),
          makeIntegerParam('type', lower = 11, upper = 12)
        )
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(brnn)$', x)) {
        pset <- makeParamSet(
          makeIntegerParam('neurons', lower = 2, upper = max(7, floor(nrow(input_std)^(1/3))))
        )
        ctrl <- makeTuneControlGrid()
      } else if (grepl('*.(RRF)$', x)){
        pset <- makeParamSet(
          makeIntegerParam('mtry', lower = floor(length(xvec)/2), upper = length(xvec)),
          makeIntegerParam("maxnodes", lower = floor(nrow(input_std)/5), upper = min(floor(nrow(input_std)*0.75), 1000))
        )
        ctrl <- makeTuneControlGrid()
      }
      res <- 'notrain'#train(learner = lx, task = lre)
      #res <- tuneParams(learner = lx, task = lre, 
      #                  resampling= lrs, #makeResampleDesc('CV', predict = 'both', iters = 10),#lrs, 
      #                  control = ctrl,
      #                  par.set = pset,
      #                  #resample.fun = rmse,
      #                  measures = list(rmse, mae, mape, expvar, spearmanrho))
      vimp <- generateFeatureImportanceData(lre, learner = lx, nmc=100L)
      lres = lre %>>% cpoRegrResiduals(learner = x)
      lres = getTaskData(lres)[,yvec] %>% unlist %>% as.vector
      parallelStop()
      res <- list(tuneResult = res, y = input_std[,yvec] %>% unlist,
                  residuals = lres, varimp = vimp)
      return(res)
    }
    )
  return(mlr_learn)
}


### Extract Model Metric #### 
## Function : triple layered lists
extr_modelres <- function(mm, yvec, modellist = calg1, neighbor, mode = 1){
  if (mode != 1){
    mm.vals <- mm %>% 
      lapply(function(x){ x %>% 
          lapply(function(y){ 
            ym <- moran.test(x = y$residuals, listw = nb2listw(neighbor, zero.policy = TRUE), zero.policy = TRUE)$estimate[1]
            ycon <- c(y$tuneResult$y[1:3], ym)
            return(ycon)}) %>% 
            do.call(c, .) %>% 
            data.frame(measure = c('RMSE', 'MAE', 'MAPE', 'res.MI'),
                        value = .) -> df.cleaned
            return(df.cleaned)
        }) %>% 
      do.call(rbind, .) %>% 
      cbind(Model = rownames(.), .) %>% 
      mutate(YVAR = str_remove(Model, '.regr.(lm|nnet|randomForest|ranger|svm).*'),
            Model = str_extract(Model, '(lm|nnet|randomForest|ranger|svm)'))
  } else {
    mm.vals <- mm %>% 
      lapply(function(x) x %>% lapply(function(y){ 
              yt <- y$y
              ye <- y$res + y$y
              ymi <- moran.test(ye, nb2listw(neighbor, zero.policy = TRUE), zero.policy = TRUE)$estimate[1]
              ym <- c(rmse(yt, ye),
                      rmse(yt, ye)/(max(yt)-min(yt)),
                      mae(yt, ye),
                      mape(yt, ye),
                      ymi)
              return(ym)}) %>% 
               do.call(c, .) %>% 
               data.frame(measure = c('RMSE', 'NRMSE', 'MAE', 'MAPE', 'res.MI'),
                          value = .)) %>% 
      do.call(rbind, .) %>% 
      cbind(Model = rownames(.), .) %>% 
      mutate(YVAR = str_remove(Model, '.regr.(lm|nnet|randomForest|ranger|svm).*'),
             Model = str_extract(Model, '(lm|nnet|randomForest|ranger|svm)'))
    
  }    
  return(mm.vals)
}



### POST_HOC MI extraction
posthoc_mi <- function(mi, sfd, yvec, xvec, 
                         mode = 'point', maxval = 30, pp = 0.1,
                         network_wm = NULL, threshold = 5e3){
  sfd <- sfd %>% mutate_at(.vars = vars(-ncol(.)),
                           .funs = list(~as.vector(scale(.))))
  if (mode != 'network'){
    sfd.nw <- 
      nb2listw(spdep::knearneigh(sfd, 
                                 floor(pp*nrow(sfd))) %>% knn2nb,
               zero.policy = TRUE)
  } else {
    if (is.null(network_wm)) stop('Error: no network distance matrix was specified.')
    if (any(class(network_wm) == 'matrix'))
      {
      sfd.nw <- network_wm %>% 
        st_distmat_to_swm(., threshold) %>% 
        spdep::mat2listw(.)
    } else {
      sfd.nw <- network_wm %>% 
        spdep::mat2listw(.)
    }
    }
  ## 1. Y-vector moran's I
  yvec = sort(yvec)
  yvec.mi <- vector('list', length = length(yvec))
  for (i in 1:length(yvec)){
    yvec.mi[[i]] <-  
      moran.test(sfd[,yvec[i]] %>% st_set_geometry(NULL) %>% unlist, 
                 listw = sfd.nw, zero.policy = TRUE)$estimate[1]
  }
  ## 2. Extract results
  # input structure
  #   - yvar[i] each = nalgs
  #     - algs[j] rep_n = nyvar
  #       - contents
  # strategy: flatten all the nested structures up to the second level
  # e.g., [[1]][[1]] [[1]][[2]] [[2]][[1]] [[2]][[2]]
  # to [[1]] [[2]] [[3]] [[4]]
  nalgs = length(mi[[1]])
  res_ys = yvec %>% split(.,.) %>%
    lapply(function(x) rep(list(sfd[,x] %>% st_drop_geometry %>% unlist), nalgs)) %>%
    do.call(c, .)
  res_resid <- lapply(mi, function(x) lapply(x, function (y) y$residuals)) %>% 
    do.call(c, .)
    #lapply(function(x) do.call(list, x))
  res_rmses = mapply(function(y, e) { rmse_direct(y, y + e)}, res_ys, res_resid, SIMPLIFY = FALSE) %>%
    Reduce(c, .)
  res_ys_mi = res_ys %>%
    lapply(function(x) moran.test(x, listw = sfd.nw, zero.policy = TRUE)$estimate[1]) %>%
    Reduce(c, .)
  res_resid_mi = res_resid %>%
    lapply(function(x) moran.test(x, listw = sfd.nw, zero.policy = TRUE)$estimate[1]) %>%
    Reduce(c, .)
  names_algs = str_replace_all(names(res_resid), '\\d{1,2}.regr.', '')  
  data.frame(dependent = rep(yvec, each = nalgs),
             Model = rep(names_algs, length(yvec)),
             RMSE = res_rmses,
             MoranI_Y = res_ys_mi,
             MoranI_resid = res_resid_mi) -> res_total

  
  model_nums = length(names(mi[[1]]))
  res_varimp <- lapply(mi, function(x) lapply(x, function (y) rank(-y[[3]]$res))) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) x %>% data.frame %>%
             cbind(rownames(.), .)) %>% 
    do.call(rbind, .) %>% 
    data.frame()
  colnames(res_varimp)[1] <- 'Model'
  res_varimp <- res_varimp %>% 
    cbind(., dependent = rep(yvec, each = model_nums)) %>% 
    pivot_longer(cols = 2:(ncol(.)-1), names_to = 'independent', values_to = 'Rank')
  
#   res_resid <- res_resid %>% 
#     lapply(function(x) 
#       lapply(x, function(y) moran.test(y, sfd.nw, zero.policy = TRUE)$estimate[1]) %>% 
#         do.call(bind_rows, .) %>% 
#         data.frame %>% 
#         cbind(rownames(.), .)
#     ) %>% 
#     do.call(bind_rows, .)
#   colnames(res_resid)[1] <- 'Resid_Y'

#   res_total <- bind_cols(res_result, res_resid) %>% 
#     mutate(dependent = rep(yvec, each = model_nums), # yvec was not sorted (020622)
#            MoranI_Y = yvec.mi %>% do.call(c, .) %>% rep(., each = model_nums)) %>% 
#     dplyr::select(-Resid_Y)
#   colnames(res_total) <- c('Model', 'RMSE', 'MoranI_resid', 'dependent', 'MoranI_Y')
  
  ## 3. sf_analysis
  res_sf <- sf_analysis(sfd, yvec, xvec, mode = 'point', maxval = maxval, pp = pp, extract.coef = FALSE)
  res_sfc <- sf_analysis(sfd, yvec, xvec, mode = 'point', maxval = maxval, pp = pp, extract.coef = TRUE)
  
  res_total <- res_sf %>% bind_rows(res_total)
  res_varimp <- res_sfc %>% bind_rows(res_varimp)
  
  posthoc <- list(results_combined = res_total, results_varimp = res_varimp, results_sf = res_sf, results_sf_rank = res_sfc)
  
  return(posthoc)
}



### XSAC


Xmc <- function(dat, xvec){
  datx <- dat %>% st_set_geometry(NULL) %>% .[,xvec]
  datcorr <- cor(datx)
  datcorr <- datcorr[upper.tri(datcorr)]
  corr <- mean(datcorr)
  return(corr)
}

Xmi <- function(dat, xvec, lags = 1){
  if (lags > 1){
    nb <- poly2nb(dat) %>% nb2listw(., zero.policy = T)
  } else {
    nb <- dat %>% poly2nb(.) %>% nblag(2) %>% nblag_cumul() %>% nb2listw(., zero.policy = T)
  }

  xvec %>% split(.,.) %>% lapply(function(x){
  xmi <- moran.test(dat %>% st_set_geometry(NULL) %>% .[,x] %>% unlist,
                    nb, zero.policy = TRUE, randomisation = F)
  return(xmi$estimate[1])
  }) %>% do.call(c,.) %>% mean -> xmi
  return(xmi)
}

Xmi_p <- function(dat, xvec, p = .1){
  nb <- knearneigh(dat, floor(nrow(dat)*p)) %>% knn2nb %>% nb2listw
  xvec %>% split(.,.) %>% lapply(function(x){
  xmi <- moran.test(dat %>% st_set_geometry(NULL) %>% .[,x] %>% unlist,
                    nb, zero.policy = TRUE, randomisation = F)
  return(xmi$estimate[1])
  }) %>% do.call(c,.) %>% mean -> xmi
  return(xmi)
}


Xmcmi <- function(dat, xvec, lags = 1, mode = 'a'){
  if (mode == 'a'){
  xmcmi <-
    data.frame(Mean_corr = Xmc(dat, xvec),
               Mean_MI = Xmi(dat, xvec, lags))
  } else {
  xmcmi <-
    data.frame(Mean_corr = Xmc(dat, xvec),
               Mean_MI = Xmi_p(dat, xvec))

  }
  return(xmcmi)
}





# Shortened version of RFSp (Hengl et al. 2018)
rfsp <- function(sf, tolerance = 0.9, covar.names, dep.name, standardize = TRUE){
  # Step 1: make a SpatialPixelsDataFrame
  if (standardize){
    sf <- sf %>% mutate_if(is.numeric, list(~scale(.) %>% as.vector))
  }
  
  sf.pt <- sf %>% as('Spatial') 
  sf.px <- sf %>% as('Spatial') %>% 
    SpatialPixelsDataFrame(.@points, data = .@data, tolerance = tolerance)
  # Step 2: Extract spatialbuffer distance
  grid.dist0 <- buffer.dist(sf.pt[dep.name], sf.px, as.factor(1:nrow(sf)))
  sf.grid <- over(sf.pt[dep.name], grid.dist0)
  dn0 <- paste(names(grid.dist0), collapse="+")
  print(dn0)
  # Step 3: spatial predictive components
  sf.covars <- paste('~', paste(covar.names, collapse = '+'))
  sf.spc <- spc(sf.px, as.formula(sf.covars))
  sf.form <- as.formula(paste(dep.name, " ~ ", dn0," + ", paste(names(sf.spc@predicted), collapse = "+")))
  sf.ov <- over(sf.pt[dep.name], sf.spc@predicted)
  sf.rm <- do.call(cbind, list(sf.pt@data[dep.name], sf.grid, sf.ov))
  nvars = length(covar.names) + length(sf.spc@predicted)
  nvars_threshold = floor(nvars / 5)
  sf.rfsp <- ranger(sf.form, sf.rm, 
                    importance="impurity", 
                    quantreg = FALSE,
                    mtry = min(30, nvars_threshold), 
                    num.trees = 500,
                    oob.error = TRUE)
  return(sf.rfsp)
}





## Post-hoc analysis of rfsp
posthoc_rfsp <- function(rfsp, sf, dep.name, mode = 'point', 
                         maxval = 30, return.residuals = FALSE,
                         network_wm = NULL, threshold = 5e3){
  if (mode != 'network'){
    sf.nei <- sf %>% knearneigh(k = floor(0.1 * nrow(.))) %>% knn2nb
  } else {
    if (is.null(network_wm)) stop('Error: no network distance matrix was specified.')
    if (any(class(network_wm) == 'matrix'))
    {
      sf.nei <- network_wm %>% 
        st_distmat_to_swm(., threshold) %>% 
        spdep::mat2listw(.) %>% 
        .$neighbours
    } else {
      sf.nei <- network_wm %>% 
        spdep::mat2listw(.) %>% 
        .$neighbours
    }
    #sf.nei <- network_wm %>% 
    #  st_distmat_to_swm(., threshold) %>% 
    #  spdep::mat2listw(.) %>% 
    #  .$neighbours
  }
  sf.neiw <- nb2listw(sf.nei, zero.policy = TRUE)
  ys <- scale(unlist(st_set_geometry(sf, NULL)[,dep.name]))
  resids <- rfsp$predictions - ys
  sf.ymt <- moran.test(ys, listw = sf.neiw, zero.policy = TRUE)$estimate[1]
  sf.mt <- moran.test(resids, listw = sf.neiw, zero.policy = TRUE)$estimate[1]
  sf.rmse <- sqrt(mean(resids^2))
  sf.adjr2 <- 1 - ((1 - (sum(resids^2)/sum(ys^2)))*((nrow(sf)-1)/(nrow(sf)-ncol(sf))))
  
  vimp <- rfsp$variable.importance #varimp(rfsp)
  posthoc <- data.frame(Rsquared_adj=sf.adjr2, MoranI_Y = sf.ymt, MoranI_resid = sf.mt, RMSE = sf.rmse)
  if (return.residuals){
    return(list(posthoc, resids = resids, varimp = vimp))
  } else {
    return(posthoc)
  }
}
##
corr_mat = function(data, yvec, xvec, title) {
  xvec_sort = sort(xvec)
  nmat = rbind(c(xvec_sort[1], yvec[1], xvec_sort[length(xvec)], yvec[length(yvec)])
               )

  corrplot::corrplot(round(cor(data[,c(yvec, xvec_sort)]), 5),
                    diag = FALSE, title = title,
                    mar = c(0, 0, 3.6, 0), number.cex = 0.64, type = "upper",
                    method = "number", tl.cex = 0.66) %>%
      corrRect(namesMat = nmat)
}


corr_mat_simp = function(data, xvec, title) {
  xvec_sort = sort(xvec)
  library(scales)

  colors <- scales::alpha(colorRampPalette(c("red", "white", "blue"))(100), alpha = 0.8)
  corrplot::corrplot(round(cor(data[,c(xvec_sort)]), 5),
                    diag = FALSE, title = title,
                    col = "black",
                    tl.col = "black",
                    cl.pos = "n",
                    mar = c(0, 0, 3.6, 0), number.cex = 1.2, type = "upper",
                    method = "number", tl.cex = 1.2)
}


### lm with y string and x strings
lmstring = function(yvec, xvec, data) {
  form = as.formula(paste(yvec, paste(xvec, collapse = "+"), sep = "~"))
  lmform = lm(form, data = data)
  lmcoefs = coef(lmform)
  return(lmcoefs)
}
