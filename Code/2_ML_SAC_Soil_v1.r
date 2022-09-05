## Machine learning and spatial autocorrelation in soil dataset
## Code written by Insang Song (sigmafelix@hotmail.com)
## Last revised 09/05/2022

# Change the directory before running code blocks below
# The directory should have
## - Data (directory): including excel files
## - Code (directory): including codes

# Sourcing the functions
set.seed(202006, "L'Ecuyer")

# Part 1: Data preparation ####
# Data import
agri <- read_excel('./Data/Agricultural.xls')
dune <- read_excel('./Data/Coastal dune.xls')
for1 <- read_excel('./Data/Mixed forest 1.xls')
for2 <- read_excel('./Data/Mixed forest 2.xls')

# Multicore use setting
#registerDoParallel(cores = 8)

### sf transformation
agri <- agri %>% st_as_sf(coords = 1:2)
dune <- dune %>% st_as_sf(coords = 1:2)
for1 <- for1 %>% st_as_sf(coords = 1:2)
for2 <- for2 %>% st_as_sf(coords = c(2,1))

### argument setting
### Names of Y- and X-variables
agri.xvec <- colnames(agri)[1:6]
agri.yvec <- colnames(agri)[7:13]
dune.xvec <- colnames(dune)[1:6]
dune.yvec <- colnames(dune)[7:16]
colnames(for1)[grep('Elevation', colnames(for1))] = 'elev'
colnames(for2)[grep('Elevation', colnames(for2))] = 'elev'
for1.xvec <- colnames(for1)[1:6]
for1.yvec <- colnames(for1)[7:17]
for2.xvec <- colnames(for2)[1:6]
colnames(for2)[16:17] <- str_replace_all(colnames(for2)[16:17], '-', '_')
for2.yvec <- colnames(for2)[7:17]

# Exploratory: correlation matrix
(corr_mat(as.data.frame(agri), agri.yvec, agri.xvec, "Agricultural field"))
(corr_mat(as.data.frame(dune), dune.yvec, dune.xvec, "Coastal dune"))
(corr_mat(as.data.frame(for1), for1.yvec, for1.xvec, "Mixed forest I"))
(corr_mat(as.data.frame(for2), for2.yvec, for2.xvec, "Mixed forest II"))



# Part 2: Run models ####
### Define the names of machine learning algorithms (compatible with mlr package)
lalg <- c('regr.svm', 'regr.ranger', 'regr.nnet')
lalg2 <- c('regr.lm', 'regr.svm', 'regr.ranger', 'regr.nnet')

# CV
# 020922: change seq.int(1, length(yvec)) into yvec for sorting
agri.micv <- agri.yvec %>% 
  split(.,.) %>% 
  lapply(function(x) mlr_learn(algs = lalg, agri, x, agri.xvec, cvmethod = 'CV', seed.num = 202006, std = TRUE))
dune.micv <- dune.yvec %>% 
  split(.,.) %>% 
  lapply(function(x) mlr_learn(algs = lalg, dune, x, dune.xvec, cvmethod = 'CV', seed.num = 202006, std = TRUE))
for1.micv <- for1.yvec %>% 
  split(.,.) %>% 
  lapply(function(x) mlr_learn(algs = lalg, for1, x, for1.xvec, cvmethod = 'CV', seed.num = 202006, std = TRUE))
for2.micv <- for2.yvec %>% 
  split(.,.) %>% 
  lapply(function(x) mlr_learn(algs = lalg, for2, x, for2.xvec, cvmethod = 'CV', seed.num = 202006, std = TRUE))

for2.micv2 <- for2.yvec[10:11] %>% 
  split(.,.) %>% 
  lapply(function(x) mlr_learn(algs = lalg2, for2, x, for2.xvec, cvmethod = 'CV', seed.num = 202006, std = TRUE))


## Run RFSp
agri.rfsp.l <- split(agri.yvec, agri.yvec) %>% 
  lapply(function(x) rfsp(agri, covar.names = colnames(agri)[1:6], dep.name = x))
dune.rfsp.l <- split(dune.yvec, dune.yvec) %>% 
  lapply(function(x) rfsp(dune, covar.names = colnames(dune)[1:6], dep.name = x))
for1.rfsp.l <- split(for1.yvec, for1.yvec) %>% 
  lapply(function(x) rfsp(for1, covar.names = colnames(for1)[1:6], dep.name = x))
for2.rfsp.l <- split(for2.yvec, for2.yvec) %>% 
  lapply(function(x) rfsp(for2, covar.names = colnames(for2)[1:6], dep.name = x))


# Part 3: extract model results ####
# 3-1: Moran's I and variable importance 
load("./Results/MLSAC_Rerun_cv10_sort_020922.RData")

# CV (this part will be use for the drafting!!!) ####
agri.phcv <- posthoc_mi(agri.micv, agri, agri.yvec, agri.xvec)
dune.phcv <- posthoc_mi(dune.micv, dune, dune.yvec, dune.xvec)
for1.phcv <- posthoc_mi(for1.micv, for1, for1.yvec, for1.xvec)
for2.phcv <- posthoc_mi(for2.micv, for2, for2.yvec, for2.xvec)

# RFSp
agri.rfsp.ph <- agri.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, agri, y, return.residuals = TRUE)[[1]],
         ., split(agri.yvec, agri.yvec), SIMPLIFY = FALSE) %>% 
  do.call(rbind, .) %>% 
  mutate(Dataset = 'Agriculture',
         Model = 'RFSp',
         dependent = sort(agri.yvec))
dune.rfsp.ph <- dune.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, dune, y, return.residuals = TRUE)[[1]],
         ., split(dune.yvec, dune.yvec), SIMPLIFY = FALSE) %>% 
  do.call(rbind, .) %>% 
  mutate(Dataset = 'Dune',
         Model = 'RFSp',
         dependent = sort(dune.yvec))
for1.rfsp.ph <- for1.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, for1, y, return.residuals = TRUE)[[1]],
         ., split(for1.yvec, for1.yvec), SIMPLIFY = FALSE) %>% 
  do.call(rbind, .) %>% 
  mutate(Dataset = 'Forest 1',
         Model = 'RFSp',
         dependent = sort(for1.yvec))
for2.rfsp.ph <- for2.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, for2, y, return.residuals = TRUE)[[1]],
         ., split(for2.yvec, for2.yvec), SIMPLIFY = FALSE) %>% 
  do.call(rbind, .) %>% 
  mutate(Dataset = 'Forest 2',
         Model = 'RFSp',
         dependent = sort(for2.yvec))




# 3-2. Extract residuals
# LM
agri.reslm <- agri.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(agri, x, agri.xvec, mode = 'point', return.residuals = TRUE))
dune.reslm <- dune.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(dune, x, dune.xvec, mode = 'point', return.residuals = TRUE))
for1.reslm <- for1.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(for1, x, for1.xvec, mode = 'point', return.residuals = TRUE))
for2.reslm <- for2.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(for2, x, for2.xvec, mode = 'point', return.residuals = TRUE))

# SF
agri.reslmr <- agri.reslm %>% 
  lapply(function(x) cbind(unlist(x$lmres), unlist(x$sfres))) %>% 
  do.call(cbind, .)
dune.reslmr <- dune.reslm %>% 
  lapply(function(x) cbind(unlist(x$lmres), unlist(x$sfres))) %>% 
  do.call(cbind, .)
for1.reslmr <- for1.reslm %>% 
  lapply(function(x) cbind(unlist(x$lmres), unlist(x$sfres))) %>% 
  do.call(cbind, .)
for2.reslmr <- for2.reslm %>% 
  lapply(function(x) cbind(unlist(x$lmres), unlist(x$sfres))) %>% 
  do.call(cbind, .)
colnames(agri.reslmr) <- str_c(rep(names(agri.reslm), each = 2), rep(c('_OLS', '_SF'), length(names(agri.reslm))))
colnames(dune.reslmr) <- str_c(rep(names(dune.reslm), each = 2), rep(c('_OLS', '_SF'), length(names(dune.reslm))))
colnames(for1.reslmr) <- str_c(rep(names(for1.reslm), each = 2), rep(c('_OLS', '_SF'), length(names(for1.reslm))))
colnames(for2.reslmr) <- str_c(rep(names(for2.reslm), each = 2), rep(c('_OLS', '_SF'), length(names(for2.reslm))))

# ML (CV)
agri.micvr <- agri.micv %>% 
  lapply(function(x) lapply(x, function(y) unlist(y$residuals)) %>% do.call(cbind, .)) %>% 
  do.call(cbind, .)
dune.micvr <- dune.micv %>% 
  lapply(function(x) lapply(x, function(y) unlist(y$residuals)) %>% do.call(cbind, .)) %>% 
  do.call(cbind, .)
for1.micvr <- for1.micv %>% 
  lapply(function(x) lapply(x, function(y) unlist(y$residuals)) %>% do.call(cbind, .)) %>% 
  do.call(cbind, .)
for2.micvr <- for2.micv %>% 
  lapply(function(x) lapply(x, function(y) unlist(y$residuals)) %>% do.call(cbind, .)) %>% 
  do.call(cbind, .)

calg_n <- c('_ANN', '_RF', '_SVM')
colnames(agri.micvr) <- str_c(rep(sort(agri.yvec), each = 3), rep(calg_n, length(agri.yvec)))
colnames(dune.micvr) <- str_c(rep(sort(dune.yvec), each = 3), rep(calg_n, length(dune.yvec)))
colnames(for1.micvr) <- str_c(rep(sort(for1.yvec), each = 3), rep(calg_n, length(for1.yvec)))
colnames(for2.micvr) <- str_c(rep(sort(for2.yvec), each = 3), rep(calg_n, length(for2.yvec)))


# RFSp 
agri.rfsp.phr <- agri.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, agri, y, return.residuals = TRUE),
         ., split(agri.yvec, agri.yvec), SIMPLIFY = FALSE) %>% 
  lapply(function(x) x$resids) %>% 
  do.call(cbind, .)
dune.rfsp.phr <- dune.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, dune, y, return.residuals = TRUE),
         ., split(dune.yvec, dune.yvec), SIMPLIFY = FALSE) %>% 
  lapply(function(x) x$resids) %>%
  do.call(cbind, .)
for1.rfsp.phr <- for1.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, for1, y, return.residuals = TRUE),
         ., split(for1.yvec, for1.yvec), SIMPLIFY = FALSE) %>% 
  lapply(function(x) x$resids) %>% 
  do.call(cbind, .)
for2.rfsp.phr <- for2.rfsp.l %>% 
  mapply(function(x, y) posthoc_rfsp(x, for2, y, return.residuals = TRUE),
         ., split(for2.yvec, for2.yvec), SIMPLIFY = FALSE) %>% 
  lapply(function(x) x$resids) %>% 
  do.call(cbind, .)

colnames(agri.rfsp.phr) <- str_c(sort(agri.yvec), '_rfsp')
colnames(dune.rfsp.phr) <- str_c(sort(dune.yvec), '_rfsp')
colnames(for1.rfsp.phr) <- str_c(sort(for1.yvec), '_rfsp')
colnames(for2.rfsp.phr) <- str_c(sort(for2.yvec), '_rfsp')


## Varimp ####
# LM
agri.vilm <- agri.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(agri, x, agri.xvec, mode = 'point', extract.coef = TRUE))
dune.vilm <- dune.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(dune, x, dune.xvec, mode = 'point', extract.coef = TRUE))
for1.vilm <- for1.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(for1, x, for1.xvec, mode = 'point', extract.coef = TRUE))
for2.vilm <- for2.yvec %>% split(.,.) %>% 
  lapply(function(x) sf_analysis(for2, x, for2.xvec, mode = 'point', extract.coef = TRUE))
watr.vilm <- watr.yvecs %>% 
           mapply(function(dat, yvec, xvec, swm) split(yvec, yvec) %>% 
                    lapply(function(yv) sf_analysis(dat, yv, xvec, 
                                                   mode = 'network', 
                                                   gdbpath = swm, 
                                                   threshold = thres,
                                                   extract.coef = TRUE)
                           ), watr, ., watr.xvecs, dms, SIMPLIFY = FALSE)

# SF
agri.resvi <- agri.vilm %>% 
  do.call(rbind, .) %>% 
  dplyr::select(-coef) %>% 
  pivot_wider(names_from = Model, values_from = Rank) %>% 
  pivot_wider(values_from = 3:4, names_from = dependent)
dune.resvi <- dune.vilm %>% 
  do.call(rbind, .) %>% 
  dplyr::select(-coef) %>% 
  pivot_wider(names_from = Model, values_from = Rank) %>% 
  pivot_wider(values_from = 3:4, names_from = dependent)
for1.resvi <- for1.vilm %>% 
  do.call(rbind, .) %>% 
  dplyr::select(-coef) %>% 
  pivot_wider(names_from = Model, values_from = Rank) %>% 
  pivot_wider(values_from = 3:4, names_from = dependent)
for2.resvi <- for2.vilm %>% 
  do.call(rbind, .) %>% 
  dplyr::select(-coef) %>% 
  pivot_wider(names_from = Model, values_from = Rank) %>% 
  pivot_wider(values_from = 3:4, names_from = dependent)

colnames(agri.resvi) <- colnames(agri.resvi) %>% str_split('_') %>% lapply(rev) %>% 
  lapply(function(x) str_c(x, collapse='_')) %>% do.call(c,.)
colnames(dune.resvi) <- colnames(dune.resvi) %>% str_split('_') %>% lapply(rev) %>% 
  lapply(function(x) str_c(x, collapse='_')) %>% do.call(c,.)
colnames(for1.resvi) <- colnames(for1.resvi) %>% str_split('_') %>% lapply(rev) %>% 
  lapply(function(x) str_c(x, collapse='_')) %>% do.call(c,.)
colnames(for2.resvi) <- colnames(for2.resvi) %>% str_split('_') %>% lapply(rev) %>% 
  lapply(function(x) str_c(x, collapse='_')) %>% do.call(c,.)


# Part 4: Combine results ####
soil_ranks <- bind_rows(
  agri.phcv[[2]] %>% mutate(Dataset = 'Agriculture'),
  dune.phcv[[2]] %>% mutate(Dataset = 'Dune'),
  for1.phcv[[2]] %>% mutate(Dataset = 'Forest 1'),
  for2.phcv[[2]] %>% mutate(Dataset = 'Forest 2')
  ) %>%
  ungroup() %>%
  mutate(independent = if_else(independent == 'Elevation', 'elev', independent),
        Model = plyr::mapvalues(Model, unique(Model), c('OLS', 'SF', 'ANN', 'RF', 'SVM')),
        Model = factor(Model, levels = c('OLS', 'SF', 'SVM', 'ANN', 'RF')))


## end of file ##
