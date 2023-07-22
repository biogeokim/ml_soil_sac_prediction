## Plots and exporting results

ggplot(data = soil_ranks,
         mapping = aes(x = Model, y = Rank, group = independent)) +
    stat_summary(geom = 'point', fun = mean, size = 4, pch = 19, lwd = 2, mapping = aes(color = independent)) +
    stat_summary(geom = 'point', fun = mean, size = 4, pch = 1, lwd = 2, mapping = aes(group = independent)) +
    stat_summary(geom = 'line', fun = mean, pch = 1, lwd = 1, mapping = aes(color = independent)) +
    scale_y_reverse() +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.title = element_blank(),
          axis.text.x = element_text(size = 14)) +
    guides(color = guide_legend(nrow = 1))

soil_ranks_export = soil_ranks %>%
  dplyr::select(Dataset, Model, dependent, independent, Rank) 
write_csv(soil_ranks_export, "./Results/Rank_Table.csv")

# Spearman rank corr
soil_ranks %>%
  dplyr::select(Dataset, dependent, independent, Model, Rank) %>%
  pivot_wider(names_from = Model, values_from = Rank) %>%
  group_by(Dataset, dependent) %>%
  nest() %>%
  mutate(rankcorr = map(data, ~cor(.x[,(ncol(.x)-4):ncol(.x)]))) %>%
  mutate(rankdata = map(rankcorr, ~data.frame(Model = rownames(.x)[2:5], spearman = .x[2:5,1]))) %>%
  dplyr::select(-data, -rankcorr) %>%
  unnest(cols = c(rankdata)) %>%
  ungroup %>%
  mutate(Model = factor(Model, levels = c('SF', 'SVM', 'ANN', 'RF'))) -> soil_ranks_nested
soil_ranks_nested %>%
  ggplot(data = ., mapping = aes(x = spearman, fill = Model)) +
    facet_wrap(~Model) +
    geom_histogram(color = 'black', binwidth = 0.2, orientation = 0.2, position = position_dodge2(width = 0.5, preserve = 'single'),
                   boundary = -0.1) +
    scale_x_continuous(breaks = c(-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(-1.1, 1.1)) +
    labs(x = "Spearman's rank correlation coefficient",
         y = "Number of occurrences") +
    ggthemes::theme_base() +
    scale_fill_manual(values = c('red', 'pink', 'blue', 'green')) +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = 8)) -> gg_soil_ranks_corrs
gg_soil_ranks_corrs
ggsave(plot = gg_soil_ranks_corrs, filename = './Results/Correlation_Distribution.tiff',
           width = 16, height = 15, units = 'cm', dpi = 300, device = 'tiff')


## Flow: CV-selected optimal hyperparameters are used to fit the entire dataset;
## Other performance measures followed.
# Combine RFSp results
rfsp_results = 
  bind_rows(agri.rfsp.ph, dune.rfsp.ph, for1.rfsp.ph, for2.rfsp.ph)

# Ranking plots
bind_rows(
  agri.phcv[[2]] %>% mutate(Dataset = 'Agriculture'),
  dune.phcv[[2]] %>% mutate(Dataset = 'Dune'),
  for1.phcv[[2]] %>% mutate(Dataset = 'Forest 1'),
  for2.phcv[[2]] %>% mutate(Dataset = 'Forest 2')
) %>%
  mutate(Model = factor(Model, levels = c('OLS', 'SF', 'regr.svm', 'regr.nnet', 'regr.ranger'),
                      labels = c('OLS', 'SF', 'SVM', 'ANN', 'RF'))) %>%
  ggplot(data = .,
         mapping = aes(x = Model, y = Rank, group = independent)) +
    stat_summary(geom = 'line', fun = mean, pch = 1, lwd = 1, mapping = aes(color = independent)) +
    stat_summary(geom = 'point', fun = mean, size = 4, pch = 19, lwd = 2, mapping = aes(color = independent)) +
    stat_summary(geom = 'point', fun = mean, size = 4, pch = 1, lwd = 2, mapping = aes(group = independent)) +
    ylab("Average ranking of\npredictor variables") +
    scale_y_reverse(breaks = c(6,5,4,3,2,1), limits = c(6.1, 0.9)) +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 14))+
    guides(color = guide_legend(nrow = 1)) -> gg_soil_rank
ggsave(plot = gg_soil_rank,
       filename = "./Results/Rank.tiff", device = 'tiff',
       height = 15, width = 16, units = 'cm', dpi = 300)

# RMI~YMI
soil_results = bind_rows(
  agri.phcv[[1]] %>% mutate(Dataset = 'Agriculture'),
  dune.phcv[[1]] %>% mutate(Dataset = 'Dune'),
  for1.phcv[[1]] %>% mutate(Dataset = 'Forest 1'),
  for2.phcv[[1]] %>% mutate(Dataset = 'Forest 2'),
  rfsp_results
) 


soil_results_export = soil_results %>%
  mutate(Model = factor(Model, levels = c('OLS', 'SF', 'svm', 'nnet', 'ranger', 'RFSp'),
                      labels = c('OLS', 'SF', 'SVM', 'ANN', 'RF', 'RFSp'))) %>%
  dplyr::select(Dataset, Model, dependent, MoranI_Y, MoranI_resid, RMSE) %>%
  unique
write_csv(soil_results_export, './Results/MoranI_RMSE_Table.csv')

soil_results_export %>%
  ggplot(data = .,
         mapping = aes(x = MoranI_Y, y = MoranI_resid, color = Model)) +
    geom_point(aes(shape = Dataset), color = 'black') +
    geom_smooth(aes(color = Model), method = 'lm') +
    theme(legend.position = 'bottom')


# boxplots + violin for RMSE
soil_results %>%
  mutate(Model = factor(Model, levels = c('OLS', 'SF', 'svm', 'nnet', 'ranger', 'RFSp'),
                      labels = c('OLS', 'SF', 'SVM', 'ANN', 'RF', 'RFSp'))) %>%
  group_by(Model) %>%
  summarize(RMSE = mean(RMSE)) %>%
  ungroup
soil_results %>%
  mutate(Model = factor(Model, levels = c('OLS', 'SF', 'svm', 'nnet', 'ranger', 'RFSp'),
                      labels = c('OLS\n(0.91)', 'SF\n(0.75)', 'SVM\n(0.90)', 'ANN\n(0.77)', 'RF\n(0.44)', 'RFSp\n(0.82)'))) %>%
  ggplot(data = .,
         mapping = aes(x = Model, y = RMSE, fill = Model)) +
    geom_boxplot() +
    scale_fill_manual(values = c('orange', 'orange', 'skyblue', 'skyblue', 'skyblue', 'skyblue')) +
    ylab('Root mean squared errors') +
    theme_bw() +
    theme(legend.position = 'none',
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text = element_text(size = 14, color = 'black')) -> gg_soil_rmse

# boxplots + violin for MoranI_resid
soil_results %>%
  mutate(Model = factor(Model, levels = c('OLS', 'SF', 'svm', 'nnet', 'ranger', 'RFSp'),
                      labels = c('OLS', 'SF', 'SVM', 'ANN', 'RF', 'RFSp'))) %>%
  group_by(Model) %>%
  summarize(MoranI_resid = mean(MoranI_resid)) %>%
  ungroup
soil_results %>%
  mutate(Model = factor(Model, levels = c('OLS', 'SF', 'svm', 'nnet', 'ranger', 'RFSp'),
                      labels = c('OLS\n(0.15)', 'SF\n(-0.04)', 'SVM\n(0.17)', 'ANN\n(0.09)', 'RF\n(0.05)', 'RFSp\n(-0.05)'))) %>%
  ggplot(data = .,
         mapping = aes(x = Model, y = MoranI_resid, fill = Model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = 'dark grey', lwd = 1.5) +
    geom_boxplot() +
    scale_fill_manual(values = c('orange', 'orange', 'skyblue', 'skyblue', 'skyblue', 'skyblue')) +
    ylab(expression(paste("Moran\'s ", italic(I) ," of residuals", sep = ""))) +
    theme_bw() +
    theme(legend.position = 'none',
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text = element_text(size = 14, color = 'black')) -> gg_soil_residI


gg_soil_rmseresi = wrap_plots(gg_soil_rmse, gg_soil_residI, ncol = 2) + 
  plot_annotation(tag_levels = 'a',
                  theme = theme(plot.title = element_text(size = 20, family = 'bold')))
ggsave(plot = gg_soil_rmseresi,
       filename = './Results/RMSE_ResidI.tiff',
       width = 32, height = 15, units = 'cm', dpi = 300, device = 'tiff')

## Fig 2
soil_results_fig2 = soil_results_export %>%
  nest()
soil_results_export_miy = soil_results_export %>%
  as_tibble %>%
  dplyr::select(Dataset, dependent, MoranI_Y) %>%
  unique
extr_lmcoefs = function(lmres) {
  lmres = summary(lmres)
  rsq = lmres$adj.r.squared
  fsig = pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
  fsig = cut(fsig, c(0, 0.001, 0.01, 0.05, 0.1, 1),
              labels = c('***', '**', '*','.', ''))
  fsig = as.character(fsig)
  return(c(rsq, fsig))
}

soil_results_fig2b = soil_results_export %>%
  as_tibble %>%
  filter(!Model %in% c('SVM','ANN')) %>%
  dplyr::select(Dataset, dependent, Model, MoranI_resid) %>%
  pivot_wider(values_from = MoranI_resid, names_from = Model) %>%
  mutate(MoranI_delta_rf = RF - RFSp,
         MoranI_delta_lm = OLS - SF) %>%
  dplyr::select(Dataset, dependent, starts_with('MoranI')) %>%
  pivot_longer(cols = 3:4, names_to = 'flag', values_to = 'MoranI_delta') %>%
  left_join(soil_results_export_miy) %>%
  group_by(flag) %>%
  nest() %>%
  mutate(rsq = map_chr(.x = data, .f = ~extr_lmcoefs(lm(formula = MoranI_delta ~ MoranI_Y, data = .))[1] %>% as.numeric),
         pval = map_chr(.x = data, .f = ~extr_lmcoefs(lm(formula = MoranI_delta ~ MoranI_Y, data = .))[2])) %>%
  unnest(data) %>%
  mutate(grouplabel = paste(ifelse(grepl('_rf$', flag), 'RF-RFsp ', 'OLS-SF ('), expression(R^{2}), "=", round(as.numeric(rsq), 3), pval, ")", sep = ''))
fig2b = ggplot(data = soil_results_fig2b,
              mapping = aes(x = MoranI_Y, y = MoranI_delta, group = grouplabel, color = grouplabel)) +
        geom_point()




# Part 5: Concatenate residuals ####
## residual combination
# reslmr 
# micvr
# rfsp.phr
agri.resids <- cbind(agri.reslmr, agri.micvr, agri.rfsp.phr) %>% 
  as.data.frame %>% 
  mutate(Dataset = 'Agriculture') %>%
  dplyr::select(Dataset, sort(colnames(.))) %>% 
  bind_cols(st_coordinates(agri) %>% data.frame) %>%
  pivot_longer(cols = 2:(6 * length(agri.yvec) + 1)) %>%
  mutate(dependent = str_extract(name, str_c(agri.yvec, collapse = '|')),
         Model = str_extract(name, str_c(c('ANN', 'SVM', 'RF', 'OLS', 'SF', 'rfsp'), collapse = '|'))) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = Model, values_from = value)
dune.resids <- cbind(dune.reslmr, dune.micvr, dune.rfsp.phr) %>% 
  as.data.frame %>% 
  mutate(Dataset = 'Dune') %>%
  dplyr::select(Dataset, sort(colnames(.))) %>% 
  bind_cols(st_coordinates(dune) %>% data.frame) %>%
  pivot_longer(cols = 2:(6 * length(dune.yvec) + 1)) %>%
  mutate(dependent = str_extract(name, str_c(dune.yvec, collapse = '|')),
         Model = str_extract(name, str_c(c('ANN', 'SVM', 'RF', 'OLS', 'SF', 'rfsp'), collapse = '|'))) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = Model, values_from = value)
for1.resids <- cbind(for1.reslmr, for1.micvr, for1.rfsp.phr) %>% 
  as.data.frame %>% 
  mutate(Dataset = 'Forest 1') %>%
  dplyr::select(Dataset, sort(colnames(.))) %>% 
  bind_cols(st_coordinates(for1) %>% data.frame) %>%
  pivot_longer(cols = 2:(6 * length(for1.yvec) + 1)) %>%
  mutate(dependent = str_extract(name, str_c(for1.yvec, collapse = '|')),
         Model = str_extract(name, str_c(c('ANN', 'SVM', 'RF', 'OLS', 'SF', 'rfsp'), collapse = '|'))) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = Model, values_from = value)
for2.resids <- cbind(for2.reslmr, for2.micvr, for2.rfsp.phr) %>% 
  as.data.frame %>% 
  mutate(Dataset = 'Forest 2') %>%
  dplyr::select(Dataset, sort(colnames(.))) %>% 
  bind_cols(st_coordinates(for2) %>% data.frame) %>%
  pivot_longer(cols = 2:(6 * length(for2.yvec) + 1)) %>%
  mutate(dependent = str_extract(name, str_c(for2.yvec, collapse = '|')),
         Model = str_extract(name, str_c(c('ANN', 'SVM', 'RF', 'OLS', 'SF', 'rfsp'), collapse = '|'))) %>%
  dplyr::select(-name) %>%
  pivot_wider(names_from = Model, values_from = value)

soil_resids = 
  bind_rows(agri.resids,
            dune.resids,
            for1.resids,
            for2.resids)
write_csv(soil_resids, './Results/Residuals_Table.csv')

agri.resids %>%
  group_by(Dataset, dependent) %>%
  nest %>%
  mutate(CorrMat = map(data, ~cor(.x[,-1:-2]))) %>%
  .$CorrMat


# Concatenate the variable importance
agri.vi <- left_join(agri.resvi, agri.mitrvi) %>% 
  dplyr::select(1, sort(colnames(.)[-1]))
dune.vi <- left_join(dune.resvi, dune.mitrvi) %>% 
  dplyr::select(1, sort(colnames(.)[-1]))
for1.vi <- left_join(for1.resvi, for1.mitrvi) %>% 
  dplyr::select(1, sort(colnames(.)[-1]))
for2.vi <- left_join(for2.resvi, for2.mitrvi) %>% 
  dplyr::select(1, sort(colnames(.)[-1]))

# .mitrr: residuals

# Part VI: extract optimized params ####
# tip: the tune results are stored in every first list element
# tip: $x then as.data.frame
extr_optparams = function(lst, yvars) {
  depvars = sort(yvars)
  algs = c("ANN", "RF", 'SVM')

  tuneres = lapply(
    lst,
    function(l) { 
      
      lst_df = 
      do.call(bind_rows,      
        lapply(
          l,
          function(l_in) {
            as.data.frame(l_in[[1]]$x)
          }
        ))
      
      lst_df = mutate(lst_df, Model = algs)
    }
  )
  tuneres_df = do.call(bind_rows, tuneres) %>%
    mutate(dependent = rep(depvars, each = 3))# %>% 

  return(tuneres_df)


}

optparams_df = 
  bind_rows(
    extr_optparams(agri.micv, agri.yvec) %>%
      mutate(Dataset = "Agricultural field") %>%
      dplyr::select(Dataset, 1:(ncol(.)-1)),
    extr_optparams(dune.micv, dune.yvec) %>%
      mutate(Dataset = "Coastal dune") %>%
      dplyr::select(Dataset, 1:(ncol(.)-1)),
    extr_optparams(for1.micv, for1.yvec) %>%
      mutate(Dataset = "Mixed forest I") %>%
      dplyr::select(Dataset, 1:(ncol(.)-1)),
    extr_optparams(for2.micv, for2.yvec) %>%
      mutate(Dataset = "Mixed forest II") %>%
      dplyr::select(Dataset, 1:(ncol(.)-1))
  )

write_csv(optparams_df, "./Results/Optimal_parameters.csv")



## Corrplots
tiff("./Results/Corrplots_by_datasets.tiff", width = 8.8, height = 8, units = 'in',
     pointsize = 12, compression = 'lzw', bg = 'white', res = 300)
par(mfrow = c(2,2))
corr_mat_simp(as.data.frame(for1), for1.xvec, "Mixed forest I")
corr_mat_simp(as.data.frame(for2), for2.xvec, "Mixed forest II")
corr_mat_simp(as.data.frame(agri), agri.xvec, "Agricultural field")
corr_mat_simp(as.data.frame(dune), dune.xvec, "Coastal dune")
par(mfrow = c(1,1))
dev.off()


## Supplementary: all lm coefficients
agri_lmstring = 
    agri.yvec %>%
    split(.,.) %>%
    lapply(function(x) lmstring(x, sort(agri.xvec), data =agri)) %>%
    Reduce(rbind, .) %>%
    data.frame(data = "Agriculture", yvec = agri.yvec, .)
dune_lmstring = 
    dune.yvec %>%
    split(.,.) %>%
    lapply(function(x) lmstring(x, sort(dune.xvec), data =dune)) %>%
    Reduce(rbind, .) %>%
    data.frame(data = "Coastal dune", yvec = dune.yvec, .)
for1_lmstring = 
    for1.yvec %>%
    split(.,.) %>%
    lapply(function(x) lmstring(x, sort(for1.xvec), data =for1)) %>%
    Reduce(rbind, .) %>%
    data.frame(data = "Mixed forest I", yvec = for1.yvec, .)
for2_lmstring = 
    for2.yvec %>%
    split(.,.) %>%
    lapply(function(x) lmstring(x, sort(for2.xvec), data =for2)) %>%
    Reduce(rbind, .) %>%
    data.frame(data = "Mixed forest II", yvec = for2.yvec, .)
results_lmstring = 
    rbind(agri_lmstring, dune_lmstring,
          for1_lmstring, for2_lmstring)
write.csv(results_lmstring, "./Output/LM_Results.csv")

## END OF FILE ####