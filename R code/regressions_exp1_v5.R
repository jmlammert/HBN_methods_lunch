##Version 5: multiclass classification

# INSTALL/CALL PACKAGES -----------------------------------------------------------------------

.libPaths(new="./other/r_libs")

library(easypackages)
libraries("tidyverse", "janitor", "naniar", "cluster", "Rtsne", "dtw", "dendextend",
          "KernSmooth", "dtwclust", "glmnet", "psych", "caret", "splitTools")

set.seed(1003)
rand_seed = 1003


# IMPORT DATA -------------------------------------------------------------

dat <- clean_names(read_csv("./data/tc_behav_clusters_k4.csv")) ### CHANGE CLUSTER DATASET HERE

# ITERABLES ---------------------------------------------------------------

network_clusters <- c("net_1_clus", "net_2_clus", "net_3_clus", "net_4_clus", "net_5_clus", "net_6_clus", "net_7_clus", "net_8_clus")
networks <- c("net_1", "net_2", "net_3", "net_4", "net_5", "net_6", "net_7", "net_8")

# REGRESSION -- ALL VARIABLES ----------------------------------------------

#standardize data
dat[4:30] <- scale(dat[4:30])
dat[4:30] <- rescale(dat[4:30], mean=100, df=TRUE)

#set minority/majority data for classification--minority = 1
##for (n in network_clusters){
  ###dat[[n]] <- factor(dat[[n]], levels = c(1,2), labels=c(1,2))
  ###print(table(dat[[n]]))
  ##if (sum(dat[[n]]==1) > sum(dat[[n]]==2)){
    ##print("Cluster 1 larger:")
    ##print(table(dat[[n]]))
    ##dat[[n]] <- factor(dat[[n]], levels = c(1,2), labels=c(2,1))
    ##print("After:")
    ##print(table(dat[[n]]))
  ##} 
  ##else {
  ##  print("ok")
  ##}
  ##dat[[n]] <- as.numeric(as.character(dat[[n]]))
  ##}

# train/test split
train_dfs <- c()
test_dfs <- c()
for (n in network_clusters){
  inds <- partition(dat[[n]], p = c(train = 0.80, test = 0.2))
  train <- dat[inds$train, ]
  test <- dat[inds$test, ]
  train_dfs <- append(train_dfs, list(train))
  test_dfs <- append(test_dfs, list(test))
}

# LASSO -- ALL VARIABLES -------------------------------------------------------------------

#create lists to save outputs from each network for later access
las_lasso_mods <- c()
las_predplots <- c()
las_cv_mods <- c()
las_lamplots <- c()
las_best_lams <- c()
las_preds <- c()
las_accs <- c()
las_conmats <- c()
las_perf_meas <- c()
las_final_mod <- c()
las_final_coefs <- c()

for (n in 1:8){ #repeat for each network clustering (y)
  net <- network_clusters[[n]]
  y <- dat[[net]]
  x <- dat[3:25]
  xtrain <- train_dfs[[n]]
  xtest <- test_dfs[[n]]
  ytrain <- as.matrix(xtrain[[net]])
  ytest <- as.numeric(xtest[[net]])
  folds <- create_folds(y=xtrain[[net]], k=10, seed=rand_seed, invert=TRUE) # create validation folds
  train_dfs[[n]]$foldids <- rep(NA, nrow(xtrain))
  for (f in 1:length(folds)){ 
    inds <- folds[[f]]
    train_dfs[[n]]$foldids[inds] <- f
    }
  xtrain <- xtrain[3:25] #leave out composites
  xtest <- xtest[3:25]
##
  fraction <- table(ytrain)/length(ytrain)
  weights <- data.frame(1 - fraction[ytrain])
  weights <- weights$Freq
##  
  all_fraction <- table(y)/length(y)
  all_weights <- data.frame(1-all_fraction[y])
  all_weights <- all_weights$Freq

  xvalid <- data.frame(train_dfs[[n]])
  xvalid$foldids <- train_dfs[[n]]$foldids
  xvalid$weights <- weights
  
  lasso_mod = glmnet(as.matrix(xtrain), ytrain, alpha = 1, family="multinomial", weights = weights)
  plot(lasso_mod, label=TRUE) #most important predictors enter earlier
  predplot <- recordPlot()
  
  cvlasso <- cv.glmnet(as.matrix(xvalid[3:25]), as.matrix(xvalid[[net]]), family="multinomial", alpha=1, type.measure = "class", nfolds=10, foldid = xvalid$foldids, weights = xvalid$weights)
  plot(cvlasso)
  lamplot <- recordPlot()
  lamlasso <-  cvlasso$lambda.min
  
  #evaluate prediction accuracy with test data
  ##predlasso <- predict(lasso_mod, newx = as.matrix(xtest), family="multinomial",type="class", s=lamlasso) #class predictions
  predlasso <- predict(cvlasso, newx = as.matrix(xtest), family="multinomial",type="class", s=lamlasso) #class predictions
  acclasso <- length(which(sapply(mapply(`%in%`, ytest, predlasso), isTRUE))) / length(ytest) #output accuracy 
 
  #conmat <- confusion.glmnet(lasso_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lamlasso)
  conmat <- confusionMatrix(factor(as.vector(predlasso)), factor(ytest), positive="1", mode="everything")

  perf <- assess.glmnet(lasso_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lamlasso)
  #roc <- roc.glmnet(lasso_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lamlasso)
  
  #fit cross-validated model on full dataset
  final_mod <- glmnet(x, y, alpha=1, weights = all_weights)
  final_coefs <- predict(final_mod, type="coefficients", s = lamlasso)[1:24,]
  
  #save outputs
  las_lasso_mods <- append(las_lasso_mods, list(lasso_mod))
  las_predplots <- append(las_predplots, list(predplot))
  las_cv_mods <- append(las_cv_mods, list(cvlasso))
  las_lamplots <- append(las_lamplots, list(lamplot))
  las_best_lams <- append(las_best_lams, list(lamlasso))
  las_preds <- append(las_preds, list(predlasso)) 
  las_accs <- append(las_accs, list(acclasso))
  las_conmats <- append(las_conmats, list(conmat))
  las_perf_meas <- append(las_perf_meas, list(perf))
  las_final_mod <- append(las_final_mod, list(final_mod))
  las_final_coefs <- append(las_final_coefs, list(final_coefs))
}

#get top coefficients for lam with top variables in each network
top_coefs <- c()
for (n in 1:8){ 
  net <- network_clusters[n]
  y <- dat[[net]]
  top_mod <- glmnet(x, y, alpha=1)
  cfs <- predict(top_mod, type="coefficients", s = 0.007)[1:24,]
  top_coefs <- append(top_coefs, list(cfs))
  print(n)
}

#extract variables from models
top_names <- c()
for (n in 1:8){
  print(paste("Network", n))
  top <- (sort(abs(top_coefs[[n]]), decreasing=TRUE))
  names <- names(top[2:6])
  print(paste(names))
  top_names <- append(top_names, list(names))
}

top_df <- data.frame(top_names[1:8])
colnames(top_df) <- 1:8

# RIDGE -- ALL VARIABLES -------------------------------------------------------

#create lists to save outputs from each network for later access
rid_ridge_mods <- c()
rid_predplots <- c()
rid_cv_mods <- c()
rid_lamplots <- c()
rid_best_lams <- c()
rid_preds <- c()
rid_accs <- c()
rid_conmats <- c()
rid_perf_meas <- c()
rid_final_mod <- c()
rid_final_coefs <- c()

for (n in 1:8){ #repeat for each network clustering (y)
  net <- network_clusters[[n]]
  y <- dat[[net]]
  x <- dat[3:25]
  xtrain <- train_dfs[[n]]
  xtest <- test_dfs[[n]]
  ytrain <- as.matrix(xtrain[[net]])
  ytest <- as.numeric(xtest[[net]])
  folds <- create_folds(y=xtrain[[net]], k=10, seed=rand_seed, invert=TRUE) # create validation folds
  train_dfs[[n]]$foldids <- rep(NA, nrow(xtrain))
  for (f in 1:length(folds)){ 
    inds <- folds[[f]]
    train_dfs[[n]]$foldids[inds] <- f
  }
  xtrain <- xtrain[3:25] #leave out composites
  xtest <- xtest[3:25]
  
  fraction <- table(ytrain)/length(ytrain)
  weights <- data.frame(1 - fraction[ytrain])
  weights <- weights$Freq
  
  test_fraction <- table(ytest)/length(ytest)
  test_weights <- data.frame(1-test_fraction[ytest])
  test_weights <- test_weights$Freq
  
  all_fraction <- table(y)/length(y)
  all_weights <- data.frame(1-all_fraction[y])
  all_weights <- all_weights$Freq
  
  xvalid <- data.frame(train_dfs[[n]])
  xvalid$weights <- weights
  
  ridge_mod <- glmnet(xtrain, ytrain, alpha=0, family = "multinomial", weights = weights)
  #coef(ridge_mod) #can get coefficients at specific lam
  plot(ridge_mod, xvar = "lambda", label = TRUE) #most important predictors shrink slower 
  predplot <- recordPlot()
  
  cv_mod <- cv.glmnet(as.matrix(xvalid[3:25]), as.matrix(xvalid[[net]]), familiy="multinomial", alpha=0, type.measure = "deviance", nfolds=10, foldid = xvalid$foldids, weights=xvalid$weights)
  plot(cv_mod)
  lamplot <- recordPlot()
  lam <- cv_mod$lambda.min #value of lambda that gives minimum mean cross-validated error
  
  #evaluate prediction accuracy with test data and best lambda
  pred <- predict(ridge_mod, newx = as.matrix(xtest), type="class", s=lam)
  acc <- length(which(sapply(mapply(`%in%`, ytest, pred), isTRUE))) / length(ytest) #output accuracy
  
  #conmat <- confusion.glmnet(ridge_mod, newx = as.matrix(xtest), newy=as.factor(ytest), family="multinomial", s=lam)
  conmat <- confusionMatrix(factor(as.vector(pred)), factor(ytest), positive="1", mode="everything")
  
  perf <- assess.glmnet(ridge_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lam)
  
  #fit cross-validated model on full dataset
  final_mod <- glmnet(x, y, alpha=0, weights=all_weights)
  final_coefs <- predict(final_mod, type="coefficients", s = lam)[1:24,]
  
  #save outputs
  rid_ridge_mods <- append(rid_ridge_mods, list(ridge_mod))
  rid_predplots <- append(rid_predplots, list(predplot))
  rid_cv_mods <- append(rid_cv_mods, list(cv_mod))
  rid_lamplots <- append(rid_lamplots, list(lamplot))
  rid_best_lams <- append(rid_best_lams, list(lam))
  rid_preds <- append(rid_preds, list(pred)) 
  rid_accs <- append(rid_accs, list(acc))
  rid_conmats <- append(rid_conmats, list(conmat))
  rid_perf_meas <- append(rid_perf_meas, list(perf))
  rid_final_mod <- append(rid_final_mod, list(final_mod))
  rid_final_coefs <- append(rid_final_coefs, list(final_coefs))
}

# REGRESSION -- COMPOSITE VARIABLES ----------------------------------------------

# LASSO -- COMPOSITE VARIABLES ----------------------------------------------

#create lists to save outputs from each network for later access
las2_lasso_mods <- c()
las2_predplots <- c()
las2_cv_mods <- c()
las2_lamplots <- c()
las2_best_lams <- c()
las2_preds <- c()
las2_accs <- c()
las2_conmats <- c()
las2_perf_meas <- c()
las2_final_mod <- c()
las2_final_coefs <- c()

for (n in 1:8){ #repeat for each network clustering (y)
  net <- network_clusters[[n]]
  y <- dat[[net]]
  x <- dat[c(13, 17, 26:30)]
  xtrain <- train_dfs[[n]]
  xtest <- test_dfs[[n]]
  ytrain <- as.matrix(xtrain[[net]])
  ytest <- as.numeric(xtest[[net]])
  folds <- create_folds(y=xtrain[[net]], k=10, seed=rand_seed, invert=TRUE) # create validation folds
  train_dfs[[n]]$foldids <- rep(NA, nrow(xtrain))
  for (f in 1:length(folds)){ 
    inds <- folds[[f]]
    train_dfs[[n]]$foldids[inds] <- f
  }
  xtrain <- xtrain[c(13, 17, 26:30)] #composites
  xtest <- xtest[c(13, 17, 26:30)]

  fraction <- table(ytrain)/length(ytrain)
  weights <- data.frame(1 - fraction[ytrain])
  weights <- weights$Freq
  
  all_fraction <- table(y)/length(y)
  all_weights <- data.frame(1-all_fraction[y])
  all_weights <- all_weights$Freq
  
  xvalid <- data.frame(train_dfs[[n]])
  xvalid$foldids <- train_dfs[[n]]$foldids
  xvalid$weights <- weights
  
  lasso_mod = glmnet(as.matrix(xtrain), as.matrix(ytrain), alpha = 1, family="multinomial", weights = weights, type.multinomial = "grouped")

  plot(lasso_mod, label=TRUE) #most important predictors enter earlier
  predplot <- recordPlot()

  cvlasso <- cv.glmnet(as.matrix(xvalid[c(13, 17, 26:30)]), as.matrix(xvalid[[net]]), family="multinomial", alpha=1, type.measure = "class", nfolds=10, foldid = xvalid$foldids, weights = xvalid$weights)
  plot(cvlasso)
  lamplot <- recordPlot()
  lamlasso <-  cvlasso$lambda.min
  
  #evaluate prediction accuracy with test data
  ##predlasso <- predict(lasso_mod, newx = as.matrix(xtest), family="multinomial",type="class", s=lamlasso) #class predictions
  predlasso <- predict(cvlasso, newx = as.matrix(xtest), family="multinomial",type="class", s=lamlasso) #class predictions
  acclasso <- length(which(sapply(mapply(`%in%`, ytest, predlasso), isTRUE))) / length(ytest) #output accuracy 
  
  #conmat <- confusion.glmnet(lasso_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lamlasso)
  conmat <- confusionMatrix(factor(as.vector(predlasso)), factor(ytest), positive="1", mode="everything")
  
  perf <- assess.glmnet(lasso_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lamlasso)
  #roc <- roc.glmnet(lasso_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lamlasso)
  
  #fit cross-validated model on full dataset
  final_mod <- glmnet(x, y, alpha=1, weights = all_weights)
  final_coefs <- predict(final_mod, type="coefficients", s = lamlasso)[1:8,]
  
  #save outputs
  las2_lasso_mods <- append(las2_lasso_mods, list(lasso_mod))
  las2_predplots <- append(las2_predplots, list(predplot))
  las2_cv_mods <- append(las2_cv_mods, list(cvlasso))
  las2_lamplots <- append(las2_lamplots, list(lamplot))
  las2_best_lams <- append(las2_best_lams, list(lamlasso))
  las2_preds <- append(las2_preds, list(predlasso)) 
  las2_accs <- append(las2_accs, list(acclasso))
  las2_conmats <- append(las2_conmats, list(conmat))
  las2_perf_meas <- append(las2_perf_meas, list(perf))
  las2_final_mod <- append(las2_final_mod, list(final_mod))
  las2_final_coefs <- append(las2_final_coefs, list(final_coefs))
}

# RIDGE -- COMPOSITE VARIABLES ----------------------------------------------

#create lists to save outputs from each network for later access
rid2_ridge_mods <- c()
rid2_predplots <- c()
rid2_cv_mods <- c()
rid2_lamplots <- c()
rid2_best_lams <- c()
rid2_preds <- c()
rid2_accs <- c()
rid2_conmats <- c()
rid2_perf_meas <- c()
rid2_final_mod <- c()
rid2_final_coefs <- c()

for (n in 1:8){ #repeat for each network clustering (y)
  net <- network_clusters[[n]]
  y <- dat[[net]]
  x <- dat[c(13, 17, 26:30)]
  xtrain <- train_dfs[[n]]
  xtest <- test_dfs[[n]]
  ytrain <- as.matrix(xtrain[[net]])
  ytest <- as.numeric(xtest[[net]])
  folds <- create_folds(y=xtrain[[net]], k=10, seed=rand_seed, invert=TRUE) # create validation folds
  train_dfs[[n]]$foldids <- rep(NA, nrow(xtrain))
  for (f in 1:length(folds)){ 
    inds <- folds[[f]]
    train_dfs[[n]]$foldids[inds] <- f
  }
  xtrain <- xtrain[c(13, 17, 26:30)] #leave out composites
  xtest <- xtest[c(13, 17, 26:30)]
  
  fraction <- table(ytrain)/length(ytrain)
  weights <- data.frame(1 - fraction[ytrain])
  weights <- weights$Freq
  
  all_fraction <- table(y)/length(y)
  all_weights <- data.frame(1-all_fraction[y])
  all_weights <- all_weights$Freq
  
  xvalid <- data.frame(train_dfs[[n]])
  xvalid$weights <- weights
  
  ridge_mod <- glmnet(xtrain, ytrain, alpha=0, family = "multinomial", weights = weights, type.multinomial = "grouped")
  #coef(ridge_mod) #can get coefficients at specific lam
  plot(ridge_mod, xvar = "lambda", label = TRUE) #most important predictors shrink slower 
  predplot <- recordPlot()
  
  cv_mod <- cv.glmnet(as.matrix(xvalid[c(13, 17, 26:30)]), as.matrix(xvalid[[net]]), familiy="multinomial", alpha=0, type.measure = "deviance", nfolds=10, foldid = xvalid$foldids, weights=xvalid$weights)
  plot(cv_mod)
  lamplot <- recordPlot()
  lam <- cv_mod$lambda.min #value of lambda that gives minimum mean cross-validated error
  
  #evaluate prediction accuracy with test data and best lambda
  pred <- predict(ridge_mod, newx = as.matrix(xtest), type="class", s=lam)
  acc <- length(which(sapply(mapply(`%in%`, ytest, pred), isTRUE))) / length(ytest) #output accuracy
  
  #conmat <- confusion.glmnet(ridge_mod, newx = as.matrix(xtest), newy=as.factor(ytest), family="multinomial", s=lam)
  conmat <- confusionMatrix(factor(as.vector(pred)), factor(ytest), positive="1", mode="everything")
  
  perf <- assess.glmnet(ridge_mod, newx = as.matrix(xtest), newy=ytest, family="multinomial", s=lam)
  
  #fit cross-validated model on full dataset
  final_mod <- glmnet(x, y, alpha=0, weights=all_weights)
  final_coefs <- predict(final_mod, type="coefficients", s = lam)[1:8,]
  
  #save outputs
  rid2_ridge_mods <- append(rid2_ridge_mods, list(ridge_mod))
  rid2_predplots <- append(rid2_predplots, list(predplot))
  rid2_cv_mods <- append(rid2_cv_mods, list(cv_mod))
  rid2_lamplots <- append(rid2_lamplots, list(lamplot))
  rid2_best_lams <- append(rid2_best_lams, list(lam))
  rid2_preds <- append(rid2_preds, list(pred)) 
  rid2_accs <- append(rid2_accs, list(acc))
  rid2_conmats <- append(rid2_conmats, list(conmat))
  rid2_perf_meas <- append(rid2_perf_meas, list(perf))
  rid2_final_mod <- append(rid2_final_mod, list(final_mod))
  rid2_final_coefs <- append(rid2_final_coefs, list(final_coefs))
}

# REGRESSION -- IMPORTANT VARIABLES ----------------------------------------------

#subset datasets in each network to important variables as determined by LASSO
xtrains <- c()
xtests <- c()
xalls <- c()
for (n in 1:8){
  train <- train_dfs[[n]][top_names[[n]]]
  test <- test_dfs[[n]][top_names[[n]]]
  all <- dat[top_names[[n]]]
  xtrains <- append(xtrains, list(train))
  xtests <- append(xtests, list(test))
  xalls <- append(xalls, list(all))
}

# RIDGE -- IMPORTANT VARIABLES --------------------------------------------

#create lists to save outputs from each network for later access
rid3_ridge_mods <- c()
rid3_predplots <- c()
rid3_cv_mods <- c()
rid3_lamplots <- c()
rid3_best_lams <- c()
rid3_preds <- c()
rid3_accs <- c()
rid3_conmats <- c()
rid3_perf_meas <- c()
rid3_final_mod <- c()
rid3_final_coefs <- c()

for (n in 1:8){ #repeat for each network clustering (y)
  ##net <- network_clusters[[n]]
  ##y = dat[[net]]
  ##ytrain <- train_dat[[net]] #var
  ##ytest <- test_dat[[net]] #var
  net <- network_clusters[[n]]
  y <- dat[[net]]
  x <- dat[3:25]
  xtrain <- xtrains[[n]]
  xtest <- xtests[[n]]
  ytrain <- as.matrix(train_dfs[[n]][[net]])
  ytest <- as.numeric(test_dfs[[n]][[net]])
  folds <- create_folds(y=train_dfs[[n]][[net]], k=10, seed=rand_seed, n_bins=2, invert=TRUE) # create validation folds
  train_dfs[[n]]$foldids <- rep(NA, nrow(xtrain))
  for (f in 1:length(folds)){ 
    inds <- folds[[f]]
    train_dfs[[n]]$foldids[inds] <- f
  }

  fraction <- table(ytrain)/length(ytrain)
  weights <- data.frame(1 - fraction[ytrain])
  weights <- weights$Freq
  
  test_fraction <- table(ytest)/length(ytest)
  test_weights <- data.frame(1-test_fraction[ytest])
  test_weights <- test_weights$Freq
  
  all_fraction <- table(y)/length(y)
  all_weights <- data.frame(1-all_fraction[y])
  all_weights <- all_weights$Freq
  
  xvalid <- data.frame(xtrain)
  xvalid$weights <- weights
  xvalid$net <- c(ytrain)
  xvalid$foldids <- train_dfs[[n]]$foldids

  ridge_mod <- glmnet(xtrain, ytrain, alpha=0, family = "multinomial", weight=weights)
  #coef(ridge_mod) #can get coefficients at specific lam
  plot(ridge_mod, xvar = "lambda", label = TRUE) #most important predictors shrink slower 
  predplot <- recordPlot()
  
  #cv_mod <- cv.glmnet(as.matrix(xvalid[[n]]), ytrain, familiy="multinomial", alpha=0, type.measure = "mse")
  cv_mod <- cv.glmnet(as.matrix(xvalid[1:5]), as.matrix(xvalid$net), familiy="multinomial", alpha=0, type.measure = "deviance", nfolds=10, foldid = xvalid$foldids, weights=xvalid$weights)
  plot(cv_mod)
  lamplot <- recordPlot()
  lam <- cv_mod$lambda.min #value of lambda that gives minimum mean cross-validated error
  
  #evaluate prediction accuracy with test data and best lambda
  pred <- predict(ridge_mod, newx = as.matrix(xtests[[n]]), type="class", s=lam)
  acc <- length(which(sapply(mapply(`%in%`, ytest, pred), isTRUE))) / length(ytest) #output accuracy
  
  conmat <- confusionMatrix(factor(as.vector(pred)), factor(ytest), positive="1", mode="everything")
  perf <- assess.glmnet(ridge_mod, newx = as.matrix(xtests[[n]]), newy=ytest, family="multinomial", s=lam)
  
  #fit cross-validated model on full dataset
  final_mod <- glmnet(xalls[[n]], y, alpha=0, weights=all_weights)
  final_coefs <- predict(final_mod, type="coefficients", s = lam)[1:5,]
  
  #save outputs
  rid3_ridge_mods <- append(rid3_ridge_mods, list(ridge_mod))
  rid3_predplots <- append(rid3_predplots, list(predplot))
  rid3_cv_mods <- append(rid3_cv_mods, list(cv_mod))
  rid3_lamplots <- append(rid3_lamplots, list(lamplot))
  rid3_best_lams <- append(rid3_best_lams, list(lam))
  rid3_preds <- append(rid3_preds, list(pred)) 
  rid3_accs <- append(rid3_accs, list(acc))
  rid3_conmats <- append(rid3_conmats, list(conmat))
  rid3_perf_meas <- append(rid3_perf_meas, list(perf))
  rid3_final_mod <- append(rid3_final_mod, list(final_mod))
  rid3_final_coefs <- append(rid3_final_coefs, list(final_coefs))
}

# COMPARE MODELS ----------------------------------------------------------

models <- list("las_all"=las_conmats, "rid_all"=rid_conmats, "las_comp"=las2_conmats, "rid_comp"=rid2_conmats, "rid_fs"=rid3_conmats)
info_dfs <- c()
  for (m in 1:length(models)){
    model <- models[[m]]
    modname <- names(models[m])
    for (n in 1:8){
      netmod <- data.frame(model[[n]]$byClass)
      netmod$model <- modname
      netmod$acc <- model[[n]]$overall["Accuracy"]
      netmod$p_acc <- model[[n]]$overall["AccuracyPValue"]
      netmod$net <- n
      netmod$cluster <- c(1,2,3,4) ##
      info_dfs <- append(info_dfs, list(netmod))
  }
}

all_info <- clean_names(Reduce(rbind, info_dfs))
all_infos <- split(all_info, all_info$net)

# WRITE DATA ------------------------------------------------------------------

##write.csv(all_info, "./data/binom_reg.csv", row.names=TRUE)
for (n in 1:8){
  fn <- paste0("4multi_reg_net_", n, ".csv")
  write.csv(all_infos[[n]], paste0("./data/LR_outputs/", fn), row.names=FALSE)
}
