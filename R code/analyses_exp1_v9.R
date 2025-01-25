# Version 9: output 3&4 cluster data for multi-class regressions

# INSTALL/CALL PACKAGES -----------------------------------------------------------------------

.libPaths(new="./other/r_libs")

library(easypackages)
libraries("tidyverse", "janitor", "naniar", "cluster", "Rtsne", "dtw", "dendextend",
          "KernSmooth", "dtwclust", "glmnet", "psych", "ggrain")

set.seed(1003)
rand_seed = 1003

# IMPORT DATA ---------------------------------------------------------------------------------

#call preprocessed behavioural data, time course, and ISC
isc <- clean_names(read_csv("./data/0724_isc_pairs_nets1-8.csv"))
behav <- clean_names(read_csv("./data/0824_MRI-behav_comp-final.csv"))
tc <- clean_names(read_csv("./data/0724_timecourse_nets1-8.csv"))
ids <- clean_names(read_csv("./data/0724_ids_final.csv"))

clusd <- clean_names(read_csv("./data/tc_behav_clusters_sdtw_nosqr.csv"))
load("./data/k2_data_sdtw_nosqr.RData")

# TIDY DATA -----------------------------------------------------------------------------------

# ITERABLES -----------------------------------------------------------------------------------

#make iterables
networks <- c("net_1", "net_2", "net_3", "net_4", "net_5", "net_6", "net_7", "net_8")
net_names <- c("Network 1", "Network 2", "Network 3", "Network 4", "Network 5", "Network 6", "Network 7", "Network 8")
net_colours <- c("#781286", "#4682b4", "#00760e", "#f43afa", "#92b574", "#e69422", "#cd3e4f", "#0808c2")
anat_names <- c("Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "DMN", "ToM")
network_clusters <- c("net_1_clus", "net_2_clus", "net_3_clus", "net_4_clus", "net_5_clus", "net_6_clus", "net_7_clus", "net_8_clus")
net_df <- data.frame(networks, net_names, net_colours, anat_names)
cldists <- c("cldist1", "cldist2", "cldist3", "cldist4", "cldist5", "cldist6", "cldist7", "cldist8")

ivars <- c("celf_total", "ctopp_elision", "ctopp_blending", "ctopp_nwrep", "ctopp_symbol",
           "gfta_siw", "towre_swe", "towre_pde", "wisc_fsiq", "wisc_vsi", "wisc_vci", 
           "wisc_fri", "wisc_wmi", "wisc_psi", "srs_awr", "srs_cog", "srs_com", "srs_mot", "srs_rrb",
           "nih_card", "nih_flanker", "iv_lang", "iv_phon", "iv_soco", "iv_asd", "iv_attn", "age", "sex", "barratt")
scores <- ivars[1:21]

# DATA STRUCTURE ------------------------------------------------------------------------------

colnames(isc)[2:9] <- networks #set isc column names

#subset data to final subjects
ids <- subset(ids, p_id %in% behav$id)##fix for adding barratt
behav <- subset(behav, id %in% ids$p_id) 
isc <- subset(isc, sub_1 %in% ids$p_id)
isc <- subset(isc, sub_2 %in% ids$p_id)
nps <- nrow(ids) #get number of participants for loops
sub_ids <- c(ids$p_id)

#merge MRI with behavioural data
isc_behav <- merge(isc, behav, by.x=c("sub_1"), by.y=c("id"))
tc_behav <- merge(tc, behav, by.x=("p_id"), by.y=c("id"))

#fix timecourse data types
tc_behav <- tc_behav %>% #make single string back into list
  mutate_if(~any(str_detect(., fixed('|'))), ~str_split(., fixed('|')))
tc_behav$net_1[1][[1]][1] #access single item in list in df column

for (n in networks){ #make time course values numeric
  print(n)
  for (row in 1:nps){
    print(row)
    a <- as.numeric(unlist(tc_behav[n][[row,1]]))
    print(typeof(a))
    tc_behav[n][[row,1]] <- a
    #tc_behav[n][[row,1]] <- list(a)
  }
}
print(typeof(tc_behav$net_1))
print(typeof(tc_behav$net_1[[1]]))

#isolate time course data for each network

data_func <- function(net){
  d <- tc_behav %>% select(c("p_id", net))
  return(d)
}

dfs <- c()
for (n in networks){
  net_data <- data_func(n)
  dfs <- append(dfs, list(net_data))
}

#make data wide
dfs_wide <- c()
for (d in dfs){
  wide_data <- d %>% unnest_wider(2, names_sep = "_") %>% select(!c("p_id"))
  dfs_wide <- append(dfs_wide, list(wide_data))
}

#make data long
dfs_long <- c()
for (d in dfs_wide){
  long_data <- data.frame(t(d))
  dfs_long <- append(dfs_long, list(long_data))
}

# prepare correlation matrix data
cm <- data.frame(isc_behav)
miss_cors <- data.frame(isc_behav) #construct inverse correlations
miss_cors <- miss_cors %>% rename(sub_1=sub_2, sub_2=sub_1)
self_cors <- data.frame(sub_ids, sub_ids) #construct self correlations
for (n in net_df$networks){
  self_cors[n] <- 1
}
self_cors <- self_cors %>% rename(sub_1=sub_ids, sub_2=sub_ids.1)
cm <- data.frame(bind_rows(cm, self_cors, miss_cors)) #combine correlations = long form all networks

#make into wide correlation matrix for each network
cm_func <- function(net){
  c <- cm %>% select(sub_1, sub_2, net) %>% #select network
    pivot_wider(names_from = sub_1, values_from = net)
  c <- data.frame(c[order(c$sub_2),])
  rownames(c) <- c$sub_2
  c <- c[,2:(nps+1)] #remove subject column
  return(c)
}

cms = c() #store correlation matrices
for (n in networks){
  net_cm <- cm_func(n)
  cms <- append(cms, list(net_cm))
}

# DATA PREPROCESSING  ---------------------------------------------------------------------------

#smooth time series
# times series must be standardized before smoothing to allow for bandwidth
bw  <- 8 ## SET SMOOTHING KERNEL (BANDWIDTH)

# standardize
dfs_stand <- c() # store standardized data
for (d in dfs_wide){
  stand_df <- data.frame(t(d))
  stand_df <- as.data.frame(scale(stand_df, center=TRUE))
  stand_df$timepoint <- seq(1,245,1)
  dfs_stand <- append(dfs_stand, list(stand_df))
}

dfs_smooth <- c()
for (d in dfs_stand){
  smooth_df <- data.frame(d)
  for (p in 1:nps){ #smooth each participant column in each network df
    smoothed <- locpoly(d$timepoint, d[[p]], bandwidth = bw, gridsize = 245)
    smooth_df[[p]] <- smoothed$y
  }
  dfs_smooth <- append(dfs_smooth, list(smooth_df))
}

#trim time courses
start_time <- 10
end_time <- 235
for (n in 1:length(networks)){
  trim_df_stand <- data.frame(dfs_stand[[n]])
  trim_df_stand <- trim_df_stand[start_time:end_time,]
  dfs_stand[[n]] <- trim_df_stand
  trim_df_smooth <- data.frame(dfs_smooth[[n]])
  trim_df_smooth <- trim_df_smooth[start_time:end_time,]
  dfs_smooth[[n]] <- trim_df_smooth
  trim_df_long <- data.frame(dfs_long[[n]])
  trim_df_long <- trim_df_long[start_time:end_time,]
  dfs_long[[n]] <- trim_df_long
}

ttps <- nrow(dfs_smooth[[1]])
print(paste("Trimmed timepoints:", ttps)) #should be 226

#make data long for plotting (raw, standardized, smooth)
dfs_plot <- c()
for(n in 1:length(networks)){
  colnames(dfs_long[[n]]) <- tc_behav$p_id
  dfs_long[[n]]$timepoint <- dfs_stand[[n]]$timepoint
  long_data <- dfs_long[[n]] %>% pivot_longer(cols=-timepoint, names_to = "id", values_to = "raw")
  long_stand <- dfs_stand[[n]] %>% pivot_longer(cols=-timepoint, names_to = "id", values_to = "stand")
  long_smooth <- dfs_smooth[[n]] %>% pivot_longer(cols=-timepoint, names_to = "id", values_to = "smooth")
  long_data$stand <- long_stand$stand
  long_data$smooth <- long_smooth$smooth
  dfs_plot <- append(dfs_plot, list(long_data))
}

# ANALYSES ------------------------------------------------------------------------------------

# CLUSTERING ----------------------------------------------------------------------------------

# TIME SERIES CLUSTERING --------------------------------------------------

# prepare clustering data (time series only)
dfs_clus <- c()
for (d in dfs_smooth){
  clus_data <- as.data.frame(t(d[1:nps]))
  dfs_clus <- append(dfs_clus, list(clus_data))
}

# cluster -- TAKES SEVERAL MINUTES EACH; 5 ks (2-6) x 100 repetitions = 500 indices
c1 <- tsclust(dfs_clus[[1]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #Visual
c2 <- tsclust(dfs_clus[[2]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #Somatomotor
c3 <- tsclust(dfs_clus[[3]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #Ventral attention
c4 <- tsclust(dfs_clus[[4]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #Dorsal attention
c5 <- tsclust(dfs_clus[[5]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #Limbic
c6 <- tsclust(dfs_clus[[6]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #Frontoparietal
c7 <- tsclust(dfs_clus[[7]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #DMN
c8 <- tsclust(dfs_clus[[8]], type = "partitional", k=c(2,3,4), distance="sdtw", centroid = "pam", control = partitional_control(iter.max=20L, nrep=100), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(gamma=0.1, sqrt.dist=TRUE))) #ToM

#c5 <- tsclust(dfs_clus[[5]], type = "partitional", k=c(2,3,4), distance="dtw_basic", centroid = "pam", control = partitional_control(iter.max=20L, nrep=10), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(norm="L2"), cent=list(window.size=2))) #Limbic
#c6 <- tsclust(dfs_clus[[6]], type = "partitional", k=c(2,3,4), distance="dtw_basic", centroid = "pam", control = partitional_control(iter.max=20L, nrep=10), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(norm="L2"), cent=list(window.size=2))) #Frontoparietal
#c7 <- tsclust(dfs_clus[[7]], type = "partitional", k=c(2,3,4), distance="dtw_basic", centroid = "pam", control = partitional_control(iter.max=20L, nrep=10), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(norm="L2"), cent=list(window.size=2))) #DMN
#c8 <- tsclust(dfs_clus[[8]], type = "partitional", k=c(2,3,4), distance="dtw_basic", centroid = "pam", control = partitional_control(iter.max=20L, nrep=10), preproc=zscore, trace = TRUE, seed=rand_seed, args=tsclust_args(dist=list(norm="L2"), cent=list(window.size=2))) #ToM

cs <- list(c1, c2, c3, c4, c5, c6, c7, c8) #

# evaluate clusters -- TAKES SEVERAL MINUTES
c_evals <- c()
for (n in cs){
  c_list <- c()
  for (c in n){
    c_eval <- cvi(c)
    c_list <- append(c_list, list(c_eval))
  }
  c_evals <- append(c_evals, list(c_list))
}

# get silhouette scores per network/cluster number
n_count = 1
sil_list = list()
for (n in c_evals){
  k2 <- n[1:100] #separate each k
  k3 <- n[101:200]
  k4 <- n[201:300]
  ks <- list(k2, k3, k4)
  k_count <- 2
  for (k in ks){
    sil_scores <- c()
    for (r in seq(1,100,1)){
      sil_scores <- append(sil_scores, k[[r]][[1]])
    }
    mean_sil <- mean(sil_scores)
    max_sil <- max(sil_scores)
    max_ind <- which(sil_scores == max_sil)[1]
    print(max_ind)
    sil_data <- c(mean_sil, max_sil, k_count, n_count, max_ind)
    print(sil_data)
    sil_list <- append(sil_list, list(sil_data))
    k_count = k_count + 1
  }
  n_count = n_count + 1
}
sil_df <- data.frame(t(data.frame(sil_list)))
row.names(sil_df) <- NULL
colnames(sil_df) <- c("mean_sil", "max_sil", "k_count", "n_count", "max_ind")
#make wide for comparisons 
sil_wide <- data.frame(sil_df[2:4])
sil_wide <- pivot_wider(sil_wide, names_from = c(k_count), values_from = c(max_sil))

View(sil_wide)
prev_sil <- sil_df


# get clustering for k=2 index with max silhouette score
inds <- subset(sil_df, k_count==2)
count = 1
k2_data = c()
for (n in cs){
  k2 <- n[1:100]
  i <- inds["max_ind"][[1]][[count]]
  print(i)
  k2_best <- k2[i]
  print(k2_best)
  k2_data <- append(k2_data, list(k2_best))
  count=count+1
}

# get clustering for k=3 index with max silhouette score
inds <- subset(sil_df, k_count==3)
count = 1
k3_data = c()
for (n in cs){
  k3 <- n[101:200]
  i <- inds["max_ind"][[1]][[count]]
  print(i)
  k3_best <- k3[i]
  print(k3_best)
  k3_data <- append(k3_data, list(k3_best))
  count=count+1
}

# get clustering for k=4 index with max silhouette score
inds <- subset(sil_df, k_count==4)
count = 1
k4_data = c()
for (n in cs){
  k4 <- n[201:300]
  i <- inds["max_ind"][[1]][[count]]
  print(i)
  k4_best <- k4[i]
  print(k4_best)
  k4_data <- append(k4_data, list(k4_best))
  count=count+1
}

##** START HERE IF IMPORTING CLUSTER DATA **## -------
# add cluster labels to data
dat <- data.frame(tc_behav)
for(n in 1:length(networks)){
  netcluscol <- network_clusters[[n]]
  dat[[netcluscol]] <- k4_data[[n]][[1]]@cluster ### CHANGE HERE TO GET K DAT
}

corrdat <- data.frame(cm)
for (n in network_clusters){
  corrdat <- merge(corrdat, dat[c("p_id", n)], by.x=("sub_1"), by.y=c("p_id"))
}

#prepare distance matrices for plotting/saving

# extract distance matrices
dmats <- c()
for (n in k2_data){
  dmat <- n[[1]]@distmat
  dmats <- append(dmats, list(dmat))
}

#prepare distance matrix for plotting
dmats_plot <- c()
dmats_save <- c()
for (n in 1:8){
  d <- as.data.frame(dmats[[n]])
  rownames(d) <- rownames(cms[[n]])
  colnames(d) <- colnames(cms[[n]])
  dmats_save <- append(dmats_save, list(d))
  d <- d %>% rownames_to_column(var="row")
  d <- d %>% pivot_longer(cols=c(-"row"), names_to = "id", values_to = "dist")
  dmats_plot <- append(dmats_plot, d["dist"])
}
dmats_plot[c("sub_1", "sub_2")] <- d[c("row", "id")]
dmats_plot <- as.data.frame(dmats_plot)
colnames(dmats_plot) <- c(networks, "sub_1", "sub_2")

clus_dist <- data.frame(unlist(k2_data[[1]][[1]]@cldist), unlist(k2_data[[2]][[1]]@cldist), unlist(k2_data[[3]][[1]]@cldist), unlist(k2_data[[4]][[1]]@cldist)
                        , unlist(k2_data[[5]][[1]]@cldist), unlist(k2_data[[6]][[1]]@cldist), unlist(k2_data[[7]][[1]]@cldist), unlist(k2_data[[8]][[1]]@cldist))
names(clus_dist) <- cldists
dmats_plot <- cbind(dmats_plot, corrdat[c(network_clusters)], clus_dist)

# ISC CLUSTERING

#hierarchical clustering (for visualization purposes)
ddist <- diana(dmats_save[[1]], stand=TRUE, keep.diss=TRUE)
hc1 <- hclust(ddist$diss, method = "ward.D2")

plot(hc1)
summary(hc1)
hc1$order
dend <- as.dendrogram(hc1)
hord1 <- cutree(dend, h=5)#max(hc1$height)-1) #get cluster labels at highest merge ## SELECT CLUSTER LABEL LEVEL
#c <- cutree(dend, h=170) #get cluster labels at highest merge


# CLUSTER INFORMATION -------------------------------------------------------------------------

# get cluster 1 and 2 members for each network
clus_1_ids <- c()
clus_2_ids <- c()

for (n in network_clusters){
  clus1 <- dat[dat[[n]]==1,]
  clus2 <- dat[dat[[n]]==2,]  
  ids1 <- clus1$p_id
  ids2 <- clus2$p_id
  clus_1_ids <- append(clus_1_ids, list(ids1))
  clus_2_ids <- append(clus_2_ids, list(ids2))
}

clusdat <- data.frame(corrdat)

#within cluster ISC
isc_wthn_1 <- c()
isc_wthn_2 <- c()
for (n in 1:length(networks)){
  cd1 <- clusdat %>% filter(clusdat$sub_1 %in% clus_1_ids[[n]] & clusdat$sub_2 %in% clus_1_ids[[n]]) #both ids in cluster
  cd2 <- clusdat %>% filter(clusdat$sub_1 %in% clus_2_ids[[n]] & clusdat$sub_2 %in% clus_2_ids[[n]])
  net <- networks[[n]]
  iscs1 <- cd1[[net]]
  iscs1 <- iscs1[iscs1 != 1] #remove self-correlations
  isc_wthn_1 <- append(isc_wthn_1, list(iscs1))
  iscs2 <- cd2[[net]]
  iscs2 <- iscs2[iscs2 != 1]
  isc_wthn_2 <- append(isc_wthn_2, list(iscs2))
  print(paste(net, "mean within cluster ISC", "cluster 1:", mean(iscs1), "cluster 2:", mean(iscs2)))
}

#between cluster ISC
isc_btwn <- c()
for (n in 1:length(networks)){
  cd_btwn <- clusdat[-which(clusdat$sub_1 %in% clus_1_ids[[n]] & clusdat$sub_2 %in% clus_1_ids[[n]]),] #remove same-cluster pairs, sub1 and sub2 same cluster
  cd_btwn <- cd_btwn[-which(clusdat$sub_1 %in% clus_2_ids[[n]] & clusdat$sub_2 %in% clus_2_ids[[n]]),]
  net <- networks[[n]]
  net_isc_btwn <- cd_btwn[[net]]
  net_isc_btwn <- net_isc_btwn[net_isc_btwn != 1] #remove self-correlations
  isc_btwn <- append(isc_btwn, list(net_isc_btwn))
  print(paste(net, "mean between cluster ISC", mean(net_isc_btwn)))
}

#compare means
for (n in 1:length(networks)){
  print((paste(networks[[n]], "mean ISC -- between:", round(mean(isc_btwn[[n]]), 4), "within 1:", round(mean(isc_wthn_1[[n]]), 4), "within 2:", round(mean(isc_wthn_2[[n]]), 4))))
}

#statistically comapre ISC values between clusters, within/between clusters -- correct for multiple comparisons?
t1 <- t.test(cd1$net_8, cd2$net_8)
t2 <- t.test(cd1$net_8, cd_btwn$net_8)
t3 <- t.test(cd2$net_8, cd_btwn$net_8)

ggplot(data=cd1, aes(x=net_8)) +
  geom_histogram(aes(fill="Cluster 1"), alpha=0.6, bins=100) +
  geom_histogram(data=cd2, aes(x=net_8, fill="Cluster 2"), alpha =0.5, bins=100)+
  #geom_density(aes(fill="Cluster 1"), alpha=0.5) +
  #geom_density(data=cd2, aes(x=net_8, fill="Cluster 2"), alpha =0.5)+
  scale_fill_manual(values=c("#f7b148", "#8164cc"), labels = c(1,2), name="Cluster")+
  xlim(-.7,.7) +
  labs(y="Count", x="Pairwise ISC (r)") +
  theme_bw()
#geom_density()# +

# REGRESSION -----------------------------------------------------------

## SEE REGRESSIONS_EXP1 SCRIPT

# FIGURES -------------------------------------------------------------------------------------

##plotting functions
tick_func = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

# CLUSTER PLOTS -----------------------------------------------------------

# plot time series by cluster
plot(k2_data[[1]][[1]], type="series") #Visual #centroids, sc
plot(k2_data[[2]][[1]], type="series") #Somatomotor
plot(k2_data[[3]][[1]], type="series") #Dorsal Attention
plot(k2_data[[4]][[1]], type="series") #Ventral Attention
plot(k2_data[[5]][[1]], type="series") #Limbic
plot(k2_data[[6]][[1]], type="series") #Frontoparietal
plot(k2_data[[7]][[1]], type="series") #DMN
plot(k2_data[[8]][[1]], type="series") #ToM

# prepare silhouette plot data
silhouettes <- c()
for (n in 1:length(networks)){
  sil <- silhouette(k2_data[[n]][[1]]@cluster, dmats[[n]])
  silhouettes <- append(silhouettes, list(sil))
}

# silhouette plots
plot(silhouettes[[1]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[1]))
plot(silhouettes[[2]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[2]))
plot(silhouettes[[3]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[3]))
plot(silhouettes[[4]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[4]))
plot(silhouettes[[5]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[5]))
plot(silhouettes[[6]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[6]))
plot(silhouettes[[7]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[7]))
plot(silhouettes[[8]], border=NA, col=c("#f7b148", "#8164cc"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[8]))

# prepare cluster plot data
tsnes <- c()
for (d in dmats){
  tsne <- Rtsne(d, is_distance=TRUE, exaggeration_factor = 10)
  tsnes <- append(tsnes, list(tsne))
}

# cluster plots
for (n in 1:8){
  lab1 <- paste0("1", " (N=", as.character(sum(clusd[[network_clusters[[n]]]] == 1)), ")")
  lab2 <- paste0("1", " (N=", as.character(sum(clusd[[network_clusters[[n]]]] == 2)), ")")
  tite <- paste0("t-SNE 2D Projections of k-means Clusters - ", anat_names[n])
  p <- tsnes[[n]]$Y %>% data.frame() %>% setNames(c("X", "Y")) %>% mutate(cluster = factor(k2_data[[n]][[1]]@cluster)) %>% 
    ggplot(aes(x = X, y = Y, colour = cluster)) +
    scale_color_manual(values=c("#f7b148", "#8164cc"), name="Cluster", labels=c(lab1, lab2)) +
    geom_point(size=3) +
    theme_light() +
    labs(title = tite)
  print(p)
}

# exploratory cluster plots
ex_clus <- c8[201:300]#network[cluster location] (1-100 = 2; 101-200 = 3, 201-300 = 4, 301-400 = 5, 401-500 = 6)
ex_clus <- ex_clus[[29]]#[[index of highest silhouette for cluster number]] View(sil_df) 

plot(ex_clus, type="series")
ex_sil <- silhouette(ex_clus@cluster, ex_clus@distmat)
plot(ex_sil, border=NA)
ex_tsne <- Rtsne(ex_clus@distmat, is_distance=TRUE)
ex_tsne$Y %>% data.frame() %>% setNames(c("X", "Y")) %>% mutate(cluster = factor(ex_clus@cluster)) %>% 
  ggplot(aes(x = X, y = Y, colour = cluster)) +
  #scale_color_manual(values=c("#f7b148", "#8164cc"), name="Cluster") +
  geom_point(size=3) +
  theme_light() +
  labs(title = 'Exploratory clustering')

# plot silhouettes + clusters for 3 & 4 clusters
# 3 clusters
silhouettes <- c()
for (n in 1:length(networks)){
  sil <- silhouette(k3_data[[n]][[1]]@cluster, dmats[[n]])
  silhouettes <- append(silhouettes, list(sil))
}
plot(silhouettes[[1]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[1]))
plot(silhouettes[[2]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[2]))
plot(silhouettes[[3]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[3]))
plot(silhouettes[[4]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[4]))
plot(silhouettes[[5]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[5]))
plot(silhouettes[[6]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[6]))
plot(silhouettes[[7]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[7]))
plot(silhouettes[[8]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[8]))
# cluster plots
for (n in 1:8){
  lab1 <- paste0("1", " (N=", as.character(sum(k3_data[[n]][[1]]@cluster == 1)), ")")
  lab2 <- paste0("1", " (N=", as.character(sum(k3_data[[n]][[1]]@cluster == 2)), ")")
  lab3 <- paste0("1", " (N=", as.character(sum(k3_data[[n]][[1]]@cluster == 3)), ")")
  tite <- paste0("t-SNE 2D Projections of k-means Clusters - ", anat_names[n])
  p <- tsnes[[n]]$Y %>% data.frame() %>% setNames(c("X", "Y")) %>% mutate(cluster = factor(k3_data[[n]][[1]]@cluster)) %>% 
    ggplot(aes(x = X, y = Y, colour = cluster)) +
    scale_color_manual(values=c("#f7b148", "#8164cc", "#80d9a5"), name="Cluster", labels=c(lab1, lab2, lab3)) +
    geom_point(size=3) +
    theme_light() +
    labs(title = tite)
  print(p)
}

# 4 clusters
silhouettes <- c()
for (n in 1:length(networks)){
  sil <- silhouette(k4_data[[n]][[1]]@cluster, dmats[[n]])
  silhouettes <- append(silhouettes, list(sil))
}
plot(silhouettes[[1]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[1]))
plot(silhouettes[[2]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[2]))
plot(silhouettes[[3]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[3]))
plot(silhouettes[[4]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[4]))
plot(silhouettes[[5]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[5]))
plot(silhouettes[[6]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[6]))
plot(silhouettes[[7]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[7]))
plot(silhouettes[[8]], border=NA, col=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), do.n.k=FALSE, do.clus.stat=FALSE, main=paste0("Silhouette plot - ", anat_names[8]))
# cluster plots
for (n in 1:8){
  lab1 <- paste0("1", " (N=", as.character(sum(k3_data[[n]][[1]]@cluster == 1)), ")")
  lab2 <- paste0("1", " (N=", as.character(sum(k3_data[[n]][[1]]@cluster == 2)), ")")
  lab3 <- paste0("1", " (N=", as.character(sum(k3_data[[n]][[1]]@cluster == 3)), ")")
  lab4 <- paste0("1", " (N=", as.character(sum(k4_data[[n]][[1]]@cluster == 4)), ")")
  tite <- paste0("t-SNE 2D Projections of k-means Clusters - ", anat_names[n])
  p <- tsnes[[n]]$Y %>% data.frame() %>% setNames(c("X", "Y")) %>% mutate(cluster = factor(k4_data[[n]][[1]]@cluster)) %>% 
    ggplot(aes(x = X, y = Y, colour = cluster)) +
    scale_color_manual(values=c("#f7b148", "#8164cc", "#80d9a5", "#d95b5b"), name="Cluster", labels=c(lab1, lab2, lab3, lab4)) +
    geom_point(size=3) +
    theme_light() +
    labs(title = tite)
  print(p)
}


# TIME SERIES PLOTS -------------------------------------------------------

##plot time courses (full)

# raw
ggplot(data=dfs_plot[[1]], aes(x=timepoint, y=raw, color=factor(id))) +
  geom_line(alpha=0.2) +
  #geom_line(aes(x=tp, y=stand), alpha=0.2) +
  #stat_smooth(geom="line", alpha=0.5, size=1, method="loess", se=FALSE, span=.1, n=245, fullrange=TRUE) +
  ylim(-5,5) + 
  guides(color="none") + #legend
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(y="BOLD Response", x="Time (volumes)", title="Raw") +
  #scale_x_discrete(breaks = tick_func(n=20), labels=seq(1,245,20)-1) +
  theme_light()

# standardized
ggplot(data=dfs_plot[[1]], aes(x=timepoint, y=stand, color=factor(id))) +
  geom_line(alpha=0.2) +
  #geom_line(aes(x=tp, y=stand), alpha=0.2) +
  #stat_smooth(geom="line", alpha=0.5, size=1, method="loess", se=FALSE, span=.1, n=245, fullrange=TRUE) +
  ylim(-5,5) + 
  guides(color="none") + #legend
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(y="BOLD Response", x="Time (volumes)", title="Standardized") +
  #scale_x_discrete(breaks = tick_func(n=20), labels=seq(1,245,20)-1) +
  theme_light()

#smoothed
ggplot(data=dfs_plot[[1]], aes(x=timepoint, y=smooth, color=factor(id))) +
  geom_line(alpha=0.5) +
  #geom_line(aes(x=timepoint, y=stand), alpha=0.05) +
  stat_smooth(geom="line", alpha=0.5, linewidth=1, method="loess", se=FALSE, span=.1, n=ttps, fullrange=TRUE, color="black") +
  ylim(-2,2) + 
  guides(color="none") + #legend
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(y="BOLD Response", x="Time (volumes)", title=paste0("Smoothed (bw = ", bw, ")")) +
  #scale_x_discrete(breaks = tick_func(n=20), labels=seq(1,245,20)-1) +
  theme_light()

# ISC PLOTS ---------------------------------------------------------------

##plot ISC heatmap

net_n = 8 ## SET NETWORK

netcluscol = network_clusters[[net_n]]
network = networks[[net_n]]
#sort heatmap data by variable
sorted_dat <- corrdat[order(corrdat[netcluscol]),] # sort variable
sorted_ids <- (unique(sorted_dat$sub_1))
#get cluster counts for dotted lines
kcount <- count(dat, dat[netcluscol])
k1count <- kcount[[2]][1]
ggplot(data=corrdat, aes(x=sub_1, y=sub_2, fill=.data[[network]])) + ## SET NETWORK ISCS
  geom_tile() +
  #scale_fill_gradientn(colours = c("#0c0078", "#2b02a6", "#0000c7", "#2855bf", "#0058fa","#3a90e0","#349cbf", "#00bfa6","#72c758", "#6fc238","#aab539","#bcc429","#d1c738","#ffe100"),
  ##scale_fill_gradientn(colours = c("#0c0078", "#0000c7", "#0058fa", "#349cbf", "#00bfa6","#3ac238","#d1c738","#ffe100", "#ffeb52","#eddf6f","#fff9e8"),
  #                     limits=c(-.6, 1), breaks=seq(-.6,1,by=.30)) +
  scale_fill_distiller(palette = "RdBu")+
  guides(fill = guide_colourbar(title= "ISC value", barwidth = 0.5, barheight = 10))+
  scale_x_discrete(limits=sorted_ids, breaks = tick_func(n=20), labels=seq(from=1,to=nps, by=20)) +
  scale_y_discrete(limits=sorted_ids, breaks = tick_func(n=20), labels=seq(from=1,to=nps, by=20)) +
  theme( 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle=55, hjust=1),
    #axis.text.y = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.ticks.y = element_blank(),
    panel.background = element_blank())+
  geom_vline(xintercept=k1count, linetype=2, linewidth=1) +
  geom_hline(yintercept=k1count, linetype=2, linewidth=1) +
  #geom_abline(intercept = 0.4, slope=1, size=1.5, color= "#fdatf32", linetype="solid") +
  coord_fixed()  #keeps tiles square

## plot isc heatmap (with distance matrix)

net_n = 6 ## SET NETWORK

netcluscol = network_clusters[[net_n]]
network = networks[[net_n]]
cldist = cldists[[net_n]]
#sort heatmap data by variable
sorted_dist <- dmats_plot[order(dmats_plot[netcluscol], dmats_plot[cldist]),]#tc_behav[order(tc_behav["age"]),]# # sort variable  
sorted_dids <- (unique(sorted_dist$sub_1))#p_id))#$sub_1))
#get cluster counts for dotted lines
kcount <- count(dat, dat[netcluscol])
k1count <- kcount[[2]][1]
distdat <- (scales::rescale(dmats_plot[[network]], to=c(1,-1)))

ggplot(data=dmats_plot, aes(x=sub_1, y=sub_2, fill=distdat)) + ## SET NETWORK ISCS
  geom_tile() +
  #scale_fill_gradientn(colours = c("#42059e", "#a878f0", "#c6a6f5", "#dfcdfa", "#eee6fa", "#f3f0f7", "#f3f0f7"))+
  ##scale_fill_gradientn(colours = c("#0c0078", "#0000c7", "#0058fa", "#349cbf", "#00bfa6","#3ac238","#d1c738","#ffe100", "#ffeb52","#eddf6f","#fff9e8"),
  #                     limits=c(-.6, 1), breaks=seq(-.6,1,by=.30)) +
  scale_fill_distiller(palette = "Spectral", direction=-1)+
  guides(fill = guide_colourbar(title= "Similarity", barwidth = 0.5, barheight = 10))+
  scale_x_discrete(limits=sorted_dids, breaks = tick_func(n=20), labels=seq(from=1,to=nps, by=20)) +
  scale_y_discrete(limits=sorted_dids, breaks = tick_func(n=20), labels=seq(from=1,to=nps, by=20)) +
  theme( 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle=55, hjust=1),
    #axis.text.y = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.ticks.y = element_blank(),
    panel.background = element_blank())+
  geom_vline(xintercept=k1count, linetype=2, linewidth=1) +
  geom_hline(yintercept=k1count, linetype=2, linewidth=1) +
  #geom_abline(intercept = 0.4, slope=1, size=1.5, color= "#fdatf32", linetype="solid") +
  coord_fixed()  #keeps tiles square

# BEHAVIOUR + CLUSTER PLOTS *** -----------------------------------------------

##UNFINISHED##

##plot behavioural violins -- CHANGE TO VIOLIN HALVES

violin_df <- scale(dat[12:34])
violin_df <- rescale(violin_df, mean=100, df=TRUE)# re-scaled to m=100 -- age not included
violin_df <- mutate(dat, violin_df)

#violin_df <- violin_df %>%
#  pivot_longer(cols=all_of(c("ctopp_elision", "nih_card", "wisc_psi", "wisc_vsi")), 
#               names_to="iv", values_to = "val")
#recombine with other data

ggplot(violin_df, aes(x=iv, y=val, fill=as.factor(net_8_clus)))+
  geom_violin(position=position_dodge(.7), alpha=0.5) +
  #geom_point(aes(color=as.factor(k8)), position = position_jitter(seed = 1, width = 0.2), alpha=0.7)+
  geom_boxplot(aes(color=as.factor(net_8_clus)), position= position_dodge(.7), width=.2, outlier.shape=NA, show.legend = FALSE)+
  scale_fill_manual(values=c("#f7b148", "#8164cc"), labels = c(1,2), name="Cluster")+
  scale_color_manual(values=c("grey15","grey15"))+
  #scale_x_discrete(labels=c("Language", "Phonology", "Social", "SRS(RRB)", "WISC(FSIQ)", "WISC(WMI)"))+
  theme_light(base_size = 15)+
  theme(axis.text.x = element_text(angle = 0))+
  labs(y="Score", x="Test")

##raincloud plots -- behavioural scores by cluster

netcluscol <- network_clusters[8]
ggplot(violin_df, aes(x=age, group=.data[[netcluscol]], fill=as.factor(.data[[netcluscol]])))+
  geom_density(adjust=1.5, alpha=0.5) + 
  scale_fill_manual(values=c("#f7b148", "#8164cc"), name="Cluster")+
  theme_bw()

# WRITE DATA ------------------------------------------------------------------

write.csv(dat[c(1, 10:56)], "./data/tc_behav_clusters_k4.csv", row.names=FALSE)
save(k2_data, file = "./data/k2_data_sdtw_g10.RData", envir=.GlobalEnv)
save(k3_data, file = "./data/k3_data_sdtw_g10.RData", envir=.GlobalEnv)
save(k4_data, file = "./data/k4_data_sdtw_g10.RData", envir=.GlobalEnv)

save(cs, file = "./data/cluster_data_sdtw.RData", envir=.GlobalEnv)

write.csv(dmats_save[[1]], "./data/dmat_1_sdtw_g10.csv")
write.csv(dmats_save[[2]], "./data/dmat_2_sdtw_g10.csv")
write.csv(dmats_save[[3]], "./data/dmat_3_sdtw_g10.csv")
write.csv(dmats_save[[4]], "./data/dmat_4_sdtw_g10.csv")
write.csv(dmats_save[[5]], "./data/dmat_5_sdtw_g10.csv")
write.csv(dmats_save[[6]], "./data/dmat_6_sdtw_g10.csv")
write.csv(dmats_save[[7]], "./data/dmat_7_sdtw_g10.csv")
write.csv(dmats_save[[8]], "./data/dmat_8_sdtw_g10.csv")

scaleddats <- c()
for (n in networks){
  scaleddat <- (scales::rescale(dmats_plot[[n]], to=c(1,-1)))
  #coerce back to matrix
  print(length(scaleddat))
  scaleddats <- append(scaleddats, list())
}

