# Version 3: updated for CBS server; final dataset

# INSTALL/CALL PACKAGES ---------------------------------------------------------------------------------

.libPaths(new="./other/r_libs")

library(easypackages)
libraries("tidyverse", "janitor", "naniar")

# IMPORT DATA ---------------------------------------------------------------------------------

filenames <- list.files("./experiments/TP-ISC/data/0724_data_pp", pattern="*.csv", full.names = TRUE)
ldf <- lapply(filenames, read.csv, header=FALSE)

# TIDY DATA -----------------------------------------------------------------------------------

rois = paste("roi", c(seq(1,200,1)), sep="_")
nets = c("net_1","net_2","net_3","net_4","net_5","net_6","net_7","net_8")
ids = lapply(filenames, gsub, pattern = "./experiments/TP-ISC/data/0724_data_pp/", replacement="")
ids = lapply(ids, gsub, pattern = "_pp.csv", replacement="")

#turn columns of signal at each time point into list of signal at all time points
rows2vec_func <- function(dat){
  for (n in seq(2,209,1)){
    print(n)
    #test <- c(dat[[n]])
    dat[[n]] <- list(dat[[n]])
  }
  dat <- filter(dat, dat$V1==1) #trim all but first row (duplicates for some reason)
  return(dat)
}

dfc <- lapply(ldf, rows2vec_func)

#combine participants into single df
d <- as.data.frame(do.call(rbind, dfc))
d$V1 <- ids
colnames(d) <- c("p_id", nets, rois)

# WRITE DATA ----------------------------------------------------------------------------------

#d = data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
d %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = '|')) %>% 
  write.csv('./data/0724_timecourse_all.csv', row.names = FALSE)
#write.csv(d, "./data/MRI-timecourse.csv", row.names=FALSE)

d[1:9] %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = '|')) %>% 
  write.csv('./data/0724_timecourse_nets1-8.csv', row.names = FALSE)

d[1] %>%
  rowwise() %>%
  mutate_if(is.list, ~paste(unlist(.), collapse = '|')) %>% 
  write.csv('./data/0724_ids_final.csv', row.names=FALSE)


