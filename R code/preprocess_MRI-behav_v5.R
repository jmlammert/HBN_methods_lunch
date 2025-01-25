## VERSION 5: celf age norm (total-criterion); added Barratt, composites

# INSTALL/CALL PACKAGES ---------------------------------------------------------------------------------

.libPaths(new="./other/r_libs")

library(easypackages)
libraries("tidyverse", "janitor", "naniar")

# IMPORT DATA ---------------------------------------------------------------------------------

# MRI tracker
MRI_track <- clean_names(read_csv("./other/MRI_track.csv"))
# demographics
raw_demos <- clean_names(read_csv("./experiments/behavioural-assessments/data/demographics/basic_demos.csv"))
raw_barratt <- clean_names(read_csv("./experiments/behavioural-assessments/data/demographics/Barratt.csv"))
# WISC
raw_WISC <- clean_names(read_csv("./experiments/behavioural-assessments/data/WISC/WISC.csv"))
# WISC-remote
raw_WISC_r <- clean_names(read_csv("./experiments/behavioural-assessments/data/WISC/WISC_remote.csv"))
# SRS
raw_SRS <- clean_names(read_csv("./experiments/behavioural-assessments/data/SRS/SRS.csv"))
# language tests
raw_CELF <- clean_names(read_csv("./experiments/behavioural-assessments/data/CELF/CELF.csv"))
raw_CTOPP <- clean_names(read_csv("./experiments/behavioural-assessments/data/CTOPP/CTOPP.csv"))
raw_GFTA <- clean_names(read_csv("./experiments/behavioural-assessments/data/GFTA/GFTA.csv"))
raw_TOWRE <- clean_names(read_csv("./experiments/behavioural-assessments/data/TOWRE/TOWRE.csv"))
# NIH
raw_NIH <- clean_names(read_csv("./experiments/behavioural-assessments/data/NIH/NIH.csv"))
# consensus diagnoses
raw_dx <- data.frame(clean_names(read_csv("./experiments/behavioural-assessments/data/diagnoses/diagnoses.csv")))

# TIDY DATA -----------------------------------------------------------------------------------

# remove incomplete from ##WISC_remote
complete_WISC_r <- filter(raw_WISC_r, wisc_remote_wisc_complete == 1)
complete_WISC <- filter(raw_WISC, wisc_wisc_complete == 1)

# combine data by id
raw_all <- raw_demos # make copy of demos for master dataframe
raw_all <- left_join(raw_all, complete_WISC, by=c('basic_demos_eid'='wisc_eid')) %>%
  left_join(., raw_barratt, by=c('basic_demos_eid'='barratt_eid')) %>%
  left_join(., raw_NIH, by=c('basic_demos_eid'='nih_scores_eid')) %>%
  left_join(., raw_SRS, by=c('basic_demos_eid'='srs_eid')) %>%
  left_join(., raw_CELF, by=c('basic_demos_eid'='celf_eid')) %>%
  left_join(., raw_CTOPP, by=c('basic_demos_eid'='ctopp_eid')) %>%
  left_join(., raw_GFTA, by=c('basic_demos_eid'='gfta_eid')) %>%
  left_join(., raw_TOWRE, by=c('basic_demos_eid'='towre_eid')) %>%
  left_join(., MRI_track, by=c('basic_demos_eid'='mri_track_eid')) %>%
  left_join(., complete_WISC_r, by=c('basic_demos_eid'='wisc_remote_eid')) ##WISC_remote (only complete)

# filter for complete/valid behavioural tests and qualified/The Present MRI
complete_all <- filter(raw_all, celf_celf_valid == 1  & ctopp_ctopp_valid == 1 
                       & gfta_gfta_siw_valid == 1 & towre_towre_valid == 1 & mri_track_movie2 == 1 
                       & mri_track_qualify == 1 & nih_scores_nih7_complete == 1, )##WISC_remote and WISC accounted for above
##SRS does not have valid field -- NAs removed below
                        

# remove participants who don't actually have MRI scan (based off data availability on server)
##MRI_IDs <- clean_names(read_csv("./data/MRI-behav_IDs-final.csv")) # import IDs of participants who have data
final_all <- complete_all
##final_all <- left_join(final_all, MRI_IDs, by=c('basic_demos_eid'='id'))
##final_all <- filter(final_all, final == 1)

# select columns/variables for analysis ##
final_tasks <- final_all %>%
  select(basic_demos_eid, basic_demos_sex, basic_demos_age, barratt_barratt_total, celf_celf_criterion_score, celf_celf_total, ctopp_ctopp_el_s, ctopp_ctopp_bw_s,
         ctopp_ctopp_nr_s, ctopp_ctopp_rsn_comp, gfta_gfta_siw_stnd, towre_towre_swe_scaled,
         towre_towre_pde_scaled, wisc_wisc_fsiq, wisc_remote_wisc_nmfsiq, wisc_wisc_vsi, wisc_remote_wisc_nmvsi,
         wisc_wisc_vci, wisc_remote_wisc_vci, wisc_wisc_fri, wisc_remote_wisc_fri, wisc_wisc_wmi, wisc_remote_wisc_wmi,
         wisc_wisc_psi, wisc_remote_wisc_nsi, srs_srs_awr_t, srs_srs_cog_t, srs_srs_com_t, srs_srs_mot_t, srs_srs_rrb_t,
         nih_scores_nih7_card, nih_scores_nih7_flanker)

##calculate celf age norm
final_tasks["celf_celf_total"] <- final_tasks["celf_celf_total"] - final_tasks["celf_celf_criterion_score"]
final_tasks <- subset(final_tasks, select=-c(celf_celf_criterion_score)) #drop criterion column

# clean data
final_clean <- data.frame(final_tasks)
# change variable names
colnames(final_clean) <- c("id", "sex", "age", "barratt", "celf_total", "ctopp_elision", "ctopp_blending", "ctopp_nwrep", 
                           "ctopp_symbol", "gfta_siw", "towre_swe", "towre_pde", "wisc_fsiq", "wisc_r_fsiq", 
                           "wisc_vsi", "wisc_r_vsi", "wisc_vci", "wisc_r_vci", "wisc_fri", "wisc_r_fri", 
                           "wisc_wmi", "wisc_r_wmi", "wisc_psi", "wisc_r_psi", "srs_awr", "srs_cog", "srs_com",
                           "srs_mot", "srs_rrb", "nih_card", "nih_flanker")
#enforce data types
final_clean$id <- as.character(final_clean$id) # ID as character
final_clean$sex <- factor(as.numeric(final_clean$sex), levels = c(1,0), labels = c("female", "male")) # sex as factor
num_cols <- c("age", "barratt", "celf_total", "ctopp_elision", "ctopp_blending", "ctopp_nwrep", "ctopp_symbol",
              "gfta_siw", "towre_swe", "towre_pde", "wisc_fsiq", "wisc_r_fsiq", "wisc_vsi",
              "wisc_r_vsi", "wisc_vci", "wisc_r_vci", "wisc_fri", "wisc_r_fri", "wisc_wmi", "wisc_r_wmi",
              "wisc_psi", "wisc_r_psi", "srs_awr", "srs_cog", "srs_com", "srs_mot", "srs_rrb", "nih_card", "nih_flanker") # list of columns to be numeric
final_clean[num_cols] <- sapply(final_clean[num_cols], as.numeric) # all other continuous variables as numeric
# check NAs
na_all <- colSums(is.na(final_clean)) # only "WISC remote" should have NAs (completed regular WISC)
print(na_all) 
## some additional NAs removed below

# remote WISC
# if valid, must combine with non-remote WISC to get single WISC score
# since none valid (10/02/2023), remove columns
final_clean <- final_clean %>%
  select(!c(wisc_r_fsiq, wisc_r_vsi, wisc_r_vci, wisc_r_fri, wisc_r_wmi, wisc_r_psi))

# remove remaining NAs
na_final <- colSums(is.na(final_clean)) # check NAs: few in ctopp, gfta, towre, wisc psi (1-11 missing)
print(na_final)
final_clean <- na.omit(final_clean) # remove rows
na_final <- colSums(is.na(final_clean)) # confirm clean
print(na_final)

# -----COMPOSITES-----
final_clean$iv_lang <- rowMeans(rescale(scale(final_clean[,c("celf_total", "ctopp_elision", "ctopp_blending", "ctopp_nwrep", "ctopp_symbol",
                                   "gfta_siw", "towre_swe", "towre_pde", "wisc_vci")]), mean=100, df=TRUE))
final_clean$iv_phon <- rowMeans(rescale(scale(final_clean[,c("ctopp_elision", "ctopp_blending", "ctopp_nwrep", "ctopp_symbol", "towre_pde")]), mean=100, df=TRUE))
final_clean$iv_soco <- rowMeans(rescale(scale(final_clean[,c("srs_awr", "srs_cog", "srs_com", "srs_mot")]), mean=100, df=TRUE))
final_clean$iv_asd <- rowMeans(rescale(scale(final_clean[,c("srs_awr", "srs_cog", "srs_com", "srs_mot", "srs_rrb")]), mean=100, df=TRUE))
final_clean$iv_attn <- rowMeans(rescale(scale(final_clean[,c("nih_card", "nih_flanker")]), mean=100, df=TRUE))

# -----DIAGNOSES-----

# tag participants with IQ < 70
final_clean["iq_below_70"] <- ifelse(final_clean$wisc_fsiq < 70, 1, 0) # assigns 1 = IQ<70, 0 = IQ>70; ~35
# remove participants with IQ < 70
#final_clean <- filter(final_clean, final_clean$wisc_fsiq > 70) # ~1554

# make long (all diagnoses in 1 column)
# diagnosis
dx_long <- pivot_longer(raw_dx, cols = c("diagnosis_clinician_consensus_dx_01", "diagnosis_clinician_consensus_dx_02",
                                         "diagnosis_clinician_consensus_dx_03", "diagnosis_clinician_consensus_dx_04",
                                         "diagnosis_clinician_consensus_dx_05", "diagnosis_clinician_consensus_dx_06",
                                         "diagnosis_clinician_consensus_dx_07", "diagnosis_clinician_consensus_dx_08",
                                         "diagnosis_clinician_consensus_dx_09", "diagnosis_clinician_consensus_dx_10"),
                        names_to = "n_dx", values_to = "dx")
#dx_long <- dx_long %>% select("n_dx", "dx")
# confirmed
dx_conf <- pivot_longer(raw_dx, cols = c("diagnosis_clinician_consensus_dx_01_confirmed", "diagnosis_clinician_consensus_dx_02_confirmed",
                                         "diagnosis_clinician_consensus_dx_03_confirmed", "diagnosis_clinician_consensus_dx_04_confirmed",
                                         "diagnosis_clinician_consensus_dx_05_confirmed", "diagnosis_clinician_consensus_dx_06_confirmed",
                                         "diagnosis_clinician_consensus_dx_07_confirmed", "diagnosis_clinician_consensus_dx_08_confirmed",
                                         "diagnosis_clinician_consensus_dx_09_confirmed", "diagnosis_clinician_consensus_dx_10_confirmed"),
                        names_to = "n_dx", values_to = "confirmed")
#dx_conf <- dx_conf %>% select("n_dx", "confirmed")

# isolate ID and diagnosis and confirmation columns
dx <- dx_long %>% 
  select(identifiers, dx)
dx_conf <- dx_conf %>%
  select(identifiers, confirmed)
# add confirm column
dx <- cbind(dx, dx_conf$confirmed)
# check all diagnoses for keywords -- reduces to only diagnoses of interest
dx <- dx %>%
  filter(if_any(everything(), ~ grepl("Autism Spectrum Disorder|ADHD|Attention-Deficit/Hyperactivity Disorder|Language Disorder|Speech Sound Disorder|Impairment in Reading", .)))
# only include confirmed diagnoses
dx <- dx %>%
  filter(`dx_conf$confirmed` == 1)
# make wide (each diagnosis in own column, 1 row per participant)
dx_wide <- pivot_wider(dx, names_from=dx, values_from=`dx_conf$confirmed`)

# clean diagnosis dataframe
dx_clean <- data.frame(sapply(dx_wide, function(x) {
  gsub(",assessment", "", x)
})) # remove ",assessment" from ID
dx_clean <- clean_names(dx_clean) # clean column names
dx_clean[1:9] <- sapply(dx_clean[1:9], as.character)
dx_clean <- dx_clean %>% replace(.=="NULL", 0) # replace nulls with 0
dx_clean[2:9] <- sapply(dx_clean[2:9], as.numeric) # make columns numeric

# subset participants with diagnoses to those who have completed tasks (~903)
dx_tasks <- dx_clean[dx_clean$identifiers %in% final_clean$id, ]

# remove participants with diagnoses who have completed tasks from all participants who completed tasks (to get TD)
final_td <- final_clean[! final_clean$id %in% dx_tasks$identifiers, ]
final_td <- final_td %>%
  filter(iq_below_70 == 0)# remove IQ<70
final_td$num_dx <- 0

# combine diagnosis data and tasks
final_dx <- left_join(final_clean, dx_tasks, by=c('id'='identifiers'))
final_dx[is.na(final_dx)] <- 0 # replace NA with 0
# create number of diagnoses column to ID participants with/without
final_dx$num_dx <- apply(final_dx[,c(31:length(final_dx))], 1, sum) #indices of diagnosis columns

# WRITE DATA ----------------------------------------------------------------------------------

write.csv(final_all, "./data/MRI-behav_raw-final.csv", row.names=FALSE)
write.csv(final_dx, "./data/0824_MRI-behav_comp.csv", row.names=FALSE)
write.csv(final_td, "./data/MRI-behav_TD.csv", row.names=FALSE)


