## Version 3: updated for composites, barratt score -- see Emiko code for APA table

# INSTALL/CALL PACKAGES ---------------------------------------------------------------------------------

.libPaths(new="./other/r_libs")

library(easypackages)
libraries("tidyverse", "janitor", "naniar", "psych", "reshape2")

# IMPORT DATA ---------------------------------------------------------------------------------

df <- clean_names(read_csv("./data/0824_MRI-behav_comp.csv")) 
ids <- clean_names(read_csv("./data/0724_ids_final.csv"))

# TIDY DATA -----------------------------------------------------------------------------------

# ensure no NAs
colSums(is.na(df))
#subset to final dataset
df <- df[df$id %in% ids$p_id, ]

# standardize 
df_z <- scale(df[4:25])# z-scale, all between +-3, m=0, sd=1
df_s <- rescale(df[4:25], mean=100, df=TRUE)# re-scaled to m=100 -- age not included
# recombine with id, sex, and age + composites
df_z_all <- cbind(df[1:3], df_z, df[26:30])
df_s_all <- cbind(df[1:3], df_s, df[26:30])

#separate dx data
dx <- df[31:40]
dx <- cbind(df$id, dx)
dx <- rename(dx, c('id'='df$id'))

#create dx composites
dx$adhd_overall <- with(dx, ifelse(adhd_inattentive_type ==1 | adhd_combined_type == 1| 
                                     other_specified_attention_deficit_hyperactivity_disorder == 1 |
                                     adhd_hyperactive_impulsive_type == 1, 1, 0))
dx$td <- with(dx, ifelse(num_dx == 0, 1, 0))
dx$mult_dx <- with(dx, ifelse(num_dx > 1, 1, 0))

#removed multiple dx from plot
dx_plot <- data.frame(dx = c("iq_below_70", "adhd_inattentive_type", "specific_learning_disorder_with_impairment_in_reading", "adhd_combined_type", "other_specified_attention_deficit_hyperactivity_disorder", "language_disorder", "speech_sound_disorder", "autism_spectrum_disorder", "adhd_hyperactive_impulsive_type", "td", "mult_dx"),
                      count = c(sum(dx$iq_below_70), sum(dx$adhd_inattentive_type), sum(dx$specific_learning_disorder_with_impairment_in_reading), sum(dx$adhd_combined_type), sum(dx$other_specified_attention_deficit_hyperactivity_disorder), sum(dx$language_disorder), sum(dx$speech_sound_disorder), sum(dx$autism_spectrum_disorder), sum(dx$adhd_hyperactive_impulsive_type), sum(dx$td), sum(dx$mult_dx)))
dx_plot$dx <- factor(dx_plot$dx, levels=c("iq_below_70", "adhd_inattentive_type", "adhd_hyperactive_impulsive_type", "adhd_combined_type", "other_specified_attention_deficit_hyperactivity_disorder", "specific_learning_disorder_with_impairment_in_reading", "language_disorder", "speech_sound_disorder", "autism_spectrum_disorder", "td", "mult_dx"),
                     labels=c("IQ < 70", "ADHD (Inattentive)", "ADHD (Hyper.-Impulsive)", "ADHD (Combined)", "ADHD (Other specified)", "Reading disorder", "Language disorder", "SSD", "ASD", "No diagnosis", "Multiple diagnoses"))
  
#create age data
df_age <- data.frame(round(df$age, digits=0))
colnames(df_age) <- "age_round"
#set age bins -- currently 4 bins of size 3
df_age <- df_age %>% mutate(age_bin=case_when(age_round == 6 | age_round == 7 | age_round == 8  ~ "6-8",
                                              age_round == 9 | age_round == 10 | age_round == 11 ~ "9-11",
                                              age_round == 12 | age_round == 13 | age_round == 14 ~ "12-14",
                                             age_round == 15 | age_round == 16 | age_round == 17 ~ "15-17",))

#age_count <- df_age %>% count(age_bin, name="count")
#age_count$age_bin <- factor(age_count$age_bin, levels=c("6-8", "9-11", "12-14", "15-17"))

#updated age counts (all ages)
age_count <- data.frame(table(df_age$age_round))
colnames(age_count) <- c("age", "count")

# DESCRIPTIVE STATS ---------------------------------------------------------------------------

summary_all <- describeBy(df)
summary_sex <- describeBy(df, df$sex) # grouping variable
summary_sex_s <- describeBy(df_s_all, df_s_all$sex) # standardized

# TESTS OF ASSUMPTIONS ------------------------------------------------------------------------

# Shapiro-Wilk test for normality if p>0.0.5 = normal distribution
norm_tests <- data.frame() # empty df for norm tests

for (col in 3:30) { # for each of subset of columns in df
  col_name <- colnames(df[col]) # get column/variable name
  if (length(unique(df[[col]])) > 1) { # Check if there is variability in the data
    sw <- shapiro.test(df[[col]]) # calculate Shapiro-Wilk test
    sw_stat <- sw$statistic # get SW statistic
    p_val <- sw$p.value # get SW p-value
    if (p_val <= 0.05) { #p-value < 0.05? code significant
      sig <- TRUE
      norm <- "non-norm"
    } else {
      sig <- FALSE
      norm <- "norm"
    }
  } else { # If all values are identical, set test results to NA
    sw_stat <- NA
    p_val <- NA
    sig <- NA
  }
  norm_test <- data.frame(variable = col_name, stat = sw_stat, p_val = p_val, sig = sig, norm = norm) # assemble into data frame
  norm_tests <- rbind(norm_tests, norm_test) # add to df
}
colnames(norm_tests) <- c("variable", "stat", "p_val", "sig", "norm")
rownames(norm_tests) <- NULL # re-index

# CORRELATIONS --------------------------------------------------------------------------------

cors <- data.frame(cor(df[3:30], method="pearson"))

# FIGURES -------------------------------------------------------------------------------------

#list of behavioural test labels
test_labs = c("age"="Age", "barratt"="SES", "celf_total"="CELF", "ctopp_elision"="CTOPP(EL)",
              "ctopp_blending"="CTOPP(BL)", "ctopp_nwrep"="CTOPP(NWR)","ctopp_symbol"="CTOPP(SYM)",
              "gfta_siw"="GFTA", "towre_swe"="TOWRE(SWE)", "towre_pde"="TOWRE(PDE)", "wisc_fsiq"="WISC(IQ)",
              "wisc_vsi"="WISC(VS)", "wisc_vci"="WISC(VC)", "wisc_fri"="WISC(FR)", "wisc_wmi"="WISC(WM)", 
              "wisc_psi"="WISC(PS)", "srs_awr"="SRS(AWR)", "srs_cog"="SRS(COG)", "srs_com"="SRS(COM)", 
              "srs_mot"="SRS(MOT)", "srs_rrb"="SRS(RRB)", "nih_card" = "DCCS", "nih_flanker"="Flanker", 
              "iv_lang"="Language", "iv_phon"="Phonology", "iv_soco"="Communication", "iv_asd"="Autistic Traits",
              "iv_attn"="Attention")

# DISTRIBUTION PLOTS ------------------------------------------------------

# distributions of continuous variables
all_dist <- multi.hist(df[3:30], global=FALSE, freq=TRUE, density = FALSE, col="#ea9999", main=test_labs, mar=c(2,2,1.5,1))  

# see how scaling changes (age remains)
all_dist_z <- multi.hist(df_z_all[3:30], global=FALSE, col="#e0d063",freq=TRUE, density = FALSE, main=test_labs, mar=c(2,2,1.5,1))
all_dist_s <- multi.hist(df_s_all[3:30], nrow=4, global=FALSE, col="#A1C7B1",freq=TRUE, density = FALSE, main=test_labs, mar=c(2,2,1.5,1))

# age histogram
ggplot(df, aes(x=age)) +
  geom_histogram(fill="royalblue2", bins=20) + 
  labs(x="Age (years)", y = "Count") +
  theme_classic()

# age ranges
ggplot(age_count, aes(x=age, y=count)) +
  geom_bar(stat="identity", fill="royalblue2") +
  labs(x="Age (years)", y = "Count") +
  #coord_flip() +
  theme_classic()

# DIAGNOSIS PLOTS ---------------------------------------------------------

#pie chart
ggplot(dx_plot[1:10,], aes(x="", y=count, fill=dx)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  theme_void()

#bar graph
ggplot(dx_plot[1:10,], aes(x=reorder(dx,-count), y=count)) + #removed multiple diagnoses
  geom_bar(stat="identity", fill="gray30") +
  labs(x="Diagnosis", y = "Count") +
  coord_flip() +
  theme_light()

# CORRELATION HEATMAP -----------------------------------------------------

#prepare data for heatmap -- individual tasks
melted_cors <- melt(as.matrix(cors[1:23, 1:23]))

#plot heatmap -- individual tasks
ggplot(data = melted_cors, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#0000FF", mid = "#FFFFCC", high = "#FF0000",
                       limits=c(-1, 1), breaks=seq(-1,1,by=0.25)) +
  #scale_fill_continuous(limits=c(-1, 1), breaks=seq(-1,1,by=0.25))+
  guides(fill = guide_colourbar(title= "Pearson  
correlation", barwidth = 0.5, barheight = 20))+ #spacing for formatting
  scale_x_discrete(labels=test_labs) +
  scale_y_discrete(labels=test_labs) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle=55, hjust=1),
    panel.background = element_blank())+
  coord_fixed()  #keeps tiles square

#prepare data for heatmap -- composites
melted_comp <- melt(as.matrix(cors[c(1:2, 11, 15, 24:28),c(1:2, 11, 15, 24:28)]))

#plot heatmap -- composites
ggplot(data = melted_comp, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#0000FF", mid = "#FFFFCC", high = "#FF0000",
                       limits=c(-1, 1), breaks=seq(-1,1,by=0.25)) +
  #scale_fill_continuous(limits=c(-1, 1), breaks=seq(-1,1,by=0.25))+
  guides(fill = guide_colourbar(title= "Pearson  
correlation", barwidth = 0.5, barheight = 15))+ #spacing for formatting
  scale_x_discrete(labels=test_labs) +
  scale_y_discrete(labels=test_labs) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle=55, hjust=1),
    panel.background = element_blank())+
  coord_fixed()  #keeps tiles square

# WRITE DATA ------------------------------------------------------------------

write.csv(df, "./data/0824_MRI-behav_comp-final.csv", row.names=FALSE)

