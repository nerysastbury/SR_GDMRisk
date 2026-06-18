library(readxl)

library(meta)

library(dplyr)

BCOFanalysisr0.8 <- read_excel("Redacted/Systematic Review for Behavioural Interventions/Analysis/BCOFanalysisr0.8xls")


##############################Triglycerides at 12 months
data_Triglycerides_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "Triglycerides (mmol/l)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_Triglycerides_12$study_label <- paste(data_Triglycerides_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_Triglycerides_12,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  studlab = data_Triglycerides_12$study_label,
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

funnel(meta_result)
metabias(meta_result, method.bias="linreg", k.min=5)


jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/tri 12 months R0.8 BCOF Random Effects HKSJ.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/l]", "95% CI [mmol/l]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm", 
       fontsize = 12)
dev.off()

##############################Total cholesterol
data_totchol_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "Total cholesterol (mmol/l)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))
data_totchol_12$study_label <- paste0(data_totchol_12$study, " (N = ", data_totchol_12$n, ")")

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_totchol_12,
  studlab = data_totchol_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)


funnel(meta_result)
metabias(meta_result, method.bias="linreg", k.min=5)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rtotal chol 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/l]", "95% CI [mmol/l]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()

##############################LDL cholesterol
data_LDL_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "Low-density lipoprotein cholesterol (mmol/l)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_LDL_12$study_label <- paste(data_LDL_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_LDL_12,
  studlab = data_LDL_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rLDL chol 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/l]", "95% CI [mmol/l]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",  
       fontsize = 12)
dev.off()
##############################SBP at 12 months
data_SBP_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "Systolic blood pressure (mmHg)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_SBP_12$study_label <- paste(data_SBP_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_SBP_12,
  studlab = data_SBP_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rSBP 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmHg]", "95% CI [mmHg]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()


##############################DBP at 12 months
data_DBP_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "Diastolic blood pressure (mmHg)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_DBP_12$study_label <- paste(data_DBP_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_DBP_12,
  studlab = data_DBP_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rDBP 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmHg]", "95% CI [mmHg]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()

##############################Weight at 12-16 months
data_weight_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "weight", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_weight_12$study_label <- paste(data_weight_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_weight_12,
  studlab = data_weight_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

funnel(meta_result)
metabias(meta_result, method.bias="linreg", k.min=5)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rWeight 12 months BCOF random effects HKSJ.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [kg]", "95% CI [kg]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()

#################################fpg
data_fpg_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "Fasting plasma glucose (mmol/l)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_fpg_12$study_label <- paste0(data_fpg_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_fpg_12,
  studlab = data_fpg_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rfpg R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/l]", "95% CI [mmol/l]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()


#################################ogtt2
data_ogtt2_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "2-hour Oral glucose tolerance test (mmol/l)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))
data_ogtt2_12$study_label <- paste0(data_ogtt2_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_ogtt2_12,
  studlab = data_ogtt2_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rogtt2 R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/l]", "95% CI [mmol/l]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()

funnel(meta_result)


#################################waist circumference at 12 months
data_wc_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "wc", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))
data_wc_12$study_label <- paste0(data_wc_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_wc_12,
  studlab = data_wc_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

funnel(meta_result)
metabias(meta_result, method.bias="linreg", k.min=5)


jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rwc 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [cm]", "95% CI [cm]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()


##########################################HbA1c (mmol/mol)
data_hba1c_12<- BCOFanalysisr0.8 %>%
  filter(bio == "HbA1c (mmol/mol)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))
data_hba1c_12$study_label <- paste0(data_hba1c_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_hba1c_12,
  studlab = data_hba1c_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/rhba1c 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/mol]", "95% CI [mmol/mol]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)

dev.off()


##########################################HbA1c (%)
data_hba1c_12<- BCOFanalysisr0.8 %>%
  filter(bio == "HbA1c (%)", months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))
data_hba1c_12$study_label <- paste0(data_hba1c_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_hba1c_12,
  studlab = data_hba1c_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/percent_hba1c 12 months R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [%]", "95% CI [%]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)

dev.off()



#################################t2dm rr at 12 months
library(glmmTMB)

library(metafor)

library(meta)

data_t2dm_12 <- BCOFanalysisr0.8 %>%
  filter(bio == "t2dm", months == 12, !study %in% c("Potzel"))
data_t2dm_12$study_label <- paste0(data_t2dm_12$study)

df_nozero <- data_t2dm_12 %>% filter(!study %in% c("Cheung, 2011"))
df_nozero$study_label <- paste0(df_nozero$study)

meta_nozero <- metabin(
  event.e = a, n.e = a + b,
  event.c = c, n.c = c + d,
  studlab = df_nozero$study_label,
  data = df_nozero,
  sm = "RR", method = "MH",  method.tau = "REML",  
  method.random.ci = "HK", 
  common = FALSE, random = TRUE, subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/t2dm_12_forest.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_nozero,
       leftcols = c("studlab","n.e","n.c", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study","N Int","N Ctrl", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [RR]", "95% CI [RR]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()


#####################now add continuity correction
meta_cc <- metabin(
  event.e = a, n.e = a + b,
  event.c = c, n.c = c + d,
  studlab = data_t2dm_12$study_label,
  data = data_t2dm_12,
  sm = "RR", method = "MH",  method.tau = "REML",  
  method.random.ci = "HK", 
  common = FALSE, random = TRUE,  subgroup = time, 
  print.subgroup.name = FALSE
)


jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/t2dm 12 months_ cc R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)

forest(meta_cc, backtransf = TRUE,
       leftcols = c("studlab","n.e","n.c", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study","N Int","N Ctrl", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [RR]", "95% CI [RR]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",   # more room for plot
       colgap.left = "0.7cm",   # more gap between left cols
       fontsize = 12)
dev.off()


####################################################COMBINING ALL LIPIDS 
data_lipids_12 <- BCOFanalysisr0.8 %>%
  filter(bio %in% c("Total cholesterol (mmol/l)", "High-density lipoprotein cholesterol (mmol/l)", "Low-density lipoprotein cholesterol (mmol/l)", "Triglycerides (mmol/l)"), months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))


data_lipids_12$study_label <- paste0(data_lipids_12$study)
data_lipids_12$bio <- factor(data_lipids_12$bio, levels = c("Total cholesterol (mmol/l)", "Low-density lipoprotein cholesterol (mmol/l)", "High-density lipoprotein cholesterol (mmol/l)", "Triglycerides (mmol/l)"))

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_lipids_12,
  studlab = data_lipids_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = bio, 
  overall = FALSE,
  print.subgroup.name = FALSE
)


jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/alllipids 12 months_ cc R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 12, units = "in", res=1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/l]", "95% CI [mmol/l]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",  
       fontsize = 12)
dev.off()

####################################################COMBINING ALL BP 
data_pressure_12 <- BCOFanalysisr0.8 %>%
  filter(bio %in% c("Systolic blood pressure (mmHg)", "Diastolic blood pressure (mmHg)"), months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_pressure_12$study_label <- paste0(data_pressure_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_pressure_12,
  studlab = data_pressure_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = bio, 
  overall = FALSE,
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/alllpressure 12 months_ cc R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res=1200)

forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmHg]", "95% CI [mmHg]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()


####################################################COMBINING ALL glycaemic outcomes
data_glyc_12 <- BCOFanalysisr0.8 %>%
  filter(bio %in% c("Fasting plasma glucose (mmol/l)", "2-hour Oral glucose tolerance test (mmol/l)", "HbA1c (mmol/mol)"), months == 12, time %in% c("Pregnancy", "Postpartum"), weightfoucsed %in% c("Weight focus", "Not weight focused"))

data_glyc_12$study_label <- paste0(data_glyc_12$study)

meta_result <- metagen(
  TE = mean_difference,
  seTE = mean_difference_se,
  data = data_glyc_12,
  studlab = data_glyc_12$study_label,
  sm = "MD",
  method.tau = "REML",  
  method.random.ci = "HK",    
  random = TRUE,
  common = FALSE,
  subgroup = bio, 
  overall = FALSE,
  print.subgroup.name = FALSE
)
jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/allglyc 12 months_ cc R0.8 BCOF Random Effects HK.jpeg", width = 15, height = 8, units = "in", res = 1200)
forest(meta_result, 
       leftcols = c("studlab","n", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study", "N", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [mmol/mol]", "95% CI [mmol/mol]", "Weight"), 
       label.right="favours control",
       label.left = "favours intervention")
dev.off()

#########################longest fup diabetes - a separate thing because I will use a different excel sheet to do this 
library(readxl)

library(meta)

library(dplyr)

longest_fup <- read_excel("Redacted/Systematic Review for Behavioural Interventions/Analysis/t2dm longest fup.xlsx")

library(glmmTMB)

library(metafor)

library(meta)

data_t2dm <- longest_fup %>%
  filter(bio == "t2dm")
data_t2dm$study_label <- paste0(data_t2dm$study)

df_nozero <- data_t2dm %>% filter(!study %in% c("Cheung, 2011"))
df_nozero$study_label <- paste0(df_nozero$study)

meta_nozero <- metabin(
  event.e = a, n.e = a + b,
  event.c = c, n.c = c + d,
  studlab = df_nozero$study_label,
  data = df_nozero,
  sm = "RR", method = "MH",  method.tau = "REML",  
  method.random.ci = "HK", 
  common = FALSE, random = TRUE, subgroup = time, 
  print.subgroup.name = FALSE
)

jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/t2dm_longestfup_forest.jpeg", width = 20, height = 6, units="in", res = 1200)
forest(meta_nozero,
       leftcols = c("studlab","n.e","n.c", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study","N Int","N Ctrl", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [RR]", "95% CI [RR]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",  
       colgap.left = "0.7cm",
       colgap.right = "0.7cm",   
       fontsize = 12)
dev.off()


#####################now add continuity correction
meta_cc <- metabin(
  event.e = a, n.e = a + b,
  event.c = c, n.c = c + d,
  studlab = data_t2dm$study_label,
  data = data_t2dm,
  sm = "RR", method = "MH",  method.tau = "REML",  
  method.random.ci = "HK", 
  common = FALSE, random = TRUE,  subgroup = time, 
  print.subgroup.name = FALSE
)


jpeg (file="Redacted/Systematic Review for Behavioural Interventions/Results/BOCF 0.8/t2dm longest fup months_ cc R0.8 BCOF Random Effects HK.jpeg", width = 20, height = 6, units="in", res = 1200)

forest(meta_cc, backtransf = TRUE,
       leftcols = c("studlab","n.e","n.c", "population", "weightfoucsed", "interventiontype"),
       leftlabs = c("Study","N Int","N Ctrl", "Risk factors", "Intervention Aims", "Intervention Type"), 
       rightcols = c("effect", "ci", "w.random"),      
       rightlabs = c("MD [RR]", "95% CI [RR]", "Weight"),
       label.right="favours control",
       label.left = "favours intervention",
       colgap.forest = "1 cm",   # more room for plot
       colgap.left = "0.7cm",   # more gap between left cols
       fontsize = 12)
dev.off()

