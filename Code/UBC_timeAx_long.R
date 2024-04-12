#load the required packaes:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/initial_settings.R")

#source the tailored functions:
source("/media/chronos/Storage/ayelet/ProjectsCode/CancerTrajectory/Code/my_functions.R")

# read GEO data:
gse_UBC_long <- getGEO("GSE128959")[[1]]
head(pData(gse_UBC_long))
gene_info <- fData(gse_UBC_long)
featureData(gse_UBC_long)

expression_data <- read.table(file = "/media/chronos/Storage/ayelet/UBC_data_GSE128959/GSE128959_Processed_data_supplement.txt", row.names = 1, header = T)
phen_data_short <- pData(gse_UBC_long)[,c("title", "patient identifier:ch1", "primary, recurrence, or progression:ch1", "tumor number:ch1", "tumor stage:ch1", "molecular subtype:ch1")]
colnames(phen_data_short) <- sub(":ch1", "", colnames(phen_data_short))
colnames(phen_data_short) <- gsub(" ", "_", colnames(phen_data_short))

#get the patients with the long time course (more than 3 time points)
long_time_course_patients <- names(table(phen_data_short$patient_identifier))[table(phen_data_short$patient_identifier) > 3]

#subset the expression and phenotypic data to include only the long time-course patients' data:
expression_data_long_time_course <- expression_data[,phen_data_short[phen_data_short$patient_identifier %in% long_time_course_patients, "title"]]
phen_data_short_long_time_course <- phen_data_short[phen_data_short$patient_identifier %in% long_time_course_patients,]
phen_data_short_long_time_course$tumor_stage <- factor(phen_data_short_long_time_course$tumor_stage, levels = c("Ta", "T1", "T2", "T3", "T4"))
  
#run TimeAx
model = modelCreation(expression_data_long_time_course, phen_data_short_long_time_course$patient_identifier)

#get pseudotime values:
pseudotimeStats = predictByConsensus(model, expression_data_long_time_course)
pseudotime = pseudotimeStats$predictions
uncertainty = pseudotimeStats$uncertainty

#check the correlation with tumor stage:
phen_data_short_long_time_course$pseudotime <- pseudotime
ggplot(phen_data_short_long_time_course, aes(x = tumor_stage, y = pseudotime)) + geom_boxplot() + geom_jitter(width = 0.1)