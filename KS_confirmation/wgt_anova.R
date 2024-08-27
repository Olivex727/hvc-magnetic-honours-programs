results_statistician_copy <- read.csv("./results_statistician_copy.csv",header = T)

names(results_statistician_copy)<- c("Name","KS_EDF","KS_EDF_Unc","Var_Sub","Var_Sub_Unc","Wgt_Mean","Wgt_Mean_Unc")

library(reshape2)

Est <- melt(results_statistician_copy[,c("Name","KS_EDF","Var_Sub","Wgt_Mean")],id.vars = "Name", measure.vars = c("KS_EDF","Var_Sub","Wgt_Mean"),value.name = "Estimate")

Uncertainity <- melt(results_statistician_copy[,c("Name","KS_EDF_Unc","Var_Sub_Unc","Wgt_Mean_Unc")],id.vars = "Name",measure.vars = c("KS_EDF_Unc","Var_Sub_Unc","Wgt_Mean_Unc"),value.name = "Weight")

data_new <- merge(Est,Uncertainity,by="Name")

base_V2_aov <- aov(data_new$Estimate ~ data_new$variable.x, weights =data_new$Weight)

summary(base_V2_aov)

tukey.test <- TukeyHSD(base_V2_aov)

tukey.test
