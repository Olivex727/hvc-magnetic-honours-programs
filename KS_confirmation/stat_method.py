import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
 
# Load the CSV file
results_statistician_copy = pd.read_csv("C:/Users/olive/OneDrive - Australian National University/Honours-Olivia/Resources/CSV/results_statistician_copy.csv")

# Rename the columns
results_statistician_copy.columns = ["Name", "KS_EDF", "KS_EDF_Unc", "Var_Sub", "Var_Sub_Unc", "Wgt_Mean", "Wgt_Mean_Unc"]
 
# Melt the DataFrame for estimates
Est = pd.melt(results_statistician_copy[["Name", "KS_EDF", "Var_Sub", "Wgt_Mean"]],
    id_vars=["Name"],
    value_vars=["KS_EDF", "Var_Sub", "Wgt_Mean"],
    var_name="Variable",
    value_name="Estimate")
 
# Melt the DataFrame for uncertainty
Uncertainty = pd.melt(results_statistician_copy[["Name", "KS_EDF_Unc", "Var_Sub_Unc", "Wgt_Mean_Unc"]],
    id_vars=["Name"],
    value_vars=["KS_EDF_Unc", "Var_Sub_Unc", "Wgt_Mean_Unc"],
    var_name="Variable_Unc",
    value_name="Weight")

# Merge the two DataFrames
data_new = pd.merge(Est, Uncertainty, on="Name")

# Perform the ANOVA
model = ols('Estimate ~ C(Variable)', data=data_new, weights=data_new['Weight']).fit()
anova_table = model.summary()
print(anova_table)

# Tukey's HSD test
tukey_test = pairwise_tukeyhsd(endog=data_new['Estimate'], groups=data_new['Variable'], alpha=0.05)
print(tukey_test)