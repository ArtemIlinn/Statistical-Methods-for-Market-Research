import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# Original data
y_T_0 = [10, 10, 9, 11]
y_T_1 = [8, 8, 9, 7]
y_C_0 = [11, 12, 13, 13, 11]
y_C_1 = [7, 7, 7, 8, 9]

# Manual DiD calculation
did = (np.mean(y_T_1) - np.mean(y_T_0)) - (np.mean(y_C_1) - np.mean(y_C_0))
print(f"Manual DiD estimate: {did:.2f}")

# Create dataframe
y = y_T_0 + y_T_1 + y_C_0 + y_C_1
x = [1]*len(y_T_0 + y_T_1) + [0]*len(y_C_0 + y_C_1)
z = [0]*len(y_T_0) + [1]*len(y_T_1) + [0]*len(y_C_0) + [1]*len(y_C_1)

df = pd.DataFrame({'y': y, 'treatment': x, 'post': z})

# Linear regression model
model = smf.ols('y ~ treatment + post + treatment:post', data=df).fit()
print("\nRegression results:")
print(model.summary())

# Bonus: Fixed effects estimation (using panel data format)
# First create entity (group) and time indicators
df['entity'] = ['Treatment']*len(y_T_0 + y_T_1) + ['Control']*len(y_C_0 + y_C_1)
df['time'] = ['Pre']*len(y_T_0) + ['Post']*len(y_T_1) + ['Pre']*len(y_C_0) + ['Post']*len(y_C_1)

# Using linearmodels (needs installation: pip install linearmodels)
from linearmodels import PanelOLS

df = df.set_index(['entity', 'time'])
fe_model = PanelOLS.from_formula('y ~ treatment + post + treatment:post', data=df).fit()
print("\nFixed effects results:")
print(fe_model.summary)