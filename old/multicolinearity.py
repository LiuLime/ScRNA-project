from statsmodels.stats.outliers_influence import variance_inflation_factor
import pandas as pd


merge_df = pd.read_csv('merge_df.csv')
#%%

x = merge_df.iloc[:, 1:-1].astype(float)
# 计算VIF因子
vif_data = pd.DataFrame()
vif_data["feature"] = x.columns
vif_data["VIF"] = [variance_inflation_factor(x.values, i) for i in range(len(x.columns))]
colinear=vif_data.query('VIF > 10')