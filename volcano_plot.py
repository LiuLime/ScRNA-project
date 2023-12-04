import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
filename='Bladder.Fibroblast_FSC10_L1'

df = pd.read_csv(f'./ML/{filename}')
df_orgin = pd.read_csv('./ML/gene-tissue-cellType_mtx.log10P_greaterSide_med_副本.txt', delimiter='\t', header=0,
                       index_col=0)
# %%
tissue = '17.human_TabulaSapiens:normal:Bladder:Fibroblast'
genes = df['gene']
pvalue = pd.DataFrame(df_orgin.loc[genes, tissue])


# 由于ling-san的old是正值，young是负值，而我的young是正值，所以她的原数据需要*-1
def reverse_label(x):
    return x * -1


pvalue = pvalue.applymap(reverse_label)
df=pd.merge(df,pvalue,how='right',on='gene')
# %%
top_pos_genes = df.nlargest(10, 'coef')
top_neg_genes = df.nsmallest(10, 'coef')
plt.figure(figsize=(10, 6))
plt.scatter(df['coef'], df[tissue])

# 设置图形属性
plt.axhline(y=-1, color='black', linestyle='--', linewidth=0.7)
plt.axvline(x=0, color='black', linestyle='--', linewidth=0.7)
plt.xlabel('coef')
plt.ylabel('log10P')
plt.title('Volcano Plot')
plt.legend()
plt.grid(True)

# 显示图形
plt.show()
