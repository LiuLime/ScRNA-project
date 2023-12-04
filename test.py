import os

import pandas as pd


path=os.getcwd()
with open (f'{path}/Lung.Cholangiocyte_merge.csv','r') as f:
    data = pd.read_csv(f)

print(data.head)


label=data['label']
count=0
for l in label:
    if l == 'Young':
        count+=1

print(count)


#%%
a = [1,2,3,4,6,8,9]
batch=(len(a)//3)+1

#%%
