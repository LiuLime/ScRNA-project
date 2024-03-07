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
import pandas as pd

# 示例数据
data = {
    '部门': ['A', 'A', 'B', 'B', 'C', 'C'],
    '产品1销售额': [100, 150, 200, 250, 300, 350],
    '产品2销售额': [100, 150, 200, 250, 300, 350]
}
df = pd.DataFrame(data)

# 自定义函数1：计算总和
def custom_sum(x):
    return x.sum()

# 自定义函数2：计算平均值
def custom_mean(x):
    return x.mean()

# 使用 groupby + agg
result = df.groupby('部门').agg({'产品1销售额': custom_sum, '产品2销售额': custom_mean}).reset_index()

print(result)

#%%
import pandas as pd

# 示例数据
data = {
    '部门': ['A', 'A', 'B', 'B', 'C', 'C'],
    '产品1销售额': [100, 150, 200, 250, 300, 350],
    '产品2销售额': [100, 150, 200, 250, 300, 350]
}
df = pd.DataFrame(data)

# 自定义函数，对两列求和
def custom_sum(x):
    return x['产品1销售额'] + x['产品2销售额']

# 使用 groupby + apply
result = df.groupby('部门').apply(lambda x: custom_sum(x)).reset_index(name='销售总额')
result2=df.apply(lambda x: custom_sum(x))
print(result2)
