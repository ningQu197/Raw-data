import numpy as np
import pandas as pd

df = pd.read_csv("./TF_AUC_score.txt", header=0, index_col=0, sep="\t")
### 去掉含有_extended的转录因子列
TF_columns = []
for item in df.columns:
	if "_extended" not in item:
		TF_columns.append(item)
df = df.loc[:, TF_columns]
### 将宽数据转变成长数据
total_arr = []
for i in range(df.shape[1]):
	sub_df = df.iloc[:, [i]]
	sub_df["TF"] = sub_df.columns[0]
	sub_df.columns = ["Score", "TF"]
	total_arr.append(sub_df)
total_df = pd.concat(total_arr, axis=0)
### 保留大于0的细胞
total_df = total_df.loc[total_df.Score>0, :]
total_df.to_csv("./TF_score_boxplot.txt", header=True, index=True, sep="\t")





