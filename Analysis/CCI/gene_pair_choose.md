data = pd.read_csv("./BC3-CCI/module/module/cellchat/cellchat_module4_v.csv", index_col = 0)
# 1. get CCI score top 100
data_filted = data.loc[data.sum(axis=1).sort_values()[-100:].index, ]
# 2. sorted by cv
data_cv = data_filted.std(axis=1) / data_filted.mean(axis=1)
data_filted = data_filted.loc[data_cv.sort_values(ascending = False).index]
data_filted.to_csv("./BC3-CCI/module/module/cellchat/filted/module4_data_filted.csv")