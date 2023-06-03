# 取基因表达
```R
scRNA = qread("...")
write.csv(as.matrix(scRNA[gene$Gene]@assays$RNA@data), "./gene_df.csv")
```

# sg_score
```python
meta = pd.read_csv("./meta.csv", index_col=0)
gene_df = pd.read_csv('./gene_df.csv,' index_col=0).T
s_celltype = pd.DataFrame(index = ['score'], columns = set(meta.sub_cluster)).fillna(0)

for i in s_celltype.columns:
    major = meta.loc[meta.sub_cluster == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Mye/batch/cca/crc_mye_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'TNK'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/TNK/batch/cca/crc_tnk_score/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Fib'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Fib/crc_fib_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'Epi'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/Epi/crc_epi_score/"+ i + "_pct.csv", index_col = 0)
        s_celltype.loc['score', i] = res.mean().values[0]
    if (major == 'End'):
        res = pd.read_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/End/percent/"+ i + "_pct.csv", index_col=0)
        s_celltype.loc['score', i] = res.mean().values[0]

s_celltype.to_csv("/public/home/yuwenqi/sc-data/selected/CRC/workspace/module/cluster_sg_scores/s_celltype.csv")

module = pd.read_csv("./module.csv", encoding = 'gbk', index_col=0)
percent = pd.read_csv("./percent.csv", index_col=0)

data = pd.DataFrame(index = percent.index, columns = percent.columns).fillna(0)
for i in percent.index:
    for j in percent.columns:
        mod_col = j.replace(' ', '_')
        data.loc[i, mod_col] = percent.loc[i, j] * s_celltype.loc['score', j]

res = pd.DataFrame(index = data.index, columns = module.columns)
for p in res.index:
    for i in module.columns:
        s_subtme = 0
        for j in module.loc[:, i].dropna():
            s_subtme += (percent.loc[p, j].mean() * s_celltype.loc['score', j]) ** 0.5
        res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())
```

# gene_expr
```python
s_celltype = pd.DataFrame(index = set(meta['orig.ident']), columns = set(meta.sub_cluster)).fillna(0)
gene_df = gene_df.loc[meta.index,]

for i in s_celltype.columns:
    major = meta.loc[meta.sub_cluster == i, 'cluster'].values[0]
    if (major == 'Mye'):
        res = pd.read_csv("./Mye/bc_mye_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'TNK'):
        res = pd.read_csv("./TNK/bc_tnk_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Fib'):
        res = pd.read_csv("./Fib/bc_fib_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'Epi'):
        res = pd.read_csv("./Epi/bc_epi_score/"+ i + "_pct.csv", index_col = 0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
    if (major == 'End'):
        res = pd.read_csv("./End/bc_end_score/"+ i + "_pct.csv", index_col=0)
        if(len(res) == 0):
            continue
        cur = gene_df.loc[meta.sub_cluster == i, res.index]
        cur['sample'] = meta.loc[cur.index, 'orig.ident']
        s_celltype.loc[:, i] = cur.groupby('sample').mean().mean(axis=1)
s_celltype = s_celltype.fillna(0)
s_celltype.to_csv("./BC3-CCI/module/cluster_sg_gene_scores/s_celltype_gene.csv")

data = pd.DataFrame(index = percent.index, columns = percent.columns).fillna(0)
for i in data.index:
    for j in data.columns:
        data.loc[i, j] = percent.loc[i, j] * s_celltype.loc[i, j]

res = pd.DataFrame(index = data.index, columns = module.columns)
for p in res.index:
    for i in module.columns:
        s_subtme = 0
        for j in module.loc[:, i].dropna():
            s_subtme += (percent.loc[p, j].mean() * s_celltype.loc[p, j]) ** 0.5
        res.loc[p, i] = s_subtme / len(module.loc[:, i].dropna())
```