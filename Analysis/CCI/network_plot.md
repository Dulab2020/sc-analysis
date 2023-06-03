# get the cnt_data
```python
idx = set(d.columns.str.split('|').str[0]) | set(d.columns.str.split('|').str[1])
res = pd.DataFrame(index = idx, columns = idx)
res = res.fillna(0)
for i in d.columns:
    c1 = i.split('|')[0]
    c2 = i.split('|')[1]
    for v in d.loc[:, i]:
        if (v > 0):
            res.loc[c1, c2] += 1  
res.to_csv("./network_counts.csv")
```

```R
counts = read.csv("./module4_cnt.csv", row.names=1)
count_inter <- as.matrix(counts)
pdf("./module4_plot.pdf", width = 7, height = 7)
netVisual_circle(count_inter,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()
```

## plot for each cell_type
```R
pdf("./each_cluster.pdf", width = 30, height = 80)
par(mfrow = c(8, 3), xpd=TRUE)
for (i in 1:nrow(count_inter)) {
  mat2 <- matrix(0, nrow = nrow(count_inter), ncol = ncol(count_inter), dimnames = dimnames(count_inter))
  mat2[i, ] <- count_inter[i, ]
  netVisual_circle(mat2, 
                   weight.scale = T, 
                   edge.weight.max = max(count_inter), 
                   title.name = rownames(count_inter)[i],
                   arrow.size=0.2)
}
dev.off()
```

