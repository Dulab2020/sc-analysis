## bubble plot
sel_means=read.csv("./bc/sorted_bymodule_v.csv", row.names=1)
sel_pval=read.csv("./bc/sorted_bymodule_v.csv", row.names=1)
selected_rows = rownames(sel_means)
selected_cols = colnames(sel_means)

df_names = expand.grid(selected_rows, selected_cols)
pval = unlist(sel_pval)
pval[] = 0.0009
#pval[pval == 0] = 0.0009 
plot.data = cbind(df_names, pval)
pr = unlist(as.data.frame(sel_means))
pr = pr + 1
plot.data = cbind(plot.data, log2(pr))
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

my_palette = colorRampPalette(c("white", "blue", 'yellow', 'red'), alpha = TRUE)(n = 399)

ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),     
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggsave("./bc/cellchat.pdf", width = 250, height = 20, device = cairo_pdf, limitsize=F, units = 'in')
