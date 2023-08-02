methods     = list.files(out_dir)
preds_f   = list.files('.')[grep('pred', list.files('.'))] 

for(i in 1){
	preds = preds_f[i] %>% readH5AD  %>% assays %>% .[['X']] %>% t
	if(i == 1){
		dists = JSD(preds) 
		dists = dists/norm(dists)
	}else{
		d  = JSD(preds) 
		d  = d/norm(d)
		dists = dists + d
	}
}

rm(d)
gc()

preds = preds_f[3] %>% readH5AD  %>% assays %>% .[['X']] %>% t
col_data = preds_f[3] %>% readH5AD %>% colData
pred_labs = apply(preds, 1, function(x) ifelse(max(x) <0.3, NA, which.max(x)))
# preds = apply(preds, 1, function(x) ifelse(x == max(x), 1, 0)) %>% t
dists = dist(preds[1:10000,]) 

================
dists1_f = 'dists_atlas_8.5.h5ad'
dists1_f %>% readH5AD -> dists1
ind = which( ((rank(-dists1$ent_01) < 5000)	+ 
	(rank(-dists1$ent_05) < 5000) +
	(rank(-dists1$ent_5) < 5000) +
	(rank(-dists1$ent_1) < 5000)	+ 
	(rank(-dists1$ent_10) < 5000)) >3)
dists1$alpha = 1
dists1$alpha[ind] = 0.5
dists2_f = 'dists_embryo3_2_atlas_8.5.h5ad'
dists2_f %>% readH5AD -> dists2
ind = which( ((rank(-dists2$ent_01) < 5000)  + 
    (rank(-dists2$ent_05) < 5000) +
    (rank(-dists2$ent_5) < 5000) +
    (rank(-dists2$ent_1) < 5000) + 
    (rank(-dists2$ent_10) < 5000)) >3)
dists2$alpha = 1
dists2$alpha[ind] = 0.5

dist_m1 = assays(dists1)$X / norm(assays(dists1)$X)
dist_m2 = assays(dists2)$X / norm(assays(dists2)$X)
map2d = Rtsne(dist_m2, is_distance=T, perplexity=60)$Y

exprs_umap_p1 = plot_2d( 
    dim_df = map2d,
    labels = factor(as.character(dists1$celltype), ordered=TRUE),
    label_cols = celltype_colours,
    hide_legend = TRUE,
    title = 'cell type',
    alpha = dists1$alpha * dists2$alpha,
    sz = 1
) + labs(x = '', y = 'aggregated space')

dists$ent_01 = dists$ent_01/sqrt(sum(dists$ent_01^2))
dists$ent_05 = dists$ent_05/sqrt(sum(dists$ent_05^2))
dists$ent_1 = dists$ent_1/sqrt(sum(dists$ent_1^2))
dists$ent_5 = dists$ent_5/sqrt(sum(dists$ent_5^2))
dists$ent_10 = dists$ent_10/sqrt(sum(dists$ent_10^2))
dists$ent = dists$ent_01 + dists$ent_05 + dists$ent_1 + dists$ent_5 + dists$ent_10



exprs_umap_p2 = plot_2d( 
    dim_df = map2d,
    labels = factor(as.character(pred_labs[1:10000], ordered=TRUE)),
    label_cols = .palette1,
    hide_legend = TRUE,
    title = 'predicted region, res = 0.1',
    sz = 1
) + labs(x = '', y = 'aggregated space')

exprs_umap_p1 + exprs_umap_p2

exprs_umap_p1 = plot_2d( 
    dim_df = map2d,
    labels = factor(as.character(dists$pred_10), ordered=TRUE),
    label_cols = c(.palette1, .palette2, .palette3),
    hide_legend = TRUE,
    title = 'predicted region, res = 1',
    sz = 1
) + labs(x = '', y = 'aggregated space')

exprs_umap_p1 = plot_2d_cont( 
    dim_df = map2d,
    labels = dists$ent,
    title = 'entropy',
    sz = 1
) + labs(x = '', y = '')




rownames(map2d) = colData(dists)$cell_id
g_out    = map2d %>% 
      get_delaunay(plot=FALSE) %>%
      .$graph
colnames(dists) = colData(dists)$cell_id
ccc = cellCellContact(
  sce    = dists,
  group  = 'cell_type',
  graph  = g_out,
  nperm  = 500,
  plot   = FALSE,
  cellID ='cell_id'
)

colData(dists) %>%
  as.data.table %>% 
  ggplot +
  aes(x = factor(cell_type, ordered=TRUE, levels=names(celltype_colours)), y = ent, fill=factor(cell_type, ordered=TRUE, levels=names(celltype_colours))) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=celltype_colours, na.value='gray') +
  labs(y='entropy', fill = 'cell type', x= 'cell type') +
  theme(legend.position='none') +
  scale_color_manual(values=celltype_colours, na.value='gray')+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))