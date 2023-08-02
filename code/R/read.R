library(SpatialExperiment)
SE_ML_I     = read10xVisium('data_raw/ML_II_A_1')


plasmid_f   = 'data_raw/decodes/Plasmid_DECODEv1.csv'
olftasvmn_f = 'data_raw/decodes/OLFTASVMN_DECODEv1.csv'
bc_f        = 'data_raw/decodes/BC_DECODEv1.csv'
plasmid_dt = fread(plasmid_f)
olftasvmn_dt = fread(olftasvmn_f)
bc_dt = fread(bc_f)
SE_ML_I = scuttle::logNormCounts(SE_ML_I) 
gene_vars = modelGeneVar(SE_ML_I)
hvg       = getTopHVGs(gene_vars) %>% .[1:500]
genes = c(plasmid_dt$ID, olftasvmn_dt$V3, bc_dt$ID, hvg)
genes = unique(genes)
genes = genes[genes!='']
exprs = counts(SE_ML_I)
genes = intersect(genes, rownames(exprs))

exprs = as.matrix(exprs[genes,]) 
exprs = exprs[rowSums(exprs) != 0,]
# exprs = log(exprs+1)
rownames(exprs) = rowData(SE_ML_I[rownames(exprs),])$symbol
# graph_all = counts2graph(
#     mtx       = exprs,
#     rho       = 0.2,
#     threshold = 0.02,
#     no_hvg    = TRUE
# )

net = get_network(t(exprs), rho = .1, threshold = .1)
net = as.undirected(net)

Isolated = which(igraph::degree(net)==0)
net = delete.vertices(net, Isolated)
adj = as_adj(net)
genes1 = colnames(adj)
# gg = analyze_network(net)


lc = cluster_louvain(net)
m = membership(lc)
cc = communities(lc)


graph_col_comm(
    graph  = net, 
    lay    = layout_nicely(net), 
    grp    = m, 
    sz     = 4, 
    title  = 'GGM on ML_I', 
    labels = names(m)
)

layout_nicely(graphs_embryo[[t]]$graph)



library(reticulate)

np <- import("numpy")
a <- np$array(adjacency_matrix)

np$save("data_tidy/test_ajd" , a)




act_all_nodules  = act_colors(t(exprs_list[[1]]), samples[[1]]$nodule)
act_all_nodules = act_all_nodules[,names(memberships[[1]])]


mtx = t(act_all_nodules)
col_cols = .palette_all[rank(rownames(act_all_nodules))]
names(col_cols) = rownames(act_all_nodules)
col_annots  = HeatmapAnnotation(
    group                = rownames(act_all_nodules), 
    col                  = list(group  = col_cols),
    annotation_name_side = 'left', 
    show_legend          = c(group=FALSE)
)

row_cols = .palette1[memberships[[1]]]
names(row_cols) = memberships[[1]]
# names(row_cols) = comm_dt$community
row_annots  = rowAnnotation(
    gene_module          = memberships[[1]], 
    col                  = list(gene_module  = row_cols),
    show_legend          = c(gene_module=FALSE)
)

grp_vals  = seq(min(mtx), max(mtx), length.out=9)
grp_cols  = circlize::colorRamp2(grp_vals, viridis::viridis(9))
if(dim(as.matrix(mtx))[2] > 1){
    seriate_obj  = seriate(as.matrix(mtx - min(mtx)), method = "BEA_TSP")
    row_order    = get_order(seriate_obj, 1)
    column_order = get_order(seriate_obj, 2)
}else{
    row_order    = 1:length(mtx)
    column_order = 1
}

# do heatmap for these genes
hm_obj      = Heatmap(
    matrix=mtx, col=grp_cols,
    row_order=row_order, 
    cluster_rows=TRUE, cluster_columns=TRUE,
    row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize = 10),
    name="Scaled\nGene Expression", 
    row_names_side="left",
    column_names_side="top",
    row_split=memberships[[1]],
    cluster_column_slices=FALSE,
    top_annotation=col_annots, left_annotation=row_annots
    )





plotSpots(samples[[1]], annotate = "nodule", palette=col_cols)




