library(STdeconvolve)

pos <- spatialCoords(samples[[2]])
colnames(pos) <- c('x', 'y')
colnames(samples[[2]]) <- colData(samples[[2]])$Barcode
cd <- as(counts(samples[[2]]), 'Matrix')
annot <- colData(samples[[2]])$leiden_sagenet
## remove pixels with too few genes
cd <- cleanCounts(cd, min.lib.size = 100)
## feature select for genes
corpus <- restrictCorpus(cd, removeAbove=1.0, removeBelow = 0.05)
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(9, 9, by = 1))
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
## visualize deconvolved cell-type proportions
pl <- vizAllTopics(deconProp, pos[rownames(deconProp),],
             groups = annot[rownames(deconProp)], 
             group_cols = rainbow(length(levels(annot[rownames(deconProp)]))),
             r=100) 

for(t in 1:dim(deconProp)[2]){
    pl <- vizTopic(theta = deconProp, pos = pos[rownames(deconProp),], topic = t, plotTitle = paste('celltype', t),
         size = 5, stroke = 1, alpha = 0.5,
         low = "white",
         high = "red")
    print(pl)

}




colnames(deconGexp) = rowData(samples[[2]])[colnames(deconGexp),]

deconGexp <- deconGexp[, intersect(colnames(deconGexp), markers)]


mtx = t(deconGexp) %>% apply(., 1, function(x) (x)/(max(x)))  %>% t

mtx=cor_mat
col_cols = .palette2[colnames(cor_mat)]
names(col_cols) = colnames(cor_mat)
col_annots  = HeatmapAnnotation(
    group                = colnames(cor_mat), 
    col                  = list(celltype  = col_cols),
    annotation_name_side = 'right', 
    show_legend          = c(group=FALSE)
)

is_genes = setdiff(rownames(cor_mat), names(cons_mem))
sel_mem = cons_mem[rownames(cor_mat)] %>% na.exclude %>% c
is_mem  = rep(max(sel_mem, na.rm = T) + 1, length(is_genes))
names(is_mem) = is_genes
sel_mem = c(sel_mem, is_mem)
row_cols = .palette1[sel_mem]
names(row_cols) = sel_mem
# names(row_cols) = comm_dt$community
row_annots  = rowAnnotation(
    gene_module          = sel_mem, 
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
    matrix=cor_mat[names(sel_mem),], col=grp_cols,
    # row_order=row_order, 
    cluster_rows=TRUE, cluster_columns=TRUE,
    row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize = 10),
    name="correlation", 
    row_names_side="left",
    column_names_side="top",
    row_split=sel_mem,
    cluster_column_slices=FALSE,
    top_annotation=col_annots, left_annotation=row_annots
    )
draw(hm_obj)


corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp,
                                        type = "t")




sample_test = samples_raw[[2]]
sample_markers = sample_test[marker_ens,]
rownames(sample_markers) = markers
cor_mat = cor(t(as.matrix(as(logcounts(samples[[2]][,rownames(deconProp)]), 'Matrix'))), (deconProp))
rownames(cor_mat) <- rowData(samples[[2]])[rownames(cor_mat),]

hm_obj      = Heatmap(
    matrix=cor_mat,
    # , col=grp_cols,
    # row_order=row_order, 
    cluster_rows=TRUE, cluster_columns=TRUE,
    row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize = 10),
    name="Scaled\nGene Expression", 
    row_names_side="left",
    column_names_side="top",
    # row_split=sel_mem,
    cluster_column_slices=FALSE,
    # top_annotation=col_annots, left_annotation=row_annots
    )
draw(hm_obj)




for(m in markers){
    pl = vis_gene(
    spe = sample_markers,
    # sampleid = "151673",
    geneid = m,
    point_size = 3,
    alpha=0.5
)
    print(pl)

}












mobProxyTheta <- model.matrix(~ 0 + annot)
rownames(mobProxyTheta) <- names(annot)
# fix names
colnames(mobProxyTheta) <- unlist(lapply(colnames(mobProxyTheta), function(x) {
  unlist(strsplit(x, "annot"))[2]
}))
rownames(mobProxyTheta) <- colnames(samples[[2]])

mobProxyGexp <- cd[rownames(corpus),] %*% mobProxyTheta[colnames(cd),]
rownames(mobProxyGexp) <- rowData(samples[[2]])[rownames(mobProxyGexp),]


corMtx_beta <- getCorrMtx(# the deconvolved cell-type `beta` (celltypes x genes)
                          m1 = as.matrix(deconGexp),
                          # the reference `beta` (celltypes x genes)
                          m2 = t(as.matrix(mobProxyGexp)),
                          # "b" = comparing beta matrices, "t" for thetas
                          type = "b")



rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))

correlationPlot(mat = corMtx_beta,
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Transcriptional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))













dt_list = list()
for(s in names(samples)){
    tt             = table(samples[[s]]$nodule, samples[[s]]$leiden_sagenet) 
    norm_props     = t(apply(tt,1, function(x) x/sum(x))) 
    nodule_unified = apply(norm_props, 1, function(x) names(x)[which.max(x)])
    purity         = apply(norm_props, 1, function(x) x[which.max(x)])
    dt_list[[s]]   = data.table(sample=s, nodule_pre=rownames(tt), nodule_unified=nodule_unified, purity=purity)
}
dt_comb = dt_list %>% purrr::reduce(rbind)

fwrite(dt_comb, 'output/nodules_unified.txt')








