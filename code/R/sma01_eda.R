#  2021-03-05 11:07 
#  elihei  [<eheidari@student.ethz.ch>]
# /Volumes/Projects/MS_lesions/code/ms10_glasso.R


library('biomaRt')
library('mltools')
# library('pika')
library('leiden')
library('ComplexHeatmap')
library('seriation')
library('viridis')
library('magrittr')


nngraph_comm <- function(nngraph, min_cells = 100, res_param=0.1){
    comm_vec = leiden(nngraph, resolution_parameter=res_param)
    comm_dt  = table(comm_vec)
    to.keep  = max(as.numeric(names(comm_dt)[comm_dt > min_cells]))
    comm_vec = ifelse(comm_vec <= to.keep, comm_vec, NA)
    names(comm_vec) = names(V(nngraph))
    comm_vec
}


counts2graph <- function(mtx, rho = 0.1, threshold = 0.01, n_hvgs=500, keep_isolated=FALSE, no_hvg=FALSE, corr=FALSE, cut_off=0.1){
    if(!no_hvg){
        gene_vars = modelGeneVar(as.matrix(mtx))
        hvg       = getTopHVGs(gene_vars) %>% .[1:n_hvgs]
        data_norm = t(mtx[hvg,]) %>% as.data.table
    }else
        data_norm = mtx %>% as.data.table
    colnames(data_norm) = strsplit2(colnames(data_norm), split='_')[,1]
    ## penalization parameter
    ggm   = ggm(data_norm, rho=rho, threshold=threshold)
    wi    = ggm$model$wi
    ## exclude isolated nodes
    gg = ggm$graph
    if(!keep_isolated){
        isolated  = which(igraph::degree(gg) == 0)
        gg = delete.vertices(gg, isolated)
        ## gg weights
        wi %<>% .[-isolated, -isolated]
    }

    colnames(wi) = rownames(wi) = names(V(gg))
    diag(wi) = 0
    ## visualize gg
    comm  = leiden(abs(wi) > threshold)
    names(comm) = colnames(wi)
    markers = unique(comm) %>% 
        map(~names(sort(colSums((abs(wi)>threshold)[names(comm)[comm==.x], names(comm)[comm==.x]]), decreasing=T))[1]) %>%
        unlist
    comm_dt = data.table(GENE=names(comm), community=comm, marker=markers[comm]) %>%
        setkey(community) %>% 
        .[, color := c(nice_cols_1, nice_cols_2)[community]] %>%
        setkey(GENE) %>%
        .[names(V(gg))]
    if(!keep_isolated)
        data_norm  %<>% as.matrix %>% .[, -isolated] 
    trans = one_hot(comm_dt[, .(ID = GENE, factor(color))]) %>%
        setkey(ID) %>%
        .[colnames(data_norm)] %>%
        .[, -'ID'] %>%
        as.matrix 
    rownames(trans) = comm_dt$GENE
    data_trans = data_norm %*% trans
    gg = graph_from_adjacency_matrix(abs(wi) > threshold, 'undirected')
    if(corr){
        S         = stats::cov.wt(as.data.table(data_norm), method='ML')$cov
        wi        = stats::cov2cor(S)
        diag(wi)  = 0
        gg = graph_from_adjacency_matrix(abs(wi) > cut_off, 'undirected')
    }

    
    lay = layout_nicely(gg)
    #TODO: turn to a class
    list(
        graph      = gg, 
        data       = data_norm, 
        data_trans = data_trans, 
        comm_dt    = comm_dt, 
        wi         = wi, 
        trans      = trans,
        lay        = lay
    )
}

# corr_graph <- function(data, threshold=0.5){
#     model     = gRim::cmod(~ . ^ ., data=data)
#     S         = stats::cov.wt(data, method='ML')$cov
#     C         = stats::cov2cor(S)
#     AM        = abs(C) > threshold
#     diag(AM)  = F
#     g.corr   = as(AM, 'graphNEL')
#     graph::nodes(g.corr) = colnames(data)
#     g.corr
    
# }

act_colors <- function(data, var_vec){
    grp    = list()
    for(l in unique(var_vec)){
        idx = which(var_vec == l)
        dn  = data[idx, ]
        if(length(idx) == 1)
            grp[[l]] = dn
        else
            grp[[l]] = colMeans(dn)
    }
    if(length(grp) > 1){
        grp = grp %>% do.call(rbind, .)
        grp = apply(grp, 2, function(x) (x)/(max(x)))
    }else{
        grp = unlist(grp) %>% as.matrix %>% t
        rownames(grp) = unique(var_vec)
    }
    grp[is.infinite(grp)] = 0
    grp
}

graph_col_comm <- function(graph, lay, grp, sz, title=NULL, labels){
    igraph::V(graph)$color <- grp
    v <-  igraph::V(graph)
    # sprintf(comm_out, title) %>% pdf()
    p = plot.igraph(
        graph,
        vertex.size = 6,
        layout = lay,
        vertex.label = labels,
        vertex.frame.color = igraph::V(graph)$color,
        vertex.label.family = 'Helvetica',
        vertex.label.dist = 0,
        vertex.label.cex = .5,
        vertex.label.font = 2,
        vertex.label.color = '#585c59',
        main=title)
    # dev.off()
    p
}

graph_col_act <- function(graph, grp, lay, title){
    # sprintf(graph_out, paste0(title, '_train')) %>% pdf()
    grp_range = c(min(grp)^(1)/sum(grp^(1)),max(grp)^(1)/sum(grp^(1)))
    grp_vals  = seq(grp_range[1],grp_range[2],length.out=9)
    grp_cols  = circlize::colorRamp2(grp_vals, viridis::viridis(9))
    igraph::V(graph)$color = grp_cols(grp^(1)/sum(grp^(1)))
    p = plot.igraph(graph,
        vertex.size = 5,
        layout = lay,
        vertex.frame.color = igraph::V(graph)$color,
        vertex.label = "",
        main=title)
    p
}

plot_2d <- function(dim_df, labels, label_cols=nice_cols, title='', label_title='label', hide_legend=FALSE){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        geom_point(size=0.5) +
        theme_bw() + 
        theme(axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
        scale_color_manual(values=label_cols, na.value='gray') +
        labs(title=title, x='', y='', color=label_title)
    if(hide_legend)
        dim_plot = dim_plot + theme(legend.position='none')
    dim_plot
}

plot_2d_cont <- function(dim_df, labels, label_cols=nice_cols, title='', label_title='label'){
    dim_dt = data.table(label=labels,
                         dim1=unlist(dim_df[,1]), dim2=unlist(dim_df[,2])) 
    dim_plot = dim_dt %>%
        ggplot +
        aes(dim1, dim2, color=label) +
        # geom_hex(bins = 30) + 
        geom_point(size=0.5) +
        # coord_fixed() +
        scale_color_viridis() +
        theme_bw() + 
        theme(legend.position = "none",
            axis.text= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
        labs(title=title, x='', y='', color=label_title)
    dim_plot
}

hm_col_act <- function(
    mtx, 
    comm_dt,
    col_dt,
    cluster_rows=FALSE,
    cluster_columns=TRUE) {
    mtx = t(mtx)
    stopifnot(dim(mtx)[1] == dim(comm_dt)[1])
    stopifnot(dim(mtx)[2] == dim(col_dt)[1])
    col_dt = col_dt[colnames(mtx)]
    col_cols = col_dt$color
    names(col_cols) = col_dt$group
    # make column annotations
    col_annots  = HeatmapAnnotation(
        group                = col_dt$group, 
        col                  = list(group  = col_cols),
        annotation_name_side = 'left', 
        show_legend          = c(group=FALSE)
    )

    row_cols = comm_dt$color
    names(row_cols) = comm_dt$community
    row_annots  = rowAnnotation(
        gene_module          = comm_dt$community, 
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
        cluster_rows=cluster_rows, cluster_columns=cluster_columns,
        row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize = 10),
        name="Scaled\nGene Expression", 
        row_names_side="left",
        column_names_side="top",
        row_split=comm_dt$community,
        cluster_column_slices=FALSE,
        top_annotation=col_annots, left_annotation=row_annots
        )

    return(hm_obj)
}
