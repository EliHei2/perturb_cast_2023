#  2021-06-30 15:17 
#  elihei  [<eheidari@student.ethz.ch>]
# /Volumes/Projects/meso2021/code/R/00_prep_mouse_data.R

# Setup
set.seed(1996)
## libraries
suppressMessages({
    # library('pika')
    library('tidyverse')
    library('data.table')
    library('SingleCellExperiment')
    library('magrittr')
    library('zellkonverter')
    # library('MouseGastrulationData')
})

## helper scripts
source('code/R/utils.R')

## Input paths
data_dir   = 'data_raw/seqfish_mouse_embryo'
exprs_f    = file.path(data_dir, 'exprs.Rds')
metadata_f = file.path(data_dir, 'metadata.Rds')
# nngraph_f  = file.path(data_dir, 'neighbourGraph_1.3.Rds')

atlas_f    = file.path(data_dir, 'MGD_atlas/atlas_SCE.Rds')
type_map_f = file.path(data_dir, 'MGD_atlas/atlas_seqFISH_celltype_map.Rds')


exprs    = exprs_f %>% readRDS
metadata = metadata_f %>% readRDS 






%>% as.data.table %>%
    setnames(
        c('uniqueID', 'celltype_mapped_refined', 'x_global_affine', 'y_global_affine'), 
        c('cell_id', 'cell_type', 'x', 'y')
    ) %>%
    setkey('cell_id') %>%
    .[cells] %>%
    .[, embryo := paste(embryo, z, sep='_')] %>%
    .[cell_type != 'Low quality'] %>%
    .[, .(cell_id, embryo, x, y, UMAP1, UMAP2, cell_type)]
embryos       = unique(metadata$embryo)
cells_list   = embryos %>%
    map(~unlist(metadata[embryo == .x, 'cell_id'])) %>%
    setNames(embryos)
nngraph_list = cells_list %>% 
    map(~induced.subgraph(nngraph, .x))
gc()

atlas = atlas_f %>% readRDS
rownames(atlas) = rowData(atlas)$SYMBOL
type_map = type_map_f %>% readRDS %>% as.data.table %>%
    setnames(c('celltype', 'cell_type')) %>%
    setkey('celltype') %>%
    .[is.na(cell_type), cell_type := celltype] %>%
    .[celltype == 'ExE endoderm', cell_type := 'ExE endoderm']%>%
    .[celltype == 'Paraxial mesoderm', cell_type := 'Paraxial mesoderm'] %>%
    unique
col_dt = atlas %>% colData %>% as.data.table %>%
    .[!is.na(cluster)] %>%
    setkey('celltype') %>%
    .[type_map] %>%
    .[!is.na(celltype)] %>%
    setnames('cell', 'cell_id') %>%
    .[, class_ := cell_type] %>%
    DataFrame 
atlas %<>% .[, col_dt$cell_id]
colData(atlas) = col_dt