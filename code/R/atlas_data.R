EmbryoAtlasData <- function(type=c("processed", "raw"), samples=NULL, get.spliced=FALSE) {
    type <- match.arg(type)
    versions <- list(base="1.0.0")
    extra_a <- NULL
    if(get.spliced){
        if(type=="raw"){
            stop("Cannot get spliced counts with the raw data")
        }
        extra_a <- list(
            spliced_counts="counts-spliced",
            unspliced_counts="counts-unspliced",
            ambiguous_counts="counts-ambig")
        versions <- c(versions, list(
            "counts-spliced"="1.4.0",
            "counts-unspliced"="1.4.0",
            "counts-ambig"="1.4.0"))
    }
    .getRNAseqData("atlas", type, versions, samples, sample.options=as.character(c(1:10, 12:37)), sample.err="1:10 or 12:37", extra_assays = extra_a)
}

AtlasSampleMetadata <- read.table(
    text = 
        "sample,stage,pool_index,seq_batch,ncells
        1,E6.5,1,1,360
        2,E7.5,2,1,356
        3,E7.5,3,1,458
        4,E7.5,4,1,276
        5,E6.5,5,1,1207
        6,E7.5,6,1,2798
        7,E6.75,7,1,2169
        8,E7.75,8,1,3254
        9,E7.75,8,1,3093
        10,E7.0,9,1,2359
        12,E7.75,11,2,5305
        13,E7.75,11,2,6068
        14,E7.0,12,2,1311
        15,E7.0,12,2,1620
        16,E8.0,13,2,6230
        17,E8.5,14,2,4483
        18,E6.5,15,2,2130
        19,E7.5,16,2,6996
        20,E7.5,16,2,1992
        21,mixed_gastrulation,17,2,4651
        22,mixed_gastrulation,17,2,4674
        23,E7.25,18,2,1429
        24,E8.25,19,2,6707
        25,E8.25,19,2,7289
        26,E7.25,20,2,6649
        27,E7.25,20,2,7216
        28,E8.25,21,2,4646
        29,E8.5,22,3,7569
        30,E7.0,23,3,3785
        31,E7.0,23,3,3778
        32,E7.0,23,3,3718
        33,E8.0,24,3,5443
        34,E8.0,24,3,5314
        35,E8.0,24,3,5072
        36,E8.5,25,3,4915
        37,E8.5,26,3,4011",
    header = TRUE,
    sep = ",",
    stringsAsFactors = FALSE)

.getData <- function(
    dataset,
    version,
    samples,
    sample.options,
    sample.err,
    names,
    object.type=c("SingleCellExperiment"),
    return.list=FALSE
){
    object.type <- match.arg(object.type)
    hub <- ExperimentHub()
    host <- file.path("MouseGastrulationData", dataset)
    #default to all samples
    if (is.null(samples)) {
        samples <- sample.options
    }
    #check for sample boundaries
    samples <- as.character(samples)
    if (!all(samples %in% sample.options) || length(samples)==0) {
        stop(sprintf("'samples' must be in %s", sample.err))
    }

    if(return.list){
        out <- lapply(samples, function(x){ .getData(dataset, version, x,
            sample.options, sample.err, names, object.type, return.list=FALSE)})
        names(out) <- samples
        return(out)
    }

    data = list()

    # Temporary function for data extraction
    EXTRACTOR <- function(target) {
        ver <- .fetch_version(version, target)
        lapply(samples, function(i){
            hub[hub$rdatapath==file.path(host, ver, sprintf("%s-sample%s.rds", target, i))][[1]]
        })
    }
    GET_ASSAYS <- function(ass){
        assay_list <- lapply(seq_along(ass), function(x){
            samp_list = EXTRACTOR(ass[[x]])
            do.call(cbind, samp_list)
        })
        names(assay_list) = names(ass)
        return(assay_list)
    }

    data$assays <- GET_ASSAYS(names$assays)
    
    if(!is.null(names$rd)){
        ver <- .fetch_version(version, "rowdata")
        data$rowData <- hub[hub$rdatapath==file.path(host, ver, paste0(names$rd, ".rds"))][[1]]
    }

    if(!is.null(names$cd)){
        #class change - atlas (at least) requires converting 
        #to "new" version of DataFrame to work
        cd <- do.call(rbind, lapply(EXTRACTOR(names$cd), DataFrame))
        #This is a patch for the Lohoff data due to SpatialExperiment changes
        #previously, sample_id was not required
        if(object.type == "SpatialExperiment" & 
            !"sample_id" %in% names(cd))
            cd$sample_id <- cd$embryo_pos_z
        data$colData <- cd
    }

    if(!is.null(names$dimred)){
        dr_list <- EXTRACTOR(names$dimred)
        dr_types <- names(dr_list[[1]])
        dr_sce <- lapply(dr_types, function(x){
            do.call(rbind, lapply(dr_list, function(y) y[[x]]))
        })
        names(dr_sce) <- dr_types
        data$reducedDims <- dr_sce
    }

    if(!is.null(names$coords)){
        data$spatialData <- do.call(rbind, EXTRACTOR(names$coords))
        coords <- c("x", "y", "z")
        data$spatialCoordsNames <- coords[coords %in% names(data$spatialData)]
    }

    command <- sprintf("%s(%s)",
        object.type,
        paste(sapply(names(data), function(x) paste0(x, "=data$", x)),
            collapse = ","))
    sce <- eval(parse(text = command))

    if(!is.null(names$sf)){
        sizeFactors(sce) <- do.call(c, EXTRACTOR(names$sf))
    }

    if("ENSEMBL" %in% names(rowData(sce))){
        rownames(sce) <- rowData(sce)$ENSEMBL
    }
    if("cell" %in% names(colData(sce))){
        colnames(sce) <- colData(sce)$cell
    }
    return(sce)
}

####
# Simpler interfaces for specific data types
####
.getRNAseqData <- function(dataset, type, version, samples, sample.options, sample.err, extra_assays=NULL){
    if(type == "processed"){ return(
        .getData(
            dataset,
            version,
            samples,
            sample.options,
            sample.err,
            names = list(
                assays=c(list("counts" = "counts-processed"), extra_assays),
                rd="rowdata",
                cd="coldata",
                sf="sizefac",
                dimred="reduced-dims"
            ),
            object.type="SingleCellExperiment"
        ))
    } else if (type == "raw"){ return(
        .getData(
            dataset,
            version,
            samples,
            sample.options,
            sample.err,
            names = list(
                assays=list("counts" = "counts-raw"),
                rd="rowdata"
            ),
            object.type="SingleCellExperiment",
            return.list=TRUE
        ))
    }
}

.fetch_version <- function(version, field) {
    opt <- version[[field]]
    if (is.null(opt)) {
        version[[1]]
    } else {
        opt
    }
}
