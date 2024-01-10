library(glue) ; g=glue ; len=length
library(matrixStats)
library(ggnewscale)
library(stringdist)
library(gridExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggrastr)
library(stringr)
library(Seurat)
library(dbscan)
library(future)
library(rlist)
library(dplyr)
library(purrr)
library(furrr)
library(rhdf5)
library(qpdf)
library(qs)

# setwd("~/pip/")

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/230317_SL-NSA_0564_AHVCHWDRX2/SI-TT-E11"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/230317_SL-NSA_0564_AHVCHWDRX2/SI-NT-C4"

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/231020_VL00181_99_AAF2KJGM5/SI-TS-A11"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/231020_VL00181_99_AAF2KJGM5/SI-TT-A1"

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/230310_SL-NVN_0914_AHCWK5DSX5/SI-TT-D5"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/230310_SL-NVN_0914_AHCWK5DSX5/SI-NT-D5"

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/231016_JW/SI-TT-C10"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/231016_JW/SI-NT-A1"

### Download files #############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript spatial.R RNApath SBpath", call. = FALSE)
}

RNApath <- args[1] ; print(g("RNApath: {RNApath}"))
SBpath <- args[2] ; print(g("SBpath: {SBpath}"))

stopifnot(!dir.exists("RNAcounts"))
stopifnot(!dir.exists("SBcounts"))

checkgsfile <- function(path) {return(system(g("gsutil ls {path}"),intern=F,ignore.stdout=T,ignore.stderr=T)==0)}
stopifnot(checkgsfile(RNApath)) ; stopifnot(checkgsfile(SBpath))

system("mkdir RNAcounts")
if (checkgsfile(file.path(RNApath,"outs/filtered_feature_bc_matrix.h5"))) {
  RNAtech = "cellranger count"
  system(g("gsutil cp {RNApath}/outs/filtered_feature_bc_matrix.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/raw_feature_bc_matrix.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/molecule_info.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/metrics_summary.csv RNAcounts"))
} else if (checkgsfile(file.path(RNApath,"outs/multi"))) {
  RNAtech = "cellranger multi"
  system(g("gsutil cp {RNApath}/outs/multi/count/raw_feature_bc_matrix.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/per_sample_outs/{basename(RNApath)}/count/sample_filtered_feature_bc_matrix.h5 RNAcounts/filtered_feature_bc_matrix.h5"))
  system(g("gsutil cp {RNApath}/outs/per_sample_outs/{basename(RNApath)}/count/sample_molecule_info.h5 RNAcounts/molecule_info.h5"))
  system(g("gsutil cp {RNApath}/outs/per_sample_outs/{basename(RNApath)}/metrics_summary.csv RNAcounts"))
} else {
  print("Unknown RNA directory structure, exiting...")
  stopifnot(F)
}

system("mkdir SBcounts")
system(g("gsutil cp {SBpath}/SBcounts.h5 SBcounts"))

# Folder requirements:
stopifnot(length(list.files("RNAcounts")) >= 1)
stopifnot(length(list.files("SBcounts")) == 1)

# Get the 10X dictionary
CBdictpath = "3M-february-2018.txt"
if (!file.exists(CBdictpath)) {
  system("gsutil cp gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/scripts/3M-february-2018.txt .")
}
stopifnot(file.exists(CBdictpath))

### Load the RNA data ##########################################################

# Load the RNA count matrix
obj <- "RNAcounts/filtered_feature_bc_matrix.h5" %>% Read10X_h5 %>% CreateSeuratObject

# Add metadata
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
obj[["logumi"]] <- log10(obj$nCount_RNA+1)
obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .))
obj[["type"]] <- "unknown"

# PCA, Cluster, and UMAP
obj %<>% Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose=F) %>%
    Seurat::FindNeighbors(dims=1:30) %>%
    Seurat::FindClusters(resolution=0.8) %>%
    Seurat::RunUMAP(dims=1:30, verbose=F, n.epochs=NULL)

# Add %intronic
if (file.exists("RNAcounts/molecule_info.h5")) {
  fetch <- function(x){return(h5read("RNAcounts/molecule_info.h5",x))}
  barcodes = fetch("barcodes")
  info = data.frame(barcode=fetch("barcode_idx")+1,umi_type=fetch("umi_type"))
  info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
  obj$pct.intronic = info$pct.intronic[match(obj$cb,barcodes[info$barcode])] * 100
  rm(info) ; rm(barcodes)
}

cb_whitelist = unname(obj$cb)
stopifnot(!duplicated(cb_whitelist))
stopifnot(len(unique(nchar(cb_whitelist)))==1)
stopifnot(map_lgl(strsplit(cb_whitelist,""),~all(.%in%c("A","C","G","T"))))

Misc(obj, "method") <- RNAtech ; rm(RNAtech)
Misc(obj, "called_cells") <- len(cb_whitelist)
#Misc(obj, "RNA_umi_total") <- sum(obj$nCount_RNA)
#Misc(obj, "RNA_umi_percell") <- sum(obj$nCount_RNA)/len(cb_whitelist)
Misc(obj, "RNA_path") <- RNApath ; rm(RNApath)
Misc(obj, "SB_path") <- SBpath ; rm(SBpath)

gc()

### Load the SB data ###########################################################

# Accessor method for the spatial .h5
f <- function(p){return(h5read("SBcounts/SBcounts.h5",p))}
cb_list = f("lists/cb_list") ; stopifnot(!any(duplicated(cb_list))) #; stopifnot(!grepl("N",cb_list))
sb_list = f("lists/sb_list") ; stopifnot(!any(duplicated(sb_list)))

# load the puck information
pucks = f("puck/puck_list") ; stopifnot(len(pucks)==1) # 2+ pucks not implemented yet
puckdf = data.frame(sb=f("puck/sb"), x=f("puck/x"), y=f("puck/y"), puck_index=f("puck/puck_index"))
stopifnot(unique(puckdf$puck_index)==1)
puckdf = puckdf[match(sb_list,puckdf$sb),]
stopifnot(puckdf$sb == sb_list)
bn = nrow(puckdf)
if (bn < 150000) {
  k = 0.73
} else if (bn < 600000) {
  k = 0.73 * 2
} else {
  k = 0.645
}
puckdf %<>% transmute(sb_index = 1:nrow(puckdf), x_um = x*k, y_um = y*k)
Misc(obj, "num_beads") <- bn
Misc(obj, "scaling_factor") <- k
rm(sb_list, bn, k)
gc()

# Load the SB count matrix
df = data.frame(cb_index=f("matrix/cb_index"),
                umi_2bit=f("matrix/umi"),
                sb_index=f("matrix/sb_index"),
                reads=f("matrix/reads"))

# Determine the whitelist remap status
cb_list_remap = f("lists/cb_list_remap")
reads_noremap = df %>% filter(cb_list[cb_index] %in% cb_whitelist) %>% pull(reads) %>% sum
reads_remap = df %>% filter(cb_list_remap[cb_index] %in% cb_whitelist) %>% pull(reads) %>% sum
remap = reads_remap > reads_noremap
rm(reads_remap) ; rm(reads_noremap) ; rm(cb_list_remap)

# Remap the whitelist
if (remap) {
  print("Remapping CB whitelist")
  cb_dict <- read.table(CBdictpath) %>% {setNames(.[[2]],.[[1]])}
  stopifnot(cb_whitelist %in% names(cb_dict))
  cb_whitelist = cb_dict[cb_whitelist] %>% unname
  stopifnot(!duplicated(cb_whitelist))
} ; rm(CBdictpath)

# Add metadata to the seurat object
Misc(obj, "puck") <- f("puck/puck_list")
Misc(obj, "R1s") <- f("metadata/R1s")
Misc(obj, "R2s") <- f("metadata/R2s")
Misc(obj, "switchR1R2") <- f("metadata/switch") %>% as.logical
Misc(obj, "remapCB") <- remap
Misc(obj, "UP_matching_type") <- f("metadata/UP_matching/type")
Misc(obj, "UP_matching_count") <- f("metadata/UP_matching/count")
Misc(obj, "SB_matching_type") <- f("metadata/SB_matching/type")
Misc(obj, "SB_matching_count") <- f("metadata/SB_matching/count")
Misc(obj, "SB_reads") <- f("metadata/num_reads")
Misc(obj, "SB_reads_filtered") <- sum(df$reads)

gc()

# Perform the HD1 matching between observed cell barcodes and the whitelist
listHD1neighbors <- function(input_string) {
  nucleotides <- c('A','C','G','T')
  result <- c()
  for (i in 1:nchar(input_string)) {
    current_char <- substr(input_string, i, i)
    for (nuc in nucleotides) {
      if (nuc != current_char) {
        new_string <- paste0(substr(input_string, 1, i - 1), nuc, substr(input_string, i + 1, nchar(input_string)))
        result <- c(result, new_string)
      }
    }
  }
  return(result)
}

print("Performing HD1 CB fuzzy matching")

# Exact matching dictionary
exact_dict = match(cb_list, cb_whitelist)

# HD1 fuzzy matching dictionary
neighbors = map(cb_whitelist,listHD1neighbors) %>% flatten_chr
originals = map(1:len(cb_whitelist),~rep(.,nchar(cb_whitelist[[1]])*3)) %>% flatten_int
m = is.na(match(neighbors, cb_list)) ; neighbors=neighbors[!m] ; originals=originals[!m] ; rm(m)
HD1ambig = unique(neighbors[duplicated(neighbors)])
m = neighbors%in%HD1ambig ; neighbors=neighbors[!m] ; originals=originals[!m] ; rm(m)
stopifnot(!any(duplicated(neighbors)))
fuzzy_dict = originals[match(cb_list, neighbors)]
rm(neighbors) ; rm(originals)
# stopifnot(table(!is.na(fuzzy_dict),!is.na(exact_dict))["TRUE","TRUE"] == 0) # exact and fuzzy matches should be distinct

# Perform matching
df %<>% mutate(exact = exact_dict[cb_index], HD1 = fuzzy_dict[cb_index])
rm(exact_dict) ; rm(fuzzy_dict)

# Write metadata
Misc(obj, "SB_reads_filtered_exact") <- df %>% filter(!is.na(exact)) %>% pull(reads) %>% sum
Misc(obj, "SB_reads_filtered_HD1") <- df %>% filter(!is.na(HD1)) %>% pull(reads) %>% sum
Misc(obj, "SB_reads_filtered_HD1ambig") <- df %>% filter(cb_index %in% match(HD1ambig,cb_list)) %>% pull(reads) %>% sum
rm(cb_list) ; rm(HD1ambig)
gc()
# stopifnot(table(!is.na(df$exact),!is.na(df$HD1))["TRUE","TRUE"] == 0) # exact and fuzzy matches should be distinct

# Remove duplicate rows introduced by fuzzy matching
df1 = df %>% filter(!is.na(exact)) %>% mutate(cb_index = exact) %>% select(1:4)
df2 = df %>% filter(!is.na(HD1)) %>% mutate(cb_index = HD1) %>% select(1:4)
df3 = df %>% filter(is.na(exact)&is.na(HD1)) %>% mutate(cb_index = -cb_index) %>% select(1:4)
df2 %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads=sum(reads)) %>% ungroup
df12 <- full_join(df1, df2, by = c("cb_index","umi_2bit","sb_index")) ; rm(df1) ; rm(df2)
df12$reads.x %<>% tidyr::replace_na(0) ; df12$reads.y %<>% tidyr::replace_na(0)
df12 %<>% mutate(reads=reads.x+reads.y) %>% select(-reads.x,-reads.y)
df = rbind(df12, df3) ; rm(df12) ; rm(df3) ; gc()

# Remove chimeric reads
print("Removing chimeras")
df %<>% arrange(cb_index,umi_2bit,desc(reads))
before_same = tidyr::replace_na(df$cb_index==lag(df$cb_index) & df$umi_2bit==lag(df$umi_2bit), FALSE) 
after_same = tidyr::replace_na(df$cb_index==lead(df$cb_index) & df$umi_2bit==lead(df$umi_2bit) & df$reads==lead(df$reads), FALSE)
chimeric = before_same | after_same

Misc(obj, "SB_reads_filtered_chimeric") <- df[chimeric,]$reads %>% sum
df = df[!chimeric,] ; rm(chimeric, before_same, after_same)

Misc(obj, "SB_reads_final") <- df %>% pull(reads) %>% sum

# remap the whitelist back
if (remap) {
  cb_whitelist = cb_dict[cb_whitelist] %>% unname
  stopifnot(!duplicated(cb_whitelist))
  rm(cb_dict)
}
rm(remap) ; gc()

count_umis <- function(df) {
  # This method implements the following logic, except faster and more memory efficient:
  # df %<>% group_by(cb_index, sb_index) %>% summarize(umi=n()) 
  
  gdf = df %>% select(cb_index,sb_index) %>% arrange(cb_index, sb_index) 
  bnds = (gdf$cb_index!=lead(gdf$cb_index) | gdf$sb_index!=lead(gdf$sb_index)) %>% tidyr::replace_na(T) %>% which
  gdf %<>% distinct()
  gdf$umi = (bnds - lag(bnds)) %>% tidyr::replace_na(bnds[[1]]) ; rm(bnds)
  
  gdf %<>% arrange(desc(umi))

  return(gdf)
}

# Make a plot of the distributions
plotSBumicurves <- function(df) {
  gdf = count_umis(df)
  
  cb.data = gdf %>% group_by(cb_index) %>% dplyr::summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="all cell barcodes")}
  cb.data2 = cb.data %>% filter(cb_index>0) %>% {mutate(.,index=1:nrow(.),filter="called cell barcodes only")}
  sb_pct_in_called_cells = round(sum(filter(cb.data,cb_index>0)$umi)/sum(cb.data$umi)*100,2)
  p1 = ggplot(mapping=aes(x=index, y=umi,col=filter))+geom_line(data=cb.data)+geom_line(data=cb.data2) +
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per cell")+ylab("SB UMI counts")+xlab("Cell barcodes") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines")) +
    annotate("text", x = Inf, y = Inf, label = g("SB UMI in called cells: {sb_pct_in_called_cells}%"), hjust = 1, vjust = 1.3)
  rm(cb.data) ; rm(cb.data2) ; gc()
  
  sb.data = gdf %>% group_by(sb_index) %>% dplyr::summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="all cell barcodes")}
  sb.data2 = gdf %>% filter(cb_index > 0) %>% group_by(sb_index) %>% dplyr::summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="called cell barcodes only")}
  p2 = ggplot(mapping=aes(x=index,y=umi,col=filter))+geom_line(data=sb.data)+geom_line(data=sb.data2)+
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per bead")+ylab("SB UMI counts")+xlab("Beads")+
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
  rm(sb.data) ; rm(sb.data2) ; gc()
  
  return(list(p1, p2))
}
umicurves = plotSBumicurves(df)

# Compute metrics
gdf = df %>% count_umis
Misc(obj, "SB_umi_final") <- gdf %>% pull(umi) %>% sum
Misc(obj, "SB_umi_final_called") <- gdf %>% filter(cb_index>0) %>% pull(umi) %>% sum
Misc(obj, "SB_umi_final_uncalled") <- gdf %>% filter(cb_index<0) %>% pull(umi) %>% sum
Misc(obj, "SB_umi_final_pctcalled") <- round(Misc(obj, "SB_umi_final_called")/Misc(obj, "SB_umi_final")*100,2)
Misc(obj, "SB_umi_filtered_downsampling") <- f("metadata/downsampling")
rm(gdf) ; gc()

# remove reads that didn't match a called cell
df %<>% filter(cb_index>0)

gc()

### Positioning methods ########################################################

chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}

# Do a grid search to find the ideal DBSCAN parameters
ncores = 20L ; plan(multisession, workers=ncores)
opt_dbscan <- function(data.list) {
  eps.vec = c(50) ; minPts.vec = c(3:42)
  res = data.frame() ; i = 0
  repeat{
    params = expand.grid(eps.vec,minPts.vec) %>% setNames(c("eps","minPts"))
    row_lists = chunk_vector(1:nrow(params), round(nrow(params)/ncores))

    params$pct = furrr::future_map(row_lists, function(v) {
      map_dbl(v, function(i) {
        m = map_lgl(data.list, ~max(dbscan::dbscan(.[c("x_um","y_um")], eps=params$eps[[i]], minPts=params$minPts[[i]], weights=.$umi)$cluster) == 1)
        return(sum(m)/length(m))
      })
    }, .options=furrr_options(seed=T)) %>% flatten_dbl
    
    res = rbind(res, params)
    
    if (which.max(res$pct)<0.9*nrow(res) || i >= 26) {
      break
    }
    
    minPts.vec = minPts.vec + 40
    i = i + 1
  }
  
  params = res ; rm(res)
  params$is.max = params$pct==max(params$pct)
  eps = params$eps[params$is.max][[1]] ; minPts = params$minPts[params$is.max][[1]]
  print(g("Optimal eps: {eps}    Optimal minPts: {minPts}    %placed: {round(max(params$pct)*100,2)}"))
  
  return(c(eps,minPts))
}

# Add the DBSCAN clusters to the dataframes
run_dbscan <- function(data.list, eps, minPts) {
  lapply(data.list, function(df){
    df$cluster <- dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$umi)$cluster
    return(df)
  })
}

# assign centroid and record metadata
create_coords <- function(data.list) {
  lapply(data.list, function(df) {
    p = c(x_um=NA,
          y_um=NA,
          DBSCAN_clusters=max(df$cluster),
          num_SBumi = sum(df$origumi),
          SNR=NA,
          SB_bin = unique(df$bin) %>% {ifelse(is.null(.), NA, .)},
          minPts = unique(df$minPts) %>% {ifelse(is.null(.), NA, .)},
          eps = unique(df$eps) %>% {ifelse(is.null(.), NA, .)})
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p[["x_um"]] = matrixStats::weightedMedian(sdf$x_um,w=sdf$umi)
      p[["y_um"]] = matrixStats::weightedMedian(sdf$y_um,w=sdf$umi)
      p[["SNR"]] = sum(sdf$origumi)/sum(df$origumi)
    }
    return(p)
  }) %>% bind_rows %>% as.data.frame %>% mutate(cb_index=as.numeric(names(data.list))) %>% select(cb_index, everything())
}

# Run these methods on the entire data
normal_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  params = opt_dbscan(data.list)
  data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
  coords <- create_coords(data.list)
  coords %<>% mutate(eps=params[[1]],minPts=params[[2]])
  return(coords)
}

# Split the data into 10 SB UMI buckets and run on each
binned_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  
  # create the deciles
  quants = c(0,quantile(map_dbl(data.list,~sum(.$umi)), probs = seq(0.1, 1, by = .1))) %>% unname
  paste("Deciles: ", paste(round(quants),collapse=" "))
  umicounts = map(data.list,~sum(.$umi))
  data.lists=map(1:(len(quants)-1),~data.list[umicounts>quants[.] & umicounts<=quants[.+1]])
  
  # run positioning on each decile
  data.list = map2(data.lists,quants[-1],function(data.list,quant) {
    if (len(data.list) == 0) {print(g({"skipping quantile, empty"}));return(list())}
    params = opt_dbscan(data.list)
    data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
    data.list %<>% map(~mutate(.,bin=quant,eps=params[[1]],minPts=params[[2]]))
    return(data.list)
  }) %>% list_flatten
  
  stopifnot(len(data.list)==len(unique(df$cb_index)))
  
  coords <- create_coords(data.list)
  return(coords)
}

### Positioning ################################################################

correct_beadloss <- function(df) {
  df %<>% group_by(sb_index) %>% 
    mutate(origumi = umi, totumi=sum(umi), umi=ifelse(totumi>256, umi/totumi*256, umi)) %>% 
    ungroup %>% 
    select(-totumi) %>% 
    arrange(desc(umi))
}

# Perform positioning at various levels of downsampling
original_df <- df
coords_list = list()
for (i in seq(0.05,1,0.05)) {
  # Downsample the reads
  print(g("Downsampling: {round(i*100)}%"))
  df = original_df %>% mutate(reads=rmultinom(n=1, size=round(sum(original_df$reads)*i), prob=original_df$reads) %>% as.vector)
  df %<>% filter(reads>0)
  df %<>% count_umis
  
  # correct for bead loss (downsample beads above 256 umi)
  df %<>% correct_beadloss
  
  # add spatial positions from puck
  df = merge(x=df,y=puckdf,all.x=T,by="sb_index")
  
  # run binned positioning (see method above)
  data.list = split(df, df$cb_index)
  quants = c(0,quantile(map_dbl(data.list,~sum(.$umi)), probs = seq(0.1, 1, by = .1))) %>% unname
  paste("Deciles: ", paste(round(quants),collapse=" "))
  umicounts = map(data.list,~sum(.$umi))
  data.lists=map(1:(len(quants)-1),~data.list[umicounts>quants[.] & umicounts<=quants[.+1]])
  data.list = map2(data.lists,quants[-1],function(data.list,quant) {
    if (len(data.list) == 0) {print(g({"skipping quantile, empty"}));return(list())}
    params = opt_dbscan(data.list)
    data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
    data.list %<>% map(~mutate(.,bin=quant,eps=params[[1]],minPts=params[[2]]))
    return(data.list)
  }) %>% list_flatten
  stopifnot(len(data.list)==len(unique(df$cb_index)))
  coords <- create_coords(data.list)
  
  coords_list %<>% list.append(coords)
  gc()
}

# clean up
# rm(original_df) ; rm(puckdf)

# merge with seurat object
coords %<>% mutate(cb=cb_whitelist[cb_index])
rownames(coords) = paste0(coords$cb,"-1")
obj = AddMetaData(obj,coords)

emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_")

coords %>% select(cb,everything()) %>% select(-cb_index) %>% write.csv("coords.csv",quote=F,row.names=F)
qsave(obj, "seurat.qs")

################################################################################
### CREATE PDF #################################################################
################################################################################

print("Making summary.pdf")

system("mkdir plots")

fetch <- function(x){return(h5read("RNAcounts/molecule_info.h5",x))}
plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}
make.pdf <- function(plots,name,w,h) {
  if (any(class(plots)=="gg")||class(plots)=="Heatmap") {plots=list(plots)}
  pdf(file=name,width=w,height=h)
  lapply(plots,function(x){print(x)})
  dev.off()
}

### Page 0: cell ranger output #################################################

if (file.exists("RNAcounts/metrics_summary.csv")) {
  plotdf = read.table("RNAcounts/metrics_summary.csv",header=F,comment.char="",sep=",")
  if (nrow(plotdf)==2) { # count
    plotdf %<>% t
  } else if (ncol(plotdf)==6) { # multi
    colnames(plotdf) = as.character(plotdf[1,])
    plotdf = plotdf[-1,c(5,6)]
  }
  rownames(plotdf) = NULL
  plot = plot_grid(ggdraw()+draw_label(""),
                   ggdraw()+draw_label("Cell Ranger Metrics Summary"),
                   plot.tab(plotdf),
                   ggdraw()+draw_label(""),
                   ncol=1,rel_heights=c(0.1,0.1,0.7,0.2))
  make.pdf(plot,"plots/0cellranger.pdf",7,8)
}

# Add metadata to seurat object
for (i in 1:nrow(plotdf)) {Misc(obj,plotdf[i,1]) <- plotdf[i,2]}

### Page 1: cell calling #######################################################

UvsI <- function(obj) {
  barcodes = fetch("barcodes")
  fetch <- function(x){return(h5read(g("RNAcounts/molecule_info.h5"),x))}
  molecule_info = data.frame(barcode=fetch("barcode_idx"),
                             umi_type=fetch("umi_type"),
                             reads=fetch("count"))
  
  # Panel 1: downsampling curve
  tab = table(molecule_info$reads)
  downsampling = map_int(seq(0,1,0.05),function(p){sum(map2_int(tab, as.numeric(names(tab)), function(v,k){length(unique(floor(sample(0:(k*v-1), round(k*v*p), replace=F)/k)))}))})
  plotdf = data.frame(x=seq(0,1,0.05)*sum(molecule_info$reads)/1000/1000,y=downsampling/1000/1000)
  p0 = ggplot(plotdf, aes(x=x,y=y))+geom_line()+theme_bw()+xlab("Millions of reads")+ylab("Millions of filtered UMIs")+ggtitle("RNA Downsampling curve")
  
  df = molecule_info %>% group_by(barcode) %>% summarize(umi=n(), pct.intronic=sum(umi_type==0)/umi) %>% arrange(desc(umi)) %>% mutate(logumi=log10(umi))
  
  # Panel 2 and 3: intronic density
  if (!all(df$pct.intronic==0)) {
    ct = 500
    if (any(df$umi>=ct)) {
      p1 = df %>% filter(umi>=ct) %>% ggplot(aes(x = logumi, y = pct.intronic)) + 
        geom_bin2d(bins=100) +
        scale_fill_viridis(trans="log", option="A", name="density") + 
        theme_minimal() +
        labs(title = g("Intronic vs. UMI droplets (>{ct} umi)"), x = "logumi", y = "%intronic") & NoLegend()
      
      max_density_x = density(filter(df,umi>=ct)$pct.intronic) %>% {.$x[which.max(.$y)]}
      p2 = df %>% filter(umi>=ct) %>% ggplot(aes(x = pct.intronic)) +
        geom_density() + 
        theme_minimal() +
        labs(title = g("Intronic density (>{ct} umi)"), x = "%intronic", y = "Density") + 
        geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
        annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
    } else {
      p1 = ggdraw()+draw_label(g("No cells with {ct}+ UMI"))
      p2 = ggdraw()+draw_label(g("No cells with {ct}+ UMI"))
    }
  } else {
    p1 = ggdraw()+draw_label("No intronic information")
    p2 = ggdraw()+draw_label("No intronic information")
  }

  # Panel 4: cell barcode knee plot
  df %<>% mutate(index=1:nrow(df), called=barcodes[barcode+1]%in%cb_whitelist)
  p3 = ggplot(df,aes(x=index,y=umi,col=called))+geom_line()+theme_bw()+scale_x_log10()+scale_y_log10()+
    ggtitle("Barcode rank plot")+xlab("Barcodes")+ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))#, legend.title=element_text(size=10), legend.text=element_text(size=8), legend.margin=margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"), legend.key.size = unit(0.5, "lines"))
  
  plot = plot_grid(p3,p1,p0,p2,ncol=2)
  
  return(plot)
}

if (file.exists("RNAcounts/molecule_info.h5")) {
  plot=UvsI(obj)
  make.pdf(plot,"plots/1cellcalling.pdf",7,8)
  gc()
}

### Page 2: UMAP + metrics #####################################################

plot = plot_grid(DimPlot(obj,label=T)+ggtitle(g("UMAP"))+NoLegend()+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())+coord_fixed(ratio=1),
                 VlnPlot(obj,"logumi")+NoLegend()+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank()),
                 FeaturePlot(obj,"percent.mt")+ggtitle("%MT")+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())+coord_fixed(ratio=1)+theme(legend.position="top",legend.justification="center",legend.key.width=unit(2, "lines")),
                 FeaturePlot(obj,"pct.intronic")+ggtitle("%Intronic")+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())+coord_fixed(ratio=1)+theme(legend.position="top",legend.justification="center",legend.key.width=unit(2,"lines")),
                 ncol=2)
make.pdf(plot,"plots/2umap.pdf",7,8)

### Page 3: FeaturePlots #####################################################

human_genes = c("MBP","PLP1","CSPG5","VCAN","AQP4","GFAP","FLT1","DCN","CX3CR1","GAD1", "GAD2", "RELN","SLC17A6","SLC17A7", "TH", "NR4A2")
mouse_genes = c("Slc17a7", "Slc17a6", "Aqp4","Plp1","Rgs5","Matn4","Resp18","Nrp2","Meis2")

if (sum(human_genes%in%rownames(obj)>0)) {
  print("Human genes detected")
  genes=human_genes
} else if (sum(mouse_genes%in%rownames(obj)>0)) {
  print("Mouse genes detected")
  genes=mouse_genes
} else {
  print("No mouse or human genes detected")
  genes = c()
}

if (len(genes) > 0) {
  plot1 = FeaturePlot(obj,genes)&theme_void()+
    theme(legend.key.width=unit(0.3,"cm"), plot.title=element_text(hjust=0.5)) & coord_fixed(ratio=1)
  plot2 = FeaturePlot(obj,genes,reduction="spatial")&theme_void()+
    theme(legend.key.width=unit(0.3,"cm"), plot.title=element_text(hjust=0.5)) & coord_fixed(ratio=1)
} else {
  plot1 = ggdraw()+draw_label("No human or mouse marker genes found")
  plot2 = ggdraw()+draw_label("No human or mouse marker genes found")
}
make.pdf(plot1,"plots/3features.pdf",7,8)
make.pdf(plot2,"plots/8sfeatures.pdf",7,8)

### Page 4: Raw spatial data ###################################################

# Panel 1+2: elbow plots
p1 = umicurves[[1]]
p2 = umicurves[[2]]

# Panel 3: Downsampling UMI
x = seq(0,1,0.05)*f("metadata/num_reads")/1000000
plot.df = data.frame(x=x,y=f("metadata/downsampling")/1000000)
p3 = ggplot(plot.df, aes(x=x,y=y))+geom_point()+theme_bw()+xlab("Reads (millions)")+ylab("SB-whitelist UMIs (millions)")+ggtitle("Downsampling curve")

# Panel 4: Downsampling placement
plot.df = map(coords_list,function(df){return(df$DBSCAN_clusters %>% {c(sum(.==0),sum(.==1),sum(.>=2))/nrow(df)*100})}) %>% {do.call(rbind,.)} %>% {rbind(c(0,0,0),.)}
plot.df %<>% as.data.frame %>% setNames(c("0","1","2+")) %>% mutate(x=x)
plot.df = tidyr::gather(plot.df, key = "column", value = "value", -x)

mastercoord = coords %>% transmute(cb=cb_index,x1=x_um,y1=y_um)
rmse = coords_list %>% head(-1) %>% map_dbl(function(coord){
  coord %<>% transmute(cb=cb_index,x2=x_um,y2=y_um)
  df = merge(mastercoord,coord,by="cb")
  df %<>% filter(!is.na(x1),!is.na(x2))
  df %<>% mutate(dist=sqrt((x2-x1)^2+(y2-y1)^2))
  return(mean(df$dist))
}) %>% {c(NA,.,NA)}
#plot.df2 = data.frame(x=x,rmse=rmse)
m1 = max(plot.df$value, na.rm=T)
m2 = max(rmse, na.rm=T)
#plot.df2 %<>% mutate(rmse=rmse/m2*m1)
plot.df %<>% rbind(data.frame(x=x,column="disp.",value=rmse/m2*m1))

p4 <- ggplot(plot.df,aes(x=x, y=value, color=column)) + geom_line() +
  theme_bw() +
  #geom_line(data=plot.df, aes(x=x, y=value, color=column)) + 
  #geom_line(data=plot.df2, aes(x=x, y=rmse, color="grey"), col="grey") +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "grey")) +
  scale_y_continuous(sec.axis = sec_axis(~ . * m2/m1, name = "average displacement (µm)")) +
  labs(title = "Downsampling Placements", x = "Reads (millions)", y = "Percent placed", color = "") +
  theme(legend.position=c(0.85, 0.85), legend.background=element_blank(), legend.key=element_blank(), legend.key.height=unit(0.75, "lines"))

plot = plot_grid(p1,p2,p3,p4,ncol=2)
make.pdf(plot,"plots/4rawspatial.pdf",7,8)

### Page 5: Beadplot ###########################################################

df = original_df %>% count_umis %>% correct_beadloss %>% merge(y=puckdf,all.x=T,by="sb_index")
sb.data = df %>% group_by(sb_index) %>% dplyr::summarize(umi=sum(origumi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]}
sb.data_weighted = df %>% group_by(sb_index) %>% dplyr::summarize(umi=sum(umi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]}

beadplot <- function(sb.data, m, text){
  sb.data = sb.data[nrow(sb.data):1,]
  ggplot(sb.data, aes(x=x_um,y=y_um,col=umi)) +
    rasterize(geom_point(size=0.1), dpi=200) +
    coord_fixed() +
    theme_classic() +
    labs(x="x (um)", y="y (um)") +
    scale_color_viridis(trans="log", option="B", name="UMI", limits = c(1, m)) + 
    ggtitle(g("SB UMI per bead ({text})"))
}
p1 = beadplot(sb.data, max(sb.data$umi), "raw")
p2 = beadplot(sb.data_weighted, max(sb.data$umi), "corrected")

plot = plot_grid(p1,p2,ncol=1)
make.pdf(plot,"plots/5beadplot.pdf",7,8)

### Page 6: DBSCAN #############################################################

# Panel 1: DBSCAN cluster distribution
d=data.frame(x=coords$DBSCAN_clusters)
p1 = ggplot(d,aes(x=x)) + geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))*100), binwidth=.5) +
  geom_text(aes(label = sprintf("%1.0f%%", after_stat(count)/sum(after_stat(count))*100), y=after_stat(count)/sum(after_stat(count))*100), stat="bin", binwidth=1, vjust=-0.5)+
  theme_classic() + xlab("Num DBSCAN clusters") + ylab("Percent") +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_continuous(breaks=min(d$x):max(d$x)) +
  ggtitle("DBSCAN cluster distribution")

# Panel 2: SNR density
max_density_x = density(obj$SNR %>% na.omit) %>% {.$x[which.max(.$y)]}
p2 = obj@meta.data %>% filter(!is.na(x_um)) %>% ggplot(aes(x = SNR)) +
  geom_density() + 
  theme_minimal() +
  labs(title = "SNR density", x = "SNR", y = "Density") + 
  geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
  annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")

# Panel 3: RNA umi vs SB umi
p3 = data.frame(x=obj$nCount_RNA,y=obj$num_SBumi,placed=!is.na(obj$x_um)) %>% {ggplot(.,aes(x=log10(x),y=log10(y),col=placed))+geom_point(size=0.2)+theme_bw()+xlab("RNA UMI")+ylab("SB UMI")+ggtitle("SB UMI vs. RNA UMI")+theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.background = element_blank(), legend.title=element_text(size=10), legend.text=element_text(size=8), legend.margin=margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"), legend.spacing.y = unit(0.1,"lines"), legend.key.size = unit(0.5, "lines"))}

# Panel 4: DBSCAN parameters
p4 = coords %>% select(SB_bin,minPts,eps) %>% distinct %>% arrange(SB_bin) %>% mutate(SB_bin=round(SB_bin,2)) %>% plot.tab

plot = plot_grid(p1,p2,p3,p4,ncol=2)
make.pdf(plot,"plots/6DBSCAN.pdf",7,8)

### Page 7: Spatial ############################################################

p1 = DimPlot(obj,reduction="spatial")+coord_fixed(ratio=1)+ggtitle(g("%placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)} ({sum(!is.na(obj$x_um))}/{ncol(obj)})")) + NoLegend() + xlab("x-position (µm)") + ylab("y-position (µm)")
p2 = DimPlot(obj, reduction="spatial",split.by="seurat_clusters",ncol=7) + theme_void() + coord_fixed(ratio=1) + NoLegend()
plot = plot_grid(p1, p2, ncol=1, rel_heights=c(1,1))
make.pdf(plot,"plots/7spatial.pdf",7,8)

### Page 9: Create metrics plot ################################################

df = list(
  c("Reads",Misc(obj,"SB_reads")),
  c("Puck file",Misc(obj,"puck")),
  c("Number of beads",Misc(obj,"num_beads")),
  c("Scaling factor",Misc(obj,"scaling_factor")),
  c("R1<->R2",Misc(obj,"switchR1R2")),
  c("Remap CB",Misc(obj,"remapCB"))
) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("metric","value"))
p1 = plot_grid(ggdraw()+draw_label("SB metrics"), plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

UP_matching <- setNames(Misc(obj,"UP_matching_count"),Misc(obj,"UP_matching_type"))
SB_matching <- setNames(Misc(obj,"SB_matching_count"),Misc(obj,"SB_matching_type"))
CB_matching <- setNames(map_int(c("SB_reads_filtered_exact","SB_reads_filtered_HD1","SB_reads_filtered_HD1ambig"),~Misc(obj,.)), c("exact","fuzzy","ambig"))
CB_matching["none"] <- Misc(obj,"SB_reads_filtered") - sum(CB_matching)

df = data.frame(a=c("exact","fuzzy", "none", "GG"),b=c(UP_matching[["-"]],sum(UP_matching[c("1D-","1D-1X","-1X","-1D","-2X")]),UP_matching[["none"]],UP_matching[["GG"]]) %>% {./sum(.)*100} %>% round(2) %>% paste0("%")) %>% arrange(desc(b)) %>% unname
p2 = plot_grid(ggdraw()+draw_label("UP matching"), plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

df = data.frame(a=c("exact","fuzzy","none","ambig"),b=SB_matching[c("exact","HD1","none","HD1ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
p3 = plot_grid(ggdraw()+draw_label("SB matching"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

df = data.frame(a=c("exact","fuzzy","none","ambig"),b=CB_matching[c("exact","fuzzy","none","ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
p4 = plot_grid(ggdraw()+draw_label("CB matching"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

df = list(
  c("UMI filter",UP_matching[["R1lowQ"]]),
  c("GG seq",UP_matching[["GG"]]),
  c("no UP",UP_matching[["none"]]),
  c("no SB",SB_matching[["none"]]),
  c("ambig SB",SB_matching[["HD1ambig"]]),
  c("chimeric",Misc(obj,"SB_reads_filtered_chimeric")),
  c("no CB",sum(CB_matching[c("none","ambig")]))
) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("filter","%removed"))
df$`%removed` = as.numeric(df$`%removed`) / (Misc(obj,"SB_reads")-cumsum(df$`%removed`))
df$`%removed`= round(df$`%removed` * 100,2) %>% paste0("%")
p5 = plot_grid(ggdraw()+draw_label("SB filtering"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

# Script input parameters (with common prefix removed)
a=Misc(obj,"RNA_path") ; b=Misc(obj,"SB_path")
m=min(which(!map_lgl(1:min(nchar(a),nchar(b)), ~str_sub(a,1,.)==str_sub(b,1,.))))
a%<>%str_sub(m-1,999) ; b%<>%str_sub(m-1,999)
df=data.frame(a=c("RNA_path","SB_path"),b=c(a,b)) %>% unname
p6 = plot_grid(ggdraw()+draw_label("Run parameters"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

plot = plot_grid(
  ggdraw()+draw_label(""), #spacer
  plot_grid(p1,p5,ncol=2),
  plot_grid(p2,p3,p4,ncol=3),
  p6,
  ggdraw()+draw_label(""), #spacer
  ncol=1,
  rel_heights = c(0.27,0.6,0.35,0.25,0.27)
)

make.pdf(plot,"plots/9metrics.pdf",7,8)

### Sample bead plots ##########################################################

cols = list(c("steelblue4","steelblue1"),c("red4","red1"),c("springgreen4","springgreen1"))
plot.sb <- function(subdf) {
  limits = range(subdf$origumi)
  nlegend = min(3,max(subdf$cluster))
  pts = lapply(1:nlegend,function(i){
    if (i <= 2) {
      list(geom_point(data=filter(subdf,cluster==i),mapping=aes(x=x_um,y=y_um,col=origumi),size=0.5), scale_color_gradient(low=cols[[i]][[1]],high=cols[[i]][[2]],limits=limits,breaks=limits), new_scale_color())
    } else {
      list(geom_point(data=filter(subdf,cluster>=3),mapping=aes(x=x_um,y=y_um,col=origumi),size=0.5), scale_color_gradient(low=cols[[i]][[1]],high=cols[[i]][[2]],limits=limits,breaks=limits), new_scale_color())
    }
  })
  Reduce(`+`,pts,init=ggplot())+coord_fixed(ratio=1,xlim=range(df$x_um),ylim=range(df$y_um))+theme_void()+
    theme(legend.key.width=unit(0.5,"lines"), legend.position="right", legend.key.height=unit(0.3,"lines"), legend.title=element_blank(), legend.spacing.y=unit(0.2,"lines"), legend.margin=margin(0,0,0,0,"lines"), legend.box.margin=margin(0,0,0,0,"pt"), legend.box.background=element_blank(), legend.background=element_blank(), legend.direction="vertical", legend.justification="left",legend.box.just="left",legend.box.spacing=unit(0,"cm"))
}

list0 = data.list %>% keep(~max(.$cluster)==0) %>% map(~mutate(.,cluster=cluster+1) %>% arrange(origumi))
list1 = data.list %>% keep(~max(.$cluster)==1) %>% map(~mutate(.,cluster=cluster+1) %>% arrange(origumi))
list2 = data.list %>% keep(~max(.$cluster)>=2) %>% map(~mutate(.,cluster=cluster+1) %>% arrange(origumi))

p1 = plot_grid(ggdraw()+draw_label("DBSCAN=0"), map(sample(list0,20,replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=4)}, ncol=1,rel_heights=c(0.1,2))
p2 = plot_grid(ggdraw()+draw_label("DBSCAN=1"), map(sample(list1,20,replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=4)}, ncol=1,rel_heights=c(0.1,2))
p3 = plot_grid(ggdraw()+draw_label("DBSCAN>=2"), map(sample(list2,20,replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=4)}, ncol=1,rel_heights=c(0.1,2))

make.pdf(list(p1,p2,p3),"plots/SB.pdf",7,7)

# legend.key.size
# legend.text
# legend.text.align

### Save output ################################################################

pdfs = c("0cellranger.pdf","1cellcalling.pdf", "2umap.pdf", "3features.pdf", "4rawspatial.pdf", "5beadplot.pdf", "6DBSCAN.pdf","7spatial.pdf","8sfeatures.pdf","9metrics.pdf","SB.pdf") %>% paste0("plots/",.)
qpdf::pdf_combine(input = pdfs, output = "summary.pdf")
qsave(obj, "seurat.qs") # we added some more metadata while plotting

# Adjust the plotting return
# Adjust CB matching storage
# Adjust origumi and umi_corrected

print("Done!")
