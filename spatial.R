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

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/230411_VL00297_168_AACMJCLM5/SI-TT-A9"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/230411_VL00297_168_AACMJCLM5/SI-NT-A9"

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
f <- function(p){return(h5read("SBcounts/SBcounts.h5",p))}

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
obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .))
obj[["logumi"]] <- log10(obj$nCount_RNA+1)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")

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
Misc(obj, "RNA_path") <- RNApath ; rm(RNApath)
Misc(obj, "SB_path") <- SBpath ; rm(SBpath)

gc()

### Load the SB data ###########################################################

# load the puck information
load_puck <- function() {
  puckdf = data.frame(sb=f("puck/sb"), x=f("puck/x"), y=f("puck/y"), puck_index=f("puck/puck_index"))
  pucks = f("puck/puck_list") ; stopifnot(len(pucks)==1) # 2+ pucks not implemented yet
  stopifnot(unique(puckdf$puck_index)==1)
  
  sb_list = f("lists/sb_list") ; stopifnot(!any(duplicated(sb_list)))
  puckdf$sb_index = match(puckdf$sb,sb_list)
  puckdf %<>% select(sb_index,x,y)
  return(puckdf)
}
puckdf <- load_puck()

# scale the puck coordinates (estimate using number of beads)
bn = nrow(puckdf)
if (bn < 150000) {
  k = 0.73
} else if (bn < 600000) {
  k = 0.73 * 2
} else {
  k = 0.645
}
puckdf %<>% transmute(sb_index=sb_index, x_um=x*k, y_um=y*k)
Misc(obj, "num_beads") <- bn ; rm(bn)
Misc(obj, "scaling_factor") <- k ; rm(k)

# Load the SB count matrix
df = data.frame(cb_index=f("matrix/cb_index"),
                umi_2bit=f("matrix/umi"),
                sb_index=f("matrix/sb_index"),
                reads=f("matrix/reads"))
cb_list = f("lists/cb_list") ; stopifnot(!any(duplicated(cb_list)))

# Remap cb_whitelist
determine_remap <- function(df, cb_list, cb_whitelist) {
  cb_list_remap = f("lists/cb_list_remap")
  reads_noremap = df %>% filter(cb_list[cb_index] %in% cb_whitelist) %>% pull(reads) %>% sum
  reads_remap = df %>% filter(cb_list_remap[cb_index] %in% cb_whitelist) %>% pull(reads) %>% sum
  remap = reads_remap > reads_noremap
  return(remap)
}
remap <- determine_remap(df, cb_list, cb_whitelist)
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

# Fuzzy match and convert cb_index from a cb_list index to a cb_whitelist index
fuzzy_matching <- function(df, cb_list, cb_whitelist) {
  print("Performing HD1 CB fuzzy matching")
  listHD1neighbors <- function(input_string) {
    nucleotides <- c('A','C','G','T','N')
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
  
  # Exact matching dictionary
  exact_dict = match(cb_list, cb_whitelist)
  
  # HD1 fuzzy matching dictionary
  neighbors = map(cb_whitelist,listHD1neighbors) %>% flatten_chr
  originals = map(1:len(cb_whitelist),~rep(.,nchar(cb_whitelist[[1]])*4)) %>% flatten_int
  m = (neighbors %in% cb_list) ; neighbors=neighbors[m] ; originals=originals[m] ; rm(m)
  m = !(neighbors %in% cb_whitelist) ; neighbors=neighbors[m] ; originals=originals[m] ; rm(m)
  HD1ambig = unique(neighbors[duplicated(neighbors)])
  m = !(neighbors %in% HD1ambig) ; neighbors=neighbors[m] ; originals=originals[m] ; rm(m)
  stopifnot(!any(duplicated(neighbors)))
  fuzzy_dict = originals[match(cb_list, neighbors)]
  HD1ambig_dict = !is.na(match(cb_list, HD1ambig))
  rm(neighbors) ; rm(originals) ; rm(HD1ambig)
  
  # Perform matching
  df %<>% mutate(exact = exact_dict[cb_index], HD1 = fuzzy_dict[cb_index], HD1ambig = HD1ambig_dict[cb_index])
  rm(exact_dict) ; rm(fuzzy_dict) ; rm(HD1ambig_dict)
  
  # Record metadata
  CB_matching_type = c("exact", "HD1", "HD1ambig","none")
  CB_matching_count = c(df %>% filter(!is.na(exact)) %>% pull(reads) %>% sum,
                        df %>% filter(!is.na(HD1)) %>% pull(reads) %>% sum,
                        df %>% filter(HD1ambig) %>% pull(reads) %>% sum,
                        df %>% filter(is.na(exact),is.na(HD1),!HD1ambig) %>% pull(reads) %>% sum)
  stopifnot(Misc(obj,"SB_reads_filtered") == sum(CB_matching_count))
  
  # Perform the cb_index conversion
  df1 = df %>% filter(!is.na(exact)) %>% mutate(cb_index = exact) %>% select(1:4)
  df2 = df %>% filter(!is.na(HD1)) %>% mutate(cb_index = HD1) %>% select(1:4)
  df3 = df %>% filter(is.na(exact)&is.na(HD1)) %>% mutate(cb_index = -cb_index) %>% select(1:4)
  df2 %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads=sum(reads)) %>% ungroup
  df12 <- full_join(df1, df2, by = c("cb_index","umi_2bit","sb_index")) ; rm(df1) ; rm(df2)
  df12$reads.x %<>% tidyr::replace_na(0) ; df12$reads.y %<>% tidyr::replace_na(0)
  df12 %<>% mutate(reads=reads.x+reads.y) %>% select(-reads.x,-reads.y)
  df = rbind(df12, df3) ; rm(df12) ; rm(df3) ; gc()
  
  res=list(df, CB_matching_type, CB_matching_count)
  return(res)
}
res <- fuzzy_matching(df, cb_list, cb_whitelist)
df <- res[[1]]
Misc(obj,"CB_matching_type") = res[[2]]
Misc(obj,"CB_matching_count") = res[[3]]
rm(cb_list) ; rm(res) ; gc()

# Remove chimeric reads
print("Removing chimeras")
df %<>% arrange(cb_index,umi_2bit,desc(reads))
before_same = tidyr::replace_na(df$cb_index==lag(df$cb_index) & df$umi_2bit==lag(df$umi_2bit), FALSE) 
after_same = tidyr::replace_na(df$cb_index==lead(df$cb_index) & df$umi_2bit==lead(df$umi_2bit) & df$reads==lead(df$reads), FALSE)
chimeric = before_same | after_same
Misc(obj, "SB_reads_filtered_chimeric") <- df[chimeric,]$reads %>% sum
df = df[!chimeric,] ; rm(chimeric, before_same, after_same)

# remap cb_whitelist back
if (remap) {
  cb_whitelist = cb_dict[cb_whitelist] %>% unname
  stopifnot(!duplicated(cb_whitelist))
  rm(cb_dict)
}
rm(remap) ; gc()

count_umis <- function(df) {
  # This method implements the following logic, except faster and more memory efficient:
  # df %<>% group_by(cb_index, sb_index) %>% summarize(umi=n()) 
  
  gdf = df %>% filter(reads>0) %>% select(cb_index,sb_index) %>% arrange(cb_index, sb_index) 
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

# remove reads that didn't match a called cell
df %<>% filter(cb_index>0)

# Compute metrics
Misc(obj, "SB_umi_filtered_downsampling") <- f("metadata/downsampling")
Misc(obj, "SB_umi_final") <- df %>% count_umis %>% pull(umi) %>% sum
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
  pct.placed = round(max(params$pct)*100,2)
  print(g("Optimal eps: {eps}    Optimal minPts: {minPts}    %placed: {pct.placed}"))
  
  return(c(eps,minPts,pct.placed))
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
  coords = lapply(data.list, function(df) {
    p = c(x_um=NA,
          y_um=NA,
          DBSCAN_clusters=max(df$cluster),
          SB_umi = sum(df$origumi),
          SNR=NA,
          SB_bin = unique(df$bin) %>% {ifelse(is.null(.), NA, .)},
          minPts = unique(df$minPts) %>% {ifelse(is.null(.), NA, .)},
          eps = unique(df$eps) %>% {ifelse(is.null(.), NA, .)},
          pct.placed = unique(df$pct.placed) %>% {ifelse(is.null(.), NA, .)})
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p[["x_um"]] = matrixStats::weightedMedian(sdf$x_um,w=sdf$umi)
      p[["y_um"]] = matrixStats::weightedMedian(sdf$y_um,w=sdf$umi)
      p[["SNR"]] = sum(sdf$origumi)/sum(df$origumi)
    }
    return(p)
  }) %>% bind_rows %>% as.data.frame %>% mutate(cb_index=as.numeric(names(data.list))) %>% select(cb_index, everything())
  return(coords)
}

# Run these methods on the entire data
normal_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  params = opt_dbscan(data.list)
  data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
  data.list %<>% map(~mutate(.,bin=quant,eps=params[[1]],minPts=params[[2]],pct.placed=params[[3]]))
  coords <- create_coords(data.list)
  return(list(coords,data.list))
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
    data.list %<>% map(~mutate(.,bin=quant,eps=params[[1]],minPts=params[[2]],pct.placed=params[[3]]))
    return(data.list)
  }) %>% list_flatten
  
  stopifnot(len(data.list)==len(unique(df$cb_index)))
  
  coords <- create_coords(data.list)
  return(list(coords,data.list))
}

correct_beadloss <- function(gdf) {
  gdf %<>% group_by(sb_index) %>% 
    mutate(origumi = umi, totumi=sum(umi), umi=ifelse(totumi>256, umi/totumi*256, umi)) %>% 
    ungroup %>% 
    select(-totumi) %>% 
    arrange(desc(umi))
}

### Positioning ################################################################

# Perform positioning at various levels of downsampling
original_df <- df
coords_list = list()
for (i in seq(0.05,1,0.05)) {
  # Downsample the reads
  print(g("Downsampling: {round(i*100)}%"))
  df = original_df
  if (i != 1) {
    df %<>% mutate(reads=rmultinom(n=1, size=round(sum(original_df$reads)*i), prob=original_df$reads) %>% as.vector)
  }
  df %<>% count_umis
  
  # correct for bead loss (downsample beads above 256 umi)
  df %<>% correct_beadloss
  
  # add spatial positions from puck
  df = merge(x=df,y=puckdf,all.x=T,by="sb_index")
  
  # run binned positioning (see method above)
  res = binned_positioning(df)
  coords = res[[1]] ; data.list = res[[2]] ; rm(res)
  
  coords_list %<>% list.append(coords)
  gc()
}

# merge with seurat object
coords %<>% mutate(cb=cb_whitelist[cb_index])
rownames(coords) = paste0(coords$cb,"-1")
obj = AddMetaData(obj,coords)
Misc(obj,"pct.placed") = round(sum(!is.na(obj$x_um))/ncol(obj)*100,2)
obj$cb_index <- NULL

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
gdraw <- function(text,s=14) {ggdraw()+draw_label(text,size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}
add.commas <- function(num){prettyNum(num,big.mark=",")}
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
Misc(obj,"RNA_metrics") <- list(plotdf[,1],plotdf[,2])

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
    ggtitle("Barcode rank plot")+xlab("Cell barcodes")+ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))#, legend.title=element_text(size=10), legend.text=element_text(size=8), legend.margin=margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"), legend.key.size = unit(0.5, "lines"))
  
  plot = plot_grid(p3,p1,p0,p2,ncol=2)
  
  return(plot)
}

if (file.exists("RNAcounts/molecule_info.h5")) {
  plot=UvsI(obj)
} else {
  plot = gdraw("No molecule_info.h5 found")
}
make.pdf(plot,"plots/1cellcalling.pdf",7,8)
gc()

### Page 2: UMAP + metrics #####################################################

plot = plot_grid(DimPlot(obj,label=T)+ggtitle(g("UMAP"))+NoLegend()+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())+coord_fixed(ratio=1),
                 VlnPlot(obj,"logumi")+NoLegend()+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank()),
                 FeaturePlot(obj,"percent.mt")+ggtitle("%MT")+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())+coord_fixed(ratio=1)+theme(legend.position="top",legend.justification="center",legend.key.width=unit(2, "lines")),
                 FeaturePlot(obj,"pct.intronic")+ggtitle("%Intronic")+theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())+coord_fixed(ratio=1)+theme(legend.position="top",legend.justification="center",legend.key.width=unit(2,"lines")),
                 ncol=2)
make.pdf(plot,"plots/2umap.pdf",7,8)

### Page 3: Raw spatial data ###################################################

# Panel 1+2: elbow plots
p1 = umicurves[[1]]
p2 = umicurves[[2]]

# Panel 3: Downsampling UMI
x = seq(0,1,0.05)*f("metadata/num_reads")/1000000
plot.df = data.frame(x=x,y=f("metadata/downsampling")/1000000)
p3 = ggplot(plot.df, aes(x=x,y=y))+geom_point()+theme_bw()+xlab("Millions of reads")+ylab("Millions of filtered SB UMIs")+ggtitle("Downsampling curve")

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
  return(median(df$dist))
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
  scale_y_continuous(sec.axis = sec_axis(~ . * m2/m1, name = "median displacement (\u00B5m)")) +
  labs(title = "Downsampling Placements", x = "Reads (millions)", y = "Percent placed", color = "") +
  theme(legend.position=c(0.85, 0.85), legend.background=element_blank(), legend.key=element_blank(), legend.key.height=unit(0.75, "lines"))

plot = plot_grid(p1,p2,p3,p4,ncol=2)
make.pdf(plot,"plots/3rawspatial.pdf",7,8)

### Page 4: Beadplot ###########################################################

df = original_df %>% count_umis %>% correct_beadloss %>% merge(y=puckdf,all.x=T,by="sb_index")
sb.data = df %>% group_by(sb_index) %>% dplyr::summarize(umi=sum(origumi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]}
sb.data_weighted = df %>% group_by(sb_index) %>% dplyr::summarize(umi=sum(umi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]}

beadplot <- function(sb.data, m, text){
  sb.data = sb.data[nrow(sb.data):1,]
  ggplot(sb.data, aes(x=x_um,y=y_um,col=umi)) +
    rasterize(geom_point(size=0.1), dpi=200) +
    coord_fixed() +
    theme_classic() +
    labs(x="x (\u00B5m)", y="y (\u00B5m)") +
    scale_color_viridis(trans="log", option="B", name="UMI", limits = c(1, m)) + 
    ggtitle(g("SB UMI per bead ({text})"))
}
p1 = beadplot(sb.data, max(sb.data$umi), "raw")
p2 = beadplot(sb.data_weighted, max(sb.data$umi), "corrected")

plot = plot_grid(p1,p2,ncol=1)
make.pdf(plot,"plots/4beadplot.pdf",7,8)

### Page 5: DBSCAN #############################################################

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
  labs(title = "SNR per cell (density)", x = "SNR", y = "Density") + 
  geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
  annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")

# Panel 3: RNA umi vs SB umi
p3 = data.frame(x=obj$nCount_RNA,y=obj$SB_umi,placed=!is.na(obj$x_um)) %>% {ggplot(.,aes(x=log10(x),y=log10(y),col=placed))+geom_point(size=0.2)+theme_bw()+xlab("RNA UMI")+ylab("SB UMI")+ggtitle("SB UMI vs. RNA UMI")+theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.background = element_blank(), legend.title=element_text(size=10), legend.text=element_text(size=8), legend.margin=margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"), legend.spacing.y = unit(0.1,"lines"), legend.key.size = unit(0.5, "lines"))}

# Panel 4: DBSCAN parameters
df = coords %>% select(SB_bin,minPts,eps,pct.placed) %>% distinct %>% arrange(SB_bin) %>% mutate(SB_bin=round(SB_bin,2),pct.placed=round(pct.placed,2) %>% paste0("%"))
rownames(df) <- NULL
p4 = plot_grid(gdraw("DBSCAN parameters"),plot.tab(df),ncol=1,rel_heights=c(1,17))

plot = plot_grid(p1,p2,p3,p4,ncol=2)
make.pdf(plot,"plots/5DBSCAN.pdf",7,8)

### Page 6: Spatial ############################################################

p1 = DimPlot(obj,reduction="spatial")+coord_fixed(ratio=1)+
  ggtitle(g("%placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)} ({sum(!is.na(obj$x_um))}/{ncol(obj)})")) + 
  NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)")
p2 = DimPlot(obj, reduction="spatial",split.by="seurat_clusters",ncol=7) + theme_void() + coord_fixed(ratio=1) + NoLegend()
plot = plot_grid(p1, p2, ncol=1, rel_heights=c(1,1))
make.pdf(plot,"plots/6spatial.pdf",7,8)

### Page 7: Create metrics plot ################################################

df = list(
  c("Reads",Misc(obj,"SB_reads") %>% add.commas),
  c("Puck file",Misc(obj,"puck")),
  c("Number of beads",Misc(obj,"num_beads") %>% add.commas),
  c("Scaling factor",Misc(obj,"scaling_factor")),
  c("R1<->R2",Misc(obj,"switchR1R2")),
  c("Remap CB",Misc(obj,"remapCB"))
) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("metric","value"))
p1 = plot_grid(ggdraw()+draw_label("SB metrics"), plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

UP_matching <- setNames(Misc(obj,"UP_matching_count"),Misc(obj,"UP_matching_type"))
SB_matching <- setNames(Misc(obj,"SB_matching_count"),Misc(obj,"SB_matching_type"))
CB_matching <- setNames(Misc(obj,"CB_matching_count"),Misc(obj,"CB_matching_type"))

df = data.frame(a=c("exact","fuzzy", "none", "GG"),b=c(UP_matching[["-"]],
                                                       sum(UP_matching[c("1D-","1D-1X","-1X","-1D","-2X")]),
                                                       UP_matching[["none"]],
                                                       UP_matching[["GG"]]) %>% {./sum(.)*100} %>% round(2) %>% paste0("%")) %>% arrange(desc(b)) %>% unname
p2 = plot_grid(ggdraw()+draw_label("UP matching"), plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

df = data.frame(a=c("exact","fuzzy","none","ambig"),b=SB_matching[c("exact","HD1","none","HD1ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
p3 = plot_grid(ggdraw()+draw_label("SB matching"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

df = data.frame(a=c("exact","fuzzy","none","ambig"),b=CB_matching[c("exact","HD1","none","HD1ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
p4 = plot_grid(ggdraw()+draw_label("CB matching"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

df = list(
  c("Valid UMI",UP_matching[["R1lowQ"]]),
  c("Valid UP",UP_matching[c("none","GG")] %>% sum),
  c("Valid SB",SB_matching[c("none","HD1ambig")] %>% sum),
  c("Valid CB",sum(CB_matching[c("none","HD1ambig")])),
  c("Chimeric",Misc(obj,"SB_reads_filtered_chimeric"))
) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("filter","percent"))
df$percent = as.numeric(df$percent) / (Misc(obj,"SB_reads") - c(0,head(cumsum(df$percent),-1)))
df[1:4,2] = 1-df[1:4,2] # convert from fraction removed to fraction retained
df$percent = round(df$percent * 100,2) %>% paste0("%")
p5 = plot_grid(ggdraw()+draw_label("SB filtering"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

Misc(obj, "SB_filtering_type") <- df$filter
Misc(obj, "SB_filtering_count") <- df$percent

# Script input parameters (with common prefix removed)
a=Misc(obj,"RNA_path") ; b=Misc(obj,"SB_path")
m=min(which(!map_lgl(1:min(nchar(a),nchar(b)), ~str_sub(a,1,.)==str_sub(b,1,.))))
a%<>%str_sub(m-1,999) ; b%<>%str_sub(m-1,999)
df=data.frame(a=c("RNA_path","SB_path"),b=c(a,b)) %>% unname
p6 = plot_grid(ggdraw()+draw_label("Run parameters"),plot.tab(df),ncol=1,rel_heights=c(0.1,0.9))

plot = plot_grid(
  gdraw("Additional metadata",18),
  plot_grid(p1,p5,ncol=2),
  plot_grid(p2,p3,p4,ncol=3),
  p6,
  ggdraw()+draw_label(""), #spacer
  ncol=1,
  rel_heights = c(0.27,0.5,0.35,0.25,0.27)
)

make.pdf(plot,"plots/7metrics.pdf",7,8)

### Sample bead plots ##########################################################

cols = list(c("steelblue4","steelblue1"),c("red4","red1"),c("springgreen4","springgreen1"))
plot.sb <- function(subdf) {
  limits = range(subdf$origumi)
  nlegend = min(3,max(subdf$cluster))
  pts = lapply(1:nlegend,function(i){
    if (i <= 2) {
      list(geom_point(data=filter(subdf,cluster==i),mapping=aes(x=x_um,y=y_um,col=origumi),size=0.3), scale_color_gradient(low=cols[[i]][[1]],high=cols[[i]][[2]],limits=limits,breaks=limits), new_scale_color())
    } else {
      list(geom_point(data=filter(subdf,cluster>=3),mapping=aes(x=x_um,y=y_um,col=origumi),size=0.3), scale_color_gradient(low=cols[[i]][[1]],high=cols[[i]][[2]],limits=limits,breaks=limits), new_scale_color())
    }
  })
  Reduce(`+`,pts,init=ggplot())+coord_fixed(ratio=1,xlim=range(puckdf$x_um),ylim=range(puckdf$y_um))+theme_void()+
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

pdfs = c("0cellranger.pdf","1cellcalling.pdf", "2umap.pdf", "3rawspatial.pdf", "4beadplot.pdf", "5DBSCAN.pdf","6spatial.pdf","7metrics.pdf","SB.pdf") %>% paste0("plots/",.)
qpdf::pdf_combine(input = pdfs, output = "summary.pdf")
qsave(obj, "seurat.qs") # we added some more metadata while plotting

print("Done!")
