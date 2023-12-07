library(glue) ; g=glue ; len=length
library(matrixStats)
library(ggnewscale)
library(stringdist)
library(gridExtra)
library(parallel)
library(magrittr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggrastr)
library(Seurat)
library(dbscan)
library(dplyr)
library(purrr)
library(rhdf5)
library(qpdf)
library(qs)

setwd("~/SBcounts/317/")

### Download files #############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript spatial.R RNApath SBpath", call. = FALSE)
}

RNApath <- args[1] ; print(g("RNApath: {RNApath}"))
SBpath <- args[2] ; print(g("SBpath: {SBpath}"))

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/230317_SL-NSA_0564_AHVCHWDRX2/SI-TT-E11"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/230317_SL-NSA_0564_AHVCHWDRX2/SI-NT-C4"

# RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/231020_VL00181_99_AAF2KJGM5/SI-TS-A11"
# SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/231020_VL00181_99_AAF2KJGM5/SI-TT-A1"

stopifnot(!file.exists("RNAcounts"))
stopifnot(!file.exists("SBcounts"))

checkgsfile <- function(path) {return(system(g("gsutil ls {path}"),intern=F,ignore.stdout=T,ignore.stderr=T)==0)}

system("mkdir RNAcounts")
if (checkgsfile(file.path(RNApath,"outs/filtered_feature_bc_matrix.h5"))) {
  RNAtech = "count"
  system(g("gsutil cp {RNApath}/outs/filtered_feature_bc_matrix.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/raw_feature_bc_matrix.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/molecule_info.h5 RNAcounts"))
} else if (checkgsfile(file.path(RNApath,"outs/multi"))) {
  RNAtech = "multi"
  system(g("gsutil cp {RNApath}/outs/multi/count/raw_feature_bc_matrix.h5 RNAcounts"))
  system(g("gsutil cp {RNApath}/outs/per_sample_outs/{basename(RNApath)}/count/sample_filtered_feature_bc_matrix.h5 RNAcounts/filtered_feature_bc_matrix.h5"))
  system(g("gsutil cp {RNApath}/outs/per_sample_outs/{basename(RNApath)}/count/sample_molecule_info.h5 RNAcounts/molecule_info.h5"))
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
  info = data.frame(barcode=fetch("barcodes")[fetch("barcode_idx")+1] %>% paste0("-1"),
                             feature=fetch("features/name")[fetch("feature_idx")+1],
                             umi=fetch("umi"),
                             umi_type=fetch("umi_type"))
  info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
  obj$pct.intronic = info$pct.intronic[match(colnames(obj),info$barcode)]
  rm(info) ; rm(fetch)
}

cb_whitelist = unname(obj$cb)

### Load the SB counts matrix ##################################################

# Accessor method for the spatial .h5
f <- function(p){return(h5read("SBcounts/SBcounts.h5",p))}

cb_list = f("lists/cb_list")
cb_list_remap = f("lists/cb_list_remap")
sb_list = f("lists/sb_list")

df = data.frame(cb_index=f("matrix/cb_index"),
                umi_2bit=f("matrix/umi"),
                sb_index=f("matrix/sb_index"),
                reads=f("matrix/reads"))

# Remap the whitelist (if necessary)
reads_noremap = df %>% filter(cb_list[cb_index] %in% cb_whitelist) %>% pull(reads) %>% sum
reads_remap = df %>% filter(cb_list_remap[cb_index] %in% cb_whitelist) %>% pull(reads) %>% sum
remap = reads_remap > reads_noremap ; rm(reads_remap) ; rm(reads_noremap)
if (remap) {
  cb_dict <- read.table(CBdictpath) %>% {setNames(.[[2]],.[[1]])}
  print("Remapping CB whitelist")
  cb_whitelist = cb_dict[cb_whitelist]
}
stopifnot(!duplicated(cb_whitelist))

reads_after_cb_filter_nofuzzy = df %>% mutate(cb=cb_list[cb_index]) %>% filter(cb %in% cb_whitelist) %>% pull(reads) %>% sum

# Perform the HD1 matching between observed cell barcodes and the whitelist
chunk_size = round(50000000/len(cb_whitelist))
chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}
res = lapply(chunk_vector(cb_list,chunk_size), function(chunk){
  dists=stringdist::stringdistmatrix(chunk, cb_whitelist, method = "hamming")
  row0s = rowSums(dists == 0)
  row1s = rowSums(dists == 1)
  m = (row0s==1) | (row1s==1)
  mins = apply(dists[m,],1,function(row){which.min(row)})
  return(list(which(m),mins))
})
matching_dict = setNames( # keys are cb_list_index, values are cb_whitelist_index
  map(res,pluck(2)) %>% flatten_int,
  map(res,pluck(1)) %>% map2(1:length(res),~.x+chunk_size*(.y-1)) %>% flatten_int
)
stopifnot(!duplicated(names(matching_dict)))
rm(res) ; rm(chunk_size)

# Convert cb_index from an index into cb_list to an index into cb_whitelist
df %<>% mutate(cb_index = matching_dict[as.character(cb_index)])
rm(matching_dict) ; rm(cb_list) ; rm(cb_list_remap)

# remove reads that didn't match a called cell
reads_before_cb_filter = sum(df$reads)
df %<>% filter(!is.na(cb_index))
reads_after_cb_filter = sum(df$reads)

# compress duplicate rows introduced by fuzzy matching
df %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads = sum(reads)) %>% ungroup

# remap the whitelist back
if (remap) {
  cb_whitelist = cb_dict[cb_whitelist] %>% unname
  rm(cb_dict)
}
stopifnot(!duplicated(cb_whitelist))

# Remove chimeric reads
df %<>% group_by(cb_index,umi_2bit) %>% mutate(ismax = reads==max(reads), toptie = sum(ismax)>1) %>% ungroup
reads_before_chimera_filter = sum(df$reads)
df %<>% filter(ismax==T & toptie == F) %>% select(-ismax,-toptie)
reads_after_chimera_filter = sum(df$reads)

# load the puck information
pucks = f("puck/puck_list") ; stopifnot(len(pucks)==1)
puckdf = data.frame(sb=f("puck/sb"), x=f("puck/x"), y=f("puck/y"), puck_index=f("puck/puck_index"))
stopifnot(unique(puckdf$puck_index)==1)
puckdf = puckdf[match(sb_list,puckdf$sb),]
stopifnot(puckdf$sb == sb_list) ; rm(sb_list)
bn = nrow(puckdf)
if (bn < 150000) {
  k = 0.73
} else if (bn < 600000) {
  k = 0.73 * 2
} else {
  k = 0.645
}
puckdf %<>% transmute(sb_index = 1:nrow(puckdf), x_um = x*k, y_um = y*k)

gc()

### Positioning methods ########################################################

chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}

# Do a grid search to find the ideal DBSCAN parameters
opt_dbscan <- function(data.list) {
  eps.vec = c(50) ; minPts.vec = c(3:500)
  params = expand.grid(eps.vec,minPts.vec) %>% setNames(c("eps","minPts"))
  row_lists = chunk_vector(1:nrow(params), round(nrow(params)/30))
  params$pct = parallel::mclapply(row_lists, function(v) {
    map_dbl(v, function(i) {
      m = map_lgl(data.list, ~max(dbscan::dbscan(.[c("x_um","y_um")], eps=params$eps[[i]], minPts=params$minPts[[i]], weights=.$umi)$cluster) == 1)
      return(sum(m)/length(m))
    })
  },mc.cores=30L) %>% flatten_dbl
  params$is.max = params$pct==max(params$pct)
  
  eps = params$eps[params$is.max][[1]] ; minPts = params$minPts[params$is.max][[1]]
  print(g("Optimal eps: {eps}    Optimal minPts: {minPts}    %placed: {round(max(params$pct)*100,2)}"))
  
  return(c(eps,minPts,max(params$pct)))
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
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p = c(x_um = matrixStats::weightedMedian(sdf$x_um,w=sdf$umi),
            y_um = matrixStats::weightedMedian(sdf$y_um,w=sdf$umi),
            DBSCAN_clusters = 1,
            num_SBumi = sum(df$origumi),
            SNR = sum(sdf$origumi)/sum(df$origumi),
            SB_bin=unique(df$bin)#,
            # a = max(df$z1),
            # b = max(df$z2),
            # c = max(sdf$z1),
            # d = max(sdf$z2)
            )
    } else {
      p = c(x_um=NA,y_um=NA,DBSCAN_clusters=max(df$cluster),num_SBumi=sum(df$origumi),SNR=NA,SB_bin=unique(df$bin))#,z1=NA,z2=NA)
      #p = NA
    }
  }) %>% {do.call(rbind,.)} %>% as.data.frame %>% mutate(cb_index=names(data.list)) %>% select(cb_index, everything())
}

# Run these methods on the entire data
normal_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  params = opt_dbscan(data.list)
  data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
  coords <- create_coords(data.list)
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
    params = opt_dbscan(data.list)
    data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
    data.list %<>% map(~mutate(.,bin=quant))
    return(data.list)
  }) %>% list_flatten
  coords <- create_coords(data.list)
  return(coords)
}

### Positioning ################################################################

# old_df = df
# old_obj = obj

# Perform positioning at various levels of downsampling
nreads = sum(df$reads)
reads_cumsum = map(1:nrow(df),~rep(.,df$reads[[.]])) %>% flatten_int
coords_list = lapply(seq(0.05,1,0.05), function(i) {
  print(g("Downsampling: {round(i*100)}%"))
  
  # downsample the dataframe
  newreads = sample(reads_cumsum,round(nreads*i),replace=F) %>% table %>% as.data.frame %>% setNames(c("row","newreads"))
  df %<>% mutate(row=1:nrow(df)) %>% merge(newreads,by="row") %>% mutate(reads=newreads) %>% select(-newreads,-row)
  
  # group by cb,sb to count umis
  df %<>% group_by(cb_index,sb_index) %>% summarize(umi=n()) %>% ungroup %>% arrange(desc(umi))
  
  # correct for bead loss (downsample beads above 256 umi)
  maxumi = 256
  df %<>% group_by(sb_index) %>% mutate(origumi = umi, totumi=sum(umi), umi=ifelse(totumi>maxumi, umi/totumi*maxumi, umi)) %>% ungroup %>% select(-totumi) %>% arrange(desc(umi))
  
  # add spatial positions from puck
  df = merge(x=df,y=puckdf,all.x=T,by="sb_index")
  
  # run positioning
  coords <- binned_positioning(df)
  return(coords)
})

# merge with seurat object
coords <- tail(coords_list,1)[[1]]
coords %<>% mutate(cb=cb_whitelist[as.numeric(cb_index)])
rownames(coords) = paste0(coords$cb,"-1")
obj = AddMetaData(obj,coords)

emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_")

coords %>% select(cb,everything()) %>% select(-cb_index) %>% write.csv("coords.csv",quote=F,row.names=F)
qsave(obj, "seurat.qs")

df %<>% group_by(cb_index,sb_index) %>% summarize(umi=n()) %>% ungroup %>% arrange(desc(umi))
maxumi = 256
df %<>% group_by(sb_index) %>% mutate(origumi = umi, totumi=sum(umi), umi=ifelse(totumi>maxumi, umi/totumi*maxumi, umi)) %>% ungroup %>% select(-totumi) %>% arrange(desc(umi))
df = merge(x=df,y=puckdf,all.x=T,by="sb_index")
################################################################################
### CREATE PDF #################################################################
################################################################################

system("mkdir plots")

plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}

make.pdf <- function(plots,name,w,h) {
  if (any(class(plots)=="gg")||class(plots)=="Heatmap") {plots=list(plots)}
  pdf(file=name,width=w,height=h)
  lapply(plots,function(x){print(x)})
  dev.off()
}

### Page 1: cell calling #######################################################

UvsI <- function() {
  fetch <- function(x){return(h5read(g("RNAcounts/molecule_info.h5"),x))}
  molecule_info = data.frame(barcode=fetch("barcodes")[fetch("barcode_idx")+1] %>% paste0("-1"),
                             feature=fetch("features/name")[fetch("feature_idx")+1],
                             umi=fetch("umi"),
                             umi_type=fetch("umi_type"),
                             count=fetch("count"))
  df = molecule_info %>% group_by(barcode) %>% summarize(umi=n(), pct.intronic=sum(umi_type==0)/umi) %>% arrange(desc(umi)) %>% mutate(logumi=log10(umi))
  
  p1 = df %>% filter(umi>50) %>% ggplot(aes(x = logumi, y = pct.intronic)) + 
    geom_bin2d(bins=100) +
    scale_fill_viridis(trans="log", option="A", name="density") + 
    theme_minimal() +
    labs(title = "Intronic vs. UMI droplets (>50 umi)", x = "logumi", y = "%intronic")
  
  max_density_x = density(filter(df,umi>500,pct.intronic>0.3)$pct.intronic) %>% {.$x[which.max(.$y)]}
  p2 = df %>% filter(umi>500) %>% ggplot(aes(x = pct.intronic)) +
    geom_density() + 
    theme_minimal() +
    labs(title = "Intronic density (>500 umi)", x = "%intronic", y = "Density") + 
    geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
    annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
  
  # p3 = data.frame(metric=c("#cells","mean logumi", "mean %MT"),value=c(ncol(obj) %>% round %>% as.integer,
  #                                                          mean(obj$logumi) %>% round(2),
  #                                                          mean(obj$percent.mt) %>% round(2))) %>% plot.tab

  plot = plot_grid(p1,p2,ncol=1)
  
  return(plot)
}

plot=UvsI()

make.pdf(plot,"plots/1cellcalling.pdf",7,8)


### Page 2: UMAP + metrics #####################################################

plot = plot_grid(DimPlot(obj,label=T)+ggtitle(g("UMAP"))+NoLegend(),
                 VlnPlot(obj,"logumi")+NoLegend(),
                 FeaturePlot(obj,"percent.mt"),
                 FeaturePlot(obj,"pct.intronic"),
                 ncol=2)

make.pdf(plot,"plots/2umap.pdf",7,8)

### Page 3: FeaturePlots #####################################################

plot = FeaturePlot(obj,c("MBP","PLP1","CSPG5","VCAN","AQP4","GFAP","FLT1","DCN","CX3CR1","GAD1", "GAD2", "RELN","SLC17A6","SLC17A7", "TH", "NR4A2"))&theme_void()+
  theme(legend.key.width = unit(0.3, "cm"), plot.title=element_text(hjust=0.5)) & coord_fixed(ratio=1)
make.pdf(plot,"plots/3features.pdf",7,8)

### Page 4: Raw spatial data ###################################################

cb.data = df %>% group_by(cb_index) %>% dplyr::summarize(unique.sb=n(),umi=sum(origumi)) %>% {.[order(.$umi,decreasing=T),]} %>% {mutate(.,index=1:nrow(.))}
sb.data = df %>% group_by(sb_index) %>% dplyr::summarize(unique.cb=n(),umi=sum(origumi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]} %>% {mutate(.,index=1:nrow(.))}
sb.data_weighted = df %>% group_by(sb_index) %>% dplyr::summarize(unique.cb=n(),umi=sum(umi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]} %>% {mutate(.,index=1:nrow(.))}

p1 = ggplot(sb.data,aes(x=log10(index),y=log10(umi)))+rasterize(geom_point(size=0.5), dpi=200)+theme_bw()+ggtitle("SB UMI per bead")
p2 = ggplot(cb.data,aes(x=log10(index),y=log10(umi)))+rasterize(geom_point(size=0.5), dpi=200)+theme_bw()+ggtitle("SB UMI per cell (called cells only)")

x = seq(0,1,0.05)*reads/1000000

plot.df = data.frame(x=x,y=f("metadata/downsampling")/1000000)
p3 = ggplot(plot.df, aes(x=x,y=y))+geom_point()+theme_bw()+xlab("millions of reads")+ylab("millions of SB-whitelist UMIs")+ggtitle("Downsampling curve")

reads = f("metadata/num_reads")
plot.df = map(coords_list,function(df){return(df$DBSCAN_clusters %>% {c(sum(.==0),sum(.==1),sum(.==2))/nrow(df)*100})}) %>% {do.call(rbind,.)} %>% {rbind(c(0,0,0),.)}
plot.df %<>% as.data.frame %>% setNames(c("0","1","2")) %>% mutate(x=x)
plot.df = tidyr::gather(plot.df, key = "column", value = "value", -x)
p4 = ggplot(plot.df, aes(x = x, y = value, color = column)) + geom_line() +
  labs(title = "Downsampling Placements", x = "millions of reads", y = "Percent", color = "Clusters") +
  theme_minimal()

plot = plot_grid(p1,p2,p3,p4,ncol=2)
make.pdf(plot,"plots/4rawspatial.pdf",7,8)

### Page 5: Beadplot ###########################################################

beadplot <- function(sb.data, m, text){
  sb.data = sb.data[nrow(sb.data):1,]
  ggplot(sb.data, aes(x=x_um,y=y_um,col=umi)) +
    rasterize(geom_point(size=0.5), dpi=200) +
    coord_fixed() +
    theme_classic() +
    labs(x="x (um)", y="y (um)") +
    scale_color_viridis(trans="log", option="B", name="UMI", limits = c(1, m)) + 
    ggtitle(g("SB UMI per bead ({text})"))
}
p1 = beadplot(sb.data, max(sb.data$umi), "raw")
p2 = beadplot(sb.data_weighted, max(sb.data$umi), "corrected")
#p = plot_grid(p1,ggdraw()+draw_label(g("Max UMI per bead: {maxumi}")),p2,ncol=1,rel_heights=c(1,0.1,1))
plot = plot_grid(p1,p2,ncol=1)
make.pdf(plot,"plots/5beadplot.pdf",7,8)

### Page 6: DBSCAN #############################################################

d=data.frame(x=coords$DBSCAN_clusters)
p1 = ggplot(d,aes(x=x))+geom_histogram(aes(y = after_stat(count)/sum(after_stat(count)) * 100))+
  geom_text(aes(label = sprintf("%1.0f%%", after_stat(count)/sum(after_stat(count))*100), y=after_stat(count)/sum(after_stat(count))*100), stat="bin", binwidth=1, vjust=-0.5)+
  theme_classic()+xlab("Num DBSCAN clusters")+ylab("Percent")+scale_y_continuous(limits=c(0,100))+scale_x_continuous(breaks=min(d$x):max(d$x))+ggtitle("DBSCAN cluster distribution")

p2 = data.frame(x=obj$nCount_RNA,y=obj$num_SBumi,placed=!is.na(obj$x_um)) %>% {ggplot(.,aes(x=log10(x),y=log10(y),col=placed))+geom_point(size=0.2)+theme_bw()+xlab("RNA UMI")+ylab("SB UMI")+ggtitle("SB UMI vs. RNA UMI")+theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.background = element_blank(), legend.title=element_text(size=10), legend.text=element_text(size=8), legend.margin=margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"), legend.spacing.y = unit(0.1,"lines"), legend.key.size = unit(0.5, "lines"))}

# res = map(data.lists,opt_dbscan) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("eps","minPts","%placement")) %>% mutate(quantile=round(quants[-1])) %>% select(c("SB UMI quantile","eps","minPts","pct"))
# res %<>% mutate(pct = round(pct*100,2))
# p2 = plot.tab(res)

plot = plot_grid(p1,p2,ncol=1)
make.pdf(plot,"plots/6DBSCAN.pdf",7,8)

### Page 7: Spatial ############################################################

plot = plot_grid(DimPlot(obj,reduction="spatial")+coord_fixed(ratio=1)+ggtitle(g("{ncol(obj)} cells                %placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)}")),
               DimPlot(obj, reduction="spatial",split.by="seurat_clusters",ncol=6)&theme_void()&coord_fixed(ratio=1)&NoLegend(),
               ncol=1, rel_heights=c(1,1))

make.pdf(plot,"plots/7spatial.pdf",7,8)

### Create metrics plot ########################################################

Rs = data.frame(R1=f("metadata/R1s") %>% map_chr(basename),
                R2=f("metadata/R2s") %>% map_chr(basename))
Rs = rbind(Rs,map(1:(4-nrow(Rs)),~c(R1="",R2="")) %>% {do.call(rbind,.)} %>% as.data.frame)

UPmatching = setNames(f("metadata/UP_matching/count"),f("metadata/UP_matching/type"))
UP_matching = data.frame(type=f("metadata/UP_matching/type"),
                         count=f("metadata/UP_matching/count")) %>% arrange(desc(count)) %>% mutate(pct=round(count/sum(count)*100,2)) %>% mutate(type=gsub("^-$","exact",type))

SBmatching = setNames(f("metadata/SB_matching/count"),f("metadata/SB_matching/type"))
SB_matching = data.frame(type=f("metadata/SB_matching/type"),
                         count=f("metadata/SB_matching/count")) %>% arrange(desc(count)) %>% mutate(pct=round(count/sum(count)*100,2))

# chimerastats = setNames(f("metadata/chimera_stats/count"),f("metadata/chimera_stats/type"))
# chimera_stats = data.frame(type=c("reads_removed","umis_removed"),
#                            count=c(chimerastats[["reads_before"]]-chimerastats[["reads_after"]],chimerastats[["umis_before"]]-chimerastats[["umis_after"]]))
# chimera_stats$pct = round(chimera_stats$count/c(chimerastats[["reads_before"]],chimerastats[["umis_before"]])*100,2)

CB_matching = data.frame(type=c("exact","HD1"),count=c(reads_before_cb_filter-reads_after_cb_filter_nofuzzy,reads_after_cb_filter-reads_after_cb_filter_nofuzzy)) %>% mutate(pct=round(count/sum(count)*100,2))

reads = f("metadata/num_reads")
metadata = c(SB_reads=reads,
             switchR1R2=f("metadata/switch") %>% as.logical,
             CBremap=remap,
             puck=pucks,
             num_beads=bn,
             scaling_factor=k
             ) %>% {data.frame(name=names(.),value=.)} ; rownames(metadata) <- NULL

filterdf = c(R1lowQ=UPmatching[["R1lowQ"]],
             noUP=UPmatching[["none"]]+UPmatching[["GG"]],
             noSB=SBmatching[["none"]]+SBmatching[["HD1ambig"]],
             chimeric = reads_before_chimera_filter-reads_after_chimera_filter,
             uncalledCB = reads_before_cb_filter-reads_after_cb_filter) %>% {data.frame(filter=names(.),removed=unname(.))}
filterdf %<>% mutate(pct.total.removed = round(removed/reads*100,2))
filterdf$pct.removed = round(filterdf$removed/(c(reads,reads-cumsum(filterdf$removed))%>%{.[1:(len(.)-1)]})*100,2)
filterdf %<>% transmute(filter=filter,pct=pct.removed,pct.total=pct.total.removed)
filterdf = rbind(filterdf, c("total",NA,sum(filterdf$pct.total)))

# plot.tab(Rs)
# plot.tab(UP_matching)
# plot.tab(SB_matching)
# plot.tab(chimera_stats)
# plot.tab(metadata)
# plot.tab(filterdf)

metricplot = plot_grid(
          plot_grid(plot_grid(ggdraw()+draw_label("Metadata"),plot.tab(metadata),ncol=1,rel_heights=c(0.1,1)),
                    plot_grid(ggdraw()+draw_label("Reads statistics"),plot.tab(filterdf),ncol=1,rel_heights=c(0.1,1)),
                    rel_widths=c(0.5,0.5), ncol=2),
          plot_grid(plot_grid(ggdraw()+draw_label("UP matching"),plot.tab(UP_matching),ncol=1,rel_heights=c(0.1,0.7)),
                    plot_grid(ggdraw()+draw_label("SB matching"),plot.tab(SB_matching),ggdraw()+draw_label("CB matching"),plot.tab(CB_matching),ncol=1,rel_heights=c(1,5,1,3)),
                    rel_widths=c(0.5,0.5), ncol=2),
          ggdraw()+draw_label("FASTQ paths"), plot.tab(Rs),
          rel_heights=c(3.5,4.5,0.4,2), ncol=1)

make.pdf(metricplot,"plots/8metrics.pdf",7,8)

### Sample bead plots ##########################################################
# 
# cols = list(c("steelblue4","steelblue1"),c("red4","red1"),c("springgreen4","springgreen1"))
# plot.sb <- function(subdf) {
#   nlegend = min(4,max(subdf$cluster))
#   pts = lapply(1:nlegend,function(i){
#     if (i <= 3) {
#       list(geom_point(data=filter(subdf,cluster==i),mapping=aes(x=x_um,y=y_um,col=umi),size=0.5), scale_color_gradient(low=cols[[i]][[1]],high=cols[[i]][[2]]), new_scale_color())
#     } else {
#       list(geom_point(data=filter(subdf,cluster>=4),mapping=aes(x=x_um,y=y_um,col=umi),size=0.5), scale_color_gradient(low="grey20",high="grey80"), new_scale_color())
#     }
#   })
#   Reduce(`+`,pts,init=ggplot())+coord_fixed(ratio=1,xlim=range(df$x_um),ylim=range(df$y_um))+theme_void()+
#     theme(legend.key.width=unit(0.5,"lines"), legend.position="right", legend.key.height=unit(1/nlegend*1.5,"lines"), legend.title=element_blank(), legend.spacing.y=unit(0.2,"lines"), legend.margin=margin(0,0,0,0,"lines"), legend.box.margin=margin(0,0,0,0,"pt"), legend.box.background=element_blank(), legend.background=element_blank(), legend.direction="vertical")
# }
# 
# list0 = data.list %>% keep(~max(.$cluster)==0) %>% map(~mutate(.,cluster=cluster+1) %>% arrange(umi))
# list1 = data.list %>% keep(~max(.$cluster)==1) %>% map(~mutate(.,cluster=cluster+1) %>% arrange(umi))
# list2 = data.list %>% keep(~max(.$cluster)>=2) %>% map(~mutate(.,cluster=cluster+1) %>% arrange(umi))
# 
# p1 = plot_grid(ggdraw()+draw_label("DBSCAN=0"), map(sample(list0,16,replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=4)}, ncol=1,rel_heights=c(0.1,2))
# p2 = plot_grid(ggdraw()+draw_label("DBSCAN=1"), map(sample(list1,16,replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=4)}, ncol=1,rel_heights=c(0.1,2))
# p3 = plot_grid(ggdraw()+draw_label("DBSCAN>=2"), map(sample(list2,16,replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=4)}, ncol=1,rel_heights=c(0.1,2))
# 
# make.pdf(list(p1,p2,p3),"plots/SB.pdf",10,10)
# 
# legend.key.size,
# legend.text,
# legend.text.align,
# legend.justification,
# legend.box,
# legend.box.just,
# legend.box.spacing
# 
# obj$numbead = unname(map_int(list1,~sum(.$cluster==2)))[match(obj$cb, obj$cb[unname(map_int(list1,~unique(.$cb_index)))] %>% unname)]
# obj$numbeadlog = log10(obj$numbead)
# obj$bnb=obj$numbead<2
# 
# pt = map(unique(obj$seurat_clusters) %>% sort, ~DimPlot(obj[,obj$seurat_clusters==.],group.by="bnb",reduction="spatial")+theme_void()) %>% {plot_grid(plotlist=.,ncol=7)}
# pt
# 
# pt = map(unique(obj$seurat_clusters) %>% sort, ~FeaturePlot(obj[,obj$seurat_clusters==.],features="numbeadlog",reduction="spatial")+theme_void()+ggtitle(.)) %>% {plot_grid(plotlist=.,ncol=7)}
# pt
# 
# obj[,obj$seurat_clusters==2] %>% DimPlot(group.by="bnb",reduction="spatial")

### Save output ################################################################

pdfs = c("1cellcalling.pdf", "2umap.pdf", "3features.pdf", "4rawspatial.pdf", "5beadplot.pdf", "6DBSCAN.pdf","7spatial.pdf","8metrics.pdf") %>% paste0("plots/",.)
qpdf::pdf_combine(input = pdfs, output = "summary.pdf")

