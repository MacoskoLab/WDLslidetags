library(glue) ; g=glue ; len=length
library(gridExtra)
library(jsonlite)
library(parallel)
library(magrittr)
library(ggplot2)
library(cowplot)
library(viridis)
library(stringr)
library(Seurat)
library(rlist)
library(dplyr)
library(purrr)
library(qpdf)
library(qs)
# library(matrixStats)
# library(dbscan)

# This script requires 3 things:
#    INDEX/ with coords.csv metrics.json params.csv
#    INTRON/
#    NOINTRON/

INDEX = list.dirs(full.names=F)[map_lgl(list.dirs(),~all(c("coords.csv","metrics.json","params.csv") %in% list.files(.)))]
stopifnot(len(INDEX)==1)

################################################################################

plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}

add.meta <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  obj[["logumi"]] <- log10(obj$nCount_RNA+1)
  obj[["cb"]] <- map_chr(colnames(obj), ~sub("-1$", "", .))
  obj[["type"]] <- "unknown"
  return(obj)
}

process <- function(obj,res=0.8, n.epochs=NULL) {
  obj <- obj %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose=F) %>%
    Seurat::FindNeighbors(dims=1:30) %>%
    Seurat::FindClusters(resolution=res) %>%
    Seurat::RunUMAP(dims=1:30, verbose=F, n.epochs=n.epochs)
}

searchplot <- function(mat) {
  ggplot(mat %>% as.data.frame, aes(x=eps,y=minPts,fill=pct)) + geom_tile() +
    geom_text(aes(label=round(pct*100,1),color=is.max)) + theme(legend.position="none") +
    scale_x_continuous(breaks = mat$eps %>% unique %>% sort) + 
    scale_y_continuous(breaks = mat$minPts %>% unique %>% sort) +
    scale_color_manual(values = c("black", "red"))
}

DBSCAN <- function(data.list, eps, minPts) {
  coords = lapply(data.list, function(df) {
    df$cluster <- dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$umi)$cluster
    m = max(df$cluster) ; umisum = sum(df$umi) ; umimax = max(df$umi)
    if (m == 1) {
      df %<>% dplyr::filter(cluster==1)
      # xmu = weighted.mean(df$x_um,w=df$umi) ; ymu = weighted.mean(df$y_um,w=df$umi)
      xmu = matrixStats::weightedMedian(df$x_um,w=df$umi) ; ymu = matrixStats::weightedMedian(df$y_um,w=df$umi)
      # df %<>% mutate(wd=umi*((xmu-x_um)^2+(ymu-y_um)^2)) ; SE = sqrt(sum(df$wd)/sum(df$umi)/sum(df$umi))
      return(c(m, xmu, ymu, sum(df$umi),umisum,umimax))
    } else {
      return(c(m, NA, NA, NA, umisum,umimax)) 
    }
  }) %>% do.call(rbind,.) %>% as.data.frame %>% setNames(c("clusters","x_um","y_um","clustumi","umisum","umimax"))
  return(coords)
}

RUN.DBSCAN <- function(data, normalize=F, rpm=F) {
  data.list = split(data, data$cb)
  if (normalize) {data.list%<>%resample.normalize}
  
  # Do a grid search
  eps.vec = c(50)
  minPts.vec = c(3:100)
  params = expand.grid(eps.vec,minPts.vec) ; colnames(params) = c("eps","minPts")
  mat = mclapply(1:nrow(params), function(i) {
    eps = params$eps[[i]] ; minPts = params$minPts[[i]]
    m = data.list %>% map_lgl(~max(dbscan::dbscan(.[c("x_um","y_um")], eps=eps, minPts=minPts, weights=.$umi)$cluster) == 1)
    return(c(eps, minPts, sum(m)/length(m)))
  },mc.cores=30L) %>% {do.call(rbind,.)} %>% as.data.frame ; colnames(mat) = c("eps","minPts","pct")
  mat$is.max = mat$pct==max(mat$pct)
  
  if(rpm){return(mat$pct[mat$is.max][[1]])}
  
  # Plot the grid
  print(searchplot(mat))
  ggsave("searchplot.pdf",width=20,height=20)
  
  # Rerun under the optimal parameter set
  eps = mat$eps[mat$is.max][[1]] ; minPts = mat$minPts[mat$is.max][[1]]
  print(g("Optimal eps: {eps}    Optimal minPts: {minPts}"))
  coords = DBSCAN(data.list,eps,minPts)
  return(coords)
}

################################################################################
#### REQUIRED: MAKE SEURAT OBJECT ##############################################
################################################################################

obj <- "INTRON" %>% list.files(full.names=T) %>% Read10X_h5 %>% CreateSeuratObject

mat1 = "INTRON" %>% list.files(full.names=T) %>% Read10X_h5
mat2 = "NOINTRON" %>% list.files(full.names=T) %>% Read10X_h5
mat2 = mat2[,colnames(mat2) %in% colnames(mat1)]
stopifnot(colnames(mat1) == colnames(mat2))
stopifnot(rownames(mat1) == rownames(mat2))
diff = mat1 - mat2
stopifnot(colnames(obj) == colnames(diff))

obj[["NOINTRON"]] <- CreateAssayObject(counts=mat2)
DefaultAssay(obj) = "RNA"
obj %<>% add.meta %>% process
obj$percent.intronic = colSums(diff) / colSums(mat1) * 100

# run DBSCAN
data <- read.csv(file.path(INDEX, "coords.csv"))
bn = list.files("PUCK",full.names=T) %>% read.csv(header=F) %>% nrow
if (bn < 150000) {
  k = 0.73
} else if (bn < 600000) {
  k = 0.73 * 2
} else {
  k = 0.645
}
print(g("Scaling factor: {k}"))
data$x_um = data$x * k
data$y_um = data$y * k
coords <- RUN.DBSCAN(data)

m = match(obj$cb,rownames(coords))
obj$x_um = coords$x_um[m]
obj$y_um = coords$y_um[m]
obj$clusters = coords$clusters[m]
# obj@meta.data[,c("x_um","y_um","cb")] <- coords[match(rownames(obj@meta.data), row.names(coords)),c("x_um","y_um","cb")]

emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_")

qsave(obj,g("{INDEX}.qs"))

################################################################################
#### OPTIONAL: MAKE SUMMARY PDF ################################################
################################################################################

# Print run parameters
params = read.csv(file.path(INDEX,"params.csv"))
p1 = plot.tab(params[1:(nrow(params)-1),] %>% setNames(c("param","value")))
fastqs = params[nrow(params),2] %>% fromJSON
p21 = as.data.frame(basename(fastqs[grep("_R1_",fastqs)])) %>% setNames("R1") %>% plot.tab
p22 = as.data.frame(basename(fastqs[grep("_R2_",fastqs)])) %>% setNames("R2") %>% plot.tab
plot_grid(p1,plot_grid(p21,p22,nrow=1),ncol=1,rel_heights=c(11,len(fastqs)/2+1))
ggsave("params.pdf",width=15,height=5)
# params[4,2] %<>% strsplit("/") %>% unlist %>% tail(4) %>% paste(collapse="/")
# params[5,2] %<>% strsplit("/") %>% unlist %>% tail(4) %>% paste(collapse="/")
# params[6,2] %<>% basename
# p0 = params %>% {.[c(1,2,3,4,5,6,10),]} %>% plot.tab ; p0
# p11 = params %>% {.[c(7,8,9,11),]} %>% plot.tab ; p11
# p12 = eval(params[12,2]) %>% fromJSON %>% data.frame(fastqs=.) %>% plot.tab ; p12
# p0,plot_grid(p11,p12,ncol=2),

# Print positioning metrics
data_list <- fromJSON(file.path(INDEX,"metrics.json"))
array_3D <- array(0, dim = c(5, 4, 5))
dimnames(array_3D) <- list(dim1=c("none","uncalled","fuzzy.many","fuzzy","exact"), dim2=c("none","allG","fuzzy","exact"), dim3=c("noUP","none","fuzzy.many","fuzzy","exact"))
for (i in dimnames(array_3D)$dim1) { for (j in dimnames(array_3D)$dim2) { for (k in dimnames(array_3D)$dim3) { array_3D[i,j,k]<-data_list[[i]][[j]][[k]] }}}
dimnames(array_3D) <- lapply(dimnames(array_3D), function(x){map_chr(x,~(if (.=="fuzzy.many") "many" else .))})
array_3D <- array_3D / sum(array_3D) * 100

p21 = apply(array_3D, c(1), sum) %>% round(2) %>% map_chr(~paste0(.,"%")) %>% as.data.frame %>% setNames("CB") %>% plot.tab
p22 = apply(array_3D, c(2), sum) %>% round(2) %>% map_chr(~paste0(.,"%")) %>% as.data.frame %>% setNames("UP") %>% plot.tab
p23 = apply(array_3D, c(3), sum) %>% round(2) %>% map_chr(~paste0(.,"%")) %>% as.data.frame %>% setNames("SB") %>% plot.tab

p3 = apply(array_3D, c(1, 3), sum) %>% round(2) %>% apply(c(1,2),function(x) paste0(x, "%")) %>% as.data.frame %>% plot.tab
p4 = apply(array_3D, c(1, 2), sum) %>% t %>% round(2) %>% apply(c(1,2),function(x) paste0(x, "%")) %>% as.data.frame %>% plot.tab
p5 = apply(array_3D, c(2, 3), sum) %>% round(2) %>% apply(c(1,2),function(x) paste0(x, "%")) %>% as.data.frame %>% plot.tab

plot_grid(plot_grid(p21,p22,p23,nrow=1,rel_widths=c(1,1,1)),plot_grid(p3,p4,p5,ncol=1),ncol=1,rel_heights=c(1,3))
ggsave("metrics.pdf",width=5,height=7)


sb.data = data %>% group_by(sb) %>% dplyr::summarize(unique.cb=n(),umi=sum(umi),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi,decreasing=T),]}
cb.data = data %>% group_by(cb) %>% dplyr::summarize(unique.sb=n(),umi=sum(umi)) %>% {.[order(.$umi,decreasing=T),]}


# Make the bead plot
beadplot <- function(sb.data){
  sb.data = sb.data[nrow(sb.data):1,]
  ggplot(sb.data,aes(x=x_um,y=y_um,col=umi)) +
    geom_point(size=1) +
    coord_fixed() +
    theme_classic() +
    labs(x="x (um)", y="y (um)") +
    scale_color_viridis(trans="log", option="B", name="UMI") + 
    ggtitle("Spatial Barcode nUMIs (CB grouped)")
}
beadplot(sb.data)
ggsave("beadplot.pdf",width=6,height=6)

# Make some RNA plots
DimPlot(obj, reduction = "umap") + coord_fixed()
ggsave("umap.pdf",width=6,height=6)

DimPlot(obj, reduction = "spatial") + coord_fixed()
ggsave("cellplot.pdf",width=6,height=6)
       
pdfs = c("searchplot.pdf","params.pdf","metrics.pdf","beadplot.pdf","umap.pdf","cellplot.pdf")
qpdf::pdf_combine(input = pdfs, output = file.path(INDEX,"summary.pdf"))

# ################################################################################

# # plot2 <- plot_grid(
# #   ggplot(cb.data,aes(x=1:nrow(cb.data),y=log10(umi))) + geom_point(size=0.8) + theme_bw() + ggtitle("SB UMI per cell"),
# #   ggplot(sb.data,aes(x=1:nrow(sb.data),y=log10(umi))) + geom_point(size=0.8) + theme_bw() + ggtitle("SB UMI per bead"),
# #   ncol=2
# # ) %>% plot_grid(beadplot(sb.data),ncol=1,rel_heights=c(2,3)) ; print(plot2)
# 
# sum(data$umi)/len(unique(data$cb))
# s = 250/(sum(data$umi)/len(unique(data$cb)))*(sum(data$umi)) ; s
# data$umi = sample(1:nrow(data), size = s, replace = T, prob = data$umi/sum(data$umi)) %>% factor(levels = 1:nrow(data)) %>% table
# sum(data$umi)/len(unique(data$cb))
# data$logumi = log10(data$umi)
# RUN.DBSCAN(data,normalize=F, rpm=T)
# 
# 
# # 
# # obj = "/home/nsachdev/spatial/c.h5" %>% Read10X_h5 %>% CreateSeuratObject %>% add.meta %>% process ; obj %<>% RenameCells(new.names=str_sub(colnames(obj),1,16))
# # obj = add.coords(obj,coords)
# # dimplot(obj)
# # dimplot(obj,split.by="seurat_clusters",ncol=8)&theme_void()
# # DimPlot(obj)+coord_fixed()
# # 

