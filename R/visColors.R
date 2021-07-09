### CAN RUN WITH R CMD BATCH  or Rscript

rm(list = ls())
INSTALL = FALSE
normalise = FALSE
log=FALSE
if(INSTALL){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("ggplot2")
	BiocManager::install("RColorBrewer")
	BiocManager::install("gplots")
	BiocManager::install("ggrepel")
		BiocManager::install("VGAM")
}

library(stats)
library(gplots)
#library(writexl)
library(ggplot2)
library(RColorBrewer)
#library(ggrepel)

hm2 = TRUE   ## whether to use heatmap2
tab_all = NULL
color_all = NULL
###############################
path ="/home/lachlan/github/japsa_coverage/R"   ## NEED TO CHANGE 
source(paste(path, "visColorsFuncts.R", sep="/"))

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
args=c("cumul.css"   ,       "sep.css" )
}



suffix="krkn"
design =read.csv("design.csv",sep=",", head=T, as.is=T,row.names=NULL)
design = design[order(design$Name),]
samples = design$Name


#grps[grep("aqip", samples)] = 2


	matrices = lapply(args, .readCSS,samples = samples, names_index =1, thresh =NULL)
samples = dimnames(matrices[[1]])[[2]]
resdir = "resdir"
dir.create(resdir)
resdir1j = resdir
indsj = 1:length(samples)
grps = rep(0, length(samples))
nmes = c("species")
if(file.exists("groups")){
groups1 = read.csv("groups", sep=".", header=F, as.is=T)
nmes = c(nmes, groups1[,1])
for(i in 1:(dim(groups1)[[1]])){
  grps[grep(groups1[i,1],samples)] = i

}
}
grps_type = sort(unique(grps))

	names(matrices) = args
	
	

#sums = attr(matrices[[grep("cumul", args, inv=T)[1]]] ,"sums")
if(normalise){
	matrices[grep("cumul", args)] = lapply(matrices[grep("cumul", args)] , norm, T) 
	matrices[grep("cumul", args, inv=T)] = lapply(matrices[grep("cumul", args, inv=T)] , norm, F) 
}
	 



	


for(ij in 1:length(grps_type)){
col_inds = NULL
ij1 = grps_type[ij]
samp_inds = which(grps==ij1)
indsj = samp_inds
matrices_sub = lapply(matrices, .subinds, samp_inds)
names(matrices_sub) = names(matrices)

for(i in 1:length(matrices_sub)){
		dn = dimnames(matrices_sub[[i]])[[2]]
		dn = sub(".fastq", "", dn)
		dn  = sub(".resistancRT.", "_",dn)
		dn  = sub("jR", "",dn)
		dn = unlist(lapply(dn, function(x) strsplit(x,"\\.")[[1]][1]))
		dimnames(matrices_sub[[i]])[[2]] = gsub(suffix,"",dn)
		}


pos = attr(matrices_sub[[1]], "pos")
levs = getlev(pos)

todo2 = sort(as.numeric(as.character(levs[levs[,2]>=1,1])))
todo = sort(as.numeric(as.character(levs[levs[,2]>=5,1])))
todo = list(todo)
	todo1 = list()
	len= length(todo)
	for(i in 1:length(todo)){
		todo1[[i]] = unlist(todo[(len-i+1):len])
	}

 resdir1j = paste(resdir, nmes[ij1+1],sep="/")
dir.create(resdir1j)

height=0.05 * dim(matrices_sub[[2]])[[1]]+2
width = 12
pdf(paste(resdir1j,"/",nmes[ij1+1],".heatmap_all_sep.pdf",sep=""), height=height, width=width)
dendro =try( .plotHeatmap(matrices_sub[[2]], main=paste("sep",nmes[ij1+1]), todo = rev(todo1)[1],log=log, addspace=T, Colv = F,grps = grps[indsj], hm2=hm2, margins = c(10,15)))
 if(inherits(dendro,"try-error")) {
	dendro = NULL
	Colv = NULL
  	col_inds = NULL
}else{

	Colv = dendro$colDendrogram
	col_inds = dendro$colInd
}
dev.off()

pdf(paste(resdir1j,"/",nmes[ij1+1],".heatmap_cumul.pdf",sep=""), height=height, width = width)
 try(.plotHeatmap(matrices_sub[[1]],main=paste("cumul",nmes[ij1+1]), todo = rev(todo1)[1], addspace=T,log=log,Colv=Colv, grps = grps[indsj], hm2=hm2,margins = c(10,15)))
dev.off()

if(FALSE){

 pdf(paste(resdir1j,"/",nmes[ij1+1],".heatmap_sep.pdf", sep=""), height = height, width =width)
 try(.plotHeatmap(matrices_sub[[1]], main=paste("sep",nmes[ij1+1]),todo = todo,log=log, Colv = Colv,grps = grps[indsj], hm2=hm2,margins = c(10,15)))
dev.off()
}




if(TRUE){
width = max(10,dim(matrices_sub[[1]])[2]*2+3)
if(FALSE){
 pdf(paste(resdir1j,"/",nmes[ij1+1],".barcharts.pdf",sep=""), width=width)
tab2 = .plotBarchart(matrices_sub[[1]], main=nmes[ij1+1], todo = todo, colInd=col_inds)
ggps = .barChart(tab2,   show.legend=T, legsize=NULL,textsize=10)
print(ggps)
for(i in 1:length(ggps)) print(ggps[[i]])
dev.off()
}
 #pdf(, width=width)
tab_all1 = .plotBarchart(matrices_sub[[2]],main=nmes[ij1+1], todo = rev(todo1)[1],  colInd=col_inds)
if(is.null(tab_all)){
	tab_all = tab_all1
	color_all = attr(tab_all1,"color")
}else{
	tab_all = cbind(tab_all, tab_all1)
	color_all =c(color_all, attr(tab_all1,"color"))
}
#for(i in 1:length(ggps1)) print(ggps1[[i]])
#dev.off()
}


}
attr(tab_all,"color") = color_all
ggps =  .barChart(tab_all,   show.legend=T, legsize=NULL,textsize=10, facet=T)
ggsave(paste(resdir,"/",".barcharts_all.png",sep=""), ggps, width = 30*length(grps_type), height = 50, units = "cm", limitsize=F)
##SAVE PNG
#if(length(ggps1)==1){
#try(ggsave(paste(resdir1j,"barcharts_all.png",sep="/"), plot=ggps1[[1]], width = width, height = 10, units = "in"))
#}
#ggsave(outfile2, plot=ml, width = 30, height = 30, units = "cm")

