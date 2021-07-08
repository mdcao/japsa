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

###############################
path ="/home/lachlan/github/japsa_coverage/R"   ## NEED TO CHANGE 
source(paste(path, "visColorsFuncts.R", sep="/"))

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
args=grep(".css",dir(),v=T)
}

if(FALSE){
	file="aqip_bracken_combined.tsv"
	matrix = norm(.readBracken(file),F)
	matrices = list(matrix, matrix)
	 pos = attr(matrices[[1]], "pos")
	todo = sort(as.numeric(levels(factor(pos))))
	todo1 = list()
	len= length(todo)
	for(i in 1:length(todo)){
		todo1[[i]] = unlist(todo[(len-i+1):len])
	}
indsj = 1:dim(matrices[[1]])[2]
	grps = rep(1, length(indsj))
	resdir = "resdir"
resdir1j = resdir
}else{

suffix="krkn"
design =read.csv("design.csv",sep=",", head=T, as.is=T,row.names=NULL)
design = design[order(design$Name),]
samples = design$Name
resdir = "resdir"
resdir1j = resdir
indsj = 1:length(samples)
grps = rep(1, length(samples))
grps[grep("aqip", samples)] = 2

dir.create(resdir)



	matrices = lapply(args, .readCSS,samples = samples, names_index =1)
	for(i in 1:length(matrices)){
		dimnames(matrices[[i]])[[2]] = unlist(lapply(dimnames(matrices[[i]])[[2]], function(x) strsplit(x,"\\.")[[1]][1]))

	}
	names(matrices) = args
	
	

#sums = attr(matrices[[grep("cumul", args, inv=T)[1]]] ,"sums")
if(normalise){
	matrices[grep("cumul", args)] = lapply(matrices[grep("cumul", args)] , norm, T) 
	matrices[grep("cumul", args, inv=T)] = lapply(matrices[grep("cumul", args, inv=T)] , norm, F) 
}
	 pos = attr(matrices[[1]], "pos")
	todo = sort(as.numeric(levels(factor(pos))))
	todo1 = list()
	len= length(todo)
	for(i in 1:length(todo)){
		todo1[[i]] = unlist(todo[(len-i+1):len])
	}


dimnames(matrices[[1]])[[2]] = gsub(suffix,"",dimnames(matrices[[1]])[[2]])
dimnames(matrices[[2]])[[2]] = gsub(suffix,"",dimnames(matrices[[2]])[[2]])

}

height=0.05 * dim(matrices[[2]])[[1]]+2
width = 12
pdf(paste(resdir1j,"heatmap_all_sep.pdf",sep="/"), height=height, width=width)
dendro =try( .plotHeatmap(matrices[[2]], todo = rev(todo1)[1],log=log, addspace=T, Colv = F,grps = grps[indsj], hm2=hm2, margins = c(10,15)))
 if(inherits(dendro,"try-error")) {
	dendro = NULL
	Colv = NULL
  	col_inds = NULL
}else{

	Colv = dendro$colDendrogram
	col_inds = dendro$colInd
}
dev.off()

pdf(paste(resdir1j,"heatmap_cumul.pdf",sep="/"), height=height, width = width)
 try(.plotHeatmap(matrices[[1]], todo = rev(todo1)[1], addspace=T,log=log,Colv=Colv, grps = grps[indsj], hm2=hm2,margins = c(10,15)))
dev.off()



 pdf(paste(resdir1j,"heatmap_sep.pdf",sep="/"), height = height, width =width)
 try(.plotHeatmap(matrices[[1]], todo = todo,log=log, Colv = Colv,grps = grps[indsj], hm2=hm2,margins = c(10,15)))
dev.off()




if(TRUE){
width = max(10,dim(matrices[[1]])[2]*2+3)

 pdf(paste(resdir1j,"barcharts.pdf",sep="/"), width=width)
ggps = .plotBarchart(matrices[[1]], todo = todo, show.legend=T, legsize=NULL,textsize=5, colInd=col_inds)
for(i in 1:length(ggps)) print(ggps[[i]])
dev.off()
 #pdf(, width=width)
ggps1 = .plotBarchart(matrices[[2]], todo = rev(todo1)[1], show.legend=T, legsize=NULL, colInd=col_inds)
ggsave(paste(resdir1j,"barcharts_all.png",sep="/"), ggps1[[1]], width = 30, height = 50, units = "cm", limitsize=F)
#for(i in 1:length(ggps1)) print(ggps1[[i]])
dev.off()
}

##SAVE PNG
#if(length(ggps1)==1){
#try(ggsave(paste(resdir1j,"barcharts_all.png",sep="/"), plot=ggps1[[1]], width = width, height = 10, units = "in"))
#}
#ggsave(outfile2, plot=ml, width = 30, height = 30, units = "cm")

