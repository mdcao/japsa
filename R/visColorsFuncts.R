

#.strfunc<-function(string) which(strsplit(string, "")[[1]]=="+")[1]
.strfunc1<-function(string,patt) regexpr(patt, string)[[1]]
.readCSS<-function(arg, samples, names_index = 1,patt="[a-zA-Z]" ){
	 treef = read.table(arg,skip=1, head=F, sep="\t", comment.char='%', fill=T, as.is=T)
	header =  read.table(arg,nrows=1, head=F, sep="\t", comment.char='%', fill=T, as.is=T)
 	treef = treef[grep('----',treef[,1], inv=T),]
	treef = treef[nchar(treef[,2])>0,]
	names(treef) = header
		pos = as.factor(unlist(lapply(treef[,1], .strfunc1, patt)))
	dim_m = dim(treef)
	if(length(names_index)>1){
		nmes = apply(treef[,names_index,drop=F],1,paste, collapse=".")
	}else{
		nmes = treef[,names_index]
	}	
	dupl = which(duplicated(nmes))
	if(length(dupl)>0) {
		print(paste("warning duplicated", nmes[dupl]))
		treef = treef[-dupl,]
		nmes = nmes[-dupl]
		pos = pos[-dupl]
}	
	dimnames(treef)[[1]] = gsub("\\s{2,}", "",nmes)


if(dim(treef)[[2]]>8){
 matrix = apply(treef[,-(1:8)],c(1,2), as.numeric)
}else{
	matrix = array(1, dim = c(dim(treef)[[1]],2))
	dimnames(matrix)[[1]] = treef[,1]
	dimnames(matrx)[[2]] = 1:2
}
indsk = match(samples, dimnames(matrix)[[2]])
matrix = matrix[,indsk]
cols = gsub("css=","",treef[,2])
attr(matrix, "cols") <-cols
attr(matrix,"pos") <-pos
attr(matrix,"taxon")<-treef$taxon
matrix
}

## if cumul then uses the first line to normalise
norm<-function(tab, cumul, sample_inds = 1:dim(tab)[[2]],  mult=1e6, sums = NULL){

	if(is.null(sums)){
top_inds = which(attr(tab,"pos")==5)  ## should include everything
	if(cumul) sums = apply(tab[top_inds,sample_inds],2,sum) else sums = apply(tab[,sample_inds],2,sum)
	}
	for(i in 1:length(sample_inds)){
	  tab[,sample_inds[i]] = mult* (tab[,sample_inds[i]] / sums[[i]])
	}
print(sums)
attr(tab,"sums_") <-sums
	tab
}

getHMCol<-function(len = 128){
  rev(c(rev(colorRampPalette(brewer.pal(9,"Reds"),bias = 0.5)(len)), colorRampPalette(brewer.pal(9,"Blues"), bias = 0.5)(len)))
}

.plotHeatmap<-function(matrix,todo,log=T, addspace=F, Colv=TRUE, sameDendro=T, hm2=F,grps = grps, epsilon=1e-4, offsetRow=0,margins=c(15,5)){
	cols = attr(matrix,"cols") 
	pos = attr(matrix, "pos")
	nr = nrow(matrix)
    	#colinds = list()
	#todo = todo[todo>0]
grpsn = as.numeric(grps)
   cols1 = brewer.pal(length(levels(grps)),"Set1")
cols2 = cols1[grpsn]

	.transformLog<-function(x) log10(x+epsilon)
	.transform1<-function(x) as.numeric(x)
	matrix = apply(matrix,c(1,2),if(log).transformLog else .transform1)
	dendro= NULL

	if(addspace){
		
		lev = levels(pos)
		convertM = cbind(lev,1:length(lev))
		pos1 = as.numeric(convertM[match(pos,lev),2])
		pos1 = pos1 - min(pos1)
		postr = unlist(lapply(pos1, function(x) paste(rep("_", x), collapse="")))
		nme = dimnames(matrix)[[1]]
		dimnames(matrix)[[1]] = apply(cbind( postr,nme),1,paste,collapse="")
	}
	for(i in 1:length(todo)){
		row_inds = which(pos %in% todo[[i]])
		if(length(row_inds)>=2){
if(hm2){
		ncol=dim(matrix)[2]
		nr  = dim(matrix)[1]
		dendro = heatmap.2(matrix[row_inds,],   offsetRow=offsetRow, margins = margins,Rowv=NA, main=paste(todo[[i]], collapse=","),Colv=Colv, RowSideColors=cols[row_inds], dendrogram="column", scale="none", trace="none",  col=getHMCol(),ColSideColors=cols2, keysize=1,
			cexRow = 0.2 + 1/log10(nr),cexCol = 0.2 + 1/log10(nr)
)
}else{
		dendro = heatmap(matrix[row_inds,],main=paste(todo[[i]], collapse=","),Colv=Colv,Rowv = NA,ColSideColors=cols2, RowSideColors=cols[row_inds], dendrogram="column", scale="none", trace="none",  cexCol = 0.5 ,cexRow = 0.5, margins = margins)
}
		if(sameDendro) Colv = dendro$colDendrogram
		#colinds[[i]] = dendro$colInd

		}
	}
	invisible(dendro)
}

.plotBarchart<-function(matrix, 
		todo, thresh = 10,add.text=T,
		show.legend=F,legsize = 10, colInd=1:(dim(matrix)[2])
			){
	ggps = list()
if(is.null(colInd)) colInd=1:(dim(matrix)[2])
	for(i in 1:length(todo)){
		print(todo[[i]])
		row_inds =  which(pos %in% todo[[i]])
		if(length(row_inds)>2){
		ggps[[i]]=try(.barChart(matrix, row_inds ,colInd =colInd,  customCols = T, legsize = legsize, show.legend=show.legend, add.text=add.text, thresh = thresh,  main = paste(todo[[i]], collapse=",")))
		}
	}
	
	invisible(ggps[!is.null(ggps)])
}




.convertForBarChart<-function(tab1_, row_inds = 1:(dim(tab1_)[[1]]), colInd=1:(dim(tab1)[2]), thresh =10){
tab1 = tab1_[row_inds,colInd,drop=F]
 taxid = attr(tab1_,"taxon")[row_inds]
dn = dimnames(tab1)[[1]]
cols = attr(tab1_,"cols") [row_inds]
rem=apply(tab1,1,sum)>0
tab1 = tab1[rem,,drop=F]
taxid = taxid[rem]
dn = dn[rem]
cols = cols[rem]
samples = dimnames(tab1)[[2]]
print(dn)
sample_nme = c()
nrow = length(dn)* length(samples)
tab2 = data.frame(matrix(nrow = nrow, ncol = 7))
start = 1

#tax_lev = tab1$taxon
nrow1  = dim(tab1)[1]
for(i in 1:length(samples)){
	print(i)
	end = length(dn)+start -1
taxid1 = taxid
	pos1 = sum(tab1[,i]) - cumsum(tab1[,i])   
	pos2 = sum(tab1[,i]) - c(0,cumsum(tab1[1:(nrow1-1),i]))
	na_ind = tab1[,i]<thresh
	taxid1[na_ind]=""
	pos1[na_ind] = NA
	pos2[na_ind] = NA
	tab2[start:end,1] = tab1[,i]
	tab2[start:end,2] = dn
	tab2[start:end,3] =  cols
	#total  = sum(tab1[,col_inds[1]])
	#tab2[start:end,1] = tab2[start:end,1]/total
	tab2[start:end,4] = rep(samples[i], length(dn))
	tab2[start:end,5] = taxid1
	tab2[start:end,6] = pos1
	tab2[start:end,7] = (pos2+pos1)/2
	start = end +1
	
}
names(tab2) = c("count", "taxon", "color", "sample","taxid","pos1","pos2")
tab2$taxon = factor(tab2$taxon, levels = dn)
#tab2$color = as.factor(tab2$color)
tab2$sample = factor(tab2$sample, levels = samples)
attr(tab2,"color")=as.character(cols)
#match(levels(tab2$taxon), tab1$taxon)
tab2
}


.barChart<-function(matrix, row_inds,  customCols = T,colInd=1:(dim(matrix)[2]), thresh = thresh, show.legend=F, add.text=T,legsize=5,main= ""){
	tab2 = .convertForBarChart(matrix,  row_inds, colInd, thresh = thresh)
	#samp1 = which(tab2$sample==tab2$sample[dim(tab2)[1]])
	#tab22  = tab2[samp1,]
  #geom_bar(aes(y = percentage, x = year, fill = product), data = charts.data, stat="identity") +
	ggp<-ggplot(tab2) +  geom_bar(aes(fill=taxon, y=count, x=sample), position="stack", stat="identity", show.legend=show.legend)
	if(add.text) ggp<-ggp+geom_text(aes(x=sample,y=pos2, label=taxid), color="black",size=2)
	if(customCols){	
		ggp<-ggp+scale_fill_manual("Legend", values=attr(tab2,"color"))
	}
	ggp<-ggp+theme(axis.text.x = element_text(size=legsize, angle=45, hjust=1), legend.text=element_text(size=legsize), legend.key.size = unit(2.5,"line"),
		legend.box.just = "top")   +ggtitle(main)
	#if(add.text) ggp<-ggp+geom_text(data=tab2, aes(x = sample, y = count, label = taxon), colour="white", family="OfficinaSanITC-Book", size=4)
	invisible(ggp)
}

