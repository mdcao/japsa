
.readBracken<-function(file="aqip_bracken_combined.tsv", ident="_num"){

tab = read.table(file, sep="\t", head=T)

tab2 = tab[,grep(ident, names(tab))]
dimnames(tab2)[[1]] = tab[,1]

cols = rep("white", dim(tab)[1])
taxon =tab[,2]
pos = rep(5, dim(tab)[1])
attr(tab2,"sums_") = apply(tab2,2,sum)
attr(tab2,"cols") = cols
attr(tab2,"pos") = factor(pos)
attr(tab2,"taxon") = taxon
tab2
}

#.strfunc<-function(string) which(strsplit(string, "")[[1]]=="+")[1]
.strfunc1<-function(string,patt) regexpr(patt, string)[[1]]
.readCSS<-function(arg, samples, names_index = 1,patt="[a-zA-Z]" ){
	 treef = read.csv(arg,skip=1, head=F, sep="\t", comment.char='%', fill=T, as.is=T)
	header =  read.table(arg,nrows=1, head=F, sep="\t", comment.char='%', fill=T, as.is=T)
 	treef = treef[grep('----',treef[,1], inv=T),]
	treef = treef[nchar(treef[,2])>0,]

	names(treef) = header
	leaves = treef$level1==0
	#last_tp = unlist(lapply(treef$taxon_parents, function(x) rev(unlist(strsplit(x,",")[[1]])))
	#parent_ind = unlist(lapply()[1]), function(x) which(x==treef$taxon)))

		pos = as.factor(unlist(lapply(treef[,1], .strfunc1, patt)))
#pos = treef$level+1
	dim_m = dim(treef)
	if(length(names_index)>1){
		nmes = apply(treef[,names_index,drop=F],1,paste, collapse=".")
	}else{
		nmes = treef[,names_index]
	}	
	nmes=gsub("\\s{2,}", "",nmes)
	dupl = which(duplicated(nmes))
	if(length(dupl)>0) {
		stop("duplication")
}	
	dimnames(treef)[[1]] = nmes

if(dim(treef)[[2]]>9){
 matrix = apply(treef[,-(1:9)],c(1,2), as.numeric)
}else{
	matrix = array(1, dim = c(dim(treef)[[1]],2))
	dimnames(matrix)[[1]] = treef[,1]
	dimnames(matrx)[[2]] = 1:2
}
indsk = match(samples, dimnames(matrix)[[2]])
matrix = matrix[,indsk]
taxon = treef$taxon

 
if(arg=="sep.css"){
leaves1 = which(leaves)
len1 = dim(matrix)[[2]]
for(leaf in leaves1){
	pj = which(treef$taxon ==rev(strsplit(treef$taxon_p[leaf],",")[[1]])[1])
	matrix[pj,] = matrix[pj,] + matrix[leaf,]
	matrix[leaf,] = rep(0, len1)
}
 
}
 
#if(args=="sep.css")
cols = gsub("css=","",treef[,2])
cols[cols=="null"] = "white"
matrix = matrix[!leaves,]
  cols = cols[!leaves]
  pos = pos[!leaves]
  taxon = taxon[!leaves]
attr(matrix, "cols") <-cols
attr(matrix,"pos") <-pos
attr(matrix,"taxon")<-taxon
#attr(matrix,"leaves")<- leaves
#attr(matrix,"parent_ind") <-parent_ind
matrix
}

## if cumul then uses the first line to normalise
norm<-function(tab, cumul, sample_inds = 1:dim(tab)[[2]],  mult=1e6, sums = NULL){

	if(is.null(sums)){
top_inds = which(attr(tab,"pos")==1)  ## should include everything
	if(cumul) sums = apply(tab[top_inds,sample_inds,drop=F],2,sum) else sums = apply(tab[,sample_inds],2,sum)
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

.plotHeatmap<-function(matrix,todo,log=T, addspace=F, Colv=TRUE,basec = 0.3, sameDendro=T, hm2=F,grps = grps, epsilon=1e-4, offsetRow=0,margins=c(15,5) ){
	cols = attr(matrix,"cols")
	cols[cols=="null"] = "white"
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
			cexRow = basec + 1/log10(nr),cexCol = basec + 1/log10(nr)
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
		todo, thresh = 10,add.text=F,
		show.legend=F,legsize = 10, textsize=legsize, colInd=1:(dim(matrix)[2])
			){
	ggps = list()
if(is.null(colInd)) colInd=1:(dim(matrix)[2])
	for(i in 1:length(todo)){
		print(todo[[i]])
		row_inds =  which(pos %in% todo[[i]])
		if(length(row_inds)>1){
		ggps[[i]]=try(.barChart(matrix, row_inds ,colInd =colInd,  customCols = T, legsize = legsize, show.legend=show.legend, add.text=add.text, textsize=textsize, thresh = thresh,  main = paste(todo[[i]], collapse=",")))
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


.barChart<-function(matrix, row_inds,  customCols = T,colInd=1:(dim(matrix)[2]), thresh = thresh, show.legend=F, add.text=F,legsize=5,main= "", textsize = 5){
	tab2 = .convertForBarChart(matrix,  row_inds, colInd, thresh = thresh)
	numcats = length(levels(tab2$taxon))
	if(is.null(legsize)) legsize = 15/log(numcats+1)
	#samp1 = which(tab2$sample==tab2$sample[dim(tab2)[1]])
	#tab22  = tab2[samp1,]
  #geom_bar(aes(y = percentage, x = year, fill = product), data = charts.data, stat="identity") +
	ggp<-ggplot(tab2) +  geom_bar(aes(fill=taxon, y=count, x=sample), position="stack", stat="identity", show.legend=show.legend)
	if(add.text) ggp<-ggp+geom_text(aes(x=sample,y=pos2, label=taxid), color="black",size=2)
	if(customCols){	
		ggp<-ggp+scale_fill_manual("Legend", values=attr(tab2,"color"))
	}
	ggp<-ggp+theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text.x = element_text(size=textsize, angle=45, hjust=1), legend.text=element_text(size=legsize), legend.position="bottom", legend.key.size = unit(legsize/2,"line"),
		legend.box.just = "bottom")   +ggtitle(main)

	#if(add.text) ggp<-ggp+geom_text(data=tab2, aes(x = sample, y = count, label = taxon), colour="white", family="OfficinaSanITC-Book", size=4)
	invisible(ggp)
}

