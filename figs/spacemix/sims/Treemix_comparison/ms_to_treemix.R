
setwd("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/")
ms.files<-dir("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/",pattern="*Robj")

###run to convert the ms file Robjects into treemix format
for(i in 1:length(ms.files)){
	ms.file<-strsplit(ms.files[i],".Robj")[[1]]
	convert.ms.treemix(ms.file)
}


####RUN & PLOT Treemix trees.
setwd("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/")
for(i in 1:length(ms.files)){
	ms.file<-strsplit(ms.files[i],".Robj")[[1]]
	treemix.file<-  paste("treemix_", ms.file,".out",sep="") 
	 system(paste("treemix -i ",treemix.file,  ".gz -root myoutgroup "," -o ","treemix_output/",ms.file,"_treemixed", sep=""))
}

# pdf(file="~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_trees.pdf")
# for(i in 1:length(ms.files)){
	# ms.file<-strsplit(ms.files[i],".Robj")[[1]]
	# plot_tree(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_output/",ms.file,"_treemixed", sep=""))
	# mtext(ms.file)
# }
# dev.off()


setwd("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/")

for(num.migs in 1:2){
	for(i in 1:length(ms.files)){
		ms.file<-strsplit(ms.files[i],".Robj")[[1]]
		cat(ms.file," ",num.migs,"\n")
		treemix.file<-  paste("treemix_", ms.file,".out",sep="") 
		###To run each tree separately
		#system(paste("treemix -i ",treemix.file,  ".gz -root myoutgroup -m ",num.migs," -o ","treemix_output/",ms.file,"_treemixed_m",num.migs, sep=""))
		if(num.migs==1) prev.treemix<-paste("treemix_output/",ms.file,"_treemixed",sep="");
		if(num.migs> 1) prev.treemix<-paste("treemix_output/",ms.file,"_continuing_treemixed_m",num.migs-1,sep="");
		##adding 1 mig at time
		system(paste("treemix -i ",treemix.file,  ".gz -root myoutgroup  -g ",prev.treemix, ".vertices.gz ",prev.treemix, ".edges.gz -m 1 -o ","treemix_output/",ms.file,"_continuing_treemixed_m",num.migs, sep=""))
	}
}


###PLOT trees & graphs for each scenario 
pdf(file=paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_trees_mig_",num.migs,".pdf",sep=""))
for(i in 1:length(ms.files)){
	ms.file<-strsplit(ms.files[i],".Robj")[[1]]
	plot_tree(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_output/",ms.file,"_treemixed_m",num.migs, sep=""))
	mtext(ms.file)
}
dev.off()


for(i in 1:length(ms.files)){
	ms.file<-strsplit(ms.files[i],".Robj")[[1]]

	pdf(file=paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/",ms.file,"treemix_fig_sequential.pdf",sep=""),width=12,height=8)
	layout(matrix(1:6,nrow=2,byrow=TRUE))
	if(ms.file == "ms_dataset_stationary_pops") migration.rate.graphic(x.pops = 5,y.pops = 6,migration.rate=1,migration.arrows=FALSE,jitter=0.25,labels=TRUE,colors=TRUE,arrow.width=1,pop.pt.cex=2.5,pop.lab.cex=2.5)
	if(ms.file == "ms_dataset_barrier_and_inland_admixture"){
		expansion.list2 <- list(list(parent = 61,
							daughters = 105,
							time.point = 1))
		migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,migration.arrows=FALSE,barrier.effect=5,jitter = 0.25,expansion.list=expansion.list2,labels=TRUE,colors=TRUE,arrow.col="green",arrow.width=0.4,pop.pt.cex=2.5,pop.lab.cex=2.5)
		abline(v=6.2,lwd=6,lty=2)
			}
	
	if(ms.file =="ms_dataset_barrier"){
			migration.rate.graphic(x.pops = 5,y.pops = 6,migration.rate=1,migration.arrows=FALSE,jitter=0.25,barrier.effect=5,labels=TRUE,colors=TRUE,pop.pt.cex=2.5,pop.lab.cex=2.5)
		abline(v=6.2,lwd=6,lty=2)
		
	}
	
	if(ms.file =="ms_dataset_corner_admixture"){
		expansion.list3 <- list(list(parent = 13,
							daughters = 131,
							time.point = 1))
		migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,migration.arrows=FALSE,jitter = 0.25,expansion.list=expansion.list3,labels=TRUE,colors=TRUE,arrow.col="green",curve=0.2,arrow.width=0.8,pop.pt.cex=2.5,pop.lab.cex=2.5)
	}
	
	if(ms.file =="ms_dataset_expansion"){
	parents <- c(78:88)
	time.points <- rep(0.07,11)
	expansion.list <- vector(mode="list")
		for(i in 1:length(parents)){
			expansion.list[[i]] <- list(parent=parents[i],
										daughters = parents[i]+c(11,22,33,44,55),
										time.point = time.points[i])
		}
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,migration.arrows=FALSE,jitter = 0.25,expansion.list=expansion.list,labels=TRUE,colors=TRUE,ylim=c(0,10.2),curve=0.3,pop.pt.cex=2.5,pop.lab.cex=2.5)	
	}
	
	if(ms.file =="ms_dataset_neighbor_admixture"){
	expansion.list2 <- list(list(parent = 61,
							daughters = 105,
							time.point = 1))
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,migration.arrows=FALSE,barrier.effect=5,jitter = 0.25,expansion.list=expansion.list2,labels=TRUE,colors=TRUE,arrow.col="green",arrow.width=0.4,pop.pt.cex=2.5,pop.lab.cex=2.5)
		abline(v=6.2,lwd=6,lty=2)	
	}
	
	##plot first tree
	plot_tree(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_output/",ms.file,"_treemixed", sep=""))
	##get ordered list of names of populations
	tmp<-read.table(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_output/",ms.file,"_treemixed.treeout.gz",sep=""),nrow=1,as.is=TRUE)[1,]
	tmp2<-strsplit(tmp,"[[:punct:]]")
	these<-sort(c(grep("outgroup",tmp2[[1]]),grep("pop",tmp2[[1]])))
	pops<-tmp2[[1]][these]
	
	write.table(file="temp.pops",pops,quote=FALSE,row.names=FALSE,col.names=FALSE)
	plot_resid(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_output/",ms.file,"_treemixed", sep=""),"temp.pops")	
	pops[pops=="myoutgroup"]<-"OG"
 	pops[pops!="OG"]<-sapply(pops[pops!="OG"],function(x){strsplit(x,"pop")[[1]][2]})
	npop<-length(pops)
 	mtext(pops, side = 2, at = 1-(1:npop)/npop+0.5/npop, las = 1, cex = 0.7)
    mtext(pops, side = 1, at =  (1:npop)/npop-0.5/npop, las = 3, cex = 0.7)
	for(num.migs in 1:3){
			###remove the continuing to plot trees generated separately.
			plot_tree(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_output/",ms.file,  "_continuing_treemixed_m",num.migs, sep=""))
	}
	dev.off()

}

convert.ms.treemix<-function(file){
	print(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/",file,".Robj",sep=""))
	## pops x SNPs 
	 show(load(paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/",file,".Robj",sep="")))
	print(dim(spacemix.dataset$allele.counts))
	allele.1<- spacemix.dataset$allele.counts 
	allele.2<-spacemix.dataset$sample.sizes - spacemix.dataset$allele.counts 
	
	treemix.data<-character()
	 for(pop in 1:nrow(allele.1)){
	 	pop.name<-paste("pop",pop,sep="")
	 	temp<-c(pop.name, paste(allele.1[pop,], allele.2[pop,],sep=","))
		treemix.data<-cbind(treemix.data,temp) 	
	 	
	 }
	
	 #write.table(file="~/Dropbox/Students/gideon/spacemix_analysis/treemix_data_grid.out", treemix.data,quote=FALSE,col.names=FALSE,row.names=FALSE)
	# system("gzip ~/Dropbox/Students/gideon/spacemix_analysis/treemix_data_grid.out")
	
	# treemix -i treemix_data_grid.out.gz -o out_stem
	
	##Add outgroup
	temp<-c("myoutgroup",paste(rep(0,ncol(allele.1)),rep(1,ncol(allele.1)) ,sep=","))
	treemix.data<-cbind(treemix.data,temp) 	
	treemix.file<-  paste("~/Dropbox/Students/gideon/spacemix_analysis/ms_spacemix/treemix_", file,".out",sep="") 
	 write.table(file=treemix.file, treemix.data,quote=FALSE,col.names=FALSE,row.names=FALSE)
	 system(paste("gzip ",treemix.file,sep=""))
}
 

 