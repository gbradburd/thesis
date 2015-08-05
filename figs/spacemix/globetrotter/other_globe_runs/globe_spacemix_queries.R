################################################################
################################################################
#	make figures for SpaceMix paper
################################################################
################################################################

procrusteez <- function(obs.locs,target.locs,k,source.locs = NULL,option){
	require(vegan)
	proc.loc <- procrustes(obs.locs,target.locs,scale=TRUE)
	if(option==1){
		proc.pop.loc <- proc.loc$scale * target.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	} else if(option==2){
		proc.pop.loc <- proc.loc$scale * source.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	}
	return(proc.pop.loc)	
}

get.procrustes.locations.posterior.list <- function(observed.coords,population.coordinates.posterior){
	target.coords.list <- vector(mode="list",length = length(population.coordinates.posterior))
	source.coords.list <- vector(mode="list",length = length(population.coordinates.posterior))
	k <- nrow(observed.coords)
	for(i in 1:length(target.coords.list)){
		target.coords.list[[i]] <- procrusteez(obs.locs = observed.coords,
												target.locs = population.coordinates.posterior[[i]][1:k,],
												k = k,
												option = 1)
		source.coords.list[[i]] <- procrusteez(obs.locs = observed.coords,
												target.locs = population.coordinates.posterior[[i]][1:k,],
												k = k,
												source.locs = population.coordinates.posterior[[i]][(k+1):(2*k),],
												option = 2)
	}
	return(list(target.coords.list= target.coords.list, source.coords.list= source.coords.list))
}

fade.admixture.source.points <- function(pop.cols,admix.proportions){
	faded.colors <- numeric(length(pop.cols))
	for(i in 1:length(pop.cols)){
		faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i])
	}
	return(faded.colors)
}

load_MCMC_output <- function(MCMC.output.file){
    tmpenv <- environment()
	tmp <- load(MCMC.output.file,envir=tmpenv)
	mcmc.output <- lapply(tmp,get,envir=tmpenv)
	names(mcmc.output) <- tmp
	return(mcmc.output)
}

get.posterior.location.matrix.from.list <- function(posterior.list,population.index){
	post.location.matrix <- matrix(unlist(
								lapply(posterior.list,
									FUN=function(elem){elem[population.index,]})),
								nrow=length(posterior.list),ncol=2,byrow=TRUE)
	return(post.location.matrix)
}

get.credible.ellipse <- function(posterior.points,quantile){
	require(MASS)
	require(cluster)
	fit <- cov.mve(posterior.points, quantile.used = nrow(posterior.points) * quantile)
	points_in_ellipse <- posterior.points[fit$best, ]
	ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
	return(ellipse_boundary)
}

plot.credible.ellipse <- function(ellipse_boundary,population.color,fading=0.3,lty=1){
	polygon(ellipse_boundary,col=adjustcolor(population.color,fading),border=1,lty=lty)
}

make.spacemix.map.list <- function(MCMC.output.file,observed.coords,name.vector,color.vector,quantile=0.95){
	MCMC.output <- load_MCMC_output(MCMC.output.file)
	best <- which.max(MCMC.output$Prob)
	admix.source.color.vector <- fade.admixture.source.points(color.vector,rowMeans(MCMC.output$admix.proportions))
	k <- MCMC.output$last.params$k
	target.coords <- procrusteez(observed.coords,MCMC.output$population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(observed.coords,MCMC.output$population.coordinates[[best]][1:k,],k,
									source.locs=MCMC.output$population.coordinates[[best]][(k+1):(2*k),],option=2)
	procrustes.coord.posterior.lists <- get.procrustes.locations.posterior.list(observed.coords= observed.coords,
																				population.coordinates.posterior=MCMC.output$population.coordinates)
	posterior.target.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$target.coords.list)
	posterior.source.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$source.coords.list)
	posterior.target.ellipses <- lapply(posterior.target.location.matrices,get.credible.ellipse,quantile)
	posterior.source.ellipses <- lapply(posterior.source.location.matrices,get.credible.ellipse,quantile)
	spacemix.map.list <- c(MCMC.output,
							list(observed.coords=observed.coords),
								list(name.vector=name.vector),list(color.vector=color.vector),
								list(quantile=quantile),list(source=source),list(best = best),
								list(admix.source.color.vector = admix.source.color.vector),
								list(k = k),list(target.coords = target.coords),list(source.coords = source.coords),
								list(procrustes.coord.posterior.lists = procrustes.coord.posterior.lists),
								list(posterior.target.location.matrices = posterior.target.location.matrices),
								list(posterior.source.location.matrices = posterior.source.location.matrices),
								list(posterior.target.ellipses = posterior.target.ellipses),
								list(posterior.source.ellipses = posterior.source.ellipses))
	return(spacemix.map.list)
}

make.spacemix.map <- function(spacemix.map.list,text=FALSE,ellipses=TRUE,source.option=TRUE,xlim=NULL,ylim=NULL){
	with(spacemix.map.list,{ 
		plot(target.coords,type='n',xlim=xlim,ylim=ylim,xlab="",ylab="")
			if(ellipses){
				lapply(1:k,FUN=function(i){plot.credible.ellipse(posterior.target.ellipses[[i]],color.vector[i])})
			}
			if(text){
				text(target.coords,col=color.vector,font=2,labels=name.vector,cex=0.7)
			}
			if(source.option){
				if(ellipses){
					lapply(1:k,FUN=function(i){plot.credible.ellipse(posterior.source.ellipses[[i]],admix.source.color.vector[i],fading=1,lty=2)})
				}
				text(source.coords,col= admix.source.color.vector,font=3,labels=name.vector,cex=0.7)
				arrows(	x0 = source.coords[,1],
						y0 = source.coords[,2],
						x1 = target.coords[,1],
						y1 = target.coords[,2],
						col= admix.source.color.vector,
						lwd= spacemix.map.list $admix.proportions[,best],
						length=0.1)
			}
				box(lwd=2)
	})
}

query.spacemix.map <- function(focal.pops,spacemix.map.list,source.option=TRUE){
	with(spacemix.map.list,{
		# browser()
		focal.indices <- match(focal.pops,name.vector)
			for(i in 1:length(focal.indices)){
				plot.credible.ellipse(posterior.target.ellipses[[focal.indices[i]]],color.vector[focal.indices[i]],fading=1)
			}
			if(source.option){
				for(i in 1:length(focal.indices)){
					plot.credible.ellipse(posterior.source.ellipses[[focal.indices[i]]], admix.source.color.vector[focal.indices[i]],fading=1,lty=2)
				}
				text(source.coords[focal.indices,,drop=FALSE],col=1,font=3,labels=name.vector[focal.indices])
				arrows(	x0 = source.coords[focal.indices,1],
						y0 = source.coords[focal.indices,2],
						x1 = target.coords[focal.indices,1],
						y1 = target.coords[focal.indices,2],
						col= admix.source.color.vector[focal.indices[i]],
						lwd=1,
						length=0.1)
			}
			text(target.coords[focal.indices,,drop=FALSE],col=1,font=2,labels=name.vector[focal.indices],cex=1)
				box(lwd=2)
	})
}

################
#	GLOBE preliminaries
################
load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/rand_prior2/globetrotter_dataset.Robj")

	globe.coords <- cbind(globetrotter.long, globetrotter.lat)
	pops <- row.names(globetrotter.counts)
	k <- length(pops)

		continent.col <- numeric(k)
			continent.col[which(globetrotter.long < -50)] <- "orange"
			continent.col[match(c("BantuKenya","BantuSouthAfrica","BiakaPygmy",
									"Egyptian","Ethiopian","EthiopianJew","Hadza","Mandenka",
									"MbutiPygmy","Moroccan","Mozabite","Sandawe","SanNamibia",
									"SanKhomani","Tunisian","Yoruba"),pops)] <- "forestgreen"
		continent.col[which(globetrotter.long > 100 &
							globetrotter.lat < 5)] <- "brown"
		continent.col[which(continent.col==0)] <- rainbow(length(continent.col[which(continent.col==0)]),
															start=4/6,end=6/6)[as.numeric(cut(globetrotter.long[which(continent.col==0)],length(which(continent.col==0))))]
		americas <- which(continent.col=="orange")
		africa <- which(continent.col=="forestgreen")
		oceania <- which(continent.col=="brown")
		east.asia <- which(globetrotter.long > 95 & 
								globetrotter.lat > 11.5)
		western.eurasia <- c(1:k)[-c(americas,africa,oceania,east.asia)]
		eurasia <- c(1:k)[-c(americas,africa,oceania)]

		# af.eff.lat <- (globetrotter.lat[africa] + abs(min(globetrotter.lat[africa])))/max(globetrotter.lat[africa] + abs(min(globetrotter.lat[africa])))
		af.loc.cols <- hsv(h = seq(0.22,0.69,length.out=length(africa))[rank(globetrotter.lat[africa])],
				s = 1,
				v = 1)
		adj.nonamaf.long <- globetrotter.long[-c(africa,americas)] + abs(min(globe.coords[western.eurasia,1]))
		eur.eff.long <- adj.nonamaf.long/max(adj.nonamaf.long)
		eur.loc.cols <- hsv(h = eur.eff.long * 0.4 + 0.6,s=1,v=1)
		am.eff.long <- (globetrotter.long[americas] + abs(min(globetrotter.long[americas])))/max(globetrotter.long[americas] + abs(min(globetrotter.long[americas])))
		am.loc.cols <- hsv(h = (am.eff.long) * 0.08 + 0.03,s=1,v=1)
		continent.col <- numeric(k)
		continent.col[americas] <- am.loc.cols
		continent.col[africa] <- af.loc.cols
		continent.col[-c(africa,americas)] <- eur.loc.cols

	pop.order <- pop.order <- c(africa[order(globe.coords[africa,2])],
					western.eurasia[order(globe.coords[western.eurasia,1])],
					east.asia[order(globe.coords[east.asia,1])],
					oceania[rev(order(globe.coords[oceania,2]))],
					americas[rev(order(globe.coords[americas,2]))])

################
#	RandPr2
################

setwd("~/Desktop/globe_spacemix_query_figs")
globe.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/rand_prior2/globe_spaceruns_randpr1_LongRun/globe_spaceruns_randpr1space_MCMC_output1.Robj",
							observed.coords = globe.coords,
							name.vector = pops,
							color.vector = continent.col,
							quantile=0.95)
for(i in 1:k){
	cat(i,"\t")
		xlim <- c(min(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,1],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,1]),
					max(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,1],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,1]))
		ylim <- c(min(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,2],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,2]),
					max(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,2],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,2]))
	plotting.width <- 9
	plotting.height <- max(5,9/(diff(range(abs(xlim)))/diff(range(abs(ylim)))))
	png(file=paste(pops[pop.order[i]],"_randpr2",".png",sep=""),res=300,width=plotting.width*300,height=plotting.height*300)
		make.spacemix.map(spacemix.map.list = globe.spacemix.map.list,
							text=TRUE,
							ellipses=FALSE,
							source.option=TRUE,
							xlim=xlim,
							ylim=ylim)
		mtext(side=3,text=pops[pop.order[i]],font=2,cex=2)
		query.spacemix.map(focal.pops=pops[pop.order[i]],spacemix.map.list= globe.spacemix.map.list)
	dev.off()
}


################
#	RealPr2
################

globe.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/real_prior2/globe_spaceruns_realpr2_LongRun/globe_spaceruns_realpr2space_MCMC_output1.Robj",
							observed.coords = globe.coords,
							name.vector = pops,
							color.vector = continent.col,
							quantile=0.95)
for(i in 1:k){
	cat(i,"\t")
		xlim <- c(min(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,1],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,1]),
					max(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,1],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,1]))
		ylim <- c(min(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,2],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,2]),
					max(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,2],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,2]))
	plotting.width <- 9
	plotting.height <- max(5,9/(diff(range(abs(xlim)))/diff(range(abs(ylim)))))
	png(file=paste(pops[pop.order[i]],"_realpr2",".png",sep=""),res=300,width=plotting.width*300,height=plotting.height*300)
		make.spacemix.map(spacemix.map.list = globe.spacemix.map.list,
							text=TRUE,
							ellipses=FALSE,
							source.option=TRUE,
							xlim=xlim,
							ylim=ylim)
		mtext(side=3,text=pops[pop.order[i]],font=2,cex=2)
		query.spacemix.map(focal.pops=pops[pop.order[i]],spacemix.map.list= globe.spacemix.map.list)
	dev.off()
}


################
#	RealPr3
################

globe.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/real_prior3/globe_spaceruns_realpr3_LongRun/globe_spaceruns_realpr3space_MCMC_output1.Robj",
							observed.coords = globe.coords,
							name.vector = pops,
							color.vector = continent.col,
							quantile=0.95)
for(i in 1:k){
	cat(i,"\t")
		xlim <- c(min(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,1],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,1]),
					max(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,1],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,1]))
		ylim <- c(min(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,2],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,2]),
					max(globe.spacemix.map.list$posterior.target.ellipses[[pop.order[i]]][,2],
						globe.spacemix.map.list$posterior.source.ellipses[[pop.order[i]]][,2]))
	plotting.width <- 9
	plotting.height <- max(5,9/(diff(range(abs(xlim)))/diff(range(abs(ylim)))))
	png(file=paste(pops[pop.order[i]],"_realpr3",".png",sep=""),res=300,width=plotting.width*300,height=plotting.height*300)
		make.spacemix.map(spacemix.map.list = globe.spacemix.map.list,
							text=TRUE,
							ellipses=FALSE,
							source.option=TRUE,
							xlim=xlim,
							ylim=ylim)
		mtext(side=3,text=pops[pop.order[i]],font=2,cex=2)
		query.spacemix.map(focal.pops=pops[pop.order[i]],spacemix.map.list= globe.spacemix.map.list)
	dev.off()
}