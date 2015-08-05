################################################################
################################################################
#	query globetrotter admix map
################################################################
################################################################

#load Robj that contains a list with all pertinents for the 
#	admixture runs on globetrotter
load("ms/figs/globe_Ad_object.Robj")

# this function draws the admixture map between all populations
#	if no xlim/ylim values are supplied, it will default to the 
#	min and max of source and target coordinates, which will probably be pretty zoomed out
make.admix.map <- function(xlim=NULL,ylim=NULL,globe_ad_obj){
	with(globe_ad_obj,{
		if(is.null(xlim)){
			xlim <- c(min(target.coords[,1],source.coords[,1]),max(target.coords[,1],source.coords[,1]))
		}
		if(is.null(ylim)){
			ylim <- c(min(target.coords[,2],source.coords[,2]),max(target.coords[,2],source.coords[,2]))
		}
		plot(target.coords,type='n',
					xlim = xlim,
					ylim = ylim,
					xlab="long",
					ylab="lat")
		text(target.coords[c(1:k),],
				labels=pops,
				col=adjustcolor(continent.col,1),
				font=2,cex=0.8)
		text(source.coords[,1],
				source.coords[,2],
				labels=pops,font=3,
				col=globe.admix.plot.cols,
				cex=0.8,family="HersheySerif")
		arrows(	x0 = source.coords[,1],
				y0 = source.coords[,2],
				x1 = target.coords[,1],
				y1 = target.coords[,2],
				col=globe.admix.plot.cols,
				lwd=MAP.admix.props,
				length=0.1)
		box(lwd=2)
		legend(x="topleft",
				lwd = c(1,0.5,0.1),
				col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
				legend = c("w = 0.5","w = 0.25","w = 0.05"),
				title = "Admixture proportions")
		})
}

# this function highlights the admixture source locations for a specified set of populations
#	the pop.name argument takes a character vector, length 1 or more, of the names of the 
#	population(s) whose source(s) of admixture you want to highlight
#
# if you already have a plot going, and you don't want to replot, you can use the add argument
#	which will just add the highlighted source and target locations to the existing map
query.admix.map <- function(pop.name,globe_ad_obj=globe_ad_obj,add=NULL,xlim=NULL,ylim=NULL){
	if(is.null(add)){
		make.admix.map(xlim,ylim,globe_ad_obj)
	}
	with(globe_ad_obj,{
		focal.pop <- match(pop.name,pops)
		text(target.coords[focal.pop,,drop=FALSE],col=1,font=2,labels=pops[focal.pop])
		text(source.coords[focal.pop,,drop=FALSE],col=1,font=3,labels=pops[focal.pop],family="HersheySerif")
		arrows(	x0 = source.coords[focal.pop,1],
				y0 = source.coords[focal.pop,2],
				x1 = target.coords[focal.pop,1],
				y1 = target.coords[focal.pop,2],
				col=globe.admix.plot.cols[focal.pop],
				lwd=MAP.admix.props[focal.pop],
				length=0.1)
	})
	return(0)
}

# all names can be found in globe_ad_obj$pops

# demonstrating using the make.admix.map() function with reasonable xlim vals
	make.admix.map(xlim=c(0,70),ylim=c(-30,50),globe_ad_obj = globe_ad_obj)

# demonstrating using the query.admix.map() function with reasonable xlim vals
	query.admix.map(pop.name = c("SanKhomani"),add=TRUE,globe_ad_obj = globe_ad_obj,xlim=c(15,70),ylim=c(-20,70))
	
	###Zoom in on Europe/North Africa/Middle East
	make.admix.map(xlim=c(42.5,46),ylim=c(34,37.5),globe_ad_obj = globe_ad_obj)
	query.admix.map(pop.name = with(globe_ad_obj, pops[africa]),add=TRUE,globe_ad_obj = globe_ad_obj)
	
	make.admix.map(xlim=c(42.5,46),ylim=c(34,37.5),globe_ad_obj = globe_ad_obj)
	query.admix.map(pop.name = with(globe_ad_obj, pops[africa]),add=TRUE,globe_ad_obj = globe_ad_obj)

	make.admix.map(xlim=c(0,70),ylim=c(-30,50),globe_ad_obj = globe_ad_obj)
	query.admix.map(pop.name = with(globe_ad_obj, pops[western.eurasia]),add=TRUE,globe_ad_obj = globe_ad_obj)
	
# demonstrating using the make.admix.map() function with reasonable xlim vals
	make.admix.map(xlim=c(40,55),ylim=c(32,40),globe_ad_obj = globe_ad_obj)

# demonstrating using the query.admix.map() function with reasonable xlim vals
	query.admix.map(pop.name = c("Tu"),add=TRUE,globe_ad_obj = globe_ad_obj,xlim=c(15,70),ylim=c(-20,70))
