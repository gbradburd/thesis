library(igraph)
iArrows <- igraph:::igraph.Arrows   #http://kbroman.wordpress.com/2012/10/11/curved-arrows-in-r/

alpha=1
par(mar=c(0,3.5,0.1,1))
loc<-c(0,0.75,2.5,3.5)

covar.d<-function(d){
exp(-alpha*sqrt(d^2)^1.3)
}
cols<-c("red","black","blue","orange")

d<-seq(0,loc[4],length=1000)
plot(d,covar.d(d-loc[1]),type="l",col=cols[1],lwd=2,ylim=c(-0.3,1),axes=FALSE,xlab="",ylab="",cex=1.5)
axis(side=2,at=c(0,0.5,1),cex=1.2)
mtext(side=2,line=2,"Covariance",cex=1.4)
abline(h=0)

lines(d,covar.d(d-loc[3]),col=cols[3],lwd=2)
lines(d,covar.d(d-loc[4]),col=cols[4],lwd=2)
w1<-0.4

Bs.covar<-c((1-w1)*covar.d(loc[1]-loc[2]) + w1*covar.d(loc[1]-loc[4]),  ##with B-A
					(1-w1)*(1-w1)*1 + w1*w1*1 + 2*(1-w1)*w1*covar.d(loc[4]-loc[2]),  ##B-B
					(1-w1)*covar.d(loc[3]-loc[2]) + w1*covar.d(loc[3]-loc[4]),  ##B-C
					(1-w1)*covar.d(loc[4]-loc[2]) + w1*covar.d(0)   ##B-D
)

text(x=rep(loc[2]-0.18,4),y=Bs.covar,col=cols[2],"B-",cex=1.5) 
text(x=rep(loc[2]-0.09,4),y=Bs.covar,col=cols,c("A","B","C","D"),cex=1.5) 

#text(x=rep(loc[2]-0.05,4),y=Bs.covar,col=cols[2],"B",cex=1.5) 
#text(x=rep(loc[2]+0.05,4),y=Bs.covar,col=cols,c("-A","-B","-C","-D"),cex=1.5) 
points(x=rep(loc[2],4),y=Bs.covar,col=cols,pch=19,cex=1.5) 

#abline(h=(1-w1)*covar.d(loc[1]-loc[2]) + w1*covar.d(loc[1]-loc[4]),col=cols[1],lty=2)
#abline(h= (1-w1)*covar.d(loc[3]-loc[2]) + w1*covar.d(loc[3]-loc[4]),col=cols[3],lty=2)
#abline(h= (1-w1)*covar.d(loc[4]-loc[2]) + w1*covar.d(0),col=cols[4],lty=2)
#abline(h= (1-w1)*(1-w1)*1 + w1*w1*1 + 2*(1-w1)*w1*covar.d(loc[4]-loc[2]),col=cols[2],lty=2)


text(loc,rep(-0.1,4),c("A","B","C",'D'),col=cols,cex=1.5)
iArrows(loc[4], -0.18, loc[2], -0.18,
          h.lwd=2, sh.lwd=2, sh.col="black",
          curve=0.15 , width=1, size=0.7)
          dev.copy2pdf(file="~/Dropbox/Students/gideon/spacemix_ms/figs/sims/Admix_covar_toy_fig.pdf")
          
