library("Sim.DiffProc", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
library(stringr)
fileList= list.files(path="dat",pattern="resRE.*");
bigTs <- seq(0,10,length.out=5000)
keyvalue <- list(); 
REs <- c();
vars <- c(); means <- c(); meanFPTs <- c()
pdfsurf <- matrix(ncol = 5000,nrow = 0)
for (file in fileList){
	load(sprintf("dat/%s",file)); 
	RE <- as.numeric(str_match(file,"\\d*\\.\\d*")); 
	
	REs <- append(REs , RE)
	vars <- append(vars,var(t(tail(res$X,n=1))))
	means <- append(means,mean(t(tail(res$X,n=1))))
	traces <- res$X
	FPs <- c()
	for (ii in 1:dim(traces)[2]){
	  trace <- traces[,ii]
	  #plot(trace)
	  FPs <- append(FPs , min(which(trace > -50.)))
	}
  meanFPTs <- append(meanFPTs , mean(bigTs[FPs[is.finite(FPs)]]))
}

p <- plot_ly(x = ~bigTs , y = ~REs ,z = ~pdfsurf) %>% add_surface()
p

results <- data.frame(REs=REs,vars=vars,means=means,meanFPTs=meanFPTs)

g1 <- ggplot(results , aes(REs , means)) + geom_line(color="grey") + geom_point(color="firebrick")
g1 <- g1+labs(x=expression(paste(lambda[e](1/s))), 
              y=expression(paste("Mean",(V))))#, title="Temperature")
g1 <- g1 + coord_trans(x="log10") + theme(axis.text.x=element_text(angle=50, vjust=0.5, size=20),axis.text.y=element_text(angle=50, size=20))
g1 <- g1 + theme(axis.title.x=element_text( size=20),axis.title.y=element_text( size=20))
g1 <- g1 + ylim(c(-60,-50))

g2 <- ggplot(results , aes(REs , vars)) + geom_line(color="grey")  + geom_point(color="firebrick")
g2 <- g2+labs(x=expression(paste(lambda[e](1/s))), 
              y=expression(paste("Variance ",(V^2))))#, title="Temperature")
g2 <- g2 + coord_trans(x="log10")  + theme(axis.text.x=element_text(angle=50, vjust=0.5, size=20),axis.text.y=element_text(angle=50, size=20))
g2 <- g2 + theme(axis.title.x=element_text( size=20),axis.title.y=element_text( size=20))

g3 <- ggplot(results , aes(REs , meanFPTs)) + geom_line(color="grey")  + geom_point(color="firebrick")
g3 <- g3+labs(x=expression(paste(lambda[e](1/s))), 
              y=expression(paste(E(T^sigma))))#, title="Temperature")
g3 <- g3 + coord_trans(x="log10")  + theme(axis.text.x=element_text(angle=50, vjust=0.5, size=20),axis.text.y=element_text(angle=50, size=20))
g3 <- g3 + theme(axis.title.x=element_text( size=20),axis.title.y=element_text( size=20))


grid.arrange(g1,g2,g3,ncol=3)
