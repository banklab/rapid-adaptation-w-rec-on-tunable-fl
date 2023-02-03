args=commandArgs(trailingOnly = TRUE)


popSize=c(100, 500, 5000)
add=c(0.01)
#add=c(0.001)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "distance2initial.pdf", sep="."), width=2, height=2)

############# plot fixation time
par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
lwidth=c(0.5, 0.6)
dots=c(17, 19)

k=1

yup=rep(1,6) # a=0.01
ylow=rep(0,6) # a=0.01

for(a in add){
    for(e in epi){
        boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(ylow[k],yup[k]))
        title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
        mtext("Recombination probability", side=1, line=1.5, cex=0.7)
        if(k==1 || k==4 || k==7){
            mtext("Diverg. from initial pop.", side=2, line=1.8, cex=0.8, las=0) }
        axis(side=1, label=c(0, 0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85,2), cex.axis=0.75)
        
        
        #if(k==3){
         #   legend(0, 1.05, legend=c("100", "500", "5000"), title="PopSize", col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15, 15, 15), cex=0.8)
          #  legend(3, 1, legend=c("L=5", "L=15"), pch=dots, bty="n", cex=0.8)
        #}
            
        k=k+1
            
            
        j=1
        i=0.4
        for(N in popSize){
            c1=1
            for(L in c(5, 15)){
                file=paste(args[3], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".popDetail", sep="")
                
                data=read.table(file, comment.char = "", header=TRUE)
                
                tmp1=tapply(data$distance, data$recombination_rate, quantile, p=c(0.05, 0.25, 0.5, 0.75, 0.95))
                boxplot(tmp1, range=100, cex.axis=0.7, boxwex = 0.12, at=seq(i,i+8,2),xlab="", ylab="", xaxt="n", yaxt="n", col=rgb(0,0,0, alpha = 0), add=TRUE, border=rgb(red[j], green[j], blue[j],  alpha=0.6), whisklty = 1, lwd=lwidth[c1])
                lines(seq(i,i+8,2), tapply(data$distance, data$recombination_rate, quantile, p=0.5), col=rgb(red[j], green[j], blue[j]), lwd=1)
                points(seq(i,i+8,2), tapply(data$distance, data$recombination_rate, quantile, p=0.5), col=rgb(red[j], green[j], blue[j]), pch=dots[c1], cex=0.8)
                
                i=i+0.6
                c1=c1+1
                
            }
                i=i-1
                j=j+1
        }
    }
}

dev.off()