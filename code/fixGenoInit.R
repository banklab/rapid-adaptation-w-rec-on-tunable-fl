args=commandArgs(trailingOnly = TRUE)


popSize=c(100, 500,  5000)
add=c(0.01)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "finalGenoIntFrq.pdf", sep="."), width=2, height=2)
par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
lwidth=c(0.5, 0.6)
dots=c(17, 19)

k=1
for(a in add){
    for(e in epi){
        boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(0,1))
        title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
        mtext("Recombination probability", side=1, line=1.5, cex=0.7)
        #if(k==1 || k==4 || k==7){
        mtext("Pro. of present finalGeno.", side=2, line=2, cex=0.8, las=0)
        #}
        axis(side=1, label=c(0, 0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85,2), cex.axis=0.75)
        timePoint=c(0.01, 0.03, 0.1, 0.2, 0.5, 1, 2, 5)

        #if(k==3){
        #    legend(0,0.5, legend=c("100", "500", "5000"), title="PopSize", col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15,15,15), cex=0.8)
         #   legend(3, 0.3, legend=c("L=5", "L=15"), pch=dots, bty="n", cex=0.8)
        #}
        
        k=k+1
        
        j=1
        for(N in popSize){
            
           
            c1=1
            d=1
            for(L in c(5, 15)){
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".finalGenoIntFrq", sep="")
                data=read.table(file)
                
                tmp1=tapply(data$V5, data$V2, function(x){length(x[x>0])/length(x)})
                points(seq(0.85, 8.85,2), tmp1, col=rgb(red[j], green[j], blue[j]), lwd=1, type="p", pch=dots[d], cex=0.8)
                points(seq(0.85, 8.85,2), tmp1, col=rgb(red[j], green[j], blue[j]), lwd=1, type="l")
                
                c1=c1+1
                d=d+1
            }
            
            j=j+1
        }
    }
}

k=1
for(a in add){
    for(e in epi){
        boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(-4,0), yaxt="n")
        title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
        mtext("Recombination probability", side=1, line=1.5, cex=0.7)
        #if(k==1 || k==4 || k==7){
        mtext("Init. Frq of finalGeno.", side=2, line=2, cex=0.8, las=0)
        #}
        axis(side=1, label=c(0,0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85,2), cex.axis=0.75)
        timePoint=c(0.01, 0.03, 0.1, 0.2, 0.5, 1, 2, 5)
        axis(side=2, at=log10(c(0.0001, 0.001, 0.01, 0.1, 1)), label=c(0.0001, 0.001, 0.01, 0.1, 1), cex.axis=0.85)
        
        
        j=1
        i=0.4
        for(N in popSize){
            k=k+1
            c1=1
            d=1
            for(L in c(5, 15)){
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".finalGenoIntFrq", sep="")
                data=read.table(file)
                frq=data$V5
                frq[frq<1]=NA
                tmp1=tapply(log10(frq/N), data$V2, quantile, p=c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
                boxplot(tmp1, range=0, cex.axis=0.7, boxwex = 0.12, at=seq(i,i+8,2),xlab="", ylab="", xaxt="n", yaxt="n",  add=TRUE, col=rgb(0,0,0, alpha = 0), border=rgb(red[j], green[j], blue[j],  alpha=0.6), whisklty = 1, lwd=lwidth[c1])
                points(seq(i,i+8,2), tapply(log10(frq/N), data$V2, quantile, p=0.5, na.rm=TRUE), col=rgb(red[j], green[j], blue[j]), lwd=1, type="l")
                points(seq(i,i+8,2), tapply(log10(frq/N), data$V2, quantile, p=0.5, na.rm=TRUE), col=rgb(red[j], green[j], blue[j]), lwd=1, type="p",pch=dots[d], cex=0.8)
                i=i+0.6
                c1=c1+1
                d=d+1
            }
            i=i-1
            j=j+1
        }
    }
}





dev.off()
