
library(wesanderson)

args=commandArgs(trailingOnly = TRUE)

N=5000  # popSize
add=c(0.01)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "oneLocEffect.pdf", sep="."), width=6, height=2)

par(mfrow=c(1,3), mai=c(0.4, 0.58, 0.2, 0.02), mgp=c(2,0.5,0), cex =0.8,  las=1)

alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)

lwidth=c(0.5, 0.6)

ltype=c("11", "solid", "solid", "solid","solid")


dots=c(4,3)



L=15

k=1

for(a in add){
    for(e in epi){
        if(a>e){
            ylow=-a
            yup=a
        }
        else{
        ylow=-e
        yup=3*e
        }
        plot(NULL, xlim=c(1,L), ylim=c(ylow*1.2, yup*1.05), ylab="", xlab="", xaxt="n", yaxt="n")
        title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.4, cex.main=0.9)
        mtext("Fixation step", side=1, line=1.4, cex=0.8)
       #text(1, yup[k]-0.05*(yup[k]-ylow[k]), labels=paste("L=", L, sep=""), adj=0, cex=1)
        abline(h=e, lty=2, col="grey")
        abline(h=-e, lty=2, col="grey")
        steps=floor(L/5)
        print(steps)
        print(c(seq(1, L-L*0.2, steps), L))
        axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))
            
        if(a>e){
            values=seq(ylow, yup, by=a)
        }else{
            values=seq(ylow, yup, e)}
        axis(2, at=values, labels = values)
            
        if(k==1){
            mtext("One-locus effect", side=2, line=2.6, las=0, cex=0.8)
        }

        legend("topleft", legend=substitute(paste("Epi. ", sigma, "=", e, sep=""), list(e=e)),  cex=0.8, bty="n", adj=0.2)
            k=k+1
            j=1
            
            for(bg in c("major", "final")){
                 file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".oneLocEffect2", bg, sep="")
                data=read.table(file)

                numCol=length(data[1,])
                
                stepEpi=t(apply(data[, 4:numCol], 1, function(x) x))

                tmpMean=matrix(, ncol=5, nrow=0)
                #tmpSD=matrix(, ncol=4, nrow=0)
                for(i in 1:(L-1)){
                    tmpMean=rbind(tmpMean,tapply(stepEpi[,i], data[,2], mean))
                   # tmpSD=rbind(tmpSD,tapply(stepEpi[,i], data[,2], sd))
                }
                
                t=1
                for(c in  1:5){
                    #rgbValues=col2rgb(colors[j])/255
                    lines(1:(L-1), tmpMean[,c], col=rgb(red=0, green=0, blue=0, alpha=alpha[t]), lty=ltype[t])
                    points(1:(L-1), tmpMean[,c], col=rgb(red=0, green=0, blue=0, alpha=alpha[t]), pch=dots[j], cex=0.4)
                    t=t+1
                #    for(step in 1:(L-1)){
                 #       lines(rep(step+N*L/300000+c*L/200-L/50,2)+0.5, c(tmpMean[step,c]-tmpSD[step,c], tmpMean[step,c]+tmpSD[step,c]), col=rgb(red[j], green[j], blue[j], alpha=alpha[j]), lwd=0.3)
                  #  }
                }
                j=j+1
        }
    }
}


dev.off()
