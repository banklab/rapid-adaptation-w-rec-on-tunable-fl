args=commandArgs(trailingOnly = TRUE)

popSize=c(5000)
add=c(0.01)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "adaptivePath.ratio.mt.pdf", sep="."), width=6, height=2)

par(mfrow=c(1,3), mai=c(0.4, 0.58, 0.2, 0.02), mgp=c(2,0.5,0), cex = 0.8,  las=1)


lwidth=c(0.5, 0.6)
dots=c(17, 19)


alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)

lwidth=c(0.5, 0.6)
ltype=c("11", "solid", "solid", "solid","solid")


k=1
L=15
for(a in add){
    for(e in epi){
        ylow=0
        yup=2
        stepSize=0.5
        if(a>e){
            yup=0.2
            stepSize=0.05
        }
        plot(NULL, xlim=c(1,L), ylim=c(ylow, yup), ylab="", xlab="", xaxt="n", yaxt="n")
        title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.4, cex.main=0.9)
        mtext("Fixation step", side=1, line=1.4, cex=0.8)
        
        abline(h=1, lty=2, col="grey")
        abline(h=1.414, lty=2, col="grey")

        if(k==1 || k==4){
            mtext(expression(Psi), side=2, line=2.6, las=0, cex=0.8, las=1)
        }
            
            steps=floor(L/5)
            axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))
            if(a>e){
                values=seq(ylow, yup, stepSize)
                axis(2, at=values, labels = c("0.00", 0.05, "0.10",0.15, "0.20"))
                abline(h=0.1, lty=2, col="grey")
            }
            else{
            values=seq(ylow, yup, stepSize)
            axis(2, at=values, labels = c("0.0", 0.5, "1.0",1.5, "2.0"))
            }
            k=k+1
            
            j=1
            
            for(N in popSize){
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".oneLocEffect2final", sep="")
                finalOne=read.table(file)
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".oneLocEffect2major", sep="")
                majorOne=read.table(file)
                
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".epiTraj2final", sep="")
                finalTwo=read.table(file)
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".epiTraj2major", sep="")
                majorTwo=read.table(file)

                finalMeanRatio=matrix(, ncol=5, nrow=0)
                majorMeanRatio=matrix(, ncol=5, nrow=0)
                
                for(i in 1:(L-1)){
                    finalMeanRatio=rbind(finalMeanRatio,tapply(abs(finalTwo[,i+3]), finalTwo[,2], mean)/tapply(abs(finalOne[,i+3]), finalTwo[,2], mean))
                    majorMeanRatio=rbind(majorMeanRatio,tapply(abs(majorTwo[,i+3]), majorTwo[,2], mean)/tapply(abs(majorOne[,i+3]), majorTwo[,2], mean))
                    
                }
                
                t=1
                for(c in  1:5){
                    lines(1:(L-1)+0.5, finalMeanRatio[,c], col=rgb(0,0,0, alpha=alpha[t]), lty=ltype[t])
                    points(1:(L-1)+0.5, finalMeanRatio[,c], col=rgb(0,0,0, alpha=alpha[t]), pch=3, cex=0.4)
                    
                    lines(1:(L-1)+0.5, majorMeanRatio[,c], col=rgb(0,0,0, alpha=alpha[t]), lty=ltype[t])
                    points(1:(L-1)+0.5, majorMeanRatio[,c], col=rgb(0,0,0, alpha=alpha[t]), pch=4, cex=0.4)
                    t=t+1

                }
                j=j+1
            
                

            }
        }
    }
    

dev.off()
