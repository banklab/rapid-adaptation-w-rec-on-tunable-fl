args=commandArgs(trailingOnly = TRUE)


popSize=c(100, 500, 5000)
add=c(0.01)
#add=c(0.001)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "initMAF.pdf", sep="."), width=2, height=2)

par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
lwidth=c(0.5, 0.6)
dots=c(17, 19)


alpha=c(0.3, 0.25, 0.5, 0.7, 0.9)

lwidth=c(0.5, 0.6)
ltype=c("11", "solid", "solid", "solid","solid")


k=1

for(L in c(5, 15)){
    for(a in add){
        for(e in epi){
            plot(NULL, xlim=c(1,L), ylim=c(0, 0.5), ylab="", xlab="", xaxt="n")
            title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
            #            mtext("Fixation step", side=1, line=1.5, cex=0.7)
			mtext(paste("Fixation step (L=", L, ")", sep=""), side=1, line=1.4, cex=0.7)
            steps=floor(L/5)
            
            if(k==1 || k==4 || k==7){
                mtext("Mean of initial MAF (fixed)", side=2, line=1.9, cex=0.7, las=0)
            }

            k=k+1
            
            axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))
            
            j=1
            
            for(N in popSize){
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
                data=read.table(file)
                
                ###landscape_id   recombination_rate      conf    alleletype      fixStep fixation_time   initialFreq
                
                l=length(data[1,])
            
                tmpMean=matrix(, ncol=5, nrow=0) ### five recombination rates
                for(i in seq(6, l, 3)){
                    frq=data[, i]
                    frq[frq>0.5]=NA
                 
                    tmpMean=rbind(tmpMean,tapply(frq, data[,2], mean, na.rm=TRUE))
                }
                #print(tmpMean)
                      
                t=1
                for(c in  1:5){
                    lines(1:L, tmpMean[,c], col=rgb(red[j], green[j], blue[j], alpha=alpha[t]), lty=ltype[t])
                    t=t+1
                }
                j=j+1
            }
        }
    }
}
                
dev.off()


for(L in c(5, 15)){
    for(a in add){
        for(e in 0.1){
            for(N in popSize){
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
                data=read.table(file)
                
                ###landscape_id   recombination_rate      conf    alleletype      fixStep fixation_time   initialFreq
                
                l=length(data[1,])
                
                tmpMean=matrix(, ncol=5, nrow=0) ### five recombination rates
                
                frqMajor=matrix(, ncol=2, nrow=0)
                print(c(L, N, a, e))
                
                for(i in seq(6, l, 3)){
                    frq=data[, i]
                    frq[frq<0.5]=1-frq[frq<0.5]
                    frqMajor=rbind(frqMajor, cbind(frq, data[,2]))
                    print(tapply(frq, data[,2], mean, na.rm=TRUE))
                }
                print(tapply(frqMajor[,1], frqMajor[,2], mean, na.rm=TRUE))
                
            }
        }
    }
}

