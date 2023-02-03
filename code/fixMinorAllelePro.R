
args=commandArgs(trailingOnly = TRUE)

popSize=c(100, 500, 5000)
add=c(0.01)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "fixMinorAllelePro.pdf", sep="."), width=2, height=2)

### plot fixation allele type

par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
dots=c(17, 19)

alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)
lwidth=c(0.5, 0.6)
ltype=c("11", "solid", "solid", "solid","solid")



k=1

for(L in c(5, 15)){
    for(a in add){
        for(e in epi){
            
             plot(NULL, xlim=c(1,L), ylim=c(0,1), ylab="", xlab="", xaxt="n")
            title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
            #mtext("Fixation step", side=1, line=1.5, cex=0.7)
			 mtext(paste("Fixation step (L=", L, ")", sep=""), side=1, line=1.4, cex=0.7)
             #if(k==1 || k==4 || k==7){
                 mtext("Prop. of fixed minor alleles", side=2, line=1.9, cex=0.7, las=0)
             #}
             
            # if(k==3 ){
             #    legend("topleft", legend=c("100", "500", "5000"), title="PopSize", col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15,15,15), cex=0.8)
              #   legend("topright", legend=c(rep(c(0.001, 0.01, 0.1, 0.5), 1)), col=rgb(red=0, green = 0, blue=0, alpha=alpha[c(2,3,4,5)]), seg.len=0.7, lwd=1.5,lty = rep(1,4), bty="n", title = "Recomb.Prob.", cex=0.8)
             #}
             
             
             k=k+1
             
             steps=floor(L/5)
            
             
             axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))
             
             j=1
        
             for(N in popSize){
                file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
                
                data=read.table(file)
                l=length(data[1,])
                
                tmpMean=matrix(, ncol=5, nrow=0)
                
                for(i in seq(5, l, 3)){
                    tmpMean=rbind(tmpMean,tapply(data[,i], data[,2], mean, na.rm=TRUE))
                }
          
                
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
