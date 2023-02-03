
library(wesanderson)

args=commandArgs(trailingOnly = TRUE)

N=5000
add=0.01
epi=c(0.001, 0.01, 0.1)
L=15

pdf(paste(args[1], add, "fixAlleleFrq.mt.pdf", sep="."), width=2, height=2)

par(mfrow=c(1,1), mai=c(0.4, 0.5, 0.2, 0.02), mgp=c(2,0.5,0), las=1, cex =0.8)

lwidth=c(0.5, 0.6)
dots=c(17, 19)

alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)

lwidth=c(0.5, 0.6)
ltype=c("11", "solid", "solid", "solid","solid")



colors=wes_palette("Zissou1", 10, type="continuous")[c(1,4,9)]



a=add

#### plot mean allele frequency
plot(NULL, xlim=c(1,L), ylim=c(0, 1), ylab="", xlab="", xaxt="n")
#title(main="Rouggeness=10", cex.main=1, line=0.2)
mtext("Fixation step", side=1, line=1.5, cex=0.8)
steps=floor(L/5)
mtext("Initila freq. of fixed alleles", side=2, line=1.9, cex=0.8, las=0)
axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))

#legend("bottomleft", legend=c(rep(c(0.001, 0.01, 0.1, 0.5), 1)), col=rgb(red=0, green = 0, blue=0, alpha=alpha[c(2,3,4, 5)]), 
 #      seg.len=0.7, lwd=1.5,lty = rep(1,4), bty="n", title = "Recomb.", cex=0.7)
w=1
for(e in epi){
    rgbValues=col2rgb(colors[w])/255
    w=w+1
  
    file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
    data=read.table(file)
    ###landscape_id   recombination_rate      conf    alleletype      fixStep fixation_time   initialFreq
    l=length(data[1,])
    tmpMeanAF=matrix(, ncol=5, nrow=0) #### mean allele frequency

    for(i in seq(6, l, 3)){
        frq=data[, i]
        tmpMeanAF=rbind(tmpMeanAF,tapply(frq, data[,2], mean, na.rm=TRUE))
        }

t=1
for(c in  1:5){
    lines(1:L, tmpMeanAF[,c], col=rgb(rgbValues[1],  rgbValues[2], rgbValues[3], alpha=alpha[t]), lwd=1.5, lty=ltype[t])
    t=t+1
}               
            
}


#
#tmpPro=rbind(tmpPro,tapply(data[,i-1], data[,2], mean, na.rm=TRUE))


#
#
#### plot proportion of minor alleles
plot(NULL, xlim=c(1,L), ylim=c(0, 1), ylab="", xlab="", xaxt="n")
#title(main="Rouggeness=10", cex.main=1, line=0.2)
mtext("Fixation step", side=1, line=1.5, cex=0.8)
steps=floor(L/5)
mtext("Prop. of fixed minor alleles", side=2, line=1.9, cex=0.8, las=0)
axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))

w=1
for(e in epi){
    rgbValues=col2rgb(colors[w])/255
    w=w+1
    
    file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
    data=read.table(file)
    ###landscape_id   recombination_rate      conf    alleletype      fixStep fixation_time   initialFreq
    l=length(data[1,])

    tmpPro=matrix(, ncol=5, nrow=0) ### proportion of minor alelles
    for(i in seq(5, l, 3)){
        tmpPro=rbind(tmpPro,tapply(data[,i], data[,2], mean, na.rm=TRUE))
    }
    
    t=1
    for(c in  1:5){
        lines(1:L, tmpPro[,c], col=rgb(rgbValues[1],  rgbValues[2], rgbValues[3], alpha=alpha[t]), lwd=1.5, lty=ltype[t])
        t=t+1
    }               
    
}


#### plot MAF
plot(NULL, xlim=c(1,L), ylim=c(0, 0.5), ylab="", xlab="", xaxt="n")
#title(main="Rouggeness=10", cex.main=1, line=0.2)
mtext("Fixation step", side=1, line=1.5, cex=0.8)
steps=floor(L/5)
mtext("Mean of initial MAF (fixed)", side=2, line=1.9, cex=0.8, las=0)
axis(1, at=c(seq(1, L-L*0.2, steps), L), labels = c(seq(1, L-L*0.2, steps), L))

w=1
for(e in epi){
    rgbValues=col2rgb(colors[w])/255
    w=w+1
    print(rgbValues)

    file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
    data=read.table(file)
    ###landscape_id   recombination_rate      conf    alleletype      fixStep fixation_time   initialFreq
    l=length(data[1,])
    
    tmpMeanMAF=matrix(, ncol=5, nrow=0) ### mean minor allele frequency

    for(i in seq(6, l, 3)){
        frq=data[, i]
        frq[frq>0.5]=NA
        tmpMeanMAF=rbind(tmpMeanMAF,tapply(frq, data[,2], mean, na.rm=TRUE))
    }
    
    t=1
    for(c in  1:5){
        lines(1:L, tmpMeanMAF[,c], col=rgb(rgbValues[1],  rgbValues[2], rgbValues[3], alpha=alpha[t]), lwd=1.5, lty=ltype[t])
        t=t+1
    }               
    
}

dev.off()
