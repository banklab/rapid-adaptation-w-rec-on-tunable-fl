

library(vioplot)
library(wesanderson)

args=commandArgs(trailingOnly = TRUE)

N=5000
add=c(0.01)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "intFixAlleleFrq.mt.pdf", sep="."), width=2, height=2)
par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)

colors=wes_palette("Zissou1", 10, type="continuous")[c(1,4,9)]

for(L in c(5, 15)){
    boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(0,1))
    title(main=paste("L=", L, sep=""), line=0.2, cex.main=0.8)
    mtext("Recombination probability", side=1, line=1.5, cex=0.7)
    mtext("Initial freq. of fixed alleles", side=2, line=1.9, cex=0.7, las=0)
    axis(side=1, label=c(0, 0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85, 2), cex.axis=0.75)
    
    for(a in add){
    	c=1
		i=0.5
        for(e in epi){
            file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
            data=read.table(file)
           
			expandData=data.frame(matrix(, ncol=2, nrow=0))
            colnames(expandData)=c("recomb", "frq")
            for(frq in seq(6, length(data[1,]), 3)){
				extract=data[, c(2, frq)]
				colnames(extract)=c("recomb", "frq")
               expandData=rbind(expandData, extract)
            }
            
			tmp=tapply(expandData[,2], expandData[,1], quantile, p=c(0.05, 0.25, 0.5, 0.75, 0.95))
            boxplot(tmp, range=100, cex.axis=0.7, boxwex = 0.12, at=seq(i,i+8,2),
                 xlab="", ylab="", xaxt="n", yaxt="n", add=TRUE, col="white", border=colors[c], whisklty = 1)

#			vioplot(frq~recomb, data=expandData, at=seq(i,i+6,2), add=TRUE, wex = 0.5, xlab="", ylab="", xaxt="n", yaxt="n", col="white", border=colors[c], colMed = colors[c], lineCol = colors[c], rectCol = colors[c], lwd=0.4)

            i=i+0.4
            c=c+1
		}

	}

    for(a in add){
    	c=1
		i=0.5
        for(e in epi){
            file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
            data=read.table(file)

			expandData=data.frame(matrix(, ncol=2, nrow=0))
            colnames(expandData)=c("recomb", "frq")
            for(frq in seq(6, length(data[1,]), 3)){
				extract=data[, c(2, frq)]
				colnames(extract)=c("recomb", "frq")
                expandData=rbind(expandData, extract)
            }
           
           lines(seq(i,i+8,2), tapply(expandData[,2], expandData[,1], quantile, p=0.5), col=colors[c])
            i=i+0.4
            c=c+1
        }
    }
}

dev.off()


for(L in c(5, 15)){
    for(a in 0.01){
        for(e in epi){
            file=paste(args[2], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".fixAlleleInitFreq.ordered", sep="")
            data=read.table(file)
            
            expandData=data.frame(matrix(, ncol=2, nrow=0))
            colnames(expandData)=c("recomb", "frq")
            for(frq in seq(6, length(data[1,]), 3)){
                extract=data[, c(2, frq)]
                colnames(extract)=c("recomb", "frq")
                expandData=rbind(expandData, extract)
            }
            
            
            values=expandData[,2]
            groups=expandData[,1]
            
            print(tapply(expandData[,2], expandData[,1], mean))
            
            #print(kruskal.test(values ~ groups) )
            #print(pairwise.wilcox.test(values, groups,p.adjust.method = "BH"))
        }
    }
}

