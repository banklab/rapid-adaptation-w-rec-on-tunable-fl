library(wesanderson)

args=commandArgs(trailingOnly = TRUE)

N=5000 # population size
add=c(0.01)
epi=c(0.001, 0.01, 0.1)

pdf(paste(args[1], add, "distance2initial.mt.pdf", sep="."), width=2, height=2)

############# plot fixation time
par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)

colors=wes_palette("Zissou1", 10, type="continuous")[c(1,4,9)]

for(L in c(5, 15)){
    boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(0,1))
    title(main=paste("L=", L, sep=""), line=0.2, cex.main=0.8)
    mtext("Recombination probability", side=1, line=1.5, cex=0.7)
    mtext("Diverg. from initial pop.", side=2, line=1.8, cex=0.7, las=0)
    axis(side=1, label=c(0, 0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85,2), cex.axis=0.75)
    
    #if(L==15){
    #    legend("topleft", legend=c("0.1", "1", "10"), title = "Ruggedness", title.adj = 0.1, pch=c(19,19,19), col=colors, bty = "n", horiz=TRUE)
    #}
    
    i=0.5
    c=1
    for(a in add){
        for(e in epi){
            file=paste(args[3], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".popDetail", sep="")
            data=read.table(file, comment.char = "", header=TRUE)
                
            tmp1=tapply(data$distance, data$recombination_rate, quantile, p=c(0.05, 0.25, 0.5, 0.75, 0.95))
            boxplot(tmp1, range=100, cex.axis=0.7, boxwex = 0.12, at=seq(i,i+8,2),
                    xlab="", ylab="", xaxt="n", yaxt="n", add=TRUE, col="white", border=colors[c], whisklty = 1)
            lines(seq(i,i+8,2), tapply(data$distance, data$recombination_rate, quantile, p=0.5), col=colors[c])
            i=i+0.4
            c=c+1
        }
    }
}

dev.off()


for(L in c(5, 15)){
    for(a in 0.01){
        for(e in 0.001){
            file=paste(args[3], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".popDetail", sep="")
            data=read.table(file, comment.char = "", header=TRUE)
            
            values=data$distance
            groups=data$recombination_rate
            
            
            print(kruskal.test(values ~ groups) )
            print(pairwise.wilcox.test(values, groups,p.adjust.method = "BH"))
        }
    }
}
