

args=commandArgs(trailingOnly = TRUE)


popSize=c(100, 500, 5000)
add=c(0.01)

epi=c(0.001, 0.01, 0.1)

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
lwidth=c(0.5, 0.6)
dots=c(17, 19)

################ plot fixed genotype position away from global maximum on FL

pdf(paste(args[1], add, "dist2Global.pdf", sep="."), width=2, height=2)
par(mfrow=c(1,1), mai=c(0.4, 0.5, 0.2,0.1), mgp=c(2,0.5,0), cex = 0.7,  las=1)
ylow=c(-0.15, -0.15, -0.32)
k=1

for(a in add){
    for(e in epi){
        boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(ylow[k],0))
        title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
        mtext("Recombination probability", side=1, line=1.5, cex=0.7)
        if(k==1 || k==4 || k==7){
            mtext("Fitness diff. (final-global)", side=2, line=2.6, cex=0.7, las=0)}
        axis(side=1, label=c(0, 0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85,2), cex.axis=0.75)
        k=k+1

        j=1
        i=0.4
        for(N in popSize){
            c1=1
            for(L in c(5, 15)){
                file=paste(args[2], "/data/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, "_.dat", sep="")
                data=read.table(file, comment.char = "", header=TRUE)
                data$tfixation=data$tfixation/N
                tmp1=tapply(-data$genetic_load, data$recombination_rate, quantile, p=c(0.05, 0.25, 0.5, 0.75, 0.95))
                boxplot(tmp1, range=100, cex.axis=0.7, boxwex = 0.12, at=seq(i,i+8,2),xlab="", ylab="", xaxt="n", yaxt="n", add=TRUE, col=rgb(0,0,0, alpha = 0), border=rgb(red[j], green[j], blue[j],  alpha=0.6), whisklty = 1, lwd=lwidth[c1])
                lines(seq(i,i+8,2), tapply(-data$genetic_load, data$recombination_rate, quantile, p=0.5), col=rgb(red[j], green[j], blue[j]), lwd=1)
                points(seq(i,i+8,2), tapply(-data$genetic_load, data$recombination_rate, quantile, p=0.5), col=rgb(red[j], green[j], blue[j]), pch=dots[c1], cex=0.8)
                i=i+0.6
                c1=c1+1
            }
                i=i-1
                j=j+1
        }
    }
}


################ plot mean population fitness difference (initial-final)/initial


pdf(paste(args[1], add, "fitnessChange.pdf", sep="."), width=2, height=2)
par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2, 0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1)
k=1
yup=c(0.15, 0.15, 0.5)
for(a in add){
        for(e in epi){
            boxplot(NULL, xlim=c(0,10), xlab="", ylab="", ylim=c(-0.01,yup[k]))
            title(main=substitute(paste("Ruggedness ", lambda,"=", x, sep=""), list(x=e/a)), line=0.5, cex.main=0.8)
            mtext("Recombination probability", side=1, line=1.5, cex=0.7)
            if(k==1 || k==4 || k==7){
                mtext("Fitness change", side=2, line=2, cex=0.8, las=0)
            }
            axis(side=1, label=c(0, 0.001, 0.01, 0.1, 0.5), at=seq(0.85, 8.85,2), cex.axis=0.75)
            
           # if(k==1){
            #    legend(0, 0.155, legend=c("100", "500", "5000"), title="PopSize", col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15,15,15), cex=0.8)
             #   legend(3, 0.15, legend=c("L=5", "L=15"), pch=dots, bty="n", cex=0.8)
            #}
            
            k=k+1
            i=0.4
            j=1
        for(N in popSize){

                c1=1
            for(L in c(5, 15)){
                file=paste(args[3], "/rmf_L", L, "_N", N, "_mua", a, "_sa0_sb", e, ".popDetail", sep="")
                data=read.table(file, comment.char = "", header=TRUE)
                tmp1=tapply((data$mean_fitness_final-data$mean_fitness_initial)/data$mean_fitness_initial, data$recombination_rate, quantile, p=c(0.05, 0.25, 0.5, 0.75, 0.95))

                boxplot(tmp1, range=100, cex.axis=0.7, boxwex = 0.12, at=seq(i,i+8,2),xlab="", ylab="", xaxt="n", yaxt="n", add=TRUE, col=rgb(0,0,0, alpha = 0), border=rgb(red[j], green[j], blue[j], alpha=0.6), whisklty = 1, lwd=lwidth[c1])
                points(seq(i,i+8,2), tapply((data$mean_fitness_final-data$mean_fitness_initial)/data$mean_fitness_initial, data$recombination_rate, quantile, p=0.5), col=rgb(red[j], green[j], blue[j]), pch=dots[c1], cex=0.8)
                lines(seq(i,i+8,2), tapply((data$mean_fitness_final-data$mean_fitness_initial)/data$mean_fitness_initial, data$recombination_rate, quantile, p=0.5), col=rgb(red[j], green[j], blue[j]))
                
                i=i+0.6
                c1=c1+1
            }
                i=i-1
                j=j+1
        }
    }
}


dev.off()