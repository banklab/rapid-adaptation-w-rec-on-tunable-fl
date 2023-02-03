library(wesanderson)

pdf("legend.pdf", width=6, height=0.4)

par(mai=rep(0,4), cex=1.5)

plot.new()
colors=wes_palette("Zissou1", 10, type="continuous")[c(1,4,9)]
text(0.4, 1.3,  substitute(bold("Ruggedness")), adj=c(1,1.85))
legend(0.42,1.3, legend=c("0.1", "1", "10"), pch=c(19,19,19), col=colors, bty = "n", horiz=TRUE)



plot.new()
colors=wes_palette("Zissou1", 10, type="continuous")[c(1,4,9)]
text(0.25, 1.2,  substitute(bold("Ruggedness")), adj=c(1,1.85), cex=0.5)
legend(0.3,1.2, legend=c("0.1", "1", "10"), pch=c(19,19,19), col=colors, bty = "n", horiz=TRUE, cex=0.5)

ltype=c("11", "solid", "solid", "solid","solid")
alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)
text(0.25, 0.7,  substitute(bold("Recombination probability")),  adj=c(1,1.9), cex=0.5)
legend(0.3,0.7, legend=c(0, 0.001, 0.01, 0.1, 0.5), lty=ltype,  col=rgb(red=0, green = 0, blue=0, alpha=alpha), lwd=3, seg.len=1.5, bty = "n", horiz=TRUE, cex=0.5, xjust=0)



plot.new()
text(0.25, 1.2,  substitute(bold("Background genotype")), adj=c(1,1.9), cex=0.5)
legend(0.3,1.2, legend=c("Final", "Major"), pch=c(3, 4), bty = "n", horiz=TRUE, cex=0.5)
ltype=c("11", "solid", "solid", "solid","solid")
alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)
text(0.25, 0.7,  substitute(bold("Recombination probability")),  adj=c(1,1.9), cex=0.5)
legend(0.3,0.7, legend=c(0, 0.001, 0.01, 0.1, 0.5), lty=ltype,  col=rgb(red=0, green = 0, blue=0, alpha=alpha), lwd=3, seg.len=1.5, bty = "n", horiz=TRUE, cex=0.5, xjust=0)

plot.new()

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
text(0.65, 1.2,  substitute(bold("Population size")), adj=c(1,1.85), cex=0.5)
legend(0.67, 1.2, legend=c("100", "500", "5000"),  col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15, 15, 15), horiz=TRUE,cex=0.5)

text(0.25, 1.2,  substitute(bold("Background genotype")), adj=c(1,1.85), cex=0.5)
legend(0.27,1.2, legend=c("Final", "Major"), pch=c(3, 4), bty = "n", horiz=TRUE, cex=0.5)

ltype=c("11", "solid", "solid", "solid","solid")
alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)
text(0.25, 0.7,  substitute(bold("Recombination probability")),  adj=c(1,1.9), cex=0.5)
legend(0.27,0.7, legend=c(0, 0.001, 0.01, 0.1, 0.5), lty=ltype,  col=rgb(red=0, green = 0, blue=0, alpha=alpha), lwd=3, seg.len=1.5, bty = "n", horiz=TRUE, cex=0.5, xjust=0)



plot.new()

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
text(0.3, 1.2,  substitute(bold("Population size")), adj=c(1,1.85), cex=0.5)
legend(0.32, 1.2, legend=c("100", "500", "5000"),  col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15, 15, 15), horiz=TRUE,cex=0.5)

ltype=c("11", "solid", "solid", "solid","solid")
alpha=c(0.5, 0.25, 0.5, 0.7, 0.9)
text(0.3, 0.7,  substitute(bold("Recombination probability")),  adj=c(1,1.9), cex=0.5)
legend(0.32,0.7, legend=c(0, 0.001, 0.01, 0.1, 0.5), lty=ltype,  col=rgb(red=0, green = 0, blue=0, alpha=alpha), lwd=3, seg.len=1.5, bty = "n", horiz=TRUE, cex=0.5, xjust=0)

dev.off()

pdf("legend1.pdf", width=6, height=0.2)

par(mai=rep(0,4), cex=1.5)

plot.new()

red=c(0.2, 0.4, 0.8)
green=rep(0.5, 3)
blue=c(0.8, 0.4, 0.2)
text(0.6, 1.2,  substitute(bold("Population size")), adj=c(1,1.85), cex=0.5)
legend(0.62, 1.2, legend=c("100", "500", "5000"),  col=rgb(red=red, green=green, blue=blue),  bty="n", pch=c(15, 15, 15), horiz=TRUE,cex=0.5)

text(0.2, 1.13,  substitute(bold("Number of Loci")),  adj=c(1,1.85), cex=0.5)
legend(0.22, 1.2, legend=c(5,15), pch=c(17, 19), bty = "n", horiz=TRUE, cex=0.5)

dev.off()