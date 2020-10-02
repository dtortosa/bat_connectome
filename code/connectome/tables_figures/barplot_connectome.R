#### code ofr plotting bar plot with different seeds (different set of random genes)
seeds = c(2,3,4,4675)

for(i in 1:length(seeds)){

    #select seed
    seed = seeds[i]

    ####load data of all the pairs across genome
    load(paste("/Users/diegosalazar/Google Drive/science/open_projects/connectome/results/connectome_results/rdata/test3_ucp_figure_", seed, ".RData", sep=""))
    
    ###Creamos una matriz a #partir de los datos.
    xtabs.data<-xtabs(proportions ~factor(gene_class)+factor(distance_class), data=data)
    xtabs.data
    
    ###plot
    cairo_pdf(paste("/Users/diegosalazar/Google Drive/science/open_projects/connectome/results/connectome_results/figures/barplot_proportions_", seed, ".pdf", sep=""))#we use cairo     because pdf() gives problems with \u2265 (simbol >=)
    
    #barplot
    xs<-barplot(xtabs.data, beside=TRUE, xpd=F, axes=FALSE, axisnames=F, ylim=c(0,1), width=3, lwd=1:2, col=c("gray", "red"))
    
    #axis texts
    mtext(text=expression(bold("")), side = 3, line=1.8, outer=FALSE, cex=0.8)
    mtext(text=expression(bold("Proportion")), side = 2, line=2.8, outer=FALSE, cex=0.8)
    mtext(text=expression(bold("Biological distance")), side = 1, line=3, outer=FALSE, cex=0.8)
    
    #axis
    axis(side=1, at=c(-10, xs[1,1],xs[1,2], xs[1,3], xs[1,4], xs[1,5], 42) ,labels=c("","", "","", "", "", ""), cex.lab=0.9, cex.axis=1, lwd=1, lwd.ticks=0, font=2, mgp=c(3, 1, 0))
    axis(side=2, cex.lab=1.6, cex.axis=1, font=2, mgp=c(3, 1, 0))
    
    #axis labels of X axis
    mtext(at=(xs[1,1]+xs[2,1])/2, text=expression(bold("<10")), side = 1, line=1, outer=FALSE, cex=0.8)
    mtext(at=(xs[1,2]+xs[2,2])/2, text=expression(bold("<20")), side = 1, line=1, outer=FALSE, cex=0.8)
    mtext(at=(xs[1,3]+xs[2,3])/2, text=expression(bold("<30")), side = 1, line=1, outer=FALSE, cex=0.8)
    mtext(at=(xs[1,4]+xs[2,4])/2, text=expression(bold("<40")), side = 1, line=1, outer=FALSE, cex=0.8)
    mtext(at=(xs[1,5]+xs[2,5])/2, text=expression(bold("\u2265"~"40")), side = 1, line=1, outer=FALSE, cex=0.8)
    
    #legend
    legend(x=xs[2,4], y=1, legend=c(expression(bold("All genes")), expression(bold("BAT genes"))), cex=0.9, xpd=T, fill=c("gray", "red"))
    dev.off()
    
    #clean the workspace (in the next iteration you will load results of other iteration with objetcs with the same name)
    rm(data)
}