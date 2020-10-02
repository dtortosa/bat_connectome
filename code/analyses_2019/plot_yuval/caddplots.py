#!/usr/bin/python

import os, sys
import string
from pandas import *
from scipy import stats
from scipy.stats import gaussian_kde

import rpy2.robjects as robjects
import pandas.rpy.common as com
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import grid
from rpy2.robjects.lib import ggplot2

r = robjects.r
grdevices = importr('grDevices')

#--------------------------------------------------------------------#
#                               Globals                              #
#--------------------------------------------------------------------#
mytheme = {
            'panel.background':ggplot2.element_rect(fill='white',colour='white'),
            'axis.text':ggplot2.element_text(colour="black",size=15),
            'axis.line':ggplot2.ggplot2.element_line(size = 1.2, colour="black"),
            'axis.title':ggplot2.element_text(colour="black",size=15),
            'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            'panel.grid.minor':ggplot2.element_blank(),
            'panel.grid.major':ggplot2.element_blank(),
            'legend.key':ggplot2.element_blank(),
            'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15),
            'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15)
            }


def fixRLevels( r_dataframe, column, tlevels ):
    replace = robjects.FactorVector( r_dataframe.rx2(column), 
                                    levels=robjects.StrVector(tlevels) )
    allcolumns = list(r_dataframe.colnames)
    allcolumns = [ x+"_old" if x == column else x for x in allcolumns]
    new_r_df = r_dataframe.cbind(replace)
    new_r_df.colnames = robjects.StrVector(allcolumns+[column])
    return new_r_df
# END fixRLevel

def randomSamplings( rawdatalist, samplesize, nboots=1000, maxx=20 ) :
    assert len(rawdatalist) > 0, "list should be non-empty"
    # Prune empty values
    rawdatalist = [x for x in rawdatalist if x is not None]
    samplings = []
    ind = np.linspace(min(rawdatalist), maxx,512)
    kde = gaussian_kde( rawdatalist )
    for boot in range(0,nboots) :
        kdesub = gaussian_kde( kde.resample(len(resdf)) )
        kdedf = DataFrame( {"subsetname":"Random%d" % boot,
                          "Density":kdesub.evaluate(ind), "value":ind} )
        samplings.append( kdedf )

    samplings = concat( samplings ).reset_index(drop=True)
    return samplings
# END randomSamplings

################################################################################
# Main
################################################################################
if __name__ == "__main__" :

    caddscoresfile = "new_p.txt"
    caddscores = [x.rstrip() for x in open(caddscoresfile).readlines()]
    caddf = DataFrame( {"cadd":caddscores,"dtype":"Expected p-values"} )
    caddf = caddf[caddf["cadd"].notnull()]
    caddf["cadd"] = caddf["cadd"].astype(float)

    resscoresfile = "predicted_p.txt"
    resscores = [x.rstrip() for x in open(resscoresfile).readlines()]
    resdf = DataFrame( {"cadd":resscores,"dtype":"Observed p-values"} )
    resdf = resdf[resdf["cadd"].notnull()]
    resdf["cadd"] = resdf["cadd"].astype(float)

    alldf = concat( [resdf, caddf] )

    maxx = resdf["cadd"].max()
    # Generate Random Sample Densities
    samplings = randomSamplings( caddf["cadd"].tolist(), len(resdf), nboots=1000, maxx=maxx )

    # Calculate Quantiles
    quants = samplings.groupby(["value"])["Density"].quantile([.025,.5,.975]).reset_index()
    quants.rename(columns={'level_1':"Quantile",0:"Density"}, inplace=True )
    quants["linetype"] = ["Mean" if x == .5 else "95% threshold" for x in quants.Quantile]

    # Calculate Case KDE
    # the loop is a bit unnecessary, but it was easier to just leave it 
    # in case more classes of data were produced
    casedf = []
    ind = np.linspace(caddf["cadd"].min(), maxx,512)
    for dtype,lofclass in alldf.groupby("dtype") :
        kde = gaussian_kde( alldf[alldf.dtype == dtype]["cadd"].tolist() )
        tmpdf = DataFrame({"dtype":dtype, "Density":kde.evaluate(ind), "cadd":ind})
        casedf.append( tmpdf )

    casedf = concat( casedf ).reset_index(drop=True)

    rsamplings = com.convert_to_r_dataframe(samplings)
    rlofgenes = com.convert_to_r_dataframe(casedf)
    rquants = com.convert_to_r_dataframe(quants)
    rquants = fixRLevels( rquants,"linetype", ["Mean","95% threshold"] )
    p = (ggplot2.ggplot(rlofgenes) +
                ggplot2.aes_string(x="cadd",y="Density") + 
                ggplot2.geom_line( ggplot2.aes_string(x="value", y="Density", 
                                                      group="factor(subsetname)"),
                                  color="grey", data=rsamplings ) +
                ggplot2.geom_line( ggplot2.aes_string(x="value",y="Density",
                                                      linetype="factor(linetype)", 
                                                      group="factor(Quantile)"),
                                  color="black", data=rquants ) +
                ggplot2.geom_line( ggplot2.aes_string(color="factor(dtype)") ) +
                ggplot2.scale_y_continuous("Density") +
                ggplot2.scale_x_continuous("P-value") +
                ggplot2.scale_linetype("Confidence Interval") +
                ggplot2.scale_colour_brewer("Variant Type",palette="Set1") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.theme(**mytheme) )
    figname = "cadd_bootstrap.pdf"
    print("Writing file:",figname)
    grdevices.pdf(figname, width=6, height=4)
    p.plot()
    grdevices.dev_off()

    figname = "cadd_bootstrap.png"
    grdevices.png(figname, width=6, height=4, units="in", res=300)
    p.plot()
    grdevices.dev_off()


# END MAIN

