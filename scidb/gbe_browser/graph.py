from __future__ import print_function
from __future__ import division
from random import shuffle

# Written by Manuel A. Rivas
# Updated 04.21.2017

import sys,re
import os
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()


def graph(**data):
    print(data.keys())
    print("HERE")
    print(data['thetainvfindict'])
    key = data['key']
    thetainvfindict = data['thetainvfindict']
    clustvalue = data['clustvalue']
    iter = data['iter']
    x = 4
    fn = 'static/images/' + str(key) + '.svg'
    ro.r('''
      graphfn <- function(fn, precmatrix){
        require(GGally)
       require(network)
       require(sna)
       require(ggplot2)
       require(corpcor)
       require(RColorBrewer)

         net <- cor2pcor(cov2cor(solve(precmatrix)))
      print(net)
     net <- as.network.matrix(net)
      print(net)
      letters <- seq(1,10,1)
      network.vertex.names(net) <- c("CD","UC")
        p <- ggnet2(net, color = network.vertex.names(net), palette = "Set2")
       print("HI")
       ggsave(fn, plot = p, width = 6.5, height = 5.5)
}
       ''')
    graphf = ro.globalenv['graphfn']
    graphf = ro.r['graphfn']
    precmatrix = thetainvfindict[iter-1,2]
    res = graphf(fn, precmatrix)
