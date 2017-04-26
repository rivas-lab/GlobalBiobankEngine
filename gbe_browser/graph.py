from __future__ import print_function
from __future__ import division
from random import shuffle

# Written by Manuel A. Rivas
# Updated 04.21.2017

import sys,re
import os
import rpy2.robjects as ro

def graph(key):
    print(key)
    x = 10
    fn = 'static/images/' + str(key) + '.svg'
    ro.r('''
       require(GGally)
       require(network)
       require(sna)
       require(ggplot2)
       require(corpcor)
      diag({})*.2*.2 -> b
        b[1,2] <- .2*.2*.9
        b[2,1] <- .2*.2*.9
         net <- cor2pcor(b)
      print(net)
print(b)
       p <- ggnet2(net)
       print("HI")
       ggsave("{}", plot = p, width = 6.5, height = 5.5)
       '''.format(x, fn));

