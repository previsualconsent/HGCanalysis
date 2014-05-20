#!/usr/bin/env python

import ROOT

"""
generate in-time pu sets
"""
def generateInTimePUSets(totalEvts,avgPU):

    #check inputs
    evRanges=[]
    if avgPU<=0 : return evRanges
    if avgPU>totalEvts:
        print '[Warning] requested %d <PU> is larger than the total number of events to mix: %d'%(avgPU,totalEvts)

    #generate the event sets with <avgPU>
    curEv=1
    while True:
        iPu=ROOT.gRandom.Poisson(avgPU)
        evRanges.append([curEv,curEv+iPu])
        curEv=curEv+iPu+1
        if curEv>=totalEvts : break
    return evRanges
