#!/usr/bin/env python
"""
===============
Bacteria movement
===============

"""
#    Copyright (C) 2006-2017
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

import math
import time
from enum import Enum

import matplotlib.pyplot as plt

import networkx as nx

from functools import reduce
from operator import or_

def merge(*dicts):
    return 
    { k: reduce(lambda d, x: x.get(k, d), dicts, None) for k in reduce(or_, map(lambda x: x.keys(), dicts), set()) }

BACT_COUNT = 8
AI_COUNT = 10
AIC_HCOUNT = 4
AIC_WCOUNT = 4

AIC_RFRAC = 1.0 / AIC_HCOUNT
AIC_CFRAC = 1.0 / AIC_WCOUNT

AI_TYPE_COUNT = 3
BACT_TYPE_COUNT = 2
class AIType:
  BACT_1, BACT_2, TCELL = range(AI_TYPE_COUNT)

class BactType:
  BACT_1, BACT_2 = range(BACT_TYPE_COUNT)



def main():
  bact_colors = ['blue', 'green']
  ai_colors = ['red', 'orange', 'yellow']
  bact_speed = [0.1, 0.05]
  bact_map = {BactType.BACT_1: AIType.BACT_1, 
              BactType.BACT_2: AIType.BACT_2}

  # initial data structures
  bactGraphs = []
  aiGraphs = []
  bactPosInfoOrig = []
  bactPosInfo = []
  for i in xrange(BACT_TYPE_COUNT):
    b = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    bactPosInfoOrig += [b]
  bactPosInfo = bactPosInfoOrig
  
  aiPosInfo = []
  for i in xrange(AI_TYPE_COUNT):
    a = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    aiPosInfo += [a]

  for b in xrange(BACT_TYPE_COUNT):
    bg = nx.empty_graph(BACT_COUNT)
    bactGraphs += [bg]

    initialBgPos = nx.random_layout(bg)
    
    for node in initialBgPos.keys():
      posArray = initialBgPos[node]
      x = posArray[0]
      y = posArray[1]
      c = int(x / AIC_CFRAC)
      r = int(y / AIC_RFRAC)
      bactPosInfo[b][r][c][node] = posArray

  for a in xrange(AI_TYPE_COUNT):
    ai = nx.empty_graph(AI_COUNT)
    aiGraphs += [ai]
    initialAiPos = nx.random_layout(ai)
    for node in initialAiPos.keys():
      posArray = initialAiPos[node]
      x = posArray[0]
      y = posArray[1]
      c = int(x / AIC_CFRAC)
      r = int(y / AIC_RFRAC)
      aiPosInfo[a][r][c][node] = posArray

  # originally we want random distribution of AI & Bacteria


  try:
    import pygraphviz
    from networkx.drawing.nx_agraph import graphviz_layout
    layout = graphviz_layout
  except ImportError:
    try:
        import pydot
        from networkx.drawing.nx_pydot import graphviz_layout
        layout = graphviz_layout
    except ImportError:
        print("PyGraphviz and pydot not found;\n"
              "drawing with spring layout;\n"
              "will be slow.")
        layout = nx.spring_layout

  first = True
  while(1):
    time.sleep(1)
    plt.cla()

    region = 110  # for pylab 2x2 subplot layout
    plt.subplots_adjust(left=0, right=1, bottom=0, top=0.95, wspace=0.01, hspace=0.01)
        
    region += 1
    plt.subplot(region)
    plt.title("p = %6.3f" % (23))

    for i in xrange(BACT_TYPE_COUNT):
      color_map = [bact_colors[i]] * BACT_COUNT
      posDicts = bactPosInfo[i]
      #print(posDicts)
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])
      
      nx.draw(bactGraphs[i], all_coords, 
        node_color=color_map, with_labels=False, node_size=25)

    for i in xrange(AI_TYPE_COUNT):
      color_map = [ai_colors[i]] * AI_COUNT
      posDicts = aiPosInfo[i]
      #print(posDicts)
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])
      
      #print(all_coords)

      nx.draw(aiGraphs[i], all_coords, 
        node_color=color_map, with_labels=False, node_size=10)


    # we have information about the previous positions.
    # so, now go through each and find the quadrant with max number
    # of ais of each type
    maxAIQuadrants = []
    for a in xrange(AI_TYPE_COUNT):
      posDicts = aiPosInfo[i]
      maxLen = 0
      (maxR, maxC) = (-1, -1)
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          if (len(posDicts[r][c]) > maxLen):
            maxLen = len(posDicts[r][c])
            maxR = r
            maxC = c
      maxAIQuadrants += [(r,c)]

    for b in xrange(BACT_TYPE_COUNT):
      ai = bact_map[b]
      (maxr, maxc) = maxAIQuadrants[ai]
      # find the mid point of this quadrant
      (x, y) = ((maxc + 0.5) * AIC_CFRAC,
                (maxr + 0.5) * AIC_RFRAC)
      posDicts = bactPosInfo[b]
      #print(posDicts)
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])

      bactPosInfo = bactPosInfoOrig
      for node in all_coords.keys():
        coords = all_coords[node]
        (dx, dy) = (x - coords[0], y - coords[1])
        (newx, newy) = (coords[0] + dx*bact_speed[b], 
                        coords[1] + dy*bact_speed[b])
        all_coords[node][0] = newx;
        all_coords[node][1] = newy;

        posArray = all_coords[node]
        x = posArray[0]
        y = posArray[1]
        c = int(x / AIC_CFRAC)
        r = int(y / AIC_RFRAC)
        bactPosInfo[b][r][c][node] = posArray

    plt.draw()

    plt.pause(1e-17)
    
  #global plot
  #plt.plot.show()

if __name__ == "__main__":
    main()
