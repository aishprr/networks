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

  # initial data structures
  bactGraphs = []
  aiGraphs = []
  bactPosInfo = []
  for i in xrange(BACT_TYPE_COUNT):
    b = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    bactPosInfo += [b]
  
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


  while(1):
    time.sleep(0.1)
    
    region = 110  # for pylab 2x2 subplot layout
    plt.subplots_adjust(left=0, right=1, bottom=0, top=0.95, wspace=0.01, hspace=0.01)
        
    print "*********************"
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
      
      #print(all_coords)

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

    
    # identify largest connected component
    # Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
    # G0 = Gcc[0]
    # nx.draw_networkx_edges(G0, pos,
    #                        with_labels=False,
    #                        edge_color='r',
    #                        width=6.0
    #                       )
    # # show other connected components
    # for Gi in Gcc[1:]:
    #     if len(Gi) > 1:
    #         nx.draw_networkx_edges(Gi, pos,
    #                                with_labels=False,
    #                                edge_color='r',
    #                                alpha=0.3,
    #                                width=5.0
    #                                 )
    plt.draw()
    plt.pause(1e-17)
    #while(True):
    #  pass

    #break

  #global plot
  #plt.plot.show()

if __name__ == "__main__":
    main()
