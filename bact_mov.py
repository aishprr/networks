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
import random

import numpy

import matplotlib.pyplot as plt

import networkx as nx

from functools import reduce
from operator import or_

def merge(*dicts):
    return 
    { k: reduce(lambda d, x: x.get(k, d), dicts, None) for k in reduce(or_, map(lambda x: x.keys(), dicts), set()) }

BACT_COUNT = 3
AI_COUNT = 30
AIC_HCOUNT = 4
AIC_WCOUNT = 4

AIC_RFRAC = 1.0 / AIC_HCOUNT
AIC_CFRAC = 1.0 / AIC_WCOUNT

AI_TYPE_COUNT = 1
BACT_TYPE_COUNT = 1
class AIType:
  BACT_1 = 0#, BACT_2, TCELL = range(AI_TYPE_COUNT)

class BactType:
  BACT_1 = 0#, BACT_2 = range(BACT_TYPE_COUNT)



def main():
  bact_colors = ['blue', 'green']
  bact_dis_thresh = [0.1, 0.05]
  ai_conv_dis_thresh = [0.05, 0.05]
  ai_colors = ['red', 'orange', 'yellow']
  bact_speed = [0.2, 0.15]
  bact_ai_map = {BactType.BACT_1: AIType.BACT_1}#, 
              #BactType.BACT_2: AIType.BACT_2}


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
      

  #return
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
  bact_all_coords_map = dict()
  bact_all_colors_map = dict()
  ai_all_coords_map = dict()


  for i in xrange(BACT_TYPE_COUNT):
    color_map = [bact_colors[i]] * BACT_COUNT
    bact_all_colors_map[i] = color_map

  step_count = 0
  while(1):
    step_count += 1
    time.sleep(0.5)
    
    plt.cla()

    region = 110  # for pylab 2x2 subplot layout
    plt.subplots_adjust(left=0, right=1, bottom=0, top=0.95, wspace=0.01, hspace=0.01)
        
    region += 1
    plt.subplot(region)
    plt.title("p = %6.3f" % (23))

    for i in xrange(BACT_TYPE_COUNT):
      posDicts = bactPosInfo[i]
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])
      num = bactGraphs[i].number_of_nodes()
      num2 = len(all_coords)
      print ("NODE COUNT %d coordCOUNT %d", num, num2)
      nx.draw(bactGraphs[i], all_coords, 
        node_color=bact_all_colors_map[i], with_labels=False, node_size=40)

    for i in xrange(AI_TYPE_COUNT):
      color_map = [ai_colors[i]] * AI_COUNT
      posDicts = aiPosInfo[i]
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])
      ai_all_coords_map[i] = all_coords
      
      nx.draw(aiGraphs[i], all_coords, 
        node_color=color_map, with_labels=False, node_size=10)

    # START OF STEP COUNT % 2 == 0
    if (step_count % 5 != 0):
    #if (True):
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
        maxAIQuadrants += [(maxR,maxC)]

      for b in xrange(BACT_TYPE_COUNT):
        ai = bact_ai_map[b]
        (maxr, maxc) = maxAIQuadrants[ai]
        # find the mid point of this quadrant
        (x, y) = ((maxc + 0.5) * AIC_CFRAC,
                  (maxr + 0.5) * AIC_RFRAC)

        posDicts = bactPosInfo[b]
        all_coords = dict();
        for r in xrange(AIC_HCOUNT):
          for c in xrange(AIC_WCOUNT):
            all_coords.update(posDicts[r][c])
        
        bact_all_coords_map[b] = all_coords

        bactPosInfo = bactPosInfoOrig
        for node in all_coords.keys():
          coords = all_coords[node]
          (dx, dy) = (x - coords[0], y - coords[1])
          # only change if x & y are not already in the same
          # otherwise move them in a random direction
          orig_r = int(coords[1] / AIC_RFRAC)
          orig_c = int(coords[0] / AIC_CFRAC)
          if (orig_r == maxr and orig_c == maxc):
            xrandstart = int(maxc * AIC_CFRAC * 10000)
            xrandend = int((maxc + 1) * AIC_CFRAC * 10000)
            yrandstart = int((maxr * AIC_RFRAC * 10000))
            yrandend = int((maxr + 1) * AIC_RFRAC * 10000)
            newx = random.randint(xrandstart, xrandend) / 10000.0
            newy = random.randint(yrandstart, yrandend) / 10000.0
            #(newx, newy) = (coords[0] + dx*bact_speed[b], 
            #                coords[1] + dy*bact_speed[b])

          else:
            #print orig_r, orig_c , maxr, maxc
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

      for b in xrange(BACT_TYPE_COUNT):
        # we need to connect the bacteria which are kind of close enough
        all_coords = bact_all_coords_map[b]
        for node1 in all_coords.keys():
          posArray1 = all_coords[node1]
          x1 = posArray1[0]
          y1 = posArray1[1]
          for node2 in all_coords.keys():
            if (node1 != node2):
              posArray2 = all_coords[node2]
              x2 = posArray2[0]
              y2 = posArray2[1]
              dis = (abs(x1 - x2)**2 + abs(y1 - y2)**2)**(0.5)
              
              if (dis < bact_dis_thresh[b]):
                if (not bactGraphs[b].has_edge(node1, node2)):
                  bactGraphs[b].add_edge(node1, node2)
        ai = bact_ai_map[b]
        aicoords = ai_all_coords_map[ai]
        removeList = []

        for nodeai in aicoords.keys():
          ai_pos = aicoords[nodeai]
          aix = ai_pos[0]
          aiy = ai_pos[1]
          
          for nodeb in all_coords.keys():
            bpos = all_coords[nodeb]
            bx = bpos[0]
            by = bpos[1]
            dis = (abs(bx - aix)**2 + abs(by - aiy)**2)**(0.5)
            if (dis < ai_conv_dis_thresh[ai]):
              if (nodeai in removeList):
                continue
              else:
                removeList += [nodeai]
              aiGraphs[ai].remove_node(nodeai)
              ai_all_coords_map[ai].pop(nodeai)
              for r in xrange(AIC_HCOUNT):
                for c in xrange(AIC_WCOUNT):
                  aiPosInfo[ai][r][c].pop(nodeai, None)
              bact_all_colors_map[b][nodeb] = 'red'
    #END OF STEP COUNT % 2 == 0
    
    else:
      for b in xrange(BACT_TYPE_COUNT):
        bg = bactGraphs[b]
        bgn = bg.number_of_nodes();
        # we want to double this number
        for i in xrange(bgn):
          bactGraphs[b].add_node(i + bgn)
        print "orignumber %d, added number ", bgn, bactGraphs[b].number_of_nodes()
        # this is a dictionary
        add_all_coords = bact_all_coords_map[b]
        this_bact_all_coords = dict()
        for node in add_all_coords.keys():
          value1 = add_all_coords[node]
          print (node, value1)
          value2 = numpy.empty_like(value1)
          value2[:] = value1
          value2[0] = random.uniform(value2[0]-0.2, value2[0]+0.2)
          value2[1] = random.uniform(value2[1]-0.2, value2[1]+0.2)
          value2[0] = min(value2[0], 0.9999999)
          value2[0] = max(0, value2[0])
          value2[1] = min(value2[1], 0.9999999)
          value2[1] = max(0, value2[0])
          new_node = node + bgn
          this_bact_all_coords[new_node] = value2
          this_bact_all_coords[node] = value1
          print "node = %d new node = %d", node, new_node

        print bact_all_coords_map[b]
        print this_bact_all_coords
        bact_all_coords_map[b] = this_bact_all_coords
        # they start off at original color
        
        bact_all_colors_map[b] = bact_all_colors_map[b] + [bact_colors[b]] * bgn
        print bact_all_colors_map[b]

        bactPosInfo = bactPosInfoOrig
        for node in bact_all_coords_map[b].keys():
          posArray = bact_all_coords_map[b][node]
          x = posArray[0]
          y = posArray[1]
          c = int(x / AIC_CFRAC)
          r = int(y / AIC_RFRAC)
          print b, r, c, node
          bactPosInfo[b][r][c][node] = posArray

    plt.draw()
    plt.pause(1e-17)
    
  #global plot
  #plt.plot.show()

if __name__ == "__main__":
    main()
