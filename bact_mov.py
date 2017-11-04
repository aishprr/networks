#!/usr/bin/env python
"""""""""""""""""""""
|===================|
| Bacteria movement |
|===================|

"""""""""""""""""""""

import math
import time
from enum import Enum
import random
import copy

import numpy as np

import matplotlib.pyplot as plt

import networkx as nx

from functools import reduce
from operator import or_

def merge(*dicts):
    return 
    { k: reduce(lambda d, x: x.get(k, d), dicts, None) for k in reduce(or_, map(lambda x: x.keys(), dicts), set()) }

BACT_INIT_COUNT = [10, 10]
BACT_COUNT_LIMIT = [500, 500]
AI_INIT_COUNT = [5, 5]
AIC_HCOUNT = 4
AIC_WCOUNT = 4
MACRO_INIT_COUNT = 1

STEP_MULTIPLE = 10

MACRO_SPEED = 0.05
BACT_DEG_THRESH = 3
BACT_CHILD_DIS = 0.05
MACRO_EAT_DIS = 0.1
MACRO_MAX_BACT_EAT = 5

BACT_DRAW_SIZE = 40
AI_DRAW_SIZE = 20

# SIZE and BACT_WITHIN_RAD are related
MACRO_DRAW_SIZE = 600
MACRO_BACT_WITHIN_RAD = 0.02

# must add to 1
BACT_STRENGTH = [0.1, 0.9]

AIC_RFRAC = 1.0 / AIC_HCOUNT
AIC_CFRAC = 1.0 / AIC_WCOUNT

MACRO_MOVE_IN_GRID = AIC_CFRAC / 10.0

TOT_BACT_TYPE_COUNT = 2
TOT_AI_TYPE_COUNT = 2

AI_TYPE_COUNT = 2
BACT_TYPE_COUNT = 2
class AIType:
  BACT_1, BACT_2 = range(TOT_AI_TYPE_COUNT)

class BactType:
  BACT_1, BACT_2 = range(TOT_BACT_TYPE_COUNT)

AI_PER_BAC = 1

MACRO_COLOR = 'pink'
GRAPH = 'graph'
GRAPHPOS = 'graphpos'
# this is a list for each type of bacteria inside
# that macro
BACTTYPELIST = 'bacttypelist'
GRAPHCOLORMAP = 'graphcolormap'

def main():
  bact_count = BACT_INIT_COUNT
  bact_colors = ['blue', 'green']
  bact_dis_thresh = [0.1, 0.05]
  ai_conv_dis_thresh = [0.05, 0.05]
  ai_colors = ['yellow', 'purple']
  bact_speed = [0.2, 0.15]
  bact_ai_map = {BactType.BACT_1: AIType.BACT_1, 
              BactType.BACT_2: AIType.BACT_2}
  inverseStrength = []            


  macroGraph = nx.empty_graph(MACRO_INIT_COUNT)
  macroPos = nx.random_layout(macroGraph)
  macroCount = MACRO_INIT_COUNT

  macroInfo = dict()
  for node in macroGraph.nodes():
    # for each node, initialize the data structures inside the 
    # dictionary for it
    macroInfo[node] = dict()
    macroInfo[node][GRAPH] = nx.empty_graph(0)
    macroInfo[node][GRAPHPOS] = dict()
    macroInfo[node][BACTTYPELIST] = dict()
    macroInfo[node][GRAPHCOLORMAP] = []
    # has no nodes inside that macro for now
    # but when they are there, then I will initialize the thing
  
  totalInverseStrength = 0
  for elem in BACT_STRENGTH:
    totalInverseStrength += 1.0 / elem
  for elem in BACT_STRENGTH:
    inverseStrength += [(1.0 / elem) / totalInverseStrength]

  # initial data structures
  bactGraphs = []
  aiGraphs = []
  bactPosInfoOrig = []
  bactPosInfo = []
  for i in xrange(BACT_TYPE_COUNT):
    b = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    bactPosInfoOrig += [b]
  bactPosInfo = copy.deepcopy(bactPosInfoOrig)
  first = True
  bact_all_coords_map = dict()
  bact_all_colors_map = dict()
  ai_all_coords_map = dict()
  
  aiPosInfo = []
  for i in xrange(AI_TYPE_COUNT):
    a = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    aiPosInfo += [a]

  for b in xrange(BACT_TYPE_COUNT):
    bg = nx.empty_graph(BACT_INIT_COUNT[b])
    bactGraphs += [bg]

    initialBgPos = nx.random_layout(bg)
    
    bact_all_coords_map[b] = initialBgPos
    print initialBgPos

  for a in xrange(AI_TYPE_COUNT):
    ai = nx.empty_graph(AI_INIT_COUNT[a])
    aiGraphs += [ai]
    initialAiPos = nx.random_layout(ai)
    for node in initialAiPos.keys():
      posArray = initialAiPos[node]
      x = posArray[0]
      y = posArray[1]
      c = int(x / AIC_CFRAC)
      r = int(y / AIC_RFRAC)

      aiPosInfo[a][r][c][node] = posArray

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

  

  for i in xrange(BACT_TYPE_COUNT):
    color_map = [bact_colors[i]] * BACT_INIT_COUNT[i]
    bact_all_colors_map[i] = color_map

  step_count = 0
  while(1):
    step_count += 1
    time.sleep(1)
    
    plt.cla()

    region = 110  # for pylab 2x2 subplot layout
    plt.subplots_adjust(left=0, right=1, bottom=0, top=0.95, wspace=0.01, hspace=0.01)
        
    region += 1
    plt.subplot(region)
    plt.title("Bacteria Network Simulation")

    # update the Bact positions in the grid map
    for b in xrange(BACT_TYPE_COUNT):
      all_coords = bact_all_coords_map[b]
      for node in all_coords.keys():
        posArray = all_coords[node]
        x = posArray[0]
        y = posArray[1]
        c = int(x / AIC_CFRAC)
        r = int(y / AIC_RFRAC)
        bactPosInfo[b][r][c][node] = posArray

    for i in xrange(BACT_TYPE_COUNT):
      all_coords = bact_all_coords_map[i]
      num = bactGraphs[i].number_of_nodes()
      bact_count[i] = num
      num2 = len(all_coords)
      print ("bact %d NODE COUNT %d coordCOUNT %d", i, num, num2)
      if (num != 0):
        nx.draw(bactGraphs[i], all_coords, edge_color=bact_colors[i], alpha=1,
          node_color=bact_all_colors_map[i], with_labels=False, 
          node_size=BACT_DRAW_SIZE)

    for i in xrange(AI_TYPE_COUNT):
      color_map = [ai_colors[i]] * aiGraphs[i].number_of_nodes()
      posDicts = aiPosInfo[i]
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])
      ai_all_coords_map[i] = all_coords
      
      nx.draw(aiGraphs[i], all_coords, 
        node_color=color_map, with_labels=False, node_size=AI_DRAW_SIZE)

    nx.draw(macroGraph, macroPos, alpha=0.5,
        node_color=MACRO_COLOR, with_labels=False, node_size=MACRO_DRAW_SIZE)

    for m in macroInfo:
      print "macro " + str(m) + " color map " + str(macroInfo[m][GRAPHCOLORMAP])
      nx.draw(macroInfo[m][GRAPH], macroInfo[m][GRAPHPOS],
        node_color=macroInfo[m][GRAPHCOLORMAP], 
        with_labels=False, node_size=BACT_DRAW_SIZE)

    # START OF STEP COUNT % 2 == 0
    if (step_count % STEP_MULTIPLE in [1,2,3,4, 5, 6, 7, 8]):
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

      #bactPosInfo = bactPosInfoOrig
      for b in xrange(BACT_TYPE_COUNT):
        ai = bact_ai_map[b]
        (maxr, maxc) = maxAIQuadrants[ai]
        # find the mid point of this quadrant
        (x, y) = ((maxc + 0.5) * AIC_CFRAC,
                  (maxr + 0.5) * AIC_RFRAC)
        for (m, mCoord) in macroPos.iteritems(): 
          macroR = int(mCoord[0] / AIC_RFRAC)
          macroC = int(mCoord[1] / AIC_CFRAC)
          if ((macroR == maxr) and (macroC == maxc)):
            # move by a small random amount only, once you get to the place
            deltaX = MACRO_MOVE_IN_GRID * np.random.random_sample()
            deltaY = MACRO_MOVE_IN_GRID * np.random.random_sample()
          else:
            (dx, dy) = (x - mCoord[0], y - mCoord[1])
            # we get the direction as above
            deltaX = dx * inverseStrength[b] * MACRO_SPEED
            deltaY = dy * inverseStrength[b] * MACRO_SPEED
          newx = mCoord[0] + deltaX
          newy = mCoord[1] + deltaY
          macroPos[m][0] = newx
          macroPos[m][1] = newy
          # fix for all the bacteria inside this macro
          for mBactNode in macroInfo[m][GRAPHPOS].keys():
            posArray = macroInfo[m][GRAPHPOS][mBactNode]
            posArray[0] = posArray[0] + deltaX
            posArray[1] = posArray[1] + deltaY
            macroInfo[m][GRAPHPOS][mBactNode] = posArray
            

        posDicts = bactPosInfo[b]
        all_coords = dict();
        for r in xrange(AIC_HCOUNT):
          for c in xrange(AIC_WCOUNT):
            all_coords.update(posDicts[r][c])
        
        bact_all_coords_map[b] = all_coords

        
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
            
          else:
            #print orig_r, orig_c , maxr, maxc
            (newx, newy) = (coords[0] + dx*bact_speed[b], 
                            coords[1] + dy*bact_speed[b])
          all_coords[node][0] = newx;
          all_coords[node][1] = newy;

      for b in xrange(BACT_TYPE_COUNT):
        # we need to connect the bacteria which are kind of close enough
        all_coords = bact_all_coords_map[b]
        #print "all_coords when connecting bacterial network"
        #print all_coords
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
              #print "nodes in this bact graph are: " + str(bactGraphs[b].nodes())
              #print "node1 = " + str(node1)
              #print "node2 = " + str(node2)
              
              if (dis < bact_dis_thresh[b] 
                and bactGraphs[b].degree(node1) < BACT_DEG_THRESH
                and bactGraphs[b].degree(node2) < BACT_DEG_THRESH):
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
              bact_all_colors_map[b][nodeb] = 'orange'
    #END OF STEP COUNT % 2 == 0
    elif (step_count % STEP_MULTIPLE in [9]):
      bactPosInfo = copy.deepcopy(bactPosInfoOrig)
      #print "orig " + str(bactPosInfoOrig)
      #print "bactPosInfo copy orig " + str(bactPosInfo)
      for b in xrange(BACT_TYPE_COUNT):
        # remove nodes from the bacteria if it is very close to something
        bg = bactGraphs[b]
        bgList = bg.nodes()
        removeBgList = []
        for bnode in bgList:
          posArray = bact_all_coords_map[b][bnode]
          bx = posArray[0]
          by = posArray[1]
          for (m, mpos) in macroPos.iteritems():
            if bnode in removeBgList:
              continue
            macroNodeCount = macroInfo[m][GRAPH].number_of_nodes()
            if (macroNodeCount >= MACRO_MAX_BACT_EAT):
              continue
            mx = mpos[0]
            my = mpos[1]
            dis = (abs(mx - bx)**2 + abs(my - by)**2)**(0.5)
            if (dis < MACRO_EAT_DIS):
              newMacBactNode = macroNodeCount + 1
              # I want to remove this particular graph node
              removeBgList += [bnode]
              # add this node to this macro's thing
              macroInfo[m][GRAPH].add_node(newMacBactNode)
              randTheta = 2 * np.pi * np.random.random_sample()
              randRadius = MACRO_BACT_WITHIN_RAD * np.random.random_sample()
              randX = mx + randRadius * np.cos(randTheta)
              randY = my + randRadius * np.sin(randTheta)
              # now clip this between 0 and 0.9999
              newNodeCoords = np.array([np.clip(randX, 0, 0.999999), 
                                  np.clip(randY, 0, 0.999999)])
              macroInfo[m][GRAPHPOS][newMacBactNode] = newNodeCoords
              try:
                macroInfo[m][BACTTYPELIST][b] += [newMacBactNode]
              except KeyError, e:
                macroInfo[m][BACTTYPELIST][b] = [newMacBactNode]
              print "adding color " + str(bact_colors[b])
              macroInfo[m][GRAPHCOLORMAP] += [bact_colors[b]]
              # add a position for this
        #print "removing " + str(removeBgList)
        for remove in removeBgList:
          try:
            bg.remove_node(remove)
            bact_all_coords_map[b].pop(remove, None)
          except nx.exception.NetworkXError:
            continue

    elif (step_count % STEP_MULTIPLE in [0]):
      for b in xrange(BACT_TYPE_COUNT):
        if (bact_count[b] >= BACT_COUNT_LIMIT[b]):
          continue
        bg = bactGraphs[b]
        bgn = bg.number_of_nodes();

        max_node_number = -1
        nodeList = bg.nodes()
        for i in nodeList:
          max_node_number = max(max_node_number, i)

        #print "max node I found uptil now is: " + str(max_node_number)

        # we want to double this number
        for i in xrange(bgn):
          #print "adding new node number " + str(i + max_node_number + 1)
          bactGraphs[b].add_node(i + max_node_number + 1)
        #print "orignumber %d, added number ", bgn, bactGraphs[b].number_of_nodes()
        # this is a dictionary
        add_all_coords = bact_all_coords_map[b]
        this_bact_all_coords = dict()
        adding_node_number = 1
        for node in add_all_coords.keys():
          value1 = add_all_coords[node]
          value2 = np.empty_like(value1)
          #print (node, value1)
          #print (node, np.add(value1[0],-BACT_CHILD_DIS), 
          #  np.add(value1[0], BACT_CHILD_DIS))
          #print (node, np.add(value1[1],-BACT_CHILD_DIS), 
          #  np.add(value1[1], BACT_CHILD_DIS))
          value2[:] = value1
          value2[0] = np.clip(np.random.uniform(np.add(value1[0],-BACT_CHILD_DIS), 
            np.add(value1[0], 0.05)), 0, 0.999999)
          value2[1] = np.clip(np.random.uniform(np.add(value1[1],-BACT_CHILD_DIS), 
            np.add(value1[1],0.05)), 0, 0.999999)
          new_node = adding_node_number + max_node_number
          adding_node_number = adding_node_number + 1
          this_bact_all_coords[new_node] = value2
          this_bact_all_coords[node] = value1
          #print "node = %d new node = %d", node, new_node
          #print (node, value1)
          #print (new_node, value2)
          
        #print bact_all_coords_map[b]
        #print this_bact_all_coords
        bact_all_coords_map[b] = this_bact_all_coords
        # they start off at original color
        
        bact_all_colors_map[b] = bact_all_colors_map[b] + [bact_colors[b]] * bgn
        #print bact_all_colors_map[b]

        bactPosInfo = copy.deepcopy(bactPosInfoOrig)
        # need to update this here, because we need this number for the 
        # AI movement immediately in this step
        for node in bact_all_coords_map[b].keys():
          posArray = bact_all_coords_map[b][node]
          x = posArray[0]
          y = posArray[1]
          c = int(x / AIC_CFRAC)
          r = int(y / AIC_RFRAC)
          bactPosInfo[b][r][c][node] = posArray

      # change the AI number, add more of them
      for b in xrange(BACT_TYPE_COUNT):
        if (bact_count[b] >= BACT_COUNT_LIMIT[b]):
          continue
        total_b = bactGraphs[b].number_of_nodes()
        ai = bact_ai_map[b]
        for r in xrange(AIC_HCOUNT):
          for c in xrange(AIC_WCOUNT):
            bCount = len(bactPosInfo[b][r][c].keys())
            # we want to add so many in this quadrant
            newAiCount = int(bCount * AI_PER_BAC)
            xrandstart = int(c * AIC_CFRAC * 10000)
            xrandend = int((c + 1) * AIC_CFRAC * 10000)
            yrandstart = int((r * AIC_RFRAC * 10000))
            yrandend = int((r + 1) * AIC_RFRAC * 10000)
            for i in xrange(newAiCount):
              newNode = i + newAiCount
              aiGraphs[ai].add_node(newNode)
              newx = random.randint(xrandstart, xrandend) / 10000.0
              newy = random.randint(yrandstart, yrandend) / 10000.0
              newCoord = np.array([newx, newy])
              aiPosInfo[ai][r][c][newNode] = newCoord
        posDicts = aiPosInfo[ai]
        all_coords = dict();
        for r in xrange(AIC_HCOUNT):
          for c in xrange(AIC_WCOUNT):
            all_coords.update(posDicts[r][c])
        ai_all_coords_map[ai] = posDicts
    
    plt.draw()
    plt.pause(1e-17)
    
  #global plot
  #plt.plot.show()

if __name__ == "__main__":
    main()
