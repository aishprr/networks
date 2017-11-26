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
import matplotlib.image as mpimg

import networkx as nx

from functools import reduce
from operator import or_

def merge(*dicts):
    return 
    { k: reduce(lambda d, x: x.get(k, d), dicts, None) for k in reduce(or_, map(lambda x: x.keys(), dicts), set()) }

BACT_INIT_COUNT = [25, 10, 5]
BACT_COUNT_LIMIT = [20, 20, 20]
AI_INIT_COUNT = [5, 5, 5]
AIC_HCOUNT = 4
AIC_WCOUNT = 4
MACRO_INIT_COUNT = 10
HELPTCELL_INIT_COUNT = 20
KILLTCELL_INIT_COUNT = 0

STEP_MULTIPLE = 10

MACRO_SPEED = 0.1
HELP_SPEED = 0.2
KILL_SPEED = 0.3

BACT_DEG_THRESH = 3
BACT_CHILD_DIS = 0.05
MACRO_EAT_DIS = 0.1
KILLTCELL_EAT_DIS = 0.05
HELPER_MACRO_DIS = 0.15
MACRO_MAX_BACT_EAT = 2

MACRO_BACT_TO_DIE = 10

HELPTCELL_DEG_THRESH = 2
HELPTCELL_EDGE_DIS = 0.08

BACT_IN_MACRO_REPR = [1, 1, 3]
BACT_IN_MACRO_REPRO_AGE = 5

BACT_DRAW_SIZE = 50
AI_DRAW_SIZE = 10
HELPTCELL_DRAW_SIZE = 50
HELPTCELL_DRAW_SHAPE = "8"
KILLTCELL_DRAW_SIZE = 50
KILLTCELL_DRAW_SHAPE = "s"

# SIZE and BACT_WITHIN_RAD are related
MACRO_DRAW_SIZE = 600
MACRO_BACT_WITHIN_RAD = 0.02

KILLTCELL_NEW_HELP = 0.02
KILL_PER_MAC = 4
KILL_DISPR_BAC_THRESH = sum(BACT_INIT_COUNT)

# must add to 1
BACT_STRENGTH =  [0.1, 0.2, 0.7]
BACT_IN_MACRO_KILL_THRESH = 0.4
BACT_IN_MACRO_REPRO_THRESH = 0.4

AIC_RFRAC = 1.0 / AIC_HCOUNT
AIC_CFRAC = 1.0 / AIC_WCOUNT

MACRO_MOVE_IN_GRID = AIC_CFRAC * 1.5
HELP_MOVE_IN_GRID = AIC_CFRAC * 1.5
KILL_MOVE_IN_GRID = AIC_CFRAC * 0.6

BACT_RAND_LIMITS = [[[0, 0.5], [0, 0.5]], [[0, 0.5], [0.5, 1]], [[0.5, 1], [0, 0.5]]] 

TOT_BACT_TYPE_COUNT = 3
TOT_AI_TYPE_COUNT = 3

AI_TYPE_COUNT = 3
BACT_TYPE_COUNT = 3
class AIType:
  BACT_1, BACT_2, BACT_3 = range(TOT_AI_TYPE_COUNT)

class BactType:
  BACT_1, BACT_2, BACT_3 = range(TOT_BACT_TYPE_COUNT)

AI_PER_BAC = 3

MACRO_COLOR = 'salmon'
HELPTCELL_COLOR = 'blue'
KILLTCELL_COLOR = 'black'
GRAPH = 'graph'
GRAPHPOS = 'graphpos'
# this is a list for each type of bacteria inside
# that macro
BACTTYPELIST = 'bacttypelist'
GRAPHCOLORMAP = 'graphcolormap'
IMAGE = 'image'
BACTEATAGE = 'bacteatage'

def main():
  bact_count = BACT_INIT_COUNT
  bact_colors = ['gold', 'orange', 'red']
  bact_dis_thresh = [0.1, 0.05, 0.05]
  ai_conv_dis_thresh = [0.05, 0.05, 0.05]
  ai_colors = ['limegreen', 'yellowgreen', 'green']
  bact_speed = [0.2, 0.15, 0.2]
  bact_ai_map = {BactType.BACT_1: AIType.BACT_1, 
              BactType.BACT_2: AIType.BACT_2, 
              BactType.BACT_3: AIType.BACT_3}
  inverseStrength = []            


  macroGraph = nx.empty_graph(MACRO_INIT_COUNT)
  macroPos = nx.random_layout(macroGraph)
  macroCount = MACRO_INIT_COUNT

  helptcellImg = mpimg.imread('tcell.png')
  helptcellGraph = nx.empty_graph(HELPTCELL_INIT_COUNT)
  helptcellPos = nx.random_layout(helptcellGraph)
  helptcellCount = HELPTCELL_INIT_COUNT

  killtcellImg = mpimg.imread('tcell.png')
  killtcellGraph = nx.empty_graph(KILLTCELL_INIT_COUNT)
  killtcellPos = nx.random_layout(killtcellGraph)
  killtcellCount = KILLTCELL_INIT_COUNT

  for t in xrange(helptcellCount): 
    helptcellGraph.node[t][IMAGE] = helptcellImg


  macroInfo = dict()
  for node in macroGraph.nodes():
    # for each node, initialize the data structures inside the 
    # dictionary for it
    macroInfo[node] = dict()
    macroInfo[node][GRAPH] = nx.empty_graph(0)
    macroInfo[node][GRAPHPOS] = dict()
    macroInfo[node][BACTTYPELIST] = dict()
    macroInfo[node][GRAPHCOLORMAP] = []
    macroInfo[node][BACTEATAGE] = dict()
    # has no nodes inside that macro for now
    # but when they are there, then I will initialize the thing
  
  totalInverseStrength = 0
  for elem in BACT_STRENGTH:
    totalInverseStrength += 1.0 / elem
  for elem in BACT_STRENGTH:
    inverseStrength += [(1.0 / elem) / totalInverseStrength]
  print inverseStrength


  # initial data structures
  bactGraphs = []
  deadMacros = []
  aiGraphs = []
  bactPosInfoOrig = []
  bactPosInfo = []
  for i in xrange(BACT_TYPE_COUNT):
    b = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    bactPosInfoOrig += [b]
  bactPosInfo = copy.deepcopy(bactPosInfoOrig)
  first = True
  bact_all_coords_map = dict()
  #bact_all_colors_map = dict()
  ai_all_coords_map = dict()
  
  aiPosInfo = []
  for i in xrange(AI_TYPE_COUNT):
    a = [[dict() for i in range(AIC_WCOUNT)] for j in range(AIC_HCOUNT)]
    aiPosInfo += [a]

  for b in xrange(BACT_TYPE_COUNT):
    bg = nx.empty_graph(BACT_INIT_COUNT[b])
    bactGraphs += [bg]
    initialBgPos = nx.random_layout(bg)
    #print type(initialBgPos)
    newBgPos = dict()
    for (node, nodeCoords) in initialBgPos.iteritems():
      nxx = nodeCoords[0]
      nyy = nodeCoords[1]
      limits = BACT_RAND_LIMITS[b]
      nxx = limits[0][0] + (limits[0][1] - limits[0][0])* nxx
      nyy = limits[1][0] + (limits[1][1] - limits[1][0])* nyy
      newBgPos[node] = np.array([nxx, nyy])

    bact_all_coords_map[b] = newBgPos
    #print newBgPos
    #print initialBgPos


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
      if (r >= 2):
        r = 1
      if (c >= 2):
        c = 1
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

  

  # for i in xrange(BACT_TYPE_COUNT):
  #   color_map = [bact_colors[i]] * BACT_INIT_COUNT[i]
  #   bact_all_colors_map[i] = color_map

  step_count = 0
  while(1):

    #print bact_all_coords_map[BactType.BACT_1]

    step_count += 1
    time.sleep(0.5)
    
    plt.cla()

    region = 110  # for pylab 2x2 subplot layout
    plt.subplots_adjust(left=0, right=1, bottom=0, top=0.95, wspace=0.01, hspace=0.01)
        
    region += 1
    plt.subplot(region)
    plt.title("Bacteria Network Simulation")

    bactPosInfo = copy.deepcopy(bactPosInfoOrig)
    # update the Bact positions in the grid map
    for b in xrange(BACT_TYPE_COUNT):
      all_coords = bact_all_coords_map[b]
      for node in all_coords.keys():
        posArray = all_coords[node]
        x = posArray[0]
        y = posArray[1]
        c = int(x / AIC_CFRAC)
        r = int(y / AIC_RFRAC)
        if (r >= 2):
          r = 1
        if (c >= 2):
          c = 1
        try:
          bactPosInfo[b][r][c][node] = posArray
        except IndexError:
          print "index error somewhere!!!"
          # print "b = " + str(b)
          # print "r = " + str(r)
          # print "c = " + str(c)
          # print "node = " + str(node)

    bactCountList = []
    for i in xrange(BACT_TYPE_COUNT):
      all_coords = bact_all_coords_map[i]
      num = bactGraphs[i].number_of_nodes()
      bact_count[i] = num
      bactCountList += [num]
      num2 = len(all_coords)
      #print ("bact %d NODE COUNT %d coordCOUNT %d", i, num, num2)
      if (num != 0):
        nx.draw(bactGraphs[i], all_coords, edge_color=bact_colors[i], alpha=1,
          node_color=bact_colors[i], with_labels=False, 
          node_size=BACT_DRAW_SIZE)
    print step_count, bactCountList

    for i in xrange(AI_TYPE_COUNT):
      # print ai_colors[i], aiGraphs[i].number_of_nodes()
      color_map = [ai_colors[i]] * aiGraphs[i].number_of_nodes()
      posDicts = aiPosInfo[i]
      all_coords = dict();
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          all_coords.update(posDicts[r][c])
      ai_all_coords_map[i] = all_coords
      
      #nx.draw(aiGraphs[i], all_coords, 
      #  node_color=color_map, with_labels=False, node_size=AI_DRAW_SIZE)

    macroCount = macroGraph.number_of_nodes()
    nx.draw(macroGraph, macroPos, alpha=0.5,
        node_color=MACRO_COLOR, with_labels=False, node_size=MACRO_DRAW_SIZE)

    helptcellCount = helptcellGraph.number_of_nodes()
    nx.draw(helptcellGraph, helptcellPos, with_labels=False, 
      node_color=HELPTCELL_COLOR, node_shape=HELPTCELL_DRAW_SHAPE, 
      node_size=HELPTCELL_DRAW_SIZE, edge_color=HELPTCELL_COLOR)

    killtcellCount = killtcellGraph.number_of_nodes()
    nx.draw(killtcellGraph, killtcellPos, with_labels=False, 
      node_color=KILLTCELL_COLOR, node_shape=KILLTCELL_DRAW_SHAPE, 
      node_size=KILLTCELL_DRAW_SIZE)

    print deadMacros
    for m in macroInfo:
      if m in deadMacros:
        continue
      #print "macro " + str(m) + " color map " + str(macroInfo[m][GRAPHCOLORMAP])
      nx.draw(macroInfo[m][GRAPH], macroInfo[m][GRAPHPOS],
        node_color=macroInfo[m][GRAPHCOLORMAP], 
        with_labels=False, node_size=BACT_DRAW_SIZE)

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

    maxBactQuadrants = []
    for b in xrange(BACT_TYPE_COUNT):
      posDicts = bactPosInfo[b]
      maxLen = 0
      (maxR, maxC) = (-1, -1)
      for r in xrange(AIC_HCOUNT):
        for c in xrange(AIC_WCOUNT):
          if (len(posDicts[r][c]) > maxLen):
            maxLen = len(posDicts[r][c])
            maxR = r
            maxC = c
      maxBactQuadrants += [(maxR,maxC)]
    #print "MAX BAXT QUADRANTS: " + str(maxBactQuadrants)

    # START OF STEP COUNT % 2 == 0
    if (step_count % STEP_MULTIPLE in [1,2,3,6,7,8]):
      # we have information about the previous positions.
      # so, now go through each and find the quadrant with max number
      # of ais of each type
      for b in xrange(BACT_TYPE_COUNT):
        ai = bact_ai_map[b]
        (maxAr, maxAc) = maxAIQuadrants[ai]
        (maxBr, maxBc) = maxBactQuadrants[b]
        # find the mid point of this quadrant
        (ax, ay) = ((maxAc + 0.5) * AIC_CFRAC,
                  (maxAr + 0.5) * AIC_RFRAC)
        (bx, by) = ((maxBc + 0.5) * AIC_CFRAC,
                  (maxBr + 0.5) * AIC_RFRAC)
        # there are 0 of this bacteria, so it shouldn't do anything for us
        if (maxBc == -1 or maxBr == -1 or bact_count[b] == 0):
          continue
        for (m, mCoord) in macroPos.iteritems(): 
          if m in deadMacros:
            continue
          macroR = int(mCoord[0] / AIC_RFRAC)
          macroC = int(mCoord[1] / AIC_CFRAC)
          # macro move towards bacteria
          if ((macroR == maxBr) and (macroC == maxBc)):
            # move by a small random amount only, once you get to the place
            deltaX = MACRO_MOVE_IN_GRID * np.random.random_sample()
            deltaY = MACRO_MOVE_IN_GRID * np.random.random_sample()
          else:
            (dx, dy) = (bx - mCoord[0], by - mCoord[1])
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

        # move the helper tcells to the popular place
        for (helpt, helpCoord) in helptcellPos.iteritems():  
          helpR = int(helpCoord[0] / AIC_RFRAC)
          helpC = int(helpCoord[1] / AIC_CFRAC)
          if ((helpR == maxBr) and (helpC == maxBc)):
            # move by a small random amount only, once you get to the place
            deltaX = HELP_MOVE_IN_GRID * np.random.random_sample()
            deltaY = HELP_MOVE_IN_GRID * np.random.random_sample()
          else:
            (dx, dy) = (bx - helpCoord[0], by - helpCoord[1])
            deltaX = dx * inverseStrength[b] * HELP_SPEED
            deltaY = dy * inverseStrength[b] * HELP_SPEED
          newx = helpCoord[0] + deltaX
          newy = helpCoord[1] + deltaY
          helptcellPos[helpt][0] = newx
          helptcellPos[helpt][1] = newy

        # move the killer tcells to the popular place
        for (killt, killCoord) in killtcellPos.iteritems():  
          killR = int(killCoord[0] / AIC_RFRAC)
          killC = int(killCoord[1] / AIC_CFRAC)
          #print str((maxBr, maxBc)) +  str((killR, killC))
          if ((killR == maxBr) and (killC == maxBc)):
            # move by a small random amount only, once you get to the place
            deltaX = KILL_MOVE_IN_GRID * np.random.random_sample()
            deltaY = KILL_MOVE_IN_GRID * np.random.random_sample()
            #print "INSIDE CORRECT ONE"
          else:
            (dx, dy) = (bx - killCoord[0], by - killCoord[1])
            deltaX = dx * (1.0/len(bact_speed)) * KILL_SPEED
            deltaY = dy * (1.0/len(bact_speed)) * KILL_SPEED
          #print deltaX, deltaY
          newx = killCoord[0] + deltaX
          newy = killCoord[1] + deltaY
          killtcellPos[killt][0] = newx
          killtcellPos[killt][1] = newy

        posDicts = bactPosInfo[b]
        all_coords = dict();
        for r in xrange(AIC_HCOUNT):
          for c in xrange(AIC_WCOUNT):
            all_coords.update(posDicts[r][c])
        
        bact_all_coords_map[b] = all_coords

        # bacteria move towards ai
        for node in all_coords.keys():
          coords = all_coords[node]
          (dx, dy) = (ax - coords[0], ay - coords[1])
          # only change if x & y are not already in the same
          # otherwise move them in a random direction
          orig_r = int(coords[1] / AIC_RFRAC)
          orig_c = int(coords[0] / AIC_CFRAC)
          if (orig_r == maxAr and orig_c == maxAc):
            xrandstart = int(maxAc * AIC_CFRAC * 10000)
            xrandend = int((maxAc + 1) * AIC_CFRAC * 10000)
            yrandstart = int((maxAr * AIC_RFRAC * 10000))
            yrandend = int((maxAr + 1) * AIC_RFRAC * 10000)
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
              #bact_all_colors_map[b][nodeb] = 'orange'
      helptcelllist = helptcellPos.items()
      for i in xrange(len(helptcelllist)):
        for j in xrange(i):
          (hc1, hCoords1) = helptcelllist[i]
          (hc2, hCoords2) = helptcelllist[j]
          if ((helptcellGraph.degree(hc1) >= HELPTCELL_DEG_THRESH) or
                (helptcellGraph.degree(hc2) >= HELPTCELL_DEG_THRESH)):
            continue
          (hc1x, hc1y) = (hCoords1[0], hCoords1[1])
          (hc2x, hc2y) = (hCoords2[0], hCoords2[1])
          dis = (abs(hc1x - hc2x)**2 + abs(hc1y - hc2y)**2)**(0.5)
          if (dis < HELPTCELL_EDGE_DIS):
            helptcellGraph.add_edge(hc1, hc2)

    #END OF STEP COUNT % 2 == 0
    elif (step_count % STEP_MULTIPLE in [4, 9]):
      for (m, mpos) in macroPos.iteritems():
        if m in deadMacros:
          continue
        addBgDict = dict()
        for b2 in xrange(BACT_TYPE_COUNT):
          addBgDict[b2] = []
        macroNodeCount = macroInfo[m][GRAPH].number_of_nodes()
        if (macroNodeCount >= MACRO_BACT_TO_DIE):
          print "KILLING MACRO!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          # a lot of people inside, so kill this macro
          deadMacros += [m]
          macroGraph.remove_node(m)
          for bactType in macroInfo[m][BACTTYPELIST].keys():
            thisBactTypeNodes = macroInfo[m][BACTTYPELIST][bactType]
            for thisBactTypeNode in thisBactTypeNodes:
              addBgDict[bactType] += [macroInfo[m][GRAPHPOS][thisBactTypeNode]]
            print "ADDING SO MANY NODES HERE of type " + str(bactType) + "!!!! " + str(len(thisBactTypeNodes))
        for bactType in addBgDict.keys():
          addBgList = addBgDict[bactType]
          for add in addBgList:
            print "adding something here"
            newbg = bactGraphs[bactType]
            if len(newbg.nodes()) == 0:
              newNodeID = 0
            else:
              newNodeID = max(newbg.nodes()) + 1
            newbg.add_node(newNodeID)
            bact_all_coords_map[bactType][newNodeID] = add
            print str(bactType) + " adding new node " + str(newNodeID)
            print str(bactType) + " new number of nodes = " + str(newbg.number_of_nodes())
            print str(bactType) + " should match with " + str(len(bact_all_coords_map[bactType].keys()))

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
          # check if there's a macrophage close by an kill the bacteria
          for (m, mpos) in macroPos.iteritems():
            if m in deadMacros:
              continue
            if bnode in removeBgList:
              continue
            macroNodeCount = macroInfo[m][GRAPH].number_of_nodes()
            if (macroNodeCount >= MACRO_MAX_BACT_EAT):
              continue
            if (macroNodeCount == 0):
              macroNewNode = 0
            else:
              macroNewNode = max(macroInfo[m][GRAPH].nodes()) + 1
            
            mx = mpos[0]
            my = mpos[1]
            dis = (abs(mx - bx)**2 + abs(my - by)**2)**(0.5)
            if (dis < MACRO_EAT_DIS):
              newMacBactNode = macroNewNode
              macroNewNode += 1
              # I want to remove this particular graph node
              removeBgList += [bnode]
              # add this node to this macro's thing
              macroInfo[m][GRAPH].add_node(newMacBactNode)
              macroInfo[m][BACTEATAGE][newMacBactNode] = step_count
              randTheta = 2 * np.pi * np.random.random_sample()
              randRadius = MACRO_BACT_WITHIN_RAD# * np.random.random_sample()
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
              #print "adding color " + str(bact_colors[b])
              macroInfo[m][GRAPHCOLORMAP] += [bact_colors[b]]
          
          # each killer only kills one in each step
          fullKillers = []
          for (killer, killpos) in killtcellPos.iteritems():
            if bnode in removeBgList:
              continue
            if killer in fullKillers:
              continue
            kx = killpos[0]
            ky = killpos[1]
            dis = (abs(kx - bx)**2 + abs(ky - by)**2)**(0.5)
            if (dis < KILLTCELL_EAT_DIS):
              # I want to remove this particular graph node
              removeBgList += [bnode]
              fullKillers += [killer]
              
        for remove in removeBgList:
          try:
            bg.remove_node(remove)
            bact_all_coords_map[b].pop(remove, None)
          except nx.exception.NetworkXError:
            continue
        

    elif (step_count % STEP_MULTIPLE in [0, 5]):
      
      usedMacros = []
      newKillTCellPos = dict()
      for (killt, killCoord) in killtcellPos.iteritems():  
        killR = int(killCoord[0] / AIC_RFRAC)
        killC = int(killCoord[1] / AIC_CFRAC)
        #print (maxBr, maxBc)
        totalBactCountInGridCell = 0
        for b in xrange(BACT_TYPE_COUNT):
          totalBactCountInGridCell += bact_count[b]
          #totalBactCountInGridCell += len(bactPosInfo[b][killR][killC].keys())
        if (totalBactCountInGridCell < KILL_DISPR_BAC_THRESH):
          killtcellGraph.remove_node(killt)
          #killtcellPos.pop(killt, None)
          print "killing the kill tcell = " + str(killt)
        else: 
          newKillTCellPos[killt] = killCoord
      killtcellPos = newKillTCellPos

      if (killtcellGraph.number_of_nodes() == 0):
        maxKillTcellNode = -1
      else:
        maxKillTcellNode = max(killtcellGraph.nodes())

      for (helptnode, helptnodeCoords) in helptcellPos.iteritems():
        # if it is close to a macro node, then I need to then
        # make KILL_PER_MAC number of killer t cells and remove
        # the weak ones from that macrophage
        (hx, hy) = (helptnodeCoords[0], helptnodeCoords[1])
        for (m, mCoord) in macroPos.iteritems():
          if m in deadMacros:
            continue
          if (m in usedMacros):
            # only use 1 macro once for 1 helper T cell
            continue
          (mx, my) = (mCoord[0], mCoord[1])
          dis = (abs(mx - hx)**2 + abs(my - hy)**2)**(0.5)
          if (dis < HELPER_MACRO_DIS):
            usedMacros += [m]
            for new in xrange(KILL_PER_MAC):
              newKillTCellNode = maxKillTcellNode + 1
              maxKillTcellNode += 1
              # means we want to make a few killer t cells close to here
              killtcellGraph.add_node(newKillTCellNode)
              killtcellx = hx + KILLTCELL_NEW_HELP * np.random.random_sample()
              killtcelly = hy + KILLTCELL_NEW_HELP * np.random.random_sample()
              killtcellPos[newKillTCellNode] = np.array([killtcellx, killtcelly])
            # you added new killer for this macro, so now remove the weak ones 
            # from this macrophage
            for macroBactType in macroInfo[m][BACTTYPELIST].keys():
              if (BACT_STRENGTH[macroBactType] < BACT_IN_MACRO_KILL_THRESH):
                bactNodeList = macroInfo[m][BACTTYPELIST][macroBactType]
                for bact in bactNodeList:
                  # remove each of these from the graph
                  macroInfo[m][GRAPH].remove_node(bact)
                  # remove the node position properly
                  macroInfo[m][GRAPHPOS].pop(bact)
                  macroInfo[m][BACTEATAGE].pop(bact)
                  # don't care about colors in the list for now
                macroInfo[m][BACTTYPELIST][macroBactType] = []

      for (m, mCoord) in macroPos.iteritems():
        if m in deadMacros:
          continue
        mx = mCoord[0]
        my = mCoord[1]
        for macroBactType in macroInfo[m][BACTTYPELIST].keys():
          # for each type of bacteria
          if (BACT_STRENGTH[macroBactType] < BACT_IN_MACRO_REPRO_THRESH):
            continue
          # else we want to actually reproduce for everything
          if (macroInfo[m][GRAPH].number_of_nodes() == 0):
            newNodeStartID = 0
          else:
            newNodeStartID = max(macroInfo[m][GRAPH].nodes()) + 1
          for oldNode in macroInfo[m][BACTTYPELIST][macroBactType]:
            # only then can you reproducen
            # print oldNode
            # print macroInfo[m][BACTEATAGE]

            if (step_count - macroInfo[m][BACTEATAGE][oldNode] > BACT_IN_MACRO_REPRO_AGE):
              for newNodePerOldNode in xrange(BACT_IN_MACRO_REPR[macroBactType]):
                newMacBactNode = newNodeStartID
                newNodeStartID += 1
                macroInfo[m][GRAPH].add_node(newMacBactNode)
                macroInfo[m][BACTEATAGE][newMacBactNode] = step_count
                randTheta = 2 * np.pi * np.random.random_sample()
                randRadius = MACRO_BACT_WITHIN_RAD #* np.random.random_sample()
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
                #print "adding color " + str(bact_colors[b])
                macroInfo[m][GRAPHCOLORMAP] += [bact_colors[macroBactType]]


      for b in xrange(BACT_TYPE_COUNT):
        if (bact_count[b] >= BACT_COUNT_LIMIT[b]):
          continue
        bg = bactGraphs[b]
        bgn = bg.number_of_nodes();

        max_node_number = -1
        nodeList = bg.nodes()
        for i in nodeList:
          max_node_number = max(max_node_number, i)

        # we want to double this number
        for i in xrange(bgn):
          bactGraphs[b].add_node(i + max_node_number + 1)
        
        # this is a dictionary
        add_all_coords = bact_all_coords_map[b]
        this_bact_all_coords = dict()
        adding_node_number = 1
        for node in add_all_coords.keys():
          value1 = add_all_coords[node]
          value2 = np.empty_like(value1)
          value2[:] = value1
          value2[0] = np.clip(np.random.uniform(np.add(value1[0],-BACT_CHILD_DIS), 
            np.add(value1[0], 0.05)), 0, 0.999999)
          value2[1] = np.clip(np.random.uniform(np.add(value1[1],-BACT_CHILD_DIS), 
            np.add(value1[1],0.05)), 0, 0.999999)
          new_node = adding_node_number + max_node_number
          adding_node_number = adding_node_number + 1
          this_bact_all_coords[new_node] = value2
          this_bact_all_coords[node] = value1
          
        bact_all_coords_map[b] = this_bact_all_coords
        # they start off at original color
        
        bactPosInfo = copy.deepcopy(bactPosInfoOrig)
        # need to update this here, because we need this number for the 
        # AI movement immediately in this step
        for node in bact_all_coords_map[b].keys():
          posArray = bact_all_coords_map[b][node]
          x = posArray[0]
          y = posArray[1]
          c = int(x / AIC_CFRAC)
          r = int(y / AIC_RFRAC)
          if (r >= 2):
            r = 1
          if (c >= 2):
            c = 1
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
