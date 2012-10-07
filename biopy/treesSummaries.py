## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

"""
Summary Tree from posterior trees
=================================
"""

__all__ = ["summaryTreeUsingTaxonPairs", "summaryTreeUsingMedianHeights", "taxaPartitionsSummaryTree"]

from numpy import median, mean
import copy

from treeutils import getTreeClades, nodeHeights, getPostOrder, toNewick
from treeMeasure import allPartitions
from parseNewick import parseNewick
from treePlotting import taxaDistance
import mau

# Should make this more efficient than O(#clades^2)

def summaryTreeUsingTaxonPairs(tree, xtrees, atz = False) :
  tree = copy.deepcopy(tree)
  
  sClades = getTreeClades(tree)
  cladeSets = [frozenset(c) for c,n in sClades]
  stats = [[0.0,0.0] for x in cladeSets]
  NS = 0

  for t in xtrees:
    nhs = nodeHeights(t, allTipsZero = atz)
    h = [0.0]*len(cladeSets)
    for c,n in getTreeClades(t) :
      c = frozenset(c)
      for k,x in enumerate(cladeSets) :
        if x.issubset(c) :
          h[k] = nhs[n.id]

    NS += 1
    for st,hk in zip(stats,h) :
      st[0] += hk
      st[1] += hk**2

  for (c,n),s in zip(sClades,stats) :
    n.data.hh = s[0]/NS

  # use t, nhs from last iteration
  for tx in tree.get_terminals():
    if atz :
      h1 = 0.0
    else :
      n2 = t.search_taxon(tree.node(tx).data.taxon)
      h1 = nhs[n2]
    tree.node(tx).data.hh = h1

  for nid in tree.all_ids() :
    if nid != tree.root :
      n = tree.node(nid)
      ph = tree.node(n.prev).data.hh
      h = n.data.hh
      assert ph >= h
      n.data.branchlength = ph - h
  return tree


def summaryTreeUsingMedianHeights(tree, xtrees) :
  tree = copy.deepcopy(tree)

  func = lambda t,(n,h) : h
  posteriorParts,rhs = allPartitions(tree, xtrees, func = func,
                                     withHeights = True, withRoot = True)
  treeParts = allPartitions(tree, [tree])

  for r,k in enumerate(treeParts) :
    # Node id
    nn = treeParts[k][0][1]
    
    if k in posteriorParts :
      tree.node(nn).data.height = median(posteriorParts[k])
    else :
      raise RuntimeError("tree incompatible with trees")

  # Assume all trees share same tip heights (not checked)
  nh = nodeHeights(xtrees[0])
  
  tree.node(tree.root).data.height = median(rhs)
  for n in getPostOrder(tree):
    if not len(n.succ) :
      n.data.height = nh[xtrees[0].search_taxon(n.data.taxon)]
    else :
      # Make sure node is heigher than descendants 
      n.data.height = max([n.data.height] + [tree.node(x).data.height for x in n.succ])

  for n in tree.all_ids() :
    node = tree.node(n)
    if node.prev is not None:
      p = tree.node(node.prev) 
      node.data.branchlength = p.data.height - node.data.height
      assert node.data.branchlength >= 0
        
  return tree


from collections import Counter, defaultdict
from combinatorics import allPairs
from numpy import mean
import itertools

# find a linear (on a line) order of tree taxa which is compatible with
# as many trees (or parts of) as possible (compatible: can be ontained by
# rotating internal nodes).

def _taxOrder(trees) :
  dtops = Counter(toNewick(tree, topologyOnly=1) for tree in trees)

  averagedPosteriorDistances = defaultdict(lambda : 0)

  for tx in dtops:
    tree = parseNewick(tx)
    distances = dict()
    taxaDistance(tree, tree.root, distances)
    for k in distances:
      averagedPosteriorDistances[k] += distances[k] * dtops[tx]/len(trees)

  taxa = sorted(trees[0].get_taxa())
  ntax = len(taxa)
  for j,k in allPairs(taxa) :
    averagedPosteriorDistances[k,j] = averagedPosteriorDistances[j,k]
    
  def groupDistance(g1,g2,dm) :
    return mean([dm[i] for i in itertools.product(g1,g2)])
  ## def groupDistance(g1,g2,dm) :
  ##   return min([dm[i] for i in itertools.product(g1,g2)])
  
  groups = [[x] for x in taxa]
  dm = averagedPosteriorDistances
  while len(groups) > 1:
    # find the two closest groups.
    dists = [(groupDistance(groups[j],groups[k],dm), (j,k)) 
             for j,k in allPairs(range(len(groups)))]
    dists = sorted(dists)
    d,(ij,ik) = dists[0]
    # 123 abc 0,0 0,1 1,0 02 11 20
    # 321 abc
    # 123 cba
    # 321 cba
    # abc 123
    # abc 321
    # cba 123
    # cba 321
    g1,g2 = groups[ij],groups[ik]
    def gid(g1,g2,dm) :
      d = []
      for n in range(len(g1)+len(g2) - 1) :
        for i in range(-1,max(-len(g1)-1,-(len(g2)-n+1)), -1) :
          #i > -(len(g2)-n+1)
          d.append(dm[g1[i], g2[n-(i+1)]])
      return d

    dis = gid(g1,g2,dm),gid(g1,list(reversed(g2)),dm),gid(g2, g1,dm), gid(list(reversed(g2)),g1,dm)
    dis = sorted(zip(dis,range(4)))  
    o = dis[0][1]
    if o & 1 :
      g2 = list(reversed(g2))
    if o > 1 :
      g1,g2 = g2,g1
    groups[ij] = g1 + g2
    del groups[ik]
  otaxa = groups[0]
  return otaxa

def _getCladeEstimates(tree, txp, hs) :
  po = getPostOrder(tree)
  for n in po:
    data = n.data
    if data.taxon :
      data.clade = ([txp[data.taxon]], True)
      data.pheight = data.branchlength
    else :
      sc = [tree.node(s).data.clade for s in n.succ]
      svalid = all([s[1] for s in sc])
      sc = [s[0] for s in sc]
      nh = max([tree.node(s).data.pheight for s in n.succ])
      mn,mx = min([min(x) for x in sc]),max([max(x) for x in sc])
      if mx - mn + 1 == sum([len(x) for x in sc]) :
        valid = True
      else :
        valid = False
      if valid and svalid:
        m1 = max(sc[0])
        if not (m1 < mx) :
          m1 = max(sc[1]) # [max(x) for x in sc[1]]
        assert (0 <= m1 < mx)
        hs[m1].append(nh)
      if n.id != tree.root :
        data.pheight = nh + data.branchlength
        data.clade = (sc[0] + sc[1], valid)
        
def taxaPartitionsSummaryTree(trees, summaryType = "median") :
  if summaryType not in ["mean", "median", "both"] :
    raise ValueError("summaryType should be one of mean, median, both")
  
  torder = _taxOrder(trees)
  hs = [list() for k in torder[:-1]]

  txp = dict(zip(torder,itertools.count()))
  for t in trees :
    _getCladeEstimates(t, txp, hs)

  if summaryType in ["median","both"] :
    md = mau.mau2Tree((torder, [median(x) if len(x) else 0  for x in hs],
                             [(0,False) for x in hs]))
  if summaryType in ["mean","both"] :
    mn = mau.mau2Tree((torder, [mean(x) if len(x) else 0 for x in hs],
                             [(0,False) for x in hs]))
    
  return (mn,md) if summaryType=="both" else mn if summaryType == "mean" else md


def taxaPartitionsSummaryTreeX(trees, summaryType = "median") :
  if summaryType not in ["mean", "median", "both"] :
    raise ValueError("summaryType should be one of mean, median, both")
  
  torder = _taxOrder(trees)
  hs = [list() for k in torder[:-1]]

  txp = dict(zip(torder,itertools.count()))
  for t in trees :
    _getCladeEstimates(t, txp, hs)

  if summaryType in ["median","both"] :
    md = mau.mau2Tree((torder, [median(x) if len(x) else 0  for x in hs],
                             [(0,False) for x in hs]))
  if summaryType in ["mean","both"] :
    mn = mau.mau2Tree((torder, [mean(x) if len(x) else 0 for x in hs],
                             [(0,False) for x in hs]))
    
  return ((mn,md) if summaryType=="both" else (mn,) if summaryType == "mean" else
          (md,)) + ((torder, hs))
