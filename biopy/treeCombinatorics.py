## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" Combinatorial utilities for trees.

Since in combinatorics we are not concerned with branch lengths, trees have a
simple representation. A tree is either a leaf - A list of length one - or a
list of trees. (Yes this prohibits internal nodes with only one descendant).

"""

from __future__ import division

from combinatorics import nPairs as c2, factorial as fac

# Various counting formulas on number of trees/forests under special conditions

def nLabeledHistories(n, k) :
  """ Total number of ranked histories of (labeled) C{k}-forests with total of
  C{n} tips.""" 
  assert n >= k
  
  return prod([long(c2(x)) for x in range(n, k, -1)])


def nRankedHistoriesCladeOnlyRoot(n, c) :
  """ Number of ranked histories of a labeled tree with a total of C{n+c} taxa, where
  the common ancestor of taxa from C{c} taxa is the root.
  """
  
  return prod([k*(k+1)//2L for k in range(2,c)]) * \
         prod([k*(k+3)//2L for k in range(c-1,n+c-1)])

def nRankedCombineWithClade(nc, no) :
  """ Number of ranked ways C{no} taxa can be combined with a C{nc} taxa tree into a
  tree.
  """
  return prod([k*(k+nc) for k in range(2,no+1)]) // 2L**(no-1) \
         if no > 0 else 1

def nForestRankedTopologies(n1, n2, k) :
  """ Number of ranked histories of a labeled forest with C{1+k} trees, the
  first tree on C{n1} taxa, the remaining C{k} on a total of C{n2} taxa.
  """
  
  return choose(n1+n2-k-2, n1-2) * prod([long(c2(x)) for x in range(3, n1+1) + range(k+1, n2+1)])

def nTreeRankedTopologies(n1, n2, k) :
  """ Number of ranked histories of a labeled tree with C{n1+n2} taxa, which
  contains a monophyletic clade on the first C{n1} taxa, and at the root of that
  clade there are exactly C{k+1} surviving lineages in the tree.
  """
  assert k <= n2
  
  return nr(k+1) * nForestRankedTopologies(n1, n2, k)

# counting for a specific tree

def _LHsize(p) :
  """ Number of labeled ranked histories of trees with the same unranked
  topology as C{p}.

  Result is returned indirectly as a pair (n,s), where n is the number of
  internal nodes and s is the number of ranked topologies which map to a single
  unranked topology.  Therefore the result is n!/s.
  """

  # leaf - no internal nodes, 1 to 1 mapping
  if len(p) == 1 :
    return (0,1)

  # get counts from all trees in forest
  q = [_LHsize(x) for x in p]
  
  #total number of internal nodes in tree
  s = sum([x[0] for x in q])+1

  #s! ways to arrange, removed multiple counting
  return s, s * prod([x[1] for x in q])

def labeledHistoriesOfTree(p) :
  """ Number of labeled ranked histories of trees with the same unranked
  topology as C{p}."""
  
  n,k = _LHsize(p)
  return fac(n)//k

def numberOfLabeledForests(forest) :
  """ Number of labeled ranked forests having the same unranked topology as
  C{forest}."""

  # number of ranked possibilities for each tree
  z = [LHsize(x) for x in forest]

  # prune single lineages 
  u = [(n,k) for n,k in z if n > 0]

  # same logic as for a single tree
  n,k = sum([x[0] for x in u]), prod([x[1] for x in u])
  return fac(n)//k


def allCompatibleLabeledTrees(forest, compat) :
  """ Iterator returning all compatible trees constructed by joining up a subset
  of the trees in C{forest}.

  Assumes all trees in forest are compatible.
  """

  t0 = forest[:1]
  yield (t0[0], forest[1:])
  
  if len(forest) > 1:
    # all trees not containing first taxon
    for t,f in allCompatibleLabeledTrees(forest[1:], compat) :
      yield (t,t0 + f)

    # all trees containing first taxon
    for t,f in allCompatibleNonTrivialWithFirstTaxaTrees(forest, compat) :
      yield (t,f)
      
  
def allCompatibleNonTrivialWithFirstTaxaTrees(forest, compat) :
  """ Iterator returning all compatible trees containing the first
  tree from C{forest} and a non empty subset of the remaining trees.

  Returns a pair each time, with the tree and the remaining taxa in the list.
  
  Assumes all trees are compatible.  
  """

  if len(forest) > 1:
    t0 = forest[:1]
    t00 = t0[0]
    # all trees not containing first taxa
    for t,f in allCompatibleLabeledTrees(forest[1:], compat) :
      t1 = [t00, t]
      if compat is None or compat(t1) :
        yield (t1,f)
        for x,f1 in allCompatibleNonTrivialWithFirstTaxaTrees([t1] + f, compat) :
          yield (x,f1)
          
def allCompatibleLabeledHistories(taxa, compat) :
  """ Iterator returning all forests whose leaves are C{taxa}. Non-compatible
  trees are pruned: C{compat} is a function which returns True for a compatible
  tree, False otherwise.

  Each taxon is a list of length 1.

  Returns forests, each forest is a list of trees whose taxa is C{taxa}.
  """
  
  n = len(taxa)
  if n == 0:
    yield []
  elif n == 1 :
    yield taxa
  else :
    # all combinations of a non-trivial tree (2 or more taxa) containing the
    # first taxon, with all possible forests on the remaining taxa
    for t,f in allCompatibleNonTrivialWithFirstTaxaTrees(taxa, compat):
      for f1 in allCompatibleLabeledHistories(f, compat) :
        yield [t] + f1

    # all forests where first taxa is in it's own tree
    for f1 in allCompatibleLabeledHistories(taxa[1:], compat):
      yield taxa[:1] + f1
      




def allBinarySplits(lst, empty=False) :
  """ Iterator returning all ways to partition elements in C{lst} into two
  groups.

  With C{empty}, allow an empty group, otherwise groups has to be non-empty.
  """
  
  assert len(lst) >= 2
  
  if empty :
    yield [lst, []]
    
  if len(lst) == 2 :
    yield [lst[0:1], lst[1:2]]
  else :
    yield [lst[:1], lst[1:]]
    for r,l in allBinarySplits(lst[1:]) :
      yield [r + lst[:1], l]
      yield [r, l + lst[:1]]      

def allTrees(taxa) :
  """  Iterator returning all labeled trees with tips labeled by taxa. """
  
  if len(taxa) == 1 :
    yield [taxa[0]]
  elif len(taxa) == 2 :
    yield [taxa[0:1], taxa[1:2]]
  else :
    for l,r in allBinarySplits(taxa) :
      for x in allTrees(l):
        for y in allTrees(r):
          yield [x,y]

def allLabeledPartitions(lst) :
  """ Iterator returning all ways to partition elements in C{lst} into 1 to
  len(lst) groups.

  Returns a list of lists each time.
  """
  
  n = len(lst)
  if n :
    yield [lst]
    
  if n > 2 :
    for l,r in allBinarySplits(lst[1:], 1) :
      for p in allLabeledPartitions(r) :
        yield [lst[:1] + l] + p
      for p in allLabeledPartitions(l) :
        yield [lst[:1] + r] + p
  elif n == 2:
    yield [lst[:1],lst[1:]]

def allForestTrees(forest) :
  if not len(forest) :
    yield []
  else :
    for t in allTrees(forest[0]) :
      for f in allForestTrees(forest[1:]):
        yield [t] + f

def allLabeledForests(taxa) :
  """ Iterator returning all labeled forests, where each taxon from taxa appears
  in exactly one tree.
  """
  
  for p in allLabeledPartitions(taxa) :
    for f in allForestTrees(p) :
      yield f

#  LocalWords:  monophyletic clade unranked combinatorics
