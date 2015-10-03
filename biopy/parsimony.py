## This file is part of biopy.
## Copyright (C) 2015 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

"""
=============================================
Parsimony stuff
=============================================
"""

from __future__ import division

from collections import Counter
from treeutils import getPostOrder

__all__ = ["parsimonyScore", "parsimonyScoreWithPolytomies", "parsimonyChanges"]

def parsimonyScore(tree, taxVals) :
  """ Parsimony score on rooted binary tree. taxVals is a
  function returning the (text) value for a tip ('?' for missing)."""
  
  for itx in tree.get_terminals() :
    n = tree.node(itx)
    v = taxVals(n.data)
    n.data.vals = set([v]) if v != '?' else set()

  pscore = 0
  pord = getPostOrder(tree)
  for n in pord :
    if n.succ :
      ch = [tree.node(x) for x in n.succ]
      v = [c.data.vals for c in ch]
      a = reduce(set.intersection, v)
      if a :
        vals = a
      else :
        vals = reduce(set.union, v)
        tot = sum([len(x) > 0 for x in v])
        if tot > 1 :
          assert len(v) == 2
          pscore += 1
      n.data.vals = vals
      
  for i in tree.all_ids() :
    n = tree.node(i)
    del n.data.vals
    
  return pscore


def cost(parValPerAssignment, parentValue) :
  if parentValue == '?' :
    if len(parValPerAssignment) == 1 :
      assert parValPerAssignment.keys() == ['?']
      v = parValPerAssignment['?']
    else :
      v = min([parValPerAssignment[x] + 1 for x in parValPerAssignment if x != '?'])
      v = min(v, parValPerAssignment['?'])
  else :
    v = parValPerAssignment[parentValue] if parentValue in parValPerAssignment else parValPerAssignment['?']
    for x in parValPerAssignment:
      if x != '?' and x != parentValue:
        v = min(v,parValPerAssignment[x] + 1)
  return v
  
def parsimonyScoreWithPolytomies(tree, taxVals, keepInternalInfo = False, parsimonyStats = False) :
  """ Parsimony score for non-binary trees. Missing data is marked with a '?'.
  taxVals is a function which returns the tip value based on node data.
  """
  for itx in tree.get_terminals() :
    n = tree.node(itx)
    v = taxVals(n.data)
    n.data.vals = {v : 0}
    n.data.character = v
    if v != '?' :
      n.data.vals['?'] = 1
  
  pord = getPostOrder(tree)
  for n in pord :
    if n.succ :
      ch = [tree.node(x) for x in n.succ]
      vch = [c.data.vals for c in ch]
      allv = reduce(set.union,[set(x.keys()) for x in vch])
      v1 = dict([(x,sum([cost(v, x) for v in vch])) for x in allv - set('?')])
      v1['?'] = sum([cost(v,'?') for v in vch])
      n.data.vals = v1
  vals = tree.node(tree.root).data.vals
  pval = min([vals[x] for x in vals if x != '?']) if len(vals) > 1 else 0

  if parsimonyStats :
    # CI #minChanges / #actual
    if pval > 0 :
      c = Counter([tree.node(itx).data.character for itx in tree.get_terminals()])
      if '?' in c:
        c.pop('?')
      minChanges = len(c) - 1;    assert x >= 0
      maxChanges = sum(c.itervalues()) - c.most_common(1)[0][1] 
      pval = (pval, (minChanges, maxChanges))
    else :
      pval = (pval, None)
      
  if not keepInternalInfo :
    for i in tree.all_ids() :
      n = tree.node(i)
      del n.data.vals
    
  return pval
  
def _parsimonyChangesSub(tree, rn, ch) :
  vals = rn.data.vals
  if not rn.succ :
    if len(vals) == 1 :
      assert vals.keys()[0] == '?'
      rn.data.ch = '?'
    else :
      rn.data.ch = tuple(set(vals.keys()) - set('?'))[0]
  else :
    minc = min(vals.values())     
    if ch in vals and vals[ch] == minc or \
      ch not in vals and vals['?'] == minc :
      nch = ch
    else :
      u = [c for c,v in vals.iteritems() if c != '?' and v == minc]
      if len(u) == 0:
        assert vals['?'] == minc, (rn.id,ch,vals)
        nch = '?'
      else :
        nch = u[0]
    rn.data.ch = nch
    
    for cd in rn.succ :
      n = tree.node(cd)
      _parsimonyChangesSub(tree, n, nch)
      
def parsimonyChanges(tree) :
  # Had parsimony run on it, with keep 
  rn = tree.node(tree.root)
  ps,ch = min([[y,x] for x,y in rn.data.vals.iteritems()])
  rn.data.ch = ch
  for cd in rn.succ :
    _parsimonyChangesSub(tree, tree.node(cd), ch)

def nonInformative(tree, vals) :
  ps = parsimonyScoreWithPolytomies(tree, lambda x : vals[x.taxon], True)
  
  c = Counter([x for x in vals.itervalues() if x != '?'])
  if len(c) <= 1 or (len(c) == 2 and '?' in c) :
    # one known value or none
    return True
  mc = c.most_common()
  mcv = mc[0][0]
  
  for n in getPostOrder(tree) :
    if not n.succ :
      n.data.nmcv = int(vals[n.data.taxon] not in (mcv,'?'))
    else :
      nvals = n.data.vals
      pv = min(nvals.values()) if len(nvals) > 0 else 0
      n.data.nmcv = sum( [tree.node(c).data.nmcv for c in n.succ] )
      # mcv gives parsimony minimum
      if n.id != tree.root:
        if not (nvals[mcv] if mcv in nvals else nvals['?']) == n.data.nmcv:
          return False
      else :
        if not (pv == 0 or (nvals.get(mcv) == pv)) :
          return False
        if not ( pv == n.data.nmcv or (mcv not in nvals and n.data.nmcv == 1) ) :
          return False
  return True
