#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2015 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import argparse, sys, os.path

parser = argparse.ArgumentParser(description= """ %(prog)s [OPTIONS] tree-nex-file""")

parser.add_argument("-r", "--ref", metavar="FILE",
                    help="""""")
parser.add_argument("-t", "--th", metavar="F",
                    default = 0.9, help="""""")
parser.add_argument("--markref", action = "store_true", default = False, help="""""")

parser.add_argument("--taxonomy-sep", dest="taxsep", metavar="CHAR",
                    default = ';', help="""(default %default)""")

parser.add_argument('treefile', metavar='FILE', help="tree as NEXUS file")

options = parser.parse_args()

th = float(options.th)
assert 0 < th < 1

markref = options.markref

from collections import Counter, defaultdict

from biopy.genericutils import flatten
from biopy.treeutils import *
from biopy.INexus import INexus

def stripSize(tx) :
  tx = tx.strip("'")
  if tx[-1] == ';' :
    tx = tx[:-1]
  p = tx.split(';')
  sz = None
  if p[-1].startswith('size=') :
    sz = int(p[-1][5:])
    p = p[:-1]
  return ';'.join(p),sz
  
taxonomy = dict()
taxonomyLength = None
if options.ref :
  for l in file(options.ref) :
    l = l.strip()
    if l and l[0] != '#' :
      nm, tx = l.split('\t')
      v = tx.split(options.taxsep)
      assert taxonomyLength is None or taxonomyLength == len(v)

      if nm in taxonomy:
        assert taxonomy[nm] == v
      else :
        taxonomy[nm] = v
        taxonomyLength = len(v)
else :
  parser.show_help()
  
for tree in INexus(simple=1).read(options.treefile) :
  for n in getPostOrder(tree):
    d = n.data
    if not n.succ :
      tx = d.taxon
      tx, sz = stripSize(tx)
      #import pdb; pdb.set_trace()
      if tx in taxonomy:
        t = [Counter({x.strip() : sz}) if x else Counter() for x in taxonomy[tx]]
        d.taxonomy = t # list(reversed(t))
        #d.needs = False
      else :
        d.taxonomy = None
        #d.needs = True
    else :
      ch = [tree.node(c).data for c in n.succ]
      if ch[0].taxonomy is None :
        d.taxonomy = ch[1].taxonomy
      elif ch[1].taxonomy is None :
        d.taxonomy = ch[0].taxonomy
      else :
        d.taxonomy = [Counter(x) for x in ch[0].taxonomy]
        for x,y in zip(d.taxonomy, ch[1].taxonomy) :
          x.update(y)
      #d.needs = ch[0].needs or ch[1].needs
  
  ca = CAhelper(tree)

  for n in getPostOrder(tree):
    d = n.data
    if not n.succ :
      tx = d.taxon
      if d.taxonomy is None :
        d.needs = [n]
        d.taxasgn = None
      else :
        d.needs = [] 
    else :
      chd = [tree.node(x).data for x in n.succ]
      needs = [x.needs for x in chd]
      d.needs = flatten(needs)
      if not(d.needs) :
        continue
      if not all([x.taxonomy for x in chd]) :
        continue
      assert d.taxonomy

      for k in range(len(d.taxonomy)-1, -1,-1) :
        c = d.taxonomy[k]
        tot = sum(c.values())
        v = max(c, key = c.get)
        if v and c[v] >= tot * th and all([v in x.taxonomy[k] for x in chd]) :
          d.taxlevel = k
          if 0 : print "assigning",k,v, "to", len(d.needs), d.taxonomy
          
          #import pdb ; pdb.set_trace()
          for x in d.needs :
            #assert x.data.taxasgn is None or x.data.taxasgn == n
            if x.data.taxasgn is None :
              x.data.taxasgn = n
            else :
              if x.data.taxasgn.data.taxlevel < k :
                #import pdb ; pdb.set_trace()
                x.data.taxasgn = n 
          if k == len(d.taxonomy)-1 :
            d.needs = []
          break
          
  getx = lambda c : max(c, key = c.get)
  def getTaxStr(d) :
    return [getx(c) for c in d.taxonomy[:d.taxlevel+1]]

  for tp in tree.get_terminals() :
    nn = tree.node(tp)
    if hasattr(nn.data, "taxasgn") and nn.data.taxasgn:
      assert nn.data.taxonomy is None
      
      d = nn.data.taxasgn.data

      if nn.data.taxon[0] != "'" :
        nn.data.taxon = "'" + nn.data.taxon + "'"
      nn.data.taxon = "'" + '; '.join(reversed(getTaxStr(d))) + "@" + nn.data.taxon[1:]
    else :
      if markref and nn.data.taxonomy is not None :
        if nn.data.taxon[0] != "'" :
          nn.data.taxon = "'" + nn.data.taxon + "'"
        nn.data.taxon = "'" + '; '.join(reversed([x.keys()[0] for x in nn.data.taxonomy if len(x)])) + "$" + nn.data.taxon[1:]
      
  labels = []
  for k,x in enumerate(tree.get_terminals()) :
    labels.append(tree.node(x).data.taxon)
    tree.node(x).data.taxon = str(len(labels))

  tlog = TreeLogger(tree.name + ".taxc.trees",
                  labels = zip(range(1,len(labels)+1), labels), overwrite=1)
  tlog.outTree(tree)
  tlog.close()

