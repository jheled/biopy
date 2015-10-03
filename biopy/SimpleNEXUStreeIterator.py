## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
from __future__ import division

" Non-memory intensive. "

import shlex
from parseNewick import parseNewick
from genericutils import fileFromName

# Assume no nexted blocks
# trees on separate lines, one per line

class SimpleNEXUSstrIter(object) :
  def __init__(self, instream) :
    self.stream = fileFromName(instream) if isinstance(instream,str) else instream
    self._grabLine()
    while self.nextline.strip() != '#NEXUS' :
      self._grabLine()

  def _grabLine(self) :
    try :
      self.nextline = next(self.stream)
    except StopIteration:
      self.iline = None
      self.nextline = None
      return
      
    # StopIteration will be passed up
    if self.nextline.startswith('tree ') :
      self.iline = None
    else :
      self.iline = shlex.shlex(self.nextline)
    
  def __iter__(self):
    return self

  def next(self):
    if self.iline :
      try :
        tok = next(self.iline)
        return (None,tok)
      except StopIteration :
        self._grabLine()
        return self.next()
    #else - a tree
    if self.nextline :
      t = self.nextline[5:]
      self._grabLine()
      return ('tree', t)

    raise StopIteration
    
def _parseOptions(txt) :
  assert txt[:2] == "[&"
  
  e = t.index(']', 2)
  opts = t[2:e]
  od = dict()
  for o in opts.split(','):
    nv = o.split('=')
    if len(nv) == 2:
      od[nv[0]] = nv[1]
    elif len(nv) == 1 :
      od[nv] = None
          
  return od, txt[e+1:].strip()
  
PUNCTUATION='()[]{}/\,;:=*\'"`+-<>'
WHITESPACE=' \t\n'
WANTQ = set(PUNCTUATION + WHITESPACE)

class SimplisticNEXUStreeIterator(object) :
  def __init__(self, instream, withAttributes=True) :
    self.itokens = SimpleNEXUSstrIter(instream)
    self.withAttributes = withAttributes
    
    while True :
      tok = ''
      while tok.lower() != 'begin' :
        tok = self._nextSimpleTok()
      tok = self._nextSimpleTok()

      if tok.lower() == 'trees' :
        break
        
      # skip block: no nested blocks allowed
      while True:
        while tok.lower() != 'end' :
          tok = self._nextSimpleTok()
        tok = self._nextSimpleTok()
        if tok == ';' :
          break
    
    if tok.lower() == 'trees' :
      tok = self._nextSimpleTok() ; assert tok == ';'
      tok = next(self.itokens)
      self.taxatable = None
      if tok[0] is None and tok[1].lower() == 'translate' :
        self.taxatable = dict()

        while True :
          tfrom = self._nextSimpleTok()
          if tfrom == ';' :
            break
            
          tto = self._nextSimpleTok()
          if set(tto).intersection(WANTQ) :
            tto = tto.replace("'","''")
            tto= "'" + tto + "'"
          
          self.taxatable[tfrom] = tto
          
          tok = self._nextSimpleTok()
          if tok == ';' :
            break
          if tok != ',' :
            raise RuntimeError("parse error")
            
        self.treetok = next(self.itokens)
      else :
         self.treetok = tok
    else :
      self.treetok = None
      
  def _nextSimpleTok(self) :
    tok = next(self.itokens)
    if not tok:
      raise RuntimeError("parse error: unexpected text")
    if tok[0] is None :
      return tok[1]

    raise RuntimeError("unexpected tree")
    
  def  __iter__(self):
    return self

  def next(self) :
    if self.treetok :
      if self.treetok[0] == 'tree' :
        t = self.treetok[1]
        self.treetok = next(self.itokens)

        t = t.strip().split('=')
        assert len(t) >= 2
        tname = t[0].strip()
        
        t = '='.join(t[1:]).strip()
        rooted=False
        weight=1.0
        
        if t[0] == '[' :
          o,t = _parseOptions(t)
          if 'R' in o :
            rooted = True
          if 'U' in o :
            rooted = False
          if 'W' in o :
             weight=float(o['W'])

        tree = parseNewick(t.strip(), weight=weight, rooted=rooted, name=tname.split()[0],
                           loadAttributes = self.withAttributes)
        if self.taxatable :
          for n in tree.get_terminals():
            data = tree.node(n).data
            try:
              data.taxon = self.taxatable[data.taxon]
            except (ValueError,KeyError):
              raise RuntimeError("translation failed")
          
        return tree
        
    raise StopIteration
  
