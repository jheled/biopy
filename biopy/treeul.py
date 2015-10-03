## This file is part of biopy.
## Copyright (C) 2013 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

"""
========================================
Non-Ultrametric/Time trees to Time Trees
========================================
"""

from __future__ import division

import scipy.optimize
import copy

from numpy import var,mean,array
from biopy.treeutils import getPostOrder, getPreOrder, nodeHeights

__all__ = ["UMTree"]

def treeHeightEstimate(t) :
  """ Distance between root to mean tip age. """
  nh = nodeHeights(t, allTipsZero=False)
  hm = mean([nh[c] for c in t.get_terminals()])
  return nh[t.root] - hm

  
def nodeRateMultiplier(t, n, hTarget, nh) :
  rt = t.node(t.root)
  #c = rt.succ[0]
  
  b,m = n.data.branchlength, mean([nh[x] for x in t.get_taxa(n.id, asIDs=1)])
  gt = b + (nh[n.id] - m)
  #print n.id,b,gt,b/gt, hTarget * b/gt
  n.data.timebranchlength = hTarget * b/gt

  hTarget -= n.data.timebranchlength
  for c in n.succ:
    nodeRateMultiplier(t, t.node(c), hTarget, nh)

def nodeRateMultiplier1(t, n, hTarget, nh) :
  rt = t.node(t.root)
  #c = rt.succ[0]
  
  b, h = n.data.branchlength, nh[n.id]
  r = mean([b/(b+(h - nh[x])) for x in t.get_taxa(n.id, asIDs=1)])
  n.data.timebranchlength = hTarget * r

  hTarget -= n.data.timebranchlength
  for c in n.succ:
    nodeRateMultiplier1(t, t.node(c), hTarget, nh)

def calcNoderateMultipliers(t, hTarget = None, method = 1) :
  nh = nodeHeights(t,allTipsZero=0)
  if hTarget is None:
    hm = mean([nh[c] for c in t.get_terminals()])
    hTarget = nh[t.root] - hm
  meth = [nodeRateMultiplier,nodeRateMultiplier1][method]
  rt = t.node(t.root)
  for c in rt.succ :
    meth(t, t.node(c), hTarget, nh)
  return hTarget

def cmbn(hpexpr, i) :
  if hpexpr :
    hpexpr = ("x+y for x,y in zip([%s], hp%d)" % (hpexpr,i))
  else :
    hpexpr = "hp%d" % i
  return hpexpr

def _targetVar(r) :
  i = [1/x for x in r[:-4]] + list(r[-4:-2])
  return var(i)

def _targetVarDerivatives(r) :
  i = [1/x for x in r[:-4]] + list(r[-4:-2])

  n = len(i)
  si = sum(i)

  p = [((1./n) * 2 * i[k] - (2 * i[k] + 2 * (si - i[k]))/n**2)
       for k in range(n)] + [0,0]
  for k in range(n-2) :
    p[k] *= -i[k]**2
  return p

def _targetNormedVar(r) :
  i = [1/x for x in r[:-4]] + list(r[-4:-2])
  m = mean(i)
  return var(i)/(m*m)

  # fm partial derivatives (after some simplifications)
  # No need for C code, those calls are insignificant compared to cycle update
  # time.
  
def _targetNormedVarDerivatives(r) :
  i = [1/x for x in r[:-4]] + list(r[-4:-2])

  n = len(i)
  si = sum(i)
  n2 = n*n

  dvdr = [2 * ((n * i[k] - si)/n2) for k in range(n)]
  dm2dr = 2 * si/n2

  ig = n2/(si*si)
  fgp = var(i) * dm2dr * (ig*ig)
  p = [fp * ig - fgp for fp in dvdr]

  for k in range(n-2) :
    p[k] *= -i[k]**2
  return p + [0,0]

def _UMTree_slsqp(tin, targetHeight = None, useDerives = True,
                 doInit = False, normed = True,
                 niter=1000, r0x = None, verb = 0, internals = False) :
  t = copy.deepcopy(tin)
  
  nmap = []
  code = []
  primecode = []
  codeNodes = []
  # Number of branches 
  nr = 2*(len(t.get_terminals()) - 1)

  if targetHeight is None :
    targetHeight = treeHeightEstimate(t)
  assert float(targetHeight) == targetHeight and targetHeight > 0
  
  for n in getPostOrder(t) :
    if n.id != t.root :
      if n.prev != t.root :
        n.data.myindx = len(nmap)
        n.data.hexpr = "%.14f * r[%d]" % (n.data.branchlength, n.data.myindx)
        c = ["0"]*(nr+2)
        c[n.data.myindx] = "%.14f" % n.data.branchlength
        n.data.hpexpr = ",".join(c)
        nmap.append(n.id)
      else :
        n.data.hexpr = ""
        n.data.hpexpr = ""
      if n.succ :
        n.data.hexpr += ' + h%d' % n.succ[0]
        n.data.hpexpr = cmbn(n.data.hpexpr, n.succ[0])
        
      if n.prev != t.root and n.data.hexpr :
        n.data.chexpr = "h%d = %s" % (n.id, n.data.hexpr)
        code.append(n.data.chexpr)
        primecode.append("hp%d = [%s]" % (n.id, n.data.hpexpr))
        
    else :
      # left,right rates, left right branch : len(nmap)-2 to len(nmap)+1
      assert nr - 2 == len(nmap)
      
      lft = t.node(n.succ[0])
      lft.data.myindx = nr-2
      lft.data.hexpr = "r[%d] " % (nr) + lft.data.hexpr
      lft.data.chexpr = "h%d = %s" % (lft.id,lft.data.hexpr)

      sss = ",".join(["0"]*(nr) + ["1","0"])
      if lft.data.hpexpr :
        lft.data.hpexpr = "x1+y1 for x1,y1 in zip([%s], %s)" % \
                          (sss,lft.data.hpexpr)
      else :
        lft.data.hpexpr = sss
      
      rht = t.node(n.succ[1])
      rht.data.myindx = nr-1
      brr = (lft.data.branchlength+rht.data.branchlength)
      rht.data.hexpr = "r[%d] " % (nr+1) + rht.data.hexpr
      rht.data.chexpr = "h%d = %s" % (rht.id,rht.data.hexpr)

      code.append(rht.data.chexpr)
      code.append("ch1 = h%d - %.14f" % (rht.id, targetHeight))
      
      sss = ",".join(["0"]*(nr) + ["0","1"])
      if rht.data.hpexpr :
        rht.data.hpexpr = "x2+y2 for x2,y2 in zip([%s], %s)" % \
                          (sss, rht.data.hpexpr)
      else :
        rht.data.hpexpr = sss
      
      nmap.extend(n.succ)
      code.append(lft.data.chexpr)
      primecode.append("cph1 = hp%d = [%s]" % (rht.id, rht.data.hpexpr))
      primecode.append("cph2 = hp%d = [%s]" % (lft.id, lft.data.hpexpr))

      code.append("ch2 = h%d - %.14f" % (lft.id, targetHeight))

      brr = rht.data.branchlength + lft.data.branchlength
      code.append("crt = (r[%d] * r[%d] + r[%d] * r[%d]) - %.15f" % (nr-2,nr,nr-1,nr+1,brr))

      sss = ",".join(["0"]*(nr-2)) + ",r[%d],r[%d],r[%d],r[%d]" % (nr,nr+1,nr-2,nr-1)
      primecode.append("cprt = [%s]" % sss)
      
    if n.succ and n.id != t.root:
      code.append("c%d = h%d - h%d" % ((n.id,) + tuple(n.succ)))
      primecode.append("cp%d = [x-y for x,y in zip(hp%d, hp%d)]" % ((n.id,) + tuple(n.succ)))
      codeNodes.append(n.id)

  cd = ["def fx(r) :"]
  cd.extend(["  " + x for x in code])
  ccs = ["c%d" % k for k in codeNodes] + ['ch1','ch2','crt']
  hhs = ["h%d" % k for k in nmap]
  cd.append("  return [" + ",".join(ccs) + "], [" + ",".join(hhs) + "]")

  exec ( "\n".join(cd) ) in globals()

  cdp = ["def fxp(r) :"]
  cdp.extend(["  " + x for x in primecode])
  ccs = ["cp%d" % k for k in codeNodes] + ['cph1','cph2','cprt']
  cdp.append("  return [" + ",".join(ccs) + "]")

  exec ( "\n".join(cdp) ) in globals()

  # optimization target

  if r0x is not None :
    r0 = r0x
  else :
    if doInit:
      calcNoderateMultipliers(t, targetHeight)    
      cl,cr = t.node(t.root).succ
      btl,btr = (t.node(cl).data.timebranchlength , t.node(cr).data.timebranchlength)
      r0 = [t.node(c).data.timebranchlength/t.node(c).data.branchlength for c in nmap] + \
           [btl, btr]
      r0[-3] = 1/r0[-3]
      r0[-4] = 1/r0[-4]
    else :
      r0 = [1]*(nr) + [targetHeight*.1]*2

  # most are constant, can improve
  #derv = array(fxp(0))
  
  assert nr == len(nmap)
  
  slsqp = scipy.optimize.slsqp.fmin_slsqp
  if useDerives :
    fm,fmp = (_targetNormedVar,_targetNormedVarDerivatives) if normed \
             else (_targetVar, _targetVarDerivatives)
    
    re = slsqp(fm, r0,
              fprime = lambda x : array(fmp(x)),
              f_eqcons = lambda x : array(fx(x)[0]),
              fprime_eqcons = lambda x : array(fxp(x)),
              bounds = [(1e-10,10)]*nr + [(0,targetHeight)]*2,
              iter = niter, iprint=verb, full_output=1)
  else :
    re = slsqp(fm, r0,
               f_eqcons = lambda x : array(fx(x)[0]),
               bounds = [(1e-10,10)]*nr + [(0,targetHeight)]*2,
              iprint=verb, full_output=1)
  r = re[0]
  if re[3] != 0 :
    if not (re[3] == 8 and (abs(targetHeight - r[-2]) < 1e-9 or
                            abs(targetHeight - r[-1]) < 1e-9)) :
      ok = re[3] in [8,9] and all([abs(x) < 1e-8 for x in fx(r)[0]])
      if not ok :
        import pdb ; pdb.set_trace()
        raise RuntimeError(re[4])
  
  hs = fx(r)[1]
  for k,i in enumerate(nmap) :
    n = t.node(i)
    n.data.ph = hs[k] 
    n.data.subsbranchlength = n.data.branchlength
    n.data.attributes = {'clockrate' : r[n.data.myindx]}
    
    if n.succ :
      n.data.branchlength = n.data.ph - t.node(n.succ[0]).data.ph
    else :
      n.data.branchlength = n.data.ph
    n.data.branchlength = max(n.data.branchlength, 0)
    
  if internals :
    return t, fm(r), r, (fx, fxp, fm, fmp)
  return t, fm(r)



def der(k, n, h) :
  if h == -1 :
    return ",".join(["0"]*k + ["h0"] + ["0"]*(n-k-1))
                    
  return ",".join((['x[%d] * d_h%d_x[%d]' % (k,h,i) for i in range(k)] +
                   ["h%d" % h] +
                   ["0"]*(n-k-1)))


def _treeBranchAssignmentExprs(tree, hRoot,
                               withDerivative = False, withInit = False) :
  
  allid = set(tree.all_ids())  # taxa
  terms = set(tree.get_terminals())
  # internal nodes
  allint = allid - terms

  #  unrooted has dated tips when rooted somewhere
  nh = nodeHeights(tree, allTipsZero = False)
  # Node before its descendants
  nInOrder = getPreOrder(tree, includeTaxa=True)
  # Reversed, child before parent
  nhInOrder = [(nh[x],x) for x in reversed(nInOrder)]
  nInOrder = [x for x in nInOrder if x not in terms]

  nParams = len(nInOrder) - 1
  
  if withInit :
    htox = []
  
  sr = []
    
  sr.append("h%d = %.15g" % (nInOrder[0], hRoot))
    
  for i,k in enumerate(nInOrder[1:]):
    h = tree.node(k).prev
    
    sr.append("h%d = x[%d] * h%d" % (k, i, h))
    if withInit:
      htox.append("x[%d] = %.15g" % (i, nh[k]/nh[h]))
      
    if withDerivative:
      if h == tree.root :
        sr.append("d_h%d_x = [%s]" % (k,der(i, nParams, -1)))
      else :
        sr.append("d_h%d_x = [%s]" % (k,der(i, nParams, h)))
        
  for h,n in nhInOrder[:-1] :
    assert n != tree.root
    p = tree.node(n).prev
    if n in terms:
      sr.append("b%d=(h%d) # %d" % (n, p, n)) 
       
      if withDerivative:
        if p == tree.root :
          sr.append("d_b%d_x = [0]*%d" % (n,nParams))
        else :
          sr.append("d_b%d_x = d_h%d_x" % (n,p))
    else :
      sr.append("b%d=(h%d - h%d) # %d" % (n, p, n, n))
      
      if withDerivative:
        if p == tree.root :
          sr.append("d_b%d_x = [-v for v in d_h%d_x]" % (n,n))
        else :
          sr.append("d_b%d_x = [u-v for u,v in zip(d_h%d_x, d_h%d_x)]" % (n,p,n))
        
  if withInit:
    return sr, htox
  return sr


if 0:
  def _varianceAndDerive(bsnr, bsr, a, variant) :
    rl = a/bsr[0][0]
    rr = ((bsr[0][1] + bsr[1][1]) - a) / bsr[1][0]
    r = [bxs[1]/bxs[0] for bxs in bsnr] + [rl,rr]
    n = len(r)
    si = sum(r)
    ai = si/n
    
    val = sum([x**2 for x in r])/n - ai*ai
    if len(bsr[0]) == 2 :
      return val
    bs = bsnr + bsr
    
    dfdir = [0.0]*(n)
    f2 = 2./(n*n)
    for k in range(n) :
      dfdir[k] = f2 * (n*r[k] - si) * -(r[k]/bs[k][0])

    a = [bs[k][2] for k in range(n)]
    m = len(a[0])
    p = [0]*m
    for j in range(m) :
      p[j] = sum([dfdir[k] * a[k][j] for k in range(n)])

    #dfda = df/drl * drl/da + df/drr * drr/da
    pa = (f2 * (n*r[-2] - si)) * (1/bsr[0][0]) + (f2 * (n*r[-1] - si)) * (-1/bsr[1][0])
    p.append(pa)
    
    return val,p
  
  def _normedVarianceAndDeriv(bsnr, bsr, a) :
    rl = a/bsr[0][0]
    rr = ((bsr[0][1] + bsr[1][1]) - a) / bsr[1][0]
    r = [bxs[1]/bxs[0] for bxs in bsnr] + [rl,rr]

    n = len(r)
    si = sum(r)
    ai = si/n
    iai2 = 1/(ai*ai)
    
    val = (sum([x**2 for x in r])/n) * iai2 - 1

    if len(bsr[0]) == 2 :
      return val
    bs = bsnr + bsr

    # dvdr = [(n * r[k] - si) for k in range(n)]

    #fgp = val * dm2dr * iai2 # var(r) * dm2dr * (ig*ig) 
    #fgp = val * si
    #cof = 2 * iai2/(n*n)
    cof = 2/(si*si)
    # dfdr = [cof * (fp - fgp) for fp in dvdr]
    a1 = cof * n
    b1 = cof * si * (1 + val)
    #dfdr = [cof * ((n * r[k] - si) - fgp) for k in range(n)]
    dfdr = [a1 * r[k] - b1 for k in range(n)]

    dfdir = [dfdr[k] * -(r[k]/bs[k][0]) for k in range(n)]

    a = [bs[k][2] for k in range(n)]
    m = len(a[0])
    p = [0]*m
    for j in range(m) :
      p[j] = sum([dfdir[k] * a[k][j] for k in range(n)])

    pa = dfdr[-2]/bsr[0][0] - (dfdr[-1]/bsr[1][0])
    p.append(pa)

    return val,p

  def _normedVarianceAndDeriv(bsnr, bsr, a, variant) :
    
    rl = a/bsr[0][0]
    rr = ((bsr[0][1] + bsr[1][1]) - a) / bsr[1][0]
    r = [bxs[1]/bxs[0] for bxs in bsnr] + [rl,rr]

    n = len(r)
    si = sum(r)
    ai = si/n

    if variant == 0 :
      iai2 = 1/(ai*ai)
      val = (sum([x**2 for x in r])/n) * iai2 - 1
    
    if variant == 1 :
      val = sum([(x/ai - 1)**2 for x in r])**0.5 + (ai - 1)**2

    if variant == 2 :
      val = sum([(x/ai - 1)**2 for x in r]) + (ai - 1)**2

    if variant == 3 :
      val = sum([(x/ai - 1)**2 for x in r]) + abs(ai - 1)
      
    if variant == 4 :
      val = sum([(x/ai - 1)**2 for x in r])**0.5 + abs(ai - 1)
    
    if variant == 5 :
      val = sum([(x/ai - 1)**2 for x in r]) + n * (ai - 1)**2

    if variant == 6 :
      val = sum([(x/ai - 1)**2 for x in r]) + n**2 * (ai - 1)**2
      
    if variant == 7 :
      val = 0
      for x in r:
        if x == 0 :
          val = float('inf')
          break
      if val == 0 :
        lx = [log(x) for x in r]
        s = sum(lx)
        av = s/n
        val = sum([(x-av)**2 for x in lx]) + n * (ai-1)**2

    if variant == 8 :
      for x in r:
        if x == 0 :
          val =  float('inf')
          break
      if val == 0 :
        lx = [log(x) for x in r]
        s = sum(lx)
        av = s/n
        val = sum([(x-av)**2 for x in lx]) + n**2 * (ai-1)**2
        
    if len(bsr[0]) == 2 :
      return val

    bs = bsnr + bsr

    if variant == 0 :
      cof = 2/(si*si)
      a1 = cof * n
      b1 = cof * si * (1 + val)
      dfdr = [a1 * r[k] - b1 for k in range(n)]

    if variant == 1 :
      m = si * ai - sum([x*x for x in r])
      ai2 = ai*ai
      sumu = (2*m)/(n * ai2*ai)

      p = 1/(2*sum([(x/ai - 1)**2 for x in r])**0.5)
      a1 = 2 * p/ai2
      c1 = (2/n) * (ai - 1)
      b1 = p * sumu - a1 * ai + c1
      dfdr = [a1 * x + b1 for x in r]

      #u = [2*(x/ai - 1) * (-x/ai**2) * (1/n) for x in r]
      #sumu = sum(u)
      #p = 1/(2*sum([(x/mean(r) - 1)**2 for x in r])**0.5)
      #dfdr = [p * (sumu - uk + 2 * (x/ai - 1) * (1/ai - (x/ai**2) * (1/n))) +
      #        2 * (ai - 1) * (1/n) for x,uk in zip(r,u)]

    if variant == 2 :
      m = si * ai - sum([x*x for x in r])
      ai2 = (ai*ai)
      sumu = (2*m)/(n * ai2*ai)

      a1 = 2/(ai2)
      c1 = (2/n) * (ai - 1) 
      b1 = sumu - a1 * ai + c1
      dfdr = [a1 * x  + b1 for x in r]

      ## u = [2*(x/ai - 1) * (-x/ai**2) * (1/n) for x in r]
      ## su = sum(u)
      ## dfdr = [su - uk + 2 * (x/ai - 1) * (1/ai - (x/ai**2) * (1/n)) + 2 * (ai - 1) * (1/n)
      ##         for x,uk in zip(r,u)]

    if variant == 3 :
      m = si * ai - sum([x*x for x in r])
      ai2 = (ai*ai)
      sumu = (2*m)/(n * ai2*ai)

      ep = 1e-6
      c1 = 0 if 1-ep < ai < 1+ep else ((1/n) if ai > 1 else (-1/n))
      
      a1 = 2/(ai2)
      b1 = sumu - a1 * ai + c1
      dfdr = [a1 * x  + b1 for x in r]
      if 0 :
        u = [2*(x/ai - 1) * (-x/ai**2) * (1/n) for x in r]
        su = sum(u)
        ep = 1e-6
        c1 = 0 if 1-ep < ai < 1+ep else ((1/n) if ai > 1 else (-1/n))
        dfdr = [su - uk + 2 * (x/ai - 1) * (1/ai - (x/ai**2) * (1/n)) + c1
                for x,uk in zip(r,u)]

    if variant == 4 :
      m = si * ai - sum([x*x for x in r])
      ai2 = ai*ai
      sumu = (2*m)/(n * ai2*ai)

      p = 1/(2*sum([(x/ai - 1)**2 for x in r])**0.5)
      a1 = 2 * p/ai2

      ep = 1e-6
      c1 = 0 if 1-ep < ai < 1+ep else ((1/n) if ai > 1 else (-1/n))

      b1 = p * sumu - a1 * ai + c1
      dfdr = [a1 * x + b1 for x in r]

    if variant == 5 :
      m = si * ai - sum([x*x for x in r])
      ai2 = (ai*ai)
      sumu = (2*m)/(n * ai2*ai)

      a1 = 2/(ai2)
      c1 = (2) * (ai - 1) 
      b1 = sumu - a1 * ai + c1
      dfdr = [a1 * x  + b1 for x in r]

    if variant == 6 :
      m = si * ai - sum([x*x for x in r])
      ai2 = (ai*ai)
      sumu = (2*m)/(n * ai2*ai)

      a1 = 2/(ai2)
      c1 = (2*n) * (ai - 1) 
      b1 = sumu - a1 * ai + c1
      dfdr = [a1 * x  + b1 for x in r]

    if variant == 7 :
      lx = [log(x) for x in r]
      av = sum(lx)/n
      c1 = 2 * (ai - 1)
      a1 = (2 * (n-2))/n
      dfdr = [((l - av) * a1)/x + c1  for l,x in zip(lx,r)]
      
    if variant == 8 :
      lx = [log(x) for x in r]
      av = sum(lx)/n
      c1 = 2 * n * (ai - 1)
      a1 = (2 * (n-2))/n
      dfdr = [((l - av) * a1)/x + c1  for l,x in zip(lx,r)]   
    
    dfdir = [dfdr[k] * -(r[k]/bs[k][0]) for k in range(n)]

    a = [bs[k][2] for k in range(n)]
    m = len(a[0])
    p = [0]*m
    for j in range(m) :
      p[j] = sum([dfdir[k] * a[k][j] for k in range(n)])

    pa = dfdr[-2]/bsr[0][0] - (dfdr[-1]/bsr[1][0])
    p.append(pa)

    return val,p
  
else :
  #from biopy.cchelp import normedVarianceAndDeriv as _normedVarianceAndDeriv
  from biopy.cchelp import varianceAndDerive as _varianceAndDerive
  from biopy.cchelp import normedVarianceAndDerivNew as _normedVarianceAndDeriv
  
def _UMTree_bfgs(tin, targetHeight = None, variant = 0, withDerivative = False,
                 doInit = True, normed = True, lbfgsb = True, limit = scipy.inf,
                 factr=1e7, warnings = False, internals = False, iprint=-1) :
  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # the copy.
  tree = copy.deepcopy(tin)
  targetHeight = calcNoderateMultipliers(tree, targetHeight)

  allids = tree.all_ids()
  allnodes = [(ni,tree.node(ni)) for ni in allids if ni != tree.root]

  nonrootids = []
  rootids = []
  for ni,n in allnodes :
    n = tree.node(ni)
    n.data.subsbranchlength = n.data.branchlength
    n.data.branchlength = n.data.timebranchlength
    del n.data.timebranchlength
    
    if n.prev == tree.root :
      rootids.append((ni,n))
    else :
      nonrootids.append((ni,n))
      
  nr = 2*len(tree.get_terminals())-2

  pee = ""

  bsubsroot = sum([n.data.subsbranchlength for ni,n in rootids])

  assert  normed or variant == 0
  
  meth = "_normedVarianceAndDeriv" if normed else "_varianceAndDerive"
  if withDerivative :
    pee += "  pnr = (" + ",".join(["(b%d,%.15g,d_b%d_x)" % (ni,n.data.subsbranchlength,ni)
                        for (ni,n) in nonrootids]) + ")\n"
    pee += "  pr = (" + ",".join(["(b%d,%.15g,d_b%d_x)" % (ni,n.data.subsbranchlength,ni)
                        for (ni,n) in rootids]) + ")\n"
    pee += "  val,valp = %s(pnr, pr, x[-1], %d)\n" % (meth, variant)
           
  else :
    pee += "  pnr = (" + ",".join(["(b%d,%.15g)" % (ni,n.data.subsbranchlength)
                        for (ni,n) in nonrootids]) + ")\n"
    pee += "  pr = (" + ",".join(["(b%d,%.15g)" % (ni,n.data.subsbranchlength)
                        for (ni,n) in rootids]) + ")\n"
    pee += "  val = %s(pnr, pr, x[-1], %d)\n" % (meth, variant)
      
  ba,htox = _treeBranchAssignmentExprs(tree, targetHeight, withDerivative = withDerivative,
                                       withInit = True)

  # Define the distance function on the fly.
  cod = ("def v1score(x):\n  " + "\n  ".join(ba) + "\n" + pee + "\n  return " +
         (("( val , array(valp) )") if withDerivative else " val "))
  exec cod in globals()

  xcod = "def htoxs():\n  x = [0]*%d\n  " % (nr//2 -1) + "\n  ".join(htox) \
         + "\n  return x + [%.14f]" % (bsubsroot/2)
  exec xcod in globals()

  # Function to obtain branch lengths from encoding
  codb = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['(%d,b%d)' % (k,k)
                              for k in allids if k != tree.root]) + ")" 
  exec codb in globals()

  x0 = htoxs()
  if not doInit :
    x0 = [.5] * (len(x0)-1) + x0[-1:]

  if lbfgsb :
    sol = scipy.optimize.fmin_l_bfgs_b(v1score, x0,
                                       approx_grad= 0 if withDerivative else 1,
                                       bounds = [[1e-3,1-1e-3]]*(len(x0)-1) + [[0,bsubsroot]],
                                       iprint=iprint, factr = factr)
    if warnings and sol[2]['warnflag'] != 0 :
      print "WARNING:", sol[2]['task']
  else :
    sol = scipy.optimize.fmin_tnc(v1score, x0, approx_grad=1,
                                  bounds = [[1e-5,1-1e-5]]*(len(x0)),
                                  messages=0)
    assert sol[2] == 1

  finaleVal = v1score(sol[0])
  if withDerivative :
    finaleVal = finaleVal[0]

  brs = code2branches(sol[0])
  for nn,br in brs:
    # numerical instability : don't permit negative branches
    n = tree.node(nn)
    n.data.branchlength = max(br, 0)
    if n.data.branchlength > 0 :
      if n.prev == tree.root :
        if n.id == rootids[0][0] :
          n.data.clockrate = sol[0][-1]/n.data.branchlength
          #print "1", sol[0][-1], n.data.clockrate
        else :
          n.data.clockrate = (bsubsroot - sol[0][-1])/n.data.branchlength
          #print "2", (bsubsroot - sol[0][-1]),n.data.clockrate
      else :
        n.data.clockrate = n.data.subsbranchlength/n.data.branchlength
      del n.data.subsbranchlength

  if internals :
    return (tree,finaleVal , (v1score,htoxs,cod,xcod,codb) )
  
  return tree,finaleVal

def UMTree(tin, targetHeight = None, method = "bgfs", normed = True,
           withDerivative = True, variant = 0) :
  if method == "slsqp" :
    tree,val = _UMTree_slsqp(tin, targetHeight = None, normed = normed,
                            useDerives = withDerivative, doInit = True)
  elif method == "bgfs" :
    assert 0 <= variant <= 8

    tree, val = _UMTree_bfgs(tin, targetHeight = None, normed = normed,
                             variant = variant,
                             withDerivative = withDerivative, doInit = True)
  else :
    raise
           
  return tree,val
    
