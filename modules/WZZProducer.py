import ROOT
from ROOT import TLorentzVector
from itertools import combinations, permutations
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os,copy
import numpy as np
from numpy import sign
from numpy import argsort

MW, MZ = 80.4, 91.2

class Particle:
  def __init__(self, p4, Id, charge):
    self.p4 = p4         
    self.Id = Id
    self.charge = charge

def METXYCorr_Met_MetPhi(uncormet,uncormet_phi,year,runera,isMC,npv):
  if npv>100: npv=100
  if isMC:
    if '2018' in year:
      METxcorr = -(0.183518*npv +0.546754)
      METycorr = -(0.192263*npv +-0.42121)
    elif '2017' in year:
      METxcorr = -(-0.300155*npv +1.90608)
      METycorr = -(0.300213*npv +-2.02232)
    elif '2016apv' in year:
      METxcorr = -(-0.188743*npv +0.136539)
      METycorr = -(0.0127927*npv +0.117747)
    else:
      METxcorr = -(-0.153497*npv +-0.231751)
      METycorr = -(0.00731978*npv +0.243323)
  else:
    if '2018' in year:
      if 'A' in runera:
        METxcorr = -(0.263733*npv +-1.91115)
        METycorr = -(0.0431304*npv +-0.112043)
      elif 'B' in runera:
        METxcorr = -(0.400466*npv +-3.05914)
        METycorr = -(0.146125*npv +-0.533233)
      elif 'C' in runera:
        METxcorr = -(0.430911*npv +-1.42865);
        METycorr = -(0.0620083*npv +-1.46021);
      elif 'D' in runera:
        METxcorr = -(0.457327*npv +-1.56856);
        METycorr = -(0.0684071*npv +-0.928372);
    elif '2017' in year:
      if 'B' in runera:
        METxcorr = -(-0.211161*npv +0.419333);
        METycorr = -(0.251789*npv +-1.28089);
      elif 'C' in runera:
        METxcorr = -(-0.185184*npv +-0.164009);
        METycorr = -(0.200941*npv +-0.56853);
      elif 'D' in runera:
        METxcorr = -(-0.201606*npv +0.426502);
        METycorr = -(0.188208*npv +-0.58313);
      elif 'E' in runera:
        METxcorr = -(-0.162472*npv +0.176329);
        METycorr = -(0.138076*npv +-0.250239);
      elif 'F' in runera:
        METxcorr = -(-0.210639*npv +0.72934);
        METycorr = -(0.198626*npv +1.028);
    elif '2016apv' in year:
      if 'B' in runera:
        METxcorr = -(-0.0214894*npv +-0.188255);
        METycorr = -(0.0876624*npv +0.812885);
      elif 'C' in runera:
        METxcorr = -(-0.032209*npv +0.067288);
        METycorr = -(0.113917*npv +0.743906);
      elif 'D' in runera:
        METxcorr = -(-0.0293663*npv +0.21106);
        METycorr = -(0.11331*npv +0.815787);
      elif 'E' in runera:
        METxcorr = -(-0.0132046*npv +0.20073);
        METycorr = -(0.134809*npv +0.679068);
      elif 'F' in runera:
        METxcorr = -(-0.0543566*npv +0.816597);
        METycorr = -(0.114225*npv +1.17266);
    else:
      if 'F' in runera:
        METxcorr = -(0.134616*npv +-0.89965);
        METycorr = -(0.0397736*npv +1.0385);
      elif 'G' in runera:
        METxcorr = -(0.121809*npv +-0.584893);
        METycorr = -(0.0558974*npv +0.891234);
      elif 'H' in runera:
        METxcorr = -(0.0868828*npv +-0.703489);
        METycorr = -(0.0888774*npv +0.902632);
    
  CorrectedMET_x = uncormet *math.cos( uncormet_phi)+METxcorr;
  CorrectedMET_y = uncormet *math.sin( uncormet_phi)+METycorr;
  CorrectedMET = math.sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y)

  if CorrectedMET_x==0 and CorrectedMET_y>0: CorrectedMETPhi = math.pi
  elif CorrectedMET_x==0 and CorrectedMET_y<0: CorrectedMETPhi = -1*math.pi
  elif CorrectedMET_x>0: CorrectedMETPhi = math.atan(CorrectedMET_y/CorrectedMET_x)
  elif CorrectedMET_x<0 and CorrectedMET_y>0: CorrectedMETPhi = math.atan(CorrectedMET_y/CorrectedMET_x) + math.pi
  elif CorrectedMET_x<0 and CorrectedMET_y<0: CorrectedMETPhi = math.atan(CorrectedMET_y/CorrectedMET_x) - math.pi
  else: CorrectedMETPhi =0;
  
  return (CorrectedMET, CorrectedMETPhi)

def w_v4(lep_v4, MET, MET_phi):
  wv4=TLorentzVector()
  vv4=TLorentzVector()
  pxl=lep_v4.Px()
  pyl=lep_v4.Py()
  pzl=lep_v4.Pz()
  El=lep_v4.E()
  pxv=MET*math.cos(MET_phi)
  pyv=MET*math.sin(MET_phi)
  Ev=MET

  a=MW*MW + 2*pxl*pxv + 2*pyl*pyv
  A=4*(El*El - pzl*pzl)
  B=-4*a*pzl
  C=4*El*El*MET*MET - a*a
  tmproot = B * B - 4.0 * A * C
  if tmproot<0:
    pzv = -B / (2 * A)
  else:
    tmpsol1 = (-B + math.sqrt(tmproot)) / (2.0 * A)
    tmpsol2 = (-B - math.sqrt(tmproot)) / (2.0 * A)
    if abs(tmpsol1) < abs(tmpsol2):
      pzv = tmpsol1
    else:
      pzv = tmpsol2

  vv4.SetPxPyPzE(pxv,pyv,pzv,math.sqrt(pxv*pxv+pyv*pyv+pzv*pzv))
  wv4=vv4+lep_v4
  return wv4

def assign_jets_from_many(jets):
  assert len(jets) >= 4
  best_score = 10000.
  best_combo = None
  for jet_indices in combinations(range(len(jets)), 4):
    selected_jets = [jets[i] for i in jet_indices]
    pairings = [
      ((0, 1), (2, 3)),
      ((0, 2), (1, 3)),
      ((0, 3), (1, 2))
    ]
    for (i1, i2), (j1, j2) in pairings:
      m1 = (selected_jets[i1]+selected_jets[i2]).M()
      m2 = (selected_jets[j1]+selected_jets[j2]).M()
      comb1=(m1 - MW)**2 + (m2 - MZ)**2
      comb2=(m1 - MZ)**2 + (m2 - MW)**2
      score = min( comb1, comb2)
      if score < best_score:
        best_score = score
        if comb1<comb2:
          best_combo = {
              'W_pair': (jet_indices[i1], jet_indices[i2]),
              'Z_pair': (jet_indices[j1], jet_indices[j2])
          }
        else:
          best_combo = {
              'W_pair': (jet_indices[j1], jet_indices[j2]),
              'Z_pair': (jet_indices[i1], jet_indices[i2])
          }
  return best_combo

def assign_2jets_from_many(jets,mass):
  assert len(jets) >= 2
  mass_boson=mass
  best_score = 10000.
  best_combo = None
  for jet_indices in combinations(range(len(jets)), 2):
    selected_jets = [jets[i] for i in jet_indices]
    mass_tmp = abs((selected_jets[0]+selected_jets[1]).M()-mass_boson)
    if mass_tmp<best_score:
      best_score=mass_tmp
      best_combo=jet_indices
  return best_combo


def assign_3lep(l1v4,l2v4,l3v4,l1charge,l2charge,l3charge,l1id,l2id,l3id):
  # return order, zl1_v4, zl2_v4, wl_v4
  assert abs(l1charge+l2charge+l3charge)==1
  if (l1charge+l2charge+l3charge)==1:
    if l1charge==-1:
      if abs((l1v4+l2v4).M()-MZ)<abs((l1v4+l3v4).M()-MZ):
        return (l1v4,l2v4,l3v4,l1id,l2id,l3id)
      else:
        return (l1v4,l3v4,l2v4,l1id,l3id,l2id)
    elif l2charge==-1:
      if abs((l2v4+l1v4).M()-MZ)<abs((l2v4+l3v4).M()-MZ):
        return (l1v4,l2v4,l3v4,l1id,l2id,l3id)
      else:
        return (l2v4,l3v4,l1v4,l2id,l3id,l1id)
    else:
      if abs((l3v4+l1v4).M()-MZ)<abs((l3v4+l2v4).M()-MZ):
        return (l1v4,l3v4,l2v4,l1id,l3id,l2id)
      else:
        return (l2v4,l3v4,l1v4,l2id,l3id,l1id)
  else:
    if l1charge==1:
      if abs((l1v4+l2v4).M()-MZ)<abs((l1v4+l3v4).M()-MZ):
        return (l1v4,l2v4,l3v4,l1id,l2id,l3id)
      else:
        return (l1v4,l3v4,l2v4,l1id,l3id,l2id)
    elif l2charge==1:
      if abs((l2v4+l1v4).M()-MZ)<abs((l2v4+l3v4).M()-MZ):
        return (l1v4,l2v4,l3v4,l1id,l2id,l3id)
      else:
        return (l2v4,l3v4,l1v4,l2id,l3id,l1id)
    else:
      if abs((l3v4+l1v4).M()-MZ)<abs((l3v4+l2v4).M()-MZ):
        return (l1v4,l3v4,l2v4,l1id,l3id,l2id)
      else:
        return (l2v4,l3v4,l1v4,l2id,l3id,l1id)

def assign_4lep(l1v4,l2v4,l3v4,l4v4,l1charge,l2charge,l3charge,l4charge,l1id,l2id,l3id,l4id):
  if l1charge+l2charge==0:
    if l1charge+l3charge==0:
      comb=[abs((l1v4+l2v4).M()-MZ), abs((l3v4+l4v4).M()-MZ), abs((l1v4+l3v4).M()-MZ), abs((l2v4+l4v4).M()-MZ)]
      if min(comb) == abs((l1v4+l2v4).M()-MZ):
        return(l1v4,l2v4,l3v4,l4v4,l1id,l2id,l3id,l4id)
      elif min(comb) == abs((l3v4+l4v4).M()-MZ):
        return(l3v4,l4v4,l1v4,l2v4,l3id,l4id,l1id,l2id)
      elif min(comb) == abs((l1v4+l3v4).M()-MZ):
        return(l1v4,l3v4,l2v4,l4v4,l1id,l3id,l2id,l4id)
      else:
        return(l2v4,l4v4,l1v4,l3v4,l2id,l4id,l1id,l3id)
    else:
      comb=[abs((l1v4+l2v4).M()-MZ), abs((l3v4+l4v4).M()-MZ), abs((l1v4+l4v4).M()-MZ), abs((l2v4+l3v4).M()-MZ)]
      if min(comb) == abs((l1v4+l2v4).M()-MZ):
        return(l1v4,l2v4,l3v4,l4v4,l1id,l2id,l3id,l4id)
      elif min(comb) == abs((l3v4+l4v4).M()-MZ):
        return(l3v4,l4v4,l1v4,l2v4,l3id,l4id,l1id,l2id)
      elif min(comb) == abs((l1v4+l4v4).M()-MZ):
        return(l1v4,l4v4,l2v4,l3v4,l1id,l4id,l2id,l3id)
      else:
        return(l2v4,l3v4,l1v4,l4v4,l2id,l3id,l1id,l4id)
  else:
    comb=[abs((l1v4+l3v4).M()-MZ), abs((l2v4+l4v4).M()-MZ), abs((l1v4+l4v4).M()-MZ), abs((l2v4+l3v4).M()-MZ)]
    if min(comb) == abs((l1v4+l3v4).M()-MZ):
      return(l1v4,l3v4,l2v4,l4v4,l1id,l3id,l2id,l4id)
    elif min(comb) == abs((l2v4+l4v4).M()-MZ):
      return(l2v4,l4v4,l1v4,l3v4,l2id,l4id,l1id,l3id)
    elif min(comb) == abs((l1v4+l4v4).M()-MZ):
      return(l1v4,l4v4,l2v4,l3v4,l1id,l4id,l2id,l3id)
    else:
      return(l2v4,l3v4,l1v4,l4v4,l2id,l3id,l1id,l4id)
     
def assign_5lep(l1,l2,l3,l4,l5):
  leptons=[l1,l2,l3,l4,l5]
  best_cost = 10000
  best_combo = None

  z_pairs = []
  for i, j in combinations(range(5), 2):
    li, lj = leptons[i], leptons[j]
    if li.charge + lj.charge != 0: continue
    z_pairs.append(((i, j), (li.p4 + lj.p4).M()))

  for (i1, j1), m1 in z_pairs:
    for (i2, j2), m2 in z_pairs:
      indices = {i1, j1, i2, j2}
      if len(indices) != 4: continue

      cost = (m1 - MZ)**2 + (m2 - MZ)**2
      if cost < best_cost:
        best_cost = cost
        used_idx = indices
        rest_idx = list(set(range(5)) - used_idx)
        if abs(m1 - MZ) < abs(m2 - MZ):
          Z1_pair = (leptons[i1], leptons[j1])
          Z2_pair = (leptons[i2], leptons[j2])
        else:
          Z1_pair = (leptons[i2], leptons[j2])
          Z2_pair = (leptons[i1], leptons[j1])
        best_combo = {
            'Z1': Z1_pair,
            'Z2': Z2_pair,
            'W': leptons[rest_idx[0]]
        }

  return best_combo

class WZZProducer(Module):
  def __init__(self , year, runera):
    self.year = year
    self.runera = runera
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("lhe_nlepton", "I")
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("met_user","F")
    self.out.branch("met_phi_user","F")
    self.out.branch("met_CorrePhi","F")
    self.out.branch("met_phi_CorrePhi","F")
    self.out.branch("TMuons_id","I",lenVar="nTMuons")
    self.out.branch("FMuons_id","I",lenVar="nFMuons")
    self.out.branch("LMuons_id","I",lenVar="nLMuons")
    self.out.branch("TElectrons_id","I",lenVar="nTElectrons")
    self.out.branch("FElectrons_id","I",lenVar="nFElectrons")
    self.out.branch("LElectrons_id","I",lenVar="nLElectrons")
    self.out.branch("TightJet_id","I",lenVar="nTightJet")
    self.out.branch("TightJet_pt","F",lenVar="nTightJet")
    self.out.branch("TightJet_eta","F",lenVar="nTightJet")
    self.out.branch("TightJet_phi","F",lenVar="nTightJet")
    self.out.branch("TightJet_mass","F",lenVar="nTightJet")

    #2L
    self.out.branch("SR_2L_FlavorBit","I")
    self.out.branch("SR_2L_FakeBit","I")
    self.out.branch("SR2L_l1_id","I")
    self.out.branch("SR2L_l2_id","I")
    self.out.branch("SR2L_l1_pt","F")
    self.out.branch("SR2L_l1_eta","F")
    self.out.branch("SR2L_l1_phi","F")
    self.out.branch("SR2L_l1_mass","F")
    self.out.branch("SR2L_l2_pt","F")
    self.out.branch("SR2L_l2_eta","F")
    self.out.branch("SR2L_l2_phi","F")
    self.out.branch("SR2L_l2_mass","F")
    self.out.branch("SR2L_zll_pt","F")
    self.out.branch("SR2L_zll_eta","F")
    self.out.branch("SR2L_zll_phi","F")
    self.out.branch("SR2L_zll_mass","F")
    self.out.branch("SR2L_zj1_id","I")
    self.out.branch("SR2L_zj1_pt","F")
    self.out.branch("SR2L_zj1_eta","F")
    self.out.branch("SR2L_zj1_phi","F")
    self.out.branch("SR2L_zj1_mass","F")
    self.out.branch("SR2L_zj2_id","I")
    self.out.branch("SR2L_zj2_pt","F")
    self.out.branch("SR2L_zj2_eta","F")
    self.out.branch("SR2L_zj2_phi","F")
    self.out.branch("SR2L_zj2_mass","F")
    self.out.branch("SR2L_wj1_id","I")
    self.out.branch("SR2L_wj1_pt","F")
    self.out.branch("SR2L_wj1_eta","F")
    self.out.branch("SR2L_wj1_phi","F")
    self.out.branch("SR2L_wj1_mass","F")
    self.out.branch("SR2L_wj2_id","I")
    self.out.branch("SR2L_wj2_pt","F")
    self.out.branch("SR2L_wj2_eta","F")
    self.out.branch("SR2L_wj2_phi","F")
    self.out.branch("SR2L_wj2_mass","F")
    self.out.branch("SR2L_zjj_pt","F")
    self.out.branch("SR2L_zjj_eta","F")
    self.out.branch("SR2L_zjj_phi","F")
    self.out.branch("SR2L_zjj_mass","F")
    self.out.branch("SR2L_wjj_pt","F")
    self.out.branch("SR2L_wjj_eta","F")
    self.out.branch("SR2L_wjj_phi","F")
    self.out.branch("SR2L_wjj_mass","F")
    self.out.branch("SR2L_dR_zll_wjj","F")
    self.out.branch("SR2L_dEta_zll_wjj","F")
    self.out.branch("SR2L_dPhi_zll_wjj","F")
    self.out.branch("SR2L_zllwjj_pt","F")
    self.out.branch("SR2L_zllwjj_eta","F")
    self.out.branch("SR2L_zllwjj_phi","F")
    self.out.branch("SR2L_zllwjj_mass","F")
    self.out.branch("SR2L_dR_zll_zjj","F")
    self.out.branch("SR2L_dEta_zll_zjj","F")
    self.out.branch("SR2L_dPhi_zll_zjj","F")
    self.out.branch("SR2L_zllzjj_pt","F")
    self.out.branch("SR2L_zllzjj_eta","F")
    self.out.branch("SR2L_zllzjj_phi","F")
    self.out.branch("SR2L_zllzjj_mass","F")
    self.out.branch("SR2L_dR_wjj_zjj","F")
    self.out.branch("SR2L_dEta_wjj_zjj","F")
    self.out.branch("SR2L_dPhi_wjj_zjj","F")
    self.out.branch("SR2L_wjjzjj_pt","F")
    self.out.branch("SR2L_wjjzjj_eta","F")
    self.out.branch("SR2L_wjjzjj_phi","F")
    self.out.branch("SR2L_wjjzjj_mass","F")
    self.out.branch("SR2L_wzz_pt","F")
    self.out.branch("SR2L_wzz_eta","F")
    self.out.branch("SR2L_wzz_phi","F")
    self.out.branch("SR2L_wzz_mass","F")
    #3L
    self.out.branch("SR_3L_FlavorBit","I")
    self.out.branch("SR_3L_FakeBit","I")
    self.out.branch("SR3L_zl1_id","I")
    self.out.branch("SR3L_zl2_id","I")
    self.out.branch("SR3L_wl_id","I")
    self.out.branch("SR3L_zl1_pt","F")
    self.out.branch("SR3L_zl1_eta","F")
    self.out.branch("SR3L_zl1_phi","F")
    self.out.branch("SR3L_zl1_mass","F")
    self.out.branch("SR3L_zl2_pt","F")
    self.out.branch("SR3L_zl2_eta","F")
    self.out.branch("SR3L_zl2_phi","F")
    self.out.branch("SR3L_zl2_mass","F")
    self.out.branch("SR3L_wl_pt","F")
    self.out.branch("SR3L_wl_eta","F")
    self.out.branch("SR3L_wl_phi","F")
    self.out.branch("SR3L_wl_mass","F")
    self.out.branch("SR3L_zll_pt","F")
    self.out.branch("SR3L_zll_eta","F")
    self.out.branch("SR3L_zll_phi","F")
    self.out.branch("SR3L_zll_mass","F")
    self.out.branch("SR3L_wlv_pt","F")
    self.out.branch("SR3L_wlv_eta","F")
    self.out.branch("SR3L_wlv_phi","F")
    self.out.branch("SR3L_wlv_mass","F")
    self.out.branch("SR3L_zj1_id","I")
    self.out.branch("SR3L_zj1_pt","F")
    self.out.branch("SR3L_zj1_eta","F")
    self.out.branch("SR3L_zj1_phi","F")
    self.out.branch("SR3L_zj1_mass","F")
    self.out.branch("SR3L_zj2_id","I")
    self.out.branch("SR3L_zj2_pt","F")
    self.out.branch("SR3L_zj2_eta","F")
    self.out.branch("SR3L_zj2_phi","F")
    self.out.branch("SR3L_zj2_mass","F")
    self.out.branch("SR3L_zjj_pt","F")
    self.out.branch("SR3L_zjj_eta","F")
    self.out.branch("SR3L_zjj_phi","F")
    self.out.branch("SR3L_zjj_mass","F")
    self.out.branch("SR3L_dR_zll_wlv","F")
    self.out.branch("SR3L_dEta_zll_wlv","F")
    self.out.branch("SR3L_dPhi_zll_wlv","F")
    self.out.branch("SR3L_zllwlv_pt","F")
    self.out.branch("SR3L_zllwlv_eta","F")
    self.out.branch("SR3L_zllwlv_phi","F")
    self.out.branch("SR3L_zllwlv_mass","F")
    self.out.branch("SR3L_dR_zll_zjj","F")
    self.out.branch("SR3L_dEta_zll_zjj","F")
    self.out.branch("SR3L_dPhi_zll_zjj","F")
    self.out.branch("SR3L_zllzjj_pt","F")
    self.out.branch("SR3L_zllzjj_eta","F")
    self.out.branch("SR3L_zllzjj_phi","F")
    self.out.branch("SR3L_zllzjj_mass","F")
    self.out.branch("SR3L_dR_wlv_zjj","F")
    self.out.branch("SR3L_dEta_wlv_zjj","F")
    self.out.branch("SR3L_dPhi_wlv_zjj","F")
    self.out.branch("SR3L_wlvzjj_pt","F")
    self.out.branch("SR3L_wlvzjj_eta","F")
    self.out.branch("SR3L_wlvzjj_phi","F")
    self.out.branch("SR3L_wlvzjj_mass","F")
    self.out.branch("SR3L_wzz_pt","F")
    self.out.branch("SR3L_wzz_eta","F")
    self.out.branch("SR3L_wzz_phi","F")
    self.out.branch("SR3L_wzz_mass","F")

    #4L
    self.out.branch("SR_4L_FlavorBit","I")
    self.out.branch("SR_4L_FakeBit","I")
    self.out.branch("SR4L_z1l1_id","I")
    self.out.branch("SR4L_z1l2_id","I")
    self.out.branch("SR4L_z2l1_id","I")
    self.out.branch("SR4L_z2l2_id","I")
    self.out.branch("SR4L_z1l1_pt","F")
    self.out.branch("SR4L_z1l1_eta","F")
    self.out.branch("SR4L_z1l1_phi","F")
    self.out.branch("SR4L_z1l1_mass","F")
    self.out.branch("SR4L_z1l2_pt","F")
    self.out.branch("SR4L_z1l2_eta","F")
    self.out.branch("SR4L_z1l2_phi","F")
    self.out.branch("SR4L_z1l2_mass","F")
    self.out.branch("SR4L_z2l1_pt","F")
    self.out.branch("SR4L_z2l1_eta","F")
    self.out.branch("SR4L_z2l1_phi","F")
    self.out.branch("SR4L_z2l1_mass","F")
    self.out.branch("SR4L_z2l2_pt","F")
    self.out.branch("SR4L_z2l2_eta","F")
    self.out.branch("SR4L_z2l2_phi","F")
    self.out.branch("SR4L_z2l2_mass","F")
    self.out.branch("SR4L_z1_pt","F")
    self.out.branch("SR4L_z1_eta","F")
    self.out.branch("SR4L_z1_phi","F")
    self.out.branch("SR4L_z1_mass","F")
    self.out.branch("SR4L_z2_pt","F")
    self.out.branch("SR4L_z2_eta","F")
    self.out.branch("SR4L_z2_phi","F")
    self.out.branch("SR4L_z2_mass","F")
    self.out.branch("SR4L_wj1_id","I")
    self.out.branch("SR4L_wj1_pt","F")
    self.out.branch("SR4L_wj1_eta","F")
    self.out.branch("SR4L_wj1_phi","F")
    self.out.branch("SR4L_wj1_mass","F")
    self.out.branch("SR4L_wj2_id","I")
    self.out.branch("SR4L_wj2_pt","F")
    self.out.branch("SR4L_wj2_eta","F")
    self.out.branch("SR4L_wj2_phi","F")
    self.out.branch("SR4L_wj2_mass","F")
    self.out.branch("SR4L_wjj_pt","F")
    self.out.branch("SR4L_wjj_eta","F")
    self.out.branch("SR4L_wjj_phi","F")
    self.out.branch("SR4L_wjj_mass","F")
    self.out.branch("SR4L_dR_z1_wjj","F")
    self.out.branch("SR4L_dEta_z1_wjj","F")
    self.out.branch("SR4L_dPhi_z1_wjj","F")
    self.out.branch("SR4L_z1wjj_pt","F")
    self.out.branch("SR4L_z1wjj_eta","F")
    self.out.branch("SR4L_z1wjj_phi","F")
    self.out.branch("SR4L_z1wjj_mass","F")
    self.out.branch("SR4L_dR_z1_z2","F")
    self.out.branch("SR4L_dEta_z1_z2","F")
    self.out.branch("SR4L_dPhi_z1_z2","F")
    self.out.branch("SR4L_z1z2_pt","F")
    self.out.branch("SR4L_z1z2_eta","F")
    self.out.branch("SR4L_z1z2_phi","F")
    self.out.branch("SR4L_z1z2_mass","F")
    self.out.branch("SR4L_dR_wjj_z2","F")
    self.out.branch("SR4L_dEta_wjj_z2","F")
    self.out.branch("SR4L_dPhi_wjj_z2","F")
    self.out.branch("SR4L_wjjz2_pt","F")
    self.out.branch("SR4L_wjjz2_eta","F")
    self.out.branch("SR4L_wjjz2_phi","F")
    self.out.branch("SR4L_wjjz2_mass","F")
    self.out.branch("SR4L_wzz_pt","F")
    self.out.branch("SR4L_wzz_eta","F")
    self.out.branch("SR4L_wzz_phi","F")
    self.out.branch("SR4L_wzz_mass","F")
 
    #5L
    self.out.branch("SR_5L_FlavorBit","I")
    self.out.branch("SR_5L_FakeBit","I")
    self.out.branch("SR5L_z1l1_id","I")
    self.out.branch("SR5L_z1l2_id","I")
    self.out.branch("SR5L_z2l1_id","I")
    self.out.branch("SR5L_z2l2_id","I")
    self.out.branch("SR5L_wl_id","I")
    self.out.branch("SR5L_z1l1_pt","F")
    self.out.branch("SR5L_z1l1_eta","F")
    self.out.branch("SR5L_z1l1_phi","F")
    self.out.branch("SR5L_z1l1_mass","F")
    self.out.branch("SR5L_z1l2_pt","F")
    self.out.branch("SR5L_z1l2_eta","F")
    self.out.branch("SR5L_z1l2_phi","F")
    self.out.branch("SR5L_z1l2_mass","F")
    self.out.branch("SR5L_z2l1_pt","F")
    self.out.branch("SR5L_z2l1_eta","F")
    self.out.branch("SR5L_z2l1_phi","F")
    self.out.branch("SR5L_z2l1_mass","F")
    self.out.branch("SR5L_z2l2_pt","F")
    self.out.branch("SR5L_z2l2_eta","F")
    self.out.branch("SR5L_z2l2_phi","F")
    self.out.branch("SR5L_z2l2_mass","F")
    self.out.branch("SR5L_wl_pt","F")
    self.out.branch("SR5L_wl_eta","F")
    self.out.branch("SR5L_wl_phi","F")
    self.out.branch("SR5L_wl_mass","F")
    self.out.branch("SR5L_z1_pt","F")
    self.out.branch("SR5L_z1_eta","F")
    self.out.branch("SR5L_z1_phi","F")
    self.out.branch("SR5L_z1_mass","F")
    self.out.branch("SR5L_z2_pt","F")
    self.out.branch("SR5L_z2_eta","F")
    self.out.branch("SR5L_z2_phi","F")
    self.out.branch("SR5L_z2_mass","F")
    self.out.branch("SR5L_w_pt","F")
    self.out.branch("SR5L_w_eta","F")
    self.out.branch("SR5L_w_phi","F")
    self.out.branch("SR5L_w_mass","F")
    self.out.branch("SR5L_dR_z1_w","F")
    self.out.branch("SR5L_dEta_z1_w","F")
    self.out.branch("SR5L_dPhi_z1_w","F")
    self.out.branch("SR5L_z1w_pt","F")
    self.out.branch("SR5L_z1w_eta","F")
    self.out.branch("SR5L_z1w_phi","F")
    self.out.branch("SR5L_z1w_mass","F")
    self.out.branch("SR5L_dR_z1_z2","F")
    self.out.branch("SR5L_dEta_z1_z2","F")
    self.out.branch("SR5L_dPhi_z1_z2","F")
    self.out.branch("SR5L_z1z2_pt","F")
    self.out.branch("SR5L_z1z2_eta","F")
    self.out.branch("SR5L_z1z2_phi","F")
    self.out.branch("SR5L_z1z2_mass","F")
    self.out.branch("SR5L_dR_w_z2","F")
    self.out.branch("SR5L_dEta_w_z2","F")
    self.out.branch("SR5L_dPhi_w_z2","F")
    self.out.branch("SR5L_wz2_pt","F")
    self.out.branch("SR5L_wz2_eta","F")
    self.out.branch("SR5L_wz2_phi","F")
    self.out.branch("SR5L_wz2_mass","F")
    self.out.branch("SR5L_wzz_pt","F")
    self.out.branch("SR5L_wzz_eta","F")
    self.out.branch("SR5L_wzz_phi","F")
    self.out.branch("SR5L_wzz_mass","F")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    if (event.PV_npvsGood<1): return False

    met_user=-99
    met_phi_user=-99
    if self.is_mc:
      met_user=event.MET_T1Smear_pt
      met_phi_user=event.MET_T1Smear_phi
    else:
      met_user=event.MET_T1_pt
      met_phi_user=event.MET_T1_phi

    self.out.fillBranch("met_user",met_user)
    self.out.fillBranch("met_phi_user",met_phi_user)

    met_CorrePhi,met_phi_CorrePhi=METXYCorr_Met_MetPhi(event.RawMET_pt, event.RawMET_phi, self.year, self.runera, self.is_mc, event.PV_npvs)
    self.out.fillBranch("met_CorrePhi",met_CorrePhi)
    self.out.fillBranch("met_phi_CorrePhi",met_phi_CorrePhi)
    
    # recover HLT
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
        for iobj in range(0,event.nTrigObj):
          if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
            HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    lhe_nlepton=0
    if self.is_lhe:
      lheparticle = Collection(event, 'LHEPart')
      for ilhe in range(0, event.nLHEPart):
        if lheparticle[ilhe].status==1 and (abs(lheparticle[ilhe].pdgId)==11 or abs(lheparticle[ilhe].pdgId)==13 or abs(lheparticle[ilhe].pdgId)==15):
          lhe_nlepton=lhe_nlepton+1

    self.out.fillBranch("lhe_nlepton", lhe_nlepton)

    TMuons_id = []
    FMuons_id = []
    LMuons_id = []
    TElectrons_id = []
    FElectrons_id = []
    LElectrons_id = []

    muons = Collection(event, 'Muon')
    for imu in range(0, event.nMuon):
      if abs(muons[imu].dxy)>0.2 or abs(muons[imu].dz)>0.5:continue
      if abs(muons[imu].eta)<2.4 and muons[imu].pt>10 and muons[imu].tightId and muons[imu].pfRelIso04_all<0.15:
        TMuons_id.append(imu)
      if abs(muons[imu].eta)<2.4 and muons[imu].pt>10 and muons[imu].tightId and muons[imu].pfRelIso04_all>0.15 and muons[imu].pfRelIso04_all<0.4:
        FMuons_id.append(imu)
      if abs(muons[imu].eta)<2.4 and muons[imu].pt>10 and muons[imu].looseId and muons[imu].pfRelIso04_all<0.4:
        LMuons_id.append(imu)
    LMuons_id = list(set(LMuons_id) - set(TMuons_id) - set(FMuons_id))
    TFMuons_id = TMuons_id+FMuons_id

    eles = Collection(event, 'Electron')
    for iele in range(0, event.nElectron):
      if not eles[iele].convVeto:continue
      if abs(eles[iele].eta)>1.4442 and abs(eles[iele].eta)<1.566:continue
      if abs(eles[iele].eta)<1.479 and (abs(eles[iele].dxy)>0.05 or abs(eles[iele].dz)>0.1):continue
      if abs(eles[iele].eta)>1.479 and (abs(eles[iele].dxy)>0.1 or abs(eles[iele].dz)>0.2):continue
      if abs(eles[iele].eta)<2.5 and eles[iele].pt>10 and eles[iele].cutBased>2:
        TElectrons_id.append(iele)
      if abs(eles[iele].eta)<2.5 and eles[iele].pt>10 and eles[iele].cutBased==2:
        FElectrons_id.append(iele)
      if abs(eles[iele].eta)<2.5 and eles[iele].pt>10 and eles[iele].cutBased==1:
        LElectrons_id.append(iele)
    TFElectrons_id=TElectrons_id+FElectrons_id 

    self.out.fillBranch("TMuons_id", TMuons_id)
    self.out.fillBranch("FMuons_id", FMuons_id)
    self.out.fillBranch("LMuons_id", LMuons_id)
    self.out.fillBranch("TElectrons_id", TElectrons_id)
    self.out.fillBranch("FElectrons_id", FElectrons_id)
    self.out.fillBranch("LElectrons_id", LElectrons_id)

    jets = Collection(event, 'Jet')
    TightJet_id = []
    TightJet_pt = []
    TightJet_eta = []
    TightJet_phi = []
    TightJet_mass = []
    TightJet_v4 = []
    jet_v4_temp=TLorentzVector()
    lep_v4_temp=TLorentzVector()
    for ijet in range(0, event.nJet):
      if abs(jets[ijet].eta)>2.4 or jets[ijet].pt_nom<30: continue
      if jets[ijet].jetId<6:continue
      jet_v4_temp.SetPtEtaPhiM(jets[ijet].pt_nom,jets[ijet].eta,jets[ijet].phi,jets[ijet].mass_nom)
      pass_mu_dr=1
      pass_ele_dr=1
      pass_jet_dr=1

      for imu in range(0,len(TFMuons_id)):
        if pass_mu_dr<1:continue
        mid_tmp=TFMuons_id[imu]
        lep_v4_temp.SetPtEtaPhiM(muons[mid_tmp].pt, muons[mid_tmp].eta, muons[mid_tmp].phi, muons[mid_tmp].mass)
        if jet_v4_temp.DeltaR(lep_v4_temp)<0.4:pass_mu_dr=0
      
      for iele in range(0,len(TFElectrons_id)):
        if pass_ele_dr<1:continue
        eid_tmp=TFElectrons_id[iele]
        lep_v4_temp.SetPtEtaPhiM(eles[eid_tmp].pt, eles[eid_tmp].eta, eles[eid_tmp].phi, eles[eid_tmp].mass)
        if jet_v4_temp.DeltaR(lep_v4_temp)<0.4:pass_ele_dr=0

      if len(TightJet_id)>0:
        for ij in range(0,len(TightJet_id)):
          if jet_v4_temp.DeltaR(TightJet_v4[ij])<0.4:pass_jet_dr=0

      if pass_mu_dr>0 and pass_ele_dr>0 and pass_jet_dr>0:
        TightJet_id.append(ijet)
        TightJet_v4.append(jet_v4_temp)
        TightJet_pt.append(jet_v4_temp.Pt())
        TightJet_eta.append(jet_v4_temp.Eta())
        TightJet_phi.append(jet_v4_temp.Phi())
        TightJet_mass.append(jet_v4_temp.M())

    self.out.fillBranch("TightJet_id", TightJet_id)
    self.out.fillBranch("TightJet_pt", TightJet_pt)
    self.out.fillBranch("TightJet_eta", TightJet_eta)
    self.out.fillBranch("TightJet_phi", TightJet_phi)
    self.out.fillBranch("TightJet_mass", TightJet_mass)

    if len(LMuons_id)>0 or len(LElectrons_id)>0:return False
    if (len(TFMuons_id)+len(TFElectrons_id))<2:return False

    # Two leptons channel
    # SR_2L_FlavorBit=26: 2mu, 
    # SR_2L_FlavorBit=22: 2ele,

    # in this region, in 1 Real lepton case, the l2 is always the Fake one
    # SR_2L_FakeBit, bit1=0 if l1 is fake, bit1=1 if l1 is prompt, bit2=0 if l2 is fake, bit2=1 if l2 is prompt
    # SR_2L_FakeBit=3: l1 and l2 prompt
    # SR_2L_FakeBit=1: l1 prompt
    # SR_2L_FakeBit=0: l1 and l2 fake

    SR_2L_FlavorBit=-1
    SR_2L_FakeBit=-1
    SR2L_l1_id=-99
    SR2L_l2_id=-99
    SR2L_l1_pt=-99
    SR2L_l1_eta=-99
    SR2L_l1_phi=-99
    SR2L_l1_mass=-99
    SR2L_l2_pt=-99
    SR2L_l2_eta=-99
    SR2L_l2_phi=-99
    SR2L_l2_mass=-99
    SR2L_zll_pt=-99
    SR2L_zll_eta=-99
    SR2L_zll_phi=-99
    SR2L_zll_mass=-99
    SR2L_zj1_id=-99
    SR2L_zj1_pt=-99
    SR2L_zj1_eta=-99
    SR2L_zj1_phi=-99
    SR2L_zj1_mass=-99
    SR2L_zj2_id=-99
    SR2L_zj2_pt=-99
    SR2L_zj2_eta=-99
    SR2L_zj2_phi=-99
    SR2L_zj2_mass=-99
    SR2L_wj1_id=-99
    SR2L_wj1_pt=-99
    SR2L_wj1_eta=-99
    SR2L_wj1_phi=-99
    SR2L_wj1_mass=-99
    SR2L_wj2_id=-99
    SR2L_wj2_pt=-99
    SR2L_wj2_eta=-99
    SR2L_wj2_phi=-99
    SR2L_wj2_mass=-99
    SR2L_zjj_pt=-99
    SR2L_zjj_eta=-99
    SR2L_zjj_phi=-99
    SR2L_zjj_mass=-99
    SR2L_wjj_pt=-99
    SR2L_wjj_eta=-99
    SR2L_wjj_phi=-99
    SR2L_wjj_mass=-99
    SR2L_dR_zll_wjj = -99
    SR2L_dEta_zll_wjj =-99
    SR2L_dPhi_zll_wjj =-99
    SR2L_zllwjj_pt = -99
    SR2L_zllwjj_eta = -99
    SR2L_zllwjj_phi = -99
    SR2L_zllwjj_mass =-99
    SR2L_dR_zll_zjj = -99
    SR2L_dEta_zll_zjj =-99
    SR2L_dPhi_zll_zjj =-99
    SR2L_zllzjj_pt = -99
    SR2L_zllzjj_eta =-99
    SR2L_zllzjj_phi =-99
    SR2L_zllzjj_mass =-99
    SR2L_dR_wjj_zjj = -99
    SR2L_dEta_wjj_zjj =-99
    SR2L_dPhi_wjj_zjj =-99
    SR2L_wjjzjj_pt = -99
    SR2L_wjjzjj_eta =-99
    SR2L_wjjzjj_phi =-99
    SR2L_wjjzjj_mass =-99
    SR2L_wzz_pt = -99
    SR2L_wzz_eta = -99
    SR2L_wzz_phi = -99
    SR2L_wzz_mass =-99

    if (len(TFMuons_id)+len(TFElectrons_id))==2:
      if len(TFMuons_id)==1:return False
      if len(TightJet_id)<4:return False
      l1_v4=TLorentzVector()
      l2_v4=TLorentzVector()
      zll_v4=TLorentzVector()
      zjj_v4=TLorentzVector()
      wjj_v4=TLorentzVector()
      # w(qq)+z(ll)+z(qq) in muon channel
      if len(TFMuons_id)==2:
        SR_2L_FlavorBit=26
        if len(TMuons_id)==2:
          SR_2L_FakeBit=3
          SR2L_l1_id=TMuons_id[0]
          SR2L_l2_id=TMuons_id[1]
        elif len(TMuons_id)==1:
          SR_2L_FakeBit=1
          SR2L_l1_id=TMuons_id[0]
          SR2L_l2_id=FMuons_id[0]
        else:
          SR_2L_FakeBit=0
          SR2L_l1_id=FMuons_id[0]
          SR2L_l2_id=FMuons_id[1]
        SR2L_l1_pt=muons[SR2L_l1_id].pt
        SR2L_l1_eta=muons[SR2L_l1_id].eta
        SR2L_l1_phi=muons[SR2L_l1_id].phi
        SR2L_l1_mass=muons[SR2L_l1_id].mass
        SR2L_l2_pt=muons[SR2L_l2_id].pt
        SR2L_l2_eta=muons[SR2L_l2_id].eta
        SR2L_l2_phi=muons[SR2L_l2_id].phi
        SR2L_l2_mass=muons[SR2L_l2_id].mass
        l1_v4.SetPtEtaPhiM(SR2L_l1_pt,SR2L_l1_eta,SR2L_l1_phi,SR2L_l1_mass)
        l2_v4.SetPtEtaPhiM(SR2L_l2_pt,SR2L_l2_eta,SR2L_l2_phi,SR2L_l2_mass)
        SR2L_zll_pt=(l1_v4+l2_v4).Pt()
        SR2L_zll_eta=(l1_v4+l2_v4).Eta()
        SR2L_zll_phi=(l1_v4+l2_v4).Phi()
        SR2L_zll_mass=(l1_v4+l2_v4).M()
      # w(qq)+z(ll)+z(qq) in ele channel
      elif len(TFElectrons_id)==2:
        SR_2L_FlavorBit=22
        if len(TElectrons_id)==2:
          SR_2L_FakeBit=3
          SR2L_l1_id=TElectrons_id[0]
          SR2L_l2_id=TElectrons_id[1]
        elif len(TElectrons_id)==1:
          SR_2L_FakeBit=1
          SR2L_l1_id=TElectrons_id[0]
          SR2L_l2_id=FElectrons_id[0]
        else:
          SR_2L_FakeBit=0
          SR2L_l1_id=FElectrons_id[0]
          SR2L_l2_id=FElectrons_id[1]
        SR2L_l1_pt=eles[SR2L_l1_id].pt
        SR2L_l1_eta=eles[SR2L_l1_id].eta
        SR2L_l1_phi=eles[SR2L_l1_id].phi
        SR2L_l1_mass=eles[SR2L_l1_id].mass
        SR2L_l2_pt=eles[SR2L_l2_id].pt
        SR2L_l2_eta=eles[SR2L_l2_id].eta
        SR2L_l2_phi=eles[SR2L_l2_id].phi
        SR2L_l2_mass=eles[SR2L_l2_id].mass
        l1_v4.SetPtEtaPhiM(SR2L_l1_pt,SR2L_l1_eta,SR2L_l1_phi,SR2L_l1_mass)
        l2_v4.SetPtEtaPhiM(SR2L_l2_pt,SR2L_l2_eta,SR2L_l2_phi,SR2L_l2_mass)
        SR2L_zll_pt=(l1_v4+l2_v4).Pt()
        SR2L_zll_eta=(l1_v4+l2_v4).Eta()
        SR2L_zll_phi=(l1_v4+l2_v4).Phi()
        SR2L_zll_mass=(l1_v4+l2_v4).M()

      zll_v4=l1_v4+l2_v4

      best_combo = assign_jets_from_many(TightJet_v4)
      w_j1, w_j2=best_combo['W_pair']
      z_j1, z_j2=best_combo['Z_pair']
      SR2L_zj1_id=TightJet_id[z_j1]
      SR2L_zj1_pt=TightJet_v4[z_j1].Pt()
      SR2L_zj1_eta=TightJet_v4[z_j1].Eta()
      SR2L_zj1_phi=TightJet_v4[z_j1].Phi()
      SR2L_zj1_mass=TightJet_v4[z_j1].M()
      SR2L_zj2_id=TightJet_id[z_j2]
      SR2L_zj2_pt=TightJet_v4[z_j2].Pt()
      SR2L_zj2_eta=TightJet_v4[z_j2].Eta()
      SR2L_zj2_phi=TightJet_v4[z_j2].Phi()
      SR2L_zj2_mass=TightJet_v4[z_j2].M()
      SR2L_wj1_id=TightJet_id[w_j1]
      SR2L_wj1_pt=TightJet_v4[w_j1].Pt()
      SR2L_wj1_eta=TightJet_v4[w_j1].Eta()
      SR2L_wj1_phi=TightJet_v4[w_j1].Phi()
      SR2L_wj1_mass=TightJet_v4[w_j1].M()
      SR2L_wj2_id=TightJet_id[w_j2]
      SR2L_wj2_pt=TightJet_v4[w_j2].Pt()
      SR2L_wj2_eta=TightJet_v4[w_j2].Eta()
      SR2L_wj2_phi=TightJet_v4[w_j2].Phi()
      SR2L_wj2_mass=TightJet_v4[w_j2].M()
      SR2L_zjj_pt=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).Pt()
      SR2L_zjj_eta=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).Eta()
      SR2L_zjj_phi=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).Phi()
      SR2L_zjj_mass=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).M()
      SR2L_wjj_pt=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).Pt()
      SR2L_wjj_eta=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).Eta()
      SR2L_wjj_phi=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).Phi()
      SR2L_wjj_mass=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).M()
      
      zjj_v4=(TightJet_v4[z_j1]+TightJet_v4[z_j2])
      wjj_v4=(TightJet_v4[w_j1]+TightJet_v4[w_j2])

      SR2L_dR_zll_wjj = zll_v4.DeltaR(wjj_v4)
      SR2L_dEta_zll_wjj = abs(zll_v4.Eta()-wjj_v4.Eta())
      SR2L_dPhi_zll_wjj = zll_v4.DeltaPhi(wjj_v4)
      SR2L_zllwjj_pt = (zll_v4+wjj_v4).Pt()
      SR2L_zllwjj_eta = (zll_v4+wjj_v4).Eta()
      SR2L_zllwjj_phi = (zll_v4+wjj_v4).Phi()
      SR2L_zllwjj_mass = (zll_v4+wjj_v4).M()
      SR2L_dR_zll_zjj = zll_v4.DeltaR(zjj_v4)
      SR2L_dEta_zll_zjj = abs(zll_v4.Eta()-zjj_v4.Eta())
      SR2L_dPhi_zll_zjj = zll_v4.DeltaPhi(zjj_v4)
      SR2L_zllzjj_pt = (zll_v4+zjj_v4).Pt()
      SR2L_zllzjj_eta = (zll_v4+zjj_v4).Eta()
      SR2L_zllzjj_phi = (zll_v4+zjj_v4).Phi()
      SR2L_zllzjj_mass = (zll_v4+zjj_v4).M()
      SR2L_dR_wjj_zjj = wjj_v4.DeltaR(zjj_v4)
      SR2L_dEta_wjj_zjj = abs(wjj_v4.Eta()-zjj_v4.Eta())
      SR2L_dPhi_wjj_zjj = wjj_v4.DeltaPhi(zjj_v4)
      SR2L_wjjzjj_pt = (wjj_v4+zjj_v4).Pt()
      SR2L_wjjzjj_eta = (wjj_v4+zjj_v4).Eta()
      SR2L_wjjzjj_phi = (wjj_v4+zjj_v4).Phi()
      SR2L_wjjzjj_mass = (wjj_v4+zjj_v4).M()
      SR2L_wzz_pt = (zll_v4+wjj_v4+zjj_v4).Pt()
      SR2L_wzz_eta = (zll_v4+wjj_v4+zjj_v4).Eta()
      SR2L_wzz_phi = (zll_v4+wjj_v4+zjj_v4).Phi()
      SR2L_wzz_mass = (zll_v4+wjj_v4+zjj_v4).M()

    self.out.fillBranch("SR_2L_FlavorBit", SR_2L_FlavorBit)
    self.out.fillBranch("SR_2L_FakeBit", SR_2L_FakeBit)
    self.out.fillBranch("SR2L_l1_id", SR2L_l1_id)
    self.out.fillBranch("SR2L_l2_id", SR2L_l2_id)
    self.out.fillBranch("SR2L_l1_pt", SR2L_l1_pt)
    self.out.fillBranch("SR2L_l1_eta", SR2L_l1_eta)
    self.out.fillBranch("SR2L_l1_phi", SR2L_l1_phi)
    self.out.fillBranch("SR2L_l1_mass", SR2L_l1_mass)
    self.out.fillBranch("SR2L_l2_pt", SR2L_l2_pt)
    self.out.fillBranch("SR2L_l2_eta", SR2L_l2_eta)
    self.out.fillBranch("SR2L_l2_phi", SR2L_l2_phi)
    self.out.fillBranch("SR2L_l2_mass", SR2L_l2_mass)
    self.out.fillBranch("SR2L_zll_pt", SR2L_zll_pt)
    self.out.fillBranch("SR2L_zll_eta", SR2L_zll_eta)
    self.out.fillBranch("SR2L_zll_phi", SR2L_zll_phi)
    self.out.fillBranch("SR2L_zll_mass", SR2L_zll_mass)
    self.out.fillBranch("SR2L_zj1_id", SR2L_zj1_id)
    self.out.fillBranch("SR2L_zj1_pt", SR2L_zj1_pt)
    self.out.fillBranch("SR2L_zj1_eta", SR2L_zj1_eta)
    self.out.fillBranch("SR2L_zj1_phi", SR2L_zj1_phi)
    self.out.fillBranch("SR2L_zj1_mass", SR2L_zj1_mass)
    self.out.fillBranch("SR2L_zj2_id", SR2L_zj2_id)
    self.out.fillBranch("SR2L_zj2_pt", SR2L_zj2_pt)
    self.out.fillBranch("SR2L_zj2_eta", SR2L_zj2_eta)
    self.out.fillBranch("SR2L_zj2_phi", SR2L_zj2_phi)
    self.out.fillBranch("SR2L_zj2_mass", SR2L_zj2_mass)
    self.out.fillBranch("SR2L_wj1_id", SR2L_wj1_id)
    self.out.fillBranch("SR2L_wj1_pt", SR2L_wj1_pt)
    self.out.fillBranch("SR2L_wj1_eta", SR2L_wj1_eta)
    self.out.fillBranch("SR2L_wj1_phi", SR2L_wj1_phi)
    self.out.fillBranch("SR2L_wj1_mass", SR2L_wj1_mass)
    self.out.fillBranch("SR2L_wj2_id", SR2L_wj2_id)
    self.out.fillBranch("SR2L_wj2_pt", SR2L_wj2_pt)
    self.out.fillBranch("SR2L_wj2_eta", SR2L_wj2_eta)
    self.out.fillBranch("SR2L_wj2_phi", SR2L_wj2_phi)
    self.out.fillBranch("SR2L_wj2_mass", SR2L_wj2_mass)
    self.out.fillBranch("SR2L_zjj_pt", SR2L_zjj_pt)
    self.out.fillBranch("SR2L_zjj_eta", SR2L_zjj_eta)
    self.out.fillBranch("SR2L_zjj_phi", SR2L_zjj_phi)
    self.out.fillBranch("SR2L_zjj_mass", SR2L_zjj_mass)
    self.out.fillBranch("SR2L_wjj_pt", SR2L_wjj_pt)
    self.out.fillBranch("SR2L_wjj_eta", SR2L_wjj_eta)
    self.out.fillBranch("SR2L_wjj_phi", SR2L_wjj_phi)
    self.out.fillBranch("SR2L_wjj_mass", SR2L_wjj_mass)
    self.out.fillBranch("SR2L_dR_zll_wjj", SR2L_dR_zll_wjj)
    self.out.fillBranch("SR2L_dEta_zll_wjj", SR2L_dEta_zll_wjj)
    self.out.fillBranch("SR2L_dPhi_zll_wjj", SR2L_dPhi_zll_wjj)
    self.out.fillBranch("SR2L_zllwjj_pt", SR2L_zllwjj_pt)
    self.out.fillBranch("SR2L_zllwjj_eta", SR2L_zllwjj_eta)
    self.out.fillBranch("SR2L_zllwjj_phi", SR2L_zllwjj_phi)
    self.out.fillBranch("SR2L_zllwjj_mass", SR2L_zllwjj_mass)
    self.out.fillBranch("SR2L_dR_zll_zjj", SR2L_dR_zll_zjj)
    self.out.fillBranch("SR2L_dEta_zll_zjj", SR2L_dEta_zll_zjj)
    self.out.fillBranch("SR2L_dPhi_zll_zjj", SR2L_dPhi_zll_zjj)
    self.out.fillBranch("SR2L_zllzjj_pt", SR2L_zllzjj_pt)
    self.out.fillBranch("SR2L_zllzjj_eta", SR2L_zllzjj_eta)
    self.out.fillBranch("SR2L_zllzjj_phi", SR2L_zllzjj_phi)
    self.out.fillBranch("SR2L_zllzjj_mass", SR2L_zllzjj_mass)
    self.out.fillBranch("SR2L_dR_wjj_zjj", SR2L_dR_wjj_zjj)
    self.out.fillBranch("SR2L_dEta_wjj_zjj", SR2L_dEta_wjj_zjj)
    self.out.fillBranch("SR2L_dPhi_wjj_zjj", SR2L_dPhi_wjj_zjj)
    self.out.fillBranch("SR2L_wjjzjj_pt", SR2L_wjjzjj_pt)
    self.out.fillBranch("SR2L_wjjzjj_eta", SR2L_wjjzjj_eta)
    self.out.fillBranch("SR2L_wjjzjj_phi", SR2L_wjjzjj_phi)
    self.out.fillBranch("SR2L_wjjzjj_mass", SR2L_wjjzjj_mass)
    self.out.fillBranch("SR2L_wzz_pt", SR2L_wzz_pt)
    self.out.fillBranch("SR2L_wzz_eta", SR2L_wzz_eta)
    self.out.fillBranch("SR2L_wzz_phi", SR2L_wzz_phi)
    self.out.fillBranch("SR2L_wzz_mass", SR2L_wzz_mass)
    # Three leptons channel
    # SR_3L_FlavorBit
    #  l3 | l2 | l1 (l1 and l2 are always considered as from Z decay.)
    # mmm: SR_3L_FlavorBit=39
    # emm: SR_3L_FlavorBit=37
    # eee: SR_3L_FlavorBit=33
    # mee: SR_3L_FlavorBit=35
    
    # SR_3L_FakeBit, bit1=0 if l1 is fake, bit1=1 if l1 is prompt, bit2=0 if l2 is fake, bit2=1 if l2 is prompt, bit3=0 if l3 is fake, bit3=1 if l3 is prompt
    # SR_3L_FakeBit=7: three prompt
    # SR_3L_FakeBit=6: l1 is fake
    # SR_3L_FakeBit=5: l2 is fake
    # SR_3L_FakeBit=4: l1 and l2 are fake
    # SR_3L_FakeBit=3: l3 is fake
    # SR_3L_FakeBit=2: l1 and l3 are fake
    # SR_3L_FakeBit=1: l2 and l3 are fake
    # SR_3L_FakeBit=0: all are fake
    SR_3L_FlavorBit=-1
    SR_3L_FakeBit=-1
    SR3L_zl1_id=-99
    SR3L_zl2_id=-99
    SR3L_wl_id=-99
    SR3L_zl1_pt=-99
    SR3L_zl1_eta=-99
    SR3L_zl1_phi=-99
    SR3L_zl1_mass=-99
    SR3L_zl2_pt=-99
    SR3L_zl2_eta=-99
    SR3L_zl2_phi=-99
    SR3L_zl2_mass=-99
    SR3L_wl_pt=-99
    SR3L_wl_eta=-99
    SR3L_wl_phi=-99
    SR3L_wl_mass=-99
    SR3L_zll_pt=-99
    SR3L_zll_eta=-99
    SR3L_zll_phi=-99
    SR3L_zll_mass=-99
    SR3L_wlv_pt=-99
    SR3L_wlv_eta=-99
    SR3L_wlv_phi=-99
    SR3L_wlv_mass=-99
    SR3L_zj1_id=-99
    SR3L_zj1_pt=-99
    SR3L_zj1_eta=-99
    SR3L_zj1_phi=-99
    SR3L_zj1_mass=-99
    SR3L_zj2_id=-99
    SR3L_zj2_pt=-99
    SR3L_zj2_eta=-99
    SR3L_zj2_phi=-99
    SR3L_zj2_mass=-99
    SR3L_zjj_pt=-99
    SR3L_zjj_eta=-99
    SR3L_zjj_phi=-99
    SR3L_zjj_mass=-99

    SR3L_dR_zll_wlv = -99
    SR3L_dEta_zll_wlv = -99
    SR3L_dPhi_zll_wlv = -99
    SR3L_zllwlv_pt = -99
    SR3L_zllwlv_eta = -99
    SR3L_zllwlv_phi = -99
    SR3L_zllwlv_mass = -99
    SR3L_dR_zll_zjj = -99
    SR3L_dEta_zll_zjj = -99
    SR3L_dPhi_zll_zjj = -99
    SR3L_zllzjj_pt = -99
    SR3L_zllzjj_eta = -99
    SR3L_zllzjj_phi = -99
    SR3L_zllzjj_mass = -99
    SR3L_dR_wlv_zjj = -99
    SR3L_dEta_wlv_zjj = -99
    SR3L_dPhi_wlv_zjj = -99
    SR3L_wlvzjj_pt = -99
    SR3L_wlvzjj_eta = -99
    SR3L_wlvzjj_phi = -99
    SR3L_wlvzjj_mass = -99
    SR3L_wzz_pt = -99
    SR3L_wzz_eta = -99
    SR3L_wzz_phi = -99
    SR3L_wzz_mass = -99

    if (len(TFMuons_id)+len(TFElectrons_id))==3:
      if len(TightJet_id)<2:return False
      l1_v4=TLorentzVector()
      l2_v4=TLorentzVector()
      l3_v4=TLorentzVector()
      zll_v4=TLorentzVector()
      zjj_v4=TLorentzVector()
      wlv_v4=TLorentzVector()

      if len(TMuons_id)==3:
        chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge
        if not abs(chargeTot_tmp)==1:return False
        SR_3L_FlavorBit=39
        SR_3L_FakeBit=7
        l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
        l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
        l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
        zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,TMuons_id[0],TMuons_id[1],TMuons_id[2])
        SR3L_zl1_pt=zl1_v4.Pt()
        SR3L_zl1_eta=zl1_v4.Eta()
        SR3L_zl1_phi=zl1_v4.Phi()
        SR3L_zl1_mass=zl1_v4.M()
        SR3L_zl2_pt=zl2_v4.Pt()
        SR3L_zl2_eta=zl2_v4.Eta()
        SR3L_zl2_phi=zl2_v4.Phi()
        SR3L_zl2_mass=zl2_v4.M()
        SR3L_wl_pt=wl_v4.Pt()
        SR3L_wl_eta=wl_v4.Eta()
        SR3L_wl_phi=wl_v4.Phi()
        SR3L_wl_mass=wl_v4.M()
        SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
        SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
        SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
        SR3L_zll_mass=(zl1_v4+zl2_v4).M()
        zll_v4=zl1_v4+zl2_v4
#        wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
        wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
      elif len(TMuons_id)==2:
        if len(FMuons_id)==1:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[FMuons_id[0]].charge
          if not abs(chargeTot_tmp)==1:return False
          SR_3L_FlavorBit=39

          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,TMuons_id[0],TMuons_id[1],FMuons_id[0])
          if SR3L_zl1_id==FMuons_id[0]:
            SR_3L_FakeBit=6
          elif SR3L_zl2_id==FMuons_id[0]:
            SR_3L_FakeBit=5
          else:
            SR_3L_FakeBit=3
          SR3L_zl1_pt=zl1_v4.Pt()
          SR3L_zl1_eta=zl1_v4.Eta()
          SR3L_zl1_phi=zl1_v4.Phi()
          SR3L_zl1_mass=zl1_v4.M()
          SR3L_zl2_pt=zl2_v4.Pt()
          SR3L_zl2_eta=zl2_v4.Eta()
          SR3L_zl2_phi=zl2_v4.Phi()
          SR3L_zl2_mass=zl2_v4.M()
          SR3L_wl_pt=wl_v4.Pt()
          SR3L_wl_eta=wl_v4.Eta()
          SR3L_wl_phi=wl_v4.Phi()
          SR3L_wl_mass=wl_v4.M()
          SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
          SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
          SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
          SR3L_zll_mass=(zl1_v4+zl2_v4).M()
          zll_v4=zl1_v4+zl2_v4
          #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
          wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
        elif len(TElectrons_id)==1:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge
          if not abs(chargeTot_tmp)==0:return False
          SR_3L_FlavorBit=37
          SR_3L_FakeBit=7
          SR3L_zl1_id=TMuons_id[0]
          SR3L_zl2_id=TMuons_id[1]
          SR3L_wl_id=TElectrons_id[0]
          SR3L_zl1_pt=muons[TMuons_id[0]].pt
          SR3L_zl1_eta=muons[TMuons_id[0]].eta
          SR3L_zl1_phi=muons[TMuons_id[0]].phi
          SR3L_zl1_mass=muons[TMuons_id[0]].mass
          SR3L_zl2_pt=muons[TMuons_id[1]].pt
          SR3L_zl2_eta=muons[TMuons_id[1]].eta
          SR3L_zl2_phi=muons[TMuons_id[1]].phi
          SR3L_zl2_mass=muons[TMuons_id[1]].mass
          SR3L_wl_pt=eles[TElectrons_id[0]].pt
          SR3L_wl_eta=eles[TElectrons_id[0]].eta
          SR3L_wl_phi=eles[TElectrons_id[0]].phi
          SR3L_wl_mass=eles[TElectrons_id[0]].mass
          l1_v4.SetPtEtaPhiM(SR3L_zl1_pt,SR3L_zl1_eta,SR3L_zl1_phi,SR3L_zl1_mass)
          l2_v4.SetPtEtaPhiM(SR3L_zl2_pt,SR3L_zl2_eta,SR3L_zl2_phi,SR3L_zl2_mass)
          l3_v4.SetPtEtaPhiM(SR3L_wl_pt,SR3L_wl_eta,SR3L_wl_phi,SR3L_wl_mass)
          SR3L_zll_pt=(l1_v4+l2_v4).Pt()
          SR3L_zll_eta=(l1_v4+l2_v4).Eta()
          SR3L_zll_phi=(l1_v4+l2_v4).Phi()
          SR3L_zll_mass=(l1_v4+l2_v4).M()
          zll_v4=l1_v4+l2_v4
          #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
          wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
        elif len(FElectrons_id)==1:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge
          if not abs(chargeTot_tmp)==0:return False
          SR_3L_FlavorBit=37
          SR_3L_FakeBit=3
          SR3L_zl1_id=TMuons_id[0]
          SR3L_zl2_id=TMuons_id[1]
          SR3L_wl_id=FElectrons_id[0]
          SR3L_zl1_pt=muons[TMuons_id[0]].pt
          SR3L_zl1_eta=muons[TMuons_id[0]].eta
          SR3L_zl1_phi=muons[TMuons_id[0]].phi
          SR3L_zl1_mass=muons[TMuons_id[0]].mass
          SR3L_zl2_pt=muons[TMuons_id[1]].pt
          SR3L_zl2_eta=muons[TMuons_id[1]].eta
          SR3L_zl2_phi=muons[TMuons_id[1]].phi
          SR3L_zl2_mass=muons[TMuons_id[1]].mass
          SR3L_wl_pt=eles[FElectrons_id[0]].pt
          SR3L_wl_eta=eles[FElectrons_id[0]].eta
          SR3L_wl_phi=eles[FElectrons_id[0]].phi
          SR3L_wl_mass=eles[FElectrons_id[0]].mass
          l1_v4.SetPtEtaPhiM(SR3L_zl1_pt,SR3L_zl1_eta,SR3L_zl1_phi,SR3L_zl1_mass)
          l2_v4.SetPtEtaPhiM(SR3L_zl2_pt,SR3L_zl2_eta,SR3L_zl2_phi,SR3L_zl2_mass)
          l3_v4.SetPtEtaPhiM(SR3L_wl_pt,SR3L_wl_eta,SR3L_wl_phi,SR3L_wl_mass)
          SR3L_zll_pt=(l1_v4+l2_v4).Pt()
          SR3L_zll_eta=(l1_v4+l2_v4).Eta()
          SR3L_zll_phi=(l1_v4+l2_v4).Phi()
          SR3L_zll_mass=(l1_v4+l2_v4).M()
          zll_v4=l1_v4+l2_v4
          #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
          wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
      elif len(TMuons_id)==1:
        if len(FMuons_id)==2:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
          if not abs(chargeTot_tmp)==1:return False
          SR_3L_FlavorBit=39
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[TMuons_id[0]].charge,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,TMuons_id[0],FMuons_id[0],FMuons_id[1])
          if SR3L_zl1_id==TMuons_id[0]:
            SR_3L_FakeBit=1
          elif SR3L_zl2_id==TMuons_id[0]:
            SR_3L_FakeBit=2
          else:
            SR_3L_FakeBit=4
          SR3L_zl1_pt=zl1_v4.Pt()
          SR3L_zl1_eta=zl1_v4.Eta()
          SR3L_zl1_phi=zl1_v4.Phi()
          SR3L_zl1_mass=zl1_v4.M()
          SR3L_zl2_pt=zl2_v4.Pt()
          SR3L_zl2_eta=zl2_v4.Eta()
          SR3L_zl2_phi=zl2_v4.Phi()
          SR3L_zl2_mass=zl2_v4.M()
          SR3L_wl_pt=wl_v4.Pt()
          SR3L_wl_eta=wl_v4.Eta()
          SR3L_wl_phi=wl_v4.Phi()
          SR3L_wl_mass=wl_v4.M()
          SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
          SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
          SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
          SR3L_zll_mass=(zl1_v4+zl2_v4).M()
          zll_v4=zl1_v4+zl2_v4
          #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
          wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
        elif len(FMuons_id)==1:
          SR_3L_FlavorBit=37
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge
          if not abs(chargeTot_tmp)==0:return False
          if len(TElectrons_id)==1:
            SR_3L_FakeBit=5
            SR3L_zl1_id=TMuons_id[0]
            SR3L_zl2_id=FMuons_id[0]
            SR3L_wl_id=TElectrons_id[0]
            l1_v4.SetPtEtaPhiM(muons[SR3L_zl1_id].pt,muons[SR3L_zl1_id].eta,muons[SR3L_zl1_id].phi,muons[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(muons[SR3L_zl2_id].pt,muons[SR3L_zl2_id].eta,muons[SR3L_zl2_id].phi,muons[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(eles[SR3L_wl_id].pt,eles[SR3L_wl_id].eta,eles[SR3L_wl_id].phi,eles[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
          else:
            SR_3L_FakeBit=1
            SR3L_zl1_id=TMuons_id[0]
            SR3L_zl2_id=FMuons_id[0]
            SR3L_wl_id=FElectrons_id[0]
            l1_v4.SetPtEtaPhiM(muons[SR3L_zl1_id].pt,muons[SR3L_zl1_id].eta,muons[SR3L_zl1_id].phi,muons[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(muons[SR3L_zl2_id].pt,muons[SR3L_zl2_id].eta,muons[SR3L_zl2_id].phi,muons[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(eles[SR3L_wl_id].pt,eles[SR3L_wl_id].eta,eles[SR3L_wl_id].phi,eles[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
        elif len(FMuons_id)==0:
          SR_3L_FlavorBit=35
          if len(TElectrons_id)==2:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_3L_FakeBit=7
            SR3L_zl1_id=TElectrons_id[0]
            SR3L_zl2_id=TElectrons_id[1]
            SR3L_wl_id=TMuons_id[0]
            l1_v4.SetPtEtaPhiM(eles[SR3L_zl1_id].pt,eles[SR3L_zl1_id].eta,eles[SR3L_zl1_id].phi,eles[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(eles[SR3L_zl2_id].pt,eles[SR3L_zl2_id].eta,eles[SR3L_zl2_id].phi,eles[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(muons[SR3L_wl_id].pt,muons[SR3L_wl_id].eta,muons[SR3L_wl_id].phi,muons[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
          elif len(TElectrons_id)==1:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_3L_FakeBit=5
            SR3L_zl1_id=TElectrons_id[0]
            SR3L_zl2_id=FElectrons_id[0]
            SR3L_wl_id=TMuons_id[0]
            l1_v4.SetPtEtaPhiM(eles[SR3L_zl1_id].pt,eles[SR3L_zl1_id].eta,eles[SR3L_zl1_id].phi,eles[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(eles[SR3L_zl2_id].pt,eles[SR3L_zl2_id].eta,eles[SR3L_zl2_id].phi,eles[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(muons[SR3L_wl_id].pt,muons[SR3L_wl_id].eta,muons[SR3L_wl_id].phi,muons[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
          else:
            chargeTot_tmp=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_3L_FakeBit=4
            SR3L_zl1_id=FElectrons_id[0]
            SR3L_zl2_id=FElectrons_id[1]
            SR3L_wl_id=TMuons_id[0]
            l1_v4.SetPtEtaPhiM(eles[SR3L_zl1_id].pt,eles[SR3L_zl1_id].eta,eles[SR3L_zl1_id].phi,eles[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(eles[SR3L_zl2_id].pt,eles[SR3L_zl2_id].eta,eles[SR3L_zl2_id].phi,eles[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(muons[SR3L_wl_id].pt,muons[SR3L_wl_id].eta,muons[SR3L_wl_id].phi,muons[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
      elif len(TMuons_id)==0:
        if len(FMuons_id)==3:
          chargeTot_tmp=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge+muons[FMuons_id[2]].charge
          if not abs(chargeTot_tmp)==1:return False
          SR_3L_FlavorBit=39
          SR_3L_FakeBit=0
          l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[2]].pt,muons[FMuons_id[2]].eta,muons[FMuons_id[2]].phi,muons[FMuons_id[2]].mass)
          zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,muons[FMuons_id[2]].charge,FMuons_id[0],FMuons_id[1],FMuons_id[2])
          SR3L_zl1_pt=zl1_v4.Pt()
          SR3L_zl1_eta=zl1_v4.Eta()
          SR3L_zl1_phi=zl1_v4.Phi()
          SR3L_zl1_mass=zl1_v4.M()
          SR3L_zl2_pt=zl2_v4.Pt()
          SR3L_zl2_eta=zl2_v4.Eta()
          SR3L_zl2_phi=zl2_v4.Phi()
          SR3L_zl2_mass=zl2_v4.M()
          SR3L_wl_pt=wl_v4.Pt()
          SR3L_wl_eta=wl_v4.Eta()
          SR3L_wl_phi=wl_v4.Phi()
          SR3L_wl_mass=wl_v4.M()
          SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
          SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
          SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
          SR3L_zll_mass=(zl1_v4+zl2_v4).M()
          zll_v4=zl1_v4+zl2_v4
          #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
          wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
        elif len(FMuons_id)==2:
          chargeTot_tmp=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
          if not abs(chargeTot_tmp)==0:return False
          SR_3L_FlavorBit=37
          if len(TElectrons_id)==1:
            SR3L_zl1_id=FMuons_id[0]
            SR3L_zl2_id=FMuons_id[1]
            SR3L_wl_id=TElectrons_id[0]
            SR_3L_FakeBit=4
            l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
          else:
            SR_3L_FakeBit=0
            SR3L_zl1_id=FMuons_id[0]
            SR3L_zl2_id=FMuons_id[1]
            SR3L_wl_id=FElectrons_id[0]
            l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
        elif len(FMuons_id)==1:
          SR_3L_FlavorBit=35
          if len(TElectrons_id)==2:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_3L_FakeBit=3
            SR3L_zl1_id=TElectrons_id[0]
            SR3L_zl2_id=TElectrons_id[1]
            SR3L_wl_id=FMuons_id[0]
            l1_v4.SetPtEtaPhiM(eles[SR3L_zl1_id].pt,eles[SR3L_zl1_id].eta,eles[SR3L_zl1_id].phi,eles[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(eles[SR3L_zl2_id].pt,eles[SR3L_zl2_id].eta,eles[SR3L_zl2_id].phi,eles[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(muons[SR3L_wl_id].pt,muons[SR3L_wl_id].eta,muons[SR3L_wl_id].phi,muons[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
          elif len(TElectrons_id)==1:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_3L_FakeBit=1
            SR3L_zl1_id=TElectrons_id[0]
            SR3L_zl2_id=FElectrons_id[0]
            SR3L_wl_id=FMuons_id[0]
            l1_v4.SetPtEtaPhiM(eles[SR3L_zl1_id].pt,eles[SR3L_zl1_id].eta,eles[SR3L_zl1_id].phi,eles[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(eles[SR3L_zl2_id].pt,eles[SR3L_zl2_id].eta,eles[SR3L_zl2_id].phi,eles[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(muons[SR3L_wl_id].pt,muons[SR3L_wl_id].eta,muons[SR3L_wl_id].phi,muons[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
          else:
            chargeTot_tmp=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_3L_FakeBit=0
            SR3L_zl1_id=FElectrons_id[0]
            SR3L_zl2_id=FElectrons_id[1]
            SR3L_wl_id=FMuons_id[0]
            l1_v4.SetPtEtaPhiM(eles[SR3L_zl1_id].pt,eles[SR3L_zl1_id].eta,eles[SR3L_zl1_id].phi,eles[SR3L_zl1_id].mass)
            l2_v4.SetPtEtaPhiM(eles[SR3L_zl2_id].pt,eles[SR3L_zl2_id].eta,eles[SR3L_zl2_id].phi,eles[SR3L_zl2_id].mass)
            l3_v4.SetPtEtaPhiM(muons[SR3L_wl_id].pt,muons[SR3L_wl_id].eta,muons[SR3L_wl_id].phi,muons[SR3L_wl_id].mass)
            SR3L_zl1_pt=l1_v4.Pt()
            SR3L_zl1_eta=l1_v4.Eta()
            SR3L_zl1_phi=l1_v4.Phi()
            SR3L_zl1_mass=l1_v4.M()
            SR3L_zl2_pt=l2_v4.Pt()
            SR3L_zl2_eta=l2_v4.Eta()
            SR3L_zl2_phi=l2_v4.Phi()
            SR3L_zl2_mass=l2_v4.M()
            SR3L_wl_pt=l3_v4.Pt()
            SR3L_wl_eta=l3_v4.Eta()
            SR3L_wl_phi=l3_v4.Phi()
            SR3L_wl_mass=l3_v4.M()
            SR3L_zll_pt=(l1_v4+l2_v4).Pt()
            SR3L_zll_eta=(l1_v4+l2_v4).Eta()
            SR3L_zll_phi=(l1_v4+l2_v4).Phi()
            SR3L_zll_mass=(l1_v4+l2_v4).M()
            zll_v4=l1_v4+l2_v4
            #wlv_v4=w_v4(l3_v4, met_user, met_phi_user)
            wlv_v4=w_v4(l3_v4, met_CorrePhi, met_phi_CorrePhi)
        else:
          SR_3L_FlavorBit=33
          if len(TElectrons_id)==3:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge
            if not abs(chargeTot_tmp)==1:return False
            SR_3L_FakeBit=7
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2])
            SR3L_zl1_pt=zl1_v4.Pt()
            SR3L_zl1_eta=zl1_v4.Eta()
            SR3L_zl1_phi=zl1_v4.Phi()
            SR3L_zl1_mass=zl1_v4.M()
            SR3L_zl2_pt=zl2_v4.Pt()
            SR3L_zl2_eta=zl2_v4.Eta()
            SR3L_zl2_phi=zl2_v4.Phi()
            SR3L_zl2_mass=zl2_v4.M()
            SR3L_wl_pt=wl_v4.Pt()
            SR3L_wl_eta=wl_v4.Eta()
            SR3L_wl_phi=wl_v4.Phi()
            SR3L_wl_mass=wl_v4.M()
            SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
            SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
            SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
            SR3L_zll_mass=(zl1_v4+zl2_v4).M()
            zll_v4=zl1_v4+zl2_v4
            #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
            wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
          elif len(TElectrons_id)==2:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp)==1:return False
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[FElectrons_id[0]].charge,TElectrons_id[0],TElectrons_id[1],FElectrons_id[0])
            if SR3L_zl1_id==FElectrons_id[0]:
              SR_3L_FakeBit=6
            elif SR3L_zl2_id==FElectrons_id[0]:
              SR_3L_FakeBit=5
            else:
              SR_3L_FakeBit=3
            SR3L_zl1_pt=zl1_v4.Pt()
            SR3L_zl1_eta=zl1_v4.Eta()
            SR3L_zl1_phi=zl1_v4.Phi()
            SR3L_zl1_mass=zl1_v4.M()
            SR3L_zl2_pt=zl2_v4.Pt()
            SR3L_zl2_eta=zl2_v4.Eta()
            SR3L_zl2_phi=zl2_v4.Phi()
            SR3L_zl2_mass=zl2_v4.M()
            SR3L_wl_pt=wl_v4.Pt()
            SR3L_wl_eta=wl_v4.Eta()
            SR3L_wl_phi=wl_v4.Phi()
            SR3L_wl_mass=wl_v4.M()
            SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
            SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
            SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
            SR3L_zll_mass=(zl1_v4+zl2_v4).M()
            zll_v4=zl1_v4+zl2_v4
            #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
            wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
          elif len(TElectrons_id)==1:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp)==1:return False
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
            zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,eles[TElectrons_id[0]].charge,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,TElectrons_id[0],FElectrons_id[0],FElectrons_id[1])
            if SR3L_zl1_id==TElectrons_id[0]:
              SR_3L_FakeBit=1
            elif SR3L_zl2_id==TElectrons_id[0]:
              SR_3L_FakeBit=2
            else:
              SR_3L_FakeBit=4
            SR3L_zl1_pt=zl1_v4.Pt()
            SR3L_zl1_eta=zl1_v4.Eta()
            SR3L_zl1_phi=zl1_v4.Phi()
            SR3L_zl1_mass=zl1_v4.M()
            SR3L_zl2_pt=zl2_v4.Pt()
            SR3L_zl2_eta=zl2_v4.Eta()
            SR3L_zl2_phi=zl2_v4.Phi()
            SR3L_zl2_mass=zl2_v4.M()
            SR3L_wl_pt=wl_v4.Pt()
            SR3L_wl_eta=wl_v4.Eta()
            SR3L_wl_phi=wl_v4.Phi()
            SR3L_wl_mass=wl_v4.M()
            SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
            SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
            SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
            SR3L_zll_mass=(zl1_v4+zl2_v4).M()
            zll_v4=zl1_v4+zl2_v4
            #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
            wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)
          else:
            SR_3L_FakeBit=0
            chargeTot_tmp=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge+eles[FElectrons_id[2]].charge
            if not abs(chargeTot_tmp)==1:return False
            l1_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[2]].pt,eles[FElectrons_id[2]].eta,eles[FElectrons_id[2]].phi,eles[FElectrons_id[2]].mass)
            zl1_v4,zl2_v4,wl_v4,SR3L_zl1_id,SR3L_zl2_id,SR3L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,eles[FElectrons_id[2]].charge,FElectrons_id[0],FElectrons_id[1],FElectrons_id[2])
            SR3L_zl1_pt=zl1_v4.Pt()
            SR3L_zl1_eta=zl1_v4.Eta()
            SR3L_zl1_phi=zl1_v4.Phi()
            SR3L_zl1_mass=zl1_v4.M()
            SR3L_zl2_pt=zl2_v4.Pt()
            SR3L_zl2_eta=zl2_v4.Eta()
            SR3L_zl2_phi=zl2_v4.Phi()
            SR3L_zl2_mass=zl2_v4.M()
            SR3L_wl_pt=wl_v4.Pt()
            SR3L_wl_eta=wl_v4.Eta()
            SR3L_wl_phi=wl_v4.Phi()
            SR3L_wl_mass=wl_v4.M()
            SR3L_zll_pt=(zl1_v4+zl2_v4).Pt()
            SR3L_zll_eta=(zl1_v4+zl2_v4).Eta()
            SR3L_zll_phi=(zl1_v4+zl2_v4).Phi()
            SR3L_zll_mass=(zl1_v4+zl2_v4).M()
            zll_v4=zl1_v4+zl2_v4
            #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
            wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)

      SR3L_wlv_pt=wlv_v4.Pt()
      SR3L_wlv_eta=wlv_v4.Eta()
      SR3L_wlv_phi=wlv_v4.Phi()
      SR3L_wlv_mass=wlv_v4.M()

      z_j1, z_j2 = assign_2jets_from_many(TightJet_v4,MZ)
      SR3L_zj1_id=TightJet_id[z_j1]
      SR3L_zj1_pt=TightJet_v4[z_j1].Pt()
      SR3L_zj1_eta=TightJet_v4[z_j1].Eta()
      SR3L_zj1_phi=TightJet_v4[z_j1].Phi()
      SR3L_zj1_mass=TightJet_v4[z_j1].M()
      SR3L_zj2_id=TightJet_id[z_j2]
      SR3L_zj2_pt=TightJet_v4[z_j2].Pt()
      SR3L_zj2_eta=TightJet_v4[z_j2].Eta()
      SR3L_zj2_phi=TightJet_v4[z_j2].Phi()
      SR3L_zj2_mass=TightJet_v4[z_j2].M()
      SR3L_zjj_pt=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).Pt()
      SR3L_zjj_eta=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).Eta()
      SR3L_zjj_phi=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).Phi()
      SR3L_zjj_mass=(TightJet_v4[z_j1]+TightJet_v4[z_j2]).M()

      zjj_v4=TightJet_v4[z_j1]+TightJet_v4[z_j2]

      SR3L_dR_zll_wlv = zll_v4.DeltaR(wlv_v4)
      SR3L_dEta_zll_wlv = abs(zll_v4.Eta()-wlv_v4.Eta())
      SR3L_dPhi_zll_wlv = zll_v4.DeltaPhi(wlv_v4)
      SR3L_zllwlv_pt = (zll_v4+wlv_v4).Pt()
      SR3L_zllwlv_eta = (zll_v4+wlv_v4).Eta()
      SR3L_zllwlv_phi = (zll_v4+wlv_v4).Phi()
      SR3L_zllwlv_mass = (zll_v4+wlv_v4).M()
      SR3L_dR_zll_zjj = zll_v4.DeltaR(zjj_v4)
      SR3L_dEta_zll_zjj = abs(zll_v4.Eta()-zjj_v4.Eta())
      SR3L_dPhi_zll_zjj = zll_v4.DeltaPhi(zjj_v4)
      SR3L_zllzjj_pt = (zll_v4+zjj_v4).Pt()
      SR3L_zllzjj_eta = (zll_v4+zjj_v4).Eta()
      SR3L_zllzjj_phi = (zll_v4+zjj_v4).Phi()
      SR3L_zllzjj_mass = (zll_v4+zjj_v4).M()
      SR3L_dR_wlv_zjj = wlv_v4.DeltaR(zjj_v4)
      SR3L_dEta_wlv_zjj = abs(wlv_v4.Eta()-zjj_v4.Eta())
      SR3L_dPhi_wlv_zjj = wlv_v4.DeltaPhi(zjj_v4)
      SR3L_wlvzjj_pt = (wlv_v4+zjj_v4).Pt()
      SR3L_wlvzjj_eta = (wlv_v4+zjj_v4).Eta()
      SR3L_wlvzjj_phi = (wlv_v4+zjj_v4).Phi()
      SR3L_wlvzjj_mass = (wlv_v4+zjj_v4).M()
      SR3L_wzz_pt = (zll_v4+wlv_v4+zjj_v4).Pt()
      SR3L_wzz_eta = (zll_v4+wlv_v4+zjj_v4).Eta()
      SR3L_wzz_phi = (zll_v4+wlv_v4+zjj_v4).Phi()
      SR3L_wzz_mass = (zll_v4+wlv_v4+zjj_v4).M()

    self.out.fillBranch("SR_3L_FlavorBit", SR_3L_FlavorBit)
    self.out.fillBranch("SR_3L_FakeBit", SR_3L_FakeBit)
    self.out.fillBranch("SR3L_zl1_id", SR3L_zl1_id)
    self.out.fillBranch("SR3L_zl2_id", SR3L_zl2_id)
    self.out.fillBranch("SR3L_wl_id", SR3L_wl_id)
    self.out.fillBranch("SR3L_zl1_pt", SR3L_zl1_pt)
    self.out.fillBranch("SR3L_zl1_eta", SR3L_zl1_eta)
    self.out.fillBranch("SR3L_zl1_phi", SR3L_zl1_phi)
    self.out.fillBranch("SR3L_zl1_mass", SR3L_zl1_mass)
    self.out.fillBranch("SR3L_zl2_pt", SR3L_zl2_pt)
    self.out.fillBranch("SR3L_zl2_eta", SR3L_zl2_eta)
    self.out.fillBranch("SR3L_zl2_phi", SR3L_zl2_phi)
    self.out.fillBranch("SR3L_zl2_mass", SR3L_zl2_mass)
    self.out.fillBranch("SR3L_wl_pt", SR3L_wl_pt)
    self.out.fillBranch("SR3L_wl_eta", SR3L_wl_eta)
    self.out.fillBranch("SR3L_wl_phi", SR3L_wl_phi)
    self.out.fillBranch("SR3L_wl_mass", SR3L_wl_mass)
    self.out.fillBranch("SR3L_zll_pt", SR3L_zll_pt)
    self.out.fillBranch("SR3L_zll_eta", SR3L_zll_eta)
    self.out.fillBranch("SR3L_zll_phi", SR3L_zll_phi)
    self.out.fillBranch("SR3L_zll_mass", SR3L_zll_mass)
    self.out.fillBranch("SR3L_wlv_pt", SR3L_wlv_pt)
    self.out.fillBranch("SR3L_wlv_eta", SR3L_wlv_eta)
    self.out.fillBranch("SR3L_wlv_phi", SR3L_wlv_phi)
    self.out.fillBranch("SR3L_wlv_mass", SR3L_wlv_mass)
    self.out.fillBranch("SR3L_zj1_id", SR3L_zj1_id)
    self.out.fillBranch("SR3L_zj1_pt", SR3L_zj1_pt)
    self.out.fillBranch("SR3L_zj1_eta", SR3L_zj1_eta)
    self.out.fillBranch("SR3L_zj1_phi", SR3L_zj1_phi)
    self.out.fillBranch("SR3L_zj1_mass", SR3L_zj1_mass)
    self.out.fillBranch("SR3L_zj2_id", SR3L_zj2_id)
    self.out.fillBranch("SR3L_zj2_pt", SR3L_zj2_pt)
    self.out.fillBranch("SR3L_zj2_eta", SR3L_zj2_eta)
    self.out.fillBranch("SR3L_zj2_phi", SR3L_zj2_phi)
    self.out.fillBranch("SR3L_zj2_mass", SR3L_zj2_mass)
    self.out.fillBranch("SR3L_zjj_pt", SR3L_zjj_pt)
    self.out.fillBranch("SR3L_zjj_eta", SR3L_zjj_eta)
    self.out.fillBranch("SR3L_zjj_phi", SR3L_zjj_phi)
    self.out.fillBranch("SR3L_zjj_mass", SR3L_zjj_mass)
    self.out.fillBranch("SR3L_dR_zll_wlv", SR3L_dR_zll_wlv)
    self.out.fillBranch("SR3L_dEta_zll_wlv", SR3L_dEta_zll_wlv)
    self.out.fillBranch("SR3L_dPhi_zll_wlv", SR3L_dPhi_zll_wlv)
    self.out.fillBranch("SR3L_zllwlv_pt", SR3L_zllwlv_pt)
    self.out.fillBranch("SR3L_zllwlv_eta", SR3L_zllwlv_eta)
    self.out.fillBranch("SR3L_zllwlv_phi", SR3L_zllwlv_phi)
    self.out.fillBranch("SR3L_zllwlv_mass", SR3L_zllwlv_mass)
    self.out.fillBranch("SR3L_dR_zll_zjj", SR3L_dR_zll_zjj)
    self.out.fillBranch("SR3L_dEta_zll_zjj", SR3L_dEta_zll_zjj)
    self.out.fillBranch("SR3L_dPhi_zll_zjj", SR3L_dPhi_zll_zjj)
    self.out.fillBranch("SR3L_zllzjj_pt", SR3L_zllzjj_pt)
    self.out.fillBranch("SR3L_zllzjj_eta", SR3L_zllzjj_eta)
    self.out.fillBranch("SR3L_zllzjj_phi", SR3L_zllzjj_phi)
    self.out.fillBranch("SR3L_zllzjj_mass", SR3L_zllzjj_mass)
    self.out.fillBranch("SR3L_dR_wlv_zjj", SR3L_dR_wlv_zjj)
    self.out.fillBranch("SR3L_dEta_wlv_zjj", SR3L_dEta_wlv_zjj)
    self.out.fillBranch("SR3L_dPhi_wlv_zjj", SR3L_dPhi_wlv_zjj)
    self.out.fillBranch("SR3L_wlvzjj_pt", SR3L_wlvzjj_pt)
    self.out.fillBranch("SR3L_wlvzjj_eta", SR3L_wlvzjj_eta)
    self.out.fillBranch("SR3L_wlvzjj_phi", SR3L_wlvzjj_phi)
    self.out.fillBranch("SR3L_wlvzjj_mass", SR3L_wlvzjj_mass)
    self.out.fillBranch("SR3L_wzz_pt", SR3L_wzz_pt)
    self.out.fillBranch("SR3L_wzz_eta", SR3L_wzz_eta)
    self.out.fillBranch("SR3L_wzz_phi", SR3L_wzz_phi)
    self.out.fillBranch("SR3L_wzz_mass", SR3L_wzz_mass)

    # Four leptons channel
    # SR_4L_FlavorBit=52: 4mu, 
    # SR_4L_FlavorBit=44: 4ele,
    # SR_4L_FlavorBit=48: 2ele2mu,

    #     z2     |     z1     |
    # bit4  bit3 | bit2  bit1 |
    #  l4    l3  |  l2    l1  |

    SR_4L_FlavorBit=-1
    SR_4L_FakeBit=-1
    SR4L_z1l1_id=-99
    SR4L_z1l2_id=-99
    SR4L_z2l1_id=-99
    SR4L_z2l2_id=-99
    SR4L_z1l1_pt=-99
    SR4L_z1l1_eta=-99
    SR4L_z1l1_phi=-99
    SR4L_z1l1_mass=-99
    SR4L_z1l2_pt=-99
    SR4L_z1l2_eta=-99
    SR4L_z1l2_phi=-99
    SR4L_z1l2_mass=-99
    SR4L_z2l1_pt=-99
    SR4L_z2l1_eta=-99
    SR4L_z2l1_phi=-99
    SR4L_z2l1_mass=-99
    SR4L_z2l2_pt=-99
    SR4L_z2l2_eta=-99
    SR4L_z2l2_phi=-99
    SR4L_z2l2_mass=-99
    SR4L_z1_pt=-99
    SR4L_z1_eta=-99
    SR4L_z1_phi=-99
    SR4L_z1_mass=-99
    SR4L_z2_pt=-99
    SR4L_z2_eta=-99
    SR4L_z2_phi=-99
    SR4L_z2_mass=-99
    SR4L_wj1_id=-99
    SR4L_wj1_pt=-99
    SR4L_wj1_eta=-99
    SR4L_wj1_phi=-99
    SR4L_wj1_mass=-99
    SR4L_wj2_id=-99
    SR4L_wj2_pt=-99
    SR4L_wj2_eta=-99
    SR4L_wj2_phi=-99
    SR4L_wj2_mass=-99
    SR4L_wjj_pt=-99
    SR4L_wjj_eta=-99
    SR4L_wjj_phi=-99
    SR4L_wjj_mass=-99
    SR4L_dR_z1_wjj = -99
    SR4L_dEta_z1_wjj =-99
    SR4L_dPhi_z1_wjj =-99
    SR4L_z1wjj_pt = -99
    SR4L_z1wjj_eta =-99
    SR4L_z1wjj_phi =-99
    SR4L_z1wjj_mass = -99
    SR4L_dR_z1_z2 = -99
    SR4L_dEta_z1_z2 = -99
    SR4L_dPhi_z1_z2 = -99
    SR4L_z1z2_pt = -99
    SR4L_z1z2_eta =-99
    SR4L_z1z2_phi =-99
    SR4L_z1z2_mass =-99
    SR4L_dR_wjj_z2 =-99
    SR4L_dEta_wjj_z2 =-99
    SR4L_dPhi_wjj_z2 =-99
    SR4L_wjjz2_pt = -99
    SR4L_wjjz2_eta =-99
    SR4L_wjjz2_phi = -99
    SR4L_wjjz2_mass =-99
    SR4L_wzz_pt = -99
    SR4L_wzz_eta =-99
    SR4L_wzz_phi =-99
    SR4L_wzz_mass = -99

    if (len(TFMuons_id)+len(TFElectrons_id))==4:
      if len(TFMuons_id)==1 or len(TFMuons_id)==3: return False
      if len(TightJet_id)<2:return False
      l1_v4=TLorentzVector()
      l2_v4=TLorentzVector()
      l3_v4=TLorentzVector()
      l4_v4=TLorentzVector()
      z1_v4=TLorentzVector()
      z2_v4=TLorentzVector()
      wjj_v4=TLorentzVector()

      if len(TFMuons_id)==4:
        SR_4L_FlavorBit=52
        if len(TMuons_id)==4:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[TMuons_id[3]].charge
          if not abs(chargeTot_tmp)==0:return False
          SR_4L_FakeBit=15
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[TMuons_id[3]].pt,muons[TMuons_id[3]].eta,muons[TMuons_id[3]].phi,muons[TMuons_id[3]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,muons[TMuons_id[3]].charge,TMuons_id[0],TMuons_id[1],TMuons_id[2],TMuons_id[3])


        elif len(TMuons_id)==3:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[FMuons_id[0]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,muons[FMuons_id[0]].charge,TMuons_id[0],TMuons_id[1],TMuons_id[2],FMuons_id[0])
          if SR4L_z1l1_id==FMuons_id[0]:
            SR_4L_FakeBit=14
          elif SR4L_z1l2_id==FMuons_id[0]:
            SR_4L_FakeBit=13
          elif SR4L_z2l1_id==FMuons_id[0]:
            SR_4L_FakeBit=11
          else:
            SR_4L_FakeBit=7

        elif len(TMuons_id)==2:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,TMuons_id[0],TMuons_id[1],FMuons_id[0],FMuons_id[1])
          id_tmp=[SR4L_z2l2_id,SR4L_z2l1_id,SR4L_z1l2_id,SR4L_z1l1_id]
          real_index1=id_tmp.index(TMuons_id[0])
          real_index2=id_tmp.index(TMuons_id[1])
          SR_4L_FakeBit = (1 << (3 - real_index1)) + (1 << (3 - real_index2))

        elif len(TMuons_id)==1:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge+muons[FMuons_id[2]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[2]].pt,muons[FMuons_id[2]].eta,muons[FMuons_id[2]].phi,muons[FMuons_id[2]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,muons[FMuons_id[2]].charge,TMuons_id[0],FMuons_id[0],FMuons_id[1],FMuons_id[2])
          if SR4L_z1l1_id==TMuons_id[0]:
            SR_4L_FakeBit=1
          elif SR4L_z1l2_id==TMuons_id[0]:
            SR_4L_FakeBit=2
          elif SR4L_z2l1_id==TMuons_id[0]:
            SR_4L_FakeBit=4
          else:
            SR_4L_FakeBit=8
        else:
          chargeTot_tmp=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge+muons[FMuons_id[2]].charge+muons[FMuons_id[3]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[2]].pt,muons[FMuons_id[2]].eta,muons[FMuons_id[2]].phi,muons[FMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[3]].pt,muons[FMuons_id[3]].eta,muons[FMuons_id[3]].phi,muons[FMuons_id[3]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,muons[FMuons_id[2]].charge,muons[FMuons_id[3]].charge,FMuons_id[0],FMuons_id[1],FMuons_id[2],FMuons_id[3])

        SR4L_z1l1_pt=z1l1_v4.Pt()
        SR4L_z1l1_eta=z1l1_v4.Eta()
        SR4L_z1l1_phi=z1l1_v4.Phi()
        SR4L_z1l1_mass=z1l1_v4.M()
        SR4L_z1l2_pt=z1l2_v4.Pt()
        SR4L_z1l2_eta=z1l2_v4.Eta()
        SR4L_z1l2_phi=z1l2_v4.Phi()
        SR4L_z1l2_mass=z1l2_v4.M()
        SR4L_z2l1_pt=z2l1_v4.Pt()
        SR4L_z2l1_eta=z2l1_v4.Eta()
        SR4L_z2l1_phi=z2l1_v4.Phi()
        SR4L_z2l1_mass=z2l1_v4.M()
        SR4L_z2l2_pt=z2l2_v4.Pt()
        SR4L_z2l2_eta=z2l2_v4.Eta()
        SR4L_z2l2_phi=z2l2_v4.Phi()
        SR4L_z2l2_mass=z2l2_v4.M()

        SR4L_z1_pt=(z1l1_v4+z1l2_v4).Pt()
        SR4L_z1_eta=(z1l1_v4+z1l2_v4).Eta()
        SR4L_z1_phi=(z1l1_v4+z1l2_v4).Phi()
        SR4L_z1_mass=(z1l1_v4+z1l2_v4).M()
        SR4L_z2_pt=(z2l1_v4+z2l2_v4).Pt()
        SR4L_z2_eta=(z2l1_v4+z2l2_v4).Eta()
        SR4L_z2_phi=(z2l1_v4+z2l2_v4).Phi()
        SR4L_z2_mass=(z2l1_v4+z2l2_v4).M()

        z1_v4 = z1l1_v4+z1l2_v4
        z2_v4 = z2l1_v4+z2l2_v4

      # in two muon case, the z1 is always decaying to 2 muons
      elif len(TFMuons_id)==2:
        SR_4L_FlavorBit=48
        if len(TMuons_id)==2:
          if len(TElectrons_id)==2:
            SR_4L_FakeBit=15
            chargeTot_tmp1=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge
            chargeTot_tmp2=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=TMuons_id[0]
            SR4L_z1l2_id=TMuons_id[1]
            SR4L_z2l1_id=TElectrons_id[0]
            SR4L_z2l2_id=TElectrons_id[1]
            l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          elif len(TElectrons_id)==1:
            SR_4L_FakeBit=7
            chargeTot_tmp1=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge
            chargeTot_tmp2=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=TMuons_id[0]
            SR4L_z1l2_id=TMuons_id[1]
            SR4L_z2l1_id=TElectrons_id[0]
            SR4L_z2l2_id=FElectrons_id[0]
            l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          else:
            SR_4L_FakeBit=3
            chargeTot_tmp1=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge
            chargeTot_tmp2=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=TMuons_id[0]
            SR4L_z1l2_id=TMuons_id[1]
            SR4L_z2l1_id=FElectrons_id[0]
            SR4L_z2l2_id=FElectrons_id[1]
            l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
        elif len(TMuons_id)==1:
          if len(TElectrons_id)==2:
            SR_4L_FakeBit=13
            chargeTot_tmp1=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge
            chargeTot_tmp2=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=TMuons_id[0]
            SR4L_z1l2_id=FMuons_id[0]
            SR4L_z2l1_id=TElectrons_id[0]
            SR4L_z2l2_id=TElectrons_id[1]
            l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          elif len(TElectrons_id)==1:
            SR_4L_FakeBit=5
            chargeTot_tmp1=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge
            chargeTot_tmp2=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=TMuons_id[0]
            SR4L_z1l2_id=FMuons_id[0]
            SR4L_z2l1_id=TElectrons_id[0]
            SR4L_z2l2_id=FElectrons_id[0]
            l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          else:
            SR_4L_FakeBit=1
            chargeTot_tmp1=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge
            chargeTot_tmp2=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=TMuons_id[0]
            SR4L_z1l2_id=FMuons_id[0]
            SR4L_z2l1_id=FElectrons_id[0]
            SR4L_z2l2_id=FElectrons_id[1]
            l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
        else:
          if len(TElectrons_id)==2:
            SR_4L_FakeBit=12
            chargeTot_tmp1=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
            chargeTot_tmp2=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=FMuons_id[0]
            SR4L_z1l2_id=FMuons_id[1]
            SR4L_z2l1_id=TElectrons_id[0]
            SR4L_z2l2_id=TElectrons_id[1]
            l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          elif len(TElectrons_id)==1:
            SR_4L_FakeBit=4
            chargeTot_tmp1=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
            chargeTot_tmp2=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=FMuons_id[0]
            SR4L_z1l2_id=FMuons_id[1]
            SR4L_z2l1_id=TElectrons_id[0]
            SR4L_z2l2_id=FElectrons_id[0]
            l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          else:
            SR_4L_FakeBit=0
            chargeTot_tmp1=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
            chargeTot_tmp2=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp1)==0:return False
            if not abs(chargeTot_tmp2)==0:return False
            SR4L_z1l1_id=FMuons_id[0]
            SR4L_z1l2_id=FMuons_id[1]
            SR4L_z2l1_id=FElectrons_id[0]
            SR4L_z2l2_id=FElectrons_id[1]
            l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)

        SR4L_z1l1_pt=l1_v4.Pt()
        SR4L_z1l1_eta=l1_v4.Eta()
        SR4L_z1l1_phi=l1_v4.Phi()
        SR4L_z1l1_mass=l1_v4.M()
        SR4L_z1l2_pt=l2_v4.Pt()
        SR4L_z1l2_eta=l2_v4.Eta()
        SR4L_z1l2_phi=l2_v4.Phi()
        SR4L_z1l2_mass=l2_v4.M()
        SR4L_z2l1_pt=l3_v4.Pt()
        SR4L_z2l1_eta=l3_v4.Eta()
        SR4L_z2l1_phi=l3_v4.Phi()
        SR4L_z2l1_mass=l3_v4.M()
        SR4L_z2l2_pt=l4_v4.Pt()
        SR4L_z2l2_eta=l4_v4.Eta()
        SR4L_z2l2_phi=l4_v4.Phi()
        SR4L_z2l2_mass=l4_v4.M()

        SR4L_z1_pt=(l1_v4+l2_v4).Pt()
        SR4L_z1_eta=(l1_v4+l2_v4).Eta()
        SR4L_z1_phi=(l1_v4+l2_v4).Phi()
        SR4L_z1_mass=(l1_v4+l2_v4).M()
        SR4L_z2_pt=(l3_v4+l4_v4).Pt()
        SR4L_z2_eta=(l3_v4+l4_v4).Eta()
        SR4L_z2_phi=(l3_v4+l4_v4).Phi()
        SR4L_z2_mass=(l3_v4+l4_v4).M()
        z1_v4=l1_v4+l2_v4
        z2_v4=l3_v4+l4_v4

      elif len(TFMuons_id)==0:
        SR_4L_FlavorBit=44
        if len(TElectrons_id)==4:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[TElectrons_id[3]].charge
          if not abs(chargeTot_tmp)==0:return False
          SR_4L_FakeBit=15
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(eles[TElectrons_id[3]].pt,eles[TElectrons_id[3]].eta,eles[TElectrons_id[3]].phi,eles[TElectrons_id[3]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,eles[TElectrons_id[3]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2],TElectrons_id[3])

        elif len(TElectrons_id)==3:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[FElectrons_id[0]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,eles[FElectrons_id[0]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2],FElectrons_id[0])
          if SR4L_z1l1_id==FElectrons_id[0]:
            SR_4L_FakeBit=14
          elif SR4L_z1l2_id==FElectrons_id[0]:
            SR_4L_FakeBit=13
          elif SR4L_z2l1_id==FElectrons_id[0]:
            SR_4L_FakeBit=11
          else:
            SR_4L_FakeBit=7

        elif len(TElectrons_id)==2:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          l4_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,TElectrons_id[0],TElectrons_id[1],FElectrons_id[0],FElectrons_id[1])
          id_tmp=[SR4L_z2l2_id,SR4L_z2l1_id,SR4L_z1l2_id,SR4L_z1l1_id]
          real_index1=id_tmp.index(TElectrons_id[0])
          real_index2=id_tmp.index(TElectrons_id[1])
          SR_4L_FakeBit = (1 << (3 - real_index1)) + (1 << (3 - real_index2))

        elif len(TElectrons_id)==1:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge+eles[FElectrons_id[2]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          l3_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
          l4_v4.SetPtEtaPhiM(eles[FElectrons_id[2]].pt,eles[FElectrons_id[2]].eta,eles[FElectrons_id[2]].phi,eles[FElectrons_id[2]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,eles[FElectrons_id[2]].charge,TElectrons_id[0],FElectrons_id[0],FElectrons_id[1],FElectrons_id[2])
          if SR4L_z1l1_id==TElectrons_id[0]:
            SR_4L_FakeBit=1
          elif SR4L_z1l2_id==TElectrons_id[0]:
            SR_4L_FakeBit=2
          elif SR4L_z2l1_id==TElectrons_id[0]:
            SR_4L_FakeBit=4
          else:
            SR_4L_FakeBit=8
        else:
          chargeTot_tmp=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge+eles[FElectrons_id[2]].charge+eles[FElectrons_id[3]].charge
          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[FElectrons_id[2]].pt,eles[FElectrons_id[2]].eta,eles[FElectrons_id[2]].phi,eles[FElectrons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(eles[FElectrons_id[3]].pt,eles[FElectrons_id[3]].eta,eles[FElectrons_id[3]].phi,eles[FElectrons_id[3]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR4L_z1l1_id,SR4L_z1l2_id,SR4L_z2l1_id,SR4L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,eles[FElectrons_id[2]].charge,eles[FElectrons_id[3]].charge,FElectrons_id[0],FElectrons_id[1],FElectrons_id[2],FElectrons_id[3])

        SR4L_z1l1_pt=z1l1_v4.Pt()
        SR4L_z1l1_eta=z1l1_v4.Eta()
        SR4L_z1l1_phi=z1l1_v4.Phi()
        SR4L_z1l1_mass=z1l1_v4.M()
        SR4L_z1l2_pt=z1l2_v4.Pt()
        SR4L_z1l2_eta=z1l2_v4.Eta()
        SR4L_z1l2_phi=z1l2_v4.Phi()
        SR4L_z1l2_mass=z1l2_v4.M()
        SR4L_z2l1_pt=z2l1_v4.Pt()
        SR4L_z2l1_eta=z2l1_v4.Eta()
        SR4L_z2l1_phi=z2l1_v4.Phi()
        SR4L_z2l1_mass=z2l1_v4.M()
        SR4L_z2l2_pt=z2l2_v4.Pt()
        SR4L_z2l2_eta=z2l2_v4.Eta()
        SR4L_z2l2_phi=z2l2_v4.Phi()
        SR4L_z2l2_mass=z2l2_v4.M()

        SR4L_z1_pt=(z1l1_v4+z1l2_v4).Pt()
        SR4L_z1_eta=(z1l1_v4+z1l2_v4).Eta()
        SR4L_z1_phi=(z1l1_v4+z1l2_v4).Phi()
        SR4L_z1_mass=(z1l1_v4+z1l2_v4).M()
        SR4L_z2_pt=(z2l1_v4+z2l2_v4).Pt()
        SR4L_z2_eta=(z2l1_v4+z2l2_v4).Eta()
        SR4L_z2_phi=(z2l1_v4+z2l2_v4).Phi()
        SR4L_z2_mass=(z2l1_v4+z2l2_v4).M()
        z1_v4=z1l1_v4+z1l2_v4
        z2_v4=z2l1_v4+z2l2_v4

      w_j1, w_j2 = assign_2jets_from_many(TightJet_v4,MW)
      SR4L_wj1_id=TightJet_id[w_j1]
      SR4L_wj1_pt=TightJet_v4[w_j1].Pt()
      SR4L_wj1_eta=TightJet_v4[w_j1].Eta()
      SR4L_wj1_phi=TightJet_v4[w_j1].Phi()
      SR4L_wj1_mass=TightJet_v4[w_j1].M()
      SR4L_wj2_id=TightJet_id[w_j2]
      SR4L_wj2_pt=TightJet_v4[w_j2].Pt()
      SR4L_wj2_eta=TightJet_v4[w_j2].Eta()
      SR4L_wj2_phi=TightJet_v4[w_j2].Phi()
      SR4L_wj2_mass=TightJet_v4[w_j2].M()
      SR4L_wjj_pt=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).Pt()
      SR4L_wjj_eta=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).Eta()
      SR4L_wjj_phi=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).Phi()
      SR4L_wjj_mass=(TightJet_v4[w_j1]+TightJet_v4[w_j2]).M()
      wjj_v4=TightJet_v4[w_j1]+TightJet_v4[w_j2]

      SR4L_dR_z1_wjj = z1_v4.DeltaR(wjj_v4)
      SR4L_dEta_z1_wjj = abs(z1_v4.Eta()-wjj_v4.Eta())
      SR4L_dPhi_z1_wjj = z1_v4.DeltaPhi(wjj_v4)
      SR4L_z1wjj_pt = (z1_v4+wjj_v4).Pt()
      SR4L_z1wjj_eta = (z1_v4+wjj_v4).Eta()
      SR4L_z1wjj_phi = (z1_v4+wjj_v4).Phi()
      SR4L_z1wjj_mass = (z1_v4+wjj_v4).M()
      SR4L_dR_z1_z2 = z1_v4.DeltaR(z2_v4)
      SR4L_dEta_z1_z2 = abs(z1_v4.Eta()-z2_v4.Eta())
      SR4L_dPhi_z1_z2 = z1_v4.DeltaPhi(z2_v4)
      SR4L_z1z2_pt = (z1_v4+z2_v4).Pt()
      SR4L_z1z2_eta = (z1_v4+z2_v4).Eta()
      SR4L_z1z2_phi = (z1_v4+z2_v4).Phi()
      SR4L_z1z2_mass = (z1_v4+z2_v4).M()
      SR4L_dR_wjj_z2 = wjj_v4.DeltaR(z2_v4)
      SR4L_dEta_wjj_z2 = abs(wjj_v4.Eta()-z2_v4.Eta())
      SR4L_dPhi_wjj_z2 = wjj_v4.DeltaPhi(z2_v4)
      SR4L_wjjz2_pt = (wjj_v4+z2_v4).Pt()
      SR4L_wjjz2_eta = (wjj_v4+z2_v4).Eta()
      SR4L_wjjz2_phi = (wjj_v4+z2_v4).Phi()
      SR4L_wjjz2_mass = (wjj_v4+z2_v4).M()
      SR4L_wzz_pt = (z1_v4+wjj_v4+z2_v4).Pt()
      SR4L_wzz_eta = (z1_v4+wjj_v4+z2_v4).Eta()
      SR4L_wzz_phi = (z1_v4+wjj_v4+z2_v4).Phi()
      SR4L_wzz_mass = (z1_v4+wjj_v4+z2_v4).M()

    self.out.fillBranch("SR_4L_FlavorBit", SR_4L_FlavorBit)
    self.out.fillBranch("SR_4L_FakeBit", SR_4L_FakeBit)
    self.out.fillBranch("SR4L_z1l1_id", SR4L_z1l1_id)
    self.out.fillBranch("SR4L_z1l2_id", SR4L_z1l2_id)
    self.out.fillBranch("SR4L_z2l1_id", SR4L_z2l1_id)
    self.out.fillBranch("SR4L_z2l2_id", SR4L_z2l2_id)
    self.out.fillBranch("SR4L_z1l1_pt", SR4L_z1l1_pt)
    self.out.fillBranch("SR4L_z1l1_eta", SR4L_z1l1_eta)
    self.out.fillBranch("SR4L_z1l1_phi", SR4L_z1l1_phi)
    self.out.fillBranch("SR4L_z1l1_mass", SR4L_z1l1_mass)
    self.out.fillBranch("SR4L_z1l2_pt", SR4L_z1l2_pt)
    self.out.fillBranch("SR4L_z1l2_eta", SR4L_z1l2_eta)
    self.out.fillBranch("SR4L_z1l2_phi", SR4L_z1l2_phi)
    self.out.fillBranch("SR4L_z1l2_mass", SR4L_z1l2_mass)
    self.out.fillBranch("SR4L_z2l1_pt", SR4L_z2l1_pt)
    self.out.fillBranch("SR4L_z2l1_eta", SR4L_z2l1_eta)
    self.out.fillBranch("SR4L_z2l1_phi", SR4L_z2l1_phi)
    self.out.fillBranch("SR4L_z2l1_mass", SR4L_z2l1_mass)
    self.out.fillBranch("SR4L_z2l2_pt", SR4L_z2l2_pt)
    self.out.fillBranch("SR4L_z2l2_eta", SR4L_z2l2_eta)
    self.out.fillBranch("SR4L_z2l2_phi", SR4L_z2l2_phi)
    self.out.fillBranch("SR4L_z2l2_mass", SR4L_z2l2_mass)
    self.out.fillBranch("SR4L_z1_pt", SR4L_z1_pt)
    self.out.fillBranch("SR4L_z1_eta", SR4L_z1_eta)
    self.out.fillBranch("SR4L_z1_phi", SR4L_z1_phi)
    self.out.fillBranch("SR4L_z1_mass", SR4L_z1_mass)
    self.out.fillBranch("SR4L_z2_pt", SR4L_z2_pt)
    self.out.fillBranch("SR4L_z2_eta", SR4L_z2_eta)
    self.out.fillBranch("SR4L_z2_phi", SR4L_z2_phi)
    self.out.fillBranch("SR4L_z2_mass", SR4L_z2_mass)
    self.out.fillBranch("SR4L_wj1_id", SR4L_wj1_id)
    self.out.fillBranch("SR4L_wj1_pt", SR4L_wj1_pt)
    self.out.fillBranch("SR4L_wj1_eta", SR4L_wj1_eta)
    self.out.fillBranch("SR4L_wj1_phi", SR4L_wj1_phi)
    self.out.fillBranch("SR4L_wj1_mass", SR4L_wj1_mass)
    self.out.fillBranch("SR4L_wj2_id", SR4L_wj2_id)
    self.out.fillBranch("SR4L_wj2_pt", SR4L_wj2_pt)
    self.out.fillBranch("SR4L_wj2_eta", SR4L_wj2_eta)
    self.out.fillBranch("SR4L_wj2_phi", SR4L_wj2_phi)
    self.out.fillBranch("SR4L_wj2_mass", SR4L_wj2_mass)
    self.out.fillBranch("SR4L_wjj_pt", SR4L_wjj_pt)
    self.out.fillBranch("SR4L_wjj_eta", SR4L_wjj_eta)
    self.out.fillBranch("SR4L_wjj_phi", SR4L_wjj_phi)
    self.out.fillBranch("SR4L_wjj_mass", SR4L_wjj_mass)
    self.out.fillBranch("SR4L_dR_z1_wjj", SR4L_dR_z1_wjj)
    self.out.fillBranch("SR4L_dEta_z1_wjj", SR4L_dEta_z1_wjj)
    self.out.fillBranch("SR4L_dPhi_z1_wjj", SR4L_dPhi_z1_wjj)
    self.out.fillBranch("SR4L_z1wjj_pt", SR4L_z1wjj_pt)
    self.out.fillBranch("SR4L_z1wjj_eta", SR4L_z1wjj_eta)
    self.out.fillBranch("SR4L_z1wjj_phi", SR4L_z1wjj_phi)
    self.out.fillBranch("SR4L_z1wjj_mass", SR4L_z1wjj_mass)
    self.out.fillBranch("SR4L_dR_z1_z2", SR4L_dR_z1_z2)
    self.out.fillBranch("SR4L_dEta_z1_z2", SR4L_dEta_z1_z2)
    self.out.fillBranch("SR4L_dPhi_z1_z2", SR4L_dPhi_z1_z2)
    self.out.fillBranch("SR4L_z1z2_pt", SR4L_z1z2_pt)
    self.out.fillBranch("SR4L_z1z2_eta", SR4L_z1z2_eta)
    self.out.fillBranch("SR4L_z1z2_phi", SR4L_z1z2_phi)
    self.out.fillBranch("SR4L_z1z2_mass", SR4L_z1z2_mass)
    self.out.fillBranch("SR4L_dR_wjj_z2", SR4L_dR_wjj_z2)
    self.out.fillBranch("SR4L_dEta_wjj_z2", SR4L_dEta_wjj_z2)
    self.out.fillBranch("SR4L_dPhi_wjj_z2", SR4L_dPhi_wjj_z2)
    self.out.fillBranch("SR4L_wjjz2_eta", SR4L_wjjz2_eta)
    self.out.fillBranch("SR4L_wjjz2_phi", SR4L_wjjz2_phi)
    self.out.fillBranch("SR4L_wjjz2_mass", SR4L_wjjz2_mass)
    self.out.fillBranch("SR4L_wzz_pt", SR4L_wzz_pt)
    self.out.fillBranch("SR4L_wzz_eta", SR4L_wzz_eta)
    self.out.fillBranch("SR4L_wzz_phi", SR4L_wzz_phi)
    self.out.fillBranch("SR4L_wzz_mass", SR4L_wzz_mass)


    # l1 is z1l1, l2 is z1l2, l3 is z2l1, l4 is z2l2, l5 is wl
    # the bit is l5 | l4 | l3 | l2 | l1 |
    # if there are Z->mm and Z->ee, the Z->mm is always Z1
    SR_5L_FlavorBit=-1
    SR_5L_FakeBit=-1
    SR5L_z1l1_id=-99
    SR5L_z1l2_id=-99
    SR5L_z2l1_id=-99
    SR5L_z2l2_id=-99
    SR5L_wl_id=-99
    SR5L_z1l1_pt=-99
    SR5L_z1l1_eta=-99
    SR5L_z1l1_phi=-99
    SR5L_z1l1_mass=-99
    SR5L_z1l2_pt=-99
    SR5L_z1l2_eta=-99
    SR5L_z1l2_phi=-99
    SR5L_z1l2_mass=-99
    SR5L_z2l1_pt=-99
    SR5L_z2l1_eta=-99
    SR5L_z2l1_phi=-99
    SR5L_z2l1_mass=-99
    SR5L_z2l2_pt=-99
    SR5L_z2l2_eta=-99
    SR5L_z2l2_phi=-99
    SR5L_z2l2_mass=-99
    SR5L_wl_pt=-99
    SR5L_wl_eta=-99
    SR5L_wl_phi=-99
    SR5L_wl_mass=-99
    SR5L_z1_pt=-99
    SR5L_z1_eta=-99
    SR5L_z1_phi=-99
    SR5L_z1_mass=-99
    SR5L_z2_pt=-99
    SR5L_z2_eta=-99
    SR5L_z2_phi=-99
    SR5L_z2_mass=-99
    SR5L_w_pt=-99
    SR5L_w_eta=-99
    SR5L_w_phi=-99
    SR5L_w_mass=-99
    SR5L_dR_z1_w = -99
    SR5L_dEta_z1_w =-99
    SR5L_dPhi_z1_w =-99
    SR5L_z1w_pt = -99
    SR5L_z1w_eta =-99
    SR5L_z1w_phi =-99
    SR5L_z1w_mass = -99
    SR5L_dR_z1_z2 = -99
    SR5L_dEta_z1_z2 = -99
    SR5L_dPhi_z1_z2 = -99
    SR5L_z1z2_pt = -99
    SR5L_z1z2_eta =-99
    SR5L_z1z2_phi =-99
    SR5L_z1z2_mass =-99
    SR5L_dR_w_z2 =-99
    SR5L_dEta_w_z2 =-99
    SR5L_dPhi_w_z2 =-99
    SR5L_wz2_pt = -99
    SR5L_wz2_eta =-99
    SR5L_wz2_phi = -99
    SR5L_wz2_mass =-99
    SR5L_wzz_pt = -99
    SR5L_wzz_eta =-99
    SR5L_wzz_phi =-99
    SR5L_wzz_mass = -99

    if (len(TFMuons_id)+len(TFElectrons_id))==5:
      l1_v4=TLorentzVector()
      l2_v4=TLorentzVector()
      l3_v4=TLorentzVector()
      l4_v4=TLorentzVector()
      l5_v4=TLorentzVector()
      z1_v4=TLorentzVector()
      z2_v4=TLorentzVector()
      wlv_v4=TLorentzVector()
      if len(TFMuons_id)==5:
        SR_5L_FlavorBit=65
        if len(TMuons_id)==5:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[TMuons_id[3]].charge+muons[TMuons_id[4]].charge

          if not abs(chargeTot_tmp)==1:return False
          SR_5L_FakeBit=31
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[TMuons_id[3]].pt,muons[TMuons_id[3]].eta,muons[TMuons_id[3]].phi,muons[TMuons_id[3]].mass)
          l5_v4.SetPtEtaPhiM(muons[TMuons_id[4]].pt,muons[TMuons_id[4]].eta,muons[TMuons_id[4]].phi,muons[TMuons_id[4]].mass)

          l1_all=Particle(p4=l1_v4, Id=TMuons_id[0], charge=muons[TMuons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TMuons_id[1], charge=muons[TMuons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=TMuons_id[2], charge=muons[TMuons_id[2]].charge)
          l4_all=Particle(p4=l4_v4, Id=TMuons_id[3], charge=muons[TMuons_id[3]].charge)
          l5_all=Particle(p4=l5_v4, Id=TMuons_id[4], charge=muons[TMuons_id[4]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']

        elif len(TMuons_id)==4:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[TMuons_id[3]].charge+muons[FMuons_id[0]].charge

          if not abs(chargeTot_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[TMuons_id[3]].pt,muons[TMuons_id[3]].eta,muons[TMuons_id[3]].phi,muons[TMuons_id[3]].mass)
          l5_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)

          l1_all=Particle(p4=l1_v4, Id=TMuons_id[0], charge=muons[TMuons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TMuons_id[1], charge=muons[TMuons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=TMuons_id[2], charge=muons[TMuons_id[2]].charge)
          l4_all=Particle(p4=l4_v4, Id=TMuons_id[3], charge=muons[TMuons_id[3]].charge)
          l5_all=Particle(p4=l5_v4, Id=FMuons_id[0], charge=muons[FMuons_id[0]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']
          id_tmp=[wl.Id,z2l2.Id,z2l1.Id,z1l2.Id,z1l1.Id]
          fake_index1=id_tmp.index(FMuons_id[0])
          SR_5L_FakeBit=31 - (1<<(4-fake_index1))

        elif len(TMuons_id)==3:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge

          if not abs(chargeTot_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l5_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)

          l1_all=Particle(p4=l1_v4, Id=TMuons_id[0], charge=muons[TMuons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TMuons_id[1], charge=muons[TMuons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=TMuons_id[2], charge=muons[TMuons_id[2]].charge)
          l4_all=Particle(p4=l4_v4, Id=FMuons_id[0], charge=muons[FMuons_id[0]].charge)
          l5_all=Particle(p4=l5_v4, Id=FMuons_id[1], charge=muons[FMuons_id[1]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']
          id_tmp=[wl.Id,z2l2.Id,z2l1.Id,z1l2.Id,z1l1.Id]
          fake_index1=id_tmp.index(FMuons_id[0])
          fake_index2=id_tmp.index(FMuons_id[1])
          SR_5L_FakeBit=31 - (1<<(4-fake_index1)) - (1<<(4-fake_index2))

        elif len(TMuons_id)==2:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge+muons[FMuons_id[2]].charge

          if not abs(chargeTot_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          l5_v4.SetPtEtaPhiM(muons[FMuons_id[2]].pt,muons[FMuons_id[2]].eta,muons[FMuons_id[2]].phi,muons[FMuons_id[2]].mass)

          l1_all=Particle(p4=l1_v4, Id=TMuons_id[0], charge=muons[TMuons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TMuons_id[1], charge=muons[TMuons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=FMuons_id[0], charge=muons[FMuons_id[0]].charge)
          l4_all=Particle(p4=l4_v4, Id=FMuons_id[1], charge=muons[FMuons_id[1]].charge)
          l5_all=Particle(p4=l5_v4, Id=FMuons_id[2], charge=muons[FMuons_id[2]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']
          id_tmp=[wl.Id,z2l2.Id,z2l1.Id,z1l2.Id,z1l1.Id]
          fake_index1=id_tmp.index(FMuons_id[0])
          fake_index2=id_tmp.index(FMuons_id[1])
          fake_index3=id_tmp.index(FMuons_id[2])
          SR_5L_FakeBit=31 - (1<<(4-fake_index1)) - (1<<(4-fake_index2)) - (1<<(4-fake_index3))

        SR5L_z1l1_id=z1l1.Id
        SR5L_z1l2_id=z1l2.Id
        SR5L_z2l1_id=z2l1.Id
        SR5L_z2l2_id=z2l2.Id
        SR5L_wl_id=wl.Id
        SR5L_z1l1_pt=z1l1.p4.Pt()
        SR5L_z1l1_eta=z1l1.p4.Eta()
        SR5L_z1l1_phi=z1l1.p4.Phi()
        SR5L_z1l1_mass=z1l1.p4.M()
        SR5L_z1l2_pt=z1l2.p4.Pt()
        SR5L_z1l2_eta=z1l2.p4.Eta()
        SR5L_z1l2_phi=z1l2.p4.Phi()
        SR5L_z1l2_mass=z1l2.p4.M()
        SR5L_z2l1_pt=z2l1.p4.Pt()
        SR5L_z2l1_eta=z2l1.p4.Eta()
        SR5L_z2l1_phi=z2l1.p4.Phi()
        SR5L_z2l1_mass=z2l1.p4.M()
        SR5L_z2l2_pt=z2l2.p4.Pt()
        SR5L_z2l2_eta=z2l2.p4.Eta()
        SR5L_z2l2_phi=z2l2.p4.Phi()
        SR5L_z2l2_mass=z2l2.p4.M()
        SR5L_wl_pt=wl.p4.Pt()
        SR5L_wl_eta=wl.p4.Eta()
        SR5L_wl_phi=wl.p4.Phi()
        SR5L_wl_mass=wl.p4.M()
        
        z1_v4=z1l1.p4+z1l2.p4
        z2_v4=z2l1.p4+z2l2.p4
        #wlv_v4=w_v4(wl.p4, met_user, met_phi_user)
        wlv_v4=w_v4(wl.p4, met_CorrePhi, met_phi_CorrePhi)

        SR5L_z1_pt=z1_v4.Pt()
        SR5L_z1_eta=z1_v4.Eta()
        SR5L_z1_phi=z1_v4.Phi()
        SR5L_z1_mass=z1_v4.M()
        SR5L_z2_pt=z2_v4.Pt()
        SR5L_z2_eta=z2_v4.Eta()
        SR5L_z2_phi=z2_v4.Phi()
        SR5L_z2_mass=z2_v4.M()
        SR5L_w_pt=wlv_v4.Pt()
        SR5L_w_eta=wlv_v4.Eta()
        SR5L_w_phi=wlv_v4.Phi()
        SR5L_w_mass=wlv_v4.M()

      elif len(TFMuons_id)==4:
        SR_5L_FlavorBit=63
        if len(TMuons_id)==4:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[TMuons_id[3]].charge

          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[TMuons_id[3]].pt,muons[TMuons_id[3]].eta,muons[TMuons_id[3]].phi,muons[TMuons_id[3]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,muons[TMuons_id[3]].charge,TMuons_id[0],TMuons_id[1],TMuons_id[2],TMuons_id[3])
          if len(TElectrons_id)==1:
            SR_5L_FakeBit=31
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            SR5L_wl_id=TElectrons_id[0]
          else:
            SR_5L_FakeBit=15
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            SR5L_wl_id=FElectrons_id[0]

        elif len(TMuons_id)==3:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge+muons[FMuons_id[0]].charge

          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,muons[FMuons_id[0]].charge,TMuons_id[0],TMuons_id[1],TMuons_id[2],FMuons_id[0])
          id_tmp=[SR5L_z2l2_id,SR5L_z2l1_id,SR5L_z1l2_id,SR5L_z1l1_id]
          fake_index=id_tmp.index(FMuons_id[0])
          SR_5L_FakeBit_ = 15 - (1<<(3-fake_index))
          if len(TElectrons_id)==1:
            SR_5L_FakeBit=16+SR_5L_FakeBit_
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            SR5L_wl_id=TElectrons_id[0]
          else:
            SR_5L_FakeBit=SR_5L_FakeBit_
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            SR5L_wl_id=FElectrons_id[0]

        elif len(TMuons_id)==2:
          chargeTot_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge

          if not abs(chargeTot_tmp)==0:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l4_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,TMuons_id[0],TMuons_id[1],FMuons_id[0],FMuons_id[1])
          id_tmp=[SR5L_z2l2_id,SR5L_z2l1_id,SR5L_z1l2_id,SR5L_z1l1_id]
          fake_index1=id_tmp.index(FMuons_id[0])
          fake_index2=id_tmp.index(FMuons_id[1])
          SR_5L_FakeBit_ = 15 - (1<<(3-fake_index1)) - (1<<(3-fake_index2))
          if len(TElectrons_id)==1:
            SR_5L_FakeBit=16+SR_5L_FakeBit_
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            SR5L_wl_id=TElectrons_id[0]
          else:
            SR_5L_FakeBit=SR_5L_FakeBit_
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            SR5L_wl_id=FElectrons_id[0]

        SR5L_z1l1_pt=z1l1_v4.Pt()
        SR5L_z1l1_eta=z1l1_v4.Eta()
        SR5L_z1l1_phi=z1l1_v4.Phi()
        SR5L_z1l1_mass=z1l1_v4.M()
        SR5L_z1l2_pt=z1l2_v4.Pt()
        SR5L_z1l2_eta=z1l2_v4.Eta()
        SR5L_z1l2_phi=z1l2_v4.Phi()
        SR5L_z1l2_mass=z1l2_v4.M()
        SR5L_z2l1_pt=z2l1_v4.Pt()
        SR5L_z2l1_eta=z2l1_v4.Eta()
        SR5L_z2l1_phi=z2l1_v4.Phi()
        SR5L_z2l1_mass=z2l1_v4.M()
        SR5L_z2l2_pt=z2l2_v4.Pt()
        SR5L_z2l2_eta=z2l2_v4.Eta()
        SR5L_z2l2_phi=z2l2_v4.Phi()
        SR5L_z2l2_mass=z2l2_v4.M()
        SR5L_wl_pt=l5_v4.Pt()
        SR5L_wl_eta=l5_v4.Eta()
        SR5L_wl_phi=l5_v4.Phi()
        SR5L_wl_mass=l5_v4.M()
        
        z1_v4=z1l1_v4+z1l2_v4
        z2_v4=z2l1_v4+z2l2_v4
        #wlv_v4=w_v4(l5_v4, met_user, met_phi_user)
        wlv_v4=w_v4(l5_v4, met_CorrePhi, met_phi_CorrePhi)

        SR5L_z1_pt=z1_v4.Pt()
        SR5L_z1_eta=z1_v4.Eta()
        SR5L_z1_phi=z1_v4.Phi()
        SR5L_z1_mass=z1_v4.M()
        SR5L_z2_pt=z2_v4.Pt()
        SR5L_z2_eta=z2_v4.Eta()
        SR5L_z2_phi=z2_v4.Phi()
        SR5L_z2_mass=z2_v4.M()
        SR5L_w_pt=wlv_v4.Pt()
        SR5L_w_eta=wlv_v4.Eta()
        SR5L_w_phi=wlv_v4.Phi()
        SR5L_w_mass=wlv_v4.M()

      elif len(TFMuons_id)==3:
        SR_5L_FlavorBit=61
        if len(TMuons_id)==3:
          chargeTotmu_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[TMuons_id[2]].charge
          if not abs(chargeTotmu_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[TMuons_id[2]].pt,muons[TMuons_id[2]].eta,muons[TMuons_id[2]].phi,muons[TMuons_id[2]].mass)
          z1l1_v4,z1l2_v4,wl_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[TMuons_id[2]].charge,TMuons_id[0],TMuons_id[1],TMuons_id[2])

          if len(TElectrons_id)==2:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTotele_tmp)==0:return False
            SR_5L_FakeBit=31
            SR5L_z2l1_id=TElectrons_id[0]
            SR5L_z2l2_id=TElectrons_id[1]
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          elif len(TElectrons_id)==1:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTotele_tmp)==0:return False
            SR_5L_FakeBit=23
            SR5L_z2l1_id=TElectrons_id[0]
            SR5L_z2l2_id=FElectrons_id[0]
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          else:
            chargeTotele_tmp=eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTotele_tmp)==0:return False
            SR_5L_FakeBit=19
            SR5L_z2l1_id=TElectrons_id[0]
            SR5L_z2l2_id=FElectrons_id[0]
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)

        elif len(TMuons_id)==2:
          chargeTotmu_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge+muons[FMuons_id[0]].charge
          if not abs(chargeTotmu_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          z1l1_v4,z1l2_v4,wl_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[TMuons_id[0]].charge,muons[TMuons_id[1]].charge,muons[FMuons_id[0]].charge,TMuons_id[0],TMuons_id[1],FMuons_id[0])
          id_tmp=[SR5L_wl_id,9,10,SR5L_z1l2_id,SR5L_z1l1_id]
          fake_index=id_tmp.index(FMuons_id[0])
          SR_5L_FakeBit_=31 - (1<<(4-fake_index)) - 12

          if len(TElectrons_id)==2:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTotele_tmp)==0:return False
            SR_5L_FakeBit=SR_5L_FakeBit_+12
            SR5L_z2l1_id=TElectrons_id[0]
            SR5L_z2l2_id=TElectrons_id[1]
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          elif len(TElectrons_id)==1:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTotele_tmp)==0:return False
            SR_5L_FakeBit=SR_5L_FakeBit_+4
            SR5L_z2l1_id=TElectrons_id[0]
            SR5L_z2l2_id=FElectrons_id[0]
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)

        elif len(TMuons_id)==1:
          chargeTotmu_tmp=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
          if not abs(chargeTotmu_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l3_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          z1l1_v4,z1l2_v4,wl_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_wl_id=assign_3lep(l1_v4,l2_v4,l3_v4,muons[TMuons_id[0]].charge,muons[FMuons_id[0]].charge,muons[FMuons_id[1]].charge,TMuons_id[0],FMuons_id[0],FMuons_id[1])
          id_tmp=[SR5L_wl_id,9,10,SR5L_z1l2_id,SR5L_z1l1_id]
          real_index=id_tmp.index(TMuons_id[0])
          SR_5L_FakeBit_= 1<<(4-real_index)

          if len(TElectrons_id)==2:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge
            if not abs(chargeTotele_tmp)==0:return False
            SR_5L_FakeBit=SR_5L_FakeBit_+12
            SR5L_z2l1_id=TElectrons_id[0]
            SR5L_z2l2_id=FElectrons_id[0]
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)

        SR5L_z1l1_pt=z1l1_v4.Pt()
        SR5L_z1l1_eta=z1l1_v4.Eta()
        SR5L_z1l1_phi=z1l1_v4.Phi()
        SR5L_z1l1_mass=z1l1_v4.M()
        SR5L_z1l2_pt=z1l2_v4.Pt()
        SR5L_z1l2_eta=z1l2_v4.Eta()
        SR5L_z1l2_phi=z1l2_v4.Phi()
        SR5L_z1l2_mass=z1l2_v4.M()
        SR5L_z2l1_pt=l4_v4.Pt()
        SR5L_z2l1_eta=l4_v4.Eta()
        SR5L_z2l1_phi=l4_v4.Phi()
        SR5L_z2l1_mass=l4_v4.M()
        SR5L_z2l2_pt=l5_v4.Pt()
        SR5L_z2l2_eta=l5_v4.Eta()
        SR5L_z2l2_phi=l5_v4.Phi()
        SR5L_z2l2_mass=l5_v4.M()
        SR5L_wl_pt=wl_v4.Pt()
        SR5L_wl_eta=wl_v4.Eta()
        SR5L_wl_phi=wl_v4.Phi()
        SR5L_wl_mass=wl_v4.M()
        
        z1_v4=z1l1_v4+z1l2_v4
        z2_v4=l4_v4+l5_v4
        #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
        wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)

        SR5L_z1_pt=z1_v4.Pt()
        SR5L_z1_eta=z1_v4.Eta()
        SR5L_z1_phi=z1_v4.Phi()
        SR5L_z1_mass=z1_v4.M()
        SR5L_z2_pt=z2_v4.Pt()
        SR5L_z2_eta=z2_v4.Eta()
        SR5L_z2_phi=z2_v4.Phi()
        SR5L_z2_mass=z2_v4.M()
        SR5L_w_pt=wlv_v4.Pt()
        SR5L_w_eta=wlv_v4.Eta()
        SR5L_w_phi=wlv_v4.Phi()
        SR5L_w_mass=wlv_v4.M()

      elif len(TFMuons_id)==2:
        SR_5L_FlavorBit=59
        if len(TMuons_id)==2:
          chargeTotmu_tmp=muons[TMuons_id[0]].charge+muons[TMuons_id[1]].charge
          if not abs(chargeTotmu_tmp)==0:return False
          SR5L_z1l1_id=TMuons_id[0]
          SR5L_z1l2_id=TMuons_id[1]
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[TMuons_id[1]].pt,muons[TMuons_id[1]].eta,muons[TMuons_id[1]].phi,muons[TMuons_id[1]].mass)
          if len(TElectrons_id)==3:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge
            if not abs(chargeTotele_tmp)==1:return False
            SR_5L_FakeBit=31
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            z2l1_v4,z2l2_v4,wl_v4,SR5L_z2l1_id,SR5L_z2l2_id,SR5L_wl_id=assign_3lep(l3_v4,l4_v4,l5_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2])
          elif len(TElectrons_id)==2:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTotele_tmp)==1:return False
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            z2l1_v4,z2l2_v4,wl_v4,SR5L_z2l1_id,SR5L_z2l2_id,SR5L_wl_id=assign_3lep(l3_v4,l4_v4,l5_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[FElectrons_id[0]].charge,TElectrons_id[0],TElectrons_id[1],FElectrons_id[0])
            id_tmp=[SR5L_wl_id,SR5L_z2l2_id,SR5L_z2l1_id,9,10]
            fake_index=id_tmp.index(FElectrons_id[0])
            SR_5L_FakeBit=31 - (1<<(4-fake_index))
          elif len(TElectrons_id)==1:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTotele_tmp)==1:return False
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
            z2l1_v4,z2l2_v4,wl_v4,SR5L_z2l1_id,SR5L_z2l2_id,SR5L_wl_id=assign_3lep(l3_v4,l4_v4,l5_v4,eles[TElectrons_id[0]].charge,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,TElectrons_id[0],FElectrons_id[0],FElectrons_id[1])
            id_tmp=[SR5L_wl_id,SR5L_z2l2_id,SR5L_z2l1_id,9,10]
            real_index=id_tmp.index(TElectrons_id[0])
            SR_5L_FakeBit=3 + (1<<(4-real_index))

        elif len(TMuons_id)==1:
          chargeTotmu_tmp=muons[TMuons_id[0]].charge+muons[FMuons_id[0]].charge
          if not abs(chargeTotmu_tmp)==0:return False
          SR5L_z1l1_id=TMuons_id[0]
          SR5L_z1l2_id=FMuons_id[0]
          l1_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          if len(TElectrons_id)==3:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge
            if not abs(chargeTotele_tmp)==1:return False
            SR_5L_FakeBit=29
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            z2l1_v4,z2l2_v4,wl_v4,SR5L_z2l1_id,SR5L_z2l2_id,SR5L_wl_id=assign_3lep(l3_v4,l4_v4,l5_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2])
          elif len(TElectrons_id)==2:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTotele_tmp)==1:return False
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            z2l1_v4,z2l2_v4,wl_v4,SR5L_z2l1_id,SR5L_z2l2_id,SR5L_wl_id=assign_3lep(l3_v4,l4_v4,l5_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[FElectrons_id[0]].charge,TElectrons_id[0],TElectrons_id[1],FElectrons_id[0])
            id_tmp=[SR5L_wl_id,SR5L_z2l2_id,SR5L_z2l1_id,9,10]
            fake_index=id_tmp.index(FElectrons_id[0])
            SR_5L_FakeBit=29 - (1<<(4-fake_index))

        else:
          chargeTotmu_tmp=muons[FMuons_id[0]].charge+muons[FMuons_id[1]].charge
          if not abs(chargeTotmu_tmp)==0:return False
          SR5L_z1l1_id=FMuons_id[0]
          SR5L_z1l2_id=FMuons_id[1]
          l1_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(muons[FMuons_id[1]].pt,muons[FMuons_id[1]].eta,muons[FMuons_id[1]].phi,muons[FMuons_id[1]].mass)
          if len(TElectrons_id)==3:
            chargeTotele_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge
            if not abs(chargeTotele_tmp)==1:return False
            SR_5L_FakeBit=28
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l5_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            z2l1_v4,z2l2_v4,wl_v4,SR5L_z2l1_id,SR5L_z2l2_id,SR5L_wl_id=assign_3lep(l3_v4,l4_v4,l5_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2])

        SR5L_z1l1_pt=l1_v4.Pt()
        SR5L_z1l1_eta=l1_v4.Eta()
        SR5L_z1l1_phi=l1_v4.Phi()
        SR5L_z1l1_mass=l1_v4.M()
        SR5L_z1l2_pt=l2_v4.Pt()
        SR5L_z1l2_eta=l2_v4.Eta()
        SR5L_z1l2_phi=l2_v4.Phi()
        SR5L_z1l2_mass=l2_v4.M()
        SR5L_z2l1_pt=z2l1_v4.Pt()
        SR5L_z2l1_eta=z2l1_v4.Eta()
        SR5L_z2l1_phi=z2l1_v4.Phi()
        SR5L_z2l1_mass=z2l1_v4.M()
        SR5L_z2l2_pt=z2l2_v4.Pt()
        SR5L_z2l2_eta=z2l2_v4.Eta()
        SR5L_z2l2_phi=z2l2_v4.Phi()
        SR5L_z2l2_mass=z2l2_v4.M()
        SR5L_wl_pt=wl_v4.Pt()
        SR5L_wl_eta=wl_v4.Eta()
        SR5L_wl_phi=wl_v4.Phi()
        SR5L_wl_mass=wl_v4.M()
        
        z1_v4=l1_v4+l2_v4
        z2_v4=z2l1_v4+z2l2_v4
        #wlv_v4=w_v4(wl_v4, met_user, met_phi_user)
        wlv_v4=w_v4(wl_v4, met_CorrePhi, met_phi_CorrePhi)

        SR5L_z1_pt=z1_v4.Pt()
        SR5L_z1_eta=z1_v4.Eta()
        SR5L_z1_phi=z1_v4.Phi()
        SR5L_z1_mass=z1_v4.M()
        SR5L_z2_pt=z2_v4.Pt()
        SR5L_z2_eta=z2_v4.Eta()
        SR5L_z2_phi=z2_v4.Phi()
        SR5L_z2_mass=z2_v4.M()
        SR5L_w_pt=wlv_v4.Pt()
        SR5L_w_eta=wlv_v4.Eta()
        SR5L_w_phi=wlv_v4.Phi()
        SR5L_w_mass=wlv_v4.M()

      elif len(TFMuons_id)==1:
        SR_5L_FlavorBit=57
        if len(TMuons_id)==1:
          SR5L_wl_id=TMuons_id[0]
          l5_v4.SetPtEtaPhiM(muons[TMuons_id[0]].pt,muons[TMuons_id[0]].eta,muons[TMuons_id[0]].phi,muons[TMuons_id[0]].mass)
          if len(TElectrons_id)==4:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[TElectrons_id[3]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_5L_FakeBit=31
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[3]].pt,eles[TElectrons_id[3]].eta,eles[TElectrons_id[3]].phi,eles[TElectrons_id[3]].mass)
            z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,eles[TElectrons_id[3]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2],TElectrons_id[3])
          elif len(TElectrons_id)==3:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp)==0:return False
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,eles[FElectrons_id[0]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2],FElectrons_id[0])
            id_tmp=[SR5L_wl_id,SR5L_z2l2_id,SR5L_z2l1_id,SR5L_z1l2_id,SR5L_z1l1_id]
            fake_index=id_tmp.index(FElectrons_id[0])
            SR_5L_FakeBit=31 - (1<<(4-fake_index))
          elif len(TElectrons_id)==2:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge
            if not abs(chargeTot_tmp)==0:return False
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
            z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[FElectrons_id[0]].charge,eles[FElectrons_id[1]].charge,TElectrons_id[0],TElectrons_id[1],FElectrons_id[0],FElectrons_id[1])
            id_tmp=[SR5L_wl_id,SR5L_z2l2_id,SR5L_z2l1_id,SR5L_z1l2_id,SR5L_z1l1_id]
            fake_index1=id_tmp.index(FElectrons_id[0])
            fake_index2=id_tmp.index(FElectrons_id[1])
            SR_5L_FakeBit=31 - (1<<(4-fake_index1)) - (1<<(4-fake_index2))
        else:
          SR5L_wl_id=FMuons_id[0]
          l5_v4.SetPtEtaPhiM(muons[FMuons_id[0]].pt,muons[FMuons_id[0]].eta,muons[FMuons_id[0]].phi,muons[FMuons_id[0]].mass)
          if len(TElectrons_id)==4:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[TElectrons_id[3]].charge
            if not abs(chargeTot_tmp)==0:return False
            SR_5L_FakeBit=15
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            l4_v4.SetPtEtaPhiM(eles[TElectrons_id[3]].pt,eles[TElectrons_id[3]].eta,eles[TElectrons_id[3]].phi,eles[TElectrons_id[3]].mass)
            z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,eles[TElectrons_id[3]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2],TElectrons_id[3])
          elif len(TElectrons_id)==3:
            chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[FElectrons_id[0]].charge
            if not abs(chargeTot_tmp)==0:return False
            l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
            l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
            l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
            l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
            z1l1_v4,z1l2_v4,z2l1_v4,z2l2_v4,SR5L_z1l1_id,SR5L_z1l2_id,SR5L_z2l1_id,SR5L_z2l2_id=assign_4lep(l1_v4,l2_v4,l3_v4,l4_v4,eles[TElectrons_id[0]].charge,eles[TElectrons_id[1]].charge,eles[TElectrons_id[2]].charge,eles[FElectrons_id[0]].charge,TElectrons_id[0],TElectrons_id[1],TElectrons_id[2],FElectrons_id[0])
            id_tmp=[SR5L_wl_id,SR5L_z2l2_id,SR5L_z2l1_id,SR5L_z1l2_id,SR5L_z1l1_id]
            fake_index=id_tmp.index(FElectrons_id[0])
            SR_5L_FakeBit=31 - (1<<(4-fake_index))

        SR5L_z1l1_pt=z1l1_v4.Pt()
        SR5L_z1l1_eta=z1l1_v4.Eta()
        SR5L_z1l1_phi=z1l1_v4.Phi()
        SR5L_z1l1_mass=z1l1_v4.M()
        SR5L_z1l2_pt=z1l2_v4.Pt()
        SR5L_z1l2_eta=z1l2_v4.Eta()
        SR5L_z1l2_phi=z1l2_v4.Phi()
        SR5L_z1l2_mass=z1l2_v4.M()
        SR5L_z2l1_pt=z2l1_v4.Pt()
        SR5L_z2l1_eta=z2l1_v4.Eta()
        SR5L_z2l1_phi=z2l1_v4.Phi()
        SR5L_z2l1_mass=z2l1_v4.M()
        SR5L_z2l2_pt=z2l2_v4.Pt()
        SR5L_z2l2_eta=z2l2_v4.Eta()
        SR5L_z2l2_phi=z2l2_v4.Phi()
        SR5L_z2l2_mass=z2l2_v4.M()
        SR5L_wl_pt=l5_v4.Pt()
        SR5L_wl_eta=l5_v4.Eta()
        SR5L_wl_phi=l5_v4.Phi()
        SR5L_wl_mass=l5_v4.M()
        
        z1_v4=z1l1_v4+z1l2_v4
        z2_v4=z2l1_v4+z2l2_v4
        #wlv_v4=w_v4(l5_v4, met_user, met_phi_user)
        wlv_v4=w_v4(l5_v4, met_CorrePhi, met_phi_CorrePhi)

        SR5L_z1_pt=z1_v4.Pt()
        SR5L_z1_eta=z1_v4.Eta()
        SR5L_z1_phi=z1_v4.Phi()
        SR5L_z1_mass=z1_v4.M()
        SR5L_z2_pt=z2_v4.Pt()
        SR5L_z2_eta=z2_v4.Eta()
        SR5L_z2_phi=z2_v4.Phi()
        SR5L_z2_mass=z2_v4.M()
        SR5L_w_pt=wlv_v4.Pt()
        SR5L_w_eta=wlv_v4.Eta()
        SR5L_w_phi=wlv_v4.Phi()
        SR5L_w_mass=wlv_v4.M()

      else:
        SR_5L_FlavorBit=55
        if len(TElectrons_id)==5:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[TElectrons_id[3]].charge+eles[TElectrons_id[4]].charge

          if not abs(chargeTot_tmp)==1:return False
          SR_5L_FakeBit=31
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(eles[TElectrons_id[3]].pt,eles[TElectrons_id[3]].eta,eles[TElectrons_id[3]].phi,eles[TElectrons_id[3]].mass)
          l5_v4.SetPtEtaPhiM(eles[TElectrons_id[4]].pt,eles[TElectrons_id[4]].eta,eles[TElectrons_id[4]].phi,eles[TElectrons_id[4]].mass)

          l1_all=Particle(p4=l1_v4, Id=TElectrons_id[0], charge=eles[TElectrons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TElectrons_id[1], charge=eles[TElectrons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=TElectrons_id[2], charge=eles[TElectrons_id[2]].charge)
          l4_all=Particle(p4=l4_v4, Id=TElectrons_id[3], charge=eles[TElectrons_id[3]].charge)
          l5_all=Particle(p4=l5_v4, Id=TElectrons_id[4], charge=eles[TElectrons_id[4]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']

        elif len(TElectrons_id)==4:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[TElectrons_id[3]].charge+eles[FElectrons_id[0]].charge

          if not abs(chargeTot_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(eles[TElectrons_id[3]].pt,eles[TElectrons_id[3]].eta,eles[TElectrons_id[3]].phi,eles[TElectrons_id[3]].mass)
          l5_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)

          l1_all=Particle(p4=l1_v4, Id=TElectrons_id[0], charge=eles[TElectrons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TElectrons_id[1], charge=eles[TElectrons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=TElectrons_id[2], charge=eles[TElectrons_id[2]].charge)
          l4_all=Particle(p4=l4_v4, Id=TElectrons_id[3], charge=eles[TElectrons_id[3]].charge)
          l5_all=Particle(p4=l5_v4, Id=FElectrons_id[0], charge=eles[FElectrons_id[0]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']
          id_tmp=[wl.Id,z2l2.Id,z2l1.Id,z1l2.Id,z1l1.Id]
          fake_index1=id_tmp.index(FElectrons_id[0])
          SR_5L_FakeBit=31 - (1<<(4-fake_index1))

        elif len(TElectrons_id)==3:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[TElectrons_id[2]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge

          if not abs(chargeTot_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[TElectrons_id[2]].pt,eles[TElectrons_id[2]].eta,eles[TElectrons_id[2]].phi,eles[TElectrons_id[2]].mass)
          l4_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          l5_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)

          l1_all=Particle(p4=l1_v4, Id=TElectrons_id[0], charge=eles[TElectrons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TElectrons_id[1], charge=eles[TElectrons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=TElectrons_id[2], charge=eles[TElectrons_id[2]].charge)
          l4_all=Particle(p4=l4_v4, Id=FElectrons_id[0], charge=eles[FElectrons_id[0]].charge)
          l5_all=Particle(p4=l5_v4, Id=FElectrons_id[1], charge=eles[FElectrons_id[1]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']
          id_tmp=[wl.Id,z2l2.Id,z2l1.Id,z1l2.Id,z1l1.Id]
          fake_index1=id_tmp.index(FElectrons_id[0])
          fake_index2=id_tmp.index(FElectrons_id[1])
          SR_5L_FakeBit=31 - (1<<(4-fake_index1)) - (1<<(4-fake_index2))

        elif len(TElectrons_id)==2:
          chargeTot_tmp=eles[TElectrons_id[0]].charge+eles[TElectrons_id[1]].charge+eles[FElectrons_id[0]].charge+eles[FElectrons_id[1]].charge+eles[FElectrons_id[2]].charge

          if not abs(chargeTot_tmp)==1:return False
          l1_v4.SetPtEtaPhiM(eles[TElectrons_id[0]].pt,eles[TElectrons_id[0]].eta,eles[TElectrons_id[0]].phi,eles[TElectrons_id[0]].mass)
          l2_v4.SetPtEtaPhiM(eles[TElectrons_id[1]].pt,eles[TElectrons_id[1]].eta,eles[TElectrons_id[1]].phi,eles[TElectrons_id[1]].mass)
          l3_v4.SetPtEtaPhiM(eles[FElectrons_id[0]].pt,eles[FElectrons_id[0]].eta,eles[FElectrons_id[0]].phi,eles[FElectrons_id[0]].mass)
          l4_v4.SetPtEtaPhiM(eles[FElectrons_id[1]].pt,eles[FElectrons_id[1]].eta,eles[FElectrons_id[1]].phi,eles[FElectrons_id[1]].mass)
          l5_v4.SetPtEtaPhiM(eles[FElectrons_id[2]].pt,eles[FElectrons_id[2]].eta,eles[FElectrons_id[2]].phi,eles[FElectrons_id[2]].mass)

          l1_all=Particle(p4=l1_v4, Id=TElectrons_id[0], charge=eles[TElectrons_id[0]].charge)
          l2_all=Particle(p4=l2_v4, Id=TElectrons_id[1], charge=eles[TElectrons_id[1]].charge)
          l3_all=Particle(p4=l3_v4, Id=FElectrons_id[0], charge=eles[FElectrons_id[0]].charge)
          l4_all=Particle(p4=l4_v4, Id=FElectrons_id[1], charge=eles[FElectrons_id[1]].charge)
          l5_all=Particle(p4=l5_v4, Id=FElectrons_id[2], charge=eles[FElectrons_id[2]].charge)
          best_comb=assign_5lep(l1_all,l2_all,l3_all,l4_all,l5_all)
          z1l1,z1l2 = best_comb['Z1']
          z2l1,z2l2 = best_comb['Z2']
          wl = best_comb['W']
          id_tmp=[wl.Id,z2l2.Id,z2l1.Id,z1l2.Id,z1l1.Id]
          fake_index1=id_tmp.index(FElectrons_id[0])
          fake_index2=id_tmp.index(FElectrons_id[1])
          fake_index3=id_tmp.index(FElectrons_id[2])
          SR_5L_FakeBit=31 - (1<<(4-fake_index1)) - (1<<(4-fake_index2)) - (1<<(4-fake_index3))

        SR5L_z1l1_id=z1l1.Id
        SR5L_z1l2_id=z1l2.Id
        SR5L_z2l1_id=z2l1.Id
        SR5L_z2l2_id=z2l2.Id
        SR5L_wl_id=wl.Id
        SR5L_z1l1_pt=z1l1.p4.Pt()
        SR5L_z1l1_eta=z1l1.p4.Eta()
        SR5L_z1l1_phi=z1l1.p4.Phi()
        SR5L_z1l1_mass=z1l1.p4.M()
        SR5L_z1l2_pt=z1l2.p4.Pt()
        SR5L_z1l2_eta=z1l2.p4.Eta()
        SR5L_z1l2_phi=z1l2.p4.Phi()
        SR5L_z1l2_mass=z1l2.p4.M()
        SR5L_z2l1_pt=z2l1.p4.Pt()
        SR5L_z2l1_eta=z2l1.p4.Eta()
        SR5L_z2l1_phi=z2l1.p4.Phi()
        SR5L_z2l1_mass=z2l1.p4.M()
        SR5L_z2l2_pt=z2l2.p4.Pt()
        SR5L_z2l2_eta=z2l2.p4.Eta()
        SR5L_z2l2_phi=z2l2.p4.Phi()
        SR5L_z2l2_mass=z2l2.p4.M()
        SR5L_wl_pt=wl.p4.Pt()
        SR5L_wl_eta=wl.p4.Eta()
        SR5L_wl_phi=wl.p4.Phi()
        SR5L_wl_mass=wl.p4.M()
        
        z1_v4=z1l1.p4+z1l2.p4
        z2_v4=z2l1.p4+z2l2.p4
        #wlv_v4=w_v4(wl.p4, met_user, met_phi_user)
        wlv_v4=w_v4(wl.p4, met_CorrePhi, met_phi_CorrePhi)

        SR5L_z1_pt=z1_v4.Pt()
        SR5L_z1_eta=z1_v4.Eta()
        SR5L_z1_phi=z1_v4.Phi()
        SR5L_z1_mass=z1_v4.M()
        SR5L_z2_pt=z2_v4.Pt()
        SR5L_z2_eta=z2_v4.Eta()
        SR5L_z2_phi=z2_v4.Phi()
        SR5L_z2_mass=z2_v4.M()
        SR5L_w_pt=wlv_v4.Pt()
        SR5L_w_eta=wlv_v4.Eta()
        SR5L_w_phi=wlv_v4.Phi()
        SR5L_w_mass=wlv_v4.M()

      SR5L_dR_z1_w = z1_v4.DeltaR(wlv_v4)
      SR5L_dEta_z1_w = abs(z1_v4.Eta()-wlv_v4.Eta())
      SR5L_dPhi_z1_w = z1_v4.DeltaPhi(wlv_v4)
      SR5L_z1w_pt = (z1_v4+wlv_v4).Pt()
      SR5L_z1w_eta = (z1_v4+wlv_v4).Eta()
      SR5L_z1w_phi = (z1_v4+wlv_v4).Phi()
      SR5L_z1w_mass = (z1_v4+wlv_v4).M()
      SR5L_dR_z1_z2 = z1_v4.DeltaR(z2_v4)
      SR5L_dEta_z1_z2 = abs(z1_v4.Eta()-z2_v4.Eta())
      SR5L_dPhi_z1_z2 = z1_v4.DeltaPhi(z2_v4)
      SR5L_z1z2_pt = (z1_v4+z2_v4).Pt()
      SR5L_z1z2_eta = (z1_v4+z2_v4).Eta()
      SR5L_z1z2_phi = (z1_v4+z2_v4).Phi()
      SR5L_z1z2_mass = (z1_v4+z2_v4).M()
      SR5L_dR_w_z2 = z2_v4.DeltaR(wlv_v4)
      SR5L_dEta_w_z2 = abs(z2_v4.Eta()-wlv_v4.Eta())
      SR5L_dPhi_w_z2 = z2_v4.DeltaPhi(wlv_v4)
      SR5L_wz2_pt = (z2_v4+wlv_v4).Pt()
      SR5L_wz2_eta = (z2_v4+wlv_v4).Eta()
      SR5L_wz2_phi = (z2_v4+wlv_v4).Phi()
      SR5L_wz2_mass = (z2_v4+wlv_v4).M()
      SR5L_wzz_pt = (z1_v4+z2_v4+wlv_v4).Pt()
      SR5L_wzz_eta = (z1_v4+z2_v4+wlv_v4).Eta()
      SR5L_wzz_phi = (z1_v4+z2_v4+wlv_v4).Phi()
      SR5L_wzz_mass = (z1_v4+z2_v4+wlv_v4).M()

    self.out.fillBranch("SR_5L_FlavorBit", SR_5L_FlavorBit)
    self.out.fillBranch("SR_5L_FakeBit", SR_5L_FakeBit)
    self.out.fillBranch("SR5L_z1l1_id", SR5L_z1l1_id)
    self.out.fillBranch("SR5L_z1l2_id", SR5L_z1l2_id)
    self.out.fillBranch("SR5L_z2l1_id", SR5L_z2l1_id)
    self.out.fillBranch("SR5L_z2l2_id", SR5L_z2l2_id)
    self.out.fillBranch("SR5L_wl_id", SR5L_wl_id)
    self.out.fillBranch("SR5L_z1l1_pt", SR5L_z1l1_pt)
    self.out.fillBranch("SR5L_z1l1_eta", SR5L_z1l1_eta)
    self.out.fillBranch("SR5L_z1l1_phi", SR5L_z1l1_phi)
    self.out.fillBranch("SR5L_z1l1_mass", SR5L_z1l1_mass)
    self.out.fillBranch("SR5L_z1l2_pt", SR5L_z1l2_pt)
    self.out.fillBranch("SR5L_z1l2_eta", SR5L_z1l2_eta)
    self.out.fillBranch("SR5L_z1l2_phi", SR5L_z1l2_phi)
    self.out.fillBranch("SR5L_z1l2_mass", SR5L_z1l2_mass)
    self.out.fillBranch("SR5L_z2l1_pt", SR5L_z2l1_pt)
    self.out.fillBranch("SR5L_z2l1_eta", SR5L_z2l1_eta)
    self.out.fillBranch("SR5L_z2l1_phi", SR5L_z2l1_phi)
    self.out.fillBranch("SR5L_z2l1_mass", SR5L_z2l1_mass)
    self.out.fillBranch("SR5L_z2l2_pt", SR5L_z2l2_pt)
    self.out.fillBranch("SR5L_z2l2_eta", SR5L_z2l2_eta)
    self.out.fillBranch("SR5L_z2l2_phi", SR5L_z2l2_phi)
    self.out.fillBranch("SR5L_z2l2_mass", SR5L_z2l2_mass)
    self.out.fillBranch("SR5L_wl_pt", SR5L_wl_pt)
    self.out.fillBranch("SR5L_wl_eta", SR5L_wl_eta)
    self.out.fillBranch("SR5L_wl_phi", SR5L_wl_phi)
    self.out.fillBranch("SR5L_wl_mass", SR5L_wl_mass)
    self.out.fillBranch("SR5L_z1_pt", SR5L_z1_pt)
    self.out.fillBranch("SR5L_z1_eta", SR5L_z1_eta)
    self.out.fillBranch("SR5L_z1_phi", SR5L_z1_phi)
    self.out.fillBranch("SR5L_z1_mass", SR5L_z1_mass)
    self.out.fillBranch("SR5L_z2_pt", SR5L_z2_pt)
    self.out.fillBranch("SR5L_z2_eta", SR5L_z2_eta)
    self.out.fillBranch("SR5L_z2_phi", SR5L_z2_phi)
    self.out.fillBranch("SR5L_z2_mass", SR5L_z2_mass)
    self.out.fillBranch("SR5L_w_pt", SR5L_w_pt)
    self.out.fillBranch("SR5L_w_eta", SR5L_w_eta)
    self.out.fillBranch("SR5L_w_phi", SR5L_w_phi)
    self.out.fillBranch("SR5L_w_mass", SR5L_w_mass)
    self.out.fillBranch("SR5L_dR_z1_w", SR5L_dR_z1_w)
    self.out.fillBranch("SR5L_dEta_z1_w", SR5L_dEta_z1_w)
    self.out.fillBranch("SR5L_dPhi_z1_w", SR5L_dPhi_z1_w)
    self.out.fillBranch("SR5L_z1w_pt", SR5L_z1w_pt)
    self.out.fillBranch("SR5L_z1w_eta", SR5L_z1w_eta)
    self.out.fillBranch("SR5L_z1w_phi", SR5L_z1w_phi)
    self.out.fillBranch("SR5L_z1w_mass", SR5L_z1w_mass)
    self.out.fillBranch("SR5L_dR_z1_z2", SR5L_dR_z1_z2)
    self.out.fillBranch("SR5L_dEta_z1_z2", SR5L_dEta_z1_z2)
    self.out.fillBranch("SR5L_dPhi_z1_z2", SR5L_dPhi_z1_z2)
    self.out.fillBranch("SR5L_z1z2_pt", SR5L_z1z2_pt)
    self.out.fillBranch("SR5L_z1z2_eta", SR5L_z1z2_eta)
    self.out.fillBranch("SR5L_z1z2_phi", SR5L_z1z2_phi)
    self.out.fillBranch("SR5L_z1z2_mass", SR5L_z1z2_mass)
    self.out.fillBranch("SR5L_dR_w_z2", SR5L_dR_w_z2)
    self.out.fillBranch("SR5L_dEta_w_z2", SR5L_dEta_w_z2)
    self.out.fillBranch("SR5L_dPhi_w_z2", SR5L_dPhi_w_z2)
    self.out.fillBranch("SR5L_wz2_pt", SR5L_wz2_pt)
    self.out.fillBranch("SR5L_wz2_eta", SR5L_wz2_eta)
    self.out.fillBranch("SR5L_wz2_phi", SR5L_wz2_phi)
    self.out.fillBranch("SR5L_wz2_mass", SR5L_wz2_mass)
    self.out.fillBranch("SR5L_wzz_pt", SR5L_wzz_pt)
    self.out.fillBranch("SR5L_wzz_eta", SR5L_wzz_eta)
    self.out.fillBranch("SR5L_wzz_phi", SR5L_wzz_phi)
    self.out.fillBranch("SR5L_wzz_mass", SR5L_wzz_mass)

WZZ2016apvMC = lambda: WZZProducer("2016apv",'MC')
WZZ2016apvB = lambda: WZZProducer("2016apv",'B')
WZZ2016apvC = lambda: WZZProducer("2016apv",'C')
WZZ2016apvD = lambda: WZZProducer("2016apv",'D')
WZZ2016apvE = lambda: WZZProducer("2016apv",'E')
WZZ2016apvF = lambda: WZZProducer("2016apv",'F')
WZZ2016MC = lambda: WZZProducer("2016",'MC')
WZZ2016F = lambda: WZZProducer("2016",'F')
WZZ2016G = lambda: WZZProducer("2016",'G')
WZZ2016H = lambda: WZZProducer("2016",'H')
WZZ2017MC = lambda: WZZProducer("2017",'MC')
WZZ2017B = lambda: WZZProducer("2017",'B')
WZZ2017C = lambda: WZZProducer("2017",'C')
WZZ2017D = lambda: WZZProducer("2017",'D')
WZZ2017E = lambda: WZZProducer("2017",'E')
WZZ2017F = lambda: WZZProducer("2017",'F')
WZZ2018MC = lambda: WZZProducer("2018",'MC')
WZZ2018A = lambda: WZZProducer("2018",'A')
WZZ2018B = lambda: WZZProducer("2018",'B')
WZZ2018C = lambda: WZZProducer("2018",'C')
WZZ2018D = lambda: WZZProducer("2018",'D')
