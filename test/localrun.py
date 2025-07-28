import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.WZZProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=True, action='store_true')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  parser.add_option('-n','--nEve', dest='nEvent', help='number of event', type='int', action='store')
  parser.add_option('-i', '--in', dest='inputs', help='input directory with files', default=None, type='string')
  parser.add_option('-o', '--out', dest='output', help='output directory with files', default=None, type='string')
  (opt, args) = parser.parse_args()

  if opt.ismc:
    if opt.year == "2016a":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2016_preAPV(),muonIDISOSF2016apv(),muonScaleRes2016a(),eleRECOSF2016apv(),eleIDSF2016apv(),jmeCorrections_UL2016APVMC(), WZZ2016apvMC()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2016_postAPV(),muonIDISOSF2016(),muonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(),jmeCorrections_UL2016MC(),WZZ2016MC()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2017(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(),jmeCorrections_UL2017MC(),WZZ2017MC()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018":
      p = PostProcessor(opt.output, [opt.inputs], modules=[countHistogramsModule(),puWeight_2018(),muonIDISOSF2018(),muonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(),jmeCorrections_UL2018MC(),WZZ2018MC()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016B(),WZZ2016apvB()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016c":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016C(),WZZ2016apvC()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016d":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016D(),WZZ2016apvD()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016e":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016E(),WZZ2016apvE()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016f_apv":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016a(),jmeCorrections_UL2016APVF(),WZZ2016apvF()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016f":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016b(),jmeCorrections_UL2016F(),WZZ2016F()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016g":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016b(),jmeCorrections_UL2016G(),WZZ2016G()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2016h":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2016b(),jmeCorrections_UL2016H(),WZZ2016H()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017B(),WZZ2017B()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017c":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017C(),WZZ2017C()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017d":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017D(),WZZ2017D()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017e":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017E(),WZZ2017E()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2017f":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2017(),jmeCorrections_UL2017F(),WZZ2017F()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018a":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018A(),WZZ2018A()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018b":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018B(),WZZ2018B()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018c":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018C(),WZZ2018C()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
    if opt.year == "2018d":
      p = PostProcessor(opt.output, [opt.inputs], modules=[muonScaleRes2018(),jmeCorrections_UL2018D(),WZZ2018D()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=opt.nEvent)
  p.run()

if __name__ == "__main__":
    sys.exit(main())
