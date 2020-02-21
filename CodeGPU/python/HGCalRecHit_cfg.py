import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *

enableGPU = False
from Configuration.ProcessModifiers.gpu_cff import gpu

from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *#HGCalRecHit
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_noise_fC, HGCAL_chargeCollectionEfficiencies

process = cms.Process("TESTgpu", gpu) if enableGPU else cms.Process("TESTnongpu")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2026D46Reco_cff')
#process.load('HeterogeneousCore.CUDAServices.CUDAService_cfi')
process.load('RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi')
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 10 ))

fNames = ['file:/afs/cern.ch/user/b/bfontana/CMSSW_11_0_0_pre11_Patatrack/src/UserCode/Samples/20495.0_CloseByParticleGun_CE_E_Front_200um+CE_E_Front_200um_2026D41_GenSimHLBeamSpotFull+DigiFullTrigger_2026D41+RecoFullGlobal_2026D41+HARVESTFullGlobal_2026D41/step3.root']
keep = 'keep *'
drop = 'drop CSCDetIdCSCALCTPreTriggerDigiMuonDigiCollection_simCscTriggerPrimitiveDigis__HLT'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fNames),
                            inputCommands = cms.untracked.vstring([keep, drop]),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"))

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( False )) #add option for edmStreams

process.HGCalRecHit = HGCalRecHit

fNameOut = 'out'
#convert this to a task!!!!!
process.path = cms.Path( process.HGCalRecHit )

process.out = cms.OutputModule("PoolOutputModule",                                                                               
                               fileName = cms.untracked.string(fNameOut+".root"))

process.outpath = cms.EndPath(process.out)
