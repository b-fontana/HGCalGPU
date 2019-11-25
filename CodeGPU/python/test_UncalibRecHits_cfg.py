import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 2 ))

fNames = ['file:/afs/cern.ch/user/b/bfontana/CMSSW_11_0_0_pre11_Patatrack/src/UserCode/Samples/20495.0_CloseByParticleGun_CE_E_Front_200um+CE_E_Front_200um_2026D41_GenSimHLBeamSpotFull+DigiFullTrigger_2026D41+RecoFullGlobal_2026D41+HARVESTFullGlobal_2026D41/step3.root']
keep = 'keep *'
drop = 'drop CSCDetIdCSCALCTPreTriggerDigiMuonDigiCollection_simCscTriggerPrimitiveDigis__HLT'
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(fNames),
                            inputCommands = cms.untracked.vstring([keep, drop]),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"))

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( False ))


process.HeterogeneousHGCalRecHits = cms.EDProducer('HeterogeneousHGCalRecHitsProd', #module label
                                                   message = cms.untracked.string('Hello world!'),
                                                   HGCEEUncalibRecHitsTok =  cms.InputTag('HGCalUncalibRecHit', 'HGCEEUncalibRecHits'),
                                                   HGCHEFUncalibRecHitsTok = cms.InputTag('HGCalUncalibRecHit', 'HGCHEFUncalibRecHits'),
                                                   HGCHEBUncalibRecHitsTok = cms.InputTag('HGCalUncalibRecHit', 'HGCHEBUncalibRecHits'))

fNameOut = 'out'
process.path = cms.Path( process.HeterogeneousHGCalRecHits )
process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(fNameOut+'.root'))
process.outpath = cms.EndPath(process.out)

