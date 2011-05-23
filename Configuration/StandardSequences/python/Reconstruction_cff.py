import FWCore.ParameterSet.Config as cms

from RecoLuminosity.LumiProducer.lumiProducer_cff import *
from RecoLocalMuon.Configuration.RecoLocalMuon_cff import *
from RecoLocalCalo.Configuration.RecoLocalCalo_cff import *
from RecoTracker.Configuration.RecoTracker_cff import *
from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
from TrackingTools.Configuration.TrackingTools_cff import *
# Global  reco
from RecoEcal.Configuration.RecoEcal_cff import *
from RecoJets.Configuration.CaloTowersRec_cff import *
from RecoMET.Configuration.RecoMET_cff import *
from RecoMuon.Configuration.RecoMuon_cff import *
# Higher level objects
from RecoVertex.Configuration.RecoVertex_cff import *
from RecoEgamma.Configuration.RecoEgamma_cff import *
from RecoPixelVertexing.Configuration.RecoPixelVertexing_cff import *


from RecoJets.Configuration.RecoJetsGlobal_cff import *
from RecoMET.Configuration.RecoPFMET_cff import *
from RecoBTag.Configuration.RecoBTag_cff import *
#
# please understand that division global,highlevel is completely fake !
#
#local reconstruction
from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
from RecoParticleFlow.Configuration.RecoParticleFlow_cff import *
#
# new tau configuration
#
from RecoTauTag.Configuration.RecoPFTauTag_cff import *
from RecoTauTag.Configuration.RecoTauTag_cff import *
# Also BeamSpot
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *

from RecoLocalCalo.CastorReco.CastorSimpleReconstructor_cfi import *

localreco = cms.Sequence(trackerlocalreco+muonlocalreco+calolocalreco+castorreco+lumiProducer)
localreco_HcalNZS = cms.Sequence(trackerlocalreco+muonlocalreco+calolocalrecoNZS+castorreco+lumiProducer)

#
# temporarily switching off recoGenJets; since this are MC and wil be moved to a proper sequence
#

from RecoLocalCalo.Castor.Castor_cff import *
from RecoLocalCalo.Configuration.hcalGlobalReco_cff import *

globalreco = cms.Sequence(offlineBeamSpot*
                          recopixelvertexing*
                          trackingGlobalReco*
                          hcalGlobalRecoSequence*
                          particleFlowCluster*
                          ecalClusters*
                          caloTowersRec*
                          vertexreco*
                          egammaGlobalReco*
                          pfTrackingGlobalReco*
                          jetGlobalReco*
                          muonrecoComplete*
                          muoncosmicreco*
                          CastorFullReco)

globalreco_plusRS = cms.Sequence(globalreco*rstracks)
globalreco_plusPL= cms.Sequence(globalreco*ctfTracksPixelLess)


reducedRecHits = cms.Sequence ( reducedEcalRecHitsSequence * reducedHcalRecHitsSequence )

highlevelreco = cms.Sequence(egammaHighLevelRecoPrePF*
                             particleFlowReco*
                             egammaHighLevelRecoPostPF*
                             jetHighLevelReco*
                             tautagging*
                             regionalCosmicTracksSeq*
                             muoncosmichighlevelreco*
                             metrecoPlusHCALNoise*
                             btagging*
                             recoPFMET*
                             PFTau*
                             reducedRecHits
                             )


from FWCore.Modules.logErrorHarvester_cfi import *

# "Export" Section
reconstruction         = cms.Sequence(localreco        *globalreco       *highlevelreco*logErrorHarvester)

#need a fully expanded sequence copy
modulesToRemove = list() # copy does not work well
modulesToRemove.append(siPixelClusters)
modulesToRemove.append(siStripZeroSuppression)
modulesToRemove.append(siStripClusters)
modulesToRemove.append(dt1DRecHits)
modulesToRemove.append(dt1DCosmicRecHits)
modulesToRemove.append(csc2DRecHits)
modulesToRemove.append(rpcRecHits)
modulesToRemove.append(ecalGlobalUncalibRecHit)
modulesToRemove.append(ecalDetIdToBeRecovered)
modulesToRemove.append(ecalRecHit)
modulesToRemove.append(ecalCompactTrigPrim)
modulesToRemove.append(ecalTPSkim)
modulesToRemove.append(ecalPreshowerRecHit)
modulesToRemove.append(selectDigi)
modulesToRemove.append(hbheprereco)
modulesToRemove.append(hbhereco)
modulesToRemove.append(hfreco)
modulesToRemove.append(horeco)
modulesToRemove.append(hcalnoise)
modulesToRemove.append(zdcreco)
modulesToRemove.append(castorreco)
reconstruction_fromRECO = reconstruction.copyAndExclude(modulesToRemove)
modulesToRemove.append(siPixelRecHits)
modulesToRemove.append(siStripMatchedRecHits)
modulesToRemove.append(pixelTracks)
modulesToRemove.append(ckftracks)
modulesToRemove.append(trackerDrivenElectronSeeds)
modulesToRemove.append(ecalDrivenElectronSeeds)
modulesToRemove.append(electronMergedSeeds)
modulesToRemove.append(electronCkfTrackCandidates)
modulesToRemove.append(electronGsfTracks)
modulesToRemove.append(generalConversionTrackProducer)
modulesToRemove.append(generalInOutOutInConversionTrackMerger)
modulesToRemove.append(gsfGeneralInOutOutInConversionTrackMerger)
modulesToRemove.append(ecalDrivenGsfElectrons)
modulesToRemove.append(gsfConversionTrackProducer)
modulesToRemove.append(allConversions)
modulesToRemove.append(gsfElectrons)
reconstruction_fromRECO_noTracking = reconstruction.copyAndExclude(modulesToRemove)


#sequences with additional stuff
reconstruction_withPixellessTk  = cms.Sequence(localreco        *globalreco_plusPL*highlevelreco*logErrorHarvester)
reconstruction_withRS  = cms.Sequence(localreco        *globalreco_plusRS*highlevelreco*logErrorHarvester)
reconstruction_HcalNZS = cms.Sequence(localreco_HcalNZS*globalreco       *highlevelreco*logErrorHarvester)

#sequences without some stuffs
#
reconstruction_woCosmicMuons = cms.Sequence(localreco*globalreco*highlevelreco       *logErrorHarvester)


# define a standard candle. please note I am picking up individual
# modules instead of sequences
#
reconstruction_standard_candle = cms.Sequence(localreco*globalreco*vertexreco*recoJetAssociations*btagging*electronSequence*photonSequence)


