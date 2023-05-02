import FWCore.ParameterSet.Config as cms

from DQM.Physics.ExoticaDQM_miniAOD_cfi import *

exoticaPhysicsminiAOD = cms.Sequence(    ExoticaDQM_miniAOD 
                                    )
