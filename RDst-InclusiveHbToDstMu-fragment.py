import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    comEnergy = cms.double(13000.0),
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table            = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            list_forced_decays     = cms.vstring(
                'MyD*+',
                'MyD*-',
            ),        
            operates_on_particles = cms.vint32(),
            convertPythiaCodes = cms.untracked.bool(False),
            user_decay_embedded= cms.vstring(
"""
Alias      MyD0        D0
Alias      Myanti-D0   anti-D0
Alias      MyD*-       D*-
Alias      MyD*+       D*+

ChargeConj MyD0   Myanti-D0
ChargeConj MyD*-  MyD*+

Decay MyD0
1.000       K-  pi+           PHSP;
Enddecay
CDecay Myanti-D0

Decay MyD*-
1.000       Myanti-D0 pi-     VSS;
Enddecay
CDecay MyD*+

End
"""
            ),
        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            'SoftQCD:nonDiffractive = on',
            'PTFilter:filter = on', # this turn on the filter
            'PTFilter:quarkToFilter = 5', # PDG id of q quark
            'PTFilter:scaleToFilter = 5.0',
#             'HardQCD:hardbbbar = on',
#             'PhaseSpace:pTHatMin = 100',
		),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters',
        ),
    ),
)

###### Filters ##########
DstFromBhadronFilter = cms.EDFilter(
    "PythiaFilterMultiAncestor",
    DaughterIDs     = cms.untracked.vint32 (  421,   211), # D0, pion+
    DaughterMaxEtas = cms.untracked.vdouble( 1.e9,  2.55),
    DaughterMaxPts  = cms.untracked.vdouble( 1.e9,  1.e9),
    DaughterMinEtas = cms.untracked.vdouble(-1.e9, -2.55),
    DaughterMinPts  = cms.untracked.vdouble( -1.0,   0.3),
    MaxEta          = cms.untracked.double ( 99.0),
    MinEta          = cms.untracked.double (-99.0),
    MinPt           = cms.untracked.double (-1.0),
    MotherIDs       = cms.untracked.vint32 (5), # a b-quark, short for any b-hadron
    ParticleID      = cms.untracked.int32  (413) # D*+
)

D0ToKpiFromDstFilter = cms.EDFilter(
    "PythiaFilterMultiAncestor",
    DaughterIDs     = cms.untracked.vint32 ( -321,   211), # K-, pion+ # this shitty piece of code is blissfully ignorant about CC https://github.com/cms-sw/cmssw/blob/1ae3be438887808e96fe9a07d72009e5db230208/GeneratorInterface/GenFilters/plugins/PythiaFilterMultiAncestor.cc#L190
    DaughterMaxEtas = cms.untracked.vdouble( 2.55,  2.55),
    DaughterMaxPts  = cms.untracked.vdouble( 1.e9,  1.e9),
    DaughterMinEtas = cms.untracked.vdouble(-2.55, -2.55),
    DaughterMinPts  = cms.untracked.vdouble(  0.5,   0.5),
    MaxEta          = cms.untracked.double ( 99.0),
    MinEta          = cms.untracked.double (-99.0),
    MinPt           = cms.untracked.double (-1.0),
    MotherIDs       = cms.untracked.vint32 (413), # D*+
    ParticleID      = cms.untracked.int32  (421) # D0
)

DstMuMaxMassFilter = cms.EDFilter(
    "MCParticlePairFilter",
    ParticleID1    = cms.untracked.vint32(413), # D*+
    ParticleID2    = cms.untracked.vint32(13), # mu
    ParticleCharge = cms.untracked.int32(0), # opposite charge
    MaxInvMass     = cms.untracked.double(5.3),
    MinPt          = cms.untracked.vdouble(-1., 7.),
    MinEta         = cms.untracked.vdouble(-1.e9, -1.55),
    MaxEta         = cms.untracked.vdouble( 1.e9,  1.55),
    Status         = cms.untracked.vint32(2, 1),
)

ProductionFilterSequence = cms.Sequence(
    generator +
    DstFromBhadronFilter +
    D0ToKpiFromDstFilter +
    DstMuMaxMassFilter
)