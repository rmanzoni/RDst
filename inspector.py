from __future__ import print_function
import ROOT
import re
import argparse
import numpy as np
import pickle
from time import time
from datetime import datetime, timedelta
from array import array
from glob import glob
from collections import OrderedDict
from scipy.constants import c as speed_of_light
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import product
# https://pypi.org/project/particle/
import particle
from particle import Particle

parser = argparse.ArgumentParser(description='')
parser.add_argument('--files_per_job', dest='files_per_job', default=2    , type=int)
parser.add_argument('--jobid'        , dest='jobid'        , default=0    , type=int)
parser.add_argument('--verbose'      , dest='verbose'      , action='store_true' )
parser.add_argument('--destination'  , dest='destination'  , default='.'  , type=str)
parser.add_argument('--maxevents'    , dest='maxevents'    , default=-1   , type=int)
args = parser.parse_args()

files_per_job = args.files_per_job
jobid         = args.jobid
verbose       = args.verbose
destination   = args.destination
maxevents     = args.maxevents

diquarks = [
    1103,
    2101,
    2103,
    2203,
    3101,
    3103,
    3201,
    3203,
    3303,
    4101,
    4103,
    4201,
    4203,
    4301,
    4303,
    4403,
    5101,
    5103,
    5201,
    5203,
    5301,
    5303,
    5401,
    5403,
    5503,
]

excitedBs = [
    513,
    523,
    533,
    543,
    # others?
]

def isAncestor(a, p):
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
    return False

def printAncestors(particle, ancestors=[], verbose=True):
    for i in xrange(0, particle.numberOfMothers()):
        mum = particle.mother(i)
#         if mum is None: import pdb ; pdb.set_trace()
        if abs(mum.pdgId())<8 or \
           abs(mum.pdgId())==21 or \
           abs(mum.pdgId()) in diquarks or\
           abs(mum.pdgId()) in excitedBs: continue
        # don't count B oscillations
        if mum.pdgId() == -particle.pdgId() and abs(particle.pdgId()) in [511, 531]:
            continue 
        if not mum.isLastCopy(): continue
        try:
            if verbose: print(' <-- ', Particle.from_pdgid(mum.pdgId()).name, end = '')
            ancestors.append(mum)
            printAncestors(mum, ancestors=ancestors, verbose=verbose)
        except:
            if verbose: print(' <-- ', 'pdgid', mum.pdgId(), end = '')
            ancestors.append(mum)
            printAncestors(mum, ancestors=ancestors, verbose=verbose)
        else:
            pass

def printOffspring(particle, indent=0, file=None):
    for idau in range(particle.numberOfDaughters()):
        daughter = particle.daughter(idau)
        try:
            print(' '*8*indent, '|---->', Particle.from_pdgid(daughter.pdgId()).name, '\t pt {0:>5}  eta {1:>5}  phi {2:>5}'.format('%.1f'%daughter.pt(), '%.2f'%daughter.eta(), '%.2f'%daughter.phi()), file=file)
        except:
            print(' '*8*indent, '|---->', 'pdgid', daughter.pdgId(), '\t pt {0:>5}  eta {1:>5}  phi {2:>5}'.format('%.1f'%daughter.pt(), '%.2f'%daughter.eta(), '%.2f'%daughter.phi()), file=file)
        printOffspring(daughter, indent+1, file=file)
    
def decayChainToLabel(particle, label='', generation=0):
    if generation==0:
        label += '%d' %particle.pdgId()
    generation += 1
    for idau in range(particle.numberOfDaughters()):
        if idau==0: 
            label += '('
        daughter = particle.daughter(idau)
        # don't count FSR photon. Fingers crossed, hope this is the correct way...        
#         if daughter.pdgId()==22 and \
#            particle.isLastCopy() and \
#            particle.isLastCopyBeforeFSR():
#             continue
        # nope... not correct
        label += '%d,' %daughter.pdgId()
        if daughter.numberOfDaughters()>0:
            label = decayChainToLabel(daughter, label, generation=generation)
        if idau==particle.numberOfDaughters()-1:
            label += '),'
    generation -= 1
    
    label = label.replace(',(', '(')
    label = label.replace(',)', ')')
    if label.endswith(','):
        label = label[:-1]

    return label
    
handles = OrderedDict()
handles['genp'   ] = ('genParticles', Handle('std::vector<reco::GenParticle>'))
handles['genInfo'] = ('generator'   , Handle('GenEventInfoProduct'           ))

files  = glob('/pnfs/psi.ch/cms/trivcat/store/user/manzoni/RDst_InclusiveHbToDstMu_no_accecptance_GEN_19apr21_v2/*root')
# events = Events(files[:20])
events = Events(files)

logfile = open('RDst_InclusiveHbToDstMu_no_acceptance_fullstat_test.txt', 'w')

start = time()
# maxevents = 2e5
maxevents = maxevents if maxevents>=0 else events.size() # total number of events in the files

alldecays = dict()

branches = [
    'run',
    'lumi',
    'event',
    
    'tmpDecayIndex',
    
    'm2_miss',
    'q2',
    'e_star_mu3',
    
    'mu_pt',
    'mu_eta',
    'mu_phi',
    'mu_E',
    'mu_mass',
    'mu_charge',
    
    'pi_dst_pt',
    'pi_dst_eta',
    'pi_dst_phi',
    'pi_dst_E',
    'pi_dst_mass',
    'pi_dst_charge',
    
    'pi_d0_pt',
    'pi_d0_eta',
    'pi_d0_phi',
    'pi_d0_E',
    'pi_d0_mass',
    'pi_d0_charge',
    
    'k_pt',
    'k_eta',
    'k_phi',
    'k_E',
    'k_mass',
    'k_charge',
]

decay_index = {}
next_decay_index = 0

fout = ROOT.TFile('decay_info.root', 'recreate')
ntuple = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
decayIndexNtuple = ROOT.TNtuple('decayIndexTree', 'decayIndexTree', 'decayIndex')
tofill = OrderedDict(zip(branches, [np.nan]*len(branches)))

for i, event in enumerate(events):

    if (i+1)>maxevents:
        break
        
    if i%100==0:
        percentage = float(i)/maxevents*100.
        speed = float(i)/(time()-start)
        eta = datetime.now() + timedelta(seconds=(maxevents-i) / max(0.1, speed))
        print('\t===> processing %d / %d event \t completed %.1f%s \t %.1f ev/s \t ETA %s s' %(i, maxevents, percentage, '%', speed, eta.strftime('%Y-%m-%d %H:%M:%S')))

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    lumi = event.eventAuxiliary().luminosityBlock()
    iev  = event.eventAuxiliary().event()
   
    event.qscale = event.genInfo.qScale()
   
    if verbose: print('=========>')
    
    dsts  = [ip for ip in event.genp if abs(ip.pdgId())==413]
#     muons = [ip for ip in event.genp if abs(ip.pdgId())==13 and ip.status()==1 and ip.pt()>7. and abs(ip.eta())<1.5]
    muons = sorted([ip for ip in event.genp if abs(ip.pdgId())==13 and ip.status()==1], key = lambda ip: ip.pt(), reverse = True)
    
    used_ancestors = set()
    for ids, imu in product(dsts, muons):
#         if iev == 2333: import pdb ; pdb.set_trace()
            
#         if (ids.p4() + imu.p4()).mass()>5.4: continue
        if (ids.p4() + imu.p4()).mass()>7.: continue
        
        if not len(getattr(imu, 'ancestors', [])):
            ancestorsmu = []
            printAncestors(imu, ancestorsmu, verbose=False)
            imu.ancestors = ancestorsmu
        
        if not len(getattr(ids, 'ancestors', [])):
            ancestorsds = []
            printAncestors(ids, ancestorsds, verbose=False)
            ids.ancestors = ancestorsds

        if len(ids.ancestors)==0 or len(imu.ancestors)==0: continue
        if (imu.charge()*ids.charge() > 0): continue #same charge, cannot be a single B0

        if ids.ancestors[-1] != imu.ancestors[-1]: continue
        
        ids.daughters = [abs(ids.daughter(idau).pdgId()) for idau in range(ids.numberOfDaughters())]
        ids.daughters = [idau for idau in ids.daughters if idau!=22] # remove photons
        ids.daughters.sort()

#         if iev == 2333: import pdb ; pdb.set_trace()
        
        if ids.daughters != [211, 421]: continue

        for idau in range(ids.numberOfDaughters()):
            if abs(ids.daughter(idau).pdgId())==421:
                id0 = ids.daughter(idau)
                id0.daughters = [abs(id0.daughter(idau).pdgId()) for idau in range(id0.numberOfDaughters())]
                id0.daughters = [idau for idau in id0.daughters if idau!=22] # remove photons
                id0.daughters.sort()
            elif abs(ids.daughter(idau).pdgId())==211:
                ipi_dst = ids.daughter(idau)

        if id0.daughters != [211, 321]: continue
        
        for idau in range(id0.numberOfDaughters()):
            if abs(id0.daughter(idau).pdgId())==321:
                ik = id0.daughter(idau)
            elif abs(id0.daughter(idau).pdgId())==211:
                ipi_d0 = id0.daughter(idau)

#         if iev == 2333: import pdb ; pdb.set_trace()
        
        
        print('\n\n'+'-'*80, file=logfile)
        print('lumi %d \t event %d' %(lumi, iev), file=logfile)

        ib = ids.ancestors[-1]
        print('first b-hadron ancestor: ', end='', file=logfile)
        try:
            print(Particle.from_pdgid(ib.pdgId()), end='', file=logfile)
        except:
            print(ib.pdgId(), end='', file=logfile)
        print('\t pt {0:>5}  eta {1:>5}  phi {2:>5}'.format('%.1f'%ib.pt(), '%.2f'%ib.eta(), '%.2f'%ib.phi()), file=logfile)
        
        printOffspring(ib, file=logfile)
#         import pdb ; pdb.set_trace()
        decayName = decayChainToLabel(ib, label='')
#         print(decayName)
        if decayName in alldecays.keys():
            alldecays[decayName] += 1
        else:
            alldecays[decayName] = 1
            decay_index[decayName] = next_decay_index
            next_decay_index += 1
        
        if ids.ancestors[-1] in used_ancestors:
            # a different pair with this ancestor was already registered in this event
            # possibly a decay with 2 muons and a Dst
            # we keep only the highest pt muon, as we sorted the muons before
            continue

        used_ancestors.add(ids.ancestors[-1])
        
        b_lab_p4 = imu.p4() + ids.p4()
        b_scaled_p4 = b_lab_p4 * ((particle.literals.B_0.mass/1000)/b_lab_p4.mass())
        
        b_scaled_p4_tlv = ROOT.TLorentzVector() ; b_scaled_p4_tlv.SetPtEtaPhiE(b_scaled_p4.pt(), b_scaled_p4.eta(), b_scaled_p4.phi(), b_scaled_p4.energy())
        imu_p4_tlv = ROOT.TLorentzVector() ; imu_p4_tlv.SetPtEtaPhiE(imu.pt(), imu.eta(), imu.phi(), imu.energy())
        
        b_scaled_p4_boost = b_scaled_p4_tlv.BoostVector()
        
        imu_p4_in_b_rf = imu_p4_tlv.Clone(); imu_p4_in_b_rf.Boost(-b_scaled_p4_boost)
        
        tofill['run'          ] = event.eventAuxiliary().run()
        tofill['lumi'         ] = event.eventAuxiliary().luminosityBlock()
        tofill['event'        ] = event.eventAuxiliary().event()
        
        tofill['m2_miss'      ] = (b_scaled_p4 - imu.p4() - ids.p4()).mass2()
        tofill['q2'           ] = (b_scaled_p4 - ids.p4()).mass2()
        tofill['e_star_mu3'   ] = imu_p4_in_b_rf.E()
        
        tofill['mu_pt        '] = imu.pt()
        tofill['mu_eta       '] = imu.eta()
        tofill['mu_phi       '] = imu.phi()
        tofill['mu_E         '] = imu.energy()
        tofill['mu_mass      '] = imu.mass()
        tofill['mu_charge    '] = imu.charge()
        
        tofill['pi_dst_pt    '] = ipi_dst.pt()
        tofill['pi_dst_eta   '] = ipi_dst.eta()
        tofill['pi_dst_phi   '] = ipi_dst.phi()
        tofill['pi_dst_E     '] = ipi_dst.energy()
        tofill['pi_dst_mass  '] = ipi_dst.mass()
        tofill['pi_dst_charge'] = ipi_dst.charge()
        
        tofill['pi_d0_pt     '] = ipi_d0.pt()
        tofill['pi_d0_eta    '] = ipi_d0.eta()
        tofill['pi_d0_phi    '] = ipi_d0.phi()
        tofill['pi_d0_E      '] = ipi_d0.energy()
        tofill['pi_d0_mass   '] = ipi_d0.mass()
        tofill['pi_d0_charge '] = ipi_d0.charge()
        
        tofill['k_pt         '] = ik.pt()
        tofill['k_eta        '] = ik.eta()
        tofill['k_phi        '] = ik.phi()
        tofill['k_E          '] = ik.energy()
        tofill['k_mass       '] = ik.mass()
        tofill['k_charge     '] = ik.charge()
        
        tofill['tmpDecayIndex'] = decay_index[decayName]
        
        ntuple.Fill(array('f', tofill.values()))
        
logfile.close()

# sorted_alldecays = {k: v for k, v in sorted(alldecays.items(), key=lambda item: item[1])}
alldecays = sorted(alldecays.items(), key=lambda x: x[1], reverse=True)

sorted_all_decays = OrderedDict(alldecays)

# merge charge conjugates and FSR
sorted_all_decays_merged = OrderedDict()
already_done = OrderedDict()


for k in sorted_all_decays.keys():
    reference = map(int, re.findall('[0-9]+', k))
    gammaless_reference = list(filter(lambda x: x != 22, reference)) # don't count FSR photons
    if repr(gammaless_reference) in already_done.keys():
        if len(k.replace('-','')) > len(already_done[repr(gammaless_reference)].replace('-','')):
            # print (k, 'longer than', already_done[repr(gammaless_reference)]) 
            continue
    already_done[repr(gammaless_reference)] = k

decay_merge_map = {}
for v in already_done.values():
    sorted_all_decays_merged[v] = 0
    reference = map(int, re.findall('[0-9]+', v))
    gammaless_reference = list(filter(lambda x: x != 22, reference)) # don't count FSR photons
    for kk, vv in sorted_all_decays.items():
        undertest = map(int, re.findall('[0-9]+', kk))
        undertest = list(filter(lambda x: x != 22, undertest)) # don't count FSR photons
        if undertest == gammaless_reference:
            sorted_all_decays_merged[v] += vv
            decay_merge_map[kk] = v


sorted_all_decays_merged = OrderedDict(sorted(sorted_all_decays_merged.items(), key=lambda x: x[1], reverse=True))

sorted_decay_index = {decay_name: index for index, decay_name in enumerate(sorted_all_decays_merged.keys())}

with open("decay_dictionary", "w") as decay_dict_out:
    for index, decay_name in enumerate(sorted_all_decays_merged.keys()):
        print(index, decay_name, file = decay_dict_out)

decay_index_ordering = {raw_decay_index: sorted_decay_index[decay_merge_map[decay_name]] for decay_name, raw_decay_index in decay_index.items()}

for event in ntuple:
    decayIndexNtuple.Fill(decay_index_ordering[event.tmpDecayIndex])

fout.cd()
ntuple.AddFriend(decayIndexNtuple)
ntuple.Write()
decayIndexNtuple.Write()
fout.Close()

with open('decay_no_acceptance_fullstat_test.pkl', 'wb') as fout:
#     pickle.dump(sorted_all_decays, fout)
    pickle.dump(sorted_all_decays_merged, fout)



