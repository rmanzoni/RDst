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
    muons = [ip for ip in event.genp if abs(ip.pdgId())==13 and ip.status()==1]
    
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

        if id0.daughters != [211, 321]: continue

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

for v in already_done.values():
    sorted_all_decays_merged[v] = 0
    reference = map(int, re.findall('[0-9]+', v))
    gammaless_reference = list(filter(lambda x: x != 22, reference)) # don't count FSR photons
    for kk, vv in sorted_all_decays.items():
        undertest = map(int, re.findall('[0-9]+', kk))
        undertest = list(filter(lambda x: x != 22, undertest)) # don't count FSR photons
        if undertest == gammaless_reference:
            sorted_all_decays_merged[v] += vv

sorted_all_decays_merged = OrderedDict(sorted(sorted_all_decays_merged.items(), key=lambda x: x[1], reverse=True))

with open('decay_no_acceptance_fullstat_test.pkl', 'wb') as fout:
#     pickle.dump(sorted_all_decays, fout)
    pickle.dump(sorted_all_decays_merged, fout)


