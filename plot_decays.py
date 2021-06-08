import re
import pickle
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from particle import Particle

def particle_to_latex(pdgid):
    pdgid = int(pdgid)
    try:
        return ' '+Particle.from_pdgid(pdgid).latex_name
    except:
        if pdgid==20423:
            return 'D_{1}(2430)^{0}'
        elif pdgid==-20423:
            return '\\bar{D}_{1}(2430)^{0}'
        elif pdgid==20413:
            return 'D_{1}(H)^{+}'
        elif pdgid==-20413:
            return 'D_{1}(H)^{-}'
        elif pdgid==10413:
            return 'D_{1}(2420)^{+}'
        elif pdgid==-10413:
            return 'D_{1}(2420)^{-}'
        elif pdgid==5214:
            return '\\Sigma*_{b}^{0}'
        elif pdgid==-5214:
            return '\\bar{\\Sigma}*_{b}^{0}'
        elif pdgid==5212:
            return '\\Sigma_{b}^{0}'
        elif pdgid==-5212:
            return '\\bar{\\Sigma}*_{b}^{0}'
        elif pdgid==5314:
            return '\\Xi*_{b}^{-}'
        elif pdgid==-5314:
            return '\\Xi*_{b}^{+}'
        elif pdgid==5312:
            return "\\Xi'_{b}^{-}"
        elif pdgid==-5312:
            return "\\Xi'_{b}^{+}"
        elif pdgid==5324:
            return '\\Xi*_{b}^{0}'
        elif pdgid==-5324:
            return '\\bar{\\Xi}*_{b}^{0}'
        elif pdgid==5322:
            return '\\Xi_{b}^{\\prime 0}'
        elif pdgid==-5322:
            return '\\bar{\\Xi}_{b}^{\\prime 0}'
        elif pdgid==543:
            return 'B*_{c}^+'
        elif pdgid==-543:
            return '\\bar{B}*_{c}^+'
        print('pdg not found', pdgid)    
        return ' %d'%pdgid
#         import pdb ; pdb.set_trace()

def relabel(label):
    newlabel = r'$'
    pdgid = ''
    for i in label:
        if i == '(':
            if len(pdgid): newlabel += particle_to_latex(pdgid)
            newlabel += '(\\rightarrow '
            pdgid = ''
        elif i == ')':
            if len(pdgid): newlabel += particle_to_latex(pdgid)
            newlabel += ')'
            pdgid = ''
        elif i == ',':
            if len(pdgid): newlabel += particle_to_latex(pdgid)
            pdgid = ''
        else:
            pdgid += i
    return newlabel+'$'
    
parser = argparse.ArgumentParser(description='')
parser.add_argument('inputFile', type=str)
parser.add_argument('--outputFile', type=str)
args = parser.parse_args()

# ff = open('decay_with_stars.pkl')
# ff = open('decay.pkl')
# ff = open('decay_no_acceptance.pkl')
# ff = open('decay_no_acceptance_fullstat.pkl')
# ff = open('decay_no_acceptance_fullstat_test.pkl')
# ff = open('test.pkl')
#ff = open('test_17may.pkl')
#ff = open('decay_test-QCD_Pt15to20.pkl')
ff = open(args.inputFile)

decays = pickle.load(ff)
ff.close()

plt.rcdefaults()
# fig, ax = plt.subplots(figsize=(100, 200))
# fig, ax = plt.subplots(figsize=(18, 8))
# fig, ax = plt.subplots(figsize=(18, 18))
# fig, ax = plt.subplots(figsize=(18, 36))
# fig, ax = plt.subplots(figsize=(18, 54))
# fig, ax = plt.subplots(figsize=(18, 72))
# fig, ax = plt.subplots(figsize=(22, 1500))
# fig, ax = plt.subplots()

alldecays = sorted(decays.items(), key=lambda x: x[1], reverse=True)
total_events = np.sum(np.array([idecay[1] for idecay in alldecays]))
# alldecays = alldecays[:25]
# alldecays = alldecays[:50]
# alldecays = alldecays[:100]

decays = [idecay[0] for idecay in alldecays] #list(alldecays.keys())
occurrences = [idecay[1] for idecay in alldecays]
y_pos = np.arange(len(decays))
frequency = np.array(occurrences).astype(np.float32)/total_events
newlabels = map(relabel, decays)
newlabels = [newlabels[ii]+'   idx %4d'%ii for ii in range(len(newlabels))]

print('total decays', len(decays))

chunk_size = 120
if len(frequency) < chunk_size:
    chunk_size = len(frequency)

counter = 1
tot = len(frequency) % chunk_size + 1
for ichunk in range(tot):
    if counter > 4: break
    print('doing chunk {0} of {1}'.format(counter, tot))
    fig, ax = plt.subplots(figsize=(18, 54))
    ini = ichunk * chunk_size
    fin = (ichunk+1) * chunk_size
    mymin = min(1.e-5, 0.5 * np.min(frequency[ini:fin]))
    mymax = min(1., 2. * np.max(frequency[ini:fin]))
    ax.barh(y_pos[ini:fin], frequency[ini:fin], align='center', xerr=np.sqrt(occurrences[ini:fin])/total_events)
    for i, v in enumerate(frequency[ini:fin]):
        vv = occurrences[ini:fin][i]
        ax.text(1.1*mymin, i, '%.3f%s - %d events'%(100.*v, '%', vv), color='lightgray', fontweight='bold')
#         ax.text(v, i, '%.3f%s'%(100.*v, '%'), color='black', fontweight='bold')
#         ax.text(v/2., i, '%d ev.'%(vv), color='white', fontweight='bold')
    print('added {0} decays'.format(chunk_size))
    ax.set_yticks(y_pos[ini:fin])
    # ax.set_yticklabels(map(relabel, decays))
    ax.set_yticklabels(newlabels[ini:fin])
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Frequency')
    # ax.set_title(r'B hadron decays giving rise to $D*(2010)^{+}(\rightarrow D^{0}(\rightarrow K^{-}\pi^{+})\pi^{+})\mu$ in the acceptance. Charge conjugation implied')
    ax.set_title(r'B hadron decays giving rise to $D*(2010)^{+}(\rightarrow D^{0}(\rightarrow K^{-}\pi^{+})\pi^{+})\mu$. Charge conjugation implied')
    ax.set_xscale('log')
    # ax.set_aspect(aspect=10)
    # ax.set_xlim(300, 400)
    # ax.set_box_aspect(10)
    ax.margins(y=0.001)
    plt.xlim(mymin, mymax)
    
    textstr = '\n'.join([
        'tot events %d' %total_events,
        'tot events in this chunk %d' %(np.sum(occurrences[ini:fin])),
        'tot fraction in this chunk %.3f%s' %(100.*np.sum(frequency[ini:fin]), '%'),
    ])

    # place a text box in upper left in axes coords
#     ax.text(0.6 * (mymax-mymin), 0.92, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top')
    print('adding summary')
    plt.text(0.05, 0.98, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top')
    
    fig.tight_layout()

    print('saving plot')
    # plt.savefig('decay_frequency_nicer.pdf')
    plt.savefig(args.outputFile + '_chunk' + str(counter) + '.pdf')
    # just pdf to halve running time
    #plt.savefig('allDecays/decay_frequency_no_acceptance_fullstat_chunk%d.png' %counter)
    
#     fig.clf()
    
    counter += 1

