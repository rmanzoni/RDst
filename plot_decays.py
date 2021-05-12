import re
import pickle
import numpy as np
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

# ff = open('decay_with_stars.pkl')
# ff = open('decay.pkl')
# ff = open('decay_no_acceptance.pkl')
# ff = open('decay_no_acceptance_fullstat.pkl')
# ff = open('decay_no_acceptance_fullstat_test.pkl')
ff = open('test.pkl')
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

counter = 1
for ichunk in range(len(frequency)%100 + 1):
    if counter > 4: break
    fig, ax = plt.subplots(figsize=(18, 54))
    print('doing chunk', counter)
    ini = ichunk*100
    fin = (ichunk+1)*100
    ax.barh(y_pos[ini:fin], frequency[ini:fin], align='center')
    for i, v in enumerate(frequency[ini:fin]):
        ax.text(v, i, '%.3f%s'%(100.*v, '%'), color='black', fontweight='bold')
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
    fig.tight_layout()

    # plt.savefig('decay_frequency_nicer.pdf')
    plt.savefig('alldecays/decay_frequency_no_acceptance_fullstat_chunk%d.pdf' %counter)
    plt.savefig('alldecays/decay_frequency_no_acceptance_fullstat_chunk%d.png' %counter)
    
#     fig.clf()
    
    counter += 1
