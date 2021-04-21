'''
Cambia il rng seed!
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMRandomNumberGeneratorService
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFastSimRandNumGen

from IOMC.RandomEngine.RandomServiceHelper import  RandomNumberServiceHelper
randHelper =  RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate()
process.RandomNumberGeneratorService.saveFileName =  cms.untracked.string("RandomEngineState.log")


'''
import os
from glob import glob

njobs = 500
# njobs = 2
# out_dir = 'RDst_InclusiveHbToDstMu_GEN_13apr21_v4'
out_dir = 'RDst_InclusiveHbToDstMu_no_accecptance_GEN_19apr21_v2'
# events_per_job = -1
events_per_job = 20000

# template_cfg = "RDst-InclusiveHbToDstMu_TEMPLATE_cfg.py"
# template_fileout = "RDst_InclusiveHbToDstMu_TEMPLATE.root"

template_cfg = "RDst-InclusiveHbToDstMu-no-acceptance_TEMPLATE_cfg.py"
template_fileout = "RDst_InclusiveHbToDstMu-no-acceptance_TEMPLATE.root"

##########################################################################################
##########################################################################################

# make output dir
if not os.path.exists(out_dir):
    try:
        os.makedirs('/pnfs/psi.ch/cms/trivcat/store/user/manzoni/'+out_dir)
    except:
        print('pnfs directory exists')
    os.makedirs(out_dir)
    os.makedirs(out_dir + '/logs')
    os.makedirs(out_dir + '/errs')
            
for ijob in range(njobs):

    tmp_cfg = template_cfg.replace('TEMPLATE', 'chunk%d' %ijob)
    tmp_fileout = template_fileout.replace('TEMPLATE', '%d'%ijob)
    
    #input file
    fin = open(template_cfg, "rt")
    #output file to write the result to
    fout = open("%s/%s" %(out_dir, tmp_cfg), "wt")
    #for each line in the input file
    for line in fin:
        #read replace the string and write to output file
        if   'HOOK_FILE_IN'    in line: fout.write(line.replace('HOOK_FILE_IN'   , files[ijob]))
        elif 'HOOK_MAX_EVENTS' in line: fout.write(line.replace('HOOK_MAX_EVENTS', '%d' %events_per_job))
        elif 'HOOK_FILE_OUT'   in line: fout.write(line.replace('HOOK_FILE_OUT'  , '/scratch/manzoni/%s/%s' %(out_dir, tmp_fileout)))
        else: fout.write(line)
    #close input and output files
    fout.close()
    fin.close()

    to_write = '\n'.join([
        '#!/bin/bash',
        'cd {dir}',
        'scramv1 runtime -sh',
        'mkdir -p /scratch/manzoni/{scratch_dir}',
        'ls /scratch/manzoni/',
        'cmsRun {cfg}',
        'xrdcp /scratch/manzoni/{scratch_dir}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/manzoni/{se_dir}/{fout}',
        'rm /scratch/manzoni/{scratch_dir}/{fout}',
        '',
    ]).format(
        dir           = '/'.join([os.getcwd(), out_dir]), 
        scratch_dir   = out_dir, 
        cfg           = tmp_cfg, 
        se_dir        = out_dir,
        fout          = tmp_fileout
        )

    with open("%s/submitter_chunk%d.sh" %(out_dir, ijob), "wt") as flauncher: 
        flauncher.write(to_write)
    
    command_sh_batch = ' '.join([
        'sbatch', 
        '-p wn', 
        '--account=t3', 
        '-o %s/logs/chunk%d.log' %(out_dir, ijob),
        '-e %s/errs/chunk%d.err' %(out_dir, ijob), 
        '--job-name=%s' %out_dir, 
        '--time=1000', 
        '%s/submitter_chunk%d.sh' %(out_dir, ijob), 
    ])

    print(command_sh_batch)
    os.system(command_sh_batch)

    
    
