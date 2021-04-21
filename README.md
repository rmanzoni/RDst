# RDst

## instructions to produce the "inclusive" D* mu MC

```
cmsrel CMSSW_10_6_24
cd CMSSW_10_6_24/src
cmsenv
git-cms-addpkg GeneratorInterface/GenFilters
git remote add riccardo git@github.com:rmanzoni/cmssw.git
git fetch riccardo
git cherry-pick fca90cefe98
scram b

# download the fragments
curl -s --insecure https://raw.githubusercontent.com/rmanzoni/RDst/main/RDst-InclusiveHbToDstMu-fragment.py --create-dirs -o Configuration/GenProduction/python/RDst-InclusiveHbToDstMu-fragment.py
curl -s --insecure https://raw.githubusercontent.com/rmanzoni/RDst/main/RDst-InclusiveHbToDstMu-no-acceptance-fragment.py --create-dirs -o Configuration/GenProduction/python/RDst-InclusiveHbToDstMu-no-acceptance-fragment.py

# mutatis mutandis for the other fragment
cmsDriver.py Configuration/GenProduction/python/RDst-InclusiveHbToDstMu-no-acceptance-fragment.py \
--fileout file:RDst-InclusiveHbToDstMu-no-acceptance.root                                         \
--mc                                                                                              \
--eventcontent RAWSIM                                                                             \
--datatier GEN                                                                                    \
--conditions 106X_upgrade2018_realistic_v11_L1v1                                                  \
--beamspot Realistic25ns13TeVEarly2018Collision                                                   \
--step GEN                                                                                        \
--geometry DB:Extended                                                                            \
--era Run2_2018                                                                                   \
--python_filename RDst-InclusiveHbToDstMu-no-acceptance_cfg.py                                    \
--no_exec                                                                                         \
--customise Configuration/DataProcessing/Utils.addMonitoring -n -1

# produce the sample
cmsRun RDst-InclusiveHbToDstMu-no-acceptance_cfg.py
```

## inspect the sample and save the list of decay that produce a D* and a muon from the same b-hadron
It also merges together charge-conjugate processes, processes that differ only by FSR and processes that differ by B oscillations

```
ipython -i inspector.py
```

## plot the results
```
ipython -i plot_decays.py
```


