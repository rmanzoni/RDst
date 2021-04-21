# RDst

instructions to produce the "inclusive" D* mu MC

```
cmsrel CMSSW_10_6_24
cd CMSSW_10_6_24/src
cmsenv
git-cms-addpkg GeneratorInterface/GenFilters
git remote add riccardo git@github.com:rmanzoni/cmssw.git
git fetch riccardo
git cherry-pick fca90cefe98
scram b
```

