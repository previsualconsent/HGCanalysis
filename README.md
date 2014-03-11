HGCanalysis
===========

HGCal analysis scripts

To add to your cmssw area

git clone git@github.com:PFCal-dev/HGCanalysis UserCode/HGCanalysis

To update in the remote

git push git@github.com:PFCal-dev/HGCanalysis

Submit local production

python scripts/submitLocalHGCalProduction.py -q 1nd -n 500 -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SingleElectron_v1 -r XXX_PID_XXX:11
python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SingleQuark_v1 -r XXX_PID_XXX:1
python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SingleMuon_v1 -r XXX_PID_XXX:13