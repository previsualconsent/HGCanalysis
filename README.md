HGCanalysis
===========

HGCal analysis scripts

To add to your cmssw area

git clone git@github.com:PFCal-dev/HGCanalysis UserCode/HGCanalysis

To update in the remote

git push git@github.com:PFCal-dev/HGCanalysis

Submit local production

python scripts/submitLocalHGCalProduction.py -q 1nd -n 100  -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SingleMuon_v14     -r XXX_PID_XXX:13
python scripts/submitLocalHGCalProduction.py -q 1nd -n 500  -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SingleElectron_v14 -r XXX_PID_XXX:11
python scripts/submitLocalHGCalProduction.py -q 1nd -n 100  -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SinglePion_v14     -r XXX_PID_XXX:211
python scripts/submitLocalHGCalProduction.py -q 1nd -n 100  -c test/runParticleGun_GEN_SIM_cfg.py -o /store/cmst3/group/hgcal/CMSSW/SingleQuark_v14    -r XXX_PID_XXX:1
python scripts/submitLocalHGCalProduction.py -q 1nd -n 2000 -c test/runMinBias_GEN_SIM_cfg.py     -o /store/cmst3/group/hgcal/CMSSW/MinBias_v14

Submit ntuples production
for i in `seq 0 50 1950`; do
	python scripts/submitLocalAnalysis_cfg.py -t MinBias_v14 -q 2nd -f ${i} -s 50;
done

#calibration studies
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SingleElectron_SLHC13_30um_SimHits
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SingleK0L_SLHC13_30um_SimHits
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SinglePion_SLHC13_30um_SimHits

#occupancy studies
pileup=(100 140 200)
cell=(5 10 25)
sdType=(0 1)
for p in ${pileup[@]}; do
   for c in ${cell[@]}; do
	for s in ${sdType[@]}; do
	   python test/analysis/runOccupancyAnalysis.py -i /data/psilva/MinBias_SimHits.root -c ${c} -p ${p} -s ${s} & 
	done
   done
done

