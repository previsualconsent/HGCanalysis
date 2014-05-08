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

Submit ntuples production
for i in `seq 0 12 1500`; do
	python scripts/submitLocalAnalysis_cfg.py -t MinBias_v4 -q 8nh -f ${i} -s 12;
done

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

