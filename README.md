HGCanalysis
===========

HGCal analysis scripts

To add to your cmssw area

git clone git@github.com:PFCal-dev/HGCanalysis UserCode/HGCanalysis

To update in the remote

git push git@github.com:PFCal-dev/HGCanalysis

GEN-SIM-RECO production

Particle gun: use the script scripts/generateParticleGun.sh to generate locally a flat eta, fixed energy particle gun

Full process: use the script scripts/generateSample.sh to generate other processes with existing cfi (to be updated)

Submit production to the batch:

energies=(5 10 20 30 50 100 250)
pids=(11 13 211)
for pid in ${pids[@]}; do
for en in ${energies[@]}; do
	python scripts/submitLocalHGCalProduction.py -q 1nd -n 3 -s generateParticleGun.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${PID}_SLHC16 -p ${pid} -n 500 -e ${en}";
done
done

#check this point forward

Submit ntuples production
for i in `seq 0 50 1950`; do
	python scripts/submitLocalAnalysis_cfg.py -t MinBias_v14 -q 2nd -f ${i} -s 50;
done

#calibration studies
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SingleElectron_SLHC13_30um_SimHits
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SingleK0L_SLHC13_30um_SimHits
python test/analysis//testSimHits.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples -t SinglePion_SLHC13_30um_SimHits

#occupancy studies
pileup=(140 200)
noise=(0) # it is generated in the summary
sdType=(0 1)
for p in ${pileup[@]}; do
for n in ${noise[@]}; do 
for s in ${sdType[@]}; do
python test/analysis/runOccupancyAnalysis.py -i /store/cmst3/group/hgcal/CMSSW/Ntuples/ -t MinBias_v16 -p ${p} -s ${s} -n ${n} &
done
done
done

for p in ${pileup[@]}; do
for n in ${noise[@]}; do
for s in ${sdType[@]}; do
python test/analysis/drawOccupancyAnalysisSummary.py -i MinBias_v16_occ_pu${p}_sd${s}.root -o pu${p}_sd${s} &
done
done
done