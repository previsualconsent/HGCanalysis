# HGCanalysis scripts

## Installation 

To add to your cmssw area

git clone git@github.com:PFCal-dev/HGCanalysis UserCode/HGCanalysis

To update in the remote

git push git@github.com:PFCal-dev/HGCanalysis

## GEN-SIM-RECO production

Use the generateEventsFromCfi.sh to steer the generation.
It will call cmsDriver.py and customize the configuration file.
The options can be inspected by calling:

generateEventsFromCfi.sh -h

A full production can be ran locally or submitted to the batch using 
the submitLocalHGCalProduction.py wrapper script. Two examples are given below:

### Particle gun 

energies=(5 10 20 30 50 100 250)
pids=(11 13 211)
for pid in ${pids[@]}; do
for en in ${energies[@]}; do
	python scripts/submitLocalHGCalProduction.py -q 1nd -n 5 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/Single${PID}_SLHC16 -p ${pid} -n 500 -e ${en}";
done
done

### Minimum bias

python scripts/submitLocalHGCalProduction.py -q 1nd -n 100 -s generateEventsFromCfi.sh -o "-o /store/cmst3/group/hgcal/CMSSW/MinBias_SLHC16 -c UserCode/HGCanalysis/python/minBias_cfi.py";

### Other processes

Can use the minimum bias example, just substitute the argument passed in the -c option to point to the new cfi snippet.


## Producing analysis ntuples

The ntuples are produced by plugins/HGCSimHitsAnalyzer.cc.
Change the code according to your needs.
To submit the production of the ntuples you can use the following script

CHECK ME!

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