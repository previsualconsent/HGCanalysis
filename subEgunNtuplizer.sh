#!/bin/bash
a=(SingleElectron_SLHC13_50um SingleElectron_SLHC13_noB_50um SingleElectron_vSLHC13 SingleElectron_vSLHC13_noB)
for k in ${a[@]}; do for i in `seq 0 50 450`; do python scripts/submitLocalAnalysis_cfg.py -t ${k} -q 8nh -f ${i} -s 50; done done