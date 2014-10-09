#!/bin/bash


#
# PARSE CONFIGURATION PARAMETERS
#
BASEDIR=/store/cmst3/group/hgcal/CMSSW
STOREDIR=/store/cmst3/group/hgcal/CMSSW/Ntuples
TAG=Single11_SLHC16
IFILE=0
NFILES=10
WORKDIR="/tmp/`whoami`/"

while getopts "hn:f:i:t:o:w:" opt; do
    case "$opt" in
    h)
        echo ""
        echo "runHGCAnalyzer.sh [OPTIONS]"
	echo "     -n      number of files to process"
	echo "     -f      first file to process"
	echo "     -i      directory with inputs"
	echo "     -t      tag to use"
	echo "     -o      output directory (local or eos)"
	echo "     -w      local work directory (by default /tmp/user)"
	echo "     -h      help"
        echo ""
	exit 0
        ;;
    n)  NFILES=$OPTARG
        ;;
    f)  IFILE=$OPTARG
        ;;
    i)  BASEDIR=$OPTARG
        ;;
    t)  TAG=$OPTARG
        ;;
    o)  STOREDIR=$OPTARG
	;;
    w)  WORKDIR=$OPTARG
	;;
    esac
done

#
# CONFIGURE JOB
#
BASEJOBNAME=HGCHits_${TAG}_${IFILE}_${NFILES}
OUTFILE=${BASEJOBNAME}.root
LOGFILE=${BASEJOBNAME}.log


#
# RUN cmsRun and at the end move the output to the required directory
#
cmsRun ${CMSSW_BASE}/src/UserCode/HGCanalysis/test/runHGCSimHitsAnalyzer_cfg.py ${BASEDIR}/${TAG} ${IFILE} ${NFILES} > ${WORKDIR}/${LOGFILE} 2>&1

#move output to storage directory
if [[ $STOREDIR =~ .*/store/cmst3.* ]]; then
    cmsMkdir ${STOREDIR}
    cmsStage -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
elif [[ $STOREDIR =~ /afs/.* ]]; then
    cmsMkdir ${STOREDIR}
    cp -f ${WORKDIR}/${OUTFILE} ${STOREDIR}/${OUTFILE}
    rm ${WORKDIR}/${OUTFILE}
fi

echo "Local output @ `hostname` stored @ ${WORKDIR}/${OUTFILE} being moved to ${STOREDIR}" 
echo "log file can be found in ${WORKDIR}/${LOGFILE}"
