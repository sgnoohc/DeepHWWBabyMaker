export SCRAM_ARCH=slc6_amd64_gcc530   # or whatever scram_arch you need for your desired CMSSW release
export CMSSW_VERSION=CMSSW_9_2_0
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
eval `scramv1 runtime -sh`
cd - > /dev/null

export COREDIR=$PWD/CORE
export LD_LIBRARY_PATH=$PWD/rooutil:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PWD/coreutil:$LD_LIBRARY_PATH
