#!/bin/env python

import time
import json

from metis.Sample import DirectorySample
from metis.CondorTask import CondorTask

from metis.StatsParser import StatsParser

import metis

import sys
import os

# Configurations
main_dir = os.path.dirname(os.path.abspath(__file__))
metis_path = os.path.dirname(os.path.dirname(metis.__file__))
job_tag = "2016_tnp_sample_v1.5"
exec_path = os.path.join(main_dir, "metis.sh")
tar_path = os.path.join(main_dir, "package.tar.gz")
hadoop_path = "metis/lepbaby/{}".format(job_tag)
args = ""

# Get into the directory where this lepmetis.py sits. So we can tar up the condor package.
os.system("tar -czf package.tar.gz setup.sh hadoopmap.txt coreutil/data rooutil/*.so coreutil/*.so libs/ processBaby CORE/Tools")

dslocs = [
    ["/DY"                                            , "/hadoop/cms/store/group/snt/run2_moriond17/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/V08-00-16/" ] ,
#    ["/TTbar"                                         , "/hadoop/cms/store/group/snt/run2_moriond17/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/"                      ] ,
#    ["/WJets"                                         , "/hadoop/cms/store/group/snt/run2_moriond17/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/"           ] ,
#    ["/WW"                                            , "/hadoop/cms/store/group/snt/run2_moriond17/WW_TuneCUETP8M1_13TeV-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/"                               ] ,
#    ["/WZ"                                            , "/hadoop/cms/store/group/snt/run2_moriond17/WZ_TuneCUETP8M1_13TeV-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/"                               ] ,
#    ["/ZZ"                                            , "/hadoop/cms/store/group/snt/run2_moriond17/EWK_LNuJJ_MJJ-120_13TeV-madgraph-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/"                    ] ,
    ["/Run2016B_SingleElectron_MINIAOD_03Feb2017_ver2-v2" , "/nfs-7/userdata/dataTuple/kludge/Run2016B_SingleElectron_MINIAOD_03Feb2017_ver2-v2/merged/V08-00-18/"                                                                         ] ,
    ["/Run2016C_SingleElectron_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016C_SingleElectron_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                              ] ,
    ["/Run2016D_SingleElectron_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016D_SingleElectron_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                              ] ,
    ["/Run2016E_SingleElectron_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016E_SingleElectron_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                              ] ,
    ["/Run2016F_SingleElectron_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016F_SingleElectron_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                              ] ,
    ["/Run2016G_SingleElectron_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016G_SingleElectron_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                              ] ,
    ["/Run2016H_SingleElectron_MINIAOD_03Feb2017_ver2-v1" , "/nfs-7/userdata/dataTuple/kludge/Run2016H_SingleElectron_MINIAOD_03Feb2017_ver2-v1/merged/V08-00-18/"                                                                         ] ,
    ["/Run2016H_SingleElectron_MINIAOD_03Feb2017_ver3-v1" , "/nfs-7/userdata/dataTuple/kludge/Run2016H_SingleElectron_MINIAOD_03Feb2017_ver3-v1/merged/V08-00-18/"                                                                         ] ,
    ["/Run2016B_SingleMuon_MINIAOD_03Feb2017_ver2-v2"     , "/nfs-7/userdata/dataTuple/kludge/Run2016B_SingleMuon_MINIAOD_03Feb2017_ver2-v2/merged/V08-00-18/"                                                                             ] ,
    ["/Run2016C_SingleMuon_MINIAOD_03Feb2017-v1"          , "/nfs-7/userdata/dataTuple/kludge/Run2016C_SingleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                  ] ,
    ["/Run2016D_SingleMuon_MINIAOD_03Feb2017-v1"          , "/nfs-7/userdata/dataTuple/kludge/Run2016D_SingleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                  ] ,
    ["/Run2016E_SingleMuon_MINIAOD_03Feb2017-v1"          , "/nfs-7/userdata/dataTuple/kludge/Run2016E_SingleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                  ] ,
    ["/Run2016F_SingleMuon_MINIAOD_03Feb2017-v1"          , "/nfs-7/userdata/dataTuple/kludge/Run2016F_SingleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                  ] ,
    ["/Run2016G_SingleMuon_MINIAOD_03Feb2017-v1"          , "/nfs-7/userdata/dataTuple/kludge/Run2016G_SingleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                  ] ,
    ["/Run2016H_SingleMuon_MINIAOD_03Feb2017_ver2-v1"     , "/nfs-7/userdata/dataTuple/kludge/Run2016H_SingleMuon_MINIAOD_03Feb2017_ver2-v1/merged/V08-00-18/"                                                                             ] ,
    ["/Run2016H_SingleMuon_MINIAOD_03Feb2017_ver3-v1"     , "/nfs-7/userdata/dataTuple/kludge/Run2016H_SingleMuon_MINIAOD_03Feb2017_ver3-v1/merged/V08-00-18/"                                                                             ] ,
#    ["/QCD_Pt_15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"         , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]         , 
#    ["/QCD_Pt_20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"         , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]         , 
#    ["/QCD_Pt_30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"         , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]         , 
#    ["/QCD_Pt_50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"         , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]         , 
#    ["/QCD_Pt_80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v3"   , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/V08-00-16/" ]   , 
#    ["/QCD_Pt_120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"       , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]       , 
#    ["/QCD_Pt_170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1"  , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/V08-00-16/" ]  , 
#    ["/QCD_Pt_300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2_v1"  , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/V08-00-16/" ]  , 
#    ["/QCD_Pt_470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1"  , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/V08-00-16/" ]  , 
#    ["/QCD_Pt_600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"       , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]       , 
#    ["/QCD_Pt_800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"      , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]      , 
#    ["/QCD_Pt_1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v3" , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/V08-00-16/" ] , 
#    ["/QCD_Pt_20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"            , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]            , 
#    ["/QCD_Pt_30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1"       , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/V08-00-16/" ]       , 
#    ["/QCD_Pt_50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1"       , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/V08-00-16/" ]       , 
#    ["/QCD_Pt_80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1"      , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/V08-00-16/" ]      , 
#    ["/QCD_Pt_120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1_v1"     , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/V08-00-16/" ]     , 
#    ["/QCD_Pt_170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"          , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]          , 
#    ["/QCD_Pt_300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"          , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]          , 
#    ["/QCD_Pt_15to20_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"                 , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt_15to20_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]                 , 
#    ["/QCD_Pt_20to30_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"                 , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt_20to30_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]                 , 
#    ["/QCD_Pt_30to80_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"                 , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt_30to80_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]                 , 
#    ["/QCD_Pt_80to170_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"         , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt_80to170_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]         , 
#    ["/QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"               , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]               , 
#    ["/QCD_Pt_250toInf_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1"               , "/hadoop/cms/store/group/snt/run2_moriond17/QCD_Pt_250toInf_bcToE_TuneCUETP8M1_13TeV_pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/" ]               , 
#    ["/Run2016B_DoubleEG_MINIAOD_03Feb2017_ver2-v2"   , "/nfs-7/userdata/dataTuple/kludge/Run2016B_DoubleEG_MINIAOD_03Feb2017_ver2-v2/merged/V08-00-18/"                                                                                                  ] ,
#    ["/Run2016C_DoubleEG_MINIAOD_03Feb2017-v1"        , "/nfs-7/userdata/dataTuple/kludge/Run2016C_DoubleEG_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                       ] ,
#    ["/Run2016D_DoubleEG_MINIAOD_03Feb2017-v1"        , "/nfs-7/userdata/dataTuple/kludge/Run2016D_DoubleEG_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                       ] ,
#    ["/Run2016E_DoubleEG_MINIAOD_03Feb2017-v1"        , "/nfs-7/userdata/dataTuple/kludge/Run2016E_DoubleEG_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                       ] ,
#    ["/Run2016F_DoubleEG_MINIAOD_03Feb2017-v1"        , "/nfs-7/userdata/dataTuple/kludge/Run2016F_DoubleEG_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                       ] ,
#    ["/Run2016G_DoubleEG_MINIAOD_03Feb2017-v1"        , "/nfs-7/userdata/dataTuple/kludge/Run2016G_DoubleEG_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                       ] ,
#    ["/Run2016H_DoubleEG_MINIAOD_03Feb2017_ver2-v1"   , "/nfs-7/userdata/dataTuple/kludge/Run2016H_DoubleEG_MINIAOD_03Feb2017_ver2-v1/merged/V08-00-18/"                                                                                                  ] ,
#    ["/Run2016H_DoubleEG_MINIAOD_03Feb2017_ver3-v1"   , "/nfs-7/userdata/dataTuple/kludge/Run2016H_DoubleEG_MINIAOD_03Feb2017_ver3-v1/merged/V08-00-18/"                                                                                                  ] ,
#    ["/Run2016B_DoubleMuon_MINIAOD_03Feb2017_ver2-v2" , "/nfs-7/userdata/dataTuple/kludge/Run2016B_DoubleMuon_MINIAOD_03Feb2017_ver2-v2/merged/V08-00-18/"                                                                                                ] ,
#    ["/Run2016C_DoubleMuon_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016C_DoubleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                     ] ,
#    ["/Run2016D_DoubleMuon_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016D_DoubleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                     ] ,
#    ["/Run2016E_DoubleMuon_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016E_DoubleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                     ] ,
#    ["/Run2016F_DoubleMuon_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016F_DoubleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                     ] ,
#    ["/Run2016G_DoubleMuon_MINIAOD_03Feb2017-v1"      , "/nfs-7/userdata/dataTuple/kludge/Run2016G_DoubleMuon_MINIAOD_03Feb2017-v1/merged/V08-00-18/"                                                                                                     ] ,
#    ["/Run2016H_DoubleMuon_MINIAOD_03Feb2017_ver2-v1" , "/nfs-7/userdata/dataTuple/kludge/Run2016H_DoubleMuon_MINIAOD_03Feb2017_ver2-v1/merged/V08-00-18/"                                                                                                ] ,
#    ["/Run2016H_DoubleMuon_MINIAOD_03Feb2017_ver3-v1" , "/nfs-7/userdata/dataTuple/kludge/Run2016H_DoubleMuon_MINIAOD_03Feb2017_ver3-v1/merged/V08-00-18/"                                                                                                ] ,
]

os.chdir(metis_path)

total_summary = {}
while True:
    allcomplete = True
    for ds,loc in dslocs:
        task = CondorTask(
                sample = DirectorySample( dataset=ds, location=loc ),
                open_dataset = False,
                flush = True,
                files_per_output = 2,
                output_name = "merged.root",
                tag = job_tag,
                cmssw_version = "CMSSW_9_2_1", # doesn't do anything
                arguments = args,
                executable = exec_path,
                tarfile = tar_path,
                condor_submit_params = {"sites" : "T2_US_UCSD"},
                special_dir = hadoop_path,
                max_jobs = 10
                )
        task.process()
        allcomplete = allcomplete and task.complete()
        # save some information for the dashboard
        total_summary[ds] = task.get_task_summary()
    # parse the total summary and write out the dashboard
    StatsParser(data=total_summary, webdir="~/public_html/dump/lepmetis/").do()
    os.system("chmod -R 755 ~/public_html/dump/lepmetis")
    if allcomplete:
        print ""
        print "Job={} finished".format(job_tag)
        print ""
        break
    if job_tag == "filescan":
        break
    print "Sleeping 300 seconds ..."
    time.sleep(300)
