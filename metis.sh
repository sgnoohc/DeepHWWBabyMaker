#!/bin/bash

PACKAGE=package.tar.gz


###################################################################################################
# ProjectMetis/CondorTask specific (Setting up some common environment)
###################################################################################################
#echo "To check whether it's on condor universe Vanilla or Local. The following variables are used."
#echo "If _CONDOR_SLOT is set, it's on Vanilla"
#echo "If X509_USER_PROXY is set, it's either on Vanilla or Local."
# if 
if [ "x${_CONDOR_JOB_AD}" == "x" ]; then
    :
    INPUTFILENAMES=$1
    shift;
else
    hostname
    uname -a
    date
    whoami
    pwd
    echo "ls'ing hadoop"
    ls -lh /hadoop/cms/store/user/phchang/
    echo "_CONDOR_SLOT" ${_CONDOR_SLOT}
    echo "X509_USER_PROXY" ${X509_USER_PROXY}
    echo "_CONDOR_SCRATCH_DIR"             ${_CONDOR_SCRATCH_DIR}
    echo "_CONDOR_SLOT"                    ${_CONDOR_SLOT}
    echo "CONDOR_VM"                       ${CONDOR_VM}
    echo "X509_USER_PROXY"                 ${X509_USER_PROXY}
    echo "_CONDOR_JOB_AD"                  ${_CONDOR_JOB_AD}
    echo "_CONDOR_MACHINE_AD"              ${_CONDOR_MACHINE_AD}
    echo "_CONDOR_JOB_IWD"                 ${_CONDOR_JOB_IWD}
    echo "_CONDOR_WRAPPER_ERROR_FILE"      ${_CONDOR_WRAPPER_ERROR_FILE}
    echo "CONDOR_IDS"                      ${CONDOR_IDS}
    echo "CONDOR_ID"                       ${CONDOR_ID}
    OUTPUTDIR=$1
    OUTPUTNAME=$2
    INPUTFILENAMES=$3
    IFILE=$4
    CMSSWVERSION=$5
    SCRAMARCH=$6
    if [ "x${_CONDOR_SLOT}" == "x" ]; then
        WORKDIR=/tmp/phchang_condor_local_${OUTPUTDIR//\//_}_${OUTPUTNAME}_${IFILE}
        mkdir -p ${WORKDIR}
        ls
        cp package.tar.gz ${WORKDIR}
        cd ${WORKDIR}
        ls
        echo "This is in Condor session with Universe=Local."
        echo "WORKDIR=${WORKDIR}"
        echo "pwd"
        pwd
    fi
    echo "OUTPUTDIR     : $1"
    echo "OUTPUTNAME    : $2"
    echo "INPUTFILENAMES: $3"
    echo "IFILE         : $4"
    echo "CMSSWVERSION  : $5"
    echo "SCRAMARCH     : $6"
    shift 6
    tar xvzf package.tar.gz
    if [ $? -eq 0 ]; then
        echo "Successfully untarred package."
        :
    else
        echo "Failed to untar package."
        exit
    fi
fi
###################################################################################################

echo ">>> source setup.sh"
source setup.sh
echo ">>> ls -l"
ls -l
echo ">>> export COREDIR=$PWD/CORE/"
export COREDIR=$PWD/CORE/
echo ">>> ./processBaby ${INPUTFILENAMES} -1"
INPUTFILES=$(echo ${INPUTFILENAMES} | tr ',' ' ')
INDEX=1
for file in $INPUTFILES; do
    ./processBaby "$file" -1 ${INDEX} &
    JOBS[${INDEX}]=$!
    INDEX=$((INDEX + 1))
done
GOODOUTPUTS=""
INDEX=1
for job in "${JOBS[@]}"; do
    wait $job
    JOBSTATUS=$?
    echo ""
    echo "Job #${INDEX} Exit Code = $JOBSTATUS"
    echo ""
    if [ $JOBSTATUS -eq 0 ]; then
        echo ""
        echo "Job #${INDEX} Success"
        echo ""
        GOODOUTPUTS="${GOODOUTPUTS} output_${INDEX}.root"
    else
        echo ""
        echo "Job #${INDEX} Failed"
        echo ""
    fi
    INDEX=$((INDEX+1))
done
hadd -f output.root ${GOODOUTPUTS}

###################################################################################################
# ProjectMetis/CondorTask specific (Copying files over to hadoop)
###################################################################################################
if [ "x${_CONDOR_JOB_AD}" == "x" ]; then
    :
else
    echo "==============================================================================="
    echo " Copying files to output directory"
    echo "==============================================================================="
    hostname
    uname -a
    date
    whoami
    pwd
    echo "ls'ing hadoop"
    ls -lh /hadoop/cms/store/user/phchang/
    if [[ ${OUTPUTDIR} == *"home/users/"* ]]; then
        mkdir -p ${OUTPUTDIR}
        INFILE=${OUTPUTROOTNAME}
        cp ${INFILE} ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root
    else
        if [ "x${X509_USER_PROXY}" == "x" ]; then
            echo "Copying outputs to Hadoop via cp."
            mkdir -p ${OUTPUTDIR}
            INFILE=${OUTPUTROOTNAME//.root/}
            INDEX=0
            #for OUTPUTFILE in $(ls ${INFILE}*.root); do
            for OUTPUTFILE in $(ls output.root); do
                if [ $INDEX -lt 1 ]; then
                    echo "cp ${OUTPUTFILE} ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root"
                    cp ${OUTPUTFILE} ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root
                else
                    echo "cp ${OUTPUTFILE} ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}_${INDEX}.root"
                    cp ${OUTPUTFILE} ${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}_${INDEX}.root
                fi
                INDEX=$((INDEX+1))
            done
        else
            echo 'ls -l'
            ls -l
            echo 'gfal-copy'
            INFILE=${OUTPUTROOTNAME//.root/}
            INDEX=0
            #for OUTPUTFILE in $(ls ${INFILE}*.root); do
            for OUTPUTFILE in $(ls output.root); do
                if [ $INDEX -lt 1 ]; then
                    echo gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUTFILE} gsiftp://gftp.t2.ucsd.edu/${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root --checksum ADLER32
                    gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUTFILE} gsiftp://gftp.t2.ucsd.edu/${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}.root --checksum ADLER32
                else
                    echo gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUTFILE} gsiftp://gftp.t2.ucsd.edu/${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}_${INDEX}.root --checksum ADLER32
                    gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUTFILE} gsiftp://gftp.t2.ucsd.edu/${OUTPUTDIR}/${OUTPUTNAME}_${IFILE}_${INDEX}.root --checksum ADLER32
                fi
                INDEX=$((INDEX+1))
            done
        fi
    fi
    if [ $? -eq 0 ]; then
        echo "Job Success"
    else
        echo "Job Failed"
    fi
    date
fi
###################################################################################################

