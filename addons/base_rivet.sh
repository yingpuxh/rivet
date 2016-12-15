#!/bin/bash

echo -e "\nJob started at "`date`" on "`hostname --fqdn`"\n"

#setup run environment 
source SUBRIVETDIR/rivetenv.sh 

# change to scratch directory
cd ${TMPDIR}

# check on how many files we need to run 
NFILES=`ls SUBINPUTDIR/SUBINPUTFILE/*.hepmc.gz | wc -l`
NREST=`echo "scale=0; ${NFILES}%SUBNFILESPERJOB " | bc`
NJOBS=`echo "scale=0; (${NFILES}-${NREST})/SUBNFILESPERJOB + 1" | bc`
# set job id 
JOBID="$(printf '%05d' "${SGE_TASK_ID}")" 
if [ ${SGE_TASK_ID} == ${NJOBS} ]; then
    HEAD=`echo "scale=0; SUBNFILESPERJOB*(${SGE_TASK_ID}-1) + ${NREST}" | bc` 
else
    HEAD=`echo "scale=0; SUBNFILESPERJOB*${SGE_TASK_ID}" | bc` 
    NREST=SUBNFILESPERJOB
fi 
echo "run rivet"
mkfifo ${TMPDIR}/fifo.hepmc
less `ls  SUBINPUTDIR/SUBINPUTFILE/*.hepmc.gz  | head -${HEAD} | tail -${NREST}` >> fifo.hepmc & 
rivet -a SUBANALYSIS --process SUBPROCESS -o output.yoda fifo.hepmc

echo "copy output"
cp output.yoda SUBOUTPUTDIR/SUBOUTPUTFILE__${JOBID}.yoda

#done 
echo -e "\nJob stoped at "`date`" on "`hostname --fqdn`"\n"

exit 0
