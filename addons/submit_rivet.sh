#!/bin/bash

### please adjust the following settings
RESUBMIT=0
NFILESPERJOB=10 
RIVETDIR=/nfs/dust/cms/user/agrohsje/exsoft/rivet 
INPUTDIR=/nfs/dust/cms/user/duerrm/data/MonoDarkHiggs/backgrounds
OUTPUTDIR=/nfs/dust/cms/user/agrohsje/samples/all/mg5_amcatnlo/yoda/${PROCESS}_${JOBTAG}
INPUTFILES=( ttbar_NLO/Wplv  ttbar_NLO/Wmlv  WmZ_NLO/Events  WpZ_NLO/Events ZZ_NLO/Events Zbb_NLO/Events Wbb )
OUTPUTFILES=(    ttbar_Wplv      ttbar_Wmlv             WmZ             WpZ            ZZ            Zbb Wbb )
PROCESSES=(           ttbar           ttbar              WZ              WZ            ZZ            Zbb Wbb )
JOBTAG=validateatlas_161208 
ANALYSIS=FatHiggsTagging


### start job submission 
JOBID=`echo "scale=0; ${#INPUTFILES[@]} -1 " | bc` 
for IJOBID in `seq 0 ${JOBID}`; do
    PROCESS=${PROCESSES[${IJOBID}]}
    INPUTFILE=${INPUTFILES[${IJOBID}]}
    OUTPUTFILE=${OUTPUTFILES[${IJOBID}]}
    # job splitting according to number of files and files per job 
    NFILES=`ls ${INPUTDIR}/${INPUTFILE}/*.hepmc.gz | wc -l`
    NREST=`echo "scale=0; ${NFILES}%${NFILESPERJOB} " | bc`
    if [ ${NREST} == 0 ] ; then
	NJOBS=`echo "scale=0; (${NFILES}-${NREST})/${NFILESPERJOB}" | bc`
    else
	NJOBS=`echo "scale=0; (${NFILES}-${NREST})/${NFILESPERJOB} + 1" | bc`
    fi
    ### resubmit jobs that failed 
    if [ ${RESUBMIT} == "1" ] ; then
	for IJOB in `seq 1 ${NJOBS}`; do
	    ID="$(printf '%05d' "${IJOB}")"
	    if [ ! -f ${OUTPUTDIR}/${PROCESS}_${JOBTAG}/${OUTPUTFILE}__${ID}.yoda ] ; then
		echo "missing file ${OUTPUTDIR}/${PROCESS}_${JOBTAG}/${OUTPUTFILE}__${ID}.yoda"
                qsub -N ${OUTPUTFILE}_${JOBTAG} -V -j y -m as -o ${OUTPUTDIR}/${PROCESS}_${JOBTAG} \
		    -l h_rt=24:00:00 -l h_vmem=8G -l distro=sld6 \
		    -t ${IJOB}:${IJOB}  ${OUTPUTDIR}/${PROCESS}_${JOBTAG}/scripts/${OUTPUTFILE}_${JOBTAG}.sh
            fi
        done
    ### initial job submission 
    else
        mkdir -p ${OUTPUTDIR}/${PROCESS}_${JOBTAG}/scripts
	### create job bash file 
	cp base_rivet.sh ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBINPUTDIR|${INPUTDIR}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBINPUTFILE|${INPUTFILE}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBOUTPUTDIR|${OUTPUTDIR}/${PROCESS}_${JOBTAG}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBOUTPUTFILE|${OUTPUTFILE}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBRIVETDIR|${RIVETDIR}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBANALYSIS|${ANALYSIS}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBPROCESS|${PROCESS}|g" ${OUTPUTFILE}_${JOBTAG}.sh 
	sed -i -e "s|SUBNFILESPERJOB|${NFILESPERJOB}|g" ${OUTPUTFILE}_${JOBTAG}.sh 

	chmod +x ${OUTPUTFILE}_${JOBTAG}.sh
	mv ${OUTPUTFILE}_${JOBTAG}.sh ${OUTPUTDIR}/${PROCESS}_${JOBTAG}/scripts/.

	qsub -N ${OUTPUTFILE}_${JOBTAG} -V -j y -m as -o ${OUTPUTDIR}/${PROCESS}_${JOBTAG} \
	    -l h_rt=12:00:00 -l h_vmem=8G  -l distro=sld6 \
	    -t 1:${NJOBS} ${OUTPUTDIR}/${PROCESS}_${JOBTAG}/scripts/${OUTPUTFILE}_${JOBTAG}.sh 
	sleep 2s
    fi
done
