#!bin/bash 

VERSION=2.5.2

cp ../Rivet-${VERSION}/rivetenv.*h ../.
source ../rivetenv.sh # need to get prefix properly filled 
if [ -f /afs/cern.ch/sw/lcg/contrib/gcc/5.2/x86_64-slc6-gcc52-opt/setup.sh ]; then 
    echo "source /afs/cern.ch/sw/lcg/contrib/gcc/5.2/x86_64-slc6-gcc52-opt/setup.sh" >> ../rivetenv.sh
fi
if [ -f /etc/profile.d/modules.sh ] ; then 
    echo "source /etc/profile.d/modules.sh" >> ../rivetenv.sh
    echo "module load root/5.34" >> ../rivetenv.sh
    echo "module load python/2.7" >> ../rivetenv.sh
fi

echo "export RIVET_ANALYSIS_PATH=${prefix}/addons" >> ../rivetenv.sh
echo "export RIVET_REF_PATH=${prefix}/addons" >> ../rivetenv.sh

source ../rivetenv.sh

for FILE in FastHiggsJets FastTopJets PartonicHiggs ; do 
    if [ -f Projections/${FILE}.cc ] ; then 	      
	cp -p Projections/${FILE}.cc  ${prefix}/Rivet-${VERSION}/src/Projections/.
    fi 
    if [ -f Projections/${FILE}.hh ] ; then	
	cp -p Projections/${FILE}.hh  ${prefix}/Rivet-${VERSION}/include/Rivet/Projections/.    
    fi
done    
cp Projections/Makefile.in  ${prefix}/Rivet-${VERSION}/src/Projections/. 
cp Rivet/Makefile.in  ${prefix}/Rivet-${VERSION}/include/Rivet/Makefile.in
cp Tools/ParticleIdUtils.hh ${prefix}/Rivet-${VERSION}/include/Rivet/Tools/ParticleIdUtils.hh

cd ${prefix}/Rivet-${VERSION}
make -j8 
make install
cd ${prefix}/addons 

cp rivet ${prefix}/bin/.
