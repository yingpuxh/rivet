. /cvmfs/sft.cern.ch/lcg/views/LCG_latest/x86_64-slc6-gcc62-opt/setup.sh
export HADOOP_CONF_DIR=/etc/hadoop/conf
chmod +x rivet-bootstrap
INSTALL_PREFIX=`pwd` MAKE="make -j4" ./rivet-bootstrap
cp rivetenv.sh.local rivetenv.sh
source rivetenv.sh
cd addons
source ./install_addons.sh
source ./make.sh
source ./run.sh
cd ..
