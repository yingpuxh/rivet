#!bin/bash
for FILE in FatHiggsTagging ; do 
    g++ -std=c++11 -o "Rivet${FILE}.so" -shared -fPIC -I${prefix}/include -I${prefix}/include -I${prefix}/include -I${prefix}/include -I${prefix}/include \
	-pedantic -Wall -Wno-long-long -Wno-format -Werror=uninitialized -Werror=delete-non-virtual-dtor -O2 -Wl,--no-as-needed \
	-L${prefix}/lib -L${prefix}/lib -L${prefix}/lib -Wl,-rpath,${prefix}/lib -lm \
	-L${prefix}/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone -lfastjetcontribfragile  ${FILE}.cc -lRivet
done
