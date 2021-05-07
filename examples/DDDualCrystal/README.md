This is a simulation of a dual readout crystal calorimeter (currently the code is work in progress).  See https://iopscience.iop.org/article/10.1088/1748-0221/15/11/P11005 for the concept.

A Very Important Thing to NOTE:

By default, Cerenkov radiation is turned off in DDG4 and cannot easily be turned on by an example.  To turn it on, please replace DDG4/python/DDSim/Helper/Physics.py with the file at: https://github.com/saraheno/DD4hep/blob/master/DDG4/python/DDSim/Helper/Physics.py and recompile.

Some instructions on how to use:

First time:

 git clone git@github.com:saraheno/DD4hep.git

 cd DD4hep

 (choose a line similar to the next one appropriate for your OS)
 source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

 mkdir build

 cd build

 cmake -DDD4HEP_USE_GEANT4=ON -DBoost_NO_BOOST_CMAKE=ON -DDD4HEP_USE_LCIO=ON -DDD4HEP_USE_HEPMC3=ON -DBUILD_TESTING=ON -DDD4HEP_USE_XERCESC=ON -DCMAKE_CXX_STANDARD=17 -DGeant4_DIR=$G4INSTALL/lib/Geant4-10.4.0 -DROOT_DIR=$ROOTSYS -DDD4HEP_BUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Release ..

 make -j 5

 make install

 cd ..

 source bin/thisdd4hep.sh

Subsequent times, in the DD4hep area:

 (choose a line similar to the next one appropriate for your OS)
   source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

  source bin/thisdd4hep.sh


To Visualize Geometry:

  go to DD4hep/examples/DDDualcrystal/SCEPCdetector

  geoDisplay DDDetectors/compact/SCEPCAL.xml
    

to make some interactive events in GEANT

  go to DD4hep/examples/DDDualcrystal/SCEPCdetector

  ddsim --compactFile=./SCEPCAL.xml --runType=vis -G -N=2  --steeringFile SCEPCALsteering.py   --outputFile=testSCEPCAL.root --gun.position "0.0 0.0 1.0*cm" --gun.direction "1.0 1.0 0.1" --gun.energy "30*GeV" --part.userParticleHandler="" --gun.particle "e-" --runType vis

and then in the geant4 shell

  /control/execute vis.mac  (save this in my root directory)

  /run/beamOn 1

  /vis/viewer/refresh

  /vis/viewer/zoomTo 10

  /vis/viewer/pan -100 200 cm



to make a few debugging events in GEANT

  go to DD4hep/examples/DDDualcrystal/SCEPCdetector

  ddsim --compactFile=./SCEPCAL.xml --runType=batch -G -N=2  --steeringFile SCEPCALsteering.py --outputFile=testSCEPCAL.root --gun.position "0.0 0.0 1.0*cm" --gun.direction "0.0 1.0 0.1" --gun.energy "30*GeV" --part.userParticleHandler=""   --gun.particle "e-"


To make some simple plots on the GEANT output

  go to DD4hep/examples/DDDualcrystal/SCEPCdetector

  root.exe

    root .L SCEPCALplot.C

    SCEPCAL(10,"testSCEPCAL.root")


  then look at the resulting histograms in hist.root

  




