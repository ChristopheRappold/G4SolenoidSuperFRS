Simple Geant4 Simulation for Solenoid setup of SuperFRS

Needed : 
* cmake (> 3.1)
* Geant4 (last version)
* ROOT6 (last version)

Build : (classic cmake build)

$ mkdir build
$ cd build
$ cmake ..

for debug flags : the following have to be added to cmake command :
-DCMAKE_BUILD_TYPE=Debug

and for release (compiled with optimization flags) :
-DCMAKE_BUILD_TYPE=Release