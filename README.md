Simple Geant4 Simulation for Solenoid setup of SuperFRS

Needed : 
* cmake (> 3.1)
* Geant4 (last version)
* ROOT6 (last version) or ROOT5 (> 5.34/34)
* boost (> 1.54.0)
Build : (classic cmake build)

$ mkdir build
$ cd build
$ cmake ..

for debug flags : the following have to be added to cmake command :
-DCMAKE_BUILD_TYPE=Debug

and for release (compiled with optimization flags) :
-DCMAKE_BUILD_TYPE=Release

Run option :
binary name : ./G4SolenoidSimple
* [-h] [--help] : provide command line help.
* [-g] [--gui] : set flags for graphical interface (if Geant4 is built with QT).
* [-i inputfile] [--input inputfile] : provide the inputfile i.e. the input file from external event generators. Be aware that this need the PrimaryGeneratorAction for your input file to be provided, included and compiled into the Geant4 code.
* [-m run.mac] [--mac run.mac] : provide Geant4 macro file for Geant4 settings (for example number of events generated).
* [-c config.par] [--config config.par] : configuration file for the code, can set different option that are implemented within the simulation code.
* Outputfile.root : mandatory ! provide name of ROOT output file. 