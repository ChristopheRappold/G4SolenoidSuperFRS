# Simple Geant4 Simulation for Solenoid setup of SuperFRS

Geant4 simulation for the solenoid-type setup at FRS and SuperFRS.

## Installation

OS X & Linux:

1.  first:

```sh
git clone git@gitlab.com:TDR-Solenoid/G4SolenoidSuperFRS.git
```

2.  The following environment variables are needed for the configuration via cmake to work: 

```sh
echo $ROOTSYS
echo $VGM_INSTALL
echo $GEANT4_DIR
```

3.  Build:
```sh 
mkdir build
cd build
cmake ..
```

4.  To create and save the Wasa geometry in a rootfile:

Inside ROOT:

```sh
.L geometry.C++
```
Then there are several options:
geometry(bool Central\_FW = false, bool Central\_Pipe = false, bool EMC = false, bool ForwardCal = false)
```sh
geometry()
```

To save the geometry:
```sh
gGeoManager->Export("name_of_your_geometry.root");
```

5. To run the simulation:

This is not needed anymore:
>>>
Go to your build directory. You would need to make symbolic link of the geometry :

cd build

ln -s ../WasaGeoRoot/"name...".root ./geometry.root

Later the geometry will be load from the configuration file, instead from a hardlink in the code.
>>>

The path for the geometry file is now set via the config file, it is not hardcoded anymore

## Needed : 

* cmake (> 3.1)
* Geant4 (last version)
* ROOT6 (last version) or ROOT5 (> 5.34/34)
* boost (> 1.54.0)
* VGM (> 4.2) 

## Build : (classic cmake build)

for debug flags : the following have to be added to cmake command :
-DCMAKE_BUILD_TYPE=Debug

and for release (compiled with optimization flags) :
-DCMAKE_BUILD_TYPE=Release

## Run option :
binary name : ./G4SolenoidSimple
* [-h] [--help] : provide command line help.
* [-g] [--gui] : set flags for graphical interface (if Geant4 is built with QT).
* [-i inputfile] [--input inputfile] : provide the inputfile i.e. the input file from external event generators. Be aware that this need the PrimaryGeneratorAction for your input file to be provided, included and compiled into the Geant4 code.
* [-m run.mac] [--mac run.mac] : provide Geant4 macro file for Geant4 settings (for example number of events generated).
* [-c config.par] [--config config.par] : configuration file for the code, can set different option that are implemented within the simulation code.
* Outputfile.root : mandatory ! provide name of ROOT output file. 

## Example :

*  Run with GUI with beam and geometry setting in config file. An output file is mandatory, in this case it is just like a dummy file since it will have only the event run in the GUI.

```sh
./G4SolenoidSimple -g -c configWasa_1.8T.par OutGUI.root
```

* Run in batch mode for full run. Config file and macro file are provided. External event generator file is set as input. And finally the output file. 

```sh
./G4SolenoidSimple -c configWasa_1.8T.par -i ../../UrQMDEvents_10_L_wAssoc.dat_div60_new.dat -m runSim.mac run_LambdaWithA_WasaFront3_Field1.8T.root
```



## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Merge Request
