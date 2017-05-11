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

## Input files :

### Configuration option :

The default option can be seen in src/G4SolConfig.cc 

Example : 
```
Target_Size 2 cm
```

The first string is the key in the configuration data structure. The second can be anything: string, integer, or double; it will be delimited by the white space. The third key is only for declaring the unit. In this example the key into the configuration data structure is "Target_Size", which stand for the half-size of target, which is set to be 2 cm. The string "Target\_Size" is use to retrieve the contain "2 cm". The unit is the unit defined within Geant4 framework. 

Any key / contain can be set in the configuration file. Please keep the convention that the key name starts with a Capital letter. It will loaded in the configuration data structure. Then it can be obtained via the config object. Mandatory key/contain do not need to be set in the configuration file since they will be set and fill by default values. If other values are needed, the default values are overridden by the values of the configuration file.

Any line in the configuration file can be made ignored with a # at the beginning of the line as a shell comment

Example of configuration file:

```
Target_Size 2 cm
#CDS_RelativePosTarget -1.5
#CD_AbsPosTarget 1. m
Wasa_Side 0
Wasa_ShiftZ 1. m
#Systematic_Shift 0 cm
Field_CDS_Bz 1.8 tesla
chambersize 300 mm
#Particle lambda
Particle H3L
Beam_Momentum 9 GeV
Geo Wasa
Geometry_Namefile ../WasaGeoRoot/WasaGeometry5_biggerEntrance.root 
SimpleGeo 0
#Physicslist G4_QGSP_FTFP_BERT
Physicslist G4Default_FTFP_BERT
HypHI_InnerTracker_Spec 1
HypHI_InnerTracker_PosZ 2. cm
HypHI_InnerTracker_Nb 4
HypHI_InnerTracker_Spacing 1. cm
HypHI_Si_minR 0. cm
HypHI_Si_maxR 15. cm
HypHI_EndCap_posZ 150. cm
HypHI_EndCap_maxR 20. cm
FRS_FMF2_posZ 3. m
HypHI_TR1_posZ 20. cm
HypHI_TR2_posZ 40. cm
#ReductionFactor 0.8
ConvertRoot GeoWasaRealSolenoidFront.root
#HypHI_Downstream_Size 0.85 m
#HypHI_Downstream_Shift -0.36 m
#HypHI_InnerTrackerBox_Visible 0
HyperNuclei_H3L_T12 0.200 ns
```

Some explanations for the most important options:

| Key name                    | possible value            | purpose |
|-----------------------------|---------------------------|-------|
| Particle                    | string                    | Name of the particle/ion for the beam definition. It is used when then Geant4 event generator is used. for ion : Z3A6 = Li6 beam, Z6A12 = C12 | 
| Beam\_SpotSizeSigma         | double + unit: 1 cm       | Size of the beam profile in the transverse plane. Primary vertex is randomize within a Gaussian dist. | 
| Beam\_Momentum              | double + unit : 10 GeV    | Beam momentum in GeV. Be careful in case of ion it is the total momentum, not per nucleon |
| Beam\_MomentumSigma         | double + unit : 1 MeV     | In case of needs of momentum randomization |
| Beam\_Momentum(X, Y, Z)     | double                    | Define the unit direction vector |
| Beam_Momentum(X, Y, Z)sigma | double                    | In case of needs of randomization of the direction vector |
| Target\_Size                | double + unit : 2 cm      | Half-size of the target (cube) |
| Field\_CDS\_Bz              | double + unit : 1.3 tesla | Set the magnetic field value |
|                             |                           | |
| SimpleGeo                   | bool                      | Set for simple geometry case for debugging only |
| Geo                         | string                    | Name of the Geometry class used for the Geant4 geometry construction. See src/G4SolGeometryController.cc |
| Wasa\_ShiftZ                | double + unit: 1 m        | Placement of the Wasa central detector from the target | 
| Wasa\_Side                  | int/bool: 0 or 1          | To flip the Wasa geometry so that the back of the central detector is now the front |
| Geometry\_Namefile          | string                    | Path of the Wasa geometry file |
| Physicslist                 | string                    | Name of the physics list option |
|                             |                           |  |
| DefaultRegionCut            | double + unit: 1 mm       | Value for region cut as defined in Geant4 | 

### Geant4 macro file:

Those files are responsible for the internal Geant4 configuration. The most important part is the /run/beamOn "Nb\_of\_Events"

Example :
```
# Macro file for G4SolSimple
# 
# Can be run in batch, without graphic
# or interactively: Idle> /control/execute run1.mac
#
# Change the default number of workers (in multi-threading mode) 
#/run/numberOfWorkers 4
#
# Initialize kernel
/run/initialize
#
#/run/verbose 1
/run/printProgress 1000
/run/beamOn 50000
```
