//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------
// Implementation of the HypHIPrimaryGeneratorAction class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "HypHIPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
// #include "G4LogicalVolumeStore.hh"
// #include "G4LogicalVolume.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HypHIPrimaryGeneratorAction::HypHIPrimaryGeneratorAction(const G4SolConfig& conf)
    : G4VUserPrimaryGeneratorAction(), status(0), fParticleGun(0), fMessenger(0), Par(conf)

{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);

  boost::optional<std::string> tempName = Par.Get<boost::optional<std::string> >("InputFile");
  if(tempName)
    {
      nameInputFile = Par.Get<std::string>("InputFile");

      fPosX = Par.Get<double>("Target_PosX");
      fPosY = Par.Get<double>("Target_PosY");
      fPosZ = Par.Get<double>("Target_PosZ");

      fSpotSizeSigma = Par.Get<double>("Beam_SpotSizeSigma");
      if(Par.IsAvailable("Beam_SpotSizeSigmaX"))
	{
	  SpotElliptical = 1;
	  fSpotSizeSigmaX = Par.Get<double>("Beam_SpotSizeSigmaX");
	}
      if(Par.IsAvailable("Beam_SpotSizeSigmaY"))
	{
	  SpotElliptical = 1;
	  fSpotSizeSigmaY = Par.Get<double>("Beam_SpotSizeSigmaY");
	}

      if(Par.IsAvailable("Target_Size"))
	{
	  fTargetSizeX    = Par.Get<double>("Target_Size");
	  fTargetSizeY    = Par.Get<double>("Target_Size");
	  fTargetSizeZ    = Par.Get<double>("Target_Size");
	}
      if(Par.IsAvailable("Target_SizeX"))
	fTargetSizeX    = Par.Get<double>("Target_SizeX");
      if(Par.IsAvailable("Target_SizeY"))
	fTargetSizeY    = Par.Get<double>("Target_SizeY");
      if(Par.IsAvailable("Target_SizeZ"))
	fTargetSizeZ    = Par.Get<double>("Target_SizeZ");


      ConstParticle = GetParticle("pi-");

      fParticleGun->SetParticlePosition(G4ThreeVector(fPosX, fPosY, fPosZ));
      fParticleGun->SetParticleDefinition(ConstParticle);

      InStream.open(nameInputFile);

      std::string firstline;
      std::getline(InStream, firstline);

      auto it_comment = firstline.find("HYPHI");
      if(it_comment == std::string::npos)
        status = -1;

      std::getline(InStream, firstline);
      std::stringstream stream(firstline);
      stream >> nEvents;

      DefineCommands();
    }
  else
    status = -2;

  // define commands for this class
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HypHIPrimaryGeneratorAction::~HypHIPrimaryGeneratorAction() { delete fParticleGun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HypHIPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  // G4double worldZHalfLength = 0;
  // G4LogicalVolume* worlLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Target");

  G4double pos_x = fPosX;          // Par.Get_Geometry_TargetPosX();
  G4double pos_y = fPosY;          // Par.Get_Geometry_TargetPosY();
  G4double pos_z = fPosZ;          // Par.Get_Geometry_TargetPosZ();
  G4double sigma = fSpotSizeSigma; // Par.Get_Beam_SpotSizeSigma();

  if(fRandomizePrimary[0])
    {
      bool not_acc_position = true;
      G4double ofpos_x, ofpos_y, ofpos_z;

      while(not_acc_position)
        {
	  if(SpotElliptical == 0)
	    {
	      G4double r0;
	      // r0 = s*std::sqrt(G4UniformRand()); // Uniformly distributed

	      r0           = sigma * std::sqrt(-std::log(G4UniformRand())); // Gaus
	      G4double phi = 2.0 * CLHEP::pi * G4UniformRand();
	      ofpos_x      = r0 * cos(phi);
	      ofpos_y      = r0 * sin(phi);
	      ofpos_z      = 2.*fTargetSizeZ * (0.5 - G4UniformRand());
	    }
	  else
	    {

	      G4double r;
	      G4double v1,v2,fac,valG1,valG2;
 	      do
		{
		  v1 = 2.0 * G4UniformRand() - 1.0;
		  v2 = 2.0 * G4UniformRand() - 1.0;
		  r = v1*v1 + v2*v2;
		}
	      while ( r > 1.0 );

	      fac = std::sqrt(-2.0*std::log(r)/r);
	      valG1 = v1*fac;
	      valG2 = v2*fac;

	      ofpos_x = fSpotSizeSigmaX * valG1;
	      ofpos_y = fSpotSizeSigmaY * valG2;
	      ofpos_z = fTargetSizeZ * (0.5 - G4UniformRand());
	    }
          if(std::abs(ofpos_x) <= fTargetSizeX && std::abs(ofpos_y) <= fTargetSizeY &&
             std::abs(ofpos_z) <= fTargetSizeZ) // Par.Get_Geometry_TargetLength()/2.0 &&
                                               // std::abs(ofpos_y)<=Par.Get_Geometry_TargetHeight()/2.0)
            not_acc_position = false;

          //  std::cout<<"r0= "<<r0<<" "<<ofpos_x<<" "<<ofpos_y<<" "<<ofpos_z<<" "<<
          //  sqrt(ofpos_x*ofpos_x+ofpos_y*ofpos_y+ofpos_z*ofpos_z)<<std::endl;
        }

      pos_x += ofpos_x;
      pos_y += ofpos_y;
      pos_z += ofpos_z;
    }

  int Multiplicity = 0;
  //	THypHi_Communicator::n_current_event++;
  // std::cout << "--------------" << THypHi_Communicator::n_current_event << std::endl;
  std::string temp_line;
  bool ret(std::getline(InStream, temp_line));
  if(ret == false)
    {
      std::cout << "!> End of the file is reached!!!\n";
      G4RunManager::GetRunManager()->AbortRun(true);
      return;
    }

  std::stringstream streamM(temp_line);
  streamM >> Multiplicity;

  if(Multiplicity <= 0)
    {
      std::cout << "!> Wrong multiplisity ";
      std::cout << Multiplicity << "\n";
      G4RunManager::GetRunManager()->AbortRun(true);
      return;
    }

  std::string particlename;
  G4double Px, Py, Pz, energy;
  G4ParticleDefinition* particle;
  G4double phi_ofset = 2.0 * CLHEP::pi * (0.5 - G4UniformRand());
  // int pnum=0;

  for(int i = 0; i < Multiplicity; i++)
    {
      std::getline(InStream, temp_line);
      std::stringstream streamP(temp_line);
      streamP >> particlename >> energy >> Px >> Py >> Pz;

      //	    std::cout<<"-I- Particle name in Generator ===>  "<<particlename<<std::endl;

      particle = GetParticle(particlename);
      //
      G4String particleName_temp = particle->GetParticleName();
      //	    std::cout<<"come here "<<particleName_temp<<std::endl;

      bool accept = true;
      //    	    if(particleName_temp!="L")
      //    	      accept = false;

      if(particleName_temp == "gamma")
        accept = false;

      // 	    if(particleName_temp=="sigma0"||
      // 	       particleName_temp=="sigma+"||
      // 	       particleName_temp=="sigma-")
      // 	      continue;

      //
      // 	    if(particleName_temp!="proton")accept=false;
      //  	    if(particleName_temp=="proton")
      //  	      {
      //  		pnum++;
      //  	      }
      //  	    if(pnum>1)accept=false;

      if(accept && energy > 0)
        {
          fParticleGun->SetParticleDefinition(particle);
          // direction
          G4ThreeVector mom(Px, Py, Pz);
          // 		std::cout<<energy<<" "<<Px<<" "<<Py<<" "<<Pz<<std::endl;
          // 		std::cout<<mom.x()<<" "<<mom.y()<<" "<<mom.z()<<std::endl;

          mom.setPhi(mom.phi() + phi_ofset); // open for rotation of repeated mom seed
          fParticleGun->SetParticleMomentumDirection(mom);
          //
          fParticleGun->SetParticleEnergy(energy);
          fParticleGun->SetParticlePosition(G4ThreeVector(pos_x, pos_y, pos_z));

          fParticleGun->GeneratePrimaryVertex(anEvent);
          //
          // Event->BeamNames.push_back(particle->GetParticleName());
          // Event->BeamMasses.push_back(particle->GetPDGMass());
          // Event->BeamCharges.push_back(particle->GetPDGCharge());
          // Event->BeamMomentums_X.push_back(mom.x());//Px);
          // Event->BeamMomentums_Y.push_back(mom.y());//Py);
          // Event->BeamMomentums_Z.push_back(mom.z());//Pz);
        }
    }

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HypHIPrimaryGeneratorAction::DefineCommands()
{
  // Define /G4SolSimple/generator command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/G4SolSimple/generator/", "Primary generator control");
  double fMomentum, fSigmaMomentum, fSigmaAngle;
  // momentum command
  G4GenericMessenger::Command& momentumCmd =
      fMessenger->DeclarePropertyWithUnit("momentum", "GeV", fMomentum, "Mean momentum of primaries.");
  momentumCmd.SetParameterName("p", true);
  momentumCmd.SetRange("p>=0.");
  momentumCmd.SetDefaultValue("1.");

  std::string nameParticle;
  // momentum command
  G4GenericMessenger::Command& ParticleCmd = fMessenger->DeclareProperty("Particle", nameParticle);
  ParticleCmd.SetGuidance("Which particle");
  ParticleCmd.SetParameterName("particle", true);
  ParticleCmd.SetDefaultValue("proton");

  // momentumCmd.SetParameterName("p", true);
  // momentumCmd.SetRange("p>=0.");

  // sigmaMomentum command
  G4GenericMessenger::Command& sigmaMomentumCmd =
      fMessenger->DeclarePropertyWithUnit("sigmaMomentum", "MeV", fSigmaMomentum, "Sigma momentum of primaries.");
  sigmaMomentumCmd.SetParameterName("sp", true);
  sigmaMomentumCmd.SetRange("sp>=0.");
  sigmaMomentumCmd.SetDefaultValue("50.");

  // sigmaAngle command
  G4GenericMessenger::Command& sigmaAngleCmd =
      fMessenger->DeclarePropertyWithUnit("sigmaAngle", "deg", fSigmaAngle, "Sigma angle divergence of primaries.");
  sigmaAngleCmd.SetParameterName("t", true);
  sigmaAngleCmd.SetRange("t>=0.");
  sigmaAngleCmd.SetDefaultValue("2.");

  // randomizePrimary command
  G4GenericMessenger::Command& randomCmd = fMessenger->DeclareProperty("randomizePrimary", fRandomizePrimary[0]);
  G4String guidance                      = "Integer flag for randomizing primary particle types.\n";
  guidance += "In case this flag is false, you can select the primary particle\n";
  guidance += "  with /gun/particle command.";
  randomCmd.SetGuidance(guidance);
  randomCmd.SetParameterName("flg", true);
  randomCmd.SetDefaultValue("1");
}

G4ParticleDefinition* HypHIPrimaryGeneratorAction::GetParticle(const G4String& particleName)
{

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();

  G4IonTable* pTableIon = G4IonTable::GetIonTable();
  pTableIon->CreateAllIon();
  //     G4ParticleDefinition* particle = pTable->FindParticle(particleName);
  //     if(particle) return particle;

  G4ParticleDefinition* particle; // = pTable->FindParticle(particleName);

  G4int Apos = particleName.rfind("A");
  G4int Zpos = particleName.rfind("Z");
  G4int Epos = particleName.rfind("E");
  G4int Lpos = particleName.rfind("L");

  if(Zpos != (G4int)G4String::npos)
    {
      //	std::cout<<"zpos >0 "<<particleName<<std::endl;
      if(particleName == "Z2A4")
        {
          particle = pTable->FindParticle("alpha");
          if(particle)
            return particle;
        }
      else if(particleName == "Z2A3")
        {
          particle = pTable->FindParticle("He3");
          // particle = pTable->FindParticle("alpha");
          if(particle)
            return particle;
        }
      else if(particleName == "Z1A2")
        {
          //	    std::cout<<" come to deuteron "<<std::endl;
          particle = pTable->FindParticle("deuteron");
          if(particle)
            return particle;
        }
      else if(particleName == "Z1A3")
        {
          particle = pTable->FindParticle("triton");
          if(particle)
            return particle;
        }
      else
        {
          // particle = pTable->FindParticle("gamma");

          //
          if(Apos == (G4int)G4String::npos || Zpos == (G4int)G4String::npos)
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Particle: " << particleName;
              std::cout << " is not found! (Primary Gen)\n";
              std::exit(1);
            }
          if(Apos < Zpos)
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Wrong ion format! should be like ";
              std::cout << " Z1A20E5.11\n";
              std::exit(1);
            }
          //
          G4double E = 0;
          G4int endl = particleName.size();
          if(Epos != (G4int)G4String::npos)
            {
              G4String Estr = particleName.substr(Epos + 1, endl - Epos - 1);
              E             = std::atof(Estr.c_str());
              endl          = Epos;
            }
          //
          G4String Zstr = particleName.substr(Zpos + 1, Apos - Zpos - 1);
          G4String Astr = particleName.substr(Apos + 1, endl - Apos - 1);
          G4int Z       = std::atoi(Zstr.c_str());
          G4int A       = std::atoi(Astr.c_str());
          // printf("A = %i; Z = %i; E = %f\n",A,Z,E);
          //
          particle = pTableIon->GetIon(Z, A, E);
          if(!particle)
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Cannot create the ion ";
              std::cout << particleName << "\n";
              std::exit(1);
            }

          //	    std::cout<<" come to something "<<std::endl;
          if(particle)
            return particle;
        }
    }
  else if(Lpos != (G4int)G4String::npos)
    {
      if(particleName == "Lambda" || particleName == "lambda" || particleName == "Lambda0" || particleName == "L")
        {
          particle = pTable->FindParticle("lambda");
          if(particle)
            return particle;
        }
      else if(particleName == "H3L")
        {
          particle = pTable->FindParticle("H3L");
          // particle = pTableIon->GetIon(1,3,1,0.);
          // particle->DumpTable();
          if(particle)
            return particle;
          else
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Cannot create the ion ";
              std::cout << particleName << "\n";
              std::exit(1);
            }
        }
      else if(particleName == "H4L")
        {
          particle = pTable->FindParticle("H4L");
          // particle = pTableIon->GetIon(1,4,1,0.);
          // particle->DumpTable();
          if(particle)
            return particle;
          else
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Cannot create the ion ";
              std::cout << particleName << "\n";
              std::exit(1);
            }
        }
      else if(particleName == "He4L")
        {
          particle = pTable->FindParticle("He4L");
          // particle = pTableIon->GetIon(1,4,1,0.);
          // particle->DumpTable();
          if(particle)
            return particle;
          else
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Cannot create the ion ";
              std::cout << particleName << "\n";
              std::exit(1);
            }
        }
      else if(particleName == "He5L")
        {
          particle = pTable->FindParticle("He5L");
          // particle = pTableIon->GetIon(1,4,1,0.);
          // particle->DumpTable();
          if(particle)
            return particle;
          else
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Cannot create the ion ";
              std::cout << particleName << "\n";
              std::exit(1);
            }
        }
      else if(particleName == "nnL")
        {
          particle = pTable->FindParticle("nnL");
          // particle = pTableIon->GetIon(1,4,1,0.);
          // particle->DumpTable();
          if(particle)
            return particle;
          else
            {
              std::cout << "!> Fatal error!\n";
              std::cout << "!>    Cannot create the ion ";
              std::cout << particleName << "\n";
              std::exit(1);
            }
        }
      else
        {
          std::cout << "!> Fatal error!\n";
          std::cout << "!>    Cannot create the ion ";
          std::cout << particleName << "\n";
          std::exit(1);
        }
    }
  else
    {
      particle = pTable->FindParticle(particleName);
      if(particle)
        return particle;
    }

  /*
  //
  if(Apos==(G4int)G4String::npos || Zpos==(G4int)G4String::npos){
      std::cout << "!> Fatal error!\n";
      std::cout << "!>    Particle: " << particleName;
      std::cout << " is not found! (Primary Gen)\n";
      std::exit(1);
  }
  if(Apos<Zpos){
      std::cout << "!> Fatal error!\n";
      std::cout << "!>    Wrong ion format! should be like ";
      std::cout << " Z1A20E5.11\n";
      std::exit(1);
  }
  //
  G4double E = 0;
  G4int endl = particleName.size();
  if(Epos!=(G4int)G4String::npos){
      G4String Estr = particleName.substr(Epos+1,endl-Epos-1);
      E = std::atof(Estr.c_str());
      endl = Epos;
  }
  //
  G4String Zstr = particleName.substr(Zpos+1,Apos-Zpos-1);
  G4String Astr = particleName.substr(Apos+1,endl-Apos-1);
  G4int Z = std::atoi(Zstr.c_str());
  G4int A = std::atoi(Astr.c_str());
  //printf("A = %i; Z = %i; E = %f\n",A,Z,E);
  //
  particle = pTable->GetIon(Z,A,E);
  if(!particle){
      std::cout << "!> Fatal error!\n";
      std::cout << "!>    Cannot create the ion ";
      std::cout << particleName << "\n";
      std::exit(1);
  }
  */
  return particle;
}
