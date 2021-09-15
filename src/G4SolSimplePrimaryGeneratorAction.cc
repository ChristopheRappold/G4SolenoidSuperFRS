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
// Implementation of the G4SolSimplePrimaryGeneratorAction class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "G4SolSimplePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolSimplePrimaryGeneratorAction::G4SolSimplePrimaryGeneratorAction(const G4SolConfig& conf)
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0), fMessenger(0), nameParticle("pi-"), fMomentum(1000. * MeV),
      fSigmaMomentum(50. * MeV), fSigmaAngle(2. * deg), Par(conf)

{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);

  nameParticle = Par.Get<std::string>("Particle");

  if(Par.IsAvailable("Beam_KineticEnergy"))
    {
      fKineticE = Par.Get<double>("Beam_KineticEnergy");
      G4cout << "!> Set Primary Generator with beam kinetic energy !\n";
    }
  else
    fKineticE = -1;

  if(Par.IsAvailable("Beam_PosRandom"))
    fRandomizePrimary[0] = Par.Get<int>("Beam_PosRandom"); //== 1 ? true : false;

  if(Par.IsAvailable("Beam_TotalMomRandom"))
    fRandomizePrimary[1] = Par.Get<int>("Beam_TotalMomRandom"); //== 1 ? true : false;

  if(Par.IsAvailable("Beam_DirMomRandom"))
    fRandomizePrimary[2] = Par.Get<int>("Beam_DirMomRandom"); //== 1 ? true : false;

  if(Par.IsAvailable("Beam_DirThetaRandom"))
    fRandomizePrimary[3] = Par.Get<int>("Beam_DirThetaRandom");

  if(Par.IsAvailable("Beam_DirPhiRandom"))
    fRandomizePrimary[4] = Par.Get<int>("Beam_DirPhiRandom");

  fMomentum      = Par.Get<double>("Beam_Momentum");
  fSigmaMomentum = Par.Get<double>("Beam_MomentumSigma");
  // fSigmaAngle;
  fDirX = Par.Get<double>("Beam_MomentumX");
  fDirY = Par.Get<double>("Beam_MomentumY");
  fDirZ = Par.Get<double>("Beam_MomentumZ");

  fDir.set(fDirX, fDirY, fDirZ);

  fDirXSigma = Par.Get<double>("Beam_MomentumXsigma");
  fDirYSigma = Par.Get<double>("Beam_MomentumYsigma");
  fDirZSigma = Par.Get<double>("Beam_MomentumZsigma");

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

  ConstParticle = GetParticle("geantino"); // nullptr;//GetParticle(nameParticle);

  fParticleGun->SetParticlePosition(G4ThreeVector(fPosX, fPosY, fPosZ));
  fParticleGun->SetParticleDefinition(ConstParticle);
  fParticleGun->SetParticleEnergy(fMomentum);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fDirX, fDirY, fDirZ));

  if(Par.IsAvailable("Beam_PosX"))
    fPosX = Par.Get<double>("Beam_PosX");
  if(Par.IsAvailable("Beam_PosY"))
    fPosY = Par.Get<double>("Beam_PosY");
  if(Par.IsAvailable("Beam_PosZ"))
    fPosZ = Par.Get<double>("Beam_PosZ");

  if(Par.IsAvailable("Beam_DirTheta_1"))
    fDirTheta_1 = Par.Get<double>("Beam_DirTheta_1");
  if(Par.IsAvailable("Beam_DirTheta_2"))
    fDirTheta_2 = Par.Get<double>("Beam_DirTheta_2");

  if(Par.IsAvailable("Beam_DirPhi_1"))
    fDirPhi_1 = Par.Get<double>("Beam_DirPhi_1");
  if(Par.IsAvailable("Beam_DirPhi_2"))
    fDirPhi_2 = Par.Get<double>("Beam_DirPhi_2");

  currentEvent = 0;
  eventPerBin  = -1;
  setBinRand   = false;
  currentBin   = std::make_tuple(0, 0);

  if(Par.IsAvailable("EventPerRandomBin"))
    {
      eventPerBin = Par.Get<double>("EventPerRandomBin");
      setBinRand  = true;
    }

  if(Par.IsAvailable("NbRandomBinPerGen"))
    {
      NbBin = Par.Get<int>("NbRandomBinPerGen");
      if(Par.IsAvailable("Beam_DirTheta_1") && Par.IsAvailable("Beam_DirTheta_2"))
        {
          double temp_bin_l = fDirTheta_1, temp_bin_h = 0;
          double binStep = (fDirTheta_2 - fDirTheta_1) / static_cast<double>(NbBin);
          for(int iBin = 0; iBin < NbBin; ++iBin)
            {
              temp_bin_h = temp_bin_l + binStep;
              ThetaBins.emplace_back(std::make_tuple(temp_bin_l, temp_bin_h));
              temp_bin_l = temp_bin_h;
            }
        }
      else
        ThetaBins.emplace_back(std::make_tuple(fDirTheta_1, fDirTheta_2));

      if(Par.IsAvailable("Beam_DirPhi_1") && Par.IsAvailable("Beam_DirPhi_2"))
        {
          double temp_bin_l = fDirPhi_1, temp_bin_h = 0;
          double binStep = (fDirPhi_2 - fDirPhi_1) / static_cast<double>(NbBin);
          for(int iBin = 0; iBin < NbBin; ++iBin)
            {
              temp_bin_h = temp_bin_l + binStep;
              PhiBins.emplace_back(std::make_tuple(temp_bin_l, temp_bin_h));
              temp_bin_l = temp_bin_h;
            }
        }
      else
        PhiBins.emplace_back(std::make_tuple(fDirPhi_1, fDirPhi_2));
    }
  else
    {
      ThetaBins.emplace_back(std::make_tuple(fDirTheta_1, fDirTheta_2));
      PhiBins.emplace_back(std::make_tuple(fDirPhi_1, fDirPhi_2));
    }

  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolSimplePrimaryGeneratorAction::~G4SolSimplePrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolSimplePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  ConstParticle = GetParticle(nameParticle); // GetParticle("Lambda");

  // ConstParticle = GetParticle(Par.Get_Beam_Type());
  // G4cout<<" PrimaryGenerator :"<<Par.Get<std::string>("Particle")<<" "<<ConstParticle<<" "<<nameParticle<<G4endl;

  G4double mom;

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
	      ofpos_z      = fTargetSizeZ * (0.5 - G4UniformRand());
	    }
	  else
	    {
	      ofpos_x = fSpotSizeSigmaX * RandTable[1](0,1);
	      ofpos_y = fSpotSizeSigmaY * RandTable[1](0,1);
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

  fParticleGun->SetParticleDefinition(ConstParticle);

  G4ThreeVector TempDir = fDir;

  if(fRandomizePrimary[2])
    {
      double xdir = fDirX; // Par.Get_Beam_MomentumDirectionX();
      double ydir = fDirY; // Par.Get_Beam_MomentumDirectionY();
      double zdir = fDirZ; // Par.Get_Beam_MomentumDirectionZ();

      xdir += fDirXSigma * std::sqrt(-std::log(G4UniformRand()));
      ydir += fDirYSigma * std::sqrt(-std::log(G4UniformRand()));
      zdir = sqrt(1.0 - ydir * ydir - xdir * xdir);
      // std::cout<<xdir<<" "<<ydir<<" "<<zdir<<std::endl;
      TempDir.set(xdir, ydir, zdir);
      // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xdir,ydir,zdir));

      // G4double angle = (G4UniformRand()-0.5)*fSigmaAngle;
      // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(std::sin(angle),0.,std::cos(angle)));
    }

  if(fRandomizePrimary[3])
    {
      double temp_theta = RandTable[fRandomizePrimary[3]](std::get<0>(ThetaBins[std::get<0>(currentBin)]),
                                                          std::get<1>(ThetaBins[std::get<0>(currentBin)]));
      TempDir.setTheta(temp_theta);
    }

  if(fRandomizePrimary[4])
    {
      double temp_phi = RandTable[fRandomizePrimary[4]](std::get<0>(PhiBins[std::get<1>(currentBin)]),
                                                        std::get<1>(PhiBins[std::get<1>(currentBin)]));
      TempDir.setPhi(temp_phi);
    }

  fParticleGun->SetParticleMomentumDirection(TempDir);

  G4double Ekin = 0;
  if(fKineticE < 0)
    {
      mom = fMomentum; // Par.Get_Beam_KineticEnergy();
      if(fRandomizePrimary[1])
        mom += std::sqrt(-std::log(G4UniformRand())) * fSigmaMomentum;

      G4double pp   = mom;
      G4double mass = ConstParticle->GetPDGMass();
      Ekin          = std::sqrt(pp * pp + mass * mass) - mass;
    }
  else
    {
      Ekin = fKineticE;
      if(fRandomizePrimary[1])
        Ekin += std::sqrt(-std::log(G4UniformRand())) * fSigmaMomentum;
    }

  fParticleGun->SetParticleEnergy(Ekin);

  fParticleGun->SetParticlePosition(G4ThreeVector(pos_x, pos_y, pos_z));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  ++currentEvent;

  if(currentEvent > eventPerBin && setBinRand == true)
    {
      ++std::get<0>(currentBin);
      if(std::get<0>(currentBin) > ThetaBins.size())
        {
          std::get<0>(currentBin) = 0;
          ++std::get<1>(currentBin);
          if(std::get<1>(currentBin) > PhiBins.size())
            std::get<1>(currentBin) = 0;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolSimplePrimaryGeneratorAction::DefineCommands()
{
  // Define /G4SolSimple/generator command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/G4SolSimple/generator/", "Primary generator control");

  // momentum command
  G4GenericMessenger::Command& momentumCmd =
      fMessenger->DeclarePropertyWithUnit("momentum", "GeV", fMomentum, "Mean momentum of primaries.");
  momentumCmd.SetParameterName("p", true);
  momentumCmd.SetRange("p>=0.");
  momentumCmd.SetDefaultValue("1.");

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

G4ParticleDefinition* G4SolSimplePrimaryGeneratorAction::GetParticle(const G4String& particleName)
{

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();

  G4IonTable* pTableIon = G4IonTable::GetIonTable();
  pTableIon->CreateAllIon();

  G4ParticleDefinition* particle = nullptr;

  G4int Apos = particleName.rfind("A");
  G4int Zpos = particleName.rfind("Z");
  G4int Epos = particleName.rfind("E");
  G4int Lpos = particleName.rfind("L");

  if(Zpos != (G4int)G4String::npos)
    {
      // std::cout<<"zpos >0 "<<particleName<<std::endl;
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
          //particle->DumpTable();
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
          //particle->DumpTable();
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
          //particle->DumpTable();
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
          //particle->DumpTable();
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
          //particle->DumpTable();
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

  return particle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
