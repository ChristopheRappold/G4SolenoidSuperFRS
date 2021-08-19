#include "G4SolWASAMapMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4SystemOfUnits.hh"
#include <memory>
#include <iostream>

G4SolWASAMapMagneticField::G4SolWASAMapMagneticField(const G4String& name):nameField(name),maxField(1.*tesla)
{
  FieldMap = std::make_unique<MField>(name);
  FieldMap->ReadParameter(nameField);
}

G4SolWASAMapMagneticField::~G4SolWASAMapMagneticField()
{
  // if(FieldMap)
  //   {
  //     delete FieldMap;
  //     FieldMap= 0;
  //   }
}

void G4SolWASAMapMagneticField::SetMaxField(double maxF)
{
  maxField =  maxF;
}

void G4SolWASAMapMagneticField::SetOriginField(double initPoint[3], double signD)
{
  for(int i=0;i<3;++i)
    originField[i] = initPoint[i];

  signDir = signD;
}

void G4SolWASAMapMagneticField::InitField()
{
  std::cout<<" SetMaxField : "<<maxField/kilogauss<<" kG"<<" "<<maxField/tesla<<" T"<<" "<<maxField<<"\n";
  FieldMap->SetScale(maxField/kilogauss);
  FieldMap->InitializeParameter();

  double checkB[3] = {0.,0.,0.};
  GetFieldValue(originField, checkB);
  std::cout<<" --> Check G4 field at : ["<<originField[0]/cm<<", "<<originField[1]/cm<<", "<<originField[2]/cm<<"] cm : Bxyz ["<<checkB[0]/tesla<<", "<<checkB[1]/tesla<<", "<<checkB[2]/tesla<<"] T \n";
}

void G4SolWASAMapMagneticField::GetFieldValue(const double Point[3],double *Bfield) const
{

  Double_t B[3];
  Double_t pos[3];
  
  //std::cout<<" MagField "<<Point[0]<<" "<<Point[1]<<" "<<Point[2]<<" ";

  pos[0]=signDir*(Point[0]-originField[0])/cm;
  pos[1]=signDir*(Point[1]-originField[1])/cm;
  pos[2]=signDir*(Point[2]-originField[2])/cm;

  FieldMap->Wfld(pos,B);
  
  Bfield[0] = B[0]*kilogauss;
  Bfield[1] = B[1]*kilogauss;
  Bfield[2] = B[2]*kilogauss;
  
  //std::cout<<"["<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"] :"<<" ["<<B[0]<<" "<<B[1]<<" "<<B[2]<<"] / "<<Bfield[2]<<" T"<<G4endl;

}
