#ifndef WASAMap_MagneticField_H
#define WASAMap_MagneticField_H

#include "globals.hh"
#include "G4MagneticField.hh"
#include "MField.hh"
#include <memory>

class G4FieldManager;

class G4SolWASAMapMagneticField : public G4MagneticField
{
    public:
  explicit G4SolWASAMapMagneticField(const G4String& namefile_field);
  ~G4SolWASAMapMagneticField();

  void GetFieldValue( const  double Point[3], double *Bfield ) const;
  void SetMaxField(double maxField );
  void SetOriginField(double initPoint[3], double signDir);
  void InitField();
private:

  std::unique_ptr<MField> FieldMap;

  double originField[3] = {0., 0., 0.};
  double signDir;
  G4String nameField;
  double maxField;

};

#endif

