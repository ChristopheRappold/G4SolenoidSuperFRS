#ifndef _MField_hh
#define _MField_hh 1
/*
 * $Id:$
 *
 * Software developement for WASA-at-COSY
 * (c) 2006 The WASA-at-COSY Collaboration
 * Created:  2006-06-09 Andrzej Kupsc 
 * Modified:
 * License: see COPYRIGHT file
 */
/*!
  \file   MField.hh
  \brief  Magnetic field maps

  \author Andrzej Kupsc
  \date   2006-06-09
*/

#include "TList.h"

#define litecr 95
#define storcr 9500
#define storcr1 storcr-1

class MFieldSet;

class  MField: public TObject
{
public:
  /** The standard c'tor should not do any registration, because Root
   *  internally calls it to create temporary objects.
   */
  MField();

  /** Sets the name and calls Register() */ 
  explicit MField(const char* path);   

  virtual ~MField();

  TList* GetEntries() const { return fEntries; }

  void AddEntry(MFieldSet* set);
 
  /* Clear() <B>must</B> reset the fProcessed flag */

  virtual void Clear(Option_t *option = "");
  virtual void Print(Option_t *option = "") const;
  
  void InitializeParameter();
  void ReadParameter(const char* pathname);
  void SetScale(double scaling = 12.9);
  void Wfld(Double_t *xyz, Double_t *Fxyz);
  void mflcdr(double *xyz,double *fret,int *iret);
  double fint(int narg,double *arg,int *nent,double *ent,double *table);
  Int_t GetNSets(){return fNSets;}

  Double_t MinMax_R[2] = {1000.,-1000.};
  Double_t MinMax_Z[2] = {1000.,-1000.};

private:

  Double_t xyz_cr[litecr];
  Double_t fxyzcr[storcr][3];
  Int_t ltstcr;
  Int_t limrcr;
  Int_t ndimcr;
  Int_t minicr[4];
  Int_t maxicr[4];
  Int_t mirrcr[4];
  Double_t eqalcr;
  Int_t mtabcr[3][3];
  Int_t numicr[3];
  Int_t xzspcr[3];
  Int_t m0tabl[3][3];
  Bool_t testcr;
  Bool_t rlimcr;

  Double_t fMFScaling; //! Factor to scale FieldMap according to user settings
  // parameter collection containing parameter sets is stored in this data member
  Int_t fNSets; 
  TList  *fEntries;  // MFieldSet
 
};


class MFieldSet : public TObject {


protected:

  Double_t       fR;     // @HASH@ R
  Double_t       fZ;     // @HASH@ Z
  Double_t       fFr;    // Fr
  Double_t       fFz;    // Fz

public:

  MFieldSet();
  MFieldSet(Double_t r,Double_t z,Double_t fr,Double_t fz);
  virtual ~MFieldSet();
  void GetAll(Double_t *a);
  Double_t GetR() {return fR;};
  Double_t GetZ() {return fZ;};
  Double_t GetFr() {return fFr;};
  Double_t GetFz() {return fFz;};
  virtual void  Print(Option_t * option="") const;

};


#endif
