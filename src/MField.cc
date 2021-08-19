
/*
 * $Id:$
 *
 * Software developement for WASA-at-COSY
 * (c) 2006 The WASA-at-COSY Collaboration
 * Created:  2006-06-09 Andrzej Kupsc
 * Modified:
 * License: see COPYRIGHT file
 */

#include <fstream>
#include <iostream>
#include <string>

//#include <stdio.h>
//#include "stdlib.h"

#include "MField.hh"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

MField::MField()
{
  fEntries   = 0;
  fNSets     = 0;
  fMFScaling = 0;
}

MField::MField(const char*)
{

  fEntries = 0;
  fNSets   = 0;

  // Factor to scale FieldMap according
  // to average value given by user
  fMFScaling = 1. / 12.9;

  ltstcr = 0;
  limrcr = 1;
  //
  ndimcr    = 2;
  mirrcr[0] = 0;
  mirrcr[1] = 0;
  mirrcr[2] = 0;
  mirrcr[3] = 0;
  minicr[0] = 0;
  minicr[1] = 0;
  minicr[2] = -20;
  minicr[3] = 0;
  maxicr[0] = 0;
  maxicr[1] = 19;
  maxicr[2] = 20;
  maxicr[3] = 0;
  //
  eqalcr = 1.e-4;

  //
  //-----preliminary code for scaling of central detector field map;
  //     flexible field map (rr 020821);
  //
  //     rectangular field map (hc 020704, update rr 020821);
  //
  //-----set the min/max indices for the rectangular field map;
  //     this overwrites the values defined in the block data;
  //
  m0tabl[0][0] = -1;
  m0tabl[1][0] = 1;
  m0tabl[2][0] = 1;
  m0tabl[0][1] = 1;
  m0tabl[1][1] = -1;
  m0tabl[2][1] = 1;
  m0tabl[0][2] = -1;
  m0tabl[1][2] = -1;
  m0tabl[2][2] = 1;
  //
  if(ltstcr > 0)
    testcr = true;
  else
    testcr = false;
  //
  if(limrcr > 0)
    rlimcr = true;
  else
    rlimcr = false;
  //
}

void MField::SetScale(double scaling) { fMFScaling = scaling / 12.9; }

MField::~MField()
{
  if(fEntries)
    delete fEntries;
}

void MField::Clear(Option_t*)
{
  if(fEntries != 0)
    fEntries->Delete();
}

void MField::Print(Option_t*) const
{

  if(fEntries)
    {
      TList* rList = fEntries;
      TIter it(rList);
      MFieldSet* entry;
      while((entry = dynamic_cast<MFieldSet*>(it())))
        {
          entry->Print();
        }
    }
}

void MField::InitializeParameter()
{

  //-----initialise pointers:;
  Int_t coord       = 1;
  minicr[0]         = 0;
  maxicr[0]         = 0;
  numicr[coord - 1] = maxicr[coord] - minicr[coord] + 1;
  xzspcr[coord - 1] = 3 + TMath::Max(abs(minicr[coord]), abs(maxicr[coord])) +
                      TMath::Max(abs(minicr[coord - 1]), abs(maxicr[coord - 1]));
  for(coord = 2; coord <= 3; coord++)
    {
      numicr[coord - 1] = maxicr[coord] - minicr[coord] + 1;
      xzspcr[coord - 1] = xzspcr[coord - 2] + 3 + TMath::Max(abs(minicr[coord]), abs(maxicr[coord])) +
                          TMath::Max(abs(minicr[coord - 1]), abs(maxicr[coord - 1]));
    }
  // Int_t nx = 1;
  Int_t ny = numicr[0] + 2;
  Int_t nz = (numicr[0] + 2) * (numicr[1] + 2);
  Int_t mx = minicr[1] - 1;
  Int_t my = minicr[2] - 1;
  Int_t mz = minicr[3] - 1;

  // Initialize lookup table
  fNSets       = 0;
  TList* rList = fEntries;
  std::cout << "MField: " << rList->GetEntries() << " sets loaded" << std::endl;
  TIter it(rList);
  MFieldSet* entry;
  while((entry = dynamic_cast<MFieldSet*>(it())))
    {

      Int_t ix = TMath::Nint(entry->GetR() / 5);
      Int_t iy = TMath::Nint(entry->GetZ() / 5);
      Int_t iz = 0;

      xyz_cr[xzspcr[0] + ix - 1]                             = entry->GetR();
      xyz_cr[xzspcr[1] + iy - 1]                             = entry->GetZ();
      fxyzcr[(ix - mx) + ny * (iy - my) + nz * (iz - mz)][0] = entry->GetFr() * fMFScaling;
      fxyzcr[(ix - mx) + ny * (iy - my) + nz * (iz - mz)][1] = entry->GetFz() * fMFScaling;

      fNSets++;
    }
  std::cout << "MField: initialized #sets:" << fNSets << std::endl;

  //
  //-----mirror fieldmap-coordinate-arrays in case of symmetry
  //
  for(Int_t coord2 = 1; coord2 <= ndimcr; coord2++)
    {
      if(mirrcr[coord2] == 1)
        {
          mirrcr[0]      = 1;
          minicr[coord2] = -maxicr[coord2];
          for(Int_t i = minicr[coord2]; i <= -1; i++)
            {
              xyz_cr[xzspcr[coord2 - 1] + i - 1] = -xyz_cr[xzspcr[coord2 - 1] - i - 1];
            }
          mtabcr[0][coord2 - 1] = m0tabl[0][coord2 - 1];
          mtabcr[1][coord2 - 1] = m0tabl[1][coord2 - 1];
          mtabcr[2][coord2 - 1] = m0tabl[2][coord2 - 1];
        }
      else
        {
          mtabcr[0][coord2 - 1] = 1;
          mtabcr[1][coord2 - 1] = 1;
          mtabcr[2][coord2 - 1] = 1;
        }
      //-----add the points which are needed for interpolation outside fieldmap
      xyz_cr[xzspcr[coord2 - 1] + minicr[coord2] - 2] = -150;
      xyz_cr[xzspcr[coord2 - 1] + maxicr[coord2] + 0] = 150;
    }

  //
  //----- test: find the field strength in the origin of coordinate system
  //
  double xyz[3], f[3];
  Int_t iret = 0;

  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  f[0]   = 0.0;
  f[1]   = 0.0;
  f[2]   = 0.0;

  mflcdr(xyz, f, &iret);

  if(iret != 0)
    {
      std::cout << "*************************************" << std::endl;
      std::cout << " error in inmfl at initial mflcdr call" << std::endl;
      std::cout << "*************************************" << std::endl;
    }
  else
    {
      std::cout << "*****************************************" << std::endl;
      std::cout << "\t  Checking FieldMap at point (0,0,0):" << std::endl;
      std::cout << "\t  Components of field are: (" << f[0] << "," << f[1] << "," << f[2] << ") kGs" << std::endl;
      std::cout << "\t  Magnitude  of field is:  " << sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]) << " kGs"
                << std::endl;
      std::cout << "*****************************************" << std::endl;
    }
}

void MField::ReadParameter(const char* filename)
{

  // basename here is name of input file

  std::ifstream in(filename);
  if(!in.good())
    {
      std::cout << "MField::ReadParameter - cannot open file: " << filename << std::endl;
      return;
    }

  MFieldSet* set  = 0;
  Int_t CountSets = 0;
  TString line;
  while(true)
    {
      line.ReadLine(in);
      if(in.eof())
        break;
      if(line.BeginsWith("#"))
        continue;
      TString delim = "\n\t ";
      TObjArray* aa = line.Tokenize(delim);
      if(aa->GetEntries() != 4)
        continue;

      TString st;
      st          = dynamic_cast<TObjString*>(aa->At(0))->String();
      Double_t r  = st.Atof();
      st          = dynamic_cast<TObjString*>(aa->At(1))->String();
      Double_t z  = st.Atof();
      st          = dynamic_cast<TObjString*>(aa->At(2))->String();
      Double_t Fr = st.Atof();
      st          = dynamic_cast<TObjString*>(aa->At(3))->String();
      Double_t Fz = st.Atof();

      set = new MFieldSet(r, z, Fr, Fz);
      AddEntry(set);
      CountSets++;
      delete aa;
    }
  in.close();
  std::cout << "*****************************************" << std::endl;
  std::cout << "MField: " << CountSets << " sets filled from " << filename << std::endl;
}

void MField::mflcdr(double* xyz, double* fret, int* iret)
{
  //-----units are cm and kgaus;
  //
  int i, ixyz[3], na[3], coord, fcoord, nx, ny, nz, mx, my, mz, npart, npart0, n, ntwo;
  int ix, iy, iz, isx, isy, isz, ixx, iyy, izz;
  double mirfac;
  double f[3], xxyyzz[9], fxyz3[3][3][3], xyzrz[3], fxyz2[3][3], xratio = 0., yratio = 0.;
  //

  na[0] = 3;
  na[1] = 3;
  na[2] = 3;
  //
  *iret     = 0;
  numicr[0] = maxicr[1] - minicr[1] + 1;
  numicr[1] = maxicr[2] - minicr[2] + 1;
  numicr[2] = maxicr[3] - minicr[3] + 1;
  nx        = 1;
  ny        = numicr[0] + 2;
  nz        = (numicr[0] + 2) * (numicr[1] + 2);
  mx        = minicr[1] - 1;
  my        = minicr[2] - 1;
  mz        = minicr[3] - 1;
  //-----transformation to cylindercoordinates if ndimcr=2
  if(ndimcr == 3)
    {
      xyzrz[0] = xyz[0];
      xyzrz[1] = xyz[1];
      xyzrz[2] = xyz[2];
    }
  else
    {
      xyzrz[0] = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
      xyzrz[1] = xyz[2];
      if(xyzrz[0] <= eqalcr)
        {
          xratio = 0.;
          yratio = 0.;
        }
      else
        {
          xratio = xyz[0] / xyzrz[0];
          yratio = xyz[1] / xyzrz[0];
        }
    }
  //-----find indices ixyz(1-ndimcr) of closest existing coordinates:;
  for(coord = 1; coord <= ndimcr; coord++)
    {
      npart0 = maxicr[coord] - minicr[coord];
      n      = minicr[coord];
      ntwo   = 2;
      do
        {
          npart = npart0 / ntwo;
          ntwo  = ntwo * 2;
          if(xyzrz[coord - 1] >= xyz_cr[xzspcr[coord - 1] + n + npart - 1])
            n += npart;
        }
      while(npart > 1);
      //
      while(((xyzrz[coord - 1] - xyz_cr[xzspcr[coord - 1] + n - 1]) >=
             (xyz_cr[xzspcr[coord - 1] + n] - xyzrz[coord - 1])) &&
            ((n + 1) <= maxicr[coord]))
        {
          n++;
          if((n + 1) > maxicr[coord])
            break;
        }
      if(n >= maxicr[coord])
        {
          if(xyzrz[coord - 1] > xyz_cr[xzspcr[coord - 1] + maxicr[coord] - 1])
            *iret = -3;
          if(xyzrz[coord - 1] > 150.)
            {
              *iret = -4;
              //	    print *, ' error 2 in mflcdr';
              return;
            }
          else
            {
              ixyz[coord - 1] = maxicr[coord];
            }
        }
      else if(n <= minicr[coord])
        {
          if(xyzrz[coord - 1] < xyz_cr[xzspcr[coord - 1] + minicr[coord] - 1])
            *iret = -3;
          if(xyzrz[coord - 1] < -150.)
            {
              //	    print *, ' error 4 in mflcdr';
              //	        std::cout << "MField: error #4 in mflcdr " <<std::endl;
              *iret = -4;
              return;
            }
          else
            {
              ixyz[coord - 1] = minicr[coord] + 1;
            }
        }
      else
        {
          ixyz[coord - 1] = n;
        }
    }
  //-----fill interpolation arrays:;
  for(coord = 1; coord <= ndimcr; coord++)
    {
      for(i = -1; i <= +1; i++)
        {
          xxyyzz[3 * coord + i - 2] = xyz_cr[xzspcr[coord - 1] + ixyz[coord - 1] + i - 1];
        }
    }
  //-----cylinder-symmetrical field map:;
  if(mirrcr[0] == 0)
    {
      if(ndimcr != 3)
        {
          for(fcoord = 1; fcoord <= 2; fcoord++)
            {
              for(iy = ixyz[1] - 1; iy <= (ixyz[1] + 1); iy++)
                {
                  for(ix = ixyz[0] - 1; ix <= (ixyz[0] + 1); ix++)
                    {
                      fxyz2[1 + iy - ixyz[1]][1 + ix - ixyz[0]] =
                          fxyzcr[(ix - mx) + ny * (iy - my) + nz * (0 - mz)][fcoord - 1];
                    }
                }
              //-------2-d linear interpolation using kern-lib s/r fint (e104);
              f[fcoord - 1] = fint(2, &xyzrz[0], na, &xxyyzz[0], &fxyz2[0][0]);
            }
          fret[0] = f[0] * xratio;
          fret[1] = f[0] * yratio;
          fret[2] = f[1];
          return;
        }
      else
        {
          //-----asymmetrical field map:;
          for(fcoord = 1; fcoord <= 3; fcoord++)
            {
              for(iz = ixyz[2] - 1; iz <= ixyz[2] + 1; iz++)
                {
                  for(iy = ixyz[1] - 1; iy <= ixyz[1] + 1; iy++)
                    {
                      for(ix = ixyz[0] - 1; ix <= ixyz[0] + 1; ix++)
                        {
                          fxyz3[1 + ix - ixyz[0]][1 + iy - ixyz[1]][1 + iz - ixyz[2]] =
                              fxyzcr[(ix - mx) + ny * (iy - my) + nz * (iz - mz)][fcoord - 1];
                        }
                    }
                }
              //-------3-d linear interpolation using kern-lib s/r fint (e104);
              f[fcoord - 1] = fint(3, &xyz[0], na, &xxyyzz[0], &fxyz3[0][0][0]);
            }
          fret[0] = f[0];
          fret[1] = f[1];
          fret[2] = f[2];
          return;
        }
    }
  else
    {
      //-----symmetrical field map:;
      for(fcoord = 1; fcoord <= 3; fcoord++)
        {
          for(iz = ixyz[2] - 1; iz <= ixyz[2] + 1; iz++)
            {
              for(iy = ixyz[1] - 1; iy <= ixyz[1] + 1; iy++)
                {
                  for(ix = ixyz[0] - 1; ix <= ixyz[0] + 1; ix++)
                    {
                      isx    = ix >= 0 ? +1 : -1;
                      isy    = iy >= 0 ? +1 : -1;
                      isz    = iz >= 0 ? +1 : -1;
                      ixx    = ix * (mirrcr[1] == 0 ? +1 : isx);
                      iyy    = iy * (mirrcr[2] == 0 ? +1 : isy);
                      izz    = iz * (mirrcr[3] == 0 ? +1 : isz);
                      mirfac = (isx == 1 ? 1 : mtabcr[fcoord - 1][0]) * (isy == 1 ? 1 : mtabcr[fcoord - 1][1]) *
                               (isz == 1 ? 1 : mtabcr[fcoord - 1][2]);
                      fxyz3[1 + ix - ixyz[0]][1 + iy - ixyz[1]][1 + iz - ixyz[2]] =
                          mirfac * fxyzcr[(ixx + 1) + ny * (iyy + 1) + nz * (izz + 1)][fcoord - 1];
                    }
                }
            }
          //-------3-d linear interpolation using kern-lib s/r fint (e104);
          f[fcoord - 1] = fint(3, &xyz[0], na, &xxyyzz[0], &fxyz3[0][0][0]);
        }
    }
  if(ndimcr == 3)
    {
      fret[0] = f[0];
      fret[1] = f[1];
      fret[2] = f[2];
    }
  else
    {
      fret[0] = f[0] * xratio;
      fret[1] = f[0] * yratio;
      fret[2] = f[1];
    }
  return;
}

void MField::AddEntry(MFieldSet* set)
{
  if(!fEntries)
    fEntries = new TList();
  fEntries->Add(set);
}

void MField::Wfld(double* xyzd, double* fd)
{
  double xyz[3], f[3];
  int iretcd = 0;
  //-----give magn. field at point xyz (options: "iflg" mfl, "w4mf");
  xyz[0] = xyzd[0];
  xyz[1] = xyzd[1];
  xyz[2] = xyzd[2];
  if(xyz[2] > 250.)
    {
      f[0] = 0.;
      f[1] = 0.;
      f[2] = 0.;
    }
  else
    {
      //------cd field;
      f[0] = 0.;
      f[1] = 0.;
      f[2] = 0.;
      mflcdr(xyz, f, &iretcd);
      //	if (iretcd<0) gLog(CLog::kError) << "MField::Wfld(): mflcdr returns " << iretcd
      //	     << " at " << Form("x=(%f,%f,%f)\n",xyz[0],xyz[1],xyz[2]);
    }
  fd[0] = f[0];
  fd[1] = f[1];
  fd[2] = f[2];
}

double MField::fint(int narg, double* arg, int* nent, double* ent, double* table)
{
  //
  //   interpolation routine. author c. letertre.;
  //   modified by b. schorr, 1.07.1982.;
  //
  int inde[32], k, istep, knots, n, ndim, loca, lmin, lmax, ishift, locb, locc, i;
  double weight[32], x, h, fin, eta;
  //
  fin = 0.;
  if(narg < 1 || narg > 5)
    return fin;
  lmax      = 0;
  istep     = 1;
  knots     = 1;
  inde[0]   = 1;
  weight[0] = 1.;
  for(n = 1; n <= narg; n++)
    {
      x    = arg[n - 1];
      ndim = nent[n - 1];
      loca = lmax;
      lmin = lmax + 1;
      lmax = lmax + ndim;
      if(ndim <= 2)
        {
          if(ndim == 1)
            continue;
          h = x - ent[lmin - 1];
          if(h == 0.)
            {
              istep = istep * ndim;
              continue;
            }
          ishift = istep;
          if(x - ent[lmin] == 0.)
            {
              for(k = 1; k <= knots; k++)
                {
                  inde[k - 1] = inde[k - 1] + ishift;
                };
              istep = istep * ndim;
              continue;
            }
          ishift = 0;
          eta    = h / (ent[lmin] - ent[lmin - 1]);
          for(k = 1; k <= knots; k++)
            {
              inde[k - 1]           = inde[k - 1] + ishift;
              inde[k + knots - 1]   = inde[k - 1] + istep;
              weight[k + knots - 1] = weight[k - 1] * eta;
              weight[k - 1]         = weight[k - 1] - weight[k + knots - 1];
            }
          knots = 2 * knots;
          istep = istep * ndim;
          continue;
        }
      locb = lmax + 1;
      do
        {
          locc = (loca + locb) / 2;
          if(x - ent[locc - 1] < 0.)
            locb = locc;
          else if(x - ent[locc - 1] > 0.)
            loca = locc;
          else
            {
              ishift = (locc - lmin) * istep;
              for(k = 1; k <= knots; k++)
                {
                  inde[k - 1] = inde[k - 1] + ishift;
                };
              istep = istep * ndim;
              goto cycle;
            }
        }
      while((locb - loca) > 1);
      loca   = loca > lmin ? loca : lmin;
      loca   = loca < lmax - 1 ? loca : lmax - 1;
      ishift = (loca - lmin) * istep;
      eta    = (x - ent[loca - 1]) / (ent[loca] - ent[loca - 1]);
      for(k = 1; k <= knots; k++)
        {
          inde[k - 1]           = inde[k - 1] + ishift;
          inde[k + knots - 1]   = inde[k - 1] + istep;
          weight[k + knots - 1] = weight[k - 1] * eta;
          weight[k - 1]         = weight[k - 1] - weight[k + knots - 1];
        }
      knots = 2 * knots;
      istep = istep * ndim;
    cycle:
      continue;
    }
  for(k = 1; k <= knots; k++)
    {
      i   = inde[k - 1];
      fin = fin + weight[k - 1] * table[i - 1];
    }
  return fin;
}

MFieldSet::MFieldSet()
{
  fR  = 0.;
  fZ  = 0.;
  fFr = 0.;
  fFz = 0.;
}

MFieldSet::MFieldSet(Double_t r, Double_t z, Double_t fr, Double_t fz)
{
  fR  = r;
  fZ  = z;
  fFr = fr;
  fFz = fz;
}

MFieldSet::~MFieldSet() {}

void MFieldSet::GetAll(Double_t* a)
{
  double* t = a;
  *t++      = fR;
  *t++      = fZ;
  *t++      = fFr;
  *t++      = fFz;
}

void MFieldSet::Print(Option_t*) const { std::cout << fR << "\t" << fZ << "\t" << fFr << "\t" << fFz << std::endl; }
