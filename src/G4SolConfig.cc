// ------------------------------------------------------- 
// Implementation of the G4SolConfig class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "G4SolConfig.hh"
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>

#include "G4SystemOfUnits.hh"

#include <getopt.h>

static struct option optlong[] = 
  {
    {"help",0,NULL,'h'},
    {"gui",0,NULL,'g'},
    {"mac",1,NULL,'m'},
    {"input",1,NULL,'i'},
    {"config",1,NULL,'c'},
  };
G4SolConfig::G4SolConfig()
{
  SetDefault();
}

G4SolConfig::G4SolConfig(int argc,char** argv)
{
  SetDefault();
  status = ParseCmd(argc,argv);
}

void G4SolConfig::ParseConfig(const std::string& namefile)
{
  std::cout<<"start reading"<<std::endl;
  std::ifstream ifs ( namefile.c_str() );
  if(ifs.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs,temp_line))
	{
	  std::stringstream stream(temp_line);
	  std::string testComment(stream.str());
	  auto it_comment = testComment.find(CommentSymbol);
	  if(it_comment!=std::string::npos)
	    {
	      //std::cout<<"!> Skip "<<test<<std::endl;
	      continue;
	    }
	  std::string key, value, unit;
	  stream >> key >> value >> unit;
	  std::cout<<"key :"<<key<<" "<<value<<" "<<unit<<std::endl;
	  if(unit.empty()==false)
	    {
	      std::string key1(key+".unit");
	      tree.put(key,value);
	      tree.put(key1,unit);
	      std::string key2(key1+"."+unit);
	      double valUnit=GetDimension(unit);
	      tree.put(key2,valUnit);
	    }
	  else
	    tree.put(key,value);
	}
    }
}
G4SolConfig::~G4SolConfig()
{
  
}

int G4SolConfig::ParseCmd(int argc,char** argv)
{

  auto print_help = [&argv] () {
    std::cout << "!> Wrong number of parameters!\n";
    std::cout << "!> Example of use:\n";
    std::cout << "!> " << argv[0];
    std::cout << "!> [-h] [--help] [-g] [--gui] [-i inputfile] [--input inputfile] [-m run.mac] [--mac run.mac] [-c config.par] [--config config.par] Outputfile.root \n";
    std::cout << " \n";
  };
  
  if (argc < 2)
    {
      print_help();
      return -1;
    }
  int option_char;
  std::string nameI, nameM, nameC("testconfig.par");
  while ((option_char = getopt_long (argc,argv,"+hgm:i:c:",optlong,NULL)) != EOF)
    switch (option_char)
      {  
      case 'h': print_help(); return -1; break;
      case 'g': std::cout<<"Gui mode "<<std::endl; tree.put("Gui",1); break; 
      case 'i': std::cout<<"Inputfile of Event :"<<optarg<<std::endl; nameI=optarg; tree.put("InputFile",nameI); break;
      case 'm': std::cout<<"Macro File :"<<optarg<<std::endl; nameM=optarg; tree.put("MacroFile",nameM); break;
      case 'c': std::cout<<"Configuration File :"<<optarg<<std::endl; nameC=optarg; tree.put("Config",nameC); break;
      case '?': print_help(); return -1;
      }

  std::string name_out;

  if(optind == argc)
    {
      print_help();
      return -1;
    }
  else
    {
      name_out = argv[optind];
      tree.put("Output_Namefile",name_out);
    }

  ParseConfig(nameC);
  
  return 0;
  
  
}

void G4SolConfig::display(const int depth, const pt::ptree& t)
{  
  BOOST_FOREACH(pt::ptree::value_type const& v, t.get_child("") )
    {  
      pt::ptree subtree = v.second;  
      std::string nodestr = t.get<std::string>(v.first);  
      
      // print current node  
      std::cout << std::string("").assign(depth*2,' ') << "* ";  
      std::cout << v.first;  
      if ( nodestr.length() > 0 )   
	std::cout << " -> \"" << t.get<std::string>(v.first) << "\"";  
      std::cout << std::endl;  
      
      // recursive go down the hierarchy  
      display(depth+1,subtree);  
    }  
}

void G4SolConfig::CheckConfig()
{
  display(0,tree);
}

double G4SolConfig::GetDimension(const std::string& dimension)
{
    // Time
    if(dimension == "ns") return ns;
    if(dimension == "ms") return ms;
    if(dimension == "s") return s;
    // Length
    if(dimension == "m") return m;
    if(dimension == "cm") return cm;
    if(dimension == "mm") return mm;
    // Fields
    if(dimension == "tesla") return tesla;
    if(dimension == "kG") return kilogauss;
    // Energy
    if(dimension == "eV") return eV;
    if(dimension == "keV") return keV;
    if(dimension == "MeV") return MeV;
    if(dimension == "GeV") return GeV;
    // Angle
    if(dimension == "degree") return deg;
    if(dimension == "rad") return rad;
    //
    std::cerr << "!> Unckown dimension " << dimension << "\n";
    std::cerr << "!> Exiting !!!\n";
    std::exit(1);

}

void G4SolConfig::SetDefault()
{
  tree.put("Gui",0);
  
  tree.put("Particle","proton");

  tree.put("Beam_SpotSizeSigma",1);
  tree.put("Beam_SpotSizeSigma.unit","cm");
  tree.put("Beam_SpotSizeSigma.unit.cm",1.*cm);

  tree.put("Beam_Momentum",2);
  tree.put("Beam_Momentum.unit","GeV");
  tree.put("Beam_Momentum.unit.GeV",1.*GeV);

  tree.put("Beam_MomentumSigma",10.);
  tree.put("Beam_MomentumSigma.unit","MeV");
  tree.put("Beam_MomentumSigma.unit.MeV",1.*MeV);
    
  tree.put("Beam_MomentumX",0.);   
  tree.put("Beam_MomentumY",0.);   
  tree.put("Beam_MomentumZ",1.);   

  tree.put("Beam_MomentumXsigma",0.1);
  tree.put("Beam_MomentumYsigma",0.1);
  tree.put("Beam_MomentumZsigma",0);
  
  tree.put("Target_Size",2);
  tree.put("Target_Size.unit","cm");
  tree.put("Target_Size.unit.cm",1.*cm);

  tree.put("Target_PosX",0);
  tree.put("Target_PosX.unit","cm");
  tree.put("Target_PosX.unit.cm",1.*cm);

  tree.put("Target_PosY",0);
  tree.put("Target_PosY.unit","cm");
  tree.put("Target_PosY.unit.cm",1.*cm);

  tree.put("Target_PosZ",0);
  tree.put("Target_PosZ.unit","cm");
  tree.put("Target_PosZ.unit.cm",1.*cm);

  tree.put("Wasa_ShiftZ",1);
  tree.put("Wasa_ShiftZ.unit","m");
  tree.put("Wasa_ShiftZ.unit.m",1.*m);

  tree.put("Systematic_Shift",0);
  tree.put("Systematic_Shift.unit","cm");
  tree.put("Systematic_Shift.unit.cm",1.*cm);
    
  tree.put("Wasa_Side",0);

  tree.put("HypHI_Si_minR",1);
  tree.put("HypHI_Si_minR.unit","cm");
  tree.put("HypHI_Si_minR.unit.cm",1.*cm);

  tree.put("HypHI_Si_maxR",3);
  tree.put("HypHI_Si_maxR.unit","cm");
  tree.put("HypHI_Si_maxR.unit.cm",1.*cm);

  tree.put("HypHI_EndCap_posZ",40);
  tree.put("HypHI_EndCap_posZ.unit","cm");
  tree.put("HypHI_EndCap_posZ.unit.cm",1.*cm);

  tree.put("HypHI_EndCap_maxR",30);
  tree.put("HypHI_EndCap_maxR.unit","cm");
  tree.put("HypHI_EndCap_maxR.unit.cm",1.*cm);

  tree.put("FRS_FMF2_posZ",2);
  tree.put("FRS_FMF2_posZ.unit","m");
  tree.put("FRS_FMF2_posZ.unit.m",1.*m);

  tree.put("FRS_TR1_posZ",40);
  tree.put("FRS_TR1_posZ.unit","cm");
  tree.put("FRS_TR1_posZ.unit.cm",1.*cm);

  tree.put("FRS_TR2_posZ",70);
  tree.put("FRS_TR2_posZ.unit","cm");
  tree.put("FRS_TR2_posZ.unit.cm",1.*cm);
  
  tree.put("CDS_RelativePosTarget",-0.7);

  tree.put("Output_Namefile","Default_Output.root");

  tree.put("SimpleGeo",1);
  tree.put("Geo","CDS");
  
  tree.put("Physicslist","G4Default_FTFP_BERT");

  tree.put("HyperNuclei_H4L_br_mode1",1.); 
  tree.put("HyperNuclei_H4L_br_mode2",0.); 
  tree.put("HyperNuclei_H4L_br_mode3",0.); 
  tree.put("HyperNuclei_H4L_br_mode4",0.); 
  tree.put("HyperNuclei_H4L_br_mode5",0.); 

  tree.put("HyperNuclei_H4LT12",0.194);
  tree.put("HyperNuclei_H4LT12.unit","ns");
  tree.put("HyperNuclei_H4LT12.unit.ns",1.*ns);
  
  tree.put("HyperNuclei_H3L_br_mode1",0); 
  tree.put("HyperNuclei_H3L_br_mode2",1); 
  tree.put("HyperNuclei_H3L_br_mode3",0); 
  tree.put("HyperNuclei_H3L_br_mode4",0); 
  tree.put("HyperNuclei_H3L_br_mode5",0);

  tree.put("HyperNuclei_H3L_T12",0.246);
  tree.put("HyperNuclei_H3L_T12.unit","ns");
  tree.put("HyperNuclei_H3L_T12.unit.ns",1.*ns);

  tree.put("HyperNuclei_He4L_br_mode1",0); 
  tree.put("HyperNuclei_He4L_br_mode2",1); 
  tree.put("HyperNuclei_He4L_br_mode3",0); 
  tree.put("HyperNuclei_He4L_br_mode4",0); 
  tree.put("HyperNuclei_He4L_br_mode5",0); 
  tree.put("HyperNuclei_He4L_br_mode6",0); 
  tree.put("HyperNuclei_He4L_br_mode7",0);

  tree.put("HyperNuclei_He4L_T12",0.256);
  tree.put("HyperNuclei_He4L_T12.unit","ns");
  tree.put("HyperNuclei_He4L_T12.unit.ns",1.*ns);

  tree.put("HyperNuclei_He5L_br_mode1",1); 
  tree.put("HyperNuclei_He5L_br_mode2",0);

  tree.put("HyperNuclei_He5L_T12",0.256);
  tree.put("HyperNuclei_He5L_T12.unit","ns");
  tree.put("HyperNuclei_He5L_T12.unit.ns",1.*ns);

  tree.put("HyperNuclei_nnL_br_mode1",1);
  tree.put("HyperNuclei_nnL_br_mode2",0);

  tree.put("HyperNuclei_nnL_T12",0.256);
  tree.put("HyperNuclei_nnL_T12.unit","ns");
  tree.put("HyperNuclei_nnL_T12.unit.ns",1.*ns);

  tree.put("DefaultRegionCut",10.);
  tree.put("DefaultRegionCut.unit","mm");
  tree.put("DefaultRegionCut.unit.mm",1.*mm);

  tree.put("DetectorRegionCut",50.);
  tree.put("DetectorRegionCut.unit","mm");
  tree.put("DetectorRegionCut.unit.mm",1.*mm);
    
  tree.put("TargetRegionCut",10.);
  tree.put("TargetRegionCut.unit","mm");
  tree.put("TargetRegionCut.unit.mm",1.*mm);

}
