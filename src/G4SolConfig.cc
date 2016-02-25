#include "G4SolConfig.hh"
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>

#include "G4SystemOfUnits.hh"


G4SolConfig::G4SolConfig(const std::string& namefile)
{
  SetDefault();
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
  tree.put("Particle","proton");

  tree.put("Beam_Momentum",2);
  tree.put("Beam_Momentum.unit","GeV");
  tree.put("Beam_Momentum.unit.GeV",1.*GeV);

  tree.put("Beam_MomentumSigma",10.);
  tree.put("Beam_MomentumSigma.unit","MeV");
  tree.put("Beam_MomentumSigma.unit.MeV",1*MeV);
    
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

  tree.put("CDS_RelativePosTarget",-0.7);

  
}
