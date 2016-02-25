#ifndef SOLCONFIG_h
#define SOLCONFIG_h

#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <complex>

namespace pt = boost::property_tree;

class G4SolConfig
{
public:
  G4SolConfig(const std::string& namefile);
  ~G4SolConfig();

  template<typename T>
  T Get(const std::string& key) const
  {
    T temp = tree.get<T>(key);
    if(boost::optional<std::string> unit = tree.get_optional<std::string>(key+".unit"))
      {
	std::string unitName (*unit);
	T unitVal = tree.get<T>(key+".unit."+unitName);
	return temp*unitVal;
      }
    else
      return temp;
  }
    
  void CheckConfig();
private:
  pt::ptree tree;

  void display(const int depth, const pt::ptree& t);
  double GetDimension(const std::string&);
  void SetDefault();
  
};



template<> inline
std::string G4SolConfig::Get(const std::string& key) const
{
  return tree.get<std::string>(key);
}


#endif
