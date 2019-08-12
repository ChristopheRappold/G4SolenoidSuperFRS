// -------------------------------------------------
// Definition of the G4SolConfig class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef SOLCONFIG_h
#define SOLCONFIG_h

#include <boost/property_tree/ptree.hpp>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

namespace pt = boost::property_tree;

class G4SolConfig
{
public:
  G4SolConfig();
  G4SolConfig(int argc, char** argv);
  ~G4SolConfig();
  void ParseConfig(const std::string& namefile);
  int ProperConf() { return status; }

  template <typename T>
  T Get(const std::string& key) const
  {
    T temp = tree.get<T>(key);
    if(boost::optional<std::string> unit = tree.get_optional<std::string>(key + ".unit"))
      {
        std::string unitName(*unit);
        T unitVal = tree.get<T>(key + ".unit." + unitName);
        return temp * unitVal;
      }
    else
      return temp;
  }
  bool IsAvailable(const std::string& key) const
  {
    size_t count = tree.count(key);
    if(count > 0)
      return true;
    else
      return false;
  }
  void CheckConfig();

private:
  pt::ptree tree;
  int status;

  int ParseCmd(int argc, char** argv);
  void ParseLine(std::stringstream& lineStream, std::vector<std::string>& res);
  void display(const int depth, const pt::ptree& t);
  double GetDimension(const std::string&);
  void SetDefault();
};

template <>
inline std::string G4SolConfig::Get(const std::string& key) const
{
  return tree.get<std::string>(key);
}
template <>
inline boost::optional<std::string> G4SolConfig::Get(const std::string& key) const
{
  return tree.get_optional<std::string>(key);
}

#endif
