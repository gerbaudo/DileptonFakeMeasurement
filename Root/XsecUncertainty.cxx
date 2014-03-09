#include "SusyTest0/XsecUncertainty.h"

#include "SusyTest0/DsidGroups.h"
#include <iostream>
#include <sstream> // std::ostringstream

using susy::wh::XsecUncertainty;

//----------------------------------------------------------
bool XsecUncertainty::determineGroup(const int dsid)
{
    float unc = 1.0;
    const McGroup group = susy::wh::groupFromDsid(dsid);
    switch(group) {
    case kUnknown   : unc = 1.00; break;
    case kTtbar     : unc = 0.5*(13.3+14.5)/252.9; break;
    case kSingleTop : unc = 1.5/22.4; break;
    case kTtbarW    : unc = 0.07/0.23; break;
    case kTtbarZ    : unc = 0.06/0.21; break;
    case kWw        : unc = 0.05; break;
    case kWz        : unc = 0.05; break;
    case kZz        : unc = 0.07; break;
    case kWwjj      : unc = 0.50; break;
    case kTriboson  : unc = 1.00; break;
    case kZjets     : unc = 0.07; break;
    case kHiggs     : unc = 0.30; break;
        // no default, so that the compiler will warn us of un-handled cases
    }
    m_uncertainty = unc;
    bool groupFound(group!=kUnknown);
    return groupFound;
}
//----------------------------------------------------------
float XsecUncertainty::fractionalUncertainty() const
{
    if(!m_keepQuiet && m_group==kUnknown)
        std::cout<<"XsecUncertainty::fractionalUncertainty() Warning! unknown group "<<std::endl;
    return m_uncertainty;
}
//----------------------------------------------------------
std::string XsecUncertainty::McGroup2str(const McGroup g)
{
    std::string name;
    switch(g) {
    case kUnknown   : name = "kUnknown"  ; break;
    case kTtbar     : name = "kTtbar"    ; break;
    case kSingleTop : name = "kSingleTop"; break;
    case kTtbarW    : name = "kTtbarW"   ; break;
    case kTtbarZ    : name = "kTtbarZ"   ; break;
    case kWw        : name = "kWw"       ; break;
    case kWz        : name = "kWz"       ; break;
    case kZz        : name = "kZz"       ; break;
    case kWwjj      : name = "kWwjj"     ; break;
    case kTriboson  : name = "kTriboson" ; break;
    case kZjets     : name = "kZjets"    ; break;
    case kHiggs     : name = "kHiggs"    ; break;
    }
    return name;
}
//----------------------------------------------------------
std::string XsecUncertainty::str() const
{
  std::ostringstream oss;
  oss<<"XsecUncertainty : group '"<<McGroup2str(m_group)<<"' uncertainty "<<m_uncertainty;
  return oss.str();
}
//----------------------------------------------------------
