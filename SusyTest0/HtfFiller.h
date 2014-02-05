// emacs -*- C++ -*-
#ifndef SUSY_WH_HTFFILLER_H
#define SUSY_WH_HTFFILLER_H

#include <string>
#include <vector>

class HistFitterTree;
namespace susy{ namespace wh { namespace kin { class DilepVars; } } }

namespace susy
{
namespace wh
{
/*!
  A class to feed the wh 2l variables to HistFitterTree
  
  davide.gerbaudo@gmail.com
  Jan 2014
*/
class HtfFiller {
public:
    HtfFiller(const std::string &mcid, float sumw);
    ~HtfFiller();
    bool init(const std::vector<std::string> &systematics);
    bool close();
    bool fill(size_t systIndex, const susy::wh::kin::DilepVars &v);
private: // rule of three 
    HtfFiller(const HtfFiller&);
    HtfFiller& operator=(const HtfFiller&);
private:
    std::string m_mcid;
    float m_sumw;
    std::vector<HistFitterTree*> m_hftTrees;
}; // end HtfFiller

} // namespace wh
} // namespace susy

#endif // end include guard
