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
class HftFiller {
public:
    HftFiller(const std::string &mcid, float sumw);
    ~HftFiller();
    bool init(const std::vector<std::string> &systematics);
    bool close();
    bool fill(size_t systIndex, const susy::wh::kin::DilepVars &v);
private: // rule of three 
    HftFiller(const HftFiller&);
    HftFiller& operator=(const HftFiller&);
private:
    std::string m_mcid;
    float m_sumw;
    std::vector<HistFitterTree*> m_hftTrees;
}; // end HftFiller

} // namespace wh
} // namespace susy

#endif // end include guard
