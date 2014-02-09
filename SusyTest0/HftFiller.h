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
    HftFiller();
    ~HftFiller();
    bool init(const std::string &mcid, const std::vector<std::string> &systematics);
    bool close(float sumw);
    bool fill(size_t systIndex, const susy::wh::kin::DilepVars &v);
    size_t nTrees() const { return m_hftTrees.size(); }
private: // rule of three 
    HftFiller(const HftFiller&);
    HftFiller& operator=(const HftFiller&);
private:
    std::vector<HistFitterTree*> m_hftTrees;
}; // end HftFiller

} // namespace wh
} // namespace susy

#endif // end include guard
