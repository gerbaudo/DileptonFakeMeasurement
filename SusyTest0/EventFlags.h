// emacs -*- C++ -*-
#ifndef SUSY_WH_EVENTFLAGS_H
#define SUSY_WH_EVENTFLAGS_H

#include <string>

namespace susy
{
namespace wh
{
/*!
  A holder to pass around the event-level flags; true = pass (including vetoes)
  
  davide.gerbaudo@gmail.com
  Feb 2014
*/
struct EventFlags {
  EventFlags() { reset(); }
  void reset() {
   grl = larErr = tileErr = ttcVeto = goodVtx = tileTrip = lAr = false;
   badJet = deadRegions = badMuon = cosmicMuon = hfor = ge2blep = eq2blep = mllMin = false;
  }
  bool grl, larErr, tileErr, ttcVeto, goodVtx, tileTrip, lAr;
  bool badJet, deadRegions, badMuon, cosmicMuon, hfor, ge2blep, eq2blep, mllMin;
  bool allTrue() {
    return  (grl && larErr && tileErr && ttcVeto && goodVtx && tileTrip && lAr &&
             badJet && deadRegions && badMuon && cosmicMuon && hfor && ge2blep && eq2blep && mllMin);
  }
  bool failAny() { return !allTrue(); }
  std::string str() const;
};

} // namespace wh
} // namespace susy

#endif // end include guard
