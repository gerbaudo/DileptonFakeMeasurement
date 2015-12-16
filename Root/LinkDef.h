#include "DileptonFakeMeasurement/MatrixPrediction.h"
#include "DileptonFakeMeasurement/SusySelection.h"
#include "DileptonFakeMeasurement/SusyPlotter.h"
#include "DileptonFakeMeasurement/TightProbability.h"
#include "DileptonFakeMeasurement/MeasureFakeRate2.h"
#include "DileptonFakeMeasurement/myHist.h"
#include "DileptonFakeMeasurement/EffObject.h"
#include "DileptonFakeMeasurement/SsPassFlags.h"
#include "DileptonFakeMeasurement/SelectionRegions.h" 
#include "DileptonFakeMeasurement/TupleMakerObjects.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#pragma link C++ class map<string,string>;

#pragma link C++ class MatrixPrediction;
#pragma link C++ class SusySelection;
#pragma link C++ class SusyPlotter;
#pragma link C++ class TightProbability;
#pragma link C++ class MeasureFakeRate2;
#pragma link C++ class myHist;
#pragma link C++ class EffObject;
#pragma link C++ class EffObject2;
#pragma link C++ struct SsPassFlags;
#pragma link C++ namespace susywh;
#pragma link C++ enum susywh::SelectionRegions;
#pragma link C++ namespace susy;
#pragma link C++ namespace susy::wh;
#pragma link C++ struct susy::wh::FourMom+;
#pragma link C++ struct susy::wh::EventParameters+;
#pragma link C++ class vector<susy::wh::FourMom>+;
#endif
