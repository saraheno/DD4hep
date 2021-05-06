//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
#ifndef EXAMPLES_DDDualCrystal_SRC_DRcaloSiPMHit_H
#define EXAMPLES_DDDualCrystal_SRC_DRcaloSiPMHit_H 1





/// Framework include files
#include "DDG4/Geant4Data.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"


typedef ROOT::Math::XYZVector Position;
typedef ROOT::Math::XYZVector Direction;



namespace ddDRcalo {

   
  class DRcaloSiPMHit : public dd4hep::sim::Geant4Calorimeter::Hit   {

  public:
    int ncerenkov,nscintillator;

  public:
    /// Default constructor
    DRcaloSiPMHit() = default;
    /// Initializing constructor
  DRcaloSiPMHit(const Position& cell_pos):dd4hep::sim::Geant4Calorimeter::Hit(cell_pos),ncerenkov(0),nscintillator(0) {}

    /// Default destructor
    virtual ~DRcaloSiPMHit() = default;
    /// Assignment operator
    //DRcaloSiPMHit& operator=(const DRcaloSiPMHit& c);
  };

  /// Helper to dump data file
  /**
   *  Usage:  
   *  $> root.exe
   *  ....
   *  root [0] gSystem->Load("libDDG4Plugins.so");
   *  root [1] gSystem->Load("libDDG4_MySensDet.so");
   *  root [2] CalVision::Dump::dumpData(<num-ebents>,<file-name>);
   *
   */
  class Dump   {
  public:


  };
}

// CINT configuration
#if defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

/// Define namespaces
#pragma link C++ namespace dd4hep;
#pragma link C++ namespace dd4hep::sim;
#pragma link C++ namespace CalVision;
#pragma link C++ class     CalVision::DRcaloSiPMHit+;
#pragma link C++ class     CalVision::DRcaloSiPMHitDump;
#endif

#endif // EXAMPLES_DDDualCrystal_SRC_DRcaloSiPMHit_H
