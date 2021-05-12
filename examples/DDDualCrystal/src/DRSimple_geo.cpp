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
//
// Specialized generic detector constructor
// 
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"


using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {


  std::cout<<"Creating DRSimple"<<std::endl;


  xml_det_t     x_det     = e;


  static double tolerance = 0e0;

  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is "<<det_name<<std::endl;



  xml_comp_t    x_dim     = x_det.dimensions();


  // look in DDCore/include/Parsers/detail/Dimension.h for arguments


  int           nphi    = x_dim.numsides();
  std::cout<<" nphi is "<<nphi<<std::endl;
  int           neta    = x_dim.nz();
  std::cout<<" neta is "<<neta<<std::endl;


  double        inner_r   = x_dim.rmin();
  std::cout<<" inner_r is "<<inner_r<<std::endl;
  double        thick   = x_dim.thickness();
  std::cout<<" thickness is "<<thick<<std::endl;
  double zmax = x_dim.z_length();
  std::cout<<" zlength is "<<zmax<<std::endl;
  
  

  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = description.pickMotherVolume(sdet);




  std::cout<<"exiting DRS creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_DRSimple,create_detector)

DECLARE_DEPRECATED_DETELEMENT(DRsimple,create_detector)
