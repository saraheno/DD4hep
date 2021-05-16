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
  // material to underly it all
  Material      air       = description.air();


  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is "<<det_name<<std::endl;


  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_dim     = x_det.dimensions();


  int           nphi    = x_dim.numsides();
  std::cout<<" nphi is "<<nphi<<std::endl;
  int           nzdiv    = x_dim.nz();
  std::cout<<" nzdiv is "<<nzdiv<<std::endl;


  double        inner_r   = x_dim.rmin();
  std::cout<<" inner_r is "<<inner_r<<std::endl;
  double        thick   = x_dim.thickness();
  double outer_r = inner_r + thick;
  double delphi=2.*M_PI/nphi;
  std::cout<<" delta phi is "<<delphi<<std::endl;
  double inphil = 2*inner_r*tan(delphi/2.);
  std::cout<<" inner phi chord is "<<inphil<<std::endl;
  double outphil = 2*outer_r*tan(delphi/2.);
  std::cout<<" outer phi chord is "<<outphil<<std::endl;
  std::cout<<" thickness is "<<thick<<std::endl;
  double zmaxb = x_dim.z_length();
  std::cout<<" zlength at bottom is "<<zmaxb<<std::endl;
  double thetamax=atan(zmaxb/inner_r);
  std::cout<<" theta range is "<<thetamax<<std::endl;
  double zmaxt = outer_r*tan(thetamax);
  std::cout<<" zlength at top is "<<zmaxt<<std::endl;
  double delzb=zmaxb/nzdiv;
  std::cout<<" z div at bottom "<<delzb<<std::endl;
  double delzt=zmaxt/nzdiv;
  std::cout<<" z div at top "<<delzt<<std::endl;

  // these refer to different fields in the xml file for this detector
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_barrel( x_det.child( _Unicode(barrel) ) );
  xml_comp_t fX_cladC( fX_struct.child( _Unicode(cladC) ) );
  xml_comp_t fX_coreC( fX_struct.child( _Unicode(coreC) ) );
  xml_comp_t fX_coreS( fX_struct.child( _Unicode(coreS) ) );


  // three structures, volumes, placedvolumes, and detelements
  // volumes need a setVisAttribute
  // DetElements. you can have volumes inherit attrivutesby setting them here
  //              instead of in the volumes
  // placed volumes need a physvolid, a setplacement in a detelement,
  //                and are created with a mother.placevolume
  
  // detector element for entire detector.  
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = description.pickMotherVolume(sdet);


  PolyhedraRegular hedra  (nphi,inner_r,outer_r+tolerance*2e0,zmaxt);
  Volume        envelope  (det_name,hedra,air);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,RotationZYX(0,0,M_PI/nphi));

  env_phv.addPhysVolID("system",det_id);
  env_phv.addPhysVolID("barrel",0);
  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element



  // detector element for a tower. need to make a volume, a placed volume, and then associate the placed volume to this detector element

  DetElement tower_det("tower0",det_id);  // detector element for a tower





  //  for(int itower=0;itower<nzdiv;itower++) {
  for(int itower=0;itower<1;itower++) {

    double aatheta = atan((itower+0.5)*delzb/inner_r);


    // since all towers in a phi ring are identical, make volumes here


  /*  ROOT:TGeoTrap
It has 11 parameters: 
the half length in z, 
the 2 polar angles from the centre of the face at low z to that at high z, 
H1 the half length in y at low z, 
LB1 the half length in x at low z and y low edge, 
LB2 the half length in x at low z and y high edge, 
TH1 the angle w.r.t. the y axis from the centre of low y edge to the centre of the high y edge, and H2, LB2, LH2, TH2, the corresponding quantities at high z
  */
    std::cout<<"making trap "<<std::endl;
    std::cout<<"half thickness "<<(thick+0.1)/2.<<std::endl;
    std::cout<<"polar angle "<<0.<<std::endl;
    std::cout<<"polar angle "<<aatheta<<" "<<std::endl;
    std::cout<<"half length in y at low z "<<inphil<<" "<<std::endl;
    std::cout<<"half length in x at low z and y low edge "<<delzb<<" "<<std::endl;
    std::cout<<"half length in x at low z and y high edge "<<delzb<<" "<<std::endl;
    std::cout<<"polar angle "<<0.<<" "<<std::endl;
    std::cout<<"half length in y at low z"<<outphil<<" "<<std::endl;
    std::cout<<"half length in x at low z and y low edge "<<delzt<<" "<<std::endl;
    std::cout<<"half length in x at low z and y high edge "<<delzt<<" "<<std::endl;
    std::cout<<"polar angle "<<0.<<std::endl;
    dd4hep::Trap tower((thick+0.1)/2.,0.,aatheta,inphil,delzb,delzb,0.,outphil,delzt,delzt,0.);
    dd4hep::Volume towerVol( "tower", tower, description.material(fX_barrel.materialStr()) );
    if(itower==0) std::cout<<"    material is "<<fX_barrel.materialStr()<<std::endl;
    towerVol.setVisAttributes(description, fX_barrel.visStr());


    double mod_x_off = 0.;             
    double mod_y_off = inner_r + thick/2;  


    for (int nPhi = 0; nPhi < nphi; nPhi++) {
    //for (int nPhi = 0; nPhi < 1; nPhi++) {
      double phi=nPhi*delphi;

      double m_pos_x = mod_x_off * std::cos(phi) - mod_y_off * std::sin(phi);
      double m_pos_y = mod_x_off * std::sin(phi) + mod_y_off * std::cos(phi);




      Transform3D tr(RotationZYX(0,phi,M_PI*0.5),Translation3D(-m_pos_x,-m_pos_y,0.));
      PlacedVolume pv = envelope.placeVolume(towerVol,tr);
      pv.addPhysVolID("system",det_id);
      pv.addPhysVolID("barrel",0);

      int imod=nPhi+itower*nphi;
      std::cout<<"making itower nphi phi imod "<<itower<<" "<<nPhi<<" "<<phi<<" "<<imod<<std::endl;
      pv.addPhysVolID("module",imod+1);

	DetElement sd = ((nPhi==0)&&(itower==0)) ? tower_det : tower_det.clone(_toString(imod,"tower0%d"));
      sd.setPlacement(pv);
      sdet.add(sd);


    }


  }



  // Set envelope volume attributes.
  envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());








  std::cout<<"exiting DRS creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_DRSimple,create_detector)

DECLARE_DEPRECATED_DETELEMENT(DRsimple,create_detector)
