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
  double zmaxm = (inner_r+(thick/2.))*tan(thetamax);
  std::cout<<" zlength at middle is "<<zmaxm<<std::endl;
  double delzb=zmaxb/nzdiv;
  std::cout<<" z div at bottom "<<delzb<<std::endl;
  double delzt=zmaxt/nzdiv;
  std::cout<<" z div at top "<<delzt<<std::endl;
  double delzm=zmaxm/nzdiv;
  std::cout<<" z div at middle "<<delzm<<std::endl;

  // these refer to different fields in the xml file for this detector
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_barrel( x_det.child( _Unicode(barrel) ) );
  xml_comp_t fX_core( fX_struct.child( _Unicode(core) ) );



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

  sens.setType("calorimeter");


  // detector element for a tower. need to make a volume, a placed volume, and then associate the placed volume to this detector element






  for(int iside=0;iside<2;iside++) {
    double aside = 1.;
    if(iside==1) aside=-1.;


  for(int itower=0;itower<nzdiv;itower++) {
    //for(int itower=nzdiv-1;itower<nzdiv;itower++) {


    std::cout<<"ITOWER is "<<itower<<std::endl;

    //    if((itower==0)||(itower==nzdiv-1)) {

    //each z division has a unique tower
    string t_name;
    if(iside==0 ) {
      t_name = _toString(itower,"towerp%d");
    }
    else {
      t_name = _toString(itower,"towerm%d");
    }
    DetElement tower_det(t_name,det_id);  // detector element for a tower

    // angle for tower at this z division with respect to x-y plane
    double aatheta = -1.*atan((itower*(delzt-delzb))/thick);
    //    aatheta=M_PI/4.;





  /*  ROOT:TGeoTrap
It has 11 parameters: 
the half length in z, 
the 2 polar angles from the centre of the face at low z to that at high z, (polar along x, polar along y)
H1 the half length in y at low z, 
LB1 the half length in x at low z and y low edge, 
LB2 the half length in x at low z and y high edge, 
TH1 the angle w.r.t. the y axis from the centre of low y edge to the centre of the high y edge, and H2, LB2, LH2, TH2, the corresponding quantities at high z
  */
    std::cout<<"making trap for eta "<<itower<<std::endl;
    std::cout<<"half thickness "<<(thick+0.1)/2.<<std::endl;
    std::cout<<"polar angle "<<0.<<std::endl;
    std::cout<<"polar angle "<<aatheta<<" "<<std::endl;
    std::cout<<"half length in y at low z "<<inphil<<" "<<std::endl;
    std::cout<<"half length in x at low z and y low edge "<<delzb<<" "<<std::endl;
    std::cout<<"half length in x at low z and y high edge "<<delzb<<" "<<std::endl;
    std::cout<<"polar angle "<<0.<<" "<<std::endl;
    std::cout<<"half length in y at low z "<<outphil<<" "<<std::endl;
    std::cout<<"half length in x at low z and y low edge "<<delzt<<" "<<std::endl;
    std::cout<<"half length in x at low z and y high edge "<<delzt<<" "<<std::endl;
    std::cout<<"polar angle "<<0.<<std::endl;
    dd4hep::Trap towertrap((thick)/2.,aatheta,0.,inphil/2.,delzb/2.,delzb/2.,0.,outphil/2.,delzt/2.,delzt/2.,0.);


    // needs a hole for fiber
    dd4hep::Tube fiberhole = dd4hep::Tube(0.,fX_core.rmax(),thick);
    // rotationzyx.  first a rotation of an angle phi around the z axis followed by a rotation of an angle theta around the y axis followed by a third rotation of an angle psi around the x axis (phi, theta,psi)
    Transform3D tra(RotationZYX(0,aatheta,0.),Translation3D(0.,0.,0.));



    //if I create a mother that is a brass trapezoid, and make the fiber a daughter, I do not need to make a hole in the brass but if I make a mother that is air and place the brass trapezoid and the fiber separately as daughters to the air mother, then I do need to make a hole in the brass
    dd4hep::SubtractionSolid tower(towertrap,fiberhole,tra);


    dd4hep::Volume towerVol( "tower", tower, description.material(fX_barrel.materialStr()) );
    if(itower==0) std::cout<<"    material is "<<fX_barrel.materialStr()<<std::endl;
    towerVol.setVisAttributes(description, fX_barrel.visStr());



    // fibers

    
    dd4hep::Tube fiber = dd4hep::Tube(0.,fX_core.rmin(),thick/2.);

    std::cout<<" making fiber from "<<fX_core.materialStr()<<std::endl;
    dd4hep::Volume coreVol("core", fiber, description.material(fX_core.materialStr()));

    std::cout<<"fX_core.isSensitive is "<<fX_core.isSensitive()<<std::endl;
    if ( fX_core.isSensitive() ) {
      std::cout<<"setting DRSimple fiber sensitive "<<std::endl;
            coreVol.setSensitiveDetector(sens);
    }


    DetElement afiber(tower_det,"fiber",det_id);
    afiber.setAttributes(description,coreVol,fX_core.regionStr(),fX_core.limitsStr(),fX_core.visStr());



    PlacedVolume fiber_phv = towerVol.placeVolume( coreVol, tra);
    //fiber_pv.addPhysVolID("system",det_id);
    // fiber_pv.addPhysVolID("barrel",0);
    //  fiber_pv.addPhysVolID("side",iside);
    //  fiber_pv.addPhysVolID("ieta",itower);
    //  fiber_pv.addPhysVolID("module",nPhi+1);

    fiber_phv.addPhysVolID("fiber",1);
    afiber.setPlacement(fiber_phv);

    

    // now do the placement

    double mod_x_off = 0.;             
    double mod_y_off = inner_r + thick/2;  
    double mod_z_off= delzm/2.;


    //for (int nPhi = 0; nPhi < 1; nPhi++) {
    std::cout<<"starting phi loop"<<std::endl;
    for (int nPhi = 0; nPhi < nphi; nPhi++) {
      //      if((nPhi%2)==0) {
      double phi=nPhi*delphi;
      //std::cout<<"placing at phi "<<phi<<std::endl;

      double m_pos_x = mod_x_off * std::cos(phi) - mod_y_off * std::sin(phi);
      double m_pos_y = mod_x_off * std::sin(phi) + mod_y_off * std::cos(phi);
      double m_pos_z = aside*(mod_z_off+1.0*(delzm*itower));


      if(nPhi==0) {
	std::cout<<"m_pos_z is "<<m_pos_z<<std::endl;
      }

      //Transform3D tr(RotationZYX(0,phi,M_PI*0.5),Translation3D(m_pos_x,m_pos_y,m_pos_z));
      double zrot=-1.*aside*M_PI*0.5;

      Transform3D tr(RotationZYX(zrot,phi,M_PI*0.5),Translation3D(-m_pos_x,-m_pos_y,m_pos_z));
      PlacedVolume pv = envelope.placeVolume(towerVol,tr);
      pv.addPhysVolID("system",det_id);
      pv.addPhysVolID("barrel",0);
      pv.addPhysVolID("side",iside);
      pv.addPhysVolID("ieta",itower);
      pv.addPhysVolID("iphi",nPhi+1);


      DetElement sd = nPhi==0 ? tower_det : tower_det.clone(t_name+_toString(nPhi,"0%d"));


      sd.setPlacement(pv);
      sdet.add(sd);
      }

      //    }

    //      }


      } // end tower loop

  } // end side loop

  // Set envelope volume attributes.
  envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());








  std::cout<<"exiting DRS creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_DRSimple,create_detector)

DECLARE_DEPRECATED_DETELEMENT(DRsimple,create_detector)
