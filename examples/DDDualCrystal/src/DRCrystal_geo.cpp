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
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {




  std::cout<<"Creating DRCrystal"<<std::endl;


  xml_det_t     x_det     = e;


  static double tol = 0.001;
  Layering      layering (e);

  // material to underly it all
  Material      air       = description.air();


  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is "<<det_name<<std::endl;


  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_towers  = x_det.staves();
  xml_comp_t    x_dim     = x_det.dimensions();


  int           nphi    = x_dim.numsides();
  std::cout<<" nphi is "<<nphi<<std::endl;
  int           nzdiv    = x_dim.nz();
  std::cout<<" nzdiv is "<<nzdiv<<std::endl;


  double        inner_r   = x_dim.rmin();
  std::cout<<" inner_r is "<<inner_r<<std::endl;
  double        thick   = x_dim.thickness();
  double outer_r = inner_r + thick;
  std::cout<<" outer_r is "<<outer_r<<std::endl;
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






  // three structures, volumes, placedvolumes, and detelements
  // volumes need a setVisAttribute
  // DetElements. you can have volumes inherit attrivutesby setting them here
  //              instead of in the volumes
  // placed volumes need a physvolid, a setplacement in a detelement,
  //                and are created with a mother.placevolume
  
  // detector element for entire detector.  
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = description.pickMotherVolume(sdet);


  //PolyhedraRegular hedra  (nphi,inner_r,outer_r+tol*2e0,zmaxt);
  PolyhedraRegular hedra  (nphi,0.,10*outer_r+tol*2e0,10*zmaxt);
  Volume        envelope  (det_name,hedra,air);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,RotationZYX(0,0,M_PI/nphi));

  env_phv.addPhysVolID("system",det_id);
  env_phv.addPhysVolID("barrel",0);
  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element

  sens.setType("calorimeter");


  // detector element for a tower. need to make a volume, a placed volume, and then associate the placed volume to this detector element






  for(int iside=0;iside<1;iside++) {  // positive and negative z parts of detector
    double aside = 1.;  // to do the reflection for negative z
    if(iside==1) aside=-1.;


    // for(int itower=0;itower<nzdiv;itower++) {
    for(int itower=0;itower<nzdiv;itower++) {
      if((itower==1)||(itower==nzdiv-1)) {

	std::cout<<"ITOWER is "<<itower<<std::endl;




    // angle for tower at this z division with respect to x-y plane

	double aatheta = atan((itower*(delzt-delzb))/thick);
	double tanaatheta = (itower*(delzt-delzb))/thick;


    //if I create a mother that is a brass trapezoid, and make the fiber a daughter, I do not need to make a hole in the brass but if I make a mother that is air and place the brass trapezoid and the fiber separately as daughters to the air mother, then I do need to make a hole in the brass



  /*  ROOT:TGeoTrap
It has 11 parameters: 
the half length in z, 
the 2 polar angles from the centre of the face at low z to that at high z, (polar along x, polar along y)
H1 the half length in y at low z, 
LB1 the half length in x at low z and y low edge, 
LB2 the half length in x at low z and y high edge, 
TH1 the angle w.r.t. the y axis from the centre of low y edge to the centre of the high y edge, and H2, LB2, LH2, TH2, the corresponding quantities at high z
  */
	std::cout<<"   tower envelope "<<std::endl;
	std::cout<<"   making trap for eta "<<itower<<std::endl;
	std::cout<<"   half thickness "<<(thick)/2.<<std::endl;
	std::cout<<"   polar angle "<<aatheta<<" "<<std::endl;
	std::cout<<"   half length in y at low z "<<inphil/2.<<" "<<std::endl;
	std::cout<<"   half length in x at low z and y low edge "<<delzb/2.<<" "<<std::endl;
	std::cout<<"   half length in x at low z and y high edge "<<delzb/2.<<" "<<std::endl;
	std::cout<<"   polar angle "<<0.<<" "<<std::endl;
	std::cout<<"   half length in y at high z "<<outphil/2.<<" "<<std::endl;
	std::cout<<"   half length in x at high z and y low edge "<<delzt/2.<<" "<<std::endl;
	std::cout<<"   half length in x at high z and y high edge "<<delzt/2.<<" "<<std::endl;
	std::cout<<"   polar angle "<<0.<<std::endl;


	double aapsi = itower==0 ? atan(thick/(delzt-delzb)) : atan(thick/(delzt+thick*tan(aatheta)-delzb));
	double tanaapsi =thick/(delzt+thick*tanaatheta-delzb); 

	std::cout<<"   psi angle psi c  at far edge "<<aapsi<<" "<<(M_PI/2.)-aapsi<<std::endl;

    // tower envelope
	dd4hep::Trap towertrap((thick)/2.,aatheta,0.,inphil/2.-0.5*tol,delzb/2.-0.5*tol,delzb/2.-0.5*tol,0.,outphil/2.-0.5*tol,delzt/2.-0.5*tol,delzt/2.-0.5*tol,0.);
	dd4hep::Volume towerVol( "tower", towertrap, air);
  std::cout<<"   tower visstr is "<<x_towers.visStr()<<std::endl;
	towerVol.setVisAttributes(description, x_towers.visStr());
	string t_name = iside==0 ? _toString(itower,"towerp%d") : _toString(itower,"towerm%d");
	DetElement tower_det(t_name,det_id);  // detector element for a tower
	towerVol.setSensitiveDetector(sens);



    // Loop over the sets of layer elements in the detector.
	double r_bottoml  = 0.;
	int l_num = 1;
	for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
	  //xml_coll_t li(x_det,_U(layer));
	  //int ihate=0;
	  // while(ihate>1) {
	  xml_comp_t x_layer = li;
	  int repeat = x_layer.repeat();
      // Loop over number of repeats for this layer.
	  for (int j=0; j<repeat; j++)    {
	    string l_name = _toString(l_num,"layer%d");
	    double l_thickness = layering.layer(l_num-1)->thickness();  // Layer's thickness.

	// find top and bottom lengths at this position and center
        // relative to tower bottom
	    double r_topl=r_bottoml + l_thickness;
	    double r_midl=r_bottoml + l_thickness/2.;


	    double inphill = 2*(r_bottoml+inner_r)*tan(delphi/2.);
	    double outphill = 2*(r_topl+inner_r)*tan(delphi/2.);


	    //double bottoml = itower==0 ? (inner_r+r_bottoml)/tan(aapsi) : delzb+(r_bottoml/tan(aapsi))-r_bottoml*tan(aatheta);

	    double bottoml = itower==0 ? (inner_r+r_bottoml)/tanaapsi : delzb-r_bottoml*tanaatheta+r_bottoml/tanaapsi;

	    double topl = itower==0 ? (inner_r+r_topl)/tanaapsi :  delzb+(r_topl/tanaapsi)-r_topl*tanaatheta;
	    double midl = itower==0 ? (inner_r+r_midl)/tanaapsi : delzb+(r_midl/tanaapsi)-r_midl*tanaatheta;


	    // I hate root.  "origin" is not at 0,0,0 but place at 0,0,0, not "origin"
	    double xmidt=delzb/2.+thick*tanaatheta/2.;
	    double xmidl=bottoml/2.+l_thickness*tanaatheta/2.;
	    double deltal=xmidt-xmidl;
	    double offsetl=r_bottoml*tanaatheta;
	    double offset2=((bottoml+l_thickness/tanaapsi)/2.-xmidl)/2.;
	    double dxmidl=-1.*(deltal-offsetl+offset2);




	    //dxmidl=0.;
	    
	    //if(r_midl>thick/2.) dxmidl=-1.*dxmidl;



	    std::cout<<"  dXMIDL is "<<dxmidl<<" "<<bottoml<<" "<<delzb<<" "<<xmidt<<" "<<xmidl<<" "<<deltal<<" "<<offsetl<<" "<<offset2<<std::endl;
	    //dxmidl=0.;

	    double ymidl=0.;

	    double zmidl=-thick/2.+r_midl+0.5*tol;
	    std::cout<<"will robinsock thick r_midl r_bottoml "<<thick<<" "<<r_midl<<" "<<r_bottoml<<std::endl;

	    Position   l_pos(dxmidl,ymidl,zmidl);      // Position of the layer.
	    

	    std::cout<<"   layer "<<l_num<<std::endl;
	    std::cout<<"      half thickness "<<(l_thickness)/2.<<std::endl;
	    std::cout<<"      polar angle "<<aatheta<<" "<<std::endl;
	    std::cout<<"      half length in y at low z "<<inphill/2.<<" "<<std::endl;
	    std::cout<<"      half length in x at low z and y low edge "<<bottoml/2.<<" "<<std::endl;
	    std::cout<<"      half length in x at low z and y high edge "<<bottoml/2.<<" "<<std::endl;
	    std::cout<<"      polar angle "<<0.<<" "<<std::endl;
	    std::cout<<"      half length in y at high z "<<outphill/2.<<" "<<std::endl;
	    std::cout<<"      half length in x at high z and y low edge "<<topl/2.<<" "<<std::endl;
	    std::cout<<"      half length in x at high z and y high edge "<<topl/2.<<" "<<std::endl;
	    std::cout<<"      polar angle "<<0.<<std::endl;
	


	    dd4hep::Trap l_box((l_thickness)/2.-tol,aatheta,0.,inphill/2.-tol,bottoml/2.-tol,bottoml/2.-tol,0.,outphill/2.-tol,topl/2.-tol,topl/2.-tol,0.);
	    dd4hep::Volume     l_vol(l_name,l_box,air);
	    DetElement layer(tower_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
	    int s_num = 1;
	    double r_bottoms=r_bottoml;
	    double r_bottoms2=0;  //relative to slice bottom
	    for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	      xml_comp_t x_slice = si;
	      string     s_name  = _toString(s_num,"slice%d");
	      double     s_thickness = x_slice.thickness();

	      // this is relative to tower bottom, not layer bottom
	      double r_tops=r_bottoms + s_thickness;
	      double r_mids=r_bottoms + s_thickness/2.;
	      double r_mids2 = r_bottoms2+s_thickness/2.;
	      


	      double inphils = 2*(r_bottoms+inner_r)*tan(delphi/2.);
	      double outphils = 2*(r_tops+inner_r)*tan(delphi/2.);


	      double bottoms = itower==0 ? (inner_r+r_bottoms)/tan(aapsi) : delzb+(r_bottoms/tan(aapsi))-r_bottoms*tan(aatheta);
	      double tops = itower==0 ? (inner_r+r_tops)/tan(aapsi) :  delzb+(r_tops/tan(aapsi))-r_tops*tan(aatheta);
	      double mids = itower==0 ? (inner_r+r_mids)/tan(aapsi) : delzb+(r_mids/tan(aapsi))-r_mids*tan(aatheta);





	      double xmids=bottoms/2.+s_thickness*tanaatheta/2.;
	      double offsets=r_bottoms2*tanaatheta;
	      double offsets2=((bottoms+l_thickness/tanaapsi)/2.-xmids)/2.;
	      double dxmids=-1.*(xmidl-xmids-offsets+offsets2)/2.;



	      double ymids=0.;

	      double zmids=-l_thickness/2.+(r_mids-r_bottoml);



	      std::cout<<"      sublayer "<<s_num<<std::endl;
	      std::cout<<"         half thickness "<<(s_thickness)/2.<<std::endl;
	      std::cout<<"         polar angle "<<aatheta<<" "<<std::endl;
	      std::cout<<"         half length in y at low z "<<inphils/2.<<" "<<std::endl;
	      std::cout<<"         half length in x at low z and y low edge "<<bottoms/2.<<" "<<std::endl;
	      std::cout<<"         half length in x at low z and y high edge "<<bottoms/2.<<" "<<std::endl;
	      std::cout<<"         polar angle "<<0.<<" "<<std::endl;
	      std::cout<<"         half length in y at high z "<<outphils/2.<<" "<<std::endl;
	      std::cout<<"         half length in x at high z and y low edge "<<tops/2.<<" "<<std::endl;
	      std::cout<<"         half length in x at high z and y high edge "<<tops/2.<<" "<<std::endl;
	      std::cout<<"         polar angle "<<0.<<std::endl;

	

	      Position   s_pos(dxmids,ymids,zmids);      // Position of the layer.
	      dd4hep::Trap s_box((s_thickness)/2.-2.*tol,aatheta,0.,inphils/2.-2.*tol,bottoms/2.-2.*tol,bottoms/2.-2.*tol,0.,outphils/2.-2.*tol,tops/2.-2.*tol,tops/2.-2.*tol,0.);


	      dd4hep::Volume     s_vol(s_name,s_box,description.material(x_slice.materialStr()));
	      DetElement slice(layer,s_name,det_id);

	      if ( x_slice.isSensitive() ) {
		s_vol.setSensitiveDetector(sens);
          }
  std::cout<<"          slice visstr is "<<x_slice.visStr()<<std::endl;
	      slice.setAttributes(description,s_vol,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());

          // Slice placement.
	      PlacedVolume slice_phv = l_vol.placeVolume(s_vol,s_pos);
	      slice_phv.addPhysVolID("slice", s_num);
	      slice.setPlacement(slice_phv);
          // Increment Z position of slice.
	      r_bottoms += s_thickness;
	      r_bottoms2 += s_thickness;

          // Increment slice number.
	      ++s_num;
	    }



        // Set region, limitset, and vis of layer.
  std::cout<<" layer visstr is "<<x_layer.visStr()<<std::endl;
	    layer.setAttributes(description,l_vol,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());

	    PlacedVolume layer_phv = towerVol.placeVolume(l_vol,l_pos);
	    layer_phv.addPhysVolID("layer", l_num);
	    layer.setPlacement(layer_phv);
        // Increment to next layer Z position.


	    r_bottoml=r_bottoml+l_thickness;

	    ++l_num;
	  }
	}
      





    

    // now do the placement

	double mod_x_off = 0.;             
	double mod_y_off = inner_r + thick/2;  
	double mod_z_off= 0.;


	//for (int nPhi = 0; nPhi < nphi; nPhi++) {
    for (int nPhi = 0; nPhi < 1; nPhi++) {

      //if((nPhi==0)||(nPhi==(nphi/2))) {
	  //if(nPhi%2==0) {
	  double phi=nPhi*delphi;


	  double m_pos_x = mod_x_off * std::cos(phi) - mod_y_off * std::sin(phi);
	  double m_pos_y = mod_x_off * std::sin(phi) + mod_y_off * std::cos(phi);
	  double m_pos_z = aside*(mod_z_off+1.0*(delzm*itower));


	  if(nPhi==0) {
	    std::cout<<"m_pos_z is "<<m_pos_z<<std::endl;
	  }


	  //double zrot=-1.*aside*M_PI*0.5;
	  double zrot=aside*M_PI*0.5;
	  //double zrot=0.;
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
	  //}
    }
      

    //          }

         }
    } // end tower loop

  } // end side loop

  // Set envelope volume attributes.
  std::cout<<" envelope visstr is "<<x_det.visStr()<<std::endl;
  envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());








  std::cout<<"exiting DRCrystal creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_DRCrystal,create_detector)

DECLARE_DEPRECATED_DETELEMENT(DRCrystal,create_detector)
