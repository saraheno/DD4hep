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

// Framework include files
#include "DualCrystalCalorimeterHit.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"




namespace CalVision {

  class DualCrystalCalorimeterSD {
  public:
    typedef DualCrystalCalorimeterHit Hit;
    // If we need special data to personalize the action, be put it here
    //int mumDeposits = 0;
    //double integratedDeposit = 0;
  };
}

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {

    using namespace CalVision;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //               Geant4SensitiveAction<MyTrackerSD>
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /** \addtogroup Geant4SDActionPlugin
     *
     * @{
     * \package DualCrystalCalorimeterSDAction
     *
     * @}
     */

    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DualCrystalCalorimeterSD>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<DualCrystalCalorimeterSD::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<DualCrystalCalorimeterSD>::process(G4Step* step,G4TouchableHistory* /*hist*/ ) {


      //std::cout<<"entering DualCrystalAction"<<std::endl;


      Geant4StepHandler h(step);
      HitContribution contrib = DualCrystalCalorimeterHit::extractContribution(step);

      Geant4HitCollection*  coll    = collection(m_collectionID);
      VolumeID cell = 0;

      try {
        cell = cellID(step);
      } catch(std::runtime_error &e) {
	std::stringstream out;
        out << std::setprecision(20) << std::scientific;
        out << "ERROR: " << e.what()  << std::endl;
        out << "Position: "
            << "Pre (" << std::setw(24) << step->GetPreStepPoint()->GetPosition() << ") "
            << "Post (" << std::setw(24) << step->GetPostStepPoint()->GetPosition() << ") "
            << std::endl;
        out << "Momentum: "
            << " Pre (" <<std::setw(24) << step->GetPreStepPoint() ->GetMomentum()  << ") "
            << " Post (" <<std::setw(24) << step->GetPostStepPoint()->GetMomentum() << ") "
            << std::endl;

	std::cout << out.str();

        return true;
      }


      DualCrystalCalorimeterHit* hit = coll->findByKey<DualCrystalCalorimeterHit>(cell);
      if ( !hit ) {
        Geant4TouchableHandler handler(step);
	DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        Position global = h.localToGlobal(pos);
        hit = new DualCrystalCalorimeterHit(global);
        hit->cellID = cell;
        coll->add(cell, hit);
        printM2("CREATE hit with deposit:%e MeV  Pos:%8.2f %8.2f %8.2f  %s",
                contrib.deposit,pos.X,pos.Y,pos.Z,handler.path().c_str());
        if ( 0 == hit->cellID )  { // for debugging only!
          hit->cellID = cellID(step);
          except("+++ Invalid CELL ID for hit!");
        }
      } else {
	//	std::cout<<"updating old hit"<<std::endl;
      }


      G4Track * track =  step->GetTrack();



      //photons
      if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
	
	//	std::cout<<"will robinson have photon "<<track->GetCreatorProcess()->G4VProcess::GetProcessName() <<std::endl;
	//std::cout<<" number of cerenkov is "<<hit->ncerenkov<<std::endl;

	
	if ( track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")  {
	  //std::cout<<" found cerenkov photon"<<std::endl;
          hit->ncerenkov+=1;
          track->SetTrackStatus(fStopAndKill);
          return false;
        }
        else {
          //      std::cout<<" why other photon?"<<std::endl;
          hit->nscintillator+=1;
          track->SetTrackStatus(fStopAndKill);
          return false;
        }
	

      }

      else {
	//        std::cout<<" not a photon"<<std::endl;


        hit->energyDeposit += contrib.deposit;
        hit->truth.emplace_back(contrib);

        mark(h.track);
        return true;
      }

	
      return true;

    }
	
      }
}

//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<DualCrystalCalorimeterSD> DualCrystalCalorimeterSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DualCrystalCalorimeterSDAction)
