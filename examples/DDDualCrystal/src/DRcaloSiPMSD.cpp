//#include "DRcaloSiPMSD.h"
#include "DRcaloSiPMHit.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"


#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4VolumeManager.h"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"
#include "DD4hep/DD4hepUnits.h"

#include "GridDRcalo.h"

namespace ddDRcalo {


  class DRcaloSiPMSD {
  public:
    typedef DRcaloSiPMHit Hit;
    // If we need special data to personalize the action, be put it here
    //int mumDeposits = 0;
    //double integratedDeposit = 0;


    dd4hep::DDSegmentation::GridDRcalo* fSeg;
    G4int fHCID;


    
    G4int fWavBin=120;
    G4int fTimeBin=600;
    G4float fWavlenStart=900;
    G4float fWavlenEnd=300;
    G4float fTimeStart=10.;
    G4float fTimeEnd=90.;
    G4float fWavlenStep= (fWavlenStart-fWavlenEnd)/(float)fWavBin;
    G4float fTimeStep= (fTimeEnd-fTimeStart)/(float)fTimeBin;
    /*
    DRsimInterface::hitRange findWavRange(G4double en){
      int i = 0;
      for ( ; i < fWavBin+1; i++) {
	if ( en < wavToE( (fWavlenStart - (float)i*fWavlenStep)*dd4hep::nm ) ) break;
      }

      if (i==0) return std::make_pair(fWavlenStart,99999.);
      else if (i==fWavBin+1) return std::make_pair(0.,fWavlenEnd);

      return std::make_pair( fWavlenStart-(float)i*fWavlenStep, fWavlenStart-(float)(i-1)*fWavlenStep );

    };
    DRsimInterface::hitRange findTimeRange(G4double stepTime){
      int i = 0;
      for ( ; i < fTimeBin+1; i++) {
	if ( stepTime < ( (fTimeStart + (float)i*fTimeStep)*CLHEP::ns ) ) break;
      }

      if (i==0) return std::make_pair(0.,fTimeStart);
      else if (i==fTimeBin+1) return std::make_pair(fTimeEnd,99999.);

      return std::make_pair( fTimeStart+(float)(i-1)*fTimeStep, fTimeStart+(float)i*fTimeStep );

    };
    G4double wavToE(G4double wav) { return dd4hep::h_Planck*dd4hep::c_light/wav; }
    */
  };
}

namespace dd4hep {
  namespace sim {
    using namespace ddDRcalo;


    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DRcaloSiPMSD>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<DRcaloSiPMSD::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<DRcaloSiPMSD>::process(G4Step* step,G4TouchableHistory* /*hist*/ ) {



      // this is Sang-Hyun's more complete code.  need to head back this way
      //
      /*
      if(step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

      auto theTouchable = step->GetPostStepPoint()->GetTouchable();
      
      dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
      dd4hep::VolumeID volID = volMgr.volumeID(theTouchable);

      G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector local = theTouchable->GetHistory()->GetTopTransform().TransformPoint( global );
      dd4hep::Position loc(local.x() * dd4hep::millimeter/CLHEP::millimeter, local.y() * dd4hep::millimeter/CLHEP::millimeter, local.z() * dd4hep::millimeter/CLHEP::millimeter);
      dd4hep::Position glob(global.x() * dd4hep::millimeter/CLHEP::millimeter, global.y() * dd4hep::millimeter/CLHEP::millimeter, global.z() * dd4hep::millimeter/CLHEP::millimeter);

      auto cID = m_userData.fSeg->cellID(loc, glob, volID);
      G4double hitTime = step->GetPostStepPoint()->GetGlobalTime();
      G4double energy = step->GetTrack()->GetTotalEnergy();


      Geant4HitCollection*  coll    = collection(m_collectionID);


      DualCrystalCalorimeterHit* hit = coll->findByKey<DualCrystalCalorimeterHit>(cell);

      if (!hit) {
	hit = new DRcaloSiPMHit(fWavBin,fTimeBin);
	hit->SetSiPMnum(cID);

	coll->add(hit);
      }


  hit->photonCount();

  DRsimInterface::hitRange wavRange = m_userData.findWavRange(energy);
  hit->CountWavlenSpectrum(wavRange);

  DRsimInterface::hitRange timeRange = m_userData.findTimeRange(hitTime);
  hit->CountTimeStruct(timeRange);


  return true;
      */



      Geant4StepHandler h(step);
      HitContribution contrib = DRcaloSiPMHit::extractContribution(step);

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


      DRcaloSiPMHit* hit = coll->findByKey<DRcaloSiPMHit>(cell);
      if ( !hit ) {
        Geant4TouchableHandler handler(step);
	DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        Position global = h.localToGlobal(pos);
        hit = new DRcaloSiPMHit(global);
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
      G4int TrPDGid = track->GetDefinition()->GetPDGEncoding();

      std::cout<<"will robinson have track "<<track->GetCreatorProcess()->G4VProcess::GetProcessName()<<" pid "<<TrPDGid <<std::endl;




      //photons
      if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
	
	std::cout<<"will robinson have photon "<<track->GetCreatorProcess()->G4VProcess::GetProcessName() <<std::endl;
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
	std::cout<<" not a photon"<<std::endl;


        hit->energyDeposit += contrib.deposit;
        hit->truth.emplace_back(contrib);

        mark(h.track);
        return true;
      }

	
      return true;




}




  } // namespace sim
}  // namespace dd4hep

//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<DRcaloSiPMSD> DRcaloSiPMSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DRcaloSiPMSDAction)


