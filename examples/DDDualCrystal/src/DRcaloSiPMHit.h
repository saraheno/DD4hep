#ifndef DRcaloSiPMHit_h
#define DRcaloSiPMHit_h 1

/*#include "DRsimInterface.h"*/

#include <vector>
#include <utility>
#include <map>
#include <tuple>


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "DD4hep/Objects.h"
#include "DD4hep/Segmentations.h"





class DRsimInterface {
 public:
  DRsimInterface();
  ~DRsimInterface();

  typedef std::pair<float,float> hitRange;
  typedef std::map<hitRange, int> DRsimTimeStruct;
  typedef std::map<hitRange, int> DRsimWavlenSpectrum;

  struct DRsimSiPMData {
    DRsimSiPMData() {};
    virtual ~DRsimSiPMData() {};

    int count;
    long long int SiPMnum;
    DRsimTimeStruct timeStruct;
    DRsimWavlenSpectrum wavlenSpectrum;
  };

  struct DRsimTowerData {
    DRsimTowerData() {};
    virtual ~DRsimTowerData() {};

    int iTheta;
    int iPhi;
    int numx;
    int numy;
    std::vector<DRsimSiPMData> SiPMs;
  };


  struct DRsimEdepFiberData {
    DRsimEdepFiberData();
    DRsimEdepFiberData(long long int fiberId64);
    DRsimEdepFiberData(long long int fiberId64, float edep, float edepEle, float edepGamma, float edepCharged);
    virtual ~DRsimEdepFiberData() {};

    void accumulate(float edep, float edepEle, float edepGamma, float edepCharged);

    long long int fiberNum;
    float Edep;
    float EdepEle;
    float EdepGamma;
    float EdepCharged;
  };



  struct DRsimEdepData {
    DRsimEdepData();
    DRsimEdepData(int theta, int phi);
    DRsimEdepData(int theta, int phi, float edep, float edepEle, float edepGamma, float edepCharged);
    virtual ~DRsimEdepData() {};

    void accumulate(float edep, float edepEle, float edepGamma, float edepCharged);

    float Edep;
    float EdepEle;
    float EdepGamma;
    float EdepCharged;
    int iTheta;
    int iPhi;
    std::vector<DRsimEdepFiberData> fibers;
  };

  struct DRsimLeakageData {
    DRsimLeakageData() {};
    virtual ~DRsimLeakageData() {};

    float E;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float vt;
    int pdgId;
  };

  struct DRsimGenData {
    DRsimGenData() {};
    virtual ~DRsimGenData() {};

    float E;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float vt;
    int pdgId;
  };

  struct DRsimEventData {
    DRsimEventData() {};
    virtual ~DRsimEventData() {};

    void clear();

    int event_number;
    std::vector<DRsimTowerData> towers;
    std::vector<DRsimEdepData> Edeps;
    std::vector<DRsimLeakageData> leaks;
    std::vector<DRsimGenData> GenPtcs;
  };
};







namespace ddDRcalo {
  class DRcaloSiPMHit : public G4VHit {
  public:

    DRcaloSiPMHit(G4int wavBin, G4int timeBin);
    DRcaloSiPMHit(const DRcaloSiPMHit &right);
    virtual ~DRcaloSiPMHit();

    const DRcaloSiPMHit& operator=(const DRcaloSiPMHit &right);
    G4bool operator==(const DRcaloSiPMHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void* aHit);

    virtual void Draw() {};
    virtual void Print() {};

    void photonCount() { fPhotons++; }
    unsigned long GetPhotonCount() const { return fPhotons; }

    void SetSiPMnum(dd4hep::DDSegmentation::CellID n) { fSiPMnum = n; }
    const dd4hep::DDSegmentation::CellID& GetSiPMnum() const { return fSiPMnum; }

    void CountWavlenSpectrum(DRsimInterface::hitRange range);
    const DRsimInterface::DRsimWavlenSpectrum& GetWavlenSpectrum() const { return fWavlenSpectrum; }

    void CountTimeStruct(DRsimInterface::hitRange range);
    const DRsimInterface::DRsimTimeStruct& GetTimeStruct() const { return fTimeStruct; }

  private:
    dd4hep::DDSegmentation::CellID fSiPMnum;
    unsigned long fPhotons;
    DRsimInterface::DRsimWavlenSpectrum fWavlenSpectrum;
    DRsimInterface::DRsimTimeStruct fTimeStruct;
    G4int fWavBin;
    G4int fTimeBin;
  };

  typedef G4THitsCollection<DRcaloSiPMHit> DRcaloSiPMHitsCollection;
  extern G4ThreadLocal G4Allocator<DRcaloSiPMHit>* DRcaloSiPMHitAllocator;

  inline void* DRcaloSiPMHit::operator new(size_t) {
    if (!DRcaloSiPMHitAllocator) DRcaloSiPMHitAllocator = new G4Allocator<DRcaloSiPMHit>;
    return (void*)DRcaloSiPMHitAllocator->MallocSingle();
  }

  inline void DRcaloSiPMHit::operator delete(void*aHit) {
    DRcaloSiPMHitAllocator->FreeSingle((DRcaloSiPMHit*) aHit);
  }
}

#endif
