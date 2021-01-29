#include "TrackEffWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "SystematicUtils.hxx"

//********************************************************************
TrackEffWeight::TrackEffWeight():EventWeightBase(), BinnedParams(std::string(getenv("PIONANALYSISROOT"))+"/data","trackEff", BinnedParams::k1D_EFF_ASSYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h TrackEffWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box){
//********************************************************************

  (void)box;


  // Get the SystBox for this event, and the appropriate selection and branch
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,box.SuccessfulBranch);


  BinnedParamsParams params;

  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Loop over all TrueParticles in the TPC
  for (Int_t itrue=0;itrue< SystBox->nRelevantTrueObjects; itrue++){      
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(SystBox->RelevantTrueObjects[itrue]);            

    //      std::cout << "- " << itrue << " " << truePart->ID << std::endl;
    // Get the tracking efficiency for bin 0 (there is a single bin);
    int index;
    if(!GetBinValues(0.5, params, index))	  continue;

    bool found=false;
    for (Int_t irec=0;irec<(int) SystBox->nRelevantRecObjects; irec++){
      AnaParticlePD* part = static_cast<AnaParticlePD*>(SystBox->RelevantRecObjects[irec]);
      if (!part->TrueObject) continue;
      //      std::cout << "   - " << irec << " " << part->TrueObject->ID << std::endl;

      if(truePart->ID==part->TrueObject->ID){            
        //        truePart->Print();
        found=true;        
        break;
      }
    }
    eventWeight *= systUtils::ComputeEffLikeWeight(found, toy, GetIndex(), index, params);

  }

  std::cout << eventWeight << std::endl;
  return eventWeight;
}


//********************************************************************
bool TrackEffWeight::IsRelevantTrueObject(const AnaEventC& event, const AnaTrueObjectC& trueObj) const{
//********************************************************************

    (void)event;

    const AnaTrueParticleB& truePart = *static_cast<const AnaTrueParticleB*>(&trueObj);

    if(truePart.Charge==0) return false;
    if(truePart.Momentum<0.4) return false;


    return true;
}

//********************************************************************
Int_t TrackEffWeight::GetRelevantTrueObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************


  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}
