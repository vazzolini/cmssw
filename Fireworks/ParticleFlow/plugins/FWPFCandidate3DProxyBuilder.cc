// -*- C++ -*-
//
// Package:     ParticleFlow
// Class  :     FWCandidate3DProxyBuilder
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Colin Bernet
//         Created:  Fri May 28 15:58:19 CEST 2010
// Edited:           sharris, Wed 9 Feb 2011, 17:34
// Edited:           Azzolini,Fri 3 oct 2014
//

// System include files
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TRandom3.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"

// User include files
#include "Fireworks/Core/interface/FWSimpleProxyBuilderTemplate.h"
#include "Fireworks/Core/interface/FWEvePtr.h"
#include "Fireworks/Core/src/CmsShowMain.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/ParticleFlow/interface/setTrackTypePF.h"
#include "Fireworks/Core/interface/FWGeometry.h"
#include "Fireworks/Core/interface/BuilderUtils.h"
#include "TEveCompound.h"
#include "TEveBoxSet.h"
//#include "TEveBox.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "Fireworks/ParticleFlow/plugins/FWPFCandidateWithHitsProxyBuilder.h"
#include "Fireworks/Core/interface/FWEventItemsManager.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Fireworks/Core/interface/fwLog.h"
#include "DataFormats/DetId/interface/DetId.h"

//-----------------------------------------------------------------------------
// FWPFCandidate3DProxyBuilder
//-----------------------------------------------------------------------------

class FWPFCandidate3DProxyBuilder : public FWSimpleProxyBuilderTemplate<reco::PFCandidate>
{
      
public:
  // ---------------- Constructor(s)/Destructor ----------------------
  FWPFCandidate3DProxyBuilder() { myRandom.SetSeed(0); }
  virtual ~FWPFCandidate3DProxyBuilder();
   
  REGISTER_PROXYBUILDER_METHODS();


protected:
  virtual void localModelChanges(const FWModelId& iId, TEveElement* iCompound,
				 FWViewType::EType viewType, const FWViewContext* vc);
private:
  TRandom3 myRandom;

  FWPFCandidate3DProxyBuilder( const FWPFCandidate3DProxyBuilder& );                    // Stop default
  const FWPFCandidate3DProxyBuilder& operator=( const FWPFCandidate3DProxyBuilder& );   // Stop default

  // --------------------- Member Functions --------------------------
  void build( const reco::PFCandidate& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext* ); //default
  void digitColor(TEveBoxSet* bs, float v, float max);
  void digitColor255(TEveBoxSet* bs, float v, float max);
};
//=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_

//______________________________________________________________________________
FWPFCandidate3DProxyBuilder::~FWPFCandidate3DProxyBuilder(){}

//______________________________________________________________________________
void 
FWPFCandidate3DProxyBuilder::build( const reco::PFCandidate& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext* ) //default
{

  //  TEveCompound* comp = createCompound(false,true);  
  TEveElement* comp = &oItemHolder;
  

  const reco::PFCandidate::ElementsInBlocks& elems = iData.elementsInBlocks();

  // --- XRays -------------------------------------
  // ref HGC collections
  edm::Handle<HGCRecHitCollection> m_collectionHGC;
  HGCRecHitCollection::const_iterator rechit;
  edm::InputTag tagEE("HGCalRecHit:HGCEERecHits");
  edm::InputTag tagHEF("HGCalRecHit:HGCHEFRecHits");
  edm::InputTag tagHEB("HGCalRecHit:HGCHEBRecHits");

  const auto* theEvent = item()->getEvent();
  //fwLog(fwlog::kError) << "1 ---------------------------------------------------------------------Event: "<< theEvent->id().event() << "-------------------"<<std::endl; 
  
  // // look for max and min energy cluster --------
  float energy_max=0;
  float energy_maxEE=0;
  float energy_maxHEF=0;
  float energy_maxHEB=0;

  //float energy_maxEvent=0;
  // ---/. XRays -------------------------------------
  
  for( unsigned i = 0 ; i < elems.size(); ++i ) {
    const reco::PFBlockElement& elem = elems[i].first->elements()[elems[i].second];

    switch( elem.type() ) {
    case reco::PFBlockElement::TRACK:
      {
	TEveRecTrack t;
	t.fBeta = 1.;
	t.fP = TEveVector( iData.px(), iData.py(), iData.pz() );
	t.fV = TEveVector( iData.vertex().x(), iData.vertex().y(), iData.vertex().z() );
	t.fSign = iData.charge();
	TEveTrack* trk = new TEveTrack(&t, context().getTrackPropagator() );      
	trk->MakeTrack();      
	fireworks::setTrackTypePF( iData, trk );    
	setupAddElement( trk, comp, true );
      }
      break;
    case reco::PFBlockElement::ECAL:
    case reco::PFBlockElement::HCAL:
    case reco::PFBlockElement::HGC_ECAL:
    case reco::PFBlockElement::HGC_HCALF:
    case reco::PFBlockElement::HGC_HCALB:
      { 
	if( elem.clusterRef().isNull() || !elem.clusterRef().isAvailable() ) {
	  TEveRecTrack t;
	  t.fBeta = 1.;
	  t.fP = TEveVector( iData.px(), iData.py(), iData.pz() );
	  t.fV = TEveVector( iData.vertex().x(), iData.vertex().y(), iData.vertex().z() );
	  t.fSign = iData.charge();
	  TEveTrack* trk = new TEveTrack(&t, context().getTrackPropagator() );      
	  trk->MakeTrack();      
	  fireworks::setTrackTypePF( iData, trk );    
	  setupAddElement( trk, comp, false );
	  continue;
	}

	const std::vector<std::pair<DetId, float> >& clusterDetIds = 
	  elem.clusterRef()->hitsAndFractions();
	TEveBoxSet* boxset = new TEveBoxSet(); // def
	boxset->Reset(TEveBoxSet::kBT_FreeBox, true, clusterDetIds.size()); // def
	boxset->SetMainTransparency(1);       // XRays 

	// --- XRays -------------------------------------
	// --- look for RecHit maximum energy within a cluster, -------------------
	//                 depending on the subdetector         -------------------
	int ta = -1;
	switch( elem.type() ) {
	case reco::PFBlockElement::HGC_ECAL:        
	  theEvent->getByLabel(tagEE,  m_collectionHGC); 
	  energy_maxEE=0; 
	  ta = 0;
	  break;
	case reco::PFBlockElement::HGC_HCALF:
	  theEvent->getByLabel(tagHEF, m_collectionHGC); 
	  energy_maxHEF=0;
	  ta = 1;
	  break;            
	case reco::PFBlockElement::HGC_HCALB:
	  theEvent->getByLabel(tagHEB, m_collectionHGC);
	  energy_maxHEB=0; 
	  ta = 2;
	  break;
	default:
	  break;
	} // elem.type
	
	if (! m_collectionHGC.product()) continue;

	boxset->SetName(Form("BOX_%d", ta));   
        
	for( std::vector<std::pair<DetId, float> >::const_iterator ip = clusterDetIds.begin(), ipEnd = clusterDetIds.end();
	     ip != ipEnd; ++ip )
	  {
	    auto iter = m_collectionHGC->find(ip->first);             
	    if( iter == m_collectionHGC->end() ) continue;
	    Float_t energy  = iter->energy();
            
	    if (energy > energy_max){
	      energy_max = energy;
	    }
            
	  }// /. ip //out of this loop energy_max and energy_min are  the maximum and min energy of the cluster 

	if(ta==0){         energy_maxEE  = energy_max;   energy_max=0;
	} else if (ta==1){ energy_maxHEF = energy_max;   energy_max=0;
	} else if (ta==2){ energy_maxHEB = energy_max;   energy_max=0;
	}

	// --- Got them, a different energy max value for each subdet -----------------------
	// --- /. XRays -------------------------------------

	for( std::vector<std::pair<DetId, float> >::const_iterator it = clusterDetIds.begin(), itEnd = clusterDetIds.end();
	     it != itEnd; ++it) 
	  {
	    const float* corners = item()->getGeom()->getCorners( (*it).first );
	    if( corners == 0 ) {
	      continue;
	    }
	    std::vector<float> pnts(24);    
	    fireworks::energyTower3DCorners(corners, (*it).second, pnts);
	    boxset->AddBox( &pnts[0]);

	    // --- XRays -------------------------------------
	    auto iterA = m_collectionHGC->find(it->first);
	    if( iterA == m_collectionHGC->end() ) continue;

	    Float_t energy2  = iterA->energy();
            
	    if(ta == 0 ){//-------------------------------------------------------------
	      if(energy2 == energy_maxEE){  
		digitColor255(boxset, energy_maxEE, energy_maxEE);
	      }else{
		digitColor(boxset, energy2, energy_maxEE);
	      }
	    }              
	    else if(ta == 1){//---------------------------------------------------------
	      if(energy2 == energy_maxHEF){ 
		digitColor255(boxset, energy_maxHEF, energy_maxHEF);
	      } else{
		digitColor( boxset, energy2, energy_maxHEF);
	      }
	    }
	    else if(ta == 2){//---------------------------------------------------------
	      if (energy2 == energy_maxHEB){
		digitColor255(boxset, energy_maxHEB, energy_maxHEB);
	      } else{
		digitColor(boxset, energy2, energy_maxHEB);
	      }
	    }
	  } // /. it
            // --- /. XRays -------------------------------------
	
	boxset->RefitPlex();                 //Instruct underlying memory allocator to regroup itself into a contiguous memory chunk.
	setupAddElement(boxset,comp, false); //(element, parents, color)

	// /.cases
      }
      break;
    default:
      break;
    } // /. of switch(elem.type) // type of PF candidate 
  } // /. elems.size         // size of the PF candidates
  
  //comp->SetMainColor((unsigned)2.0*myRandom.Uniform(50));
  //setupAddElement( comp, &oItemHolder, false );
}


//______________________________________________________________________________

void FWPFCandidate3DProxyBuilder::digitColor(TEveBoxSet* bs, float v, float max) {

  bs->SetMainTransparency(1);
  // last argument is opacitiy [0-255]
  if (item())
    //  bs->DigitColor(item()->defaultDisplayProperties().color(), 255*v/max);             // def
    //  bs->DigitColor(item()->defaultDisplayProperties().color(), 255*exp( -v/max + 1 )); // idea 1
    //  bs->DigitColor(item()->defaultDisplayProperties().color(), 255*exp( -v/max  ));    // idea 2
    bs->DigitColor(item()->defaultDisplayProperties().color(), 255*(max - v)/ max);        // idea 3 <------
  //  bs->DigitColor(item()->defaultDisplayProperties().color(), 255*exp(1-max/v));      // idea 4
  //  bs->DigitColor(item()->defaultDisplayProperties().color(), 255*pow(v/max,2));      // idea 5
}

//______________________________________________________________________________
//______________________________________________________________________________

void FWPFCandidate3DProxyBuilder::digitColor255(TEveBoxSet* bs, float v, float max) {

  bs->SetMainTransparency(1);
  // last argument is opacitiy [0-255]
  if (item())
    bs->DigitColor(item()->defaultDisplayProperties().color(), 255);                   // def <---
}

//______________________________________________________________________________





void
FWPFCandidate3DProxyBuilder::localModelChanges(const FWModelId& iId, TEveElement* iCompound,
                                               FWViewType::EType viewType, const FWViewContext* vc)
{
  UChar_t rgba[4];
  TEveUtil::ColorFromIdx(item()->defaultDisplayProperties().color() + context().colorManager()->offsetOfLimitedColors(), rgba, item()->defaultDisplayProperties().transparency()); 
  // TEveUtil::ColorFromIdx(item()->defaultDisplayProperties().color(), rgba, item()->defaultDisplayProperties().transparency()); // def
  // TEveUtil::ColorFromIdx(Color_t ci, UChar_t col[4], UChar_t transparency) 
  // this takes color-index (Short_t, can be translated to RGB) and transparency 
  // and converts them to RGBA representation (A ~ alpha value as accepted by GL: 0 - fully transparent, 255 - fully opaque)
  for (TEveElement::List_i it = iCompound->BeginChildren(); it != iCompound->EndChildren(); ++it)
    {
      TEveBoxSet* boxset = dynamic_cast<TEveBoxSet*>(*it);
      if (!boxset) continue;
      boxset->SetMainTransparency(1);
      TEveChunkManager* plex = boxset->GetPlex();
      if (plex->N())
	{
	  for (int atomIdx=0; atomIdx < plex->Size(); ++atomIdx)
	    {
	      TEveDigitSet::DigitBase_t*  digit = boxset->GetDigit(atomIdx);
	      UChar_t* x = (UChar_t*) & digit->fValue;
	      x[0] = rgba[0]; x[1] = rgba[1]; x[2] = rgba[2];
	      //printf("AMT [%d] %d %d %d alpha = %d \n",atomIdx, (int) x[0],(int) x[1], (int)x[2], (int)x[3]);
	    }
	}
    }

}
//    bs->DigitColor(item()->defaultDisplayProperties().color() + context().colorManager()->offsetOfLimitedColors() , 255);
//comp->SetMainColor((unsigned)2.0*myRandom.Uniform(50));

//______________________________________________________________________________
REGISTER_FWPROXYBUILDER(FWPFCandidate3DProxyBuilder, reco::PFCandidate,"PF Candidates", FWViewType::kAll3DBits | FWViewType::kAllRPZBits );
