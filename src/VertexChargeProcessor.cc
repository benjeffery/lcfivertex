#include "VertexChargeProcessor.h"
#include <iostream>
#include <sstream>

#include <EVENT/LCParameters.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/ParticleIDImpl.h>
#include <EVENT/LCFloatVec.h>

#include <inc/event.h>
#include <util/inc/memorymanager.h>
#include <inc/jet.h>
#include <inc/track.h>
#include <inc/algo.h>
#include <inc/vertex.h>
#include <algo/inc/vertexcharge.h>
#include <algo/inc/trackattach.h>
#include <util/inc/helixrep.h>
#include <util/inc/projection.h>
#include <inc/track.h>
#include <inc/lciointerface.h>
#include <algo/inc/twotrackpid.h>

#include <vector>
#include <string>
#include <map>

using namespace marlin ;
using namespace lcio;
using namespace vertex_lcfi;
using vertex_lcfi::util::Projection;
using std::vector;
using std::string;

VertexChargeProcessor aVertexChargeProcessor ;

VertexChargeProcessor::VertexChargeProcessor() : Processor("VertexChargeProcessor") {
    // modify processor description
  _description = "VertexChargeProcessor - takes a set of vertices as a decay chain with its associated jet and calculates the vertex charge." ;
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "JetRPCollection" , 
			      "Name of the ReconstructedParticle collection that represents jets",
			      _JetRPColName ,
			      std::string("Jets") ) ;
  registerInputCollection( lcio::LCIO::VERTEX,
			      "IPVertexCollection" , 
			      "Name of the Vertex collection that contains the primary vertex (Optional)"  ,
			      _IPVertexCollectionName ,
			      std::string("IPVertex") ) ;
  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			   "DecayChainRPCollection" , 
			   "Name of the ReconstructedParticle collection that represents decay chains"  ,
			   _DecayChainRPColName ,
			   std::string("DecayChains") ) ;
  registerOutputCollection( lcio::LCIO::LCFLOATVEC,
			    "VertexChargeCollection" , 
			    "Name of the LCFloatVec Collection that will be created to contain the flavour tag inputs"  ,
			    _VertexChargeCollectionName,
			    "BCharge" ) ;
  registerOptionalParameter( "ChargeAllSecondaryTracks",
			     "Parameter determining whether all tracks from secondary are included in the B-Charge"  ,
			     _ChargeAddAllTracksFromSecondary,
			     bool(true));
  registerOptionalParameter( "ChargeLoDCutmin",
			     "Cut determining the minimum L/D for the B-Charge"  ,
			    _ChargeLoDCutmin,
			    double(0.18 )) ;
  registerOptionalParameter( "ChargeLoDCutmax",
			     "Cut determining the maximum L/D for the B-Charge"  ,
			     _ChargeLoDCutmax,
			     double(2.5)) ;
  registerOptionalParameter( "ChargeCloseapproachCut",
			     "Upper cut on track distance of closest approach to the seed axis for the B-Charge "   ,
			     _ChargeCloseapproachCut,
			     double(1.0)) ;

}

void VertexChargeProcessor::init() 
{ 
	// usually a good idea to
	printParameters() ;
	
	_nRun = 0 ;
	_nEvt = 0 ;

	_Attach = new TrackAttach();
	MemoryManager<Algo<DecayChain*, DecayChain*> > ::Run()->registerObject(_Attach);
	_Attach->setDoubleParameter("AddAllTracksFromSecondary",_ChargeAddAllTracksFromSecondary );
	_Attach->setDoubleParameter("LoDCutmin",_ChargeLoDCutmin );
	_Attach->setDoubleParameter("LoDCutmax",_ChargeLoDCutmax );
	_Attach->setDoubleParameter("CloseapproachCut",_ChargeCloseapproachCut );

	_VertexCharge = new VertexCharge();
	MemoryManager<Algo<DecayChain*,double> >::Run()->registerObject(_VertexCharge);


}

void VertexChargeProcessor::processRunHeader( LCRunHeader* run) { 

	_JetVariableNames.clear();
	_JetVariableNames.push_back("Charge");
	run->parameters().setValues(_VertexChargeCollectionName, _JetVariableNames);
	
	_nRun++ ;
} 

void VertexChargeProcessor::processEvent( LCEvent * evt ) { 


  // this gets called for every event 
  // usually the working horse ...
if( isFirstEvent() ) 
{
	std::cout << "Collections in Event" << std::endl;
	const std::vector<std::string>* names = evt->getCollectionNames();
	for (unsigned int i = 0 ; i<names->size() ; ++i) std::cout << (*names)[i] << std::endl;
}

	LCCollection* JetRPCol = evt->getCollection( _JetRPColName );
	LCCollection* DecayChainRPCol = evt->getCollection( _DecayChainRPColName );
	
	//Find the primary vertex in the event
	LCCollection* VertexCol;
	VertexCol = evt->getCollection( _IPVertexCollectionName );
	
	//Search throught the vertices in this colection to find the primary
	Vector3 IPPos;
	Matrix3x3 IPErr;
	int nVerts = VertexCol->getNumberOfElements()  ;
	bool done = 0;
	for(int i=0; i< nVerts ; i++)
	{
		lcio::Vertex* iVertex = dynamic_cast<lcio::Vertex*>(VertexCol->getElementAt(i));
		if (iVertex->isPrimary())
		{
			IPPos.x() = iVertex->getPosition()[0];
			IPPos.y() = iVertex->getPosition()[1];
			IPPos.z() = iVertex->getPosition()[2];
			IPErr(0,0) = iVertex->getCovMatrix()[0];
			IPErr(1,0) = iVertex->getCovMatrix()[1];
			IPErr(1,1) = iVertex->getCovMatrix()[2];
			IPErr(2,0) = iVertex->getCovMatrix()[3];
			IPErr(2,1) = iVertex->getCovMatrix()[4];
			IPErr(2,2) = iVertex->getCovMatrix()[5];
			done = 1;
		}
		if (done) break;
	}
	if (!done) done = 1;//TODO Throw something
		
	Event* MyEvent = new Event(IPPos,IPErr);
	MemoryManager<Event>::Event()->registerObject(MyEvent);
	
	std::map<Jet*,DecayChain*> DecayChainOf;
	std::map<Jet*,ReconstructedParticle*> LCIORPOf;
	//LCRelationNavigator RelNav(evt->getCollection(_RelationColName));
	//std::cout << evt->getCollection(_RelationColName)->getNumberOfElements() << std::endl;
	//Jet Loop, for each jet add it to the event and its corresponding decay chain to the map 
	int nRCP = JetRPCol->getNumberOfElements()  ;
	
	if(nRCP ==0 ) std::cerr<<"Warning: VertexChargeProcessor.cc:336 : NO jets present "<<std::endl; 
	

	for(int i=0; i< nRCP ; i++)
	{
		ReconstructedParticle* JetRP = dynamic_cast<ReconstructedParticle*>(JetRPCol->getElementAt(i)); 
		Jet* ThisJet = jetFromLCIORP(MyEvent,JetRP);
		LCIORPOf[ThisJet] = JetRP;
		//Assume Jets and DecayChains in same order in LCIO
		DecayChainOf[ThisJet] = decayChainFromLCIORP(ThisJet,dynamic_cast<ReconstructedParticle*>(DecayChainRPCol->getElementAt(i)));
		//Commented Out as we rely on the order of decay chains and jets being the same
		/*//Find the Decay chain RP associated with this jet
		std::cout << JetRPCol->getElementAt(i) <<std::endl;
		LCObjectVec RelatedDecayChains = RelNav.getRelatedFromObjects(JetRPCol->getElementAt(i));
		
		LCRelation* rel = dynamic_cast<LCRelation*>( evt->getCollection(_RelationColName)->getElementAt(0) )  ;
		std::cout << rel->getFrom() << " "<< rel->getTo() << " " << rel->getWeight() << std::endl;
		LCRelation* rel2 = dynamic_cast<LCRelation*>( evt->getCollection(_RelationColName)->getElementAt(1) )  ;
		std::cout << rel2->getFrom() << " "<< rel2->getTo() << " " << rel2->getWeight() << std::endl;
		if (RelatedDecayChains.size() != 1)
		{
		std::cout << RelatedDecayChains.size() << " relations found!!" << std::endl;
		//TODO Throw Something
		}
		*/ 			
	}
	
	//Create the collection to store the result
		LCCollectionVec* OutCollection = new LCCollectionVec("LCFloatVec");
		evt->addCollection(OutCollection,_VertexChargeCollectionName);
		

		//Loop over the jets
	for (vector<Jet*>::const_iterator iJet=MyEvent->jets().begin();iJet != MyEvent->jets().end();++iJet)
	{
		
		LCFloatVec VertexCharge;

				
		VertexCharge.push_back(_VertexCharge->calculateFor(_Attach->calculateFor(DecayChainOf[*iJet])));
		LCFloatVec* OutVec = new LCFloatVec(VertexCharge);
		OutCollection->addElement(OutVec);
		
	}//End iJet Loop
	
	//std::cout << ",";std::cout.flush();
	//Clear all objects created for this event
	MetaMemoryManager::Event()->delAllObjects();
	_nEvt ++ ;
}



void VertexChargeProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void VertexChargeProcessor::end(){ 
	
	MetaMemoryManager::Run()->delAllObjects();
   	std::cout << "VertexChargeProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}
