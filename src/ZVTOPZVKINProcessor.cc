#include "ZVTOPZVKINProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>

#include <util/inc/memorymanager.h>
#include <algo/inc/zvkin.h>
#include <util/inc/matrix.h>
#include <inc/lciointerface.h>

#include <vector>
#include <string>

using namespace marlin ;
using namespace lcio;
using namespace vertex_lcfi;

ZVTOPZVKINProcessor aZVTOPZVKINProcessor ;

ZVTOPZVKINProcessor::ZVTOPZVKINProcessor() : Processor("ZVTOP_ZVKINProcessor") {
  
  // modify processor description
  _description = "ZVTOP_ZVKIN - Kinematic Vertex Reconstruction Algorithm" ;
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "JetRPCollection" , 
			      "Name of the ReconstructedParticle collection that represents jets"  ,
			      _JetRPCollectionName ,
			      std::string("Jets") ) ;
  registerInputCollection( lcio::LCIO::VERTEX,
			      "IPVertexCollection" , 
			      "Name of the Vertex collection that contains the primary vertex (Optional)"  ,
			      _IPVertexCollectionName ,
			      std::string("IPVertex") ) ;
  registerOutputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "DecayChainRPTracksCollectionName" , 
			      "Name of the ReconstructedParticle collection that represents tracks in output decay chains"  ,
			      _DecayChainRPTracksCollectionName,
			      std::string("ZVKINDecayChainRPTracks") ) ;
//  registerOutputCollection( lcio::LCIO::LCRELATION,
//			      "RelationCollection" , 
//			      "Name of the LCRelation collection where the relation between jets and decay chains will be stored"  ,
//			      _RelationCollectionName ,
//			      std::string("JetDecayChainRelations") ) ;
  registerOutputCollection( lcio::LCIO::VERTEX,
			      "VertexCollection" , 
			      "Name of the Vertex collection that contains found vertices"  ,
			      _VertexCollectionName ,
			      std::string("ZVKINVertices") ) ;
  registerOutputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "DecayChainCollectionName" , 
			      "Name of the ReconstructedParticle collection that holds RPs representing output decay chains"  ,
			      _DecayChainCollectionName,
			      std::string("ZVKINDecayChains") ) ;
  registerProcessorParameter( "ManualIPVertex" , 
			      "If false then the primary vertex from VertexCollection is used"  ,
			      _ManualPrimaryVertex ,
			      bool(1));
  FloatVec DefaultPos;
  DefaultPos.push_back(0.0);
  DefaultPos.push_back(0.0);
  DefaultPos.push_back(0.0);
  //_ManualPrimaryVertexPos = DefaultPos;
  registerOptionalParameter( "ManualIPVertexPosition" , 
			      "Manually set position of the primary vertex (cm) - non origin IP not yet fully supported"  ,
			       _ManualPrimaryVertexPos ,
			      DefaultPos,
			      DefaultPos.size()) ;
  FloatVec DefaultErr;
  DefaultErr.push_back(pow(5.0/1000.0,2.0)); //5 micron err
  DefaultErr.push_back(0.0);
  DefaultErr.push_back(pow(5.0/1000.0,2.0)); //5 micron err
  DefaultErr.push_back(0.0);
  DefaultErr.push_back(0.0);
  DefaultErr.push_back(pow(20.0/1000.0,2.0)); //20 micron err
  //_ManualPrimaryVertexErr = DefaultErr;
  registerOptionalParameter( "ManualIPVertexError" , 
			      "Manually set error matrix of the primary vertex (cm) (lower symmetric)"  ,
			      _ManualPrimaryVertexErr,
			      DefaultErr,
			      DefaultErr.size()) ;
  registerOptionalParameter( "MinimumProbability" , 
			      "If a vertex candidate has a probability below this it will not be considered - lower value results in more merging and lower vertex multiplicity"  ,
			      _MinimumProbability,
			      double(1.0/100.0)) ;
  registerOptionalParameter( "InitialGhostWidth" , 
			      "Width in cm of the ghost inital ghosttrack, also the smallest width it is allowed to have"  ,
			      _InitialGhostWidth,
			      double(25.0/1000.0)) ;
  registerOptionalParameter( "MaxChi2Allowed" , 
			      "The ghost track is widened until all forward jet tracks have a chi squared lower than this value"  ,
			      _MaxChi2Allowed,
			      double(1.0)) ;
  registerOptionalParameter( "OutputTrackChi2" , 
			      "If true the chi squared contributions of tracks to vertices is written to LCIO"  ,
			      _OutputTrackChi2,
			      false) ;
}


void ZVTOPZVKINProcessor::init() { 

  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  //Make the ZVKIN algorithm object and set its parameters
  _ZVKIN = new ZVKIN();
  MemoryManager<Algo<Jet*,DecayChain*> >::Run()->registerObject(_ZVKIN);
  
  _ZVKIN->setDoubleParameter("MinimumProbability",_MinimumProbability);
  _ZVKIN->setDoubleParameter("InitialGhostWidth",_InitialGhostWidth);
  _ZVKIN->setDoubleParameter("MaxChi2Allowed",_MaxChi2Allowed);
  _ZVKIN->setStringParameter("AutoJetAxis","TRUE");
  _ZVKIN->setStringParameter("UseEventIP","TRUE");
	
}

void ZVTOPZVKINProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void ZVTOPZVKINProcessor::processEvent( LCEvent * evt ) { 
	//Make Event from 
	LCCollection* JetCollection;
	JetCollection = evt->getCollection( _JetRPCollectionName );
	
	//Create an Event with an IP determined by the parameters or a vertex
	Vector3 IPPos;
	SymMatrix3x3 IPErr;
	if (_ManualPrimaryVertex)
	{
		//Make the manulal ip from the params
		//TODO Check length of vectors and throw if wrong
		IPPos.x() = _ManualPrimaryVertexPos[0];
		IPPos.y() = _ManualPrimaryVertexPos[1];
		IPPos.z() = _ManualPrimaryVertexPos[2];
		IPErr(0,0) = _ManualPrimaryVertexErr[0];
		IPErr(1,0) = _ManualPrimaryVertexErr[1];
		IPErr(1,1) = _ManualPrimaryVertexErr[2];
		IPErr(2,0) = _ManualPrimaryVertexErr[3];
		IPErr(2,1) = _ManualPrimaryVertexErr[4];
		IPErr(2,2) = _ManualPrimaryVertexErr[5];
	}
	else
	{
		//Find the primary vertex in the event
		LCCollection* VertexCol;
		VertexCol = evt->getCollection( _IPVertexCollectionName );
		
		//Search throught the vertices in this colection to find the primary
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
	}
	//Create an event with this IP
	vertex_lcfi::Event* MyEvent = new vertex_lcfi::Event(IPPos,IPErr);
	
	//Create jets from LCIO and add them to the event
	std::vector<std::string>::const_iterator it = find(evt->getCollectionNames()->begin(),evt->getCollectionNames()->end(),_DecayChainCollectionName);
		if (it == evt->getCollectionNames()->end())
		{
			//Doesn't exist so make collection and add
			LCCollection* MyCollection = new LCCollectionVec("ReconstructedParticle");
			evt->addCollection(MyCollection,_DecayChainCollectionName);
		}
	int nRCP = JetCollection->getNumberOfElements()  ;
	for(int i=0; i< nRCP ; i++)
	{
		Jet* MyJet = jetFromLCIORP(MyEvent,dynamic_cast<ReconstructedParticle*>(JetCollection->getElementAt(i)));
	
		//Set any jet depandant parameters
		
		//Run ZVTOP-ZVKIN
		DecayChain* ZVTOPResult = _ZVKIN->calculateFor(MyJet);
		
		//Store resulting decay chain in the LCIO file
		ReconstructedParticle* LCIOZVTOPResult = addDecayChainToLCIOEvent(evt, ZVTOPResult,_VertexCollectionName, _DecayChainRPTracksCollectionName, _OutputTrackChi2);
		
		//Store the RP that holds all the vertexed tracks in LCIO
		evt->getCollection(_DecayChainCollectionName)->addElement(LCIOZVTOPResult);
		
		//Commented out as we just rely on order
		/*//LC Relate DecayChain and Jet
		LCRelation* NewRelation = new LCRelationImpl(LCIOZVTOPResult,JetCollection->getElementAt(i));
		it = find(evt->getCollectionNames()->begin(),evt->getCollectionNames()->end(),_RelationCollectionName);
		if (it == evt->getCollectionNames()->end())
		{
			//Doesn't exist so make collection and add
			LCCollection* MyCollection = new LCCollectionVec("LC_RELATION");
			evt->addCollection(MyCollection,_RelationCollectionName);
		}
		evt->getCollection(_RelationCollectionName)->addElement(NewRelation);
		*/
	}
	//Clear all objects created for this event
	MetaMemoryManager::Event()->delAllObjects();
	_nEvt ++ ;
}



void ZVTOPZVKINProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZVTOPZVKINProcessor::end(){ 
  
	MetaMemoryManager::Run()->delAllObjects();
   	std::cout << "ZVTOPZVKINProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

