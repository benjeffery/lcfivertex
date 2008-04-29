#include "PerEventIPFitter.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>

#include <inc/lciointerface.h>
#include <algo/inc/pereventipfitter.h>
#include <util/inc/memorymanager.h>

#include <vector>
#include <string>

using namespace marlin ;
using namespace lcio;
using namespace vertex_lcfi;

PerEventIPFitterProcessor aPerEventIPFitterProcessor ;

PerEventIPFitterProcessor::PerEventIPFitterProcessor() : Processor("PerEventIPFitterProcessor") {
  
  // modify processor description
  _description = "Per Event IP fitter - trims tracks to reach probabililty threshold" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "InputRPCollection" , 
			      "Name of the ReconstructedParticle collection contains tracks to fit"  ,
			      _InputRPCollectionName ,
			      std::string("EventTracks") ) ;
  registerOutputCollection( lcio::LCIO::VERTEX,
			      "OutputVertexCollection" , 
			      "Name of the Vertex collection of the output ip vertex"  ,
			      _VertexCollectionName ,
			      std::string("IPVertex") ) ;
  FloatVec DefaultPos;
  DefaultPos.push_back(0.0);
  DefaultPos.push_back(0.0);
  DefaultPos.push_back(0.0);
  registerProcessorParameter( "DefaultIPPosition" , 
			      "Manually set default position of the IP vertex (cm)"  ,
			       _DefaultIPPos ,
			      DefaultPos,
			      DefaultPos.size()) ;
  FloatVec DefaultErr;
  DefaultErr.push_back(pow(5.0/1000.0,2.0)); //5 micron err
  DefaultErr.push_back(0.0);
  DefaultErr.push_back(pow(5.0/1000.0,2.0)); //5 micron err
  DefaultErr.push_back(0.0);
  DefaultErr.push_back(0.0);
  DefaultErr.push_back(pow(20.0/1000.0,2.0)); //20 micron err
  registerProcessorParameter( "DefaultIPError" , 
			      "Manually set default error matrix of the primary vertex (cm) (lower symmetric)"  ,
			      _DefaultIPErr,
			      DefaultErr,
			      DefaultErr.size()) ;
  registerProcessorParameter( "ProbabilityThreshold" , 
			      "Tracks are removed until this threshold is reached"  ,
			      _ProbThreshold ,
			      double(0.01)) ;
}


void PerEventIPFitterProcessor::init() { 

  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  //Make the fitter algorithm object and set its parameters
  _IPFitter = new PerEventIPFitter();
  MemoryManager<Algo<Event*,vertex_lcfi::Vertex*> >::Run()->registerObject(_IPFitter);
  
  _IPFitter->setDoubleParameter("ProbThreshold",_ProbThreshold);
}

void PerEventIPFitterProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void PerEventIPFitterProcessor::processEvent( LCEvent * evt ) { 
	
	LCCollection* RPCollection;
	RPCollection = evt->getCollection( _InputRPCollectionName );
	
	//Create an Event with an IP with default parameters
	Vector3 IPPos;
	SymMatrix3x3 IPErr;
	//Make the manulal ip from the params
	//TODO Check length of vectors and throw if wrong
	IPPos.x() = _DefaultIPPos[0];
	IPPos.y() = _DefaultIPPos[1];
	IPPos.z() = _DefaultIPPos[2];
	IPErr(0,0) = _DefaultIPErr[0];
	IPErr(1,0) = _DefaultIPErr[1];
	IPErr(1,1) = _DefaultIPErr[2];
	IPErr(2,0) = _DefaultIPErr[3];
	IPErr(2,1) = _DefaultIPErr[4];
	IPErr(2,2) = _DefaultIPErr[5];
	
	//Create an event with this IP
	vertex_lcfi::Event* MyEvent = new vertex_lcfi::Event(IPPos,IPErr);
	MemoryManager<Event>::Event()->registerObject(MyEvent);
		
	//Create jets from LCIO and add them to the event
	int nRCP = RPCollection->getNumberOfElements()  ;
	//std::cout << nRCP << std::endl;
	for(int i=0; i< nRCP ; i++)
	{
		MyEvent->addTrack(trackFromLCIORP(MyEvent,dynamic_cast<ReconstructedParticle*>(RPCollection->getElementAt(i))) );
	}
	
	//Run IP Fitter
	vertex_lcfi::Vertex* IPResult = _IPFitter->calculateFor(MyEvent);
	
	//Store resulting vertex in the LCIO file
	lcio::Vertex* LCIOIPResult = vertexFromLCFIVertex(IPResult);
	
	std::vector<std::string>::const_iterator it = find(evt->getCollectionNames()->begin(),evt->getCollectionNames()->end(),_VertexCollectionName);
	if (it == evt->getCollectionNames()->end())
	{
		//Doesn't exist so make collection and add
		LCCollection* MyCollection = new LCCollectionVec("Vertex");
		evt->addCollection(MyCollection,_VertexCollectionName);
	}
	evt->getCollection(_VertexCollectionName)->addElement(LCIOIPResult);

	//Clear all objects created for this event
	MetaMemoryManager::Event()->delAllObjects();
	_nEvt ++ ;
}



void PerEventIPFitterProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void PerEventIPFitterProcessor::end(){ 
  
	MetaMemoryManager::Run()->delAllObjects();
   	std::cout << "PerEventIPFitterProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

