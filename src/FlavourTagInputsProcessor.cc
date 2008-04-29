#include "FlavourTagInputsProcessor.h"
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
#include <algo/inc/vertexmomentum.h>
#include <algo/inc/vertexmass.h>
#include <algo/inc/vertexmultiplicity.h>
#include <algo/inc/decaysignificance.h>
#include <algo/inc/paramsignificance.h>
#include <algo/inc/trackattach.h>
#include <algo/inc/jointprob.h>
#include <algo/inc/secondvertexprob.h>
#include <util/inc/vector3.h>
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
using vertex_lcfi::SignificanceType;
using vertex_lcfi::DecaySignificanceType;
using std::vector;
using std::string;

FlavourTagInputsProcessor aFlavourTagInputsProcessor ;

FlavourTagInputsProcessor::FlavourTagInputsProcessor() : Processor("FlavourTagInputsProcessor") {
    // modify processor description
  _description = "FlavourTagInputsProcessor - takes a set of vertices as a decay chain with its associated jet and calculates flavour tag inputs stroring them in the Jet RP's pid" ;
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "JetRPCollection" , 
			      "Name of the ReconstructedParticle collection that represents jets",
			      _JetRPColName ,
			      std::string("Jets") ) ;
  registerInputCollection( lcio::LCIO::VERTEX,
			      "IPVertexCollection" , 
			      "Name of the Vertex collection that contains the primary vertex"  ,
			      _IPVertexCollectionName ,
			      std::string("IPVertex") ) ;
  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			   "DecayChainRPCollection" , 
			   "Name of the ReconstructedParticle collection that represents decay chains"  ,
			   _DecayChainRPColName ,
			   std::string("DecayChains") ) ;
  registerOutputCollection( lcio::LCIO::LCFLOATVEC,
			    "FlavourTagInputsCollection" , 
			    "Name of the LCFloatVec Collection that will be created to contain the flavour tag inputs"  ,
			    _FlavourTagInputsCollectionName,
			    "FlavourTagInputs" ) ;

  registerOptionalParameter( "VertexMassMaxMomentumAngle",
			     "Upper cut on angle between momentum of vertex and the vertex axis",
			     _VertexMassMaxMomentumAngle,
			     double(3)) ;
  registerOptionalParameter( "VertexMassMaxKinematicCorrectionSigma",
			     "Maximum Sigma (based on error matrix) that the vertex axis can move when kinematic correction is applied"  ,
			     _VertexMassMaxKinematicCorrectionSigma,
			     double(2)) ;
  registerOptionalParameter( "VertexMassMaxMomentumCorrection",
			     "Maximum factor, by which vertex mass can be corrected"  ,
			     _VertexMassMaxMomentumCorrection,
			     double(2)) ;


 registerOptionalParameter( "TrackAttachAllSecondaryTracks",
			    "Parameter determining whether all tracks from secondary are included in the track attachment"  ,
			    _TrackAttachAddAllTracksFromSecondary,
			    bool(false)) ;

 registerOptionalParameter( "TrackAttachLoDCutmin",
			     "Cut determining the minimum L/D for the track attachment"  ,
			    _TrackAttachLoDCutmin,
			     double(0.18)) ;
 registerOptionalParameter( "TrackAttachLoDCutmax",
			     "Cut determining the maximum L/D for the track attachment"  ,
			    _TrackAttachLoDCutmax,
			     double(2.5)) ;
 registerOptionalParameter( "TrackAttachCloseapproachCut",
			     "Upper cut on track distance of closest approach to the seed axis for the track attachment"  ,
			    _TrackAttachCloseapproachCut,
			     double(1.0)) ;

 registerOptionalParameter( "SecondVertexProbChisquarecut",
			    "Cut on the Chi Squared of the seed vertex",
			    _SecondVertexProbChisquarecut,
			    double(20.0)) ;
 registerOptionalParameter( "SecondVertexNtrackscut",
			    "Cut on the minimum number of tracks in the seed vertex.", 
			    _SecondVertexNtrackscut,
			    double(1.0)); 

 registerOptionalParameter( "LayersHit",
			    "Momentum cuts will be applied on number of LayersHit and LayersHit minus one",
			    _LayersHit,
			    double(5.0));
 registerOptionalParameter( "AllbutOneLayersMomentumCut",
			    "Cut on the minimum momentum if track hits LayersHit minus one", 
			    _AllbutOneLayersMomentumCut,
			    double(2.0)); 
 registerOptionalParameter( "AllLayersMomentumCut",
			    "Cut on the minimum momentum if track hits LayersHit", 
			    _AllLayersMomentumCut,
			    double(1.0));    

 registerOptionalParameter( "PIDMaxGammaMass",
			    "Cut on the upper limit of the photon candidate mass",
			    _PIDMaxGammaMass,
			    double(0.02));
 registerOptionalParameter( "PIDMinKsMass",
			    "Cut on the lower limit of the Ks candidate mass",
			    _PIDMinKsMass,
			    double(0.475));
 registerOptionalParameter( "PIDMaxKsMass",
			    "Cut on the upper limit of the Ks candidate mass",
			    _PIDMaxKsMass,
			    double(0.525));
 registerOptionalParameter( "PIDChi2Cut",
			    "Cut on the Chi squared of the two tracks beinig in the same vertex.",
			    _PIDChi2Cut,
			    double(6.63));
 registerOptionalParameter( "PIDRPhiCut",
			    "Cut on the maximum RPhi of the Ks/gamma decay vertex candidate",
			    _PIDRPhiCut,
			    double(20));
 registerOptionalParameter( "PIDSignificanceCut",
			    "Cut on the minimum RPhi significance of the tracks", 
			    _PIDSignificanceCut,
			    double(3));

   registerOptionalParameter( "JProbMaxD0Significance",
			      "Upper Cut on the maximum value of d0 significance",
			      _JProbMaxD0Significance,
			      double(200));
   registerOptionalParameter( "JProbMaxD0andZ0",
			      "Upper Cut on the maximum value of d0 and of z0", 
			      _JProbMaxD0andZ0,
			      double(5));

   FloatVec temp;
   temp.push_back(1.01313412);
   temp.push_back(0.0246350896);
   temp.push_back(0.102197811);
   temp.push_back(0.0411203019);
   temp.push_back(0.0157710761);
   registerOptionalParameter( "JProbResolutionParameterRphi",
			      "Standard deviations of the impact parameters in RPhi plane", 
			      _JProbResolutionParameterRphi,
			      temp,
			      temp.size());
   temp.clear();
   temp.push_back(1.01629865);
   temp.push_back(0.0271386635);
   temp.push_back(0.0948112309);
   temp.push_back(0.0410759225);
   temp.push_back(0.0148685882);
   registerOptionalParameter( "JProbResolutionParameterZ",
			      "Standard deviations of the impact parameters in Z direction", 
			      _JProbResolutionParameterZ,
			      temp,
			      temp.size());
}

void FlavourTagInputsProcessor::init() 
{ 
	// usually a good idea to
	printParameters() ;
	
	_nRun = 0 ;
	_nEvt = 0 ;
	
	_VertexMomentum = new VertexMomentum();
	MemoryManager<Algo<DecayChain*,double> >::Run()->registerObject(_VertexMomentum);
	

	_VertexMass = new VertexMass();
	MemoryManager<Algo<DecayChain*,double> >::Run()->registerObject(_VertexMass);  
	_VertexMass->setDoubleParameter("MaxMomentumAngle",_VertexMassMaxMomentumAngle);
	_VertexMass->setDoubleParameter("MaxKinematicCorrectionSigma",_VertexMassMaxKinematicCorrectionSigma);
	_VertexMass->setDoubleParameter("MaxMomentumCorrection",_VertexMassMaxMomentumCorrection);

	_VerticesTrackMultiplicity = new VertexMultiplicity();
	MemoryManager<Algo<DecayChain*,int> >::Run()->registerObject(_VerticesTrackMultiplicity);
	
	_VertexDecaySignificance = new VertexDecaySignificance();
	MemoryManager<Algo<DecayChain*, std::map<DecaySignificanceType,double > > >::Run()->registerObject(_VertexDecaySignificance);
	
	_TrackAttach = new TrackAttach();
	MemoryManager<Algo<DecayChain*, DecayChain*> > ::Run()->registerObject(_TrackAttach);
	_TrackAttach->setDoubleParameter("AddAllTracksFromSecondary",_TrackAttachAddAllTracksFromSecondary );
	_TrackAttach->setDoubleParameter("LoDCutmin",_TrackAttachLoDCutmin );
	_TrackAttach->setDoubleParameter("LoDCutmax",_TrackAttachLoDCutmax );
	_TrackAttach->setDoubleParameter("CloseapproachCut",_TrackAttachCloseapproachCut );
		
	_SecVertexProb = new SecVertexProb();
	MemoryManager<Algo<DecayChain*, double > > ::Run()->registerObject(_SecVertexProb);
	_SecVertexProb->setDoubleParameter("Chisquarecut",_SecondVertexProbChisquarecut);
	_SecVertexProb->setDoubleParameter("Ntrackscut",_SecondVertexNtrackscut);
			
	_ParameterSignificance = new  ParameterSignificance();
	MemoryManager<Algo<Jet*, std::map<SignificanceType,double > > >::Run()->registerObject( _ParameterSignificance);
	_ParameterSignificance->setDoubleParameter("LayersHit",_LayersHit);
	_ParameterSignificance->setDoubleParameter("AllbutOneLayersMomentumCut",_AllbutOneLayersMomentumCut);
	_ParameterSignificance->setDoubleParameter("AllLayersMomentumCut",_AllLayersMomentumCut);

	_JointProb = new  JointProb();
	MemoryManager<Algo<Jet*, std::map<Projection,double > > >::Run()->registerObject( _JointProb);
	_JointProb->setDoubleParameter("MaxD0Significance",_JProbMaxD0Significance);
	_JointProb->setDoubleParameter("MaxD0andZ0",_JProbMaxD0andZ0);
	std::vector<double> temp;  
	/*temp.push_back(_JProbResolutionParameterRphi[0]);
	temp.push_back(_JProbResolutionParameterRphi[1]);
	temp.push_back(_JProbResolutionParameterRphi[2]);
	temp.push_back(_JProbResolutionParameterRphi[3]);
	temp.push_back(_JProbResolutionParameterRphi[4]);
	_JointProb->setPointerParameter("ResolutionParameterRphi", &temp);
	std::vector<double> temp2;  
	temp2.push_back(_JProbResolutionParameterZ[0]);
	temp2.push_back(_JProbResolutionParameterZ[1]);
	temp2.push_back(_JProbResolutionParameterZ[2]);
	temp2.push_back(_JProbResolutionParameterZ[3]);
	temp2.push_back(_JProbResolutionParameterZ[4]);
	_JointProb->setPointerParameter("ResolutionParameterZ",  &temp2);
	*/
	_TwoTrackPID = new TwoTrackPid();
 	MemoryManager<Algo<Jet*, std::map<PidCutType,vector< vertex_lcfi::Track* > > > >::Run()->registerObject( _TwoTrackPID);
	_TwoTrackPID->setDoubleParameter("MaxGammaMass",_PIDMaxGammaMass);
	_TwoTrackPID->setDoubleParameter("MinKsMass",_PIDMinKsMass);
	_TwoTrackPID->setDoubleParameter("MaxKsMass",_PIDMaxKsMass);
	_TwoTrackPID->setDoubleParameter("Chi2Cut",_PIDChi2Cut);
	_TwoTrackPID->setDoubleParameter("RPhiCut",_PIDRPhiCut);
	_TwoTrackPID->setDoubleParameter("SignificanceCut",_PIDSignificanceCut);

}

void FlavourTagInputsProcessor::processRunHeader( LCRunHeader* run) { 

	_JetVariableNames.clear();
	_JetVariableNames.push_back("JointProbRPhi");
	_JetVariableNames.push_back("JointProbZ");
	//	_JetVariableNames.push_back("JointProb3D");
	_JetVariableNames.push_back("D0Significance1");
	_JetVariableNames.push_back("D0Significance2");
	_JetVariableNames.push_back("Z0Significance1");
	_JetVariableNames.push_back("Z0Significance2");
	_JetVariableNames.push_back("Momentum1");
	_JetVariableNames.push_back("Momentum2");
	_JetVariableNames.push_back("NumTracksInVertices");
	_JetVariableNames.push_back("DecayLength");
	_JetVariableNames.push_back("DecayLengthSignificance");
	_JetVariableNames.push_back("RawMomentum");
	_JetVariableNames.push_back("PTCorrectedMass");
	_JetVariableNames.push_back("SecondaryVertexProbability");
	_JetVariableNames.push_back("NumVertices");
	_JetVariableNames.push_back("DecayLength(SeedToIP)");

	
	run->parameters().setValues(_FlavourTagInputsCollectionName, _JetVariableNames);
	
	_nRun++ ;
} 

void FlavourTagInputsProcessor::processEvent( LCEvent * evt ) { 

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
	
	if(nRCP ==0 ) std::cerr<<"Warning: FlavourTagInputsProcessor.cc:336 : NO jets present "<<std::endl; 
	
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
		evt->addCollection(OutCollection,_FlavourTagInputsCollectionName);
	
	//Loop over the jets
	for (vector<Jet*>::const_iterator iJet=MyEvent->jets().begin();iJet != MyEvent->jets().end();++iJet)
	{
		LCFloatVec FlavourTagInputs;
		
		//Probability that all tracks consistant with IP
		std::map<Projection,double> JointProb;
		
		JointProb  = _JointProb->calculateFor(*iJet);
		FlavourTagInputs.push_back(JointProb[RPhi]);
		FlavourTagInputs.push_back(JointProb[Z]);
		//FlavourTagInputs.push_back(JointProb[ThreeD]);
		
		//D0, Z0 significances and momentum of the two most D0 significant tracks
		//First make a cut based on particle pid for this input
		//TODO - Clean up
		std::map<PidCutType, vector<vertex_lcfi::Track*> >* PIDCutTracks = new std::map<PidCutType	, vector<vertex_lcfi::Track*> >();
		MemoryManager<std::map<PidCutType, vector<vertex_lcfi::Track*> > > ::Event()->registerObject(PIDCutTracks);
		*PIDCutTracks = _TwoTrackPID->calculateFor(*iJet);
		std::map<SignificanceType,double> ParSignificance;
		_ParameterSignificance->setPointerParameter( "TwoTrackPidCut", PIDCutTracks);
		ParSignificance  = _ParameterSignificance->calculateFor(*iJet);
		FlavourTagInputs.push_back(ParSignificance[D0SigTrack1]);
		FlavourTagInputs.push_back(ParSignificance[D0SigTrack2]);
		FlavourTagInputs.push_back(ParSignificance[Z0SigTrack1]);
		FlavourTagInputs.push_back(ParSignificance[Z0SigTrack2]);
		FlavourTagInputs.push_back(ParSignificance[MomentumTrack1]);
		FlavourTagInputs.push_back(ParSignificance[MomentumTrack2]);
		
		//Num Tracks in secondary and upwards vertices
		FlavourTagInputs.push_back(_VerticesTrackMultiplicity->calculateFor(DecayChainOf[*iJet]))  ;
		
		//Decay Length and Significance of most significant vertex
		std::map<DecaySignificanceType,double> DecaySignificance;
		DecaySignificance  = _VertexDecaySignificance->calculateFor(DecayChainOf[*iJet]);
		FlavourTagInputs.push_back(DecaySignificance[Distance]);
		FlavourTagInputs.push_back(DecaySignificance[Significance]);
		//Using cuts attach tracks that were not associated to the decay by vertexing
	

		DecayChain* AttachedTracksChain = _TrackAttach->calculateFor(DecayChainOf[*iJet]);
		
		        //Sum momentum of all tracks in decay chain (vertexed and attached)
			FlavourTagInputs.push_back(_VertexMomentum->calculateFor(AttachedTracksChain));
			//Vertex momentum corrected mass
			FlavourTagInputs.push_back(_VertexMass->calculateFor(AttachedTracksChain));
			//Probability of all tracks in decay chain belonging to one vertex
			FlavourTagInputs.push_back(_SecVertexProb->calculateFor(AttachedTracksChain));
			
		//Num Vertices in the vertexing result 
		FlavourTagInputs.push_back(DecayChainOf[*iJet]->vertices().size());
		//std::cout << DecayChainOf[*iJet]->vertices().size();
		//Extra Decay length from seed vertex (last vertex) to IP (IP at Origin for now)
		//TODO De-obfuscate and upgrade to moveable IP
		FlavourTagInputs.push_back((*(--(DecayChainOf[*iJet]->vertices().end())))->position().mag());
		
		LCFloatVec* OutVec = new LCFloatVec(FlavourTagInputs);
		OutCollection->addElement(OutVec);
		
	}//End iJet Loop
	
	//std::cout << ",";std::cout.flush();
	//Clear all objects created for this event
	MetaMemoryManager::Event()->delAllObjects();
	_nEvt ++ ;
}



void FlavourTagInputsProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FlavourTagInputsProcessor::end(){ 
	
	MetaMemoryManager::Run()->delAllObjects();
   	std::cout << "FlavourTagInputsProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}
