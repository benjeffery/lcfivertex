#include "RPCutProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCObject.h>
#include "../vertex_lcfi/inc/lciointerface.h"
#include "util/inc/memorymanager.h"

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gearimpl/Vector3D.h>

#include <vector>
#include <string>

using namespace marlin ;
using namespace lcio;

RPCutProcessor aRPCutProcessor ;

RPCutProcessor::RPCutProcessor() : Processor("RPCutProcessor") {
  
  // modify processor description
  _description = "RPCutProcessor - cuts RPs based on several criteria removing those that have no tracks" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "InputRCPCollection" , 
			      "Name of the ReconstructedParticle collection which will be cut"  ,
			      _InRCPColName ,
			      std::string("Jets") ) ;
  registerOutputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "OutputRCPCollection" , 
			      "Name of the output ReconstructedParticle collection when WriteNewCollection is true"  ,
			      _OutRCPColName ,
			      std::string("CutJets") ) ;
  registerProcessorParameter( "WriteNewCollection" , 
			      "If true writes the cut jet to a new collection (OutputJetRCPCollection) if false overwrites input"  ,
			      _WriteNewCollection,
			      bool(1) ) ;
  registerProcessorParameter( "SubParticleLists" , 
			      "If true cuts tracks from the particle lists of particles in InputRCPCollection, if false just cuts particles from InputRCPCollection"  ,
			      _SubParticleLists,
			      bool(1) ) ;
  
  registerOptionalParameter( "a1_Chi2OverDOFEnable" , 
			      "Enable a cut on the value of each tracks chi squared over degrees of freedom"  ,
			      _Chi2OverDOFEnable,
			      bool(0) ) ;
  registerOptionalParameter( "a2_Chi2OverDOFCutLowerThan" , 
			      "If true values lower than the cut value will be cut, if false values higher will be cut"  ,
			      _Chi2OverDOFCutLowerThan,
			      bool(0) ) ;
  registerOptionalParameter( "a3_Chi2OverDOFCutValue" , 
			      "Cut Value"  ,
			      _Chi2OverDOFCutValue,
			      float(10.0) ) ;

  registerOptionalParameter( "b1_D0Enable" , 
			      "Enable a cut on the value of each tracks d0 (no correction for ref point position)"  ,
			      _D0Enable,
			      bool(0) ) ;
  registerOptionalParameter( "b2_D0CutLowerThan" , 
			      "If true values lower than the cut value will be cut, if false values higher will be cut"  ,
			      _D0CutLowerThan,
			      bool(0) ) ;
  registerOptionalParameter( "b3_D0CutValue" , 
			      "Cut Value"  ,
			      _D0CutValue,
			      float(20.0) ) ;

  registerOptionalParameter( "c1_D0ErrEnable" , 
			      "Enable a cut on the value of each tracks d0 Error (sqrt(covariance(d0,d0))"  ,
			      _D0ErrEnable,
			      bool(0) ) ;
  registerOptionalParameter( "c2_D0ErrCutLowerThan" , 
			      "If true values lower than the cut value will be cut, if false values higher will be cut"  ,
			      _D0ErrCutLowerThan,
			      bool(0) ) ;
  registerOptionalParameter( "c3_D0ErrCutValue" , 
			      "Cut Value"  ,
			      _D0ErrCutValue,
			      float(250.0/1000.0) ) ;

  registerOptionalParameter( "d1_Z0Enable" , 
			      "Enable a cut on the value of each tracks Z0 (no correction for ref point position)"  ,
			      _Z0Enable,
			      bool(0) ) ;
  registerOptionalParameter( "d2_Z0CutLowerThan" , 
			      "If true values lower than the cut value will be cut, if false values higher will be cut"  ,
			      _Z0CutLowerThan,
			      bool(0) ) ;
  registerOptionalParameter( "d3_Z0CutValue" , 
			      "Cut Value"  ,
			      _Z0CutValue,
			      float(20.0) ) ;

  registerOptionalParameter( "e1_Z0ErrEnable" , 
			      "Enable a cut on the value of each tracks z0 Error (sqrt(covariance(z0,z0))"  ,
			      _Z0ErrEnable,
			      bool(0) ) ;
  registerOptionalParameter( "e2_Z0ErrCutLowerThan" , 
			      "If true values lower than the cut value will be cut, if false values higher will be cut"  ,
			      _Z0ErrCutLowerThan,
			      bool(0) ) ;
  registerOptionalParameter( "e3_Z0ErrCutValue" , 
			      "Cut Value"  ,
			      _Z0ErrCutValue,
			      float(250.0/1000.0) ) ;
			      
  registerOptionalParameter( "f1_PTEnable" , 
			      "Enable a cut on the value of each tracks PT (radial magnitude of ReconstructedParticle->momentum())"  ,
			      _PTEnable,
			      bool(0) ) ;
  registerOptionalParameter( "f2_PTCutLowerThan" , 
			      "If true values lower than the cut value will be cut, if false values higher will be cut"  ,
			      _PTCutLowerThan,
			      bool(1) ) ;
  registerOptionalParameter( "f3_PTCutValue" , 
			      "Cut Value"  ,
			      _PTCutValue,
			      float(0.1) ) ;

  registerOptionalParameter( "g1_DetectorHitsEnable" , 
			      "Enable a cut on the number seen hits in sub detectors - for more details see documentation"  ,
			      _DetectorHitsEnable,
			      bool(0) ) ;
  std::vector<std::string> SubDetectorNames;
  SubDetectorNames.push_back("VTX");
  SubDetectorNames.push_back("FTD");
  SubDetectorNames.push_back("SIT");
  SubDetectorNames.push_back("TPC");
  registerOptionalParameter( "g2_SubDetectorNames" , 
			      "Sub detector names in same order as result of Track->getSubdetectorHitNumbers()"  ,
			      _DetectorNames,
			      SubDetectorNames ) ;
  std::vector<std::string> BoundaryDetectorNames;
  BoundaryDetectorNames.push_back("TPC");
  BoundaryDetectorNames.push_back("FTD");
  registerOptionalParameter( "g3_DetectorHitsBoundaryDetectorNames" , 
			      "Sub detector names of detectors defining the boundary between region 1 and region 2"  ,
			      _DetectorHitsBoundaryDetectorNames,
			      BoundaryDetectorNames ) ;
  std::vector<int> BoundaryCuts;
  BoundaryCuts.push_back(20);
  BoundaryCuts.push_back(3);
  registerOptionalParameter( "g4_DetectorHitsBoundaryCuts" , 
			      "Corresponding to the order of DetectorHitsBoundaryDetectorNames the max number of hits for each detector for region 1, if any of the sub detectors has a higher number of hits then region 2 is used"  ,
			      _DetectorHitsBoundaryCuts,
			      BoundaryCuts ) ;
  
  std::vector<std::string> Region1DetectorNames;
  Region1DetectorNames.push_back("VTX");
  registerOptionalParameter( "g5_DetectorHitsRegion1DetectorNames" , 
			      "Sub detector names of detectors used for cutting in region 1"  ,
			      _DetectorHitsRegion1DetectorNames,
			      Region1DetectorNames ) ;
  std::vector<int> Region1Cuts;
  Region1Cuts.push_back(3);
  registerOptionalParameter( "g6_DetectorHitsRegion1Cuts" , 
			      "Corresponding to the order of DetectorHitsRegion1DetectorNames the minimum number of hits for region 1 for a track to pass"  ,
			      _DetectorHitsRegion1Cuts,
			      Region1Cuts ) ;

  std::vector<std::string> Region2DetectorNames;
  Region2DetectorNames.push_back("VTX");
  registerOptionalParameter( "g7_DetectorHitsRegion2DetectorNames" , 
			      "Sub detector names of detectors used for cutting in region 2"  ,
			      _DetectorHitsRegion2DetectorNames,
			      Region2DetectorNames );
  std::vector<int> Region2Cuts;
  Region2Cuts.push_back(0);
  registerOptionalParameter( "g8_DetectorHitsRegion2Cuts" , 
			      "Corresponding to the order of DetectorHitsRegion2DetectorNames the minimum number of hits for region 2 for a track to pass"  ,
			      _DetectorHitsRegion2Cuts,
			      Region2Cuts );

  // The Processor parameters that will cut on Monte Carlo PDG code of a particles parent.  A boolean to enable this cut and an integer
  // vector of the PDG codes to cut on. Also need the name of the LCRelation collection which refers back to the
  // Monte Carlo data.
  registerOptionalParameter( "h1_MCPIDEnable" , 
			      "Enable a cut on the PDG code of the parent of Monte Carlo particle associated with the track. Set the codes to cut on in CutPIDS and the LCRelation collection name to the MC data in MonteCarloLCRelationCollection"  ,
			      _MonteCarloPIDEnable,
			      bool(0) ) ;
  std::vector<int> PDGCodeCuts; //just set up a blank vector for the default input
  PDGCodeCuts.push_back(0);
  registerOptionalParameter( "h2_CutPIDS" , 
			      "A list of all the PDG codes of the parent Monte Carlo particle to cut"  ,
			      _MonteCarloPIDsToCut,
			      PDGCodeCuts );
  registerInputCollection( lcio::LCIO::LCRELATION, 
			      "h3_MonteCarloLCRelationCollection" , 
			      "Name of the LCRelation collection which links InputRCPCollection to the Monte Carlo data. Required only if MCPIDEnable is true"  ,
			      _MonteCarloRelationColName ,
			      std::string("Relations") ) ;
  registerOptionalParameter( "i1_BadParametersEnable" , 
			      "Enable a cut on tracks with NaN parameters or covariances"  ,
			      _BadParametersEnable,
			      bool(0) ) ;
  registerOptionalParameter( "j1_MCVertexCut" , 
			      "Enable a cut on tracks with MC Production Vertices in material"  ,
			      _MCVertexEnable,
			      bool(0) ) ;
  			      
}


void RPCutProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  if (_MCVertexEnable) {
    _VxdPar = &(Global::GEAR->getVXDParameters());
    const gear::GearParameters& BeamPipePar
      = Global::GEAR->getGearParameters("BeamPipe");
    _BeamPipeInnerRadius=BeamPipePar.getDoubleVal("BeamPipeRadius");
    _BeamPipeOuterRadius=_BeamPipeInnerRadius
      +BeamPipePar.getDoubleVal("BeamPipeThickness");
    _BeamPipeHalfZ=BeamPipePar.getDoubleVal("BeamPipeHalfZ");

    // tracks originating from anywhere closer than the following value
    // to a ladder or to the beam pipe will be removed if the user selects
    // hadronic interaction suppression. The unit is mm.
    // I chose a distance of 20 micron based on a histogram of the distribution
    // of track origin distances from ladder material, which shows a clear
    // excess between 0 and 20 micron, probably due to the GEANT4 step size.
    _CutDistance=0.02;

    // widen beam pipe size parameters accordingly
    _BeamPipeInnerRadius-=_CutDistance;
    _BeamPipeOuterRadius+=_CutDistance;
    _BeamPipeHalfZ+=_CutDistance;

#ifdef MCFAIL_DIAGNOSTICS
    _diaghist_vxmat_xy = new TH2F("vxmathist_xy","xy distribution of track origins supposedly in VXD material",1000,-60,60,1000,-60,60);
    _diaghist_bpmat_xy = new TH2F("bpmathist_xy","xy distribution of track origins supposedly in beam pipe material",1000,-60,60,1000,-60,60);
    _diaghist_nomat_xy = new TH2F("nomathist_xy","xy distribution of track origins supposedly outside of material",1000,-60,60,1000,-60,60);
    _diaghist_vxmat_rz = new TH2F("vxmathist_rz","rz distribution of track origins supposedly in VXD material",1000,-300,300,1000,0,60);
    _diaghist_bpmat_rz = new TH2F("bpmathist_rz","rz distribution of track origins supposedly in beam pipe material",1000,-300,300,1000,0,60);
    _diaghist_nomat_rz = new TH2F("nomathist_rz","rz distribution of track origins supposedly outside of material",1000,-300,300,1000,0,60);
    _diaghist_dist = new TH1F("disthist","distance of track origins from nearest ladder",100,0.0,2);
    _diaghist_dist_vxmat = new TH1F("disthist_vxmat","distance of track origins from nearest ladder, inside VXD material",100,0.0,2);
    _diaghist_dist_nomat = new TH1F("disthist_nomat","distance of track origins from nearest ladder, outside VXD material",100,0.0,0.2);
#endif
  }
}

void RPCutProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void RPCutProcessor::processEvent( LCEvent * evt ) { 

  // this gets called for every event 
  // usually the working horse ...

	LCCollection* InCol = evt->getCollection( _InRCPColName );

	//Get the collection of associated Monte Carlo particles if the cut on MC PDG code is enabled.
	//Use a pointer for the LCRelationNavigator, because I only want it created under certain conditions but it needs to be refered to later.
	UTIL::LCRelationNavigator* pMCRelationNavigator=0;
	if(_MonteCarloPIDEnable)
	{
		lcio::LCCollection* RelCol=0;
		try
		{
			RelCol = evt->getCollection( _MonteCarloRelationColName );
		}
		catch( lcio::Exception exception )
		{
			// Don't want to quit, so don't do much here. Just ignore MC PID cuts and apply the other cuts.
			// Print a warning though.
			std::cerr << "Unable to get the Monte Carlo data for event " << _nEvt << ". The requested cuts on the true PID will not be applied." << std::endl;
			RelCol=0;
		}
		if(RelCol!=0) // the above "try" was successful
		{
			pMCRelationNavigator=new UTIL::LCRelationNavigator(RelCol);
			vertex_lcfi::MemoryManager<UTIL::LCRelationNavigator>::Event()->registerObject(pMCRelationNavigator);
		}
		else // the above "try" was unsuccessful
		{
			pMCRelationNavigator=0;
		}
	}
	
	LCCollection* OutRPCollection;
	if (_WriteNewCollection)
	{
		std::vector<std::string>::const_iterator it = find(evt->getCollectionNames()->begin(),evt->getCollectionNames()->end(),_OutRCPColName);
		if (it == evt->getCollectionNames()->end())
		{
			//Not found do add - TODO do these need memory managment?
			LCCollection* MyCollection = new LCCollectionVec("ReconstructedParticle");
			evt->addCollection(MyCollection,_OutRCPColName);
		}
		
	}
	OutRPCollection = evt->getCollection(_OutRCPColName);
	if (_WriteNewCollection && !_SubParticleLists)
	{
		dynamic_cast<LCCollectionVec*>(OutRPCollection)->setSubset(true);
	}
	
	
	//TODO Only do this if cut enables and check for data existance
	//Retrive the list of detector names and create a map so that we know at what index to find them 
	std::map<std::string,int> SubdetectorIndex;
	if (_DetectorHitsEnable)
	{
		if (_DetectorNames.empty())
		{
			std::cerr << "Subdetector Names not set" << std::endl;
		}
		int i=0;
		for (StringVec::iterator iDetector = _DetectorNames.begin(); iDetector != _DetectorNames.end();++iDetector)
		{
			SubdetectorIndex[*iDetector] = i;
			++i;
		}
	}
	//RP Loop
	int nRCP = InCol->getNumberOfElements()  ;
	for(int i=0; i< nRCP ; i++)
	{
		ReconstructedParticle*  InputRP = dynamic_cast<ReconstructedParticle*>( InCol->getElementAt( i ) );
		
		if (_SubParticleLists)
		{
			//We are cutting particles listed in the InputRP
			ReconstructedParticle* OutputRP;
			if (_WriteNewCollection)
			{
				//Get a copy of the input to cut RPs from
				OutputRP = new ReconstructedParticleImpl(*dynamic_cast<ReconstructedParticleImpl*>(InputRP));
				//To do a proper copy we need to copy the PID objects seperatly and remove the read only
				//TODO Get around nasty cast
				((ReconstructedParticleLCFI*)OutputRP)->makeWritable();
				((ReconstructedParticleLCFI*)OutputRP)->wipePIDs();
				((ReconstructedParticleLCFI*)OutputRP)->copyPIDsFrom(InputRP);
			}
			else
			{
				OutputRP = InputRP;
			}
				
			ReconstructedParticleVec RPTracks = OutputRP->getParticles();
			for (ReconstructedParticleVec::iterator iRPTrack = RPTracks.begin();iRPTrack != RPTracks.end();++iRPTrack)
			{
				if ((*iRPTrack)->getTracks().size() == 0 
				    ||((_D0Enable && _D0Fail(*iRPTrack)) ||
				       (_D0ErrEnable && _D0ErrFail(*iRPTrack)) ||
				       (_Z0Enable && _Z0Fail(*iRPTrack)) ||
				       (_Z0ErrEnable && _Z0ErrFail(*iRPTrack)) ||
				       (_PTEnable && _PTFail(*iRPTrack)) ||
				       (_Chi2OverDOFEnable && _Chi2OverDOFFail(*iRPTrack)) ||
				       (_DetectorHitsEnable && _DetectorHitsFail(*iRPTrack,SubdetectorIndex)) ||
				       (_MonteCarloPIDEnable && _MCPIDFail(*iRPTrack,pMCRelationNavigator)) ||
				       (_BadParametersEnable && _BadParametersFail(*iRPTrack)) ||
				       (_MCVertexEnable &&_MCVertexFail((*iRPTrack),pMCRelationNavigator)) ))
				{
					//TODO Remove nasty cast
					((ReconstructedParticleLCFI*)OutputRP)->removeParticle(*iRPTrack);
				}
			}
			if (_WriteNewCollection)
			{
				OutRPCollection->addElement(OutputRP);
			}
			else
			{
				//NO OP//
			}	
		}
		else //NOT SubParticleLists
		{
			//We are deciding weather to cut InputRP itself
			if ((InputRP)->getTracks().size() == 0
			    ||((_D0Enable && _D0Fail(InputRP)) ||
			       (_D0ErrEnable && _D0ErrFail(InputRP)) ||
			       (_Z0Enable && _Z0Fail(InputRP)) ||
			       (_Z0ErrEnable && _Z0ErrFail(InputRP)) ||
			       (_PTEnable && _PTFail(InputRP)) ||
			       (_Chi2OverDOFEnable && _Chi2OverDOFFail(InputRP)) ||
			       (_DetectorHitsEnable && _DetectorHitsFail(InputRP,SubdetectorIndex)) ||
			       (_MonteCarloPIDEnable && _MCPIDFail(InputRP,pMCRelationNavigator)) ||
			       (_BadParametersEnable && _BadParametersFail(InputRP)) ||
			       (_MCVertexEnable &&_MCVertexFail(InputRP,pMCRelationNavigator)) ))
			{
				//Cut InputRP either by removing it from this collection or not adding it to a new one
				if (_WriteNewCollection)
				{
				}
				else
				{
					InCol->removeElementAt(i);						
				}
			}
			else
			{
				//Keep InputRP either by not removing it from this collection or adding it to a new one
				if (_WriteNewCollection)
				{
					OutRPCollection->addElement(InputRP);
				}
				else
				{					
				}
			}
		}				     
	}
	//Clear all objects created
	vertex_lcfi::MetaMemoryManager::Event()->delAllObjects();
	_nEvt ++ ;
}



void RPCutProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void RPCutProcessor::end(){ 

#ifdef MCFAIL_DIAGNOSTICS
  if (_MCVertexEnable) {
        TFile dumpfile((name()+std::string(".root")).c_str(),"RECREATE");
	_diaghist_vxmat_xy->Write();
	_diaghist_bpmat_xy->Write();
	_diaghist_nomat_xy->Write();
	_diaghist_vxmat_rz->Write();
	_diaghist_bpmat_rz->Write();
	_diaghist_nomat_rz->Write();
	_diaghist_dist->Write();
	_diaghist_dist_vxmat->Write();
	_diaghist_dist_nomat->Write();
	dumpfile.Close();
  }
#endif

	std::cout << "RPCutProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

bool  RPCutProcessor::_Chi2OverDOFFail(ReconstructedParticle* RPTrack)
{
	lcio::Track* Track = RPTrack->getTracks()[0];
	if (_Chi2OverDOFCutLowerThan)
		return (Track->getChi2()/(float)Track->getNdf() < _Chi2OverDOFCutValue);
	else
		return (Track->getChi2()/(float)Track->getNdf() > _Chi2OverDOFCutValue);
}
bool  RPCutProcessor::_D0Fail(ReconstructedParticle* RPTrack)
{
	lcio::Track* Track = RPTrack->getTracks()[0];
	if (_D0CutLowerThan)
		return (fabs(Track->getD0()) < _D0CutValue);
	else
		return (fabs(Track->getD0()) > _D0CutValue);
}
bool  RPCutProcessor::_D0ErrFail(ReconstructedParticle* RPTrack)
{
	lcio::Track* Track = RPTrack->getTracks()[0];
	if (_D0ErrCutLowerThan)
		return (Track->getCovMatrix()[0] < _D0ErrCutValue);
	else
		return (Track->getCovMatrix()[0] > _D0ErrCutValue);
}
bool  RPCutProcessor::_Z0Fail(ReconstructedParticle* RPTrack)
{
	lcio::Track* Track = RPTrack->getTracks()[0];
	if (_Z0CutLowerThan)
		return (fabs(Track->getZ0()) < _Z0CutValue);
	else
		return (fabs(Track->getZ0()) > _Z0CutValue);
}
bool  RPCutProcessor::_Z0ErrFail(ReconstructedParticle* RPTrack)
{
	lcio::Track* Track = RPTrack->getTracks()[0];
	if (_Z0ErrCutLowerThan)
		return (Track->getCovMatrix()[9] < _Z0ErrCutValue);
	else
		return (Track->getCovMatrix()[9] > _Z0ErrCutValue);
}
bool  RPCutProcessor::_PTFail(ReconstructedParticle* RPTrack)
{
	const double* p = RPTrack->getMomentum();
	double pt = sqrt(p[0]*p[0] + p[1]*p[1]); 
	if (_PTCutLowerThan)
		return (pt < _PTCutValue);
	else
		return (pt > _PTCutValue);
}
bool  RPCutProcessor::_DetectorHitsFail(ReconstructedParticle* RPTrack, std::map<std::string,int> SubdetectorIndex)
{
	//First find out if this track is in region 1 or 2
	bool Region2 = 0;
	std::vector<int>::const_iterator iCut = _DetectorHitsBoundaryCuts.begin();
	for (StringVec::iterator iDetector = _DetectorHitsBoundaryDetectorNames.begin(); iDetector != _DetectorHitsBoundaryDetectorNames.end();++iDetector)
	{
		//TODO Check for exisance of data
		if (RPTrack->getTracks()[0]->getSubdetectorHitNumbers()[SubdetectorIndex[*iDetector]] >= (*iCut))
		{
			Region2 = 1;
			break;
		}
		++iCut;
	}
	//Cut accordingly  1 = fail
	if (Region2)
	{
		iCut = _DetectorHitsRegion2Cuts.begin();
		for (StringVec::iterator iDetector = _DetectorHitsRegion2DetectorNames.begin(); iDetector != _DetectorHitsRegion2DetectorNames.end();++iDetector)
		{
			if (RPTrack->getTracks()[0]->getSubdetectorHitNumbers()[SubdetectorIndex[*iDetector]] < (*iCut))
			{
				return 1;
			}
			++iCut;
		}
	}
	else
	{
		iCut = _DetectorHitsRegion1Cuts.begin();
		for (StringVec::iterator iDetector = _DetectorHitsRegion1DetectorNames.begin(); iDetector != _DetectorHitsRegion1DetectorNames.end();++iDetector)
		{
			if (RPTrack->getTracks()[0]->getSubdetectorHitNumbers()[SubdetectorIndex[*iDetector]] < (*iCut))
			{
				return 1;
			}
			++iCut;
		}
	}
	return 0;
}

bool RPCutProcessor::_MCPIDFail( lcio::ReconstructedParticle* RPTrack, UTIL::LCRelationNavigator* pMCRelationNavigator )
{
	if( pMCRelationNavigator==0 )
	{
		// Unable to get the MC relation data earlier, so just return false because MC PID cuts are
		// being ignored.  A warning to cerr will have already been printed so no need to do so here.
		return false;
	}

	lcio::Track* Track = RPTrack->getTracks()[0];
	std::vector<lcio::LCObject*> RelatedMCParticles = pMCRelationNavigator->getRelatedToObjects(Track);

	//Not sure what to do if it has more than one related MC particle.  For now just print a warning and ignore.
	if( RelatedMCParticles.size()!=1 )
	{
		std::cerr << RelatedMCParticles.size() << " MCParticles related to this track! No cuts performed on MC PDG code." << std::endl;
		return false;
	}
	else
	{
		std::vector<lcio::MCParticle*> Parents = dynamic_cast<lcio::MCParticle*>(RelatedMCParticles[0])->getParents();
		
		if (Parents.empty())
			return false;
		else
		{
			int truePDGCode=Parents[0]->getPDG();
			//search for this code in the integer vector of codes to cut on
			std::vector<int>::iterator iSearchResult=std::find( _MonteCarloPIDsToCut.begin(), _MonteCarloPIDsToCut.end(), truePDGCode );
			if( iSearchResult==_MonteCarloPIDsToCut.end() ) 
				return false; // the PDG code was not found in the vector
			else 
			{
				return true; // the PDG code was found
			}
		}
	}
}

bool RPCutProcessor::_BadParametersFail(lcio::ReconstructedParticle* RPTrack)
{
	lcio::Track* Track = RPTrack->getTracks()[0];
	bool fail=false;
	for (short i=0;!fail && i<15;++i)
	{
		fail = std::isnan(Track->getCovMatrix()[i]);
	}
	return (fail
	     || std::isnan(Track->getD0())
	     || std::isnan(Track->getPhi())
	     || std::isnan(Track->getZ0())
	     || std::isnan(Track->getOmega())
	     || std::isnan(Track->getTanLambda()) );
}

bool  RPCutProcessor::_MCVertexFail(lcio::ReconstructedParticle* RPTrack, UTIL::LCRelationNavigator* pMCRelationNavigator )
{
	LCObjectVec objectVec = pMCRelationNavigator->getRelatedToObjects(RPTrack->getTracks()[0]);
	
	if (objectVec.size() > 0)
	{
		MCParticle* mcp = dynamic_cast<MCParticle*> (objectVec[0]);
		const double * Vert = mcp->getVertex();
		gear::Vector3D VertVec(Vert[0],Vert[1],Vert[2]);

		double radius = sqrt((Vert[0] * Vert[0])+(Vert[1] * Vert[1]));
		if (radius <_BeamPipeInnerRadius) {
		  // this track originates from inside the beam pipe.
		  // we assume that there is no material interaction there.
#ifdef MCFAIL_DIAGNOSTICS
		  _diaghist_nomat_xy->Fill(Vert[0],Vert[1]);
		  _diaghist_nomat_rz->Fill(Vert[2],radius);
#endif
		  return false;
		}else if (radius >=_BeamPipeInnerRadius
		    && radius <=_BeamPipeOuterRadius
		    && fabs(Vert[2])<=_BeamPipeHalfZ) {
		  // this track originates from within the beam pipe wall.
		  // thus it is likely a hadronic interaction product
#ifdef MCFAIL_DIAGNOSTICS
		  _diaghist_bpmat_xy->Fill(Vert[0],Vert[1]);
		  _diaghist_bpmat_rz->Fill(Vert[2],radius);
#endif
		  return true;
		} else {
		  // next check: are we in VXD ladder material?
		  double dist=_VxdPar->distanceToNearestLadder(VertVec).r();
#ifdef MCFAIL_DIAGNOSTICS
		  _diaghist_dist->Fill(dist);
#endif
		  if (dist<=_CutDistance) {
#ifdef MCFAIL_DIAGNOSTICS
       	            _diaghist_dist_vxmat->Fill(dist);
		    _diaghist_vxmat_xy->Fill(Vert[0],Vert[1]);
		    _diaghist_vxmat_rz->Fill(Vert[2],radius);
#endif
		    return true;
		  } else {
#ifdef MCFAIL_DIAGNOSTICS
       	            _diaghist_dist_nomat->Fill(dist);
		    _diaghist_nomat_xy->Fill(Vert[0],Vert[1]);
		    _diaghist_nomat_rz->Fill(Vert[2],radius);
#endif
		    return false;
		  }
		}
	}
	else
		std::cerr << "Warning RPCutProcessor::637 No MCParticle information" << std::endl;
	return false;	
}

