// First of all, make sure that AIDA was enabled by setting the environment
// variable MARLIN_USE_AIDA when Marlin was compiled. The makefile will then
// have done the setup and defined this macro.

#ifndef MARLIN_USE_AIDA

#warning "--------------------------------------------------------------------------------"
#warning "- LCFIAIDAPlotProcessor requires MARLIN_USE_AIDA to be defined. Did you enable -"
#warning "- AIDA when compiling Marlin? LCFIAIDAPlotProcessor will not be compiled.      -"
#warning "--------------------------------------------------------------------------------"

// Can't do anything else.
#else


#include "LCFIAIDAPlotProcessor.h" 
 
// standard library includes 
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <math.h>

// LCIO includes... 
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/Vertex.h"

// Marlin includes
#include <marlin/Exceptions.h>

// AIDA includes...
#include <marlin/AIDAProcessor.h>
#include <AIDA/IDataPointSet.h>
#include <AIDA/IDataPointSetFactory.h>
#include <AIDA/IDataPoint.h>
#include <AIDA/IMeasurement.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IAxis.h>
#include <AIDA/ITree.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ICloud2D.h>

#include "TypesafeCollection.h"

// There needs to be at least one instantiation for the base constructor to register the processor with 
// the Marlin processor manager. This is it. 
LCFIAIDAPlotProcessor aLCFIAIDAPlotProcessor; 

LCFIAIDAPlotProcessor::LCFIAIDAPlotProcessor() : marlin::Processor( "LCFIAIDAPlotProcessor" ) 
{ 

  _description="Creates an AIDA plot of the LCFIVertex tagging efficiency-purity values and various other things.  Make sure that MarlinAIDAProcessor is run before this.";
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" , 
			   "Name of the collection of ReconstructedParticles that is the jet"  ,
			   _JetCollectionName ,
			   std::string("FTSelectedJets") );
  
  std::vector<std::string> FlavourTagCollectionNamesDefault;
  FlavourTagCollectionNamesDefault.push_back("FlavourTag");
  registerProcessorParameter("FlavourTagCollections" , 
			     "Names of the LCFloatVec Collections that contain the flavour tags (one purity efficiency plot per tag) (in same order as jet collection)"  ,
			     _FlavourTagCollectionNames,
			     FlavourTagCollectionNamesDefault) ;
  
  registerInputCollection( LCIO::LCINTVEC,
			   "TrueJetFlavourCollection" , 
			   "Name of the LCIntVec collection containing the true flavour of the jets (same order as jets)"  ,
			   _TrueJetFlavourColName ,
			   std::string("TrueJetFlavour") ) ;

  registerInputCollection( LCIO::LCFLOATVEC,
			   "TrueJetHadronChargeCollection",
			   "Name of the LCFloatVec collection containing the true hadron charge of the jets (same order as jets)"  ,
			   _TrueJetHadronChargeColName ,
			   std::string("TrueJetHadronCharge") ) ;  

  registerInputCollection( LCIO::LCINTVEC,
			   "TrueJetPDGCodeCollection" , 
			   "Name of the LCIntVec collection containing the true PDG code of the jets (same order as jets)"  ,
			   _TrueJetPDGCodeColName,
			   std::string("TrueJetPDGCode") ) ;

  registerInputCollection( LCIO::LCFLOATVEC,
			   "TrueJetPartonChargeCollection",
			   "Name of the LCFloatVec collection containing the true parton charge of the jets (same order as jets)"  ,
			   _TrueJetPartonChargeColName ,
			   std::string("TrueJetPartonCharge") ) ;    
  
  registerInputCollection( lcio::LCIO::VERTEX,
			   "VertexCollectionName",
			   "Name of the collection that holds the Vertices",
			   _VertexColName,
			   std::string("ZVRESVertices") ) ;

  registerInputCollection(LCIO::LCFLOATVEC,
			  "CVertexChargeCollection",
			  "Name of collection containing the vertex charge of the jets, assuming they are C-jets",
			  _CVertexChargeCollection,
			  std::string("CCharge") );
  
  registerInputCollection( LCIO::LCFLOATVEC,
			   "BVertexChargeCollection",
			   "Name of collection containing the vertex charge of the jets, assuming they are B-jets",
			   _BVertexChargeCollection,
			   std::string("BCharge") ) ;
  

  FlavourTagCollectionNamesDefault.clear();
  FlavourTagCollectionNamesDefault.push_back("FlavourTagInputs");
  registerProcessorParameter("TagInputsCollections" , 
			     "Names of the LCFloatVec Collections that contain the flavour tag inputs (in same order as jet collection)"  ,
			     _FlavourTagInputsCollectionNames,
			     FlavourTagCollectionNamesDefault) ;

  registerOptionalParameter( "CosThetaJetMax",
			     "Cut determining the maximum cos(theta) of the jet.  Default: |cos(theta)|<0.9"  ,
			     _CosThetaJetMax,
			     double(0.9)) ;
   
  registerOptionalParameter( "CosThetaJetMin",
			     "Cut determining the minimum cos(theta) of the jet.  Default: no lower cut."  ,
			     _CosThetaJetMin,
			     double(0.0)) ;

  registerOptionalParameter("PJetMax",
			    "Cut determining the maximum momentum of the jet.  Default: 10000 GeV/c"  ,
			    _PJetMax,
			    double(10000.)) ;
   
  registerOptionalParameter( "PJetMin",
			     "Cut determining the minimum momentum of the jet.  Default: no lower cut."  ,
			     _PJetMin,
			     double(0.0)) ;

  registerOptionalParameter( "PrintNeuralNetOutput",
			     "Set true if you want a print-out of the NN values (output) for the various flavour tags",
			     _PrintNeuralNetOutput,
			     bool(false));

  registerOptionalParameter( "NeuralNetOutputFile" , 
			     "Output filename for the NN values (output).  Only used if PrintNeuralNetOutput parameter is true.  If left blank, output will be directed to standard out.",
			     _NeuralNetOutputFile,
			     std::string("") ) ;

  registerOptionalParameter( "MakeTuple",
			     "Set true to make a tuple of the flavour tag input variables.  Default is true.",
			     _MakeTuple,
			     bool(true));

  registerOptionalParameter( "CTagNNCut",
			     "Cut determining the Neural Net cut used to select C-Jets",
			     _CTagNNCut,
			     double(0.7));

  registerOptionalParameter( "BTagNNCut",
			     "Cut determining the Neural Net cut used to select B-Jets",
			     _BTagNNCut,
			     double(0.7));

  registerOptionalParameter( "UseFlavourTagCollectionForVertexCharge",
			     "Integer parameter determing which FlavourTag Collection to use the determine C-Jets and B-Jets in Vertex Charge Plots",
			     _iVertexChargeTagCollection,
			     int(0));

} 

LCFIAIDAPlotProcessor::~LCFIAIDAPlotProcessor() 
{ 
} 

void LCFIAIDAPlotProcessor::init()
{

  
  if (_iVertexChargeTagCollection >=  int(_FlavourTagCollectionNames.size()) || _iVertexChargeTagCollection < 0) {
    std::cerr << " In " << __FILE__ << "(" << __LINE__ << "): Invalid parameter for UseFlavourTagCollectionForVertexCharge.  Setting to 0." << std::endl;
    _myVertexChargeTagCollection = 0;
  } else {
    _myVertexChargeTagCollection = uint(_iVertexChargeTagCollection);
  }

  _ZoomedVarNames.push_back("D0Significance1"); 
  _ZoomedVarNames.push_back("D0Significance2");
  _ZoomedVarNames.push_back("Z0Significance1");
  _ZoomedVarNames.push_back("Z0Significance2");

  _VertexCatNames.resize(N_VERTEX_CATEGORIES+1);
  _VertexCatNames[0]="AnyNumberOfVertices";
  _VertexCatNames[1]="OneVertex";
  _VertexCatNames[2]="TwoVertices";
  _VertexCatNames[3]="ThreeOrMoreVertices";

  
  _NumVertexCatDir.resize(N_VERTEX_CATEGORIES+1);
  _NumVertexCatDir[1]="OneVertex";
  _NumVertexCatDir[2]="TwoVertices";
  _NumVertexCatDir[3]="ThreeOrMoreVertices";
  _NumVertexCatDir[0]="AnyNumberOfVertices";


  _numberOfPoints=100;

  
  _pBJetBTag.resize( _FlavourTagCollectionNames.size() );
  _pBJetCTag.resize( _FlavourTagCollectionNames.size() );
  _pBJetBCTag.resize( _FlavourTagCollectionNames.size() );
  
  _pCJetCTag.resize( _FlavourTagCollectionNames.size() );
  _pCJetBTag.resize( _FlavourTagCollectionNames.size() );
  _pCJetBCTag.resize( _FlavourTagCollectionNames.size() );
  
  _pLightJetBTag.resize( _FlavourTagCollectionNames.size() ); 
  _pLightJetCTag.resize( _FlavourTagCollectionNames.size() ); 
  _pLightJetBCTag.resize( _FlavourTagCollectionNames.size() ); 
  
  _pBTagBackgroundValues.resize( _FlavourTagCollectionNames.size() );
  _pCTagBackgroundValues.resize( _FlavourTagCollectionNames.size() );
  _pBCTagBackgroundValues.resize( _FlavourTagCollectionNames.size() );

 
  
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    { 
      for (unsigned int iVertexCat=0;  iVertexCat <  N_VERTEX_CATEGORIES+1; ++iVertexCat ){
	_pLightJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
	_pLightJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	_pLightJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	
	_pBJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;	 
	_pBJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	_pBJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	
	_pCJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;	 
	_pCJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;  
	_pCJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;  
	
	_pBTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
	_pCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
	_pBCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
      }
    }

 

  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  AIDA::ITupleFactory* pTupleFactory=marlin::AIDAProcessor::tupleFactory( this );

  if(  pHistogramFactory!=0 )
    {
      bool ableToMakeAllHistograms=true;
      
      if (!(pTree->cd( "/" + name() + "/"))) {
	pTree->mkdir( "/" + name() + "/" );
	pTree->cd( "/" + name() + "/");
      }
      
 
      //some plots of vertex charge
      if (!pTree->cd("/" + name() + "/VertexChargePlots/")) {
	pTree->mkdirs("/" + name() + "/VertexChargePlots/");
	pTree->cd( "/" + name() + "/VertexChargePlots/");
      }
      
      _pBJetCharge2D = pHistogramFactory->createHistogram2D( "B Jets: Reconstructed Vertex Charge vs True Jet Charge",7,-3.5,+3.5,7,-3.5,+3.5);
      _pCJetCharge2D = pHistogramFactory->createHistogram2D( "C Jets: Reconstructed Vertex Charge vs True Jet Charge",7,-3.5,+3.5,7,-3.5,+3.5);
      
      _pBJetVertexCharge = pHistogramFactory->createHistogram1D( "B Jets: Reconstructed Vertex Charge",9,-4.5,+4.5);
      _pCJetVertexCharge = pHistogramFactory->createHistogram1D( "C Jets: Reconstructed Vertex Charge",9,-4.5,+4.5);
      
      _pCJetLeakageRate = pHistogramFactory->createHistogram1D("C Jets: Charged Leakage Rate  (DON'T TRUST ERRORS)", N_JETANGLE_BINS,0.,1.);
      _pBJetLeakageRate = pHistogramFactory->createHistogram1D("B Jets: Charged Leakage Rate  (DON'T TRUST ERRORS)", N_JETANGLE_BINS,0.,1.);
      
      
      
      for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
	{
	  for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ){
	    
	    std::string nvname = _VertexCatNames[iVertexCat];
	    
	    if (!pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection] + "/" + _NumVertexCatDir[iVertexCat])) {
	      pTree->mkdirs( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection] + "/" + _NumVertexCatDir[iVertexCat]);
	      pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection] + "/" + _NumVertexCatDir[iVertexCat]);
	    }
	    
	    _pLightJetBTag[iTagCollection][nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by B-tag NN value. ("+ nvname +")",_numberOfPoints , 0, 1.0 );
	    _pLightJetCTag[iTagCollection][nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pLightJetBCTag[iTagCollection][nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pBJetBTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of B jets by B-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
	    _pBJetCTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of B jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
	    _pBJetBCTag[iTagCollection][nvname]    = pHistogramFactory->createHistogram1D( "Numbers of B jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pCJetBTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of C jets by B-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
	    _pCJetCTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of C jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pCJetBCTag[iTagCollection][nvname]    = pHistogramFactory->createHistogram1D( "Numbers of C jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );    
	
	  }
	  
	  if (_MakeTuple) {
	    pTree->cd( "/" + name());
	    
	    //make the ntuple
	    //this breaks the paradigm of reading these in from the flavour tag collections themselves
	    std::string columnNames="int TrueJetFlavour=-1,  int NumberOfVertices=-1, int NumberOfTracksInVertices=-1, float D0Significance1 = -999., float D0Significance2 = -999., float DecayLength = -999., float DecayLength_SeedToIP= -999., float DecayLengthSignificance= -999., float JointProbRPhi= -999., float JointProbZ= -999., float Momentum1= -999.,float Momentum2= -999., float PTCorrectedMass= -999., float RawMomentum= -999., float SecondaryVertexProbability= -999., float Z0Significance1= -999., float Z0Significance2= -999., int BQVtx=-10, int CQVtx=-10";
	    
	    if( !pTree->cd(  "/"  + name() + "/" + _FlavourTagInputsCollectionNames[iTagCollection]  + "/TupleDir/"))
	      {	 
		pTree->mkdirs( "/"  + name() + "/" + _FlavourTagInputsCollectionNames[iTagCollection]  + "/TupleDir/" ) ;
		pTree->cd(  "/"  + name() + "/" + _FlavourTagInputsCollectionNames[iTagCollection]  + "/TupleDir/");
	      }
	    
	    _pMyTuple=pTupleFactory->create( "FlavourTagInputsTuple","FlavourTagInputsTuple", columnNames);
	  }
	  
	}
      
      if (!pTree->cd( "/"  + name() + "/VertexPlots/")) {
 	pTree->mkdirs( "/"  + name() + "/VertexPlots/");
	pTree->cd( "/"  + name() + "/VertexPlots/");
      }
      
      _pVertexDistanceFromIP = pHistogramFactory->createHistogram1D( "Reconstructed Vertex distance from IP",100, 0., 10.);
      _pVertexPositionX = pHistogramFactory->createHistogram1D( "Non-primary vertex: x-position", 100, -10., 10.) ;
      _pVertexPositionY = pHistogramFactory->createHistogram1D( "Non-primary vertex: y-position", 100, -10., 10.);
      _pVertexPositionZ = pHistogramFactory->createHistogram1D( "Non-primary vertex: z-position", 100, -10., 10.);
      
      _pPrimaryVertexPullX = pHistogramFactory->createHistogram1D( "Non-primary vertex: x-pull", 100, -10., 10.);
      _pPrimaryVertexPullY = pHistogramFactory->createHistogram1D( "Non-primary vertex: y-pull", 100, -10., 10.);
      _pPrimaryVertexPullZ = pHistogramFactory->createHistogram1D( "Non-primary vertex: z-pull", 100, -10., 10.);
      _pPrimaryVertexPositionX = pHistogramFactory->createHistogram1D( "Primary vertex: x-position", 100, -10., 10.);
      _pPrimaryVertexPositionY = pHistogramFactory->createHistogram1D( "Primary vertex: y-position", 100, -10., 10.);
      _pPrimaryVertexPositionZ = pHistogramFactory->createHistogram1D( "Primary vertex: z-position", 100, -10., 10.);
    
      
      pTree->cd(  "/"  + name() + "/");
      if (!pTree->cd( "ZVRESInputPlots" )) {
	pTree->mkdir( "ZVRESInputPlots/" ) ;
	pTree->cd(  "ZVRESInputPlots/" ) ;
      }
      



      if( !ableToMakeAllHistograms )
	{
	  std::cerr << "### " << __FILE__ << "(" << __LINE__ << "): Unable to create some or all of the histograms for the flavour tag values!" << std::endl;
	 
	}
    }
  else
    {
      std::cerr  << "### " << __FILE__ << "(" << __LINE__ << "): Unable to get the histogram factory! No histograms will be made."<< std::endl;
    }
  
  _lastRunHeaderProcessed=-1;
  _suppressOutputForRun=-1;
  
  _inputsHistogramsBJets.resize( _FlavourTagInputsCollectionNames.size() );
  _inputsHistogramsCJets.resize( _FlavourTagInputsCollectionNames.size() );
  _inputsHistogramsUDSJets.resize( _FlavourTagInputsCollectionNames.size() );

  _zoomedInputsHistogramsBJets.resize( _FlavourTagInputsCollectionNames.size() );
  _zoomedInputsHistogramsCJets.resize( _FlavourTagInputsCollectionNames.size() );
  _zoomedInputsHistogramsUDSJets.resize( _FlavourTagInputsCollectionNames.size() );
  

  InternalVectorInitialisation();
  
}

void LCFIAIDAPlotProcessor::processRunHeader( LCRunHeader* pRun ) 
{

	// Marlin doesn't necessarily process the run header, e.g. if you use the "SkipNEvents"
	// parameter in the steering file. The flavour tag variable/tag value names are held in
	// the run header though, so this processor has to have access to it. Set this variable
	// so that "processEvent" can tell if "processRunHeader" has been called for the run
	// it's in.
	_lastRunHeaderProcessed=pRun->getRunNumber();


	//
	// Perform a check to see if the variable names we need are here
	//
	for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag) // Loop over the different tag collection names given in the steering
	{
	  std::vector<std::string> VarNames;
	  (pRun->parameters()).getStringVals(_FlavourTagCollectionNames[iTag],VarNames);
	  
	  
	  //Fill a map so that we can get the array index from just the string
	  std::set<std::string> AvailableNames;
	  std::map<std::string,unsigned int> IndexOf;
	  
	  for (size_t i = 0;i < VarNames.size();++i)
	    {
	      AvailableNames.insert(VarNames[i]);
	      IndexOf[VarNames[i]] = i;
	    }
	  
	  //Add the index to the list
	  _IndexOfForEachTag.push_back(IndexOf);
	  
	  //Check the required information is in the LCFloatVec
	  std::set<std::string> RequiredNames;
	  RequiredNames.insert("BTag");
	  RequiredNames.insert("CTag");
	  RequiredNames.insert("BCTag");
	  
	  if (!std::includes(AvailableNames.begin(),AvailableNames.end(),RequiredNames.begin(),RequiredNames.end()))
	    {
	      std::cerr << __FILE__ << "(" << __LINE__ << "): The collection \"" << _FlavourTagCollectionNames[iTag]
			<< "\" (if it exists) does not contain the tag values required by " << type() << "." << std::endl;
	      std::cerr <<   __FILE__ << "(" << __LINE__ << "): The collection \"" << _FlavourTagCollectionNames[iTag]
			<< "\" (if it exists) does not contain the tag values required by " << type() << "." << std::endl;
	    }
	}
	
	
	
	AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
	AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
		
	for (unsigned int iInputCollection=0; iInputCollection < _FlavourTagInputsCollectionNames.size(); ++iInputCollection)
	  {
	    	  

	    std::vector<std::string> VarNames;
	    (pRun->parameters()).getStringVals(_FlavourTagInputsCollectionNames[iInputCollection],VarNames);
	    
	    //Fill the map relating names and indexes
	    std::map<std::string,unsigned int> IndexOf;
	    for (size_t i = 0;i < VarNames.size();++i)
	      {
		
		IndexOf[VarNames[i]] = i;
		
		// If there is no histogram for this name then create one
		if( _inputsHistogramsBJets[iInputCollection][VarNames[i]]==0 )
		  {
		    if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] ) )
		      {
			pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] ) ; 
			pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] ) ;
		      }
		    
		    int numberOfPoints=_numberOfPoints/4;
		    double lowerBin=-1;
		    double higerBin=1;
		   
		    //binning variables: if the name is not listed here it will use the default above
		    if( VarNames[i]=="BQVtx" || VarNames[i]=="CQVtx" )
		      {
			numberOfPoints=9;
			lowerBin=-4.5;
			higerBin=4.5;
		      }
		    else if( VarNames[i]=="NumVertices" )
		      {
			numberOfPoints=5;
			lowerBin=0.5;
			higerBin=5.5;
		      }
		    else if( VarNames[i]=="NumTracksInVertices" )
		      {
			numberOfPoints=16;
			lowerBin=-0.5;
			higerBin=15.5;
		      }
		    else if ( VarNames[i]=="D0Significance1" || VarNames[i]=="D0Significance2" )
		      {
			numberOfPoints=120;
			lowerBin=-20.;
			higerBin=100.;
		      }
		    else if ( VarNames[i]=="Z0Significance1" || VarNames[i]=="Z0Significance2")
		      {
			numberOfPoints=100;
			lowerBin=-50.;
			higerBin=50.;
		      }
		    else if (VarNames[i]=="DecayLengthSignificance") 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=100.;
		      }
		    else if (VarNames[i]=="DecayLength" || VarNames[i]=="DecayLength(SeedToIP)" ) 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=10.;
		      }
		    else if (VarNames[i]=="JointProbRPhi" || VarNames[i]=="JointProbZ"|| VarNames[i]=="SecondaryVertexProbability") 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=1.0;
		      }
		    else if (VarNames[i]=="Momentum1" || VarNames[i]=="Momentum2" ||  VarNames[i]=="RawMomentum" ) 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=50.;
		      }
		    else if (VarNames[i]=="PTCorrectedMass" ) 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=10.;
		      }
		    
		    if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/bJets" ) )
		      {
			pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/bJets" ) ; 
			pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/bJets" ) ;
		      }
		    
		    _inputsHistogramsBJets[iInputCollection][VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		    
		    if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/cJets" ) )
		      {
			pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/cJets" ) ; 
			pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/cJets" ) ;
		      }
		    _inputsHistogramsCJets[iInputCollection][VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		    
		    
		    if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/udsJets" ) )
		      {
			pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/udsJets" ) ; 
			pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/udsJets" ) ;
		      }
		    _inputsHistogramsUDSJets[iInputCollection][VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		  
		  }//end of histogram creation
	      }
	    
	    _InputsIndex.push_back(IndexOf);
	    

	    if (isFirstEvent()) {
	      //We'd like to make zoomed histograms of some of the flavour tag inputs too
	      for (size_t i = 0;i < _ZoomedVarNames.size();++i) {
		
		std::string zoomed_name = _ZoomedVarNames[i] + " (zoomed)";
		
		if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/bJets" ) )
		  {
		    pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/bJets" ) ; 
		    pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/bJets" ) ;
		    
		  }			 			  
		_zoomedInputsHistogramsBJets[iInputCollection][zoomed_name] = pHistogramFactory->createHistogram1D( zoomed_name, 100, -10., 20.);
		
		if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/cJets" ) )
		  {
		    pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/cJets" ) ; 
		    pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/cJets" ) ;
		  }
		
		_zoomedInputsHistogramsCJets[iInputCollection][zoomed_name] = pHistogramFactory->createHistogram1D( zoomed_name, 100, -10., 20.);
		
		
		if( !pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/udsJets" ) )
		  {
		    pTree->mkdirs( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/udsJets" ) ; 
		    pTree->cd(    "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection] + "/udsJets" ) ;
		  }
		_zoomedInputsHistogramsUDSJets[iInputCollection][zoomed_name] = pHistogramFactory->createHistogram1D( zoomed_name, 100, -10., 20.);
	      }
	    }
	  }
}

 
void LCFIAIDAPlotProcessor::processEvent( LCEvent* pEvent ) 
{ 
  
  // Make sure that "processRunHeader" has been called for this run (see the comment in that method).
  if( (_lastRunHeaderProcessed != pEvent->getRunNumber()) && (_suppressOutputForRun != pEvent->getRunNumber()) )
    {
      std::cerr << __FILE__ << "(" << __LINE__ << "): processRunHeader() was not called for run " << pEvent->getRunNumber()
		<< " (did you use \"SkipNEvents\"?). The order of the information in the flavour tag collection(s) is going to be guessed." << std::endl;
      
      //Only want to do this once for this run, so set a marker that this run has been done
      _suppressOutputForRun=pEvent->getRunNumber();
      
      // Just assume that the elements are in the order "BTag", "CTag", "BCTag"
      std::map<std::string,unsigned int> guessedOrder;
      guessedOrder["BTag"]=0; guessedOrder["CTag"]=1; guessedOrder["BCTag"]=2;
      _IndexOfForEachTag.clear();
      for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag) _IndexOfForEachTag.push_back( guessedOrder );
    }
  
  
  // Try and get the jet collection. If unable, show why. No need to worry about quitting
  // because getNumberOfElements will return zero, so flow will never go into the loop.
  // TypesafeCollection is just a wrapper around LCCollection with more error checking and
  // things. Just makes the code a bit easier to read (in my opinion).

  TypesafeCollection<lcio::ReconstructedParticle> jetCollection( pEvent, _JetCollectionName );
  //	if( !jetCollection.is_valid() ) _log->message<marlin::ERROR>( jetCollection.last_error() );
  
  
  
  //apply any cuts on the event here
  if( PassesEventCuts(pEvent) )
    {

      ReconstructedParticle* pJet;
      //loop over the jets
      for( int jetNumber=0; jetNumber<jetCollection.getNumberOfElements(); ++jetNumber )
	{
	  pJet=jetCollection.getElementAt(jetNumber);
	  
	  //only do anything if the jet passes the jet cuts
	  if( PassesJetCuts(pJet) )
	    {
	      FillTagPlots( pEvent, jetNumber );
	      FillInputsPlots( pEvent, jetNumber );
	    }
	}
   
      

      LCCollection* vertexCol = pEvent->getCollection(_VertexColName);

      float primaryVertexPostion[] = {0.,0.,0.};

      for (int ii=0 ; ii < vertexCol->getNumberOfElements() ; ii++){
	
	Vertex* myVertex =  dynamic_cast<Vertex*>(vertexCol->getElementAt(ii));

	//	std::cout << "vertex is primary: " << myVertex->isPrimary() << " asoc particle: " << myVertex->getAssociatedParticle() << std::endl;
	
	const float* vertexPosition = myVertex->getPosition();
	if (myVertex->isPrimary()) {
	  primaryVertexPostion[0] = vertexPosition[0];
	  primaryVertexPostion[1] = vertexPosition[1];
	  primaryVertexPostion[2] = vertexPosition[2];
	}
      }

      for (int ii=0 ; ii < vertexCol->getNumberOfElements() ; ii++){

	Vertex* myVertex =  dynamic_cast<Vertex*>(vertexCol->getElementAt(ii));
	const float* vertexPosition = myVertex->getPosition();

	double px =  double(vertexPosition[0]);
	double py =  double(vertexPosition[1]);
	double pz =  double(vertexPosition[2]);
	double ex = sqrt((double)myVertex->getCovMatrix()[0]);	  
	double ey = sqrt((double)myVertex->getCovMatrix()[2]);
	double ez = sqrt((double)myVertex->getCovMatrix()[5]);     

	
	if (! myVertex->isPrimary() ) {

	  double distanceIP = sqrt(pow(px-primaryVertexPostion[0],2)+pow(py-primaryVertexPostion[1],2)+pow(pz-primaryVertexPostion[1],2));
	  	  
	  _pVertexDistanceFromIP->fill(distanceIP);
	  _pVertexPositionX->fill(px);
	  _pVertexPositionY->fill(py);
	  _pVertexPositionZ->fill(pz);
	  
	} else { //it is the primary vertex
	  
	  _pPrimaryVertexPositionX->fill(px);
	  _pPrimaryVertexPositionY->fill(py);
	  _pPrimaryVertexPositionZ->fill(pz);

	  _pPrimaryVertexPullX->fill(px/ex);
	  _pPrimaryVertexPullY->fill(py/ey);
	  _pPrimaryVertexPullZ->fill(pz/ez);
	}

      }

    }
}

 
 
void LCFIAIDAPlotProcessor::check( LCEvent* pEvent ) 
{
}

void LCFIAIDAPlotProcessor::end() 
{
  
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::IDataPointSetFactory* pDataPointSetFactory=marlin::AIDAProcessor::dataPointSetFactory(this);
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  
  _pBJetBTagIntegral.resize( _FlavourTagCollectionNames.size() );     
  _pCJetBTagIntegral.resize( _FlavourTagCollectionNames.size() );	  
  _pLightJetBTagIntegral.resize( _FlavourTagCollectionNames.size() ); 
  _pBJetCTagIntegral.resize( _FlavourTagCollectionNames.size() );  	  
  _pCJetCTagIntegral.resize( _FlavourTagCollectionNames.size() );	  
  _pLightJetCTagIntegral.resize( _FlavourTagCollectionNames.size() ); 
  _pBJetBCTagIntegral.resize( _FlavourTagCollectionNames.size() ); 	  
  _pCJetBCTagIntegral.resize( _FlavourTagCollectionNames.size() );	  
  _pLightJetBCTagIntegral.resize( _FlavourTagCollectionNames.size() );
  

  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    {

      for (unsigned int iVertexCat=1;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	//sum over the different vertex catagories, this information goes into the "AnyNumberOfVertices" directory
      
	_pBJetBTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pBJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pBJetCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pBJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pBJetBCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pBJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCJetBTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pCJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCJetCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pCJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCJetBCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pCJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pLightJetBTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pLightJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pLightJetCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pLightJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pLightJetBCTag[iTagCollection][_VertexCatNames[0]] -> add(*_pLightJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
      }


      for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	//add up all the background values

	pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection]+ "/" +_NumVertexCatDir[iVertexCat]);
	
	_pBTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-B jets by B-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetBTag[iTagCollection][_VertexCatNames[iVertexCat]],*_pCJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-C jets by C-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetCTag[iTagCollection][_VertexCatNames[iVertexCat]],*_pBJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]); 
	_pBCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-C jets by BC-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]],*_pBJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
      }
       
    }
      

  //now calculate the efficiencies, leakage rate and purity
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    {
      
      for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {

	pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection]+ "/" +_NumVertexCatDir[iVertexCat] );
	
	std::string nvname = _VertexCatNames[iVertexCat];

	AIDA::IDataPointSet* _pBJetBTagEfficiency = CreateEfficiencyPlot( _pBJetBTag[iTagCollection][nvname] , pDataPointSetFactory->create("B-Tag efficiency  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetCTagEfficiency = CreateEfficiencyPlot( _pCJetCTag[iTagCollection][nvname] , pDataPointSetFactory->create("C-Tag efficiency  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetBCTagEfficiency = CreateEfficiencyPlot( _pCJetBCTag[iTagCollection][nvname] , pDataPointSetFactory->create("BC-Tag efficiency  ("+ nvname +")",2));
	
  
	
	_pBJetBTagIntegral[iTagCollection][nvname] =   	
	  CreateIntegralHistogram( _pBJetBTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pBJetBTag[iTagCollection][nvname]->axis().bins(),_pBJetBTag[iTagCollection][nvname]->axis().lowerEdge(),_pBJetBTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pCJetBTagIntegral[iTagCollection][nvname] =     
	  CreateIntegralHistogram( _pCJetBTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pCJetBTag[iTagCollection][nvname]->axis().bins(),_pCJetBTag[iTagCollection][nvname]->axis().lowerEdge(),_pCJetBTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pLightJetBTagIntegral[iTagCollection][nvname] = 
	  CreateIntegralHistogram( _pLightJetBTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pLightJetBTag[iTagCollection][nvname]->axis().bins(),_pLightJetBTag[iTagCollection][nvname]->axis().lowerEdge(),_pLightJetBTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pBJetCTagIntegral[iTagCollection][nvname] = 
	CreateIntegralHistogram( _pBJetCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pBJetCTag[iTagCollection][nvname]->axis().bins(),_pBJetCTag[iTagCollection][nvname]->axis().lowerEdge(),_pBJetCTag[iTagCollection][nvname]->axis().upperEdge()));
      
	_pCJetCTagIntegral[iTagCollection][nvname] =     
	CreateIntegralHistogram( _pCJetCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pCJetCTag[iTagCollection][nvname]->axis().bins(),_pCJetCTag[iTagCollection][nvname]->axis().lowerEdge(),_pCJetCTag[iTagCollection][nvname]->axis().upperEdge()));

	_pLightJetCTagIntegral[iTagCollection][nvname] = 
      CreateIntegralHistogram( _pLightJetCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pLightJetCTag[iTagCollection][nvname]->axis().bins(),_pLightJetCTag[iTagCollection][nvname]->axis().lowerEdge(),_pLightJetCTag[iTagCollection][nvname]->axis().upperEdge()));
  
	_pBJetBCTagIntegral[iTagCollection][nvname] =    
	  CreateIntegralHistogram( _pBJetBCTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pBJetBCTag[iTagCollection][nvname]->axis().bins(),_pBJetBCTag[iTagCollection][nvname]->axis().lowerEdge(),_pBJetBCTag[iTagCollection][nvname]->axis().upperEdge()));

	_pCJetBCTagIntegral[iTagCollection][nvname] =   
	  CreateIntegralHistogram( _pCJetBCTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pCJetBCTag[iTagCollection][nvname]->axis().bins(),_pCJetBCTag[iTagCollection][nvname]->axis().lowerEdge(),_pCJetBCTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pLightJetBCTagIntegral[iTagCollection][nvname] = 
	  CreateIntegralHistogram( _pLightJetBCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",
								    _pLightJetBCTag[iTagCollection][nvname]->axis().bins(),_pLightJetBCTag[iTagCollection][nvname]->axis().lowerEdge(),_pLightJetBCTag[iTagCollection][nvname]->axis().upperEdge()));
	
	//Examples of the integral plots - instead of histograms - the histogram calculate the errors wrongly
	
	//integralplots	_pBJetBTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pBJetBTag[iTagCollection][nvname], pDataPointSetFactory->create("B-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pCJetBTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pCJetBTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pLightJetBTagIntegral[iTagCollection][nvname] = 
	//integralplots//	  CreateIntegralPlot( _pLightJetBTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pBJetCTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pBJetCTag[iTagCollection][nvname], pDataPointSetFactory->create("B-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pCJetCTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pCJetCTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pLightJetCTagIntegral[iTagCollection][nvname] = 
	//integralplots//	  CreateIntegralPlot( _pLightJetCTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pBJetBCTagIntegral[iTagCollection][nvname] =    
	//integralplots//	  CreateIntegralPlot( _pBJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("B-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pCJetBCTagIntegral[iTagCollection][nvname] =    
	//integralplots//	  CreateIntegralPlot( _pCJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pLightJetBCTagIntegral[iTagCollection][nvname] = 
	//integralplots//	  CreateIntegralPlot( _pLightJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",2));
	

	AIDA::IDataPointSet* _pBJetBTagPurity =  CreatePurityPlot( _pBJetBTag[iTagCollection][nvname],  _pBTagBackgroundValues[iTagCollection][nvname] , pDataPointSetFactory->create("B-Jet purity for B-Tag  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetCTagPurity =  CreatePurityPlot( _pCJetCTag[iTagCollection][nvname],  _pCTagBackgroundValues[iTagCollection][nvname] , pDataPointSetFactory->create("C-Jet purity for C-Tag  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetBCTagPurity = CreatePurityPlot( _pCJetBCTag[iTagCollection][nvname], _pBJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jet purity for BC-Tag  ("+ nvname +")",2));      
	
	AIDA::IDataPointSet* _pCJetBTagLeakage =      CreateLeakageRatePlot( _pCJetBTag[iTagCollection][nvname],      pDataPointSetFactory->create("C-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pLightJetBTagLeakage =  CreateLeakageRatePlot( _pLightJetBTag[iTagCollection][nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pBJetCTagLeakage =      CreateLeakageRatePlot( _pBJetCTag[iTagCollection][nvname],      pDataPointSetFactory->create("B-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pLightJetCTagLeakage =  CreateLeakageRatePlot( _pLightJetCTag[iTagCollection][nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pBJetBCTagLeakage =     CreateLeakageRatePlot( _pBJetBCTag[iTagCollection][nvname],     pDataPointSetFactory->create("B-Jets: Leakage Rate into BC-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pLightJetBCTagLeakage = CreateLeakageRatePlot( _pLightJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Leakage Rate into BC-Tag Sample  ("+ nvname +")",2));     
	AIDA::IDataPointSet* _pNonBJetBTagLeakage =   CreateLeakageRatePlot( _pBTagBackgroundValues[iTagCollection][nvname],      pDataPointSetFactory->create("C-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pNonCJetCTagLeakage =   CreateLeakageRatePlot( _pCTagBackgroundValues[iTagCollection][nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pNonCJetBCTagLeakage =  CreateLeakageRatePlot( _pBCTagBackgroundValues[iTagCollection][nvname],      pDataPointSetFactory->create("B-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
	
	CreateXYPlot(_pBJetBTagEfficiency, _pBJetBTagPurity, pDataPointSetFactory->create("Purity-Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pCJetCTagPurity, pDataPointSetFactory->create("Purity-Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetBCTagEfficiency, _pCJetBCTagPurity, pDataPointSetFactory->create("Purity-Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	
	CreateXYPlot(_pBJetBTagEfficiency, _pNonBJetBTagLeakage,  pDataPointSetFactory->create("Leakage Rate-Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pNonCJetCTagLeakage,  pDataPointSetFactory->create("Leakage Rate-Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetBCTagEfficiency,_pNonCJetBCTagLeakage, pDataPointSetFactory->create("Leakage Rate-Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	
	CreateXYPlot(_pBJetBTagEfficiency, _pCJetBTagLeakage,  pDataPointSetFactory->create("Leakage Rate (C into B) vs Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pBJetBTagEfficiency, _pLightJetBTagLeakage,  pDataPointSetFactory->create("Leakage Rate (Light into B) Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pBJetCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (B into C) vs Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pLightJetCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (Light into C) Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetBCTagEfficiency, _pBJetBCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (B into BC) vs Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pLightJetBCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (Light into BC) Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	
      }
    }
      //create vertex charge leakage rate plots
  pTree->cd( "/" + name() + "/VertexChargePlots/");
  AIDA::IDataPointSet* pBJetVtxChargeDPS = pDataPointSetFactory->create("B-Jets: Vertex Charge Leakage",2);	
  AIDA::IDataPointSet* pCJetVtxChargeDPS = pDataPointSetFactory->create("C-Jets: Vertex Charge Leakage",2);
  
  CreateVertexChargeLeakagePlot(pBJetVtxChargeDPS, pCJetVtxChargeDPS);
  
  
  if (_PrintNeuralNetOutput) PrintNNOutput();

}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently applies no cuts at all*/
bool LCFIAIDAPlotProcessor::PassesEventCuts( LCEvent* pEvent )
{
  //
  // No event cuts at present
  //
  
  return true;

}
// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs

bool LCFIAIDAPlotProcessor::PassesJetCuts( ReconstructedParticle* pJet )
{
  //
  // This cut added on the suggestion of Sonja Hillert 12/Jan/07.
  //
  // Selects jets for which the cosine of the jet polar
  // angle theta for all jets is not too large.
  //
  // Make something that's easy to search for to track down erroneous cuts:
  // (too many bad experiences of long forgotten about cuts hiding somewhere)
  // GREPTAG_CUT : Jet cut on abs(cos(theta)) of jet axis
  //
  
  
  const double* jetMomentum=pJet->getMomentum(); 
  
  double momentumMagnitude = sqrt(pow(jetMomentum[0],2)+pow(jetMomentum[1],2)+pow(jetMomentum[2],2));

  double cosTheta = jetMomentum[2] / momentumMagnitude;
  if( fabs(cosTheta) < _CosThetaJetMin || fabs(cosTheta) > _CosThetaJetMax ) return false;
  
  if (momentumMagnitude > _PJetMax ||  momentumMagnitude < _PJetMin) return false;

  
  // If control gets to this point then the jet has passed
  return true;
}

void LCFIAIDAPlotProcessor::FillInputsPlots( LCEvent* pEvent, unsigned int jetNumber )
{
   
	int jetType=FindJetType( pEvent, jetNumber );
	if( jetType==0 ) return;

	int CQVtx =  FindCQVtx(pEvent, jetNumber);
	int BQVtx =  FindBQVtx(pEvent, jetNumber);

	for (unsigned int iInputsCollection=0; iInputsCollection < _FlavourTagInputsCollectionNames.size(); ++iInputsCollection )
	{
	  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _FlavourTagInputsCollectionNames[iInputsCollection] );
	  if( !inputsCollection.is_valid() )
	    {
	      //should send out an error here
	    }
	  else
	    {
	      //Do stuff...
	      lcio::LCFloatVec* pInputs=inputsCollection.getElementAt( jetNumber );
	      if( !pInputs )
		{
		}
	      else
		{
		  
		  //ok everything is okay with the data

		  if (_MakeTuple && iInputsCollection==0) {
		    
		    //this could probably be done automatically
		    
		    int  NumVertices = int((*pInputs)[_InputsIndex[iInputsCollection]["NumVertices"]]);
		    int  NumTracksInVertices = int((*pInputs)[_InputsIndex[iInputsCollection] ["NumTracksInVertices"]]);
		    float D0Significance1=(*pInputs)[_InputsIndex[iInputsCollection]     ["D0Significance1"]];
		    float D0Significance2=(*pInputs)[_InputsIndex[iInputsCollection]     ["D0Significance2"]];
		    float DecayLength=(*pInputs)[_InputsIndex[iInputsCollection]	       ["DecayLength"]];
		    float DecayLength_SeedToIP=(*pInputs)[_InputsIndex[iInputsCollection]["DecayLength(SeedToIP)"]];
		    float DecayLengthSignificance=(*pInputs)[_InputsIndex[iInputsCollection]  ["DecayLengthSignificance"]];
		    float JointProbRPhi=(*pInputs)[_InputsIndex[iInputsCollection]	    ["JointProbRPhi"]];
		    float JointProbZ=(*pInputs)[_InputsIndex[iInputsCollection]	       ["JointProbZ"]];
		    float Momentum1=(*pInputs)[_InputsIndex[iInputsCollection]	       ["Momentum1"]];
		    float Momentum2=(*pInputs)[_InputsIndex[iInputsCollection]	       ["Momentum2"]];
		    float PTCorrectedMass=(*pInputs)[_InputsIndex[iInputsCollection]    ["PTCorrectedMass"]];
		    float RawMomentum=(*pInputs)[_InputsIndex[iInputsCollection]	       ["RawMomentum"]];
		    float SecondaryVertexProbability=(*pInputs)[_InputsIndex[iInputsCollection]["SecondaryVertexProbability"]];
		    float Z0Significance1=(*pInputs)[_InputsIndex[iInputsCollection]     ["Z0Significance1"]];
		    float Z0Significance2=(*pInputs)[_InputsIndex[iInputsCollection]     ["Z0Significance2"]];
		    
		    _pMyTuple->fill( 0, jetType );
		    _pMyTuple->fill( 1, NumVertices );
		    _pMyTuple->fill( 2, NumTracksInVertices );
		    _pMyTuple->fill( 3, D0Significance1);
		    _pMyTuple->fill( 4, D0Significance2);
		    _pMyTuple->fill( 5, DecayLength);
		    _pMyTuple->fill( 6, DecayLength_SeedToIP);
		    _pMyTuple->fill( 7, DecayLengthSignificance);
		    _pMyTuple->fill( 8, JointProbRPhi);
		    _pMyTuple->fill( 9, JointProbZ);
		    _pMyTuple->fill( 10, Momentum1);
		    _pMyTuple->fill( 11, Momentum2);
		    _pMyTuple->fill( 12, PTCorrectedMass);
		    _pMyTuple->fill( 13, RawMomentum);
		    _pMyTuple->fill( 14, SecondaryVertexProbability);
		    _pMyTuple->fill( 15, Z0Significance1);
		    _pMyTuple->fill( 16 ,Z0Significance2);

		    _pMyTuple->fill( 17, BQVtx );
		    _pMyTuple->fill( 18, CQVtx );
		    
		    _pMyTuple->addRow();
		  }
		  
		  
		  for( std::map<std::string,unsigned int>::iterator iTagNames=_InputsIndex[iInputsCollection].begin(); 
		       iTagNames!=_InputsIndex[iInputsCollection].end(); ++iTagNames ) {

		    double input=(*pInputs)[(*iTagNames).second]; 

		    //if the quantity relates to the second vertex, and there is no second vertex, then don't plot it		    
		    if (! ((*pInputs)[_InputsIndex[iInputsCollection]["NumVertices"]] < 2 && 
			((*iTagNames).first == "DecayLength" || (*iTagNames).first == "RawMomentum"  ||
			 (*iTagNames).first == "SecondaryVertexProbability" || (*iTagNames).first == "PTCorrectedMass" ||
			 (*iTagNames).first == "DecayLength(SeedToIP)" || (*iTagNames).first == "DecayLengthSignificance") )) {
		      
		      if( jetType==B_JET ) _inputsHistogramsBJets[iInputsCollection][(*iTagNames).first]->fill(input);
		      else if( jetType==C_JET ) _inputsHistogramsCJets[iInputsCollection][(*iTagNames).first]->fill(input);
		      else _inputsHistogramsUDSJets[iInputsCollection][(*iTagNames).first]->fill(input);
		    }
		    
		    //fill a few extra histograms created by hand
		    for (std::vector<std::string>::const_iterator iter = _ZoomedVarNames.begin() ; iter != _ZoomedVarNames.end(); iter++){
		      if ((*iTagNames).first == *iter) {
			std::string zoomed_name = (*iTagNames).first + " (zoomed)";
			if( jetType==B_JET ) _zoomedInputsHistogramsBJets[iInputsCollection][zoomed_name]->fill(input);
			else if( jetType==C_JET ) _zoomedInputsHistogramsCJets[iInputsCollection][zoomed_name]->fill(input);
			else _zoomedInputsHistogramsUDSJets[iInputsCollection][zoomed_name]->fill(input);
			
		      }
		    }
		    
		  }
		  
		}
	    }
	}
}

void LCFIAIDAPlotProcessor::FillTagPlots( LCEvent* pEvent, unsigned int jetNumber)
{
  int jetType=FindJetType( pEvent, jetNumber );
  if( jetType==0 ) return;
  
  TypesafeCollection<lcio::ReconstructedParticle> jetCollection( pEvent, _JetCollectionName );
  ReconstructedParticle* pJet;
  pJet=jetCollection.getElementAt(jetNumber);
  
  
  
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection)
    {
      TypesafeCollection<lcio::LCFloatVec> tagCollection( pEvent, _FlavourTagCollectionNames[iTagCollection] );
      if( !tagCollection.is_valid() )
	{
	}
      else
	{
	  lcio::LCFloatVec* pJetFlavourTags=tagCollection.getElementAt( jetNumber );
	  if( !pJetFlavourTags )
	    {
	    }
	  else
	    {
	      
	      const double* jetMomentum = pJet->getMomentum();
	      double cosTheta = jetMomentum[2] / sqrt(pow(jetMomentum[0],2)+pow(jetMomentum[1],2)+pow(jetMomentum[2],2));
	      
	      double bTag= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["BTag"]];
	      double cTag= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["CTag"]];
	      double cTagBBack= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["BCTag"]];
	      unsigned int NumVertices = FindNumVertex(pEvent, jetNumber, iTagCollection);
	      int CQVtx =  FindCQVtx(pEvent, jetNumber);
	      int BQVtx =  FindBQVtx(pEvent, jetNumber);
	      int trueJetCharge = int(FindJetHadronCharge(pEvent,jetNumber));
	      
	      std::string nvname = _VertexCatNames[ (NumVertices>=N_VERTEX_CATEGORIES) ? (N_VERTEX_CATEGORIES) : (NumVertices)];
	      
	      if( jetType==B_JET )  {

		if( bTag<=1 && bTag>=0 )
		  {
		      _pBJetBTag[iTagCollection][nvname]->fill( bTag );
		    } 
		  else 
		    {
		      _pBJetBTag[iTagCollection][nvname]->fill( -0.005 );
		    }
		  
		  if( cTag<=1 && cTag>=0 )
		    {
		      _pBJetCTag[iTagCollection][nvname]->fill( cTag );
		    }
		  else 
		    {
		      _pBJetCTag[iTagCollection][nvname]->fill( -0.005 );
		    }
		  if( cTagBBack<=1 && cTagBBack>=0 ) 
		    {
		      _pBJetBCTag[iTagCollection][nvname]->fill( cTagBBack );
		    }
		  else 
		    {
		      _pBJetBCTag[iTagCollection][nvname]->fill( -0.005 );
		    }
		
	      } else if( jetType==C_JET ) {
		
		if( bTag<=1 && bTag>=0 )
		  {
		    _pCJetBTag[iTagCollection][nvname]->fill( bTag );
		  } 
		else 
		  {
		    _pCJetBTag[iTagCollection][nvname]->fill( -0.005 );
		  }
		
		if( cTag<=1 && cTag>=0 )
		  {
		      _pCJetCTag[iTagCollection][nvname]->fill( cTag );
		  }
		else 
		  {
		    _pCJetCTag[iTagCollection][nvname]->fill( -0.005 );
		  }
		
		if( cTagBBack<=1 && cTagBBack>=0 ) 
		  {
		    _pCJetBCTag[iTagCollection][nvname]->fill( cTagBBack );
		  }
		else 
		  {
		    _pCJetBCTag[iTagCollection][nvname]->fill( -0.005 );
		  }
	      } else {
		if( bTag<=1 && bTag>=0 )
		  {
		    _pLightJetBTag[iTagCollection][nvname]->fill( bTag );
		  } 
		else 
		  {
		    _pLightJetBTag[iTagCollection][nvname]->fill( -0.005 );
		  }
		
		if( cTag<=1 && cTag>=0 )
		  {
		    _pLightJetCTag[iTagCollection][nvname]->fill( cTag );
		  }
		else 
		  {
		    _pLightJetCTag[iTagCollection][nvname]->fill( -0.005 );
		  } 
		
		if( cTagBBack<=1 && cTagBBack>=0 ) 
		  {
		    _pLightJetBCTag[iTagCollection][nvname]->fill( cTagBBack );
		  }
		else 
		  {
		    _pLightJetBCTag[iTagCollection][nvname]->fill( -0.005 );
		  }
	      }
	      
	      if (iTagCollection == _myVertexChargeTagCollection) {

		//vertex charge plots
		if( jetType==C_JET && cTag > _CTagNNCut) {
		  
		  int bin = _pCJetLeakageRate->coordToIndex(fabs(cosTheta));
		  
		  if (trueJetCharge==+2)   _cJet_truePlus2++;
		  if (trueJetCharge==+1)   _cJet_truePlus++;
		  if (trueJetCharge==0)    _cJet_trueNeut++;
		  if (trueJetCharge==-1)   _cJet_trueMinus++;
		  if (trueJetCharge==-2)   _cJet_trueMinus2++;
		  
		  if (trueJetCharge==+2){ 
		    if(CQVtx>0)  _cJet_truePlus2_recoPlus++; 
		    if(CQVtx==0) _cJet_truePlus2_recoNeut++;
		    if(CQVtx<0)  _cJet_truePlus2_recoMinus++;
		  }
		  if (trueJetCharge==+1){ 
		    if(CQVtx>0)  _cJet_truePlus_recoPlus++; 
		    if(CQVtx==0) _cJet_truePlus_recoNeut++;
		    if(CQVtx<0)  _cJet_truePlus_recoMinus++;
		  }
		  if (trueJetCharge==0) { 
		    if(CQVtx>0)  _cJet_trueNeut_recoPlus++; 
		    if(CQVtx==0) _cJet_trueNeut_recoNeut++;
		    if(CQVtx<0)  _cJet_trueNeut_recoMinus++;
		  }
		  if (trueJetCharge==-1)  { 
		    if(CQVtx>0)  _cJet_trueMinus_recoPlus++; 
		    if(CQVtx==0) _cJet_trueMinus_recoNeut++;
		    if(CQVtx<0)  _cJet_trueMinus_recoMinus++;
		  }
		  if (trueJetCharge==-2) { 
		    if(CQVtx>0)  _cJet_trueMinus2_recoPlus++; 
		    if(CQVtx==0) _cJet_trueMinus2_recoNeut++;
		    if(CQVtx<0)  _cJet_trueMinus2_recoMinus++;
		  }
		  
		  _pCJetVertexCharge->fill(CQVtx);
		  _pCJetCharge2D->fill(trueJetCharge,CQVtx);

		  if (trueJetCharge==+2)   _cJet_truePlus2_angle[bin]++;
		  if (trueJetCharge==+1)   _cJet_truePlus_angle[bin]++;
		  if (trueJetCharge==0)    _cJet_trueNeut_angle[bin]++;
		  if (trueJetCharge==-1)   _cJet_trueMinus_angle[bin]++;
		  if (trueJetCharge==-2)   _cJet_trueMinus2_angle[bin]++;
		  
		  if (trueJetCharge==+2){ 
		    if(CQVtx>0)  _cJet_truePlus2_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_truePlus2_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_truePlus2_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==+1){ 
		    if(CQVtx>0)  _cJet_truePlus_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_truePlus_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_truePlus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==0) { 
		    if(CQVtx>0)  _cJet_trueNeut_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_trueNeut_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_trueNeut_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==-1)  { 
		    if(CQVtx>0)  _cJet_trueMinus_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_trueMinus_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_trueMinus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==-2) { 
		    if(CQVtx>0)  _cJet_trueMinus2_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_trueMinus2_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_trueMinus2_recoMinus_angle[bin]++;
		  }    
		  
		} else if ( jetType==B_JET && bTag > _BTagNNCut) {
		  
		  int bin = _pBJetLeakageRate->coordToIndex(fabs(cosTheta));
		  
		  if (trueJetCharge==+2)   _bJet_truePlus2++;
		  if (trueJetCharge==+1)   _bJet_truePlus++;
		  if (trueJetCharge==0)    _bJet_trueNeut++;
		  if (trueJetCharge==-1)   _bJet_trueMinus++;
		  if (trueJetCharge==-2)   _bJet_trueMinus2++;
		  
		  if (trueJetCharge==+2){ 
		    if(BQVtx>0)  _bJet_truePlus2_recoPlus++; 
		    if(BQVtx==0) _bJet_truePlus2_recoNeut++;
		    if(BQVtx<0)  _bJet_truePlus2_recoMinus++;
		  }
		  if (trueJetCharge==+1){ 
		    if(BQVtx>0)  _bJet_truePlus_recoPlus++; 
		    if(BQVtx==0) _bJet_truePlus_recoNeut++;
		    if(BQVtx<0)  _bJet_truePlus_recoMinus++;
		  }
		  if (trueJetCharge==0) { 
		    if(BQVtx>0)  _bJet_trueNeut_recoPlus++; 
		    if(BQVtx==0) _bJet_trueNeut_recoNeut++;
		    if(BQVtx<0)  _bJet_trueNeut_recoMinus++;
		  }
		  if (trueJetCharge==-1)  { 
		    if(BQVtx>0)  _bJet_trueMinus_recoPlus++; 
		    if(BQVtx==0) _bJet_trueMinus_recoNeut++;
		    if(BQVtx<0)  _bJet_trueMinus_recoMinus++;
		  }
		  if (trueJetCharge==-2) { 
		    if(BQVtx>0)  _bJet_trueMinus2_recoPlus++; 
		    if(BQVtx==0) _bJet_trueMinus2_recoNeut++;
		    if(BQVtx<0)  _bJet_trueMinus2_recoMinus++;
		  }
		  
		  
		  if (trueJetCharge==+2) _bJet_truePlus2_angle[bin]++;
		  if (trueJetCharge==+1) _bJet_truePlus_angle[bin]++;	
		  if (trueJetCharge==0)  _bJet_trueNeut_angle[bin]++;	
		  if (trueJetCharge==-1) _bJet_trueMinus_angle[bin]++;	
		  if (trueJetCharge==-2) _bJet_trueMinus2_angle[bin]++;
		  
		  if (trueJetCharge==+2){ 
		    if(BQVtx>0)  _bJet_truePlus2_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_truePlus2_recoNeut_angle[bin]++;
		    if(BQVtx<0)  _bJet_truePlus2_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==+1){ 
		    if(BQVtx>0)  _bJet_truePlus_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_truePlus_recoNeut_angle[bin]++;
		    if(BQVtx<0)  _bJet_truePlus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==0) { 
		    if(BQVtx>0) _bJet_trueNeut_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_trueNeut_recoNeut_angle[bin]++;
		    if(BQVtx<0) _bJet_trueNeut_recoMinus_angle[bin]++;
		}
		  if (trueJetCharge==-1)  { 
		    if(BQVtx>0) _bJet_trueMinus_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_trueMinus_recoNeut_angle[bin]++;
		    if(BQVtx<0) _bJet_trueMinus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==-2) { 
		    if(BQVtx>0) _bJet_trueMinus2_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_trueMinus2_recoNeut_angle[bin]++;
		    if(BQVtx<0) _bJet_trueMinus2_recoMinus_angle[bin]++;
		  }
		  
		  _pBJetVertexCharge->fill(BQVtx);
		  _pBJetCharge2D->fill(trueJetCharge,BQVtx);
		}
	      }
	    }
	}
    }
}





int LCFIAIDAPlotProcessor::FindJetPDGCode( LCEvent* pEvent, unsigned int jetNumber )
{
	TypesafeCollection<lcio::LCIntVec> trueJetPDGCodeCollection( pEvent, _TrueJetPDGCodeColName );
	if( !trueJetPDGCodeCollection.is_valid() )
	{
	  std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " <<  _TrueJetPDGCodeColName << " is not valid " << std::endl;
	  return 0; //can't do anything without this collection
	}

	int pdgCode;
	lcio::LCIntVec* pTruePDGCodeVector=trueJetPDGCodeCollection.getElementAt( jetNumber );
	if( pTruePDGCodeVector )
	{
		if( pTruePDGCodeVector->size()==1 ) pdgCode=pTruePDGCodeVector->back();
		else
		{
		  std::cerr << __FILE__ << "(" << __LINE__ << "): The LCIntVec for jet " << jetNumber << " from the collection "
			    << _TrueJetFlavourColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
			    << " is not of size 1." << std::endl;
		  return 0; //can't fill any plots if we don't know the true flavour
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCIntVec for jet " << jetNumber << " from the collection " << _TrueJetFlavourColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return pdgCode;
}

float LCFIAIDAPlotProcessor::FindJetPartonCharge( LCEvent* pEvent, unsigned int jetNumber )
{
	TypesafeCollection<lcio::LCFloatVec> trueJetPartonChargeCollection( pEvent, _TrueJetPartonChargeColName );
	if( !trueJetPartonChargeCollection.is_valid() )
	{
	  std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " <<  _TrueJetPartonChargeColName << " is not valid " << std::endl;
		return 0; //can't do anything without this collection
	}

	float partonCharge;
	lcio::LCFloatVec* pTruePartonChargeVector=trueJetPartonChargeCollection.getElementAt( jetNumber );
	if( pTruePartonChargeVector )
	{
		if( pTruePartonChargeVector->size()==1 ) partonCharge=pTruePartonChargeVector->back();
		else
		{
			std::cerr << __FILE__ << "(" << __LINE__ << "): The LCFloatVec for jet " << jetNumber << " from the collection "
					<< _TrueJetPartonChargeColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
					<< " is not of size 1." << std::endl;
			return 0; //can't fill any plots if we don't know the true flavour
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCFloatVec for jet " << jetNumber << " from the collection " << _TrueJetPartonChargeColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return partonCharge;
}


int LCFIAIDAPlotProcessor::FindJetType( LCEvent* pEvent, unsigned int jetNumber )
{
	TypesafeCollection<lcio::LCIntVec> trueJetFlavourCollection( pEvent, _TrueJetFlavourColName );
	if( !trueJetFlavourCollection.is_valid() )
	  {
	    std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " <<  _TrueJetFlavourColName << " is not valid " << std::endl;
		return 0; //can't do anything without this collection
	}

	int jetType;
	lcio::LCIntVec* pTrueJetTypeVector=trueJetFlavourCollection.getElementAt( jetNumber );
	if( pTrueJetTypeVector )
	{
		if( pTrueJetTypeVector->size()==1 ) jetType=pTrueJetTypeVector->back();
		else
		{
			std::cerr << __FILE__ << "(" << __LINE__ << "): The LCIntVec for jet " << jetNumber << " from the collection "
					<< _TrueJetFlavourColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
					<< " is not of size 1." << std::endl;
			return 0; //can't fill any plots if we don't know the true flavour
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCIntVec for jet " << jetNumber << " from the collection " << _TrueJetFlavourColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return jetType;
}

float LCFIAIDAPlotProcessor::FindJetHadronCharge( LCEvent* pEvent, unsigned int jetNumber )
{
  TypesafeCollection<lcio::LCFloatVec> trueJetHadronChargeCollection( pEvent, _TrueJetHadronChargeColName );
  if( !trueJetHadronChargeCollection.is_valid() )
    {
      std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " << _TrueJetHadronChargeColName << " is not valid " << std::endl;
      return -99; //can't do anything without this collection
    }

	float hadronCharge;
	lcio::LCFloatVec* pTrueJetChargeVector=trueJetHadronChargeCollection.getElementAt( jetNumber );
	if( pTrueJetChargeVector )
	{
		if( pTrueJetChargeVector->size()==1 ) hadronCharge=pTrueJetChargeVector->back();
		else
		{
			std::cerr << __FILE__ << "(" << __LINE__ << "): The LCFloatVec for jet " << jetNumber << " from the collection "
					<< _TrueJetHadronChargeColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
					<< " is not of size 1." << std::endl;
			return -99; 
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCFloatVec for jet " << jetNumber << " from the collection " << _TrueJetHadronChargeColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return hadronCharge;
}

int LCFIAIDAPlotProcessor::FindNumVertex( LCEvent* pEvent, unsigned int jetNumber, unsigned int iInputsCollection)
{

  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _FlavourTagInputsCollectionNames[iInputsCollection] );
  if( !inputsCollection.is_valid() )
    {
      //		_log->message<marlin::ERROR>( trueJetFlavourCollection.last_error() );
      return 0; //can't do anything without this collection 
    }
  else
    {
      //Do stuff...
      lcio::LCFloatVec* pInputs=inputsCollection.getElementAt( jetNumber );
      if( !pInputs )
	{
	}
      else
	{
	 return  int((*pInputs)[_InputsIndex[iInputsCollection]["NumVertices"]]);
	}
      
    }
  return 0;
}

int LCFIAIDAPlotProcessor::FindBQVtx( LCEvent* pEvent, unsigned int jetNumber) 
{
  
  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _BVertexChargeCollection);
  
  if( !inputsCollection.is_valid() )  {
    
    std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection " << _BVertexChargeCollection << "  for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " BQVtx will be invalid" << std::endl;
    return -99;
    
  } else { //inputsCollection.is_valid() 
    
    float bqvtx;
    lcio::LCFloatVec* pBVtxChargeVector =inputsCollection.getElementAt( jetNumber );
    
    //bool evaluation is done left to right...
    if( pBVtxChargeVector && pBVtxChargeVector->size() == 1) {
      
      bqvtx = pBVtxChargeVector->back();
      
    } else {
      
      std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection element in  " << _BVertexChargeCollection << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " corresponding to jet number: " << jetNumber << " BQVtx will be invalid" << std::endl;
      return -99;
    }
    
    return int(bqvtx);
  }
  //should never get here
  return -99;
}


int LCFIAIDAPlotProcessor::FindCQVtx( LCEvent* pEvent, unsigned int jetNumber) 
{
  
  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _CVertexChargeCollection);
  
  if( !inputsCollection.is_valid() )  {
    
    std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection " << _CVertexChargeCollection << "  for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " CQVtx will be invalid" << std::endl;
    return -99;
    
  } else { //inputsCollection.is_valid() 
    
    float cqvtx;
    lcio::LCFloatVec* pCVtxChargeVector =inputsCollection.getElementAt( jetNumber );
    
    //bool evaluation is done left to right...
    if( pCVtxChargeVector && pCVtxChargeVector->size() == 1) {
      
      cqvtx = pCVtxChargeVector->back();
      
    } else {
      
      std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection element in  " << _CVertexChargeCollection << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " corresponding to jet number: " << jetNumber << " CQVtx will be invalid" << std::endl;
      return -99;
    }
    
    return int(cqvtx);
  }
  //should never get here
  return -99;
}



AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateEfficiencyPlot(const AIDA::IHistogram1D* pSignal, AIDA::IDataPointSet* pDataPointSet)
{
  
  double totalSignal=pSignal->sumBinHeights();
  double signalPassedCut=0;
    
  const int numberOfBins=pSignal->axis().bins();
  int iPoint=0;
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber, iPoint++ )
    {
      signalPassedCut+=pSignal->binHeight( binNumber );
      
      double efficiency = signalPassedCut/totalSignal;
      
      double efficiencyError = efficiency * (1. - efficiency) / totalSignal;
      if (efficiencyError>0) efficiencyError = sqrt(efficiencyError);
      
      
      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pSignal->axis().binLowerEdge(binNumber)+pSignal->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue( efficiency );
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(efficiencyError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(efficiencyError);
    }
  
  return pDataPointSet;
}  


AIDA::IHistogram1D* LCFIAIDAPlotProcessor::CreateIntegralHistogram(const AIDA::IHistogram1D* pNN, AIDA::IHistogram1D* pIntegral)
{
  //the errors on these entries are wrong...
  
  const int numberOfBins=pNN->axis().bins();

  double integral=pNN->binHeight(AIDA::IAxis::OVERFLOW_BIN);
  pIntegral->fill(pNN->axis().binLowerEdge(AIDA::IAxis::OVERFLOW_BIN)+pNN->axis().binWidth(numberOfBins-1)/2.,integral);
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber )
    {
      integral+= pNN->binHeight( binNumber );
      pIntegral->fill( pNN->axis().binLowerEdge(binNumber)+pNN->axis().binWidth(binNumber)/2.,integral);
    }
  
  integral+= pNN->binHeight(AIDA::IAxis::UNDERFLOW_BIN);
  pIntegral->fill(pNN->axis().binUpperEdge(AIDA::IAxis::UNDERFLOW_BIN)-pNN->axis().binWidth(0)/2.,integral);

  return pIntegral;
}



AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateIntegralPlot(const AIDA::IHistogram1D* pNN, AIDA::IDataPointSet* pDataPointSet )
{
  //it might make more sense to make this a histogram, but then the error on the entries are wrong

  const int numberOfBins=pNN->axis().bins();

  double integral=0; 
  
  for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ )  
    {
      
      integral+= pNN->binHeight( binNumber );
      pDataPointSet->addPoint();
      pDataPointSet->point(binNumber)->coordinate(0)->setValue(pNN->axis().binLowerEdge(binNumber)+pNN->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(binNumber)->coordinate(1)->setValue(integral);
      pDataPointSet->point(binNumber)->coordinate(1)->setErrorPlus(sqrt(integral));
      pDataPointSet->point(binNumber)->coordinate(1)->setErrorMinus(sqrt(integral));
    }
  return pDataPointSet;
}

AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreatePurityPlot(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet)  
{
  const int numberOfBins=pSignal->axis().bins();
  int iPoint=0;
  
  double signalPassedCut=0;
  double backgroundPassedCut=0;

  for (int binNumber = numberOfBins-1; binNumber >= 0 ; --binNumber, iPoint++ )  
    {
      
      signalPassedCut+=pSignal->binHeight( binNumber );
      backgroundPassedCut+=pBackground->binHeight( binNumber );
      
      double purity = signalPassedCut/(signalPassedCut+backgroundPassedCut);
      double purityError = purity * (1. - purity) / (signalPassedCut+backgroundPassedCut);
      if (purityError>0) purityError = sqrt(purityError);
 
      
      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pSignal->axis().binLowerEdge(binNumber)+pSignal->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue(purity);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(purityError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(purityError);

      
    }  

  return pDataPointSet;
}

AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateLeakageRatePlot(const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet)
{
  
  double totalBackground = pBackground->sumBinHeights();
  double backgroundPassedCut=0;
  
  const int numberOfBins=pBackground->axis().bins();
  int iPoint=0;
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber , iPoint++ )
    {

      backgroundPassedCut+=pBackground->binHeight( binNumber );
          
      double leakageRate = backgroundPassedCut/totalBackground;
    
      double leakageRateError = leakageRate * (1. - leakageRate) / totalBackground;
      if (leakageRateError>0) leakageRateError = sqrt(leakageRateError);
    

      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pBackground->axis().binLowerEdge(binNumber)+pBackground->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue(leakageRate);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(leakageRateError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(leakageRateError);
    }
  
  return pDataPointSet;
}



AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateXYPlot(const AIDA::IDataPointSet* pDataPointSet0, const AIDA::IDataPointSet* pDataPointSet1, AIDA::IDataPointSet* xyPointSet, const int dim0, const int dim1 )
{
  
  //need to do some comparision here
  if (pDataPointSet0->size() == pDataPointSet1->size()) {

    for (int iPoint = 0 ; iPoint != pDataPointSet1->size(); iPoint++) 
      {
	xyPointSet->addPoint();
	xyPointSet->point(iPoint)->coordinate(0)->setValue(pDataPointSet0->point(iPoint)->coordinate(dim0)->value());
	xyPointSet->point(iPoint)->coordinate(1)->setValue(pDataPointSet1->point(iPoint)->coordinate(dim1)->value());
	xyPointSet->point(iPoint)->coordinate(0)->setErrorPlus(pDataPointSet0->point(iPoint)->coordinate(dim0)->errorPlus());
	xyPointSet->point(iPoint)->coordinate(1)->setErrorPlus(pDataPointSet1->point(iPoint)->coordinate(dim1)->errorPlus());
	xyPointSet->point(iPoint)->coordinate(0)->setErrorMinus(pDataPointSet0->point(iPoint)->coordinate(dim0)->errorMinus());
	xyPointSet->point(iPoint)->coordinate(1)->setErrorMinus(pDataPointSet1->point(iPoint)->coordinate(dim1)->errorMinus());
      } 
  } else {
    
    //some error message here!
  }
  return xyPointSet;
}

void LCFIAIDAPlotProcessor::CreateVertexChargeLeakagePlot(AIDA::IDataPointSet* pBJetVtxChargeDPS, AIDA::IDataPointSet* pCJetVtxChargeDPS)
{
  
  for (int j = 0 ; j < N_JETANGLE_BINS ; j++) {
    
    double c_numerator   = _cJet_truePlus2_recoPlus_angle[j]+_cJet_truePlus_recoPlus_angle[j]
      +_cJet_trueMinus2_recoMinus_angle[j]+_cJet_trueMinus_recoMinus_angle[j];
    double c_domininator = _cJet_truePlus2_angle[j]         +_cJet_truePlus_angle[j]         
      +_cJet_trueMinus2_angle[j]          +_cJet_trueMinus_angle[j];
    
    double b_numerator   = _bJet_truePlus2_recoPlus_angle[j]+_bJet_truePlus_recoPlus_angle[j]
      +_bJet_trueMinus2_recoMinus_angle[j]+_bJet_trueMinus_recoMinus_angle[j];
    double b_domininator = _bJet_truePlus2_angle[j]         +_bJet_truePlus_angle[j]         
      +_bJet_trueMinus2_angle[j]          +_bJet_trueMinus_angle[j];
    
    double b_leakage = 1. - b_numerator/b_domininator;
    double c_leakage = 1. - c_numerator/c_domininator;
    
    double b_leakage_error = sqrt( b_leakage * (1. -  b_leakage) / b_domininator);
    double c_leakage_error = sqrt( c_leakage * (1. -  c_leakage) / c_domininator);
    
    
    pBJetVtxChargeDPS->addPoint();
    pBJetVtxChargeDPS->point(j)->coordinate(0)->setValue(_pBJetLeakageRate->axis().binLowerEdge(j)+_pBJetLeakageRate->axis().binWidth(j)/2.);
    pBJetVtxChargeDPS->point(j)->coordinate(1)->setValue(b_leakage);
    pBJetVtxChargeDPS->point(j)->coordinate(1)->setErrorPlus(b_leakage_error);
    pBJetVtxChargeDPS->point(j)->coordinate(1)->setErrorMinus(b_leakage_error);
    
    
    pCJetVtxChargeDPS->addPoint();
    pCJetVtxChargeDPS->point(j)->coordinate(0)->setValue(_pCJetLeakageRate->axis().binLowerEdge(j)+_pCJetLeakageRate->axis().binWidth(j)/2.);
    pCJetVtxChargeDPS->point(j)->coordinate(1)->setValue(c_leakage);
    pCJetVtxChargeDPS->point(j)->coordinate(1)->setErrorPlus(c_leakage_error);
    pCJetVtxChargeDPS->point(j)->coordinate(1)->setErrorMinus(c_leakage_error);
    
    _pCJetLeakageRate->fill(_pCJetLeakageRate->axis().binLowerEdge(j)+_pCJetLeakageRate->axis().binWidth(j)/2.,c_leakage); 
    _pBJetLeakageRate->fill(_pBJetLeakageRate->axis().binLowerEdge(j)+_pBJetLeakageRate->axis().binWidth(j)/2.,b_leakage);
  }
}


float  LCFIAIDAPlotProcessor::CalculateDistance(const float* pos1, const float* pos2){
  return sqrt(pow((pos1[0]-pos2[0]),2)+pow((pos1[1]-pos2[1]),2)+pow((pos1[2]-pos2[2]),2));
}



void LCFIAIDAPlotProcessor::PrintNNOutput(){
  
  //if there is a _NeuralNetOutputFile string defined use that as the output stream, if not use std::cout
  std::filebuf* fb = new std::filebuf;  
  
  std::ostream outputFile( (!_NeuralNetOutputFile.empty()) ?                                  
		       fb->open(_NeuralNetOutputFile.c_str(),
				std::ios_base::out|std::ios_base::trunc):  
			std::cout.rdbuf());

  if (outputFile.rdbuf() != fb)
      {
	delete fb;
	std::cerr << "Unable to open file " <<  _NeuralNetOutputFile << "!  Redirecting output to standard out." << std::endl;
	outputFile << std::endl;
      }

 for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection)
    {
      outputFile << "\n\nRESULTS for " << iTagCollection << "-th Flavour Tag Collection " << _FlavourTagCollectionNames[iTagCollection] << std::endl;
      outputFile << "---------------------------------------------------------------------------------------\n\n";
            
      outputFile << "1vtx N(b) = " ;
      outputFile.width(10);
      outputFile << _pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->sumAllBinHeights();
      outputFile << "   N(c) = " ;
      outputFile.width(10);
      outputFile << _pCJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->sumAllBinHeights();
      outputFile << "   N(light) = ";
      outputFile.width(10);
      outputFile << _pLightJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->sumAllBinHeights();
      outputFile << std::endl << std::endl;
  
      
      outputFile.precision(5);  
      outputFile.setf(std::ios::fixed,std::ios::floatfield);  

	  
      for (unsigned int iVertexCat=1;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	
	std::string nvname = _VertexCatNames[iVertexCat];
	
	//want to integrate!
	
	outputFile << "numbers of jets in cuts for ";
	outputFile.width(2);
	outputFile << iVertexCat;
	outputFile << " ZVTOP vertices found" << std::endl;
	outputFile << "cut    b-tag b    b-tag other    c-tag c   c-tag other    c-tagbb c   c-tagbb other" << std::endl;
	
	int numberOfBins=_pBJetBTagIntegral[iTagCollection][nvname]->axis().bins();
	
	for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ ) {
	  
	  outputFile.width(5);
	  outputFile.precision(3);
	  outputFile << _pBJetBTagIntegral[iTagCollection][nvname]->axis().binUpperEdge(binNumber) << "   ";
	  outputFile.width(10);
	  outputFile << int(_pBJetBTagIntegral[iTagCollection][nvname]->binHeight(binNumber))  << "   ";
	  outputFile.width(10);
	  outputFile << int(_pCJetBTagIntegral[iTagCollection][nvname]->binHeight(binNumber)+_pLightJetBTagIntegral[iTagCollection][nvname]->binHeight(binNumber))<< "   ";
	  outputFile.width(10);
	  outputFile << int(_pCJetCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))  << "   ";
	  outputFile.width(10);
	  outputFile << int(_pBJetCTagIntegral[iTagCollection][nvname]->binHeight(binNumber)+_pLightJetCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))<< "   ";
	  outputFile.width(10);
	  outputFile << int(_pCJetBCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))  << "   ";
	  outputFile.width(10);
	  outputFile << int(_pBJetBCTagIntegral[iTagCollection][nvname]->binHeight(binNumber)+_pLightJetBCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))<< "   ";
	  outputFile << std::endl;

	}
	outputFile << std::endl;
      }
      
       outputFile << "numbers of jets in cuts summed" << std::endl;
       outputFile << "cut    b-tag b    b-tag other    c-tag c   c-tag other    c-tagbb c   c-tagbb other" << std::endl;
  	
       int numberOfBins=_pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->axis().bins();
	
       for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ ) {
       
	 outputFile.width(5);
	 outputFile.precision(3);
	 outputFile << _pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->axis().binUpperEdge(binNumber) << "   ";
	 outputFile.width(10);
	 outputFile << int(_pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))  << "   ";
	 outputFile.width(10);
	 outputFile << int(_pCJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber)+_pLightJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))<< "   ";
	 outputFile.width(10);
	 outputFile << int(_pCJetCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))  << "   ";
	 outputFile.width(10);
	 outputFile << int(_pBJetCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber)+_pLightJetCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))<< "   ";
	 outputFile.width(10);
	 outputFile << int(_pCJetBCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))  << "   ";
	 outputFile.width(10);
	 outputFile << int(_pBJetBCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber)+_pLightJetBCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))<< "   ";
	 outputFile << std::endl; 
       }      
       outputFile << std::endl;
       outputFile << std::endl;
      
    }
}

void LCFIAIDAPlotProcessor::InternalVectorInitialisation()
{
  //sets the size of some vectors and initialises some counters to zero

  _cJet_truePlus2_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus2_angle.resize(N_JETANGLE_BINS);
  
  _cJet_truePlus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_truePlus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus2_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_truePlus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueNeut_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_trueNeut_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueNeut_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_trueMinus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_trueMinus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus2_recoMinus_angle.resize(N_JETANGLE_BINS);
  
  _bJet_truePlus2_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus_angle.resize(N_JETANGLE_BINS);	
  _bJet_trueNeut_angle.resize(N_JETANGLE_BINS);	
  _bJet_trueMinus_angle.resize(N_JETANGLE_BINS);	
  _bJet_trueMinus2_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_truePlus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus2_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_truePlus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_trueNeut_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_trueNeut_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_trueNeut_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_trueMinus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_trueMinus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus2_recoMinus_angle.resize(N_JETANGLE_BINS);
 

  for (int j = 0 ; j < N_JETANGLE_BINS ; j++) {
    _cJet_truePlus2_angle[j]=0;	       
    _cJet_truePlus_angle[j]=0;	       
    _cJet_trueNeut_angle[j]=0;	       
    _cJet_trueMinus_angle[j]=0;	       
    _cJet_trueMinus2_angle[j]=0;	       
    
    _cJet_truePlus2_recoPlus_angle[j]=0;  
    _cJet_truePlus2_recoNeut_angle[j]=0;  
    _cJet_truePlus2_recoMinus_angle[j]=0; 
    _cJet_truePlus_recoPlus_angle[j]=0;   
    _cJet_truePlus_recoNeut_angle[j]=0;   
    _cJet_truePlus_recoMinus_angle[j]=0;  
    _cJet_trueNeut_recoPlus_angle[j]=0;   
    _cJet_trueNeut_recoNeut_angle[j]=0;   
    _cJet_trueNeut_recoMinus_angle[j]=0;  
    _cJet_trueMinus_recoPlus_angle[j]=0;  
    _cJet_trueMinus_recoNeut_angle[j]=0;  
    _cJet_trueMinus_recoMinus_angle[j]=0; 
    _cJet_trueMinus2_recoPlus_angle[j]=0; 
    _cJet_trueMinus2_recoNeut_angle[j]=0; 
    _cJet_trueMinus2_recoMinus_angle[j]=0;
    
    _bJet_truePlus2_angle[j]=0;	       
    _bJet_truePlus_angle[j]=0;	       
    _bJet_trueNeut_angle[j]=0;	       
    _bJet_trueMinus_angle[j]=0;	       
    _bJet_trueMinus2_angle[j]=0;	       
    _bJet_truePlus2_recoPlus_angle[j]=0;  
    _bJet_truePlus2_recoNeut_angle[j]=0;  
    _bJet_truePlus2_recoMinus_angle[j]=0; 
    _bJet_truePlus_recoPlus_angle[j]=0;   
    _bJet_truePlus_recoNeut_angle[j]=0;   
    _bJet_truePlus_recoMinus_angle[j]=0;  
    _bJet_trueNeut_recoPlus_angle[j]=0;   
    _bJet_trueNeut_recoNeut_angle[j]=0;   
    _bJet_trueNeut_recoMinus_angle[j]=0;  
    _bJet_trueMinus_recoPlus_angle[j]=0;  
    _bJet_trueMinus_recoNeut_angle[j]=0;  
    _bJet_trueMinus_recoMinus_angle[j]=0; 
    _bJet_trueMinus2_recoPlus_angle[j]=0; 
    _bJet_trueMinus2_recoNeut_angle[j]=0; 
    _bJet_trueMinus2_recoMinus_angle[j]=0;
  }
  
  _cJet_truePlus2=0;
  _cJet_truePlus=0;
  _cJet_trueNeut=0;
  _cJet_trueMinus=0;
  _cJet_trueMinus2=0;
  
  _cJet_truePlus2_recoPlus=0; 
  _cJet_truePlus2_recoNeut=0;
  _cJet_truePlus2_recoMinus=0;
  _cJet_truePlus_recoPlus=0; 
  _cJet_truePlus_recoNeut=0;
  _cJet_truePlus_recoMinus=0;
  _cJet_trueNeut_recoPlus=0; 
  _cJet_trueNeut_recoNeut=0;
  _cJet_trueNeut_recoMinus=0;
  _cJet_trueMinus_recoPlus=0; 
  _cJet_trueMinus_recoNeut=0;
  _cJet_trueMinus_recoMinus=0;
  _cJet_trueMinus2_recoPlus=0; 
  _cJet_trueMinus2_recoNeut=0;
  _cJet_trueMinus2_recoMinus=0;
  
  _bJet_truePlus2=0;
  _bJet_truePlus=0;	
  _bJet_trueNeut=0;	
  _bJet_trueMinus=0;	
  _bJet_trueMinus2=0;
  _bJet_truePlus2_recoPlus=0; 
  _bJet_truePlus2_recoNeut=0;
  _bJet_truePlus2_recoMinus=0;
  _bJet_truePlus_recoPlus=0; 
  _bJet_truePlus_recoNeut=0;
  _bJet_truePlus_recoMinus=0;
  _bJet_trueNeut_recoPlus=0; 
  _bJet_trueNeut_recoNeut=0;
  _bJet_trueNeut_recoMinus=0;
  _bJet_trueMinus_recoPlus=0; 
  _bJet_trueMinus_recoNeut=0;
  _bJet_trueMinus_recoMinus=0;
  _bJet_trueMinus2_recoPlus=0; 
  _bJet_trueMinus2_recoNeut=0;
  _bJet_trueMinus2_recoMinus=0;
  
  return;
}



#endif // endif of the check to see if MARLIN_USE_AIDA has been defined
