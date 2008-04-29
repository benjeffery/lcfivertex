#ifndef LCFIAIDAPlotProcessor_h 
#define LCFIAIDAPlotProcessor_h 


//===============================================================================================================
// LCFIAIDAPlotProcessor Class - make plots of the LCFI flavour tag and vertex charge input and output variables
//=============================================================================================================== 


/** LCFIAIDAPlotProcessor Class - make plots of the LCFI flavour tag and vertex charge code.
 * <b> Please note that LCFIAIDAPlotProcessor will not compile with RAIDA v01-03 
 * To use LCFIAIDAPlotProcessor please use AIDAJNI for the implementation of AIDA
 * run  "cmake -DBUILD_WITH="ROOT;AIDAJNI" -DAIDAJNI_HOME=${AIDAJNI_HOME}"  </b>
 *
 * This sorry states of affairs comes about becayse not all AIDA functions are defined in RAIDA v01-03 
 *
 * In order to make a histogram file, LCFIAIDAPlotProcessor must be run with AIDAProcessor.<br>
 *
 * LCFIAIDAPlotProcessor reads in one (or more) FlavourTagCollections, e.g. from FlavourTag and one (or more) TagInputCollections.
 * Histograms/plots are made of the neural net outputs, the purity and leakage rate of the flavour tag.
 * These are split into sub-samples based on the number of vertices found in the jets.
 * Plots are also made of the inputs to the FlavourTagCollections - split into sub-samples based on
 * the true (MC) flavour of the jet.<br>
 *
 * Options are given to make a tuple of the flavour tag inputs and to print out a text file
 * of the different flavour tag neural net outputs. 
 * 
 * (When providing more than one FlavourTagCollection and/or TagInputCollection plots for each collection will be made in different directories.)
 *
 * In addition LCFIAIDAPlotProcessor also requires a jet collection, and the following collections, which should refer to the <i>same</i> jet collection.
 *<p>
 *  BVertexChargeCollection  -- calculated in VertexChargeProcessor<br>
 *  CVertexChargeCollection  -- calculated in VertexChargeProcessor<br>
 *  TrueJetFlavourCollection  -- calculated in TrueAngularJetFlavourProcessor<br>
 *  TrueJetHadronChargeCollection  -- calculated in TrueAngularJetFlavourProcessor<br>
 *  TrueJetPDGCodeCollection -- calculated in TrueAngularJetFlavourProcessor<br>
 *  TrueJetPartonChargeCollection -- calculated in TrueAngularJetFlavourProcessor<br>
 *     
 *
 * <H4>Input</H4>
 * 
 * The following collections must be available:
 *
 * @param FlavourTagCollections  StringVec of LCFloatVec names representing the flavour tag inputs collections - may be more than one collection.
 * @param TagInputsCollections   StringVec of LCFloatVec names the flavour tag input collections - may be more than one collection.
 * @param JetCollectionName   Name of ReconstructedParticleCollection representing the jets.
 * @param VertexCollectionName   Name of VertexCollection representing the Vertex collection of the jets.
 * @param BVertexChargeCollection Name of LCFloatVector of the vertex charge of the jet collection, assuming the jets are b-jets  (calculated in VertexChargeProcessor)
 * @param CVertexChargeCollection  Name of LCFloatVector of the vertex charge of the jet collection, assuming the jets are c-jets  (calculated in VertexChargeProcessor)
 * @param TrueJetFlavourCollection  Name of LCIntVector of the true (MC) flavour of the jets (calculated in TrueAngularJetFlavourProcessor)
 * @param TrueJetHadronChargeCollection  Name of LCFloatVector of the true (MC) charge of the hadron initiating the jets (calculated in TrueAngularJetFlavourProcessor)
 * @param TrueJetPDGCodeCollection   Name of LCIntVector of the true (MC) PDG code of the hadron initiating the jets (calculated in TrueAngularJetFlavourProcessor)
 * @param TrueJetPartonChargeCollection  Name of LCFloatVector of the true (MC) charge of the parton (quark) initiating the jets (calculated in TrueAngularJetFlavourProcessor)
 *
 * @param VertexCollectionName  Name of VertexCollection representing the vertices.
 * @param BTagNNCut Double reprsenting the lower cut on the b-tag NN value for some of the plots.
 * @param CTagNNCut Double reprsenting the lower cut on the c-tag NN value for some of the plots.
 * @param CosThetaJetMax  Double representing upper cut on cos(theta) of the jets for the plots.
 * @param CosThetaJetMin  Double representing lower cut on cos(theta) of the jets for the plots.
 * @param PJetMax   Double representing upper cut on momentum of the jet for the plots.
 * @param PJetMin   Double representing lower cut on momentum of the jet for the plots.
 * @param MakeTuple   Bool set true if you want to make a tuple of the TagInputCollection variables.
 * @param NeuralNetOutputFile  String representing name of text file of neural net values to.  
 *    Only used if PrintNeuralNetOutput parameter is true.  If left blank, output will be directed to standard out
 * @param PrintNeuralNetOutput  Bool set true if you want to make a text file of the neural net values (useful for some scripts).
 * @param UseFlavourTagCollectionForVertexCharge For vertex charge plots we demand the cTag>CTagNNCut and bTag>BTagNNCut.  This integer is used if there is more than one tag collection, to determine which of the collections should be used to apply this cut.
 * 
 * <H4>Output</H4>
 * - An aida (or root??) file containing the histograms, plots and tuples.
 * - (Optionally) a text file containing some of the neural net tagging output
 *  
 * @author Victoria Martin (victoria.martin@ed.ac.uk)
*/


//Marlin and LCIO includes 
#include "marlin/Processor.h" 
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/LCFloatVec.h"

#include <iostream>
#include <fstream>

//AIDA includes...
#include <AIDA/IHistogram1D.h>
#include <AIDA/IDataPointSet.h>
#include <AIDA/ITuple.h>




class LCFIAIDAPlotProcessor : public marlin::Processor
{ 
public: 
	//this bit has to be done for all Marlin processors 
	virtual Processor* newProcessor() { return new LCFIAIDAPlotProcessor; } 
 
	LCFIAIDAPlotProcessor(); 
	virtual ~LCFIAIDAPlotProcessor(); 
 
	virtual void init(); 
 
	virtual void processRunHeader( LCRunHeader* pRun ); 
 
	virtual void processEvent( LCEvent* pEvent ); 
 
	virtual void check( LCEvent* pEvent ); 
 
	virtual void end(); 
protected: 
	std::vector<std::string> _FlavourTagCollectionNames;
	std::vector<std::string> _FlavourTagInputsCollectionNames;
	std::string _TrueJetFlavourColName;
	std::string _TrueJetHadronChargeColName;
	std::string _TrueJetPDGCodeColName;
	std::string _TrueJetPartonChargeColName;
	std::string _JetCollectionName;
	std::string _VertexColName;
	std::string _CVertexChargeCollection;
	std::string _BVertexChargeCollection;

	//! cuts on all jets 
	double _CosThetaJetMax;
	//! cuts on all jets 
	double _CosThetaJetMin;
	//! cuts on all jets 
	double _PJetMin;
	//! cuts on all jets 
	double _PJetMax;

	//!Cut on the NN output variables - applied in vertex charge plots
	double _BTagNNCut;
	//!Cut on the NN output variables - applied in vertex charge plots
	double _CTagNNCut;

	//!optional parameters to make an ntuple of the neural net inputs; and print out the tagging ouputs (useful for scripts)
	bool _PrintNeuralNetOutput;
	bool _MakeTuple;
	std::string _NeuralNetOutputFile;
 
	int _iVertexChargeTagCollection;
	unsigned int _myVertexChargeTagCollection;

	std::vector<std::string> _VertexCatNames;
	std::vector<std::string>  _NumVertexCatDir;
	std::vector<std::string> _ZoomedVarNames;

	//!True B-jets - vertex charge vs true charge
	AIDA::IHistogram2D* _pBJetCharge2D;
	//!True C-jets - vertex charge vs true charge
	AIDA::IHistogram2D* _pCJetCharge2D;

	//!True B-jets - vertex charge leakage rate
	AIDA::IHistogram1D* _pBJetLeakageRate;
	//!True C-jets - vertex charge leakage rate
	AIDA::IHistogram1D* _pCJetLeakageRate;
	//!True B-jets - vertex charge
	AIDA::IHistogram1D* _pBJetVertexCharge;
	//!True C-jets - vertex charge
	AIDA::IHistogram1D* _pCJetVertexCharge;

	std::vector< std::map<std::string,unsigned int> > _IndexOfForEachTag;
	std::vector< std::map<std::string,unsigned int> > _InputsIndex;
	std::vector< std::map<std::string,unsigned int> > _ZoomedInputsIndex;

	//!Histograms of the neural net inputs for true B-jets
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _inputsHistogramsBJets;
	//!Histograms of the neural net inputs for true C-jets
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _inputsHistogramsCJets;
	//!Histograms of the neural net inputs for light B-jets
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _inputsHistogramsUDSJets;
	
	//!Zoomed-in histograms of some of the neural net inputs for true B-jets
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _zoomedInputsHistogramsBJets;
	//!Zoomed-in histograms of some of the neural net inputs for true C-jets
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _zoomedInputsHistogramsCJets;
	//!Zoomed-in histograms of some of the neural net inputs for true light-jets
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _zoomedInputsHistogramsUDSJets;

	//!Histograms of the neural net B-tag outputs for true light-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBTag;
	//!Histograms of the neural net C-tag outputs for true light-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetCTag;
	//!Histograms of the neural net B-tag outputs for true B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBTag;
	//!Histograms of the neural net C-tag outputs for true B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetCTag;
	//!Histograms of the neural net B-tag outputs for true C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBTag;
	//!Histograms of the neural net C-tag outputs for true C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetCTag; 
	//!Histograms of the neural net BC-tag outputs for true B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)   
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBCTag;
	//!Histograms of the neural net BC-tag outputs for true C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBCTag;
	//!Histograms of the neural net BC-tag outputs for true light-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBCTag;
	//!Histograms of the neural net B-tag outputs for non B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBTagBackgroundValues;
	//!Histograms of the neural net C-tag outputs for non C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCTagBackgroundValues;
	//!Histograms of the neural net BC-tag outputs for non C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBCTagBackgroundValues;
	
	//!Histograms of the neural net tags - number of events that pass a given cut: jet NN value > given NN value for the three tags - B-tag, C-tag, BC-tag
	//!  - separately for true B jets, true C jets & true light jets and different number of vertices in the jets, 1, 2 or >=3 & any (sum of previous three)
	//! See comments above
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBTagIntegral;    
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetCTagIntegral;  
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetCTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetCTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBCTagIntegral; 
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBCTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBCTagIntegral; 
	
	//!Number of bins used for neural nets plots
	int _numberOfPoints;

	//!number of different vertex categories we want to look at: 1 vertex, 2 vertices, >=3 vertices
	static const unsigned int N_VERTEX_CATEGORIES=3;  
			
	//!Tuple of the input variables - only filled for one inpur collection - selected with UseFlavourTagCollectionForVertexCharge
	AIDA::ITuple* _pMyTuple;

	
	int _lastRunHeaderProcessed;
	int _suppressOutputForRun;

	bool PassesEventCuts( LCEvent* pEvent );
	bool PassesJetCuts( ReconstructedParticle* pJet );
	void FillInputsPlots( LCEvent* pEvent, unsigned int jetNumber );
	void FillTagPlots( LCEvent* pEvent, unsigned int jetNumber );
	
	static const int B_JET=5;
	static const int C_JET=4;
	
	//!number of bins used in vertex charge leakage plots
	// this is something that could potentially become a user defined parameter
	static const int N_JETANGLE_BINS=10; 


	
	float CalculateDistance(const float* pos1, const float* pos2);

	//!Finds the true flavour of a jet (uses TrueJetFlavourCollection)
	int FindJetType( LCEvent* pEvent, unsigned int jetNumber );
	//!Finds the true charge of the hadron producing a jet (uses TrueJetHadronChargeCollection)
	float FindJetHadronCharge(LCEvent* pEvent, unsigned int jetNumber);
	//!Finds the PDG code of the hadron producing a jet (uses TrueJetPDGCodeCollection)
	int FindJetPDGCode( LCEvent* pEvent, unsigned int jetNumber );
	//!Finds the true charge of the parton producing a jet (uses TrueJetPartonChargeCollection)
	float FindJetPartonCharge(LCEvent* pEvent, unsigned int jetNumber);
	
	//!Finds the number of vertices in an event (from the flavour tag inputs)
	int FindNumVertex( LCEvent* pEvent, unsigned int jetNumber, unsigned int iInputsCollection);
	//!Finds the vertex charge of the jet - using cuts tuned to find vertex charge for C-jets (from CVertexChargeCollection)
	int FindCQVtx( LCEvent* pEvent, unsigned int jetNumber);
	//!Finds the vertex charge of the jet - using cuts tuned to find vertex charge for B-jets (from BVertexChargeCollection)
	int FindBQVtx( LCEvent* pEvent, unsigned int jetNumber);
	
	
	void PrintNNOutput();
	
	void InternalVectorInitialisation();
	
	//!Makes a DataPointSet of the tag efficiency e.g number of B-jets passing a given B-tag NN cut, as a function of NN
	AIDA::IDataPointSet* CreateEfficiencyPlot(const AIDA::IHistogram1D* pSignal, AIDA::IDataPointSet* pDataPointSet);
	//!Makes a DataPointSet integrating a histogram from the first bin to the last bin  -- NOT USED
	AIDA::IDataPointSet* CreateIntegralPlot(const AIDA::IHistogram1D* pNN, AIDA::IDataPointSet* pIntegral);
	//!Makes a DataPointSet of the tag purity e.g. N(B-jets passing NN cut)/N(all-jets passing NN cut) for a given B-tag NN cut, as a function of NN 
 	AIDA::IDataPointSet* CreatePurityPlot(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet);
	//!Makes a DataPointSet showing the tagging leakage e.g. the number of non-B-jets passing a given B-tag NN cut, as a function of NN 
	AIDA::IDataPointSet* CreateLeakageRatePlot(const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet);
	//!Plots two DataPointSets against each other
	AIDA::IDataPointSet* CreateXYPlot(const AIDA::IDataPointSet* pDataPointSet0, const AIDA::IDataPointSet* pDataPointSet1, AIDA::IDataPointSet* xyPointSet, const int dim0=0, const int dim1=0 );
	//!Makes a histogram integrating a histogram from the first bin to the last bin - THE ERRORS RETURNED ARE WRONG!
	AIDA::IHistogram1D* CreateIntegralHistogram(const AIDA::IHistogram1D* pNN, AIDA::IHistogram1D* pIntegral);

	//!Makes DataPointSets for the number of 
	void CreateVertexChargeLeakagePlot(AIDA::IDataPointSet* pBJetVtxChargeDPS, AIDA::IDataPointSet* pCJetVtxChargeDPS);

	//vertex plots
	AIDA::IHistogram1D* _pVertexDistanceFromIP;
	AIDA::IHistogram1D* _pVertexPositionX;
	AIDA::IHistogram1D* _pVertexPositionY;
	AIDA::IHistogram1D* _pVertexPositionZ;	
	AIDA::IHistogram1D* _pPrimaryVertexPullX;
	AIDA::IHistogram1D* _pPrimaryVertexPullY;
	AIDA::IHistogram1D* _pPrimaryVertexPullZ;
	AIDA::IHistogram1D* _pPrimaryVertexPositionX;
	AIDA::IHistogram1D* _pPrimaryVertexPositionY;
	AIDA::IHistogram1D* _pPrimaryVertexPositionZ;
	
	//!numbers of true C-jets with true charge ++
	int _cJet_truePlus2;
	//!numbers of true C-jets with true charge +
	int _cJet_truePlus;
	//!numbers of true C-jets with true charge 0
	int _cJet_trueNeut;
	//!numbers of true C-jets with true charge -
	int _cJet_trueMinus;
	//!numbers of true C-jets with true charge --
	int _cJet_trueMinus2;
	//!numbers of true C-jets with true charge ++; reconstructed vertex charge >0
	int _cJet_truePlus2_recoPlus; 
	//!numbers of true C-jets with true charge ++; reconstructed vertex charge =0	
	int _cJet_truePlus2_recoNeut;
	//!numbers of true C-jets with true charge ++; reconstructed vertex charge <0
	int _cJet_truePlus2_recoMinus;
	//!numbers of true C-jets with true charge +; reconstructed vertex charge >0
	int _cJet_truePlus_recoPlus; 
	//!numbers of true C-jets with true charge +; reconstructed vertex charge =0
	int _cJet_truePlus_recoNeut;
	//!numbers of true C-jets with true charge +; reconstructed vertex charge <0
	int _cJet_truePlus_recoMinus;
	//!numbers of true C-jets with true charge 0; reconstructed vertex charge >0
	int _cJet_trueNeut_recoPlus; 
	//!numbers of true C-jets with true charge 0; reconstructed vertex charge =0
	int _cJet_trueNeut_recoNeut;
	//!numbers of true C-jets with true charge 0; reconstructed vertex charge <0
	int _cJet_trueNeut_recoMinus;
	//!numbers of true C-jets with true charge -; reconstructed vertex charge >0
	int _cJet_trueMinus_recoPlus; 
	//!numbers of true C-jets with true charge -; reconstructed vertex charge =0
	int _cJet_trueMinus_recoNeut;
	//!numbers of true C-jets with true charge -; reconstructed vertex charge <0
	int _cJet_trueMinus_recoMinus;
	//!numbers of true C-jets with true charge --; reconstructed vertex charge >0
	int _cJet_trueMinus2_recoPlus; 
	//!numbers of true C-jets with true charge --; reconstructed vertex charge =0
	int _cJet_trueMinus2_recoNeut;
	//!numbers of true C-jets with true charge --; reconstructed vertex charge <0
	int _cJet_trueMinus2_recoMinus;

	//!numbers of true B-jets with true charge ++
	int _bJet_truePlus2;
	//!numbers of true B-jets with true charge +
	int _bJet_truePlus;	
	//!numbers of true B-jets with true charge 0
	int _bJet_trueNeut;	
	//!numbers of true B-jets with true charge -
	int _bJet_trueMinus;	
	//!numbers of true B-jets with true charge --
	int _bJet_trueMinus2;
	//!numbers of true B-jets with true charge ++; reconstructed vertex charge >0
	int _bJet_truePlus2_recoPlus; 
	//!numbers of true B-jets with true charge ++; reconstructed vertex charge =0
	int _bJet_truePlus2_recoNeut;
	//!numbers of true B-jets with true charge ++; reconstructed vertex charge <0
	int _bJet_truePlus2_recoMinus;
	//!numbers of true B-jets with true charge +; reconstructed vertex charge >0
	int _bJet_truePlus_recoPlus; 
	//!numbers of true B-jets with true charge +; reconstructed vertex charge =0
	int _bJet_truePlus_recoNeut;
	//!numbers of true B-jets with true charge +; reconstructed vertex charge <0
	int _bJet_truePlus_recoMinus;
	//!numbers of true B-jets with true charge 0; reconstructed vertex charge >0
	int _bJet_trueNeut_recoPlus; 
	//!numbers of true B-jets with true charge 0; reconstructed vertex charge =0
	int _bJet_trueNeut_recoNeut;
	//!numbers of true B-jets with true charge 0; reconstructed vertex charge <0
	int _bJet_trueNeut_recoMinus;
	//!numbers of true B-jets with true charge -; reconstructed vertex charge >0
	int _bJet_trueMinus_recoPlus; 
	//!numbers of true B-jets with true charge -; reconstructed vertex charge =0
	int _bJet_trueMinus_recoNeut;
	//!numbers of true B-jets with true charge -; reconstructed vertex charge <0
	int _bJet_trueMinus_recoMinus;
	//!numbers of true B-jets with true charge --; reconstructed vertex charge >0
	int _bJet_trueMinus2_recoPlus; 
	//!numbers of true B-jets with true charge --; reconstructed vertex charge =0
	int _bJet_trueMinus2_recoNeut;
	//!numbers of true B-jets with true charge --; reconstructed vertex charge <0
	int _bJet_trueMinus2_recoMinus;

	//! Vector of numbers of true C-jets with true charge ++
	//! See above for details
	std::vector< unsigned int>  _cJet_truePlus2_angle;
	std::vector< unsigned int>  _cJet_truePlus_angle;
	std::vector< unsigned int>  _cJet_trueNeut_angle;
	std::vector< unsigned int>  _cJet_trueMinus_angle;
	std::vector< unsigned int>  _cJet_trueMinus2_angle;
		     
	std::vector< unsigned int>  _cJet_truePlus2_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_truePlus2_recoNeut_angle;
	std::vector< unsigned int>  _cJet_truePlus2_recoMinus_angle;
	std::vector< unsigned int>  _cJet_truePlus_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_truePlus_recoNeut_angle;
	std::vector< unsigned int>  _cJet_truePlus_recoMinus_angle;
	std::vector< unsigned int>  _cJet_trueNeut_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_trueNeut_recoNeut_angle;
	std::vector< unsigned int>  _cJet_trueNeut_recoMinus_angle;
	std::vector< unsigned int>  _cJet_trueMinus_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_trueMinus_recoNeut_angle;
	std::vector< unsigned int>  _cJet_trueMinus_recoMinus_angle;
	std::vector< unsigned int>  _cJet_trueMinus2_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_trueMinus2_recoNeut_angle;
	std::vector< unsigned int>  _cJet_trueMinus2_recoMinus_angle;
		     
	std::vector< unsigned int>  _bJet_truePlus2_angle;
	std::vector< unsigned int>  _bJet_truePlus_angle;	
	std::vector< unsigned int>  _bJet_trueNeut_angle;	
	std::vector< unsigned int>  _bJet_trueMinus_angle;	
	std::vector< unsigned int>  _bJet_trueMinus2_angle;
	std::vector< unsigned int>  _bJet_truePlus2_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_truePlus2_recoNeut_angle;
	std::vector< unsigned int>  _bJet_truePlus2_recoMinus_angle;
	std::vector< unsigned int>  _bJet_truePlus_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_truePlus_recoNeut_angle;
	std::vector< unsigned int>  _bJet_truePlus_recoMinus_angle;
	std::vector< unsigned int>  _bJet_trueNeut_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_trueNeut_recoNeut_angle;
	std::vector< unsigned int>  _bJet_trueNeut_recoMinus_angle;
	std::vector< unsigned int>  _bJet_trueMinus_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_trueMinus_recoNeut_angle;
	std::vector< unsigned int>  _bJet_trueMinus_recoMinus_angle;
	std::vector< unsigned int>  _bJet_trueMinus2_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_trueMinus2_recoNeut_angle;
	std::vector< unsigned int>  _bJet_trueMinus2_recoMinus_angle;
	 
}; 


#endif // endif for "ifndef LCFIAIDAPlotProcessor_h"
