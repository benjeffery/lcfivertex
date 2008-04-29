
#include "PlotProcessor.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>
#ifdef USEROOT
/////////////////////
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
/////////////////////
#endif // of #ifdef USEROOT

#include "EVENT/LCCollection.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "EVENT/LCParameters.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/LCIntVec.h"

#include "util/inc/vector3.h"
#include "util/inc/util.h"

using std::includes;
using std::map;
using std::string;
using std::endl;
using std::cout;
using std::vector;
using std::stringstream;
using std::abs;
using std::pair;
using std::ofstream;
using std::cerr;
using std::set;
using namespace lcio;

//Needs to be instantiated for Marlin to know about it (I think)
PlotProcessor aPlotProcessor;

PlotProcessor::PlotProcessor() : marlin::Processor("Plot")
{
	_description = "Plots various outputs from the flavour tag" ;

	// register steering parameters: name, description, class-variable, default value
	//The name of the collection of ReconstructedParticles that is the jet
	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" , 
				"Name of the collection of ReconstructedParticles that is the jet"  ,
				_JetCollectionName ,
				string("SGVJets") ) ;
	vector<string> FlavourTagCollectionNamesDefault;
	FlavourTagCollectionNamesDefault.push_back("FlavourTag");
	registerProcessorParameter("FlavourTagCollections" , 
			      "Name of the LCFloatVec Collections that contain the flavour tags (one purity efficiency plot per tag) (in same order as jet collection)"  ,
			      _FlavourTagCollectionNames,
			      FlavourTagCollectionNamesDefault) ;
       registerInputCollection( LCIO::LCINTVEC,
  			      "TrueJetFlavourCollection" , 
			      "Name of the output collection of LCIntVec (same order as jets)"  ,
			      _TrueJetFlavourColName ,
			      std::string("TrueJetFlavour") ) ;
	//The output filename
	registerProcessorParameter( "OutputFilename" , 
				"Filename for the output"  ,
				_OutputFilename ,
				string("PlotProcessorOutput") ) ;
}

PlotProcessor::~PlotProcessor()
{
}

void PlotProcessor::init()
{
	printParameters();
	cout << _description << endl
		<< "-------------------------------------------------" << endl
		<< endl;
		
	_nRun=0;
	_jetEMax=0;
	//Make as many binning containers as we need
	for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
	{
		_BTagEfficiencyPurity.push_back(efficiency_purity<double>());
		_CTagEfficiencyPurity.push_back(efficiency_purity<double>());
		_BCTagEfficiencyPurity.push_back(efficiency_purity<double>());
	}
}

void PlotProcessor::processRunHeader( LCRunHeader* pRun )
{
	//
	// Perform a check to see if the variable names we need are here
	//
	for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
	{
		vector<string> VarNames;
		(pRun->parameters()).getStringVals(_FlavourTagCollectionNames[iTag],VarNames);
		//Fill the map realting names and indexes
		set<string> AvailableNames;
		map<string,unsigned int> IndexOf;
		for (size_t i = 0;i < VarNames.size();++i)
		{
			AvailableNames.insert(VarNames[i]);
			IndexOf[VarNames[i]] = i;
		}
		
		//Add the index to the list
		_IndexOfForEachTag.push_back(IndexOf);
		
		//Check the required information is in the LCFloatVec
		set<string> RequiredNames;
		RequiredNames.insert("BTag");
		RequiredNames.insert("CTag");
		RequiredNames.insert("BCTag");
		
		if (!includes(AvailableNames.begin(),AvailableNames.end(),RequiredNames.begin(),RequiredNames.end()))
			cerr << _FlavourTagCollectionNames[iTag] << " does not contain information required by PlotProcessor";
	}
}

void PlotProcessor::processEvent( LCEvent* pEvent )
{
	//Get the collection of jets. Can't do anything if the collection isn't there
	//so don't bother catching the exception and terminate.
	LCCollection* pJetCollection=pEvent->getCollection( _JetCollectionName );
	
	//make sure the collection is of the right type
	if( pJetCollection->getTypeName()!=LCIO::RECONSTRUCTEDPARTICLE )
	{
		stringstream message;
		message << endl
			<< "########################################################################################\n"
			<< "# PlotProcessor -                                                                      #\n"
			<< "#   The jet collection requested (\"" << _JetCollectionName << "\") is not of the type \"" << LCIO::RECONSTRUCTEDPARTICLE << "\"  #\n"
			<< "########################################################################################" << endl;
		throw EventException( message.str() );
	}
	
	//apply any cuts on the event here
	if( _passesEventCuts(pEvent) )
	{
		ReconstructedParticle* pJet;
		//loop over the jets
		for( int a=0; a<pJetCollection->getNumberOfElements(); ++a )
		{
			//Dynamic casts are not the best programming practice in the world, but I can't see another way of doing this
			//in the LCIO framework.  This cast should be safe though because we've already tested the type.
			pJet=dynamic_cast<ReconstructedParticle*>( pJetCollection->getElementAt(a) );
			
			if( _passesJetCuts(pJet) )
			{
				_fillPlots(pEvent, a );
				_jetEnergy.add(pJet->getEnergy());
				if (pJet->getEnergy() > _jetEMax) _jetEMax = pJet->getEnergy();
			}
		}
	}
}


void PlotProcessor::end()
{
	_outputDataToFile( _OutputFilename );

	cout << "The largest jet energy was " << _jetEMax << endl;
}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently applies no cuts at all*/
bool PlotProcessor::_passesEventCuts( LCEvent* pEvent )
{
	//
	// No event cuts at present
	//

	return true;
}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently selects jets for which the jet polar angle theta is -0.95< cos(theta) <0.95. 
*/
bool PlotProcessor::_passesJetCuts( ReconstructedParticle* pJet )
{
	//
	// This cut added on the suggestion of Sonja Hillert 12/Jan/07.
	//
	// Selects jets for which the cosine of the jet polar
	// angle theta for all jets is not too large.
	//
	// Make something that's easy to search for to track down erroneous cuts:
	// GREPTAG_CUT : Jet cut on abs(cos(theta)) of jet axis
	//
	
	double CThJ_lower=0;	//lower cut on cos(theta) of the jet axis
	double CThJ_upper=0.95;	//upper cut (Sho Ryu Ken!)

	vertex_lcfi::util::Vector3 zAxis( 0, 0, 1 ); //work out theta from dot product with jet axis

	const double* mom=pJet->getMomentum();
	vertex_lcfi::util::Vector3 jetMomentum( mom[0], mom[1], mom[2] );

	jetMomentum.makeUnit();
	double cosTheta=jetMomentum.dot( zAxis );
	if( fabs(cosTheta)<=CThJ_lower || fabs(cosTheta)>=CThJ_upper ) return false;


	// If control gets to this point then the jet has passed
	return true;
}

void PlotProcessor::_fillPlots( LCEvent* pEvent, unsigned int jet)
{
	LCCollection* pTrueCollection=pEvent->getCollection( _TrueJetFlavourColName );
	int jetType = dynamic_cast<LCIntVec*>(pTrueCollection->getElementAt(jet))->back();
	//LCCollection* pTrueCollection=pEvent->getCollection( "SGVFlavourTagInputs" );
	//int jetType=static_cast<int>((*dynamic_cast<LCFloatVec*>(pTrueCollection->getElementAt(jet)))[19]);
	for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
	{
		LCCollection* pTagCollection=pEvent->getCollection( _FlavourTagCollectionNames[iTag] );
		double bTag= (*dynamic_cast<LCFloatVec*>(pTagCollection->getElementAt(jet)))[_IndexOfForEachTag[iTag]["BTag"]];
		double cTag= (*dynamic_cast<LCFloatVec*>(pTagCollection->getElementAt(jet)))[_IndexOfForEachTag[iTag]["CTag"]];
		double cTagBBack= (*dynamic_cast<LCFloatVec*>(pTagCollection->getElementAt(jet)))[_IndexOfForEachTag[iTag]["BCTag"]];
	
		if( jetType==B_JET )
		{
			if( bTag<=1 && bTag>=0 ) _BTagEfficiencyPurity[iTag].add_signal( bTag );
			if( cTag<=1 && cTag>=0 ) _CTagEfficiencyPurity[iTag].add_background( cTag );
			if( cTagBBack<=1 && cTagBBack>=0 ) _BCTagEfficiencyPurity[iTag].add_background( cTagBBack );
		}
		else if( jetType==C_JET )
		{
			if( bTag<=1 && bTag>=0 ) _BTagEfficiencyPurity[iTag].add_background( bTag );
			if( cTag<=1 && cTag>=0 ) _CTagEfficiencyPurity[iTag].add_signal( cTag );
			if( cTagBBack<=1 && cTagBBack>=0 ) _BCTagEfficiencyPurity[iTag].add_signal( cTagBBack );
		}
		else
		{
			if( bTag<=1 && bTag>=0 ) _BTagEfficiencyPurity[iTag].add_background( bTag );
			if( cTag<=1 && cTag>=0 ) _CTagEfficiencyPurity[iTag].add_background( cTag );
			//don't add background for the bcnet result because we're only considering the b background
		}
	}
}

void PlotProcessor::_outputDataToFile( string filename )
{
#ifdef USEROOT
	stringstream filenameStream;
	filenameStream << filename << ".root";
	TFile rootFile( filenameStream.str().c_str(), "RECREATE", "Various plots from the plot processor" );

	//have to put all the results into normal c arrays for root to deal with
	//there may be an easier way of doing this conversion but I don't know what it is
	
	for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
	{
		int numberOfPoints=100;
	
		//get the results in STL form
		vector< pair<double,double> > BTagefficiencyPurityPairs=_BTagEfficiencyPurity[iTag].eff_pur( numberOfPoints );
		vector< pair<double,double> > CTagefficiencyPurityPairs=_CTagEfficiencyPurity[iTag].eff_pur( numberOfPoints );
		vector< pair<double,double> > BCTagefficiencyPurityPairs=_BCTagEfficiencyPurity[iTag].eff_pur( numberOfPoints );
	
		//have to put all the results into normal c arrays for root to deal with
		//there may be an easier way of doing this conversion but I don't know what it is
		double BNetEfficiency[100];
		double BNetPurity[100];
		double CNetEfficiency[100];
		double CNetPurity[100];
		double BCNetEfficiency[100];
		double BCNetPurity[100];
		
		for( int a=0; a<numberOfPoints; ++a )
		{
			BNetEfficiency[a]=BTagefficiencyPurityPairs[a].first;
			BNetPurity[a]=BTagefficiencyPurityPairs[a].second;
			CNetEfficiency[a]=CTagefficiencyPurityPairs[a].first;
			CNetPurity[a]=CTagefficiencyPurityPairs[a].second;
			BCNetEfficiency[a]=BCTagefficiencyPurityPairs[a].first;
			BCNetPurity[a]=BCTagefficiencyPurityPairs[a].second;
		}
	
		TGraph BNetGraph( 99, BNetEfficiency, BNetPurity );
		BNetGraph.SetName( (_FlavourTagCollectionNames[iTag] + " B Tag").c_str() );
		TGraph CNetGraph( 99, CNetEfficiency, CNetPurity );
		CNetGraph.SetName( (_FlavourTagCollectionNames[iTag] + " C Tag").c_str() );
		TGraph BCNetGraph( 99, BCNetEfficiency, BCNetPurity );
		BCNetGraph.SetName( (_FlavourTagCollectionNames[iTag] + " BC Tag").c_str() );
	
		TMultiGraph TagGraphs;
		TagGraphs.SetName( (_FlavourTagCollectionNames[iTag]).c_str() );
		TagGraphs.Add( &BNetGraph );
		TagGraphs.Add( &CNetGraph );
		TagGraphs.Add( &BCNetGraph );
	
		TAxis* pAxis;
		TagGraphs.SetTitle( (string("Efficiency-purity for ") + _FlavourTagCollectionNames[iTag]).c_str() );
	
		//why the hell doesn't this work?
		pAxis=TagGraphs.GetXaxis();
		if(pAxis) pAxis->SetTitle( "Efficiency" );
	
		pAxis=TagGraphs.GetYaxis();
		if(pAxis) pAxis->SetTitle( "Purity" );
	
		rootFile.Add( &TagGraphs );
		rootFile.Write();
	}
	// now add in the jet energies
	TH1F jetEnergyHistogram( "jetEnergyHistogram", "Jet energies", 200, 0, 120 );
	const vector<double> jetEnergyData=_jetEnergy.sorted_data();
	for( vector<double>::const_iterator i=jetEnergyData.begin(); i<jetEnergyData.end(); ++i ) jetEnergyHistogram.Fill( (*i) );
	
	rootFile.Write();
	rootFile.Close();

#else	//of ifdef USEROOT

	//
	// Not using root so output everything to comma seperated value files (csv)
	//
	ofstream ofile;
	ofile.open( (filename + ".csv").c_str() );
	if( ofile.is_open() )
	{
		for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
		{
			string n = _FlavourTagCollectionNames[iTag];
			ofile << "E," << n << " B Tag,E," << n << " C Tag,E," << n << " BC Tag,";
		}
		ofile << endl;

		vector<vector<pair<double,double> > > BTagefficiencyPurityPairs;
		vector<vector<pair<double,double> > > CTagefficiencyPurityPairs;
		vector<vector<pair<double,double> > > BCTagefficiencyPurityPairs;

		int numberOfPoints=100;
		for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
		{
			BTagefficiencyPurityPairs.push_back(_BTagEfficiencyPurity[iTag].eff_pur( numberOfPoints ));
			CTagefficiencyPurityPairs.push_back(_CTagEfficiencyPurity[iTag].eff_pur( numberOfPoints ));
			BCTagefficiencyPurityPairs.push_back(_BCTagEfficiencyPurity[iTag].eff_pur( numberOfPoints ));
	//		std::cout << _BTagEfficiencyPurity[iTag].number_of_signal() << " " << _BTagEfficiencyPurity[iTag].number_of_background() << std::endl;
	//		std::cout << _CTagEfficiencyPurity[iTag].number_of_signal() << " " << _CTagEfficiencyPurity[iTag].number_of_background() << std::endl;
	//		std::cout << _BCTagEfficiencyPurity[iTag].number_of_signal() << " " << _BCTagEfficiencyPurity[iTag].number_of_background() << std::endl;
		}
		for( int a=0; a<numberOfPoints; ++a )
		{
			for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag)
			{
				ofile << BTagefficiencyPurityPairs[iTag][a].first << "," << BTagefficiencyPurityPairs[iTag][a].second << ",";
				ofile << CTagefficiencyPurityPairs[iTag][a].first << "," << CTagefficiencyPurityPairs[iTag][a].second << ",";
				ofile << BCTagefficiencyPurityPairs[iTag][a].first << "," << BCTagefficiencyPurityPairs[iTag][a].second << ",";
			}
				ofile << endl;
		}
		ofile.close();
	}
	else 
	{
		cerr << "########################################################################################\n"
			<< "# PlotProcessor -                                                                      #\n"
			<< "#   Unable to open file \"" << filename << "\" for output   #\n"
			<< "########################################################################################" << endl;
	}

	stringstream filenameStream;
	filenameStream << filename << "-JetEnergies.csv";

	ofile.open( filenameStream.str().c_str() );
	if( ofile.is_open() )
	{
		vector< vertex_lcfi::util::bin<double> > binned=_jetEnergy.binned_data( 200, 0, 120 );
		ofile << "Bin low,Bin High,Frequency" << endl;
		for( vector< vertex_lcfi::util::bin<double> >::iterator i=binned.begin(); i<binned.end(); ++i )
			ofile << (*i).region_low() << "," << (*i).region_high() << "," << (*i).contents() << endl;
		ofile.close();
	}
	else
	{
		cerr << "########################################################################################\n"
			<< "# PlotProcessor -                                                                      #\n"
			<< "#   Unable to open file \"" << filenameStream.str() << "\" for output   #\n"
			<< "########################################################################################" << endl;
	}
#endif //of ifdef USEROOT
}

void PlotProcessor::_displayCollectionNames( LCEvent* pEvent )
{
	const vector<string>* pCollectionNames=pEvent->getCollectionNames();
	
	cout << "The available collections are: (name - type)" << endl;
	for( vector<string>::const_iterator i=pCollectionNames->begin(); i<pCollectionNames->end(); ++i )
	{
		LCCollection* pCollection=pEvent->getCollection( (*i) );
		const string typeName=pCollection->getTypeName();
		cout << "  " << (*i) << " - " << typeName << endl;
	}
	cout << endl;
}
