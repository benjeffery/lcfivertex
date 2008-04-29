#include "../include/FlavourTag.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>

#include "EVENT/LCCollection.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ParticleIDImpl.h"
#include "EVENT/LCParameters.h"
#include "EVENT/LCFloatVec.h"
#include "util/inc/memorymanager.h"
#include "util/inc/vector3.h"

#include "nnet/inc/NeuralNet.h"

using std::set;
using std::string;

FlavourTagProcessor aFlavourTagProcessor;

FlavourTagProcessor::FlavourTagProcessor() : marlin::Processor("FlavourTag")
{
	_description = "Performs a flavour tag using previously trained neural nets" ;

	// register steering parameters: name, description, class-variable, default value
	//The name of the collection of ReconstructedParticles that is the jet
	registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" , 
				"Name of the collection of ReconstructedParticles that is the jet"  ,
				_JetCollectionName ,
				std::string("FTSelectedJets") ) ;
        registerInputCollection( lcio::LCIO::LCFLOATVEC,
  			      "FlavourTagInputsCollection" , 
			      "Name of the LCFloatVec Collection that contains the flavour tag inputs (in same order as jet collection)"  ,
			      _FlavourTagInputsCollectionName,
			      "FlavourTagInputs" ) ;
        registerOutputCollection( lcio::LCIO::LCFLOATVEC,
  			      "FlavourTagCollection" , 
			      "Name of the LCFloatVec Collection that will be created to contain the flavour tag result"  ,
			      _FlavourTagCollectionName,
			      "FlavourTag" ) ;
	//These are the variables for the input filenames of the previously trained nets
	registerProcessorParameter( "Filename-b_net-1vtx" , 
				"Filename of the previously trained 1 vertex b-tag net"  ,
				_filename["b_net-1vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-c_net-1vtx" , 
				"Filename of the previously trained 1 vertex c-tag net"  ,
				_filename["c_net-1vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-bc_net-1vtx" , 
				"Filename of the previously trained 1 vertex c-tag (b background only) net"  ,
				_filename["bc_net-1vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-b_net-2vtx" , 
				"Filename of the previously trained 2 vertex b-tag net"  ,
				_filename["b_net-2vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-c_net-2vtx" , 
				"Filename of the previously trained 2 vertex c-tag net"  ,
				_filename["c_net-2vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-bc_net-2vtx" , 
				"Filename of the previously trained 2 vertex c-tag (b background only) net"  ,
				_filename["bc_net-2vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-b_net-3plusvtx" , 
				"Filename of the previously trained 3 (or more) vertex b-tag net"  ,
				_filename["b_net-3vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-c_net-3plusvtx" , 
				"Filename of the previously trained 3 (or more) vertex c-tag net"  ,
				_filename["c_net-3vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-bc_net-3plusvtx" , 
				"Filename of the previously trained 3 (or more) vertex c-tag (b background only) net"  ,
				_filename["bc_net-3vtx"],
				std::string("") ) ;
}

FlavourTagProcessor::~FlavourTagProcessor()
{
}

void FlavourTagProcessor::init()
{
	//ofile.open("FTInputs.txt", ofstream::out);
		
	
	printParameters();
	std::cout << _description << std::endl
		<< "-------------------------------------------------" << std::endl
		<< std::endl;
		
	_nRun=0;
	_evt=0;
	//Run through the filenames given, try and open them and also see if they appear to be XML or plain text,
	//because the neural net code causes a segmentation fault if you try and open one type as another.
	//Remember "(*i).second" is the filename and "(*i).first" is the string key that identifies the net.
	for( std::map<std::string,std::string>::iterator iPair=_filename.begin(); iPair!=_filename.end(); ++iPair )
	{
		std::ifstream inputFile( (*iPair).second.c_str() );
		if( inputFile.is_open() )
		{
			nnet::NeuralNet::SerialisationMode fileFormat;//whether the file is XML or PlainText
			
			//read the first line
			std::string firstLine;//the first line (up until the first space) read from the file
			inputFile >> firstLine;//check to see if this is how we expect an XML file to start

			if ( firstLine=="<?xml" ) fileFormat=nnet::NeuralNet::XML;
			else fileFormat=nnet::NeuralNet::PlainText;

			inputFile.close();

			//Print what we're trying to do so the user knows what's happened if this goes wrong
			std::cout << "FlavourTag: Attempting to load the " << (*iPair).first << " network as ";
			if( fileFormat==nnet::NeuralNet::XML ) std::cout << "XML";
			else std::cout << "plain text";
			std::cout << " from file " << (*iPair).second << " ..." << std::endl;

			//N.B. If fileFormat is wrong could get a segmentation fault!
			_NeuralNet[ (*iPair).first ]=new nnet::NeuralNet( (*iPair).second, fileFormat );
			vertex_lcfi::MemoryManager<nnet::NeuralNet>::Run()->registerObject( _NeuralNet[ (*iPair).first ] );
		}
		else
		{
			std::stringstream message;
			message << std::endl
				<< "########################################################################################\n"
				<< "# FlavourTagProcessor -                                                                #\n"
				<< "#   Unable to open file " << (*iPair).second << " for the " << (*iPair).first << " neural net.    #\n"
				<< "########################################################################################" << std::endl;
			throw lcio::Exception( message.str() );
		}
	}

	// If control gets to here then all the files loaded okay. Tell the user in case something goes wrong later and this gets blamed.
	std::cout << "FlavourTag: All networks loaded okay." << std::endl;
}

void FlavourTagProcessor::processRunHeader( LCRunHeader* pRun )
{
	++_evt;
	//Get the current list of variable names
	std::vector<std::string> VarNames;
	(pRun->parameters()).getStringVals(_FlavourTagInputsCollectionName,VarNames);
	
	set<string> AvailableNames;
	//Fill the map realting names and indexes
	for (size_t i = 0;i < VarNames.size();++i)
	{
		AvailableNames.insert(VarNames[i]);
		_IndexOf[VarNames[i]] = i;
	}
	
	//Check the required information is in the LCFloatVec
	set<string> RequiredNames;
	RequiredNames.insert( "D0Significance1");
	RequiredNames.insert( "D0Significance2");
	RequiredNames.insert( "Z0Significance1" );
	RequiredNames.insert( "Z0Significance2" );
	RequiredNames.insert( "JointProbRPhi");
	RequiredNames.insert( "JointProbZ");
	RequiredNames.insert( "Momentum1");
	RequiredNames.insert( "Momentum2");	
	RequiredNames.insert( "DecayLengthSignificance");
	RequiredNames.insert( "DecayLength");
	RequiredNames.insert( "PTCorrectedMass");
	RequiredNames.insert( "RawMomentum");
	RequiredNames.insert( "NumTracksInVertices" );
	RequiredNames.insert( "SecondaryVertexProbability");
	
	if (!includes(AvailableNames.begin(),AvailableNames.end(),RequiredNames.begin(),RequiredNames.end()))
		std::cerr << _FlavourTagInputsCollectionName << " does not contain information required by FlavourTagProcessor";
	
	VarNames.clear();
	//Set the names of the output collection
	VarNames.push_back( "BTag" );
	VarNames.push_back( "CTag" );
	VarNames.push_back( "BCTag" );
	pRun->parameters().setValues(_FlavourTagCollectionName, VarNames);
	
}

void FlavourTagProcessor::processEvent( lcio::LCEvent* pEvent )
{
	//Get the collection of jets. Can't do anything if the collection isn't there
	//so don't bother catching the exception and terminate.
	lcio::LCCollection* pJetCollection=pEvent->getCollection( _JetCollectionName );
	
	//make sure the collection is of the right type
	if( pJetCollection->getTypeName()!=lcio::LCIO::RECONSTRUCTEDPARTICLE )
	{
		std::stringstream message;
		message << std::endl
			<< "########################################################################################\n"
			<< "# FlavourTagProcessor -                                                                #\n"
			<< "#   The jet collection requested (\"" << _JetCollectionName << "\") is not of the type \"" << lcio::LCIO::RECONSTRUCTEDPARTICLE << "\"  #\n"
			<< "########################################################################################" << std::endl;
		throw lcio::EventException( message.str() );
	}

	//Create the collection to store the result
	LCCollectionVec* OutCollection = new LCCollectionVec("LCFloatVec");
	pEvent->addCollection(OutCollection,_FlavourTagCollectionName);
	
	//loop over the jets
	for( int a=0; a<pJetCollection->getNumberOfElements(); ++a )
	{
		lcio::ReconstructedParticle* pJet;
		//Dynamic casts are not the best programming practice in the world, but I can't see another way of doing this
		//in the LCIO framework.  This cast should be safe though because we've already tested the type.
		pJet=dynamic_cast<lcio::ReconstructedParticle*>( pJetCollection->getElementAt(a) );
		
		// Find out the jet energy to work out the correct normalisation constants
		double jetEnergy=pJet->getEnergy();
		if( 0==jetEnergy )
		{
			jetEnergy=45.5;
			if( isFirstEvent() ) std::cerr << "*** FlavourTag - Warning: Jet energy undefined, assuming 45.5GeV ***" << std::cout;
		}
	
		// Variables for the normalisation of the inputs
		double Norm_D0Significance		= 100.0;
		double Norm_Z0Significance		= 100.0;
		double Norm_Momentum			= jetEnergy/3.0;
		double Norm_DecayLengthSignificance	= 6.0*jetEnergy;
		double Norm_DecayLength			= 1.0;
		double Norm_PTMassCorrection		= 5.0;
		double Norm_RawMomentum			= jetEnergy;
		double Norm_NumTracksInVertices		= 10.0;
	
		//
		// See if we can get the required info from the file
		//
		
		lcio::LCCollection* pInputs=pEvent->getCollection( _FlavourTagInputsCollectionName );
		
		//make sure the collection is of the right type
		if( pInputs->getTypeName()!=lcio::LCIO::LCFLOATVEC )
		{
			std::stringstream message;
			message << std::endl
				<< "########################################################################################\n"
				<< "# FlavourTagProcessor -                                                                #\n"
				<< "#   The jet collection requested (\"" << _FlavourTagInputsCollectionName << "\") is not of the type \"" << lcio::LCIO::LCFLOATVEC << "\"  #\n"
				<< "########################################################################################" << std::endl;
			throw lcio::EventException( message.str() );
		}
		
		LCFloatVec* FTInputs = dynamic_cast<lcio::LCFloatVec*>( pInputs->getElementAt(a) );
		double NumVertices = (*FTInputs)[_IndexOf["NumVertices"]];
		
		std::vector<double> inputs;
		if( NumVertices==1 )
		{
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["D0Significance1"]]/Norm_D0Significance) );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["D0Significance2"]]/Norm_D0Significance) );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["Z0Significance1"]]/Norm_Z0Significance) );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["Z0Significance2"]]/Norm_Z0Significance) );
			inputs.push_back( (*FTInputs)[_IndexOf["JointProbRPhi"]] );
			inputs.push_back( (*FTInputs)[_IndexOf["JointProbZ"]] );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["Momentum1"]]/Norm_Momentum) );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["Momentum2"]]/Norm_Momentum) );
			
		}
		else
		{
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["DecayLengthSignificance"]]/Norm_DecayLengthSignificance) );
			inputs.push_back( std::tanh(((*FTInputs)[_IndexOf["DecayLength"]]/10.0)/Norm_DecayLength));
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["PTCorrectedMass"]]/Norm_PTMassCorrection) );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["RawMomentum"]]/Norm_RawMomentum) );
			inputs.push_back( (*FTInputs)[_IndexOf["JointProbRPhi"]] );
			inputs.push_back( (*FTInputs)[_IndexOf["JointProbZ"]] );
			inputs.push_back( std::tanh((*FTInputs)[_IndexOf["NumTracksInVertices"]]/Norm_NumTracksInVertices) );
			inputs.push_back( (*FTInputs)[_IndexOf["SecondaryVertexProbability"]] );
		}
			
	 	// Perform the tag.
		
		std::vector<double> bTagOutput;
		std::vector<double> cTagOutput;
		std::vector<double> cTagbBackgroundOutput;
	
		if( NumVertices==1 )
		{
			bTagOutput=_NeuralNet["b_net-1vtx"]->output( inputs );
			cTagOutput=_NeuralNet["c_net-1vtx"]->output( inputs );
			cTagbBackgroundOutput=_NeuralNet["bc_net-1vtx"]->output( inputs );
		}
		else if( NumVertices==2 )
		{
			bTagOutput=_NeuralNet["b_net-2vtx"]->output( inputs );
			cTagOutput=_NeuralNet["c_net-2vtx"]->output( inputs );
			cTagbBackgroundOutput=_NeuralNet["bc_net-2vtx"]->output( inputs );
		}
		else if( NumVertices>=3 )
		{
			bTagOutput=_NeuralNet["b_net-3vtx"]->output( inputs );
			cTagOutput=_NeuralNet["c_net-3vtx"]->output( inputs );
			cTagbBackgroundOutput=_NeuralNet["bc_net-3vtx"]->output( inputs );
		}
		else
		{
			//Don't know why this happens, I had thought ZVTop always returned at least the IP...
			
			std::cerr << "FlavourTagProcessor - Warning: Unexpected multiplicity of " << NumVertices << "!" << std::endl;
			
			// The output vector sizes will be zero, so the next few lines will add in the invalid
			// output value to keep the particle ID parameters the same size as expected from the
			// names in the run header.
		}
	
		//
		// Now store the data in the file
		//
		LCFloatVec* OutVec = new LCFloatVec();
		
		if( bTagOutput.size()==1 ) // If NumVertices is 0, as sometimes happens, then this will be empty.
		{
			OutVec->push_back(bTagOutput.back());
		}
		else
		{
			// Something has gone wrong, but some data needs to be in the file otherwise it could
			// screw up the other parameters by putting them in a different order than expected.
			std::cerr << "FlavourTagProcessor - Warning: b tag output has size " << bTagOutput.size()
					<< ". Putting invalid output value of -1 in the LCIO file."<< std::endl;
			OutVec->push_back(-1);
		}
		
		if( cTagOutput.size()==1 )
		{
			OutVec->push_back(cTagOutput.back());
		}
		else
		{
			std::cerr << "FlavourTagProcessor - Warning: c tag output has size " << cTagOutput.size()
					<< ". Putting invalid output value of -1 in the LCIO file."<< std::endl;
			OutVec->push_back(-1);
		}
	
		if( cTagbBackgroundOutput.size()==1 )
		{
			OutVec->push_back(cTagbBackgroundOutput.back());
		}
		else
		{
			std::cerr << "FlavourTagProcessor - Warning: c tag (b background only) output has size " << cTagbBackgroundOutput.size()
					<< ". Putting invalid output value of -1 in the LCIO file."<< std::endl;
			OutVec->push_back(-1);
		}
		
		OutCollection->addElement(OutVec);
	//ofile << (*OutVec)[0] << "\t" << (*OutVec)[1] << "\t" << (*OutVec)[2] << std::endl;
	//ofile << "----------------------" << std::endl;
	}
	//Clear anything that may have been allocated during this event
	vertex_lcfi::MetaMemoryManager::Event()->delAllObjects();
}

void FlavourTagProcessor::end()
{
	//ofile.close();
	//free up stuff
	vertex_lcfi::MetaMemoryManager::Run()->delAllObjects();
}

void FlavourTagProcessor::_displayCollectionNames( lcio::LCEvent* pEvent )
{
	const std::vector<std::string>* pCollectionNames=pEvent->getCollectionNames();
	
	std::cout << "The available collections are: (name - type)" << std::endl;
	for( std::vector<std::string>::const_iterator i=pCollectionNames->begin(); i<pCollectionNames->end(); ++i )
	{
		lcio::LCCollection* pCollection=pEvent->getCollection( (*i) );
		const std::string typeName=pCollection->getTypeName();
		std::cout << "  " << (*i) << " - " << typeName << std::endl;
	}
	std::cout << std::endl;
}
