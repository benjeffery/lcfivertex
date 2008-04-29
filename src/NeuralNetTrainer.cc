#include "NeuralNetTrainer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>


#include "EVENT/LCCollection.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/LCParameters.h"

#include "util/inc/memorymanager.h"
#include "util/inc/vector3.h"

#include "nnet/inc/NeuralNet.h"
#include "nnet/inc/NeuralNetDataSet.h"
#include "nnet/inc/SigmoidNeuronBuilder.h"
#include "nnet/inc/BackPropagationCGAlgorithm.h"

//Needs to be instantiated for Marlin to know about it (I think)
NeuralNetTrainerProcessor aNeuralNetTrainerProcessor;

NeuralNetTrainerProcessor::NeuralNetTrainerProcessor() : marlin::Processor("NeuralNetTrainer")
{
	_description = "Trains a neural net from the lcio file" ;

	// register steering parameters: name, description, class-variable, default value

	//The name of the collection of ReconstructedParticles that is the jet
	registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" , 
				"Name of the collection of ReconstructedParticles that is the jet"  ,
				_JetCollectionName ,
				std::string("SGVJets") ) ;
	registerInputCollection( lcio::LCIO::LCFLOATVEC,
  			      "FlavourTagInputsCollection" , 
			      "Name of the LCFloatVec Collection that contains the flavour tag inputs (in same order as jet collection)"  ,
			      _FlavourTagInputsCollectionName,
			      "FlavourTagInputs" ) ;
        registerInputCollection( lcio::LCIO::LCINTVEC,
  			      "TrueJetFlavourCollection" , 
			      "Name of the LCIntVec Collection that contains the true jet flavours"  ,
			      _TrueJetFlavourCollectionName,
			      "TrueJetFlavour" ) ;

	registerProcessorParameter( "SaveAsXML",
				"Set this to 0 to output the neural nets in plain text format (default), or 1 (or anything non zero) to save in XML format",
				_serialiseAsXML,
				0 );

	//These are the variables for the output filenames of the trained nets
	//Default is "" which will switch off training for that net.
	registerProcessorParameter( "Filename-b_net-1vtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["b_net-1vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-c_net-1vtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["c_net-1vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-bc_net-1vtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["bc_net-1vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-b_net-2vtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["b_net-2vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-c_net-2vtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["c_net-2vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-bc_net-2vtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["bc_net-2vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-b_net-3plusvtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["b_net-3vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-c_net-3plusvtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["c_net-3vtx"],
				std::string("") ) ;
	registerProcessorParameter( "Filename-bc_net-3plusvtx" , 
				"Output filename for the trained net. If it is blank (default) then this net is not trained"  ,
				_filename["bc_net-3vtx"],
				std::string("") ) ;
}

NeuralNetTrainerProcessor::~NeuralNetTrainerProcessor()
{
}

void NeuralNetTrainerProcessor::init()
{
	printParameters();
	std::cout << _description << std::endl
		<< "-------------------------------------------------" << std::endl
		<< std::endl;
		
	_nRun=0;

	//Have a look through all the net output filenames supplied in the steering file.
	//If any of them are blank (the default if one isn't supplied) disable training for
	//that net.  Here "(*i).second" is the filename and "(*i).first" is the string key that
	//identifies each net in all of the maps.
	for( std::map<std::string,std::string>::iterator i=_filename.begin(); i!=_filename.end(); ++i )
	{
		if( (*i).second!="" )
		{
			_listOfSelectedNetNames.push_back( (*i).first );//make a list of the selected map names to make looping over them easier later on
			_trainThisNet[ (*i).first ]=true;//turn on training of this net
		}
		else _trainThisNet[ (*i).first ]=false;//turn off training of this net
	}
	
	//Just check that the user hasn't accidently disabled training of all the nets
	if( _listOfSelectedNetNames.size()==0 )
	{
		std::stringstream message;
		message << std::endl
			<< "############################################################################################\n"
			<< "# NeuralNetTrainerProcessor:                                                               #\n"
			<< "#   No nets have been enabled for training in the steering file!                           #\n"
			<< "#   Supply at least one output filename in the steering file. For example, put the line    #\n"
			<< "#      <parameter name=\"Filename-bc_net-1vtx\" type=\"string\"> bc_net-1vtx.xml </parameter>  #\n"
			<< "#   In with the parameters for this processor.                                             #\n"
			<< "############################################################################################" << std::endl;
		throw lcio::Exception( message.str() );
	}

	//Decide which format (plain text or XML) to save the files as
	if( _serialiseAsXML==0 ) _outputFormat=nnet::NeuralNet::PlainText;
	else _outputFormat=nnet::NeuralNet::XML;

	//Allocate a data set for each of the enabled nets.
	for( std::vector<std::string>::iterator iName=_listOfSelectedNetNames.begin(); iName<_listOfSelectedNetNames.end(); ++iName )
	{
		_dataSet[*iName]=new nnet::NeuralNetDataSet;
		vertex_lcfi::MemoryManager<nnet::NeuralNetDataSet>::Run()->registerObject( _dataSet[*iName] );
		
		//also set all of the signal/background counters to 0
		_numSignal[*iName]=0;
		_numBackground[*iName]=0;
	}

	// Initialise the event counters
	_nEvent=0;
	_nAcceptedEvents=0;
}

void NeuralNetTrainerProcessor::processRunHeader( LCRunHeader* pRun )
{
	_nRun++;
	
	//Get the list of flavour tag inputs Available
	std::vector<std::string> VarNames;
	(pRun->parameters()).getStringVals(_FlavourTagInputsCollectionName,VarNames);
	
	std::set<std::string> AvailableNames;
	for (size_t i = 0;i < VarNames.size();++i)
	{
		AvailableNames.insert(VarNames[i]);
		_IndexOf[VarNames[i]] = i;
	}
	
	//Check the required information is in the LCFloatVec
	std::set<std::string> RequiredNames;
	RequiredNames.insert( "D0Significance1");
	RequiredNames.insert( "D0Significance2");
	RequiredNames.insert( "Z0Significance1" );
	RequiredNames.insert( "D0Significance2" );
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
		std::cerr << _FlavourTagInputsCollectionName << " does not contain information required by NeuralNetTrainerProcessor";
	
}

void NeuralNetTrainerProcessor::processEvent( lcio::LCEvent* pEvent )
{
	//Output the collection names for debugging
	if( isFirstEvent() ) _displayCollectionNames( pEvent );

	//Get the collection of jets. Can't do anything if the collection isn't there
	//so don't bother catching the exception and terminate.
	lcio::LCCollection* pJetCollection=pEvent->getCollection( _JetCollectionName );
	
	//make sure the collection is of the right type
	if( pJetCollection->getTypeName()!=lcio::LCIO::RECONSTRUCTEDPARTICLE )
	{
		std::stringstream message;
		message << std::endl
			<< "########################################################################################\n"
			<< "# NeuralNetTrainerProcessor:                                                           #\n"
			<< "#   The jet collection requested (\"" << _JetCollectionName << "\") is not of the type \"" << lcio::LCIO::RECONSTRUCTEDPARTICLE << "\"  #\n"
			<< "########################################################################################" << std::endl;
		throw lcio::EventException( message.str() );
	}

	lcio::ReconstructedParticle* pJet;
	int numJets=pJetCollection->getNumberOfElements();
	
	//apply any cuts on the event here
	if( _passesCuts(pEvent) )
	{
		//loop over the jets
		for( int a=0; a<numJets; ++a )
		{
			//Dynamic casts are not the best programming practice in the world, but I can't see another way of doing this
			//in the LCIO framework.  This cast should be safe though because we've already tested the type.
			pJet=dynamic_cast<lcio::ReconstructedParticle*>( pJetCollection->getElementAt(a) );

			// Find out the jet energy to work out the correct normalisation constants
			double jetEnergy=pJet->getEnergy();
			if( 0==jetEnergy )
			{
				jetEnergy=45.5;
				if( isFirstEvent() ) std::cout << "*** NeuralNetTrainer - Warning: Jet energy undefined, assuming 45.5GeV ***" << std::cout;
			}

/*
	-----------------------------------
	-------------IMPORTANT-------------
	-----------------------------------
	If any of these normalisation constants are changed, update in the documentation by modifying
	the main class description in NeuralNetTrainer.h (at around line 20). 
*/
			// Variables for the normalisation of the inputs
			double Norm_D0Significance		= 100.0;
			double Norm_Z0Significance		= 100.0;
			double Norm_Momentum			= jetEnergy/3.0;
			double Norm_DecayLengthSignificance	= 6.0*jetEnergy;
			double Norm_DecayLength			= 1.0;
			double Norm_PTMassCorrection		= 5.0;
			double Norm_RawMomentum			= jetEnergy;
			double Norm_NumTracksInVertices		= 10.0;
		
			//Get the MC Jet type
			lcio::LCCollection* pTrueJet=pEvent->getCollection( _TrueJetFlavourCollectionName );
			//make sure the collection is of the right type
			if( pTrueJet->getTypeName()!=lcio::LCIO::LCINTVEC )
			{
				std::stringstream message;
				message << std::endl
					<< "########################################################################################\n"
					<< "# FlavourTagProcessor -                                                                #\n"
					<< "#   The jet collection requested (\"" << _TrueJetFlavourCollectionName << "\") is not of the type \"" << lcio::LCIO::LCINTVEC << "\"  #\n"
					<< "########################################################################################" << std::endl;
				throw lcio::EventException( message.str() );
			}
			int jetType = *((dynamic_cast<lcio::LCIntVec*>( pTrueJet->getElementAt(a))->begin()));
			
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
			LCFloatVec Inputs = *(dynamic_cast<lcio::LCFloatVec*>( pInputs->getElementAt(a) ));
			
			std::vector<double> inputs;
			std::vector<double> target;
			
			double NumVertices = Inputs[_IndexOf["NumVertices"]];
			//TODO Check that the inputs exist in the index
			if( NumVertices==1 )
			{
				inputs.push_back( std::tanh(Inputs[_IndexOf["D0Significance1"]]/Norm_D0Significance) );
				inputs.push_back( std::tanh(Inputs[_IndexOf["D0Significance2"]]/Norm_D0Significance) );
				inputs.push_back( std::tanh(Inputs[_IndexOf["Z0Significance1"]]/Norm_Z0Significance) );
				inputs.push_back( std::tanh(Inputs[_IndexOf["D0Significance2"]]/Norm_Z0Significance) );
				inputs.push_back( Inputs[_IndexOf["JointProbRPhi"]] );
				inputs.push_back( Inputs[_IndexOf["JointProbZ"]] );
				inputs.push_back( std::tanh(Inputs[_IndexOf["Momentum1"]]/Norm_Momentum) );
				inputs.push_back( std::tanh(Inputs[_IndexOf["Momentum2"]]/Norm_Momentum) );
			}
			else
			{
				inputs.push_back( std::tanh(Inputs[_IndexOf["DecayLengthSignificance"]]/Norm_DecayLengthSignificance) );
				inputs.push_back( std::tanh((Inputs[_IndexOf["DecayLength"]]/10.0)/Norm_DecayLength));
				inputs.push_back( std::tanh(Inputs[_IndexOf["PTCorrectedMass"]]/Norm_PTMassCorrection) );
				inputs.push_back( std::tanh(Inputs[_IndexOf["RawMomentum"]]/Norm_RawMomentum) );
				inputs.push_back( Inputs[_IndexOf["JointProbRPhi"]] );
				inputs.push_back( Inputs[_IndexOf["JointProbZ"]] );
				inputs.push_back( std::tanh(Inputs[_IndexOf["NumTracksInVertices"]]/Norm_NumTracksInVertices) );
				inputs.push_back( Inputs[_IndexOf["SecondaryVertexProbability"]] );
			}
		
			if( jetType==B_JET )
			{
				target.clear();
				target.push_back( 1.0 );
				if( _trainThisNet["b_net-1vtx"] && NumVertices==1 ){ _dataSet["b_net-1vtx"]->addDataItem( inputs, target );_numSignal["b_net-1vtx"]+=1;}
				if( _trainThisNet["b_net-2vtx"] && NumVertices==2 ){ _dataSet["b_net-2vtx"]->addDataItem( inputs, target );_numSignal["b_net-2vtx"]+=1;}
				if( _trainThisNet["b_net-3vtx"] && NumVertices>=3 ){ _dataSet["b_net-3vtx"]->addDataItem( inputs, target );_numSignal["b_net-3vtx"]+=1;}
				target.clear();
				target.push_back( 0.0 );
				if( _trainThisNet["c_net-1vtx"] && NumVertices==1 ){ _dataSet["c_net-1vtx"]->addDataItem( inputs, target );_numBackground["c_net-1vtx"]+=1;}
				if( _trainThisNet["c_net-2vtx"] && NumVertices==2 ){ _dataSet["c_net-2vtx"]->addDataItem( inputs, target );_numBackground["c_net-2vtx"]+=1;}
				if( _trainThisNet["c_net-3vtx"] && NumVertices>=3 ){ _dataSet["c_net-3vtx"]->addDataItem( inputs, target );_numBackground["c_net-3vtx"]+=1;}
				if( _trainThisNet["bc_net-1vtx"] && NumVertices==1 ){ _dataSet["bc_net-1vtx"]->addDataItem( inputs, target );_numBackground["bc_net-1vtx"]+=1;}
				if( _trainThisNet["bc_net-2vtx"] && NumVertices==2 ){ _dataSet["bc_net-2vtx"]->addDataItem( inputs, target );_numBackground["bc_net-2vtx"]+=1;}
				if( _trainThisNet["bc_net-3vtx"] && NumVertices>=3 ){ _dataSet["bc_net-3vtx"]->addDataItem( inputs, target );_numBackground["bc_net-3vtx"]+=1;}
			}
			else if( jetType==C_JET )
			{
				target.clear();
				target.push_back( 0.0 );
				if( _trainThisNet["b_net-1vtx"] && NumVertices==1 ){ _dataSet["b_net-1vtx"]->addDataItem( inputs, target );_numBackground["b_net-1vtx"]+=1;}
				if( _trainThisNet["b_net-2vtx"] && NumVertices==2 ){ _dataSet["b_net-2vtx"]->addDataItem( inputs, target );_numBackground["b_net-2vtx"]+=1;}
				if( _trainThisNet["b_net-3vtx"] && NumVertices>=3 ){ _dataSet["b_net-3vtx"]->addDataItem( inputs, target );_numBackground["b_net-3vtx"]+=1;}
				target.clear();
				target.push_back( 1.0 );
				if( _trainThisNet["c_net-1vtx"] && NumVertices==1 ){ _dataSet["c_net-1vtx"]->addDataItem( inputs, target );_numSignal["c_net-1vtx"]+=1;}
				if( _trainThisNet["c_net-2vtx"] && NumVertices==2 ){ _dataSet["c_net-2vtx"]->addDataItem( inputs, target );_numSignal["c_net-2vtx"]+=1;}
				if( _trainThisNet["c_net-3vtx"] && NumVertices>=3 ){ _dataSet["c_net-3vtx"]->addDataItem( inputs, target );_numSignal["c_net-3vtx"]+=1;}
				if( _trainThisNet["bc_net-1vtx"] && NumVertices==1 ){ _dataSet["bc_net-1vtx"]->addDataItem( inputs, target );_numSignal["bc_net-1vtx"]+=1;}
				if( _trainThisNet["bc_net-2vtx"] && NumVertices==2 ){ _dataSet["bc_net-2vtx"]->addDataItem( inputs, target );_numSignal["bc_net-2vtx"]+=1;}
				if( _trainThisNet["bc_net-3vtx"] && NumVertices>=3 ){ _dataSet["bc_net-3vtx"]->addDataItem( inputs, target );_numSignal["bc_net-3vtx"]+=1;}
			}
			else
			{
				target.clear();
				target.push_back( 0.0 );
				if( _trainThisNet["b_net-1vtx"] && NumVertices==1 ){ _dataSet["b_net-1vtx"]->addDataItem( inputs, target );_numBackground["b_net-1vtx"]+=1;}
				if( _trainThisNet["b_net-2vtx"] && NumVertices==2 ){ _dataSet["b_net-2vtx"]->addDataItem( inputs, target );_numBackground["b_net-2vtx"]+=1;}
				if( _trainThisNet["b_net-3vtx"] && NumVertices>=3 ){ _dataSet["b_net-3vtx"]->addDataItem( inputs, target );_numBackground["b_net-3vtx"]+=1;}
				if( _trainThisNet["c_net-1vtx"] && NumVertices==1 ){ _dataSet["c_net-1vtx"]->addDataItem( inputs, target );_numBackground["c_net-1vtx"]+=1;}
				if( _trainThisNet["c_net-2vtx"] && NumVertices==2 ){ _dataSet["c_net-2vtx"]->addDataItem( inputs, target );_numBackground["c_net-2vtx"]+=1;}
				if( _trainThisNet["c_net-3vtx"] && NumVertices>=3 ){ _dataSet["c_net-3vtx"]->addDataItem( inputs, target );_numBackground["c_net-3vtx"]+=1;}
				//don't fill anything for the bc net because this isn't a b or a c jet
			}

		}

		++_nAcceptedEvents;
	}
	++_nEvent;

	//Clear anything that may have been allocated during this event
	vertex_lcfi::MetaMemoryManager::Event()->delAllObjects();
}

/*
-----------------------------------
-------------IMPORTANT-------------
-----------------------------------
If you change the cuts make sure you change the line below to show the changes in the docs*/
/*! Currently selects jets for which the jet polar angle theta is -0.866<= cos(theta) <=0.866. 
*/
bool NeuralNetTrainerProcessor::_passesCuts( lcio::LCEvent* pEvent )
{
	//Any cuts would go in here.  Do a test and then return false if the event fails, let control carry on to the other
	//tests if it passes.  If the event has passed all the cuts then there is a return true at the end.
	//Currently only a cut on the jet momentum theta.

	std::vector<vertex_lcfi::util::Vector3> jetMomentums;

	//Don't want to flood the screen if the data is not available, so count how many times a warning has been given.
	static int numWarningsNoMomentum=0;//The number of times a warning has been printed that momentum data is not available.

	try
	{
		lcio::LCCollection* pJetCollection=pEvent->getCollection( _JetCollectionName );
		for( int i=0; i<pJetCollection->getNumberOfElements(); ++i )
		{
			const double* mom=(dynamic_cast<lcio::ReconstructedParticle*>( pJetCollection->getElementAt(i) ))->getMomentum();
			if( mom[0]==0 && mom[1]==0 && mom[2]==0 ) throw lcio::Exception( "Jet momentum not defined" );

			jetMomentums.push_back( vertex_lcfi::util::Vector3( mom[0], mom[1], mom[2] ) );
		}
	}
	catch( lcio::Exception exception )
	{
		//Just print a warning and proceed with the other cuts.
		if( numWarningsNoMomentum<=2 )
		{
			std::cerr << "############################################################################\n"
				<< "#   NeuralNetTrainerProcessor:                                             #\n"
				<< "#      Unable to get the data for the jet momentum because -               #\n"
				<< "#                                                                          #\n"
				<< exception.what() << std::endl
				<< "#                                                                          #\n"
				<< "#      Training will proceed with no cut on the jet momentum theta.        #\n";
			if ( numWarningsNoMomentum==2 ) std::cerr << "#   NO FURTHER WARNINGS WILL BE DISPLAYED                                  #\n";
			std::cerr << "############################################################################" << std::endl;
			++numWarningsNoMomentum;
		}
	}

	//
	// Fail this event if any of the jets have 30deg < theta < 150deg
	//
	
	vertex_lcfi::util::Vector3 zAxis( 0, 0, 1 );

	//If the above try block failed then this will be empty and control will move on to the other cuts.
	for( std::vector<vertex_lcfi::util::Vector3>::iterator iMom=jetMomentums.begin(); iMom<jetMomentums.end(); ++iMom )
	{
		(*iMom).makeUnit();
		double cosTheta=(*iMom).dot( zAxis );
		//equivalent to "if( theta>150degrees || theta<30degrees )"
		//Should maybe have this value as a steering file parameter?
		if( cosTheta>0.866 || cosTheta<-0.866 ) return false;
	}


	//other cuts would go in here, returning false if the event fails


	//If control gets to here then the event has passed all the cuts.
	return true;
}

void NeuralNetTrainerProcessor::end()
{
	//
	//The data sets should all be filled, so train the nets
	//

	//Just make one neuron builder and reuse it
	nnet::SigmoidNeuronBuilder neuronBuilder;

	//All the nets use the same amount of nodes, so just create one of these and reuse it
	int nInputs=8;
	std::vector<int> nodes;
	nodes.push_back(2 * nInputs - 2);
	nodes.push_back(1);

	//print out info on how many events passed the cuts
	std::cout << "NeuralNetTrainer: " << _nAcceptedEvents << " of " << _nEvent << " events passed the cuts. See the documentation of NeuralNetTrainerProcessor::_passesCuts() for details of the cuts applied." << std::endl;

	//Train and save any nets that have been selected for training
	for( std::vector<std::string>::iterator iName=_listOfSelectedNetNames.begin(); iName<_listOfSelectedNetNames.end(); ++iName )
	{
		//Not going to need these once the net is trained and saved so have them local to this 'for' loop
		nnet::NeuralNet thisNeuralNet( nInputs, nodes, &neuronBuilder, 1 );
		nnet::BackPropagationCGAlgorithm myAlgorithm( thisNeuralNet );

		std::cout << std::endl << "Training neural net " << *iName << " with " << _dataSet[*iName]->numberOfDataItems()
				<< " jets " << "(" << _numBackground[*iName] << " background, " << _numSignal[*iName] << " signal)..." << std::endl;

		// Make sure we can open the file before training.  Nothing worse than waiting ages to train and then losing the result!
		std::ofstream outputFile( _filename[*iName].c_str() );
		if( outputFile.is_open() )
		{
			//do the training
			_trainNet( myAlgorithm, *_dataSet[*iName] );

			//Set the output format to the one requested in the steering file
			thisNeuralNet.setSerialisationMode( _outputFormat );
			thisNeuralNet.serialise( outputFile );
			outputFile.close();
		}
		else
		{
			std::cerr << "Unable to open file " << _filename[*iName] << "! Skipping training for this net." << std::endl;
		}
	}
	

	std::cout << "Finished training all selected nets" << std::endl;
	
	//free up stuff
	vertex_lcfi::MetaMemoryManager::Run()->delAllObjects();
}

void NeuralNetTrainerProcessor::_trainNet( nnet::BackPropagationCGAlgorithm& backPropCGAlgo, nnet::NeuralNetDataSet& dataSet )
{
	//This function pretty much just calls backPropCGAlgo.train(...) at the moment, although code can easily be added
	//to check the errors after each iteration 
	double PrevErr,CurrErr;//Training errors
	int i=0;
	bool breakLoop=false;//not actually used at the moment, but will be used to cut the loop early if required

	while( i<50 && breakLoop==false )
	{
		// for CG algorithm with 10 epochs per loop index
		//A bit silly calling 10 epochs 50 times instead of just 500 epochs, but I'll
		//hopefully put in some code to cut out early depending on how the training is
		//going.
		backPropCGAlgo.train( 10, dataSet );

		//Have a look at the errors.  Not a lot of point at the moment but you could put
		//something in here to cut out early if the errors aren't getting significantly
		//smaller (e.g. "if( (CurrErr-PrevErr)/PrevErr < 0.02 ) breakLoop=true")
		std::vector<double> epochErrors=backPropCGAlgo.getTrainingErrorValuesPerEpoch();
		CurrErr=epochErrors.back();
		PrevErr = CurrErr;

		// 26/Apr/07 - Been having problems with the net not training and getting a NAN
		//error under certain conditions which are still being looked into. This takes
		//ages, and the loop keeps trying to do the same thing over and over, so check
		//for this and cut the loop to save running time. Maybe dump the data set somewhere?
		if( std::isnan( CurrErr ) )
		{
			std::cerr 	<< "NeuralNetTrainer.cc (line 523): Training the net gave an error of NaN! That's not good. Still looking\n"
					<< "into why this happens, most likely there's not enough difference between the tag variables of your\n"
					<< "signal and background data sets. Your net is going to be gibberish. Sorry." << std::endl;
			breakLoop=true;
		}
		
		++i;
	}
}

void NeuralNetTrainerProcessor::_displayCollectionNames( lcio::LCEvent* pEvent )
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
