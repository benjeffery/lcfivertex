#ifndef NeuralNetTrainer_h
#define NeuralNetTrainer_h

//Marlin and LCIO includes
#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"

//Neural Net includes
#include "nnet/inc/NeuralNet.h"
#include "nnet/inc/NeuralNetDataSet.h"
#include "nnet/inc/BackPropagationCGAlgorithm.h"

/** Trains neural networks to be used for jet flavour tagging.
 *
 * Trains flavour tagging networks using the BackPropagationCGAlgorithm (see the \ref NeuralNet page) with
 * 500 epochs. The networks are trained on the following data:<br>
 * If only 1 vertex is found (i.e. only the interaction point)
 * \f[tanh\left( \frac{D0Significance1}{100} \right)\f]
 * \f[tanh\left( \frac{D0Significance2}{100} \right)\f]
 * \f[tanh\left( \frac{Z0Significance1}{100} \right)\f]
 * \f[tanh\left( \frac{Z0Significance2}{100} \right)\f]
 * \f[JointProbRPhi\f]
 * \f[JointProbZ\f]
 * \f[tanh\left( \frac{3\times Momentum1}{E} \right)\f]
 * \f[tanh\left( \frac{3\times Momentum2}{E} \right)\f]
 *
 * If 2 or more vertices are found
 * \f[tanh\left( \frac{DecayLengthSignificance}{6\times E} \right)\f]
 * \f[tanh\left( \frac{DecayLength}{10} \right)\f]
 * \f[tanh\left( \frac{PTCorrectedMass}{5} \right)\f]
 * \f[tanh\left( \frac{RawMomentum}{E} \right)\f]
 * \f[JointProbRPhi\f]
 * \f[JointProbZ\f]
 * \f[tanh\left( \frac{NumTracksInVertices}{10} \right)\f]
 * \f[SecondaryVertexProbability\f]
 *
 * Where E is the jet energy, everything else is the data calculated by FlavourTagInputsProcessor.<br>
 * Note that <b>the processor applies it's own hard coded cuts</b>. These are documented under the full
 * description of NeuralNetTrainerProcessor::_passesCuts().
 *
 * <H4>Input</H4>
 * - A collection of ReconstructedParticles that represents the jets in the event (put in by your jet
 * finder, say SatoruJetFinderProcessor).
 * - A LCFloatVec collection that holds the true jet flavours, in the same order as the jets (put in
 * by TrueAngularJetFlavourProcessor).
 * - A LCFloatVec collection that holds the flavour tag inputs (put in by FlavourTagInputsProcessor)
 * - Between 1 and 9 filenames for the generated networks. If an output filename is left blank then
 * that network is not trained, but if none are supplied then a lcio::Exception is thrown.
 *
 * <H4>Output</H4>
 * Trained neural networks to the filenames supplied, in the format requested. The LCIO file is not
 * modified at all.
 *
 * @param JetCollectionName Name of the ReconstructedParticle collection that represents jets.
 * @param FlavourTagInputsCollection Name of the LCFloatVec collection that holds the flavour tag inputs.
 * @param TrueJetFlavourCollection Name of the LCIntVec Collection that contains the true jet flavours.
 * @param Filename-b_net-1vtx Output filename for the trained 1 vertex b-tag net.
 * @param Filename-c_net-1vtx Output filename for the trained 1 vertex c-tag net.
 * @param Filename-bc_net-1vtx Output filename for the trained 1 vertex c-tag (with only b background) net.
 * @param Filename-b_net-2vtx Output filename for the trained 2 vertex b-tag net.
 * @param Filename-c_net-2vtx Output filename for the trained 2 vertex c-tag net.
 * @param Filename-bc_net-2vtx Output filename for the trained 2 vertex c-tag (with only b background) net.
 * @param Filename-b_net-3plusvtx Output filename for the trained 3 or more vertices b-tag net.
 * @param Filename-c_net-3plusvtx Output filename for the trained 3 or more vertices c-tag net.
 * @param Filename-bc_net-3plusvtx Output filename for the trained  3 or more vertices c-tag (with only b background) net.
 *
 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
*/

class NeuralNetTrainerProcessor : public marlin::Processor
{
public:
	//The usual Marlin processor methods
	virtual Processor* newProcessor() { return new NeuralNetTrainerProcessor; }
	NeuralNetTrainerProcessor();
	virtual ~NeuralNetTrainerProcessor();
	virtual void init();
	virtual void processRunHeader( LCRunHeader* pRun );
	virtual void processEvent( LCEvent* pEvent );
	//don't need this
	//virtual void check( LCEvent* pEvent );
	virtual void end();
protected:
	//variables for the steering file options
	std::string _JetCollectionName;	//The name of the collection of ReconstructedParticles that is the jet (comes from the steering file)
	std::string _FlavourTagInputsCollectionName;
	std::string _TrueJetFlavourCollectionName;
	int _serialiseAsXML;
	nnet::NeuralNet::SerialisationMode _outputFormat;

	//These maps all use the same string keys to distinguish between the different nets.
	//The strings are of the form "c_net-2vtx", "bc_net-3vtx" etcetera.
	std::map<std::string,std::string> _filename;/**< @internal A map of the output filenames for the trained net. All of these maps use strings of the form "c_net-2vtx", "bc_net-3vtx" etcetera as the key*/
	std::map<std::string,bool> _trainThisNet;	/**< @internal Map containing true or false depending on whether this net has been selected for
							training (if the user supplies an output filename in the steering file).*/
	std::map<std::string,nnet::NeuralNetDataSet*> _dataSet;/**< @internal Map of the data sets for each of the selected nets */
	std::map<std::string,int> _numSignal;  /**< @internal Map of the number of signal events for each net. Not really used for anything, just printed
							as info when the net starts training.*/
	std::map<std::string,int> _numBackground;  /**< @internal Dito for the number of backgrounds.*/

	//List of strings for all the nets enabled for training
	//Strings are the same as the keys used in the maps above
	std::vector<std::string> _listOfSelectedNetNames; /**< @internal A list of the nets that have been selected for training, in the form of the strings
								used in the map keys above.*/

	//This map holds the position of the Inputs in the LCFloatVec
	std::map<std::string,unsigned int> _IndexOf; /**< @internal Holds the positions of the inputs, so that, for example, you can get the position of
							D0Significance1 in the inputs LCFloatVec just with _IndexOf["D0Significance1"]*/
	
	int _nRun; /**< @internal The run number.*/
	int _nEvent;/**< @internal The event number.*/
	int _nAcceptedEvents; /**< @internal The number of events that have passed the cuts so far.*/
	
	//useful constants
	static const int C_JET=4; /**< @internal Just a useful constant for testing true jet flavour.*/
	static const int B_JET=5; /**< @internal Ditto.*/

	//The following functions are just code that has been split off so that the code doesn't look quite so cluttered.
	void _displayCollectionNames( lcio::LCEvent* pEvent );/**< @internal Displays all of the available collections in the file.*/
	void _trainNet( nnet::BackPropagationCGAlgorithm& pBackPropCGAlgo, nnet::NeuralNetDataSet& dataSet );/**< @internal The training code, split off to make the code a bit more manageable*/
	bool _passesCuts( lcio::LCEvent* pEvent );///< @internal All the code for the cuts should be put in here; returns false if the event fails any of the cuts.
};

#endif //ifndef NeuralNetTrainer_h
