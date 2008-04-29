#ifndef FlavourTag_h
#define FlavourTag_h

//Marlin and LCIO includes
#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"

//Neural Net includes
#include "nnet/inc/NeuralNet.h"
#include "nnet/inc/NeuralNetDataSet.h"
#include "nnet/inc/BackPropagationCGAlgorithm.h"



/** Performs a neural net based flavour tag using data calculated by the LCFI vertex package.
*
* Loads in previously trained neural networks from the filenames provided in the steering file, and
* performs a flavour tag with them using the data previously stored in the file by the
* FlavourTagInputsProcessor. The networks can be trained using the NeuralNetTrainerProcessor.<br>
* This processor requires 9 neural networks, which are 3 for each of the 1 vertex, 2 vertices and 3 or
* more vertices cases. These 3 are a b jet tagging network, a c jet tagging network and a c jet with
* only b background tagging network. If any of these saved neural network files are not present the
* processor will throw a lcio::Exception. The networks can be in either text or XML format; the
* processor checks to see if the file starts with "<?xml" and decides whether to load as text or XML.<br>
* N.B. The code that loads the XML networks is currently a little shaky. <b>If the XML is not properly
* formed then you may get a segmentation fault or runaway memory allocation leading to Marlin crashing.</b>
* This is still being looked into.<br>
* For more information on the tagging variables used as input, have a look at the documentation for 
* FlavourTagInputsProcessor. The flavour tag result will be in the range 0 to 1; so to
* select tagged jets apply a cut on this value (e.g. the b-tag value to tag b-jets).  If anything goes
* wrong (that doesn't produce an exception) then a -1 will be stored instead.
*
* <H4>Input</H4>
* - 9 previously trained neural networks (trained by NeuralNetTrainerProcessor).<br>
* - A collection of LCFloatVec that hold the flavour tag variables (put in the lcio file by
* FlavourTagInputsProcessor).<br>
* - A collection of ReconstructedParticles that represents the jets in the event (put in by your jet
* finder, say SatoruJetFinderProcessor).<br>
*
* <H4>Output</H4>
* - A collection of LCFloatVec that contains the 3 tag results for each jet (b tag, c tag and c only b
* background tag). The collection will have a float vector for each jet in the same order as the jets; so
* for example, the tags for "pJetCollection->getElementAt(3)" will be in "pFlavourTagCollection->getElementAt(3)".
* For more details see \ref LCIO "the interface documentation"<br>
*
* @param JetCollectionName The name of the collection of ReconstructedParticles representing the jets.
* @param FlavourTagInputsCollection The name of the collection of LCFloatVec that is the flavour tag inputs.
* @param FlavourTagCollection The name of the collection of the flavour tag results that will be created.
* @param Filename-b_net-1vtx Filename for the 1 vertex b tag network.
* @param Filename-c_net-1vtx Filename for the 1 vertex c tag network.
* @param Filename-bc_net-1vtx Filename for the 1 vertex c tag (only b background) network.
* @param Filename-b_net-2vtx Filename for the 2 vertex b tag network.
* @param Filename-c_net-2vtx Filename for the 2 vertex c tag network.
* @param Filename-bc_net-2vtx Filename for the 2 vertex c (only b background) tag network.
* @param Filename-b_net-3plusvtx Filename for the 3 or more vertices b tag network.
* @param Filename-c_net-3plusvtx Filename for the 3 or more vertices c tag network.
* @param Filename-bc_net-3plusvtx Filename for the 3 or more vertices c tag (only b background) network.
*
* @author Mark Grimes (mark.grimes@bristol.ac.uk)
*/
class FlavourTagProcessor : public marlin::Processor
{
public:
	//The usual Marlin processor methods
	virtual Processor* newProcessor() { return new FlavourTagProcessor; }
	FlavourTagProcessor();
	virtual ~FlavourTagProcessor();
	virtual void init();
	virtual void processRunHeader( LCRunHeader* pRun );
	virtual void processEvent( LCEvent* pEvent );
	//don't need this
	//virtual void check( LCEvent* pEvent );
	virtual void end();
protected:
	std::string _JetCollectionName;	//The name of the collection of ReconstructedParticles that is the jet (comes from the steering file)
	std::string _FlavourTagInputsCollectionName;
	std::string _FlavourTagCollectionName;
	int _nRun;
	int _evt;
	std::map<std::string,std::string> _filename;//The input filenames for the nets.
	std::map<std::string,nnet::NeuralNet*> _NeuralNet;//Pointers to the neural nets
	//ofstream ofile;
	//This map holds the position of the Inputs in the LCFloatVec
	std::map<std::string,unsigned int> _IndexOf;
	
	void _displayCollectionNames( lcio::LCEvent* pEvent );
	
};

#endif //ifndef NeuralNetTrainer_h
