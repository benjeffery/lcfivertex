#ifndef PlotProcessor_h
#define PlotProcessor_h

#include <vector>
#include <algorithm>

#include "../vertex_lcfi/util/inc/util.h"

//Marlin and LCIO includes
#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"

using vertex_lcfi::util::efficiency_purity;
using vertex_lcfi::util::histogram_data;

/** Creates some sample plots from the data calculated by the LCFI vertex package.
 *
 * An example of getting the flavour tag results from the LCIO file and plotting an efficiency
 * purity graph with them.  Also plots a graph of jet energies for good measure.<BR>
 * The processor checks the specifed LCFloatVec collections for the flavour tag values "BTag", "CTag" and
 * "BCTag" which are the names that FlavourTagProcessor stores its b tag, c tag and c tag
 * (only b background) values in respectively. <BR>
 * These values are checked against the true jet flavour (from the TrueJetFlavour LCIntVec)
 * and efficiency-purity values calculated for a range of cuts.<BR>
 * The jet energy is taken from the energy of the reconstructed particle used to represent the jet.
 *
 * <H5>Getting Root output</H5>
 * To output to a Root file instead of CSV files the processor has to be compiled with the USEROOT preprocessor
 * flag defined. You could add "#define USEROOT" to the code, or more easily add the line
 *
 * USERINCLUDES += -D USEROOT
 *
 * to the userlib.gmk file that is in the Marlin directory. If Marlin is not already set up to use
 * Root then you will also need to add the following lines (this assumes a fully working root installation):
 *
 * USERINCLUDES += `root-config --cflags`<br>
 * USERLIBS += `root-config --libs`
 *
 * <H4>Input</H4>
 * From the LCIO file, flavour tag variable values of:
 * 
 * "BTag"     "CTag"      "BCTag"
 *
 * And
 *
 * "JetType"
 *
 * <H4>Output</H4>
 * If the USEROOT preprocessor flag was defined when this processor was compiled, then the output
 * will be a root file with the filename specified in the steering file.
 * Otherwise, the efficiency-purity values will be output as comma separated values to the file
 * <filename>+".csv", and the jet energies to <filename>+"-JetEnergies.csv".
 *
 * @param JetCollectionName Name of the ReconstructedParticle collection that represents jets.
 * @param FlavourTagCollections Names of the LCFloatVec collections holding the Flavour tags, all tags
 * in this list will be produced in one file for comparison
 * @param TrueJetFlavourCollection LCIntVec that contains the MC Jet flavour (from TrueJetFlavourProcessor)
 * @param OutputFilename The name of the file that will hold the output.
*/
class PlotProcessor : public marlin::Processor
{
public:
	//The usual Marlin processor methods
	virtual Processor* newProcessor() { return new PlotProcessor; }
	PlotProcessor();
	virtual ~PlotProcessor();
	virtual void init();
	virtual void processRunHeader( LCRunHeader* pRun );
	virtual void processEvent( LCEvent* pEvent );
	//don't need this
	//virtual void check( LCEvent* pEvent );
	virtual void end();
protected:
	std::string _JetCollectionName;	/**< @internal The name of the collection of ReconstructedParticles that is the jet (comes from the steering file).*/
	std::vector<std::string> _FlavourTagCollectionNames;	/**< @internal The names of the collection of LCFloatVec that are the flavour tags (a set of purity effiency plots will be made for each tag) (comes from the steering file).*/
	std::string _TrueJetFlavourColName; /**< @internal The name of the collection of LCIntVec that is the true jet flavour (comes from the steering file).*/
	std::string _OutputFilename; /**< @internal The filename of the output root file if using root, otherwise the directory and the first part of the filename of the comma seperated value files.*/
	std::vector<std::map<std::string,unsigned int> > _IndexOfForEachTag;
	int _nRun; /**< @internal The current run number.*/

	histogram_data<double> _jetEnergy; /**< @internal Custom storage class that holds all of the jet energies.*/

	std::vector<efficiency_purity<double> > _BTagEfficiencyPurity;/**< @internal Custom storage class that holds all the efficiency/purity data for the b tag calculated by the FlavourTagProcessor.*/
	std::vector<efficiency_purity<double> > _CTagEfficiencyPurity;/**< @internal Custom storage class that holds all the efficiency/purity data for the b tag calculated by the FlavourTagProcessor.*/
	std::vector<efficiency_purity<double> > _BCTagEfficiencyPurity;/**< @internal Custom storage class that holds all the efficiency/purity data for the b tag (only b background) calculated by the FlavourTagProcessor.*/

	//useful constants
	static const int C_JET=4;/**< @internal Useful constant for the jet flavour*/
	static const int B_JET=5;/**< @internal Useful constant for the jet flavour*/

	void _displayCollectionNames( lcio::LCEvent* pEvent );/**< @internal Just prints out the available collections in the LCIO file to standard output.*/
	bool _passesEventCuts( lcio::LCEvent* pEvent );	///< @internal A function that contains all the event cuts - returns true if the event passes all of the cuts, false otherwise.
	bool _passesJetCuts( lcio::ReconstructedParticle* pJet ); ///< @internal A function that contains all the jet cuts - returns true if the event passes all of the cuts, false otherwise.

	void _fillPlots( LCEvent* pEvent, unsigned int jet);/**< @internal Internal function that is just code split off from processEvent() to simplify it - fills the container classes with the data from the file.*/
	void _outputDataToFile( std::string filename );/**< @internal Internal function that is just code split off from end() to simplify it - writes the required data from the container classes to the output file.*/

	double _jetEMax;/**< @internal Keeps a record of the highest jet energy - gets printed to standard output at the end as a sanity check.*/
};

#endif //ifndef PlotProcessor_h
