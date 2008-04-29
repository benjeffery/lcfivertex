#ifndef VertexChargeProcessor_h
#define VertexChargeProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>

#include "inc/algo.h"
#include "inc/decaychain.h"
#include "inc/jet.h"
#include "util/inc/projection.h"


using namespace lcio ;
using namespace marlin ;
using vertex_lcfi::DecayChain;
using vertex_lcfi::Jet;
using vertex_lcfi::util::Projection;

/** Calculates the Vertex Charge. 
 *  The processor calculated the vertex charge of a Decay Chain by using the tracks and the verteces present in the chain. Two logically slightly
 * different algorithms are used depending on the hypothesis of a B or C quark Vertex. In the B hypothesis we include the inner verteces, in the C we do not include them. This choice is controlled by the parameter ChargeAllSecondaryTracks.  
 * 
 * <H4>Input</H4>
 * - A collection of ReconstructedParticles that represents the jets in the event (obtained from a jet
 * finder, say SatoruJetFinderProcessor, although in order not to break the reconstruction chain we suggest you run this after the flavour tagging. In this way the LCFI chain remains intact).
 * - A collection of vertices that contains the per event primary vertices; one for each event. (optional) This collection is filled in the vertex_lcfi::PerEventIPFitter processor.
 * - A collection of decay chains as filled by the the ZVTOPZVRESProcessor or ZVTOPZVKINProcessor. 
 *
 * <H4>Output</H4>
 * The processor writes into the selected lcio output file. All the values calculated by the processor are saved in the same LCFloatVec collection.
 * The default name of the output collection is Charge.Although this is changed in the steering files to something more appropriate, depending on B or C calculation. For more details see \ref LCIO "the interface documentation".
 *
   @param VertexChargeCollection collection where results will be stored. 
 * @param ChargeAllSecondaryTracks include or exclude tracks in the inner vertices for the track attachment.
*  @param ChargeCloseapproachCut upper cut on track distance of closest approach to the seed axis in the calculation of the vertex charge variable, used by vertex_lcfi::TrackAttach.
 *  @param ChargeLoDCutmax Cut determining the maximum L/D for the Charge, used by vertex_lcfi::TrackAttach (when calculating C-Charge)
 *  @param ChargeLoDCutmin Cut determining the minimum L/D for the Charge, used by vertex_lcfi::TrackAttach (when calculating C-Charge)
 *
 *
 *  @author Erik Devetak(erik.devetak1@physics.ox.ac.uk)
*/
class VertexChargeProcessor : public Processor {
  
 public:
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new VertexChargeProcessor ; }
  VertexChargeProcessor() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  protected:
  std::string _JetRPColName;
  std::string _DecayChainRPColName;
  std::string _RelationColName;
  std::string _IPVertexCollectionName;
  std::string _VertexChargeCollectionName;
  
  std::vector<std::string> _JetVariableNames;
  
  vertex_lcfi::Algo<DecayChain*,double>* _VertexCharge;
  vertex_lcfi::Algo<DecayChain*,DecayChain* >* _Attach;

  bool _ChargeAddAllTracksFromSecondary;
  double _ChargeLoDCutmin;
  double _ChargeLoDCutmax;
  double _ChargeCloseapproachCut;

  int _nRun ;
  int _nEvt ;
} ;

#endif



