#ifndef FlavourTagInputsProcessor_h
#define FlavourTagInputsProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>

#include "inc/algo.h"
#include "inc/decaychain.h"
#include "inc/jet.h"
#include "util/inc/projection.h"
#include "algo/inc/paramsignificance.h"
#include "algo/inc/twotrackpid.h"
#include "algo/inc/decaysignificance.h"

using namespace lcio ;
using namespace marlin ;
using vertex_lcfi::DecayChain;
using vertex_lcfi::Jet;
using vertex_lcfi::util::Projection;
using vertex_lcfi::SignificanceType;
using vertex_lcfi::PidCutType;
using vertex_lcfi::DecaySignificanceType;
/** Calculates the Flavour tag input variables for flavour tagging. 
 *
 * The aim of the processor is to calculate a series of highly discriminating tagging variables.
 * At present the default variables are the one defined in the R. Hawking LC note LC-PHSM-2000-021.
 * All the variables are calculated inside independent classes that inherit from the vertex_lcfi::Algo template and not in the main 
 * processor file. This makes the processor file extremely flexible and new variables easy to add. 
 * Similary it is also very simple to remove undesired variables. 
 * The following variables are presently calculated (variables depending on the vertex_lcfi::TrackAttach procedure are marked by *, variables depending on vertex_lcfi::TwoTrackPid are marked by ^)
 * <br> D0Significance1 - calculated in vertex_lcfi::ParameterSignificance
 * <br> D0Significance2 - calculated in vertex_lcfi::ParameterSignificance
 * <br> Z0Significance1 - calculated in vertex_lcfi::ParameterSignificance
 * <br> Z0Significance2 - calculated in vertex_lcfi::ParameterSignificance
 * <br> Momentum1 - calculated in vertex_lcfi::ParameterSignificance
 * <br> Momentum2 - calculated in vertex_lcfi::ParameterSignificance
 * <br> JointProbRPhi - ^ calculated in vertex_lcfi::JointProb
 * <br> JointProbZ -  ^ calculated in vertex_lcfi::JointProb 
 * <br> DecayLengthSignificance - calculated in vertex_lcfi::VertexDecaySignificance
 * <br> DecayLength - calculated in vertex_lcfi::VertexDecaySignificance
 * <br> PTCorrectedMass - * calculated in vertex_lcfi::VertexMass
 * <br> RawMomentum - * calculated in vertex_lcfi::VertexMomentum
 * <br> NumTracksInVertices - calculated in vertex_lcfi::VertexMultiplicity
 * <br> SecondaryVertexProbability - * calculated in vertex_lcfi::SecVertexProb 
 *
 * <br> For more information about the algorithms themselves please consult the specific algorithm documentation pages.  
 * The processor also uses the following algorithms: 
 * <br> vertex_lcfi::TwoTrackPid - algorithm that calulates the id of two charged tracks by using mass considerations. This algorithm removes 
 * tracks consistent with the hypothesis that they have been generated from Ks and gamma. This algorithm is not used in ZVTOPZVRESProcessor or ZVTOPZVKINProcessor
 * <br> vertex_lcfi::TrackAttach - algorithm that adds tracks close to the seed vertex. 
 *
 * <H4>Input</H4>
 * - A collection of ReconstructedParticles that represents the jets in the event (obtained from a jet
 * finder, say SatoruJetFinderProcessor).
 * - A collection of vertices that contains the per event primary vertices; one for each event. (optional) This collection is filled in the vertex_lcfi::PerEventIPFitter processor.
 * - A collection of decay chains as filled by the the ZVTOPZVRESProcessor or ZVTOPZVKINProcessor. 
 *
 *
 * <H4>Output</H4>
 * The processor writes into the selected lcio output file. All the values calculated by the processor are saved in the same LCFloatVec collection.
 * The default name of the output collection is FlavourTagInputs.For more details see \ref LCIO "the interface documentation".
 *
 * @param JetRPCollection Name of the ReconstructedParticle collection that represents jets.
 * @param IPVertexCollection Name of the Vertex collection that contains the primary vertices (optional)
 * @param FlavourTagInputsCollection Name of the LCFloatVec Collection that will be created to contain the flavour tag inputs
 *
 * <br>  The following parameters are parameters for the algorithms used by the processor.These parameters are all optional.
 *  @param LayersHit Momentum cuts will be applied on number of LayersHit and LayersHit minus one, used by vertex_lcfi::ParameterSignificance
 *  @param AllLayersMomentumCut Cut on the minimum momentum if track hits LayersHit, used by  vertex_lcfi::ParameterSignificance
 *  @param AllbutOneLayersMomentumCut Cut on the minimum momentum if track hits LayersHit minus one, used by  vertex_lcfi::ParameterSignificance
 *  @param JProbMaxD0Significance Maximum d0 significance of tracks used to calculate the joint probability, used in vertex_lcfi::JointProb
 *  @param JProbMaxD0andZ0 Maximum d0 and z0 of tracks used to calculate the joint probability, used in vertex_lcfi::JointProb
 *  @param PIDChi2Cut Cut on the Chi squared of two tracks being in the same vertex, used by vertex_lcfi::TwoTrackPid
 *  @param PIDMaxGammaMass Cut on the upper limit of the photon candidate mass, used by vertex_lcfi::TwoTrackPid
 *  @param PIDMaxKsMass Cut on the upper limit of the Ks candidate mass, used by vertex_lcfi::TwoTrackPid
 *  @param PIDMinKsMass Cut on the lower limit of the Ks candidate mass, used by vertex_lcfi::TwoTrackPid
 *  @param PIDRPhiCut Cut on the maximum RPhi of the Ks/gamma decay vertex candidate, used by vertex_lcfi::TwoTrackPid
 *  @param PIDSignificanceCut Cut on the minimum RPhi significance of the tracks, used by vertex_lcfi::TwoTrackPid
 *  @param SecondVertexNtrackscut Cut on the minimum number of tracks in the seed vertex, used by vertex_lcfi::SecVertexProb
 *  @param SecondVertexProbChisquarecut Cut on the Chi Squared of the seed vertex, used by vertex_lcfi::SecVertexProb
 *  @param TrackAttachAllSecondaryTracks include or exclude tracks in the inner vertices for the track attachment.
 *  @param TrackAttachCloseapproachCut upper cut on track distance of closest approach to the seed axis used by vertex_lcfi::TrackAttach (when used for * flagged variables)
 *  @param TrackAttachLoDCutmax Cut determining the maximum L/D for the track attachment, used by vertex_lcfi::TrackAttach (when used for * flagged variables)
 *  @param TrackAttachLoDCutmin Cut determining the minimum L/D for the track attachment, used by vertex_lcfi::TrackAttach (when used for * flagged variables)
 *  @param VertexMassMaxKinematicCorrectionSigma Maximum Sigma (based on error matrix) by which the vertex axis can move when kinematic correction is applied, used by  vertex_lcfi::VertexMass
 *  @param VertexMassMaxMomentumAngleCut Upper cut on angle between momentum of vertex and the vertex axis, used by  vertex_lcfi::VertexMass
 *  @param VertexMassMaxMomentumCorrection Maximum factor, by which vertex mass can be corrected, used by  vertex_lcfi::VertexMass
 *
 *  As a final remark one should notice that two additional values are stored in the Output LC Collection. 
 *  These are: 
 *  <br> NumVertices - number of vertices in the jet; used to determine what variables to use in the following flavour tag processor. Calculated in the processor.
 *  <br> DecayLength(SeedToIP)- distance from the vertex seed in the trackattach processor and IP. This veriable can be used for further analysis, but it is not used in flavour tagging. Calculated in the processor.
 *  @author Erik Devetak(erik.devetak1@physics.ox.ac.uk), 
 *  <br> interface by Ben Jeffery (ben.jeffery1@physics.ox.ac.uk)
*/
class FlavourTagInputsProcessor : public Processor {
  
 public:
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new FlavourTagInputsProcessor ; }
  FlavourTagInputsProcessor() ;
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
  std::string _FlavourTagInputsCollectionName;
  
  std::vector<std::string> _JetVariableNames;
  
  vertex_lcfi::Algo<DecayChain*,double>* _VertexMomentum;
  vertex_lcfi::Algo<DecayChain*,double>* _VertexMass;
  vertex_lcfi::Algo<DecayChain*,int>* _VerticesTrackMultiplicity;
  vertex_lcfi::Algo<DecayChain*,std::map<DecaySignificanceType, double> >* _VertexDecaySignificance;
  vertex_lcfi::Algo<DecayChain*,DecayChain* >* _TrackAttach;
  vertex_lcfi::Algo<DecayChain*,DecayChain* >* _BAttach;
  vertex_lcfi::Algo<DecayChain*,DecayChain* >* _CAttach;
  vertex_lcfi::Algo<DecayChain*,double >* _SecVertexProb;
  vertex_lcfi::Algo<Jet*,std::map<SignificanceType, double> >* _ParameterSignificance;
  vertex_lcfi::Algo<Jet*,std::map<Projection, double> >* _JointProb;
  vertex_lcfi::Algo<Jet*,std::map<PidCutType,std::vector< vertex_lcfi::Track* > > >* _TwoTrackPID; 
  double _VertexMassMaxMomentumAngle;					     
  double _VertexMassMaxKinematicCorrectionSigma;			     
  double _VertexMassMaxMomentumCorrection;			    
  bool  _TrackAttachAddAllTracksFromSecondary;
  double _TrackAttachLoDCutmin;
  double _TrackAttachLoDCutmax;
  double _TrackAttachCloseapproachCut;
  double _SecondVertexProbChisquarecut;
  double _SecondVertexNtrackscut;
  double _LayersHit;
  double _AllbutOneLayersMomentumCut;
  double _AllLayersMomentumCut;
  double _PIDMaxGammaMass;
  double _PIDMinKsMass;
  double _PIDMaxKsMass;
  double _PIDChi2Cut;
  double _PIDRPhiCut;
  double _PIDSignificanceCut;
  double _JProbMaxD0Significance;
  double _JProbMaxD0andZ0;
  FloatVec _JProbResolutionParameterRphi;
  FloatVec _JProbResolutionParameterZ;
  int _nRun ;
  int _nEvt ;
} ;

#endif



