#ifndef ZVTOPZVKINProcessor_h
#define ZVTOPZVKINProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>
#include <fstream>

#include "inc/algo.h"
#include "inc/decaychain.h"
#include "inc/jet.h"
using namespace lcio ;
using namespace marlin ;
using vertex_lcfi::DecayChain;
using vertex_lcfi::Jet;

//!Find vertices in a jet using kinematic ZVTOP-ZVKIN algorithm
/*!
<h4>Input</h4>
<table>
<tr><td>Name</td>                   <td>Type</td>                   <td>Represents</td></tr>
<tr><td>JetRPCollectionName</td>    <td>ReconstructedParticle</td>  <td>Jets to be Vertexed (eg from SatoruJetFinderProcessor</td></tr>
<tr><td>IPVertexCollectionName</td> <td>Vertex</td>                 <td>Event Interaction Point (eg from PerEventIPFitterProcessor) - optional can be manually specified</td></tr>
</table>
<h4>Output</h4>
<table>
<tr><td>Name</td>                             <td>Type</td>                  <td>Represents</td></tr>
<tr><td>DecayChainCollectionName</td>         <td>ReconstructedParticle</td> <td>Decay Chains (set of found vertices)</td></tr>
<tr><td>VertexCollection</td>                 <td>Vertex</td>                <td>Found vertices</td></tr>
<tr><td>DecayChainRPTracksCollectionName</td> <td>ReconstructedParticle</td> <td>Tracks used in Decay Chains and found vertices</td></tr>
</table>
<h4>Description</h4>
This processor finds vertices in a set of LCIO ReconstructedParticles (usually a jet)
using the algorithm ZVTOP(ZVKIN)
Also see (INSERT LINK TO ZVTOP DOC)
To be used each ReconstructedParticle must have an attached LCIO Track.
Note it is imperative that the tracks have well formed and preferably accurate
covariance matrices in d0 and z0. If the covariances are too small fake or no vertices
may be found. Too large and vertices will be combined.<br>
The algorithm also requires an interaction point in the form of an LCIO Vertex or a
manually set position and covariance. -  NOTE Only an ip at the origin (0,0,0) is supported
as the ghosttrack has that origin, this will hopefully be upgraded in a future release.<br>
The set of vertices forming a decay chain as output as set of LCIO Vertices and LCIO ReconstructedParticles
for details see \ref LCIO "the interface documentation"<br>
For more details on algorihmic parameters see the ZVKIN paper "zvkin.ps" in the doc directory. 

\author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)

\param JetRPCollection  Name of the ReconstructedParticle collection that represents jets
\param IPVertexCollection  Name of the Vertex collection that contains the primary vertex (Optional)  
\param DecayChainRPTracksCollectionName  Name of the ReconstructedParticle collection that represents tracks in output decay chains  
\param VertexCollection  Name of the Vertex collection that contains found vertices  
\param DecayChainCollectionName  Name of the ReconstructedParticle collection that holds RPs representing output decay chains  
\param ManualIPVertex  If false then the primary vertex from VertexCollection is used  
\param ManualIPVertexPosition Manually set position of the primary vertex (cm) - non origin IP not yet fully supported  
\param ManualIPVertexError  Manually set error matrix of the primary vertex (cm) (lower symmetric)  
\param MinimumProbability  If a vertex candidate has a probability below this it will not be considered - lower value results in more merging and lower vertex multiplicity  
\param InitialGhostWidth  Width in cm of the ghost inital ghosttrack also the smallest width it is allowed to have  
\param MaxChi2Allowed  The ghost track is widened until all forward jet tracks have a chi squared lower than this value  
\param OutputTrackChi2  If true the chi squared contributions of tracks to vertices is written to LCIO  
*/
class ZVTOPZVKINProcessor : public Processor {
  
 public:
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new ZVTOPZVKINProcessor ; }
  ZVTOPZVKINProcessor() ;  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  
 protected:
  std::string _JetRPCollectionName ;
  std::string _DecayChainRPTracksCollectionName;
  std::string _VertexCollectionName ;
  std::string _IPVertexCollectionName ;
  std::string _DecayChainCollectionName;
  //std::string _RelationCollectionName;
  vertex_lcfi::Algo<Jet*,DecayChain*>* _ZVKIN;
  bool _ManualPrimaryVertex;
  FloatVec _ManualPrimaryVertexPos;
  FloatVec _ManualPrimaryVertexErr;
  double _MinimumProbability;
  double _InitialGhostWidth;
  double _MaxChi2Allowed;
  bool _OutputTrackChi2;
  int _nRun ;
  int _nEvt ;
} ;

#endif



