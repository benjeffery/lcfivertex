#ifndef PerEventIPFitterProcessor_h
#define PerEventIPFitterProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>
#include <fstream>

#include "inc/algo.h"
#include "inc/event.h"
#include "inc/vertex.h"
using namespace lcio ;
using namespace marlin ;
using vertex_lcfi::Vertex;
using vertex_lcfi::Event;

//!Determine IP position and error from the tracks in an event by simple fit.
/*!
<h4>Inputs</h4>
<table>
<tr><td>Name</td>                   <td>Type</td>                   <td>Represents</td></tr>
<tr><td>InputRPCollectionName</td>  <td>ReconstructedParticle</td>  <td>Tracks to be fit</td></tr>
</table>
<h4>Outputs</h4>
<table>
<tr><td>Name</td>                 <td>Type</td>    <td>Represents</td></tr>
<tr><td>VertexCollectionName</td> <td>Vertex</td>  <td>Fitted IP</td></tr>
</table>
<h4>Description</h4>
This processor fits a set of LCIO ReconstructedParticles (which must have an LCIO Track
attached to be used) to a common point. This is performed by iterative fitting, with
removal of the track with highest chi-squared at each iteration until the fit
reaches the probability threshold. If only one track remains then the default IP 
position and error are used. The result is stored as an LCIO Vertex.
<br>This processor is highly unoptimised and untuned, and may take a long time to execute
on a large set of tracks.
<br>Currently uses VertexFitterLSM (from ZVTOP) to perform fitting.

\param InputRPCollection Name of the ReconstructedParticle collection to be fit
\param VertexCollectionName Name of the output Vertex collection
\param DefaultIPPos Length 3 Float Vector of position (x,y,z) returned (as LCIO Vertex) if no fit is found
\param DefaultIPErr Length 6 Float Vector of covariance (lower symmetric) returned (as LCIO Vertex) if no fit is found 
\param ProbabilityThreshold Once the vertex is above this probability it is returned

\author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
*/

class PerEventIPFitterProcessor : public Processor {
  
 public:
  
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new PerEventIPFitterProcessor ; }
  PerEventIPFitterProcessor() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  
 protected:
std::string _InputRPCollectionName ;
  std::string _VertexCollectionName;
  vertex_lcfi::Algo<Event*,Vertex*>* _IPFitter;
  FloatVec _DefaultIPPos;
  FloatVec _DefaultIPErr;
  double _ProbThreshold;
  int _nRun ;
  int _nEvt ;
} ;

#endif



