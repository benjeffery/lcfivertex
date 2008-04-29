#ifndef RPCutProcessor_h
#define RPCutProcessor_h 1

#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>

#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#ifdef MCFAIL_DIAGNOSTICS
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#endif

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
using namespace lcio ;
using namespace marlin ;


//!Cuts ReconstuctedParticles(RPs) from a collection (or from a list of RPs held by another RP) based on several cut criteria
/*!
<h4>Input</h4>
<table>
<tr><td>Name</td>                   <td>Type</td>                   <td>Represents</td></tr>
<tr><td>InputRCPCollection</td>     <td>ReconstructedParticle</td>  <td>Collection to be cut</td></tr>
</table>
<h4>Output</h4>
<table>
<tr><td>Name</td>                 <td>Type</td>                  <td>Represents</td></tr>
<tr><td>OutputRCPCollection</td>  <td>ReconstructedParticle</td> <td>If WriteNewCollection=true contains the RPs that passed the cuts.</td></tr>
</table>
<h4>Description</h4>
Based on several criteria this processor removes RPs from a collection,
or if SubParticleLists = true then it removes RPs held by RPs in a collection.
<br>
Depending on WriteNewCollection, the Output is either the original collection 
with the RPs removed, or a new collection with the input collection remaining untouched.
<br>
NOTE - A track is cut if its ReconstructedParticle has no Track objects attached.
<br>Most cuts follow a standard set of parameters:
<table border=0 margin=0>
<tr><td><strong>a1_{CutName}Enable</strong></td> <td>If true the cut is enabled</td></tr>
<tr><td><strong>a2_{CutName}CutLowerThan</strong></td> <td>If true RPs with a value lower than the cut value are cut, if false those higher than the cut value are cut.</td></tr>
<tr><td><strong>a3_{CutName}CutValue</strong></td> <td>The value of the cut</td></tr>
</table>
(The letter and number index prefixed to each parameter are to ensure they stay together in the
output of Marlin -x)
The cuts that follow this scheme are:
<table border=0>
<tr><td>Name</td><td>Description</td></tr>
<tr><td></td><td></td></tr>
<tr><td>Chi2OverDOF</td><td>Chi squared over degrees of freedom (Track::gethi2())</td></tr>
<tr><td>D0</td><td>Track D0 (Track::getD0())</td></tr>
<tr><td>D0Err</td><td>D0 Covariance (Track::getCovMatrix()[0])</td></tr>
<tr><td>Z0</td><td>Track Z0 (Track::getZ0())</td></tr>
<tr><td>Z0Err</td><td>Z0 Covariance (Track::getCovMatrix()[9])</td></tr>
<tr><td>PT</td><td>Track Pt (rPhi projection of Track::getMomentum())</td></tr>
</table>
<h4>Subdetector hits</h4>
The cut on subdetector hits relies on information in Track::getSubdetectorHitNumbers() this is
an array. The processor is told what order the detectors are in this array by parameter "g2_SubDetectorNames" which is the sting names of the detectors in the same order. The other parameters then use these names.
<br>
<h4>MC PID of Parents</h4>
This cut uses MC information provided by the MCParticleRelation collection to cut particles whose parents have a PID in the list provided by parameter "h2_CutPIDS"
<br>
<h4>Bad parameters cut</h4>
If enabled by "i1_BadParametersEnable" this cut removes tracks with nan covariances and parameters.
<br>
<h4>MC Vertex cut</h4>
Experimental MC Cut - most likely removed in next release
\author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
*/
class RPCutProcessor : public Processor {
  
 public:
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new RPCutProcessor ; }
  RPCutProcessor() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  protected:

  bool  _Chi2OverDOFFail(ReconstructedParticle* RPTrack);
  bool  _D0Fail(ReconstructedParticle* RPTrack);
  bool  _D0ErrFail(ReconstructedParticle* RPTrack);
  bool  _Z0Fail(ReconstructedParticle* RPTrack);
  bool  _Z0ErrFail(ReconstructedParticle* RPTrack);
  bool  _PTFail(ReconstructedParticle* RPTrack);
  bool  _DetectorHitsFail(ReconstructedParticle* RPTrack,std::map<std::string,int> SubdetectorIndex);
  bool  _MCPIDFail( lcio::ReconstructedParticle* RPTrack, UTIL::LCRelationNavigator* pMCRelationNavigator );
  bool  _BadParametersFail(lcio::ReconstructedParticle* RPTrack);
  bool  _MCVertexFail(lcio::ReconstructedParticle* RPTrack, UTIL::LCRelationNavigator* pMCRelationNavigator );
  
  std::string _InRCPColName ;
  std::string _TrackColName;
  std::string _OutRCPColName ;
  bool _WriteNewCollection;
  bool _SubParticleLists;
  
  bool _Chi2OverDOFEnable;
  bool _Chi2OverDOFCutLowerThan;
  float _Chi2OverDOFCutValue;
  
  bool _D0Enable;
  bool _D0CutLowerThan;
  float _D0CutValue;
  
  bool _D0ErrEnable;
  bool _D0ErrCutLowerThan;
  float _D0ErrCutValue;
  
  bool _Z0Enable;
  bool _Z0CutLowerThan;
  float _Z0CutValue;

  bool _Z0ErrEnable;
  bool _Z0ErrCutLowerThan;
  float _Z0ErrCutValue;
  
  bool _PTEnable;
  bool _PTCutLowerThan;
  float _PTCutValue;

  bool _MonteCarloPIDEnable;
  std::vector<int> _MonteCarloPIDsToCut;
  std::string _MonteCarloRelationColName;
  
  bool _BadParametersEnable;
  
  bool _DetectorHitsEnable;
  std::vector<std::string> _DetectorHitsBoundaryDetectorNames;
  std::vector<int> _DetectorHitsBoundaryCuts;
  std::vector<std::string> _DetectorHitsRegion1DetectorNames;
  std::vector<int> _DetectorHitsRegion1Cuts;
  std::vector<std::string> _DetectorHitsRegion2DetectorNames;
  std::vector<int> _DetectorHitsRegion2Cuts;
  std::vector<std::string> _DetectorNames;
  
  bool _MCVertexEnable;
  const gear::VXDParameters* _VxdPar;
  double _BeamPipeInnerRadius;
  double _BeamPipeOuterRadius;
  double _BeamPipeHalfZ;
  double _CutDistance;

#ifdef MCFAIL_DIAGNOSTICS
  TH2F *_diaghist_bpmat_xy;
  TH2F *_diaghist_vxmat_xy;
  TH2F *_diaghist_nomat_xy;
  TH2F *_diaghist_bpmat_rz;
  TH2F *_diaghist_vxmat_rz;
  TH2F *_diaghist_nomat_rz;
  TH1F *_diaghist_dist;
  TH1F *_diaghist_dist_vxmat;
  TH1F *_diaghist_dist_nomat;
#endif
  
  int _nRun ;
  int _nEvt ;
} ;

#endif



