#ifndef ConversionTagger_h
#define ConversionTagger_h 1

#include "HistMap.h"

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>


using namespace lcio;
using namespace marlin;
using namespace std;


class ConversionTagger : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ConversionTagger ; }
  ConversionTagger() ;
  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void end() ;
  
  
 private:

  void tagger( LCEvent *evt, const string collectionName);
  ReconstructedParticle* CreateRecoPart(Track* trk);
  double diParticleMass(float* mom1, float* mom2,
			double mass1, double mass2);

  HistMap* histos;

  double _BField;
  double _twopi;
} ;

#endif



