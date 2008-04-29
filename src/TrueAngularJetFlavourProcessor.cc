#include "TrueAngularJetFlavourProcessor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <EVENT/LCParameters.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCIntVec.h>
#include "EVENT/LCFloatVec.h"

#include <UTIL/LCRelationNavigator.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <vector>
#include <string>
#include <map>

using namespace marlin ;
using namespace lcio;
using std::vector;
using std::string;
using std::map;
using EVENT::Track;

TrueAngularJetFlavourProcessor aTrueAngularJetFlavourProcessor ;

TrueAngularJetFlavourProcessor::TrueAngularJetFlavourProcessor() : Processor("TrueAngularJetFlavourProcessor") {
  
  // modify processor description
  _description = "TrueAngularJetFlavourProcessor - Determines the true flavour of a jet from the MC Paticles associated to the Jets RP" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "JetRPCollection" , 
			      "Name of the ReconstructedParticle collection that represents jets"  ,
			      _JetRPColName ,
			      std::string("Jets") ) ;
  registerInputCollection( lcio::LCIO::MCPARTICLE,
  			      "MCParticleCollection" , 
			      "Name of the collection that holds all MC particles. "  ,
			      _MCParticleColName ,
			      std::string("MCParticle") ) ;
  registerOutputCollection( lcio::LCIO::LCINTVEC,
  			      "TrueJetFlavourCollection" , 
			      "Name of the output collection of LCIntVec (same order as jets)"  ,
			      _TrueJetFlavourColName ,
			      std::string("TrueJetFlavour") ) ;
 registerOutputCollection( lcio::LCIO::LCINTVEC,
  			      "TrueJetPDGCodeCollection" , 
			      "Name of the output collection of LCIntVec (same order as jets)"  ,
			      _TruePDGCodeColName ,
			      std::string("TrueJetPDGCode") ) ;
 registerOutputCollection( lcio::LCIO::LCFLOATVEC,
  			      "TrueJetHadronChargeCollection" , 
			      "Name of the output collection of LCIntVec (same order as jets)"  ,
			      _TrueHChargeColName ,
			      std::string("TrueJetHadronCharge") ) ;
 registerOutputCollection( lcio::LCIO::LCFLOATVEC,
  			      "TrueJetPartonChargeCollection" , 
			      "Name of the output collection of LCIntVec (same order as jets)"  ,
			      _TruePChargeColName ,
			      std::string("TrueJetPartonCharge") ) ; 
 registerOptionalParameter( "MaximumAngle" ,
   			      "Maximum value allowed between MCParticle and jet momentum expressed in degrees"  ,
			     _MaximumAngle,
			      double(180));
  }

void TrueAngularJetFlavourProcessor::init() 
{ 
	// usually a good idea to
	printParameters() ;
	
	_nRun = 0 ;
	_nEvt = 0 ;
}

void TrueAngularJetFlavourProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void TrueAngularJetFlavourProcessor::processEvent( LCEvent * evt ) { 

	LCCollection* JetRPCol = evt->getCollection( _JetRPColName );
       	LCCollection* MCParticleCol = evt->getCollection( _MCParticleColName );

	std::vector<ReconstructedParticle*> MyJets;
	std::vector<int> MCflavour;
	std::vector<float> MCCharge;
	std::vector<float> MCPCharge;
	std::vector<int> MCCode;


	LCCollectionVec* OutCollection = new LCCollectionVec("LCIntVec");
	LCCollectionVec* OutCollectionPDG = new LCCollectionVec("LCIntVec");
	LCCollectionVec* OutCollectionPCharge = new LCCollectionVec("LCFloatVec");
	LCCollectionVec* OutCollectionHCharge = new LCCollectionVec("LCFloatVec");

	evt->addCollection(OutCollection,_TrueJetFlavourColName);	
	evt->addCollection(OutCollectionPDG,_TruePDGCodeColName);	
	evt->addCollection(OutCollectionPCharge,_TruePChargeColName);	
	evt->addCollection(OutCollectionHCharge,_TrueHChargeColName);	
       	
	int nRCP = JetRPCol->getNumberOfElements()  ;

	if(nRCP ==0 ) std::cerr<<"Warning: TrueAngularFetFlavourProcessor.cc:88 : NO jets present "<<std::endl; 

	for(int i=0; i< nRCP ; i++)
	{
	  MyJets.push_back(dynamic_cast<ReconstructedParticle*>(JetRPCol->getElementAt(i))); 
	}
	
	vector<int> jetflavourdata(MyJets.size(),1);
	vector<int> jetcodedata(MyJets.size(),0);
	vector<double> jetangledata(MyJets.size(),1000);
	
	// Set default charge for light jets to -100 
	vector<float> hadronchargedata(MyJets.size(),-100);
	vector<float> partonchargedata(MyJets.size(),-100);

	std::vector<std::vector<int> > Ntrkjet;
	bool originalhadron = 0;
	int code = 0;

	float charge =-100;
	float pcharge =-100;

	std::vector<MCParticle*> HeavyMCs;
	std::vector<MCParticle*> JetMCs;
	std::vector<double>  TrueFlavour;
	
	std::vector<int> typepart;
	std::vector<MCParticle*> daughtercharge;


	if( MCParticleCol->getNumberOfElements() == 0 ) std::cerr<<"Warning: TrueAngularFetFlavourProcessor.cc:107 : NO MC data presentjets present "<<std::endl; 


	for(int nMCpart = 0; nMCpart< MCParticleCol->getNumberOfElements();nMCpart++ )
	  {

	    pcharge = -100;
	    MCParticle* MCused = dynamic_cast<MCParticle*> (MCParticleCol->getElementAt(nMCpart));
	    
	    //note this whole process relies on 551/100 = 5 in integer cast!!! (which is true here)
	      
	    code = abs( MCused->getPDG() );
	    charge = MCused->getCharge();
	    
	    //NOTE: this works only with mokka-06-04-p02 or above. 
	    // In case you are looking at older versions of mokka please remove the if statement. 
	    if(charge == -1000  && code >100)
	      {
		charge = chargefromPDG(MCused->getPDG());
	      }
						
	    //this little addition will take care of wierder mesons 
	    if( code  > 10000   )
	      {
		code  = code%1000;
	      }
	    
	    if( code  > 1000   )
	      {
		
		if( code/ 1000 == 4   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(4);
		    MCCharge.push_back(charge);
		    MCCode.push_back(MCused->getPDG());
		    if( MCused->getPDG() > 0  )
		      { pcharge = 1; }
		    else
		      { pcharge = -1;}
		    MCPCharge.push_back(pcharge);
		  }
		
		if( code/ 1000 == 5   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(5);
		    MCCharge.push_back(charge);
		    MCCode.push_back(MCused->getPDG());
		    if( MCused->getPDG() > 0 )
		      {pcharge = 1; }
		    else
		      {pcharge = -1;}
		    MCPCharge.push_back(pcharge);
		  }		
	      }
	    
	    if( code  > 100   )
	      {
		if( code / 100  == 4   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(4);
		    MCCharge.push_back(charge);
		    MCCode.push_back(MCused->getPDG());
		  
		    if ( (code/10) %11 == 0)
		      {
			pcharge = 0;
		      }
		    else
		      {
			if( MCused->getPDG() > 0  )
			  { pcharge = 1; }
			else
			  { pcharge = -1;}
		      }
		    MCPCharge.push_back(pcharge);
		  }    
		
		if( code/100 == 5   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(5);
		    MCCharge.push_back(charge);
		    MCCode.push_back(MCused->getPDG());
		    
		    if ( code %11 == 0)
		      {
			pcharge = 0;
		      }
		    else
		      {
			if( MCused->getPDG() < 0  )
			  { pcharge = 1; }
			else
			  { pcharge = -1;}
		      }
		    MCPCharge.push_back(pcharge);
		  }
		
	      }
	    //std::cout<< "hadron,  "<<MCused->getPDG() << ",    "<< pcharge<<std::endl;
	  }


	for (std::vector<MCParticle*>::const_iterator iParticle = HeavyMCs.begin(); iParticle != HeavyMCs.end() ;++iParticle)	
	  {
	    originalhadron =0;
	    MCParticle* MCused = *iParticle;
	    while(!originalhadron)
	      {
		// at hadron level particles can have only 1 parent only at parton level 
		// there might be more than 1 parent. if we are at parton level we are too far anyway.
		// refer to LCIO manual for details.
		if( MCused->getParents().size() == 1) 
		  {
		    //		    std::cout<<" ParentalPDG    "<<MCused->getParents()[0]->getPDG()<<std::endl;
		    if (abs( MCused->getParents()[0]->getPDG()) >100  || (abs(MCused->getParents()[0]->getPDG())>10 && abs(MCused->getParents()[0]->getPDG())<81))
		      {
			MCused = MCused->getParents()[0];
		      }
		    else
		      {

			//check if already in the vector

			std::vector<MCParticle*>::const_iterator Doubles = find( JetMCs.begin() ,JetMCs.end() , MCused);
			if(Doubles == JetMCs.end()|| JetMCs.size() ==0 )
			  {
			    //	    std::cout<< MCused->getPDG() <<std::endl;
			    JetMCs.push_back(MCused);
			  }  
			originalhadron = 1;
		      }
		  }
		else
		  {
		    //check if already in the vector
			//unique does not work
		    std::vector<MCParticle*>::const_iterator Doubles = find(JetMCs.begin(),JetMCs.end(), MCused);
		    if(Doubles == JetMCs.end()|| JetMCs.size() ==0)
		      {
			//			std::cout<<"NUMBER    "<<MCused->getPDG()<<std::endl;
			JetMCs.push_back(MCused);
		      }
			originalhadron = 1;
		  }
	      }
	  }

	//	std::cout<<JetMCs.size()<<std::endl;

	if(MyJets.size()>0)
	  {


		int  MCcounter =0;	    
	    for(std::vector<MCParticle*>::const_iterator iParticle3 = JetMCs.begin(); iParticle3 != JetMCs.end() ;++iParticle3)
	      {
		double dist = -999999999; 
		int  JetCounter =0;
		int JetValue =-1;
		double angle = 0;
		double anglestored = 100000;
		vector<double> normMomMC;
		const double* MomentumMC = (*iParticle3)->getMomentum();
		double length = sqrt(MomentumMC[0]*MomentumMC[0]+MomentumMC[1]*MomentumMC[1]+MomentumMC[2]*MomentumMC[2]);
		normMomMC.push_back(MomentumMC[0]/length);
		normMomMC.push_back(MomentumMC[1]/length);
		normMomMC.push_back(MomentumMC[2]/length);
		
		//		std::cout<<"numberUSED    "<< (*iParticle3)->getPDG() <<std::endl;

		for(std::vector<ReconstructedParticle*>::const_iterator iParticle2 = MyJets.begin(); iParticle2 != MyJets.end() ;++iParticle2)
		  {
		    vector<double> normMom;
		    const double* Momentum = (*iParticle2)->getMomentum();
		    double length = sqrt(Momentum[0]*Momentum[0]+Momentum[1]*Momentum[1]+Momentum[2]*Momentum[2]);
		    normMom.push_back(Momentum[0]/length);
		    normMom.push_back(Momentum[1]/length);
		    normMom.push_back(Momentum[2]/length);
		    
		    double disttemp= sqrt(((normMomMC[0]-normMom[0])*(normMomMC[0]-normMom[0]))+
					  ((normMomMC[1]-normMom[1])*(normMomMC[1]-normMom[1]))+
					  ((normMomMC[2]-normMom[2])*(normMomMC[2]-normMom[2])));

		    //		    std::cout<<"Disttemp    "<<JetCounter<< "      "<<disttemp<<std::endl;

		   

		    if (disttemp < fabs(dist))
		      {
			 angle = 180*asin( disttemp/2 )*2/3.14159265;
	
			if(angle<_MaximumAngle)
			  {
			    dist = disttemp;
			    JetValue = JetCounter;
			    anglestored = angle;
			    //			    std::cout<<"Inside    "<<(*iParticle3)->getPDG()<<"           "<<JetValue<<"     "<<anglestored<<std::endl;
			  }
			else
			  {
      			    std::cout<< "Heavy MC Particle not assigned due to angle cut   "<<std::endl;
			  }
		      }
		    
		    JetCounter++;
		  }

		if(dist>0 && JetValue >=0)
		  {
		    if(jetflavourdata[JetValue]>3)
		      {
			if( jetflavourdata[JetValue] == MCflavour[MCcounter] )
			  {
			    if(fabs(jetangledata[JetValue]) > fabs(anglestored))
			      {			    

				hadronchargedata[JetValue] =  MCCharge[MCcounter];
				partonchargedata[JetValue] =  MCPCharge[MCcounter];
				jetcodedata[JetValue] =  int(MCCode[MCcounter]);
				jetangledata[JetValue] = anglestored;

			      }
			  }
			if( jetflavourdata[JetValue] < MCflavour[MCcounter])
			  {
			    jetflavourdata[JetValue] =  MCflavour[MCcounter];
			    hadronchargedata[JetValue] =  MCCharge[MCcounter];
			    partonchargedata[JetValue] =  MCPCharge[MCcounter];
			    jetcodedata[JetValue] =  int(MCCode[MCcounter]);
			    jetangledata[JetValue] = anglestored;
			  }
		      }
		    else
		      {
			jetflavourdata[JetValue] =  MCflavour[MCcounter];
			hadronchargedata[JetValue] =  MCCharge[MCcounter];
			partonchargedata[JetValue] =  MCPCharge[MCcounter];
			jetcodedata[JetValue] =  int(MCCode[MCcounter]);
			jetangledata[JetValue] = anglestored;
		      }
		  }

		MCcounter++;

	      }
	  }
	
	for(unsigned int iii=0; iii<jetflavourdata.size(); iii++)
	  {
	    LCIntVec *OutVec = new LCIntVec();    
	    LCIntVec *OutVecPDG = new LCIntVec();    
	    LCFloatVec *OutVecPCharge = new LCFloatVec();    
	    LCFloatVec *OutVecHCharge = new LCFloatVec();    
	    OutVec->push_back(jetflavourdata[iii]);
	    OutVecPDG->push_back(jetcodedata[iii]);
	    OutVecPCharge->push_back(partonchargedata[iii]);
	    OutVecHCharge->push_back(hadronchargedata[iii]);
	    OutCollection->addElement(OutVec);
	    OutCollectionPDG->addElement(OutVecPDG);
	    OutCollectionPCharge->addElement(OutVecPCharge);
	    OutCollectionHCharge->addElement(OutVecHCharge);
	    std::cout << "TJ:  Jet # "<< iii<<  "    PDG:   "; 
       	    std::cout <<jetcodedata[iii]<< "   Flavour:   ";
       	    std::cout <<jetflavourdata[iii]<< " Parton Charge:    ";
       	    std::cout <<partonchargedata[iii] << "   Hadron Charge:    ";
       	    std::cout <<hadronchargedata[iii] <<std::endl;
	    };
	//Create the collection to store the results

	_nEvt ++ ;
}



void TrueAngularJetFlavourProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrueAngularJetFlavourProcessor::end(){ 
	
	std::cout << "TrueAngularJetFlavourProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

//reconstruct quark content from pdg code
//I am surprised I could not find this anywhere already coded
//This is needed because sometimes mokka gives us the default -1000 and not the correct charge
float TrueAngularJetFlavourProcessor::chargefromPDG(int code)
{
  double first = 0; 
  double second = 0;
  double third = 0;
  double charge = -1000;

  if( abs(code)  > 1000000000  )
    {
      std::cout<<" Note: TrueAngularFetFlavourProcessor.cc: Hevy nucleus found and ignored. "<<std::endl; 
      return -1000;
    }

  if( abs(code)  > 10000   )
    {
      code  = code%1000;
    }
  if(fabs(code) >1000)
    {
      first = (code/1000) %2;
      second = (code/100) %2;
      third = (code/10) %2;  

      if( third == 0)
	{third = third+2;}
      else
	{third = -1; }
      if (second  == 0)
	{second = second+2;}
      else
	{second = -1; }
      if (first ==0)
	{first = first+2;}
      else
	{first = -1; }
      charge = (first+second+third)/3;
      if(charge ==0)
	{return 0;}
      else if (code>0)
	{return charge;}
      else
	{return (-charge);}
    }
  else  if (abs(code)>100)
    {
      first = (code/100) % 2;
      second = (code/10) %2;
      if (second  ==0)
	{second = second+2;}
      if (first ==0)
	{first = first+2;}		 
      charge = (first-second)/3;
      if (charge != 0)
	{
	  if(code>0)
	    {return 1;}
	  if(code<0)
	    {return -1;}		   
	}
      else 
	{return 0;}
    }
  else
    {
      std::cerr<<"Warning: TrueAngularFetFlavourProcessor.cc: A non hadronic particle has been choosen for vertex calculations "<<std::endl; 
      return -1000;
    }  
  // if we are here something probably went wrong. 
  return -1000;
}
