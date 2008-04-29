#ifndef TypesafeCollection_h
#define TypesafeCollection_h

//======================================================================
//TypesafeCollection class - checks input collections for correct type
//======================================================================

#include "marlin/Processor.h" 
#include "lcio.h"
#include <iostream>
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/LCFloatVec.h"

template<class T>
class TypesafeCollection
{
public:
	TypesafeCollection( lcio::LCEvent* pEvent, std::string collectionName )
	{
		// Try and get the collection
		try
		{
			_pCollection=pEvent->getCollection( collectionName );
		}
		catch(...)
		{
			_error.str(""); // clear whatever was there beforehand
			_error << __FILE__ << "(" << __LINE__ << "): Unable to get the collection \"" << collectionName << "\" for event "
					<< pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << ".";

			_pCollection=0;
			return;
		}
		
		if( !_checkCollectionType() )  //The collection is not of the expected type
		{
			_error.str("");  //clear whatever was there beforehand
			_error << __FILE__ << "(" << __LINE__ << "): The jet collection \"" << collectionName << "\" for event " << pEvent->getEventNumber()
					<< " in run " << pEvent->getRunNumber() << " is not the correct type.";

			_pCollection=0;
		}
	}

	virtual ~TypesafeCollection(){}

	bool is_valid()
	{
		if( _pCollection==0 ) return false;
		else return true;
	}

	std::string last_error()
	{
		return _error.str();
	}

	int getNumberOfElements()
	{
		if( _pCollection==0 ) return 0;
		else return _pCollection->getNumberOfElements();
	}

	T* getElementAt( int element )
	{
		if( _pCollection==0 ) return 0;
		//TODO Modify the last error message to say way these are returning 0
		if( element<0 || element>=_pCollection->getNumberOfElements() ) return 0;

		return dynamic_cast<T*>( _pCollection->getElementAt( element ) );
	}
private:
	lcio::LCCollection* _pCollection;
	std::stringstream _error; ///< Holds a string explaining the last error encountered.
	
//	template<class T>
	bool _checkCollectionType();
};



#endif
