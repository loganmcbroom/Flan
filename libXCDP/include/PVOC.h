#pragma once

#include "Types.h"
#include "PVOCBuffer.h"

namespace xcdp {

class Audio;

class PVOC : public PVOCBuffer
{
public:

	PVOC( const PVOCBuffer::Format & other ) : PVOCBuffer( other ) {}

	//======================================================
	//	Conversions
	//======================================================

	Audio convertToAudio() const;
	const PVOC & getSpectrograph( const std::string & fileName = std::string("tempSpectrograph.tga") ) const;

	//======================================================
	//	Procs
	//======================================================

	PVOC timeSelect( double length, RealFunc selector ) const;
	PVOC holdAtTimes( const std::vector<double> & timesToHold  ) const;
	PVOC constantTimeInterpolate( int factor ) const;
	PVOC linearTimeInterpolate( int factor ) const;
//	PVOC linearTimeExtrapolate( int factor ) const;
	PVOC repitch( double factor ) const;
	//PVOC gate( RealFunc cutoff ) const;
	//PVOC gate( double cutoff ) const;
	//PVOC invert( int lowerBound, int upperBound ) const;
};

} // End namespace xcdp

