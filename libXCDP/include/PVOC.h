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

	//Note this is currently saving to tga, not the easiest to open image format, just easy to make
	const PVOC & getSpectrograph( const std::string & fileName ) const;

	//======================================================
	//	Procs
	//======================================================

	/* 
	 * Frame Selectors
	 * These choose some frame from the input at given times in the output
	*/
	PVOC timeSelect( double length, RealFunc selector ) const;
	PVOC holdAtTimes( const std::vector<double> & timesToHold  ) const;

	/* 
	 * Time Interpolators
	 * These functions choose every (decimation) frames and discards the rest, filling the space with 
	 * interpolated data
	*/
	PVOC timeInterpolate_Constant( size_t decimation ) const;
	PVOC timeInterpolate_Linear( size_t decimation ) const;
	PVOC timeInterpolate_Spline( size_t decimation ) const;

	/* 
	 * Time Extrapolators
	 * These functions attempt to extend the spectral content of a PVOC beyond its
	 * normal domain.
	*/
	PVOC timeExtrapolate_Linear( int factor ) const;

	/* 
	 * Spectral Repitch
	 * This multiplies all frequencies by factor
	 * A lot of people seem to think this algorithm is the big reason to use a phase vocoder
	 * but I just don't think it sounds that good
	*/
	PVOC repitch( double factor ) const;

	//PVOC gate( RealFunc cutoff ) const;
	//PVOC gate( double cutoff ) const;
	//PVOC invert( int lowerBound, int upperBound ) const;
};

} // End namespace xcdp

