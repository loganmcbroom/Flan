#pragma once

#include <string>
#include <functional>

#include "PVOCBuffer.h"
#include "RealFunc.h"
#include "Utility.h"

namespace xcdp {

class Audio;

class PVOC : public PVOCBuffer
{
public:

	typedef const std::vector<const PVOC &> & Vec;

	PVOC( const PVOCBuffer::Format & other ) : PVOCBuffer( other ) {}

	//========================================================================
	//	Conversions
	//========================================================================

	Audio convertToAudio() const;

	//Note this is currently saving to tga, not the easiest to open image format, just easy to make
	const PVOC & getSpectrograph( const std::string & fileName ) const;

	//========================================================================
	// Selection
	//========================================================================

	PVOC getFrame( double time ) const;
	PVOC timeSelect( double length, Surface selector ) const;
	PVOC holdAtTimes( const std::vector<double> & timesToHold  ) const;

	//========================================================================
	//	Resampling
	//========================================================================

	PVOC decimate( RealFunc decimation ) const;
	PVOC interpolate( RealFunc interpolation, Interpolator = Interpolators::linear ) const;
	PVOC interpolate_spline( RealFunc interpolation ) const; //Spline interpolation requires custom handling
	MFPair getBinInterpolated( size_t channel, double frame, size_t bin, Interpolator = Interpolators::linear ) const;
	MFPair getBinInterpolated( size_t channel, size_t frame, double bin, Interpolator = Interpolators::linear ) const;
	PVOC resample( RealFunc decimation, RealFunc interpolation, Interpolator = Interpolators::linear ) const;
	PVOC desample( RealFunc decimation, RealFunc interpolation, Interpolator = Interpolators::linear ) const;
	PVOC timeExtrapolate( double startTime, double endTime, double extrapTime, Interpolator = Interpolators::linear ) const;

	//========================================================================
	// Combinations
	//========================================================================

	PVOC replaceAmplitudes( const PVOC & ampSource, RealFunc amount = 1 ) const;
	PVOC subtractAmplitudes( const PVOC & other, RealFunc amount = 1 ) const;

	//========================================================================
	// Uncategorized (aka new stuff)
	//========================================================================

	PVOC perturb( RealFunc magAmount, RealFunc frqAmount, 
		Perturber magPerturber = Perturbers::identity, Perturber frqPerturber = Perturbers::normalDist ) const;

	PVOC retainNLoudestPartials( RealFunc numPartials ) const;
	PVOC removeNLoudestPartials( RealFunc numPartials ) const;

	PVOC resonate( double length, Surface decay ) const;

	PVOC repitch( Surface factor ) const;

	//UNFINISHED
	PVOC stretch( Surface factor ) const;
	PVOC convolve( const PVOC & other ) const;
	PVOC accumulate( double length, Surface decay ) const;

	//========================================================================
	// CDP Map
	//========================================================================

	// Blur
	PVOC blur_blur( RealFunc blurring ) const;
	PVOC blur_chorus( RealFunc fspread ) const; //mode 2

	// Combine
	PVOC combine_cross( const PVOC & infile2, RealFunc interp = 1 ) const;

};

} // End namespace xcdp

