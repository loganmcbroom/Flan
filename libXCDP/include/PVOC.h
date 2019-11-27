#pragma once

#include <string>
#include <functional>

#include "PVOCBuffer.h"
#include "RealFunc.h"

namespace xcdp {

class Audio;
struct RealFunc;

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

	PVOC timeSelect( double length, RealFunc selector ) const;
	PVOC holdAtTimes( const std::vector<double> & timesToHold  ) const;

	//========================================================================
	//	Resampling
	//========================================================================

	typedef std::function< double ( double mix, double a, double b ) > Interpolator;
	static const Interpolator constantInterpolator, linearInterpolator, sineInterpolator;

	PVOC decimate( RealFunc decimation ) const;
	PVOC interpolate( RealFunc interpolation, Interpolator interpolator = linearInterpolator ) const;
	PVOC interpolate_spline( RealFunc interpolation ) const; //Spline interpolation requires custom handling
	PVOC resample( RealFunc decimation, RealFunc interpolation, Interpolator interpolator = linearInterpolator ) const;
	PVOC desample( RealFunc decimation, RealFunc interpolation, Interpolator interpolator = linearInterpolator ) const;
	PVOC timeExtrapolate( double startTime, double endTime, double extrapTime, Interpolator interpolator = linearInterpolator ) const;

	//========================================================================
	// Combinations
	//========================================================================

	PVOC replaceAmplitudes( const PVOC & ampSource, RealFunc amount = 1 ) const;
	PVOC subtractAmplitudes( const PVOC & other, RealFunc amount = 1 ) const;

	//========================================================================
	// Uncategorized
	//========================================================================

	PVOC repitch( double factor ) const;

	typedef std::function< double ( double, double ) > Perturber;
	static const Perturber identityPerturber, normalDistPerturber, normalDistUpPerturber, 
		normalDistDownPerturber;
	PVOC perturb( RealFunc magAmount, RealFunc frqAmount, 
		Perturber magPerturber = identityPerturber, Perturber frqPerturber = normalDistPerturber ) const;

	PVOC retainLoudestPartials( RealFunc numBins ) const;

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

