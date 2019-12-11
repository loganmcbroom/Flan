#include "PVOC.h"

#include <algorithm>
#include <iostream>
#include <complex>
#include <array>
#include <fstream>

#include <fftw3.h>
#include "spline.h"

#include "RealFunc.h"
#include "WindowFunctions.h"
#include "Audio.h"

const double pi = std::acos( -1.0 );

/** TODO:
	blur scatter
	stretch time
	formants get
	formants put
	morph morph
	superaccu
	focus exag

Boilerplate:
	auto format = getFormat();
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			for( size_t bin = 0; bin < out.getNumBins(); ++bin )

	return out;

*/

namespace xcdp {

//========================================================================
//	Conversions
//========================================================================

//Phase Vocoder resynthesis. Comments are light, Audio::convertToPVOC has most of the info.
Audio PVOC::convertToAudio() const
	{
	std::cout << "Getting audio ... \n";
	const size_t numBins = getNumBins();
	const size_t frameSize = getFrameSize();
	const size_t hopSize = frameSize / getOverlaps();
	const size_t numHops = getNumFrames();
	const double binWidthD = getBinWidth();
	
	AudioBuffer::Format audioFormat;
	audioFormat.numChannels = getNumChannels();
	audioFormat.numSamples = numHops * hopSize;
	audioFormat.sampleRate = getSampleRate();
	Audio out( audioFormat );

	//Scale window for faster normalization later
	std::vector<double> window = window::Hann( frameSize );
	for( size_t i = 0; i < window.size(); ++i )
		window[i] /= frameSize * getOverlaps();

	std::vector<double> phaseBuffer( numBins );
	std::complex<double> *fftIn = (std::complex<double>*) fftw_alloc_complex( numBins );
	double * fftOut = fftw_alloc_real( frameSize );
	fftw_plan plan = fftw_plan_dft_c2r_1d( int(frameSize), (fftw_complex*) fftIn, fftOut, FFTW_MEASURE );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		///initial phase?
		std::fill( phaseBuffer.begin(), phaseBuffer.end(), 0 );

		for( size_t hop = 0; hop < numHops; ++hop )
			{
			//reconstruct phase and convert to rectangular, store in fftIn
			//It's the same as in the analysis but backwards so I'm not writing up every step
			for( size_t bin = 0; bin < numBins; ++bin )
				{
				const double magnitude = getBin( channel, hop, bin ).magnitude;

				const double trueFreq = getBin( channel, hop, bin ).frequency;
				const double binDeviation = trueFreq / binWidthD - double(bin);
				const double deltaPhase = 2.0 * pi * binDeviation / double(getOverlaps());
				const double expectedPhaseDiff = double(bin) * ( 2.0 * pi / double(getOverlaps()) );
				const double phaseDiff = deltaPhase + expectedPhaseDiff;
				phaseBuffer[bin] += phaseDiff;
				
				fftIn[bin] = std::polar( magnitude, phaseBuffer[bin] );
				}

			fftw_execute( plan );

			//Accumulate ifft output into audio buffer
			for( size_t relativeSample = 0; relativeSample < frameSize; ++relativeSample )
				{
				// - overlaps/2 is a consequency of marginal frames
				const int actualSample = (int(hop)-int(getOverlaps()/2)) * int(hopSize) + int(relativeSample);
				if( 0 <= actualSample && actualSample < out.getNumSamples() )
					{ 
					out.setSample( channel, actualSample, out.getSample( channel, actualSample ) +
						window[relativeSample] * fftOut[relativeSample] );
					}
				}
			}
		}
	fftw_destroy_plan( plan );
	fftw_free( fftIn );
	fftw_free( fftOut );

	return out;
	}

const PVOC & PVOC::getSpectrograph( const std::string & fileName ) const
	{
	std::cout << "Getting Spectrograph ... \n";

	size_t numFrames = getNumFrames();
	double maxMag = getMaxPartialMagnitude();

	//bin -> frame to match write order
	std::vector<std::vector<std::array<uint8_t,3>>> data( getNumFrames(), std::vector( getNumBins(), std::array<uint8_t,3>( {0,0,0} ) ) );
	for( size_t bin = 0; bin < getNumBins(); ++bin )
		for( size_t frame = 0; frame < numFrames; ++frame )	
			{
			/* The choice of value and hue here is an aesthetic one
			 * I've tried to make both low and high frequency information relatively visible
			 * without oversaturating the highs
			 * Feel free to change the initial hue below if you don't like matrix green
			*/

			const double initialHueAngle = 110;

			const double normMag = std::abs(getBin( 0, frame, bin ).magnitude) / maxMag;
			const double value = std::pow( std::clamp( log2( double(bin) + 3.0 )*normMag, 0.0, 1.0 ), 0.4 );

			//bin deviation from bin center on range [-.5,.5]
			const double deviation = ( getBin(0, frame, bin).frequency / double(getBinWidth()) - double(bin)) / double( getOverlaps() );
			const int hueAngle = int( deviation * 1800.0 + initialHueAngle );
			const int wrappedHueAngle = hueAngle - int(360.0*floor( double(hueAngle) / 360.0 ));

			data[frame][bin] = HSVtoRGB( wrappedHueAngle, 1.0, value );
			}

	writeBMP( fileName, data );
	return *this;
	}

//========================================================================
//	Selection
//========================================================================

PVOC PVOC::getFrame( double time ) const
	{
	const double selectedFrame = std::clamp( timeToFrame_d( time ), 0.0, double( getNumFrames() - 1 ) );

	auto format = getFormat();
	format.numFrames = 1;
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t bin = 0; bin < out.getNumBins(); ++bin )
			out.setBin( channel, 0, bin, getBinInterpolated( channel, selectedFrame, bin ) );

	return out;
	}

PVOC PVOC::timeSelect( double length, Surface selector ) const
	{
	std::cout << "Time Select ... \n";
	auto format = getFormat();
	format.numFrames = timeToFrame( length );
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			{
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const double currentTimeSelection = selector( frameToTime( frame ), binToFrequency( bin ) );
				const double selectedFrame = std::clamp( timeToFrame_d( currentTimeSelection ), 0.0, double( getNumFrames() - 1 ) );
				out.setBin( channel, frame, bin,  
					getBinInterpolated( channel, selectedFrame, bin ) );
				}
			}

	return out;
	}

PVOC PVOC::holdAtTimes( const std::vector<double> & timesToHold  ) const
	{
	std::cout << "Hold Frames ... \n";
	PVOC out( getFormat() );

	std::vector<size_t> framesToHold( timesToHold.size() );

	std::transform( timesToHold.begin(), timesToHold.end(), framesToHold.begin(), [this]( double t ){ return timeToFrame( t ); } );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		bool holding = false;
		size_t currentFrameVecTestingPosition = 0;
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			if( frame == framesToHold[currentFrameVecTestingPosition] )
				{
				holding = true;
				++currentFrameVecTestingPosition;
				}

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				if( holding )
					out.setBin( channel, frame, bin,  
						getBin( channel, framesToHold[currentFrameVecTestingPosition-1], bin ) );
				else
					out.setBin( channel, frame, bin, 
						getBin( channel, frame, bin ) );
				}
			}
		}

	return out;
	}

//========================================================================
//	Resampling
//========================================================================

PVOC PVOC::modifyFrequency( Surface outFreqFunc, Interpolator interp ) const
	{
	std::cout << "Modify Frequency ... \n";
	PVOC out( getFormat() );
	out.clearBuffer();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			double prevOutBin = frequencyToBin_d( outFreqFunc( frameToTime( frame ), 0.0 ) );
			for( size_t bin = 1; bin < getNumBins(); ++bin )
				{
				const double outBin = frequencyToBin_d( outFreqFunc( frameToTime( frame ), binToFrequency( bin ) ) );
				const bool forward = outBin > prevOutBin;

				const long long prevOutBinRound = forward? std::ceil( prevOutBin ) : std::floor( prevOutBin );
				const long long outBinRound     = forward? std::ceil( outBin     ) : std::floor( outBin     );

				const size_t startBin = std::clamp( prevOutBinRound, 0ll, long long(out.getNumBins() - 1) );
				const size_t endBin   = std::clamp( outBinRound    , 0ll, long long(out.getNumBins() - 1) );

				for( size_t interpBin = startBin; interpBin != endBin; forward? ++interpBin : --interpBin )
					{
					if( interpBin < 0 || getNumBins() <= interpBin ) continue;

					const double mix = ( double( interpBin ) - prevOutBin ) / ( outBin - prevOutBin );
					const double interpMagnitude = interp( mix, getBin( channel, frame, bin-1 ).magnitude, getBin( channel, frame, bin ).magnitude );

					auto & outMF = out.getBin( channel, frame, interpBin );

					if( interpMagnitude > outMF.magnitude )
						{
						//Set frequency to whatever frequency is more dominant
						const double pureInterp = interp( mix, 0.0, 1.0 );

						if( getBin( channel, frame, bin-1 ).magnitude * ( 1.0 - pureInterp ) > 
						    getBin( channel, frame, bin-0 ).magnitude * ( 0.0 + pureInterp ) )
							outMF.frequency = outFreqFunc( frameToTime( frame ), getBin( channel, frame, bin-1 ).frequency );
						else
							outMF.frequency = outFreqFunc( frameToTime( frame ), getBin( channel, frame, bin-0 ).frequency );
						}
					outMF.magnitude += interpMagnitude;
					}

				prevOutBin = outBin;
				}
			}

	return out;
	}

PVOC PVOC::modifyTime( Surface outPosFunc, Interpolator interp ) const
	{
	std::cout << "Modify Time ... \n";

	//Find farthest output frame
	size_t lastOutputFrame = 0;
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const double outTime = outPosFunc( frameToTime( frame ), binToFrequency( bin ) );
				if( outTime <= 0.0 ) continue;
				lastOutputFrame = std::max( timeToFrame( outTime ), lastOutputFrame );
				}

	auto format = getFormat();
	format.numFrames = lastOutputFrame;
	PVOC out( format );
	out.clearBuffer();

	std::vector< double > prevOutPos( getNumBins() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Compute initial frame output positions
		for( size_t bin = 0; bin < getNumBins(); ++bin )
			prevOutPos[bin] = timeToFrame_d( outPosFunc( 0.0, binToFrequency( bin ) ) );

		for( size_t frame = 1; frame < getNumFrames(); ++frame ) // Start at 1 because for each frame we will access previous frame
			{
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const double outPos = timeToFrame_d( outPosFunc( frameToTime( frame ), binToFrequency( bin ) ) );
				const bool forward = outPos > prevOutPos[bin];

				const long long prevOutPosRound = forward? std::ceil( prevOutPos[bin] ) : std::floor( prevOutPos[bin] );
				const long long outPosRound     = forward? std::ceil( outPos          ) : std::floor( outPos          );

				const size_t startFrame = std::clamp( prevOutPosRound, 0ll, long long(out.getNumFrames() - 1) );
				const size_t endFrame   = std::clamp( outPosRound    , 0ll, long long(out.getNumFrames() - 1) );

				for( size_t outFrame = startFrame; outFrame != endFrame; forward? ++outFrame : --outFrame )
					{
					const double mix = ( double( outFrame ) - prevOutPos[bin] ) / ( outPos - prevOutPos[bin] );
					const double interpMagnitude = interp( mix, getBin( channel, frame - 1, bin ).magnitude, getBin( channel, frame, bin ).magnitude );

					auto & outMF = out.getBin( channel, outFrame, bin );
					if( interpMagnitude > outMF.magnitude )
						outMF.frequency = interp( mix, getBin( channel, frame - 1, bin ).frequency, getBin( channel, frame, bin ).frequency );
					outMF.magnitude += interpMagnitude;
					}

				prevOutPos[bin] = outPos;
				}
			}
		}

	return out;
	}

PVOC PVOC::repitch( RealFunc factor, Interpolator interp ) const
	{
	return modifyFrequency( [&factor]( double t, double f ){ return factor(t)*f; }, interp );
	}

PVOC PVOC::stretch( RealFunc factor, Interpolator interp ) const
	{
	return modifyTime( [&factor]( double t, double f ){ return factor(t)*t; }, interp );
	}

PVOC PVOC::interpolate_spline( RealFunc interpolation ) const
	{
	std::cout << "Interpolate_spline ... \n";
	const auto safeInterpolation = [&interpolation, this]( size_t f )
		{
		return std::max( size_t(interpolation( frameToTime( f ) ) ), 1ull );
		};
	
	//Set up output and spline x coordinates
	auto format = getFormat();
	format.numFrames = 0;
	std::vector<double> Xs( getNumFrames() );
	for( size_t frame = 0; frame < getNumFrames()-1; ++frame )
		{
		Xs[frame] = format.numFrames;
		format.numFrames += safeInterpolation( frame );
		}
	Xs[getNumFrames()-1] = format.numFrames;
	PVOC out( format );

	//Allocate Y coordinate vectors
	std::vector<double> magnitudeYs( Xs.size() );
	std::vector<double> frequencyYs( Xs.size() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Channel -> bin is correct here, despite non-sequential access to pvoc buffer
		for( size_t bin = 0; bin < getNumBins(); ++bin )
			{
			//For each x coordinate get corresponding frequency and magnitude
			for( size_t frame = 0; frame < getNumFrames(); ++frame )
				{
				magnitudeYs[frame] = getBin( channel, frame, bin ).magnitude;
				frequencyYs[frame] = getBin( channel, frame, bin ).frequency;
				}

			//Do the splines
			tk::spline magnitudeSpline, frequencySpline;
			magnitudeSpline.set_points(Xs,magnitudeYs);
			frequencySpline.set_points(Xs,frequencyYs);

			//Evaluate splines into out
			for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
				{
				out.setBin( channel, frame, bin, 
					{ 
					magnitudeSpline(frame),
					frequencySpline(frame)
					} );
				}
			}
		}
			
	return out;
	}

PVOC::MFPair PVOC::getBinInterpolated( size_t channel, double frame, size_t bin, Interpolator interpolator ) const
	{
	//Subtract epsilon for correct final frame rounding
	frame = std::clamp( frame, 0.0, double( getNumFrames() - 1 ) - 0.00001 );

	const size_t lowFrame = floor( frame );
	const double rem = frame - double( lowFrame );
	const auto l = getBin( channel, lowFrame + 0, bin );
	const auto h = getBin( channel, lowFrame + 1, bin );

	return  { 
			interpolator( rem, l.magnitude, h.magnitude ),
			interpolator( rem, l.frequency, h.frequency )
			};
	}

PVOC::MFPair PVOC::getBinInterpolated( size_t channel, size_t frame, double bin, Interpolator interpolator ) const
	{
	//Subtract epsilon for correct final frame rounding
	bin = std::clamp( bin, 0.0, double( getNumBins() - 1 ) - 0.00001 );

	const size_t lowBin = floor( bin );
	const double rem = bin - double( lowBin );
	const auto l = getBin( channel, frame, lowBin + 0 );
	const auto h = getBin( channel, frame, lowBin + 1 );

	return  { 
			interpolator( rem, l.magnitude, h.magnitude ),
			interpolator( rem, l.frequency, h.frequency )
			};
	}

PVOC PVOC::desample( RealFunc factor, Interpolator interpolator ) const
	{
	return stretch( [&factor]( double t ){ return 1.0 / factor(t); }, interpolator )
		  .stretch( factor,                                           interpolator );
	}

PVOC PVOC::timeExtrapolate( double startTime, double endTime, double extrapolationTime, Interpolator interpolator ) const
	{
	//Linear extrapolation via lines matching original PVOC at x_1 and x_2
	const size_t startFrame = timeToFrame( startTime		 ); //x_1
	const size_t endFrame	= timeToFrame( endTime			 ); //x_2
	const size_t extFrames	= timeToFrame( extrapolationTime ); //Extrapolated frames to generate

	auto format = getFormat();
	format.numFrames = endFrame + extFrames;
	PVOC out( format );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Copy up to start frame
		for( size_t frame = 0; frame < startFrame; ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				out.setBin( channel, frame, bin, getBin( channel, frame, bin ) );
			
		//Interpolate and continue beyond endFrame to extrapolate
		for( size_t frame = startFrame; frame < out.getNumFrames(); ++frame )
			{
			const double interpolation = double( frame - startFrame ) / double( endFrame - startFrame ); 

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const auto leftBin  =  getBin( channel, startFrame, bin );
				const auto rightBin =  getBin( channel, endFrame, bin );

				out.setBin( channel, frame, bin, 
					{ 
					std::abs( interpolator( interpolation, leftBin.magnitude, rightBin.magnitude ) ), 
					interpolator( interpolation, leftBin.frequency, rightBin.frequency ) 
					} );
				}
			}
		}
	return out;
	}

//========================================================================
// Combinations
//========================================================================

PVOC PVOC::replaceAmplitudes( const PVOC & ampSource, RealFunc amount ) const
	{
	PVOC out( getFormat() );

	if( ampSource.getNumChannels() < getNumChannels()
	 || ampSource.getNumFrames()   < getNumFrames() 
	 || ampSource.getNumBins()	   < getNumBins() )
		out.clearBuffer();

	const size_t numChannels = std::min( ampSource.getNumChannels(), getNumChannels() );
	const size_t numFrames	 = std::min( ampSource.getNumFrames(),   getNumFrames()   );
	const size_t numBins     = std::min( ampSource.getNumBins(),     getNumBins()     );

	for( size_t channel = 0; channel < numChannels; ++channel )
		for( size_t frame = 0; frame < numFrames; ++frame )
			{
			const double currentAmount = amount( frameToTime( frame ) );
			for( size_t bin = 0; bin < numBins; ++bin )
				{
				const auto currentBin = getBin( channel, frame, bin );
				out.setBin( channel, frame, bin,
					{
					ampSource.getBin( channel, frame, bin ).magnitude * currentAmount + currentBin.magnitude * ( 1.0 - currentAmount ),
					currentBin.frequency
					} );
				}
			}

	return out;
	}

PVOC PVOC::subtractAmplitudes( const PVOC & other, RealFunc amount ) const
	{
	PVOC out( *this );

	const size_t numChannels = std::min( other.getNumChannels(), getNumChannels() );
	const size_t numFrames	 = std::min( other.getNumFrames(),   getNumFrames()   );
	const size_t numBins     = std::min( other.getNumBins(),     getNumBins()     );

	for( size_t channel = 0; channel < numChannels; ++channel )
		for( size_t frame = 0; frame < numFrames; ++frame )
			{
			const double currentAmount = amount( frameToTime( frame ) );
			for( size_t bin = 0; bin < numBins; ++bin )
				out.getBin( channel, frame, bin ).magnitude -= other.getBin( channel, frame, bin ).magnitude * currentAmount;
			}

	return out;
	}

//========================================================================
// Uncategorized
//========================================================================

PVOC PVOC::perturb( RealFunc magAmount, RealFunc frqAmount, 
					Perturber magPerturber, Perturber frqPerturber ) const
	{
	PVOC out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			const double currentMagAmount = magAmount( frameToTime( frame ) );
			const double currentFrqAmount = frqAmount( frameToTime( frame ) );

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const auto currentBin = getBin( channel, frame, bin );
				out.setBin( channel, frame, bin,
					{
					magPerturber( currentBin.magnitude, currentMagAmount / log2( 1 + bin ) ),
					frqPerturber( currentBin.frequency, currentFrqAmount )
					} );
				}
			}

	return out;
	}

PVOC predicateNLoudestPartials( const PVOC & me, RealFunc numBins, std::function< bool (size_t, size_t) > predicate )
	{
	PVOC out( me.getFormat() );

	auto safeNumBins = [&numBins, me]( size_t frame )
		{
		return size_t( std::clamp( numBins( me.frameToTime( frame ) ), 0.0, double( me.getNumBins() ) ) );
		};

	typedef std::pair<size_t, double> indexVolume;
	std::vector<indexVolume> indexAndVolumes( me.getNumBins() );

	for( size_t channel = 0; channel < me.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < me.getNumFrames(); ++frame )
			{
			for( size_t bin = 0; bin < me.getNumBins(); ++bin )
				{
				indexAndVolumes[bin].first = bin;
				indexAndVolumes[bin].second = me.getBin( channel, frame, bin ).magnitude;
				}
			std::sort( indexAndVolumes.begin(), indexAndVolumes.end(), []( indexVolume& a, indexVolume& b ){ return abs(a.second) > abs(b.second); } );
			for( size_t bin = 0; bin < me.getNumBins(); ++bin )
				{
				const size_t actualBin = indexAndVolumes[bin].first;
				if( predicate( bin, safeNumBins( frame ) ) )
					out.setBin( channel, frame, actualBin, me.getBin( channel, frame, actualBin ) );
				else
					out.setBin( channel, frame, actualBin, { 0.0, me.getBin( channel, frame, actualBin ).frequency } );
				}
			}

	return out;
	}

PVOC PVOC::retainNLoudestPartials( RealFunc numBins ) const
	{
	return predicateNLoudestPartials( *this, numBins, []( size_t a, size_t b ){ return a < b; } );
	}

PVOC PVOC::removeNLoudestPartials( RealFunc numBins ) const
	{
	return predicateNLoudestPartials( *this, numBins, []( size_t a, size_t b ){ return a >= b; } );
	}

PVOC PVOC::resonate( double length, Surface decay ) const
	{
	length = std::max( length, 0.0 );
	auto format = getFormat();
	format.numFrames = std::max( timeToFrame( length ), 1ull );
	PVOC out( format );

	//Copy first frame into out
	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t bin = 0; bin < out.getNumBins(); ++bin )
			out.setBin( channel, 0, bin, getBin( channel, 0, bin ) );
	
	const double secondsPerFrame = 1.0 / timeToFrame( 1.0 );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 1; frame < out.getNumFrames(); ++frame )
			for( size_t bin = 0; bin < out.getNumBins(); ++bin )
				{
				const double currentDecay = pow( decay( frameToTime( frame ), binToFrequency( bin ) ), secondsPerFrame );
				const double decayedAmp = out.getBin( channel, frame - 1, bin ).magnitude * currentDecay;
				if( frame < getNumFrames() && getBin( channel, frame, bin ).magnitude > decayedAmp )
					out.setBin( channel, frame, bin, getBin( channel, frame, bin ) );
				else
					out.setBin( channel, frame, bin, { decayedAmp, out.getBin( channel, frame - 1, bin ).frequency } );
				}

	return out;
	}

} //End namespace xcdp