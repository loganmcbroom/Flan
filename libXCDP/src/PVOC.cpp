#include "PVOC.h"

#include <algorithm>
#include <iostream>
#include <complex>
#include <array>
#include <fstream>
#include <random>

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
	hilite trace
	morph morph
	superaccu
	focus exag
	combine diff
*/

//	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
//		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
//			for( size_t bin = 0; bin < out.getNumBins(); ++bin )

namespace xcdp {

//========================================================================
//	Conversions
//========================================================================

//Phase Vocoder resynthesis. Comments are light, Audio::getPVOC has most of the info.
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

std::array<char,3> HSVtoRGB( int H, double S, double V ) 
	{
	const double C = S * V;
	const double X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	const double m = V - C;
	double Rs, Gs, Bs;

		 if( H >= 0   && H < 60  ) { Rs = C; Gs = X; Bs = 0; }
	else if( H >= 60  && H < 120 ) { Rs = X; Gs = C; Bs = 0; }
	else if( H >= 120 && H < 180 ) { Rs = 0; Gs = C; Bs = X; }
	else if( H >= 180 && H < 240 ) { Rs = 0; Gs = X; Bs = C; }
	else if( H >= 240 && H < 300 ) { Rs = X; Gs = 0; Bs = C; }
	else						   { Rs = C; Gs = 0; Bs = X; }

	return std::array<char,3>({ char((Rs + m) * 255),  char((Gs + m) * 255),  char((Bs + m) * 255) });
	}

const PVOC & PVOC::getSpectrograph( const std::string & fileName ) const
	{
	std::cout << "Getting Spectrograph ... \n";
	std::ofstream file( fileName, std::ios::binary );
	if( !file )
		{
		std::cout << "Error opening " << fileName << " to save spectrograph.\n";
		return *this;
		}

	//The image header
	unsigned char header[ 18 ] = { 0 };
	header[  2 ] = 1;  //Plain rgb
	header[ 12 ] = ( getNumFrames()		     )	& 0xFF;
	header[ 13 ] = ( getNumFrames()     >> 8 )	& 0xFF;
	header[ 14 ] = ( getFrameSize()/2		 )	& 0xFF;
	header[ 15 ] = ( getFrameSize()/2	>> 8 )	& 0xFF;
	header[ 16 ] = 24;  //Bits per pixel

	file.write( (const char *) header, 18 );

	size_t numFrames = getNumFrames();
	double maxMag = getMaxPartialMagnitude();

	//bin -> frame to match tga write order
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

			const std::array<char,3> rgb = HSVtoRGB( wrappedHueAngle, 1.0, value );
			file.write( rgb.data(), 3 );
			}

	file.close();
	return *this;
	}

//========================================================================
//	Selection
//========================================================================

PVOC PVOC::timeSelect( double length, RealFunc selector ) const
	{
	std::cout << "Time Select ... \n";
	auto format = getFormat();
	format.numFrames = timeToFrame( length );
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			{
			const size_t selectedFrame = std::clamp( timeToFrame( selector( frameToTime( frame ) ) ), 0ull, getNumFrames()-1 );
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				out.setBin( channel, frame, bin,  
					getBin( channel, selectedFrame, bin ) );
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

const PVOC::Interpolator PVOC::constantInterpolator = []( double mix, double a, double b ) 
	{ 
	return a;
	};

const PVOC::Interpolator PVOC::linearInterpolator = []( double mix, double a, double b ) 
	{ 
	return ( 1.0 - mix ) * a + mix * b;
	};

const PVOC::Interpolator PVOC::sineInterpolator = []( double mix, double a, double b ) 
	{ 
	return ( 1.0 - cos( pi * mix ) ) / 2.0 * ( b - a ) + a; 
	};

PVOC PVOC::decimate( RealFunc decimation ) const
	{
	std::cout << "Decimatate ... \n";
	auto safeDecimation = [&decimation, this]( size_t f )
		{
		return std::max( timeToFrame( decimation( frameToTime( f ) ) ), 1ull );
		};

	auto format = getFormat();
	format.numFrames = 0;
	for( size_t inFrame = 0; inFrame < getNumFrames(); inFrame += safeDecimation( inFrame ) )
		++format.numFrames;
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		{
		size_t inFrame = 0;
		for( size_t outFrame = 0; outFrame < out.getNumFrames(); ++outFrame )
			{
			for( size_t bin = 0; bin < out.getNumBins(); ++bin )
				out.setBin( channel, outFrame, bin, getBin( channel, inFrame, bin ) );
			inFrame += safeDecimation( inFrame );
			}
		}
	return out;
	}

PVOC PVOC::interpolate( RealFunc interpolation, PVOC::Interpolator interpolator ) const
	{
	std::cout << "Interpolate ... \n";
	const auto safeInterpolation = [&interpolation, this]( size_t frame )
		{
		return std::max( timeToFrame( interpolation( frameToTime( frame ) ) ), 1ull );
		};

	auto format = getFormat();
	format.numFrames = 0;
	for( size_t frame = 0; frame < getNumFrames(); ++frame )
		format.numFrames += safeInterpolation( frame );
	PVOC out( format );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Copy frames into output framed spaced interpolation apart
		size_t outFrame = 0;
		for( size_t inFrame = 0; inFrame < getNumFrames(); ++inFrame )
			{
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				out.setBin( channel, outFrame, bin, getBin( channel, inFrame, bin ) );
			outFrame += safeInterpolation(inFrame);
			}

		//For remaining out frames, interpolate
		size_t startFrame = 0, endFrame = 0;
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			{
			if( frame == endFrame ) 
				{
				startFrame = endFrame;
				endFrame   = std::min( out.getNumFrames()-1,  endFrame + safeInterpolation( endFrame ) );
				continue;
				}

			const double currentInterpPos = double( frame - startFrame ) / double( endFrame - startFrame );

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const auto leftBin  = out.getBin( channel, startFrame, bin );
				const auto rightBin = out.getBin( channel, endFrame,   bin );

				out.setBin( channel, frame, bin, 
					{ 
					interpolator( currentInterpPos, leftBin.magnitude, rightBin.magnitude ), 
					interpolator( currentInterpPos, leftBin.frequency, rightBin.frequency ) 
					} );
				}
			}
		}

	return out;
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

PVOC PVOC::resample( RealFunc decimation, RealFunc interpolation, PVOC::Interpolator interpolator ) const
	{
	return interpolate( interpolation, interpolator ).decimate( decimation );
	}

PVOC PVOC::desample( RealFunc decimation, RealFunc interpolation, PVOC::Interpolator interpolator ) const
	{
	return decimate( decimation ).interpolate( interpolation, interpolator );
	}

PVOC PVOC::timeExtrapolate( double startTime, double endTime, double extrapolationTime, PVOC::Interpolator interpolator ) const
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
					interpolator( interpolation, leftBin.magnitude, rightBin.magnitude ), 
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

//========================================================================
// Uncategorized
//========================================================================

PVOC PVOC::repitch( double factor ) const
	{
	std::cout << "Repitching ... \n";
	PVOC out( getFormat() );
	out.clearBuffer();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				size_t index = size_t( bin * factor );
				if( index < getNumBins() )
					{
					auto & outMF = out.getBin( channel, frame, index );
					auto & inMF  =     getBin( channel, frame, bin   );
					outMF.magnitude += inMF.magnitude;
					outMF.frequency =  inMF.frequency * factor;
					}
				}
	return out;
	}

const PVOC::Perturber PVOC::normalDistPerturber = []( double x, double amount )
	{
	///Seeding?
	static std::default_random_engine rng;
	static std::normal_distribution<double> dist( 0, 1.0 );

	return dist(rng) * amount + x;
	};

const PVOC::Perturber PVOC::identityPerturber = []( double x, double amount ){ return x; };

PVOC PVOC::perturb( RealFunc magAmount, RealFunc frqAmount, 
					PVOC::Perturber magPerturber, PVOC::Perturber frqPerturber ) const
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


} //End namespace xcdp