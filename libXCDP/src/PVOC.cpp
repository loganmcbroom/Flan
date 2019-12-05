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
				const double currentTimeSelection = selector( frameToTime( frame ), getBinFrequency( bin ) );
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

PVOC PVOC::interpolate( RealFunc interpolation, Interpolator interpolator ) const
	{
	std::cout << "Interpolate ... \n";
	const auto safeInterpolation = [&interpolation, this]( size_t frame )
		{
		return std::clamp( timeToFrame( interpolation( frameToTime( frame ) ) ), 1ull, 64ull );
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
			outFrame += safeInterpolation( inFrame );
			}

		//Copy final input frame with zeroed amplitudes into final output frame for fading out of that frame
		for( size_t bin = 0; bin < getNumBins(); ++bin )
			out.setBin( channel, out.getNumFrames()-1, bin, { 0, getBin( channel, getNumFrames()-1, bin ).frequency } );

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

PVOC PVOC::resample( RealFunc decimation, RealFunc interpolation, Interpolator interpolator ) const
	{
	return interpolate( interpolation, interpolator ).decimate( decimation );
	}

PVOC PVOC::desample( RealFunc decimation, RealFunc interpolation, Interpolator interpolator ) const
	{
	return decimate( decimation ).interpolate( interpolation, interpolator );
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

PVOC PVOC::repitch( Surface factor ) const
	{
	std::cout << "Repitching ... \n";
	PVOC out( getFormat() );
	out.clearBuffer();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const double currentFactor = std::max( 0.0, factor( frameToTime( frame ), getBinFrequency( bin ) ) );

				const double outBin_d = bin * currentFactor;
				const size_t outBin_lo = floor( outBin_d );
				const size_t outBin_hi = outBin_lo + 1;
				const double interp = outBin_d - double( outBin_lo );

				if( outBin_hi < out.getNumBins() )
					{
					const auto inMF = getBin( channel, frame, bin );
					auto & outMF_lo = out.getBin( channel, frame, outBin_lo  );
					auto & outMF_hi = out.getBin( channel, frame, outBin_hi );

					const double energy_lo = inMF.magnitude * ( 1.0 - interp );
					const double energy_hi = inMF.magnitude * ( 0.0 + interp );

					if( energy_lo > outMF_lo.magnitude ) outMF_lo.frequency = inMF.frequency * currentFactor;
					if( energy_hi > outMF_hi.magnitude ) outMF_hi.frequency = inMF.frequency * currentFactor;
					outMF_lo.magnitude += energy_lo;
					outMF_hi.magnitude += energy_hi;
					}
				}

	return out;
	}

PVOC PVOC::stretch( Surface factor ) const
	{
	//Find farthest output frame
	size_t lastOutputFrame = 0;
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				lastOutputFrame = std::max( size_t( frame * factor( frameToTime( frame ), getBinFrequency( bin ) ) ), lastOutputFrame );

	auto format = getFormat();
	format.numFrames = lastOutputFrame;
	PVOC out( format );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const double currentFactor = factor( frameToTime( frame ), getBinFrequency( bin ) );

				const double outFrame_d = frame * currentFactor;
				const size_t outFrame_lo = floor( outFrame_d );
				const size_t outFrame_hi = outFrame_lo + 1;
				const double interp = outFrame_d - double( outFrame_lo );

				if( outFrame_hi < out.getNumFrames() )
					{
					const auto inMF  = getBin( channel, frame, bin );
					auto & outMF_lo = out.getBin( channel, outFrame_lo, bin );
					auto & outMF_hi = out.getBin( channel, outFrame_hi, bin );

					const double energy_lo = inMF.magnitude * ( 1.0 - interp );
					const double energy_hi = inMF.magnitude * ( 0.0 + interp );
				
					outMF_lo.frequency = outMF_lo.frequency * outMF_lo.magnitude + inMF.frequency * energy_lo;
					outMF_hi.frequency = outMF_hi.frequency * outMF_hi.magnitude + inMF.frequency * energy_hi;

					outMF_lo.magnitude += energy_lo;
					outMF_hi.magnitude += energy_hi;

					if( outMF_lo.magnitude != 0.0 ) outMF_lo.frequency /= outMF_lo.magnitude;
					if( outMF_hi.magnitude != 0.0 ) outMF_hi.frequency /= outMF_hi.magnitude;
					}
				}

	return out;
	}

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

PVOC PVOC::convolve( const PVOC & other ) const
	{
	auto format = getFormat();
	format.numChannels = std::min( getNumChannels(), other.getNumChannels() );
	format.numFrames = getNumFrames() + other.getNumFrames();
	format.numBins = std::min( getNumBins(), other.getNumBins() );
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			{
			for( size_t bin = 0; bin < out.getNumBins(); ++bin )
				{
				MFPair convolutionAccumulator = {0,0};
				for( size_t i = 0; i < other.getNumFrames(); ++i )
					{
					if( frame - i < 0 || getNumFrames() <= frame - i ) continue;

					const auto thisBin  = getBin( channel, frame - i, bin );
					const auto otherBin = getBin( channel,	       i, bin );

					convolutionAccumulator.magnitude += thisBin.magnitude * otherBin.magnitude;
					convolutionAccumulator.frequency += thisBin.frequency * otherBin.frequency;
					}
				out.setBin( channel, frame, bin, convolutionAccumulator );
				}
			}

	return out;
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
				const double currentDecay = pow( decay( frameToTime( frame ), getBinFrequency( bin ) ), secondsPerFrame );
				const double decayedAmp = out.getBin( channel, frame - 1, bin ).magnitude * currentDecay;
				if( frame < getNumFrames() && getBin( channel, frame, bin ).magnitude > decayedAmp )
					out.setBin( channel, frame, bin, getBin( channel, frame, bin ) );
				else
					out.setBin( channel, frame, bin, { decayedAmp, out.getBin( channel, frame - 1, bin ).frequency } );
				}

	return out;
	}

} //End namespace xcdp