#include "PVOC.h"

#include <Windows.h>

#include <algorithm>
#include <iostream>
#include <complex>
#include <array>
#include <fstream>

#include <fftw3.h>

#include "WindowFunctions.h"
#include "Audio.h"

/** TODO:

	blur blur
	blur chorus
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

namespace xcdp {

//Phase Vocoder resynthesis. Comments are light, Audio::getPVOC has most of the info.
Audio PVOC::convertToAudio() const
	{
	std::cout << "Getting audio ... ";
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

	std::cout << "Done\n";
	return out;
	}

PVOC PVOC::timeSelect( double length, RealFunc selector ) const
	{
	std::cout << "Time Select ... ";
	auto format = getFormat();
	format.numFrames = timeToFrame( length );
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			{
			const size_t selectedFrame = std::clamp( timeToFrame( selector( frameToTime( frame ) ) ), 0ll, long long(getNumFrames()-1) );
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				out.setBin( channel, frame, bin,  
					getBin( channel, selectedFrame, bin ) );
				}
			}

	std::cout << "Done\n";
	return out;
	}

PVOC PVOC::holdAtTimes( const std::vector<double> & timesToHold  ) const
	{
	std::cout << "Hold Frames ... ";
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

	std::cout << "Done\n";
	return out;
	}

//should just call holdAtFrames
PVOC PVOC::constantTimeInterpolate( int factor ) const
	{
	std::cout << "Constant Time Interpolation ... ";
	PVOC out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			const size_t startFrame = size_t(floor( frame / factor )) * factor;

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				out.setBin( channel, frame, bin,  
					getBin( channel, startFrame, bin ) );

				}
			}

	std::cout << "Done\n";
	return out;
	}

PVOC PVOC::linearTimeInterpolate( int factor ) const
	{
	std::cout << "Linear Time Interpolation ... ";
	PVOC out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			const size_t startFrame = size_t(floor( frame / factor )) * factor;
			const size_t endFrame	  = min( getNumFrames()-1, startFrame + factor );
			const size_t clampedFactor = endFrame-startFrame;
			const double interpolation = double(frame % factor) / double(clampedFactor);

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				out.setBin( channel, frame, bin, 
					{ 
					  getBin( channel, startFrame, bin ).magnitude*(1.0-interpolation) 
					+ getBin( channel, endFrame  , bin ).magnitude*interpolation,
					  getBin( channel, startFrame, bin ).frequency*(1.0-interpolation) 
					+ getBin( channel, endFrame  , bin ).frequency*interpolation
					} );
				}
			}

	std::cout << "Done\n";
	return out;
	}

PVOC PVOC::repitch( double factor ) const
	{
	std::cout << "Repitching ... ";
	PVOC out( getFormat() );
	out.clearBuffer();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				size_t index = size_t( bin * factor );
				if( index < getNumBins() )
					{
					out.setBin( channel, frame, index,
						{
						out.getBin( channel, frame, index ).magnitude + getBin( channel, frame, bin ).magnitude,
						getBin( channel, frame, bin ).frequency * factor
						});
					}
				}

	std::cout << "Done\n";
	return out;
	}

//PVOC PVOC::gate( RealFunc cutoff ) const
//	{
//	SpectralBuffer out;
//	out.copyFormat( in );
//	out.setBufferSize( in );
//
//	for( int channel = 0; channel < in.getNumChannels(); ++channel )
//		for( int bin = 0; bin < out.getNumBins(); ++bin )
//			out[channel][bin] = abs(in[channel][bin]) > cutoff( in.getBinFreq( bin ) ) ? in[channel][bin] : 0;
//	
//
//	return out;
//	}
//PVOC PVOC::gate( double cutoff ) const
//	{
//	return gate( in, [cutoff]( double f ){ return cutoff; } );
//	}
//
//PVOC PVOC::invert( int lowerBound, int upperBound  ) const
//	{
//	lowerBound = std::clamp( lowerBound, 0, int( getNumBins() - 1 ) );
//	upperBound = std::clamp( upperBound, 0, int( getNumBins() - 1 ) );
//	upperBound = std::max( lowerBound, upperBound );
//
//	PVOC out( in );
//	for( int channel = 0; channel < out.getNumChannels(); ++channel )
//		std::reverse( out.buffer[channel].begin() + lowerBound, out.buffer[channel].begin() + upperBound );
//
//	return out;
//	}

/*
 * H(Hue): 0 - 360 degree (integer)
 * S(Saturation): 0 - 1.00 (double)
 * V(Value): 0 - 1.00 (double)
 */
std::array<char,3> HSVtoRGB( int H, double S, double V ) 
	{
	double C = S * V;
	double X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	double m = V - C;
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
	std::cout << "Getting Spectrograph ... ";
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
			const double normMag = getBin( 0, frame, bin ).magnitude / maxMag;
			const double value = std::pow( std::clamp( log( double(bin) + 3.0 )*normMag, 0.0, 1.0 ), 0.4 );

			//bin deviation from bin center on range [-.5,.5]
			const double deviation = ( getBin(0, frame, bin).frequency / double(getBinWidth()) - double(bin)) / double( getOverlaps() );
			const int hueAngle = int( deviation * 3000.0 + 110.0 );
			const int wrappedHueAngle = hueAngle - int(360.0*floor( double(hueAngle) / 360.0 ));

			const std::array<char,3> rgb = HSVtoRGB( wrappedHueAngle, 1.0, value );
			file.write( rgb.data(), 3 );
			}

	file.close();

	system( (std::string("start \"C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe\" ") + fileName).c_str() );

	std::cout << "Done\n";
	return *this;
	}

} //End namespace xcdp