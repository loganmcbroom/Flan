#include "PVOC.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <array>
#include <complex>

#include <fftw3.h>

#include "Audio.h"
#include "WindowFunctions.h"

namespace xcdp {

//PVOC::PVOC( const Audio & in )
//	{
//
//	}

//======================================================
//	Getters
//======================================================

size_t PVOC::getNumChannels() const
	{
	return buffer.size();
	}
size_t PVOC::getNumFrames() const
	{
	return buffer[0].size();
	}
size_t PVOC::getNumBins() const
	{
	return buffer[0][0].size();
	}
size_t PVOC::getFrameSize() const
	{
	return getNumBins();
	}
size_t PVOC::getSampleRate() const
	{
	return sampleRate;
	}
double PVOC::getBinWidth() const
	{
	return double( getSampleRate() ) / double( getFrameSize() );
	}
double PVOC::getMaxPartialMagnitude() const
	{
	double maxMagnitude = 0;
	for( auto & channel : buffer )
		for( auto & frame : channel  )
			for( auto bin : frame  )
				if( bin.magnitude > maxMagnitude )
					maxMagnitude = bin.magnitude;
	return maxMagnitude;
	}

//======================================================
//	Setters
//======================================================

void PVOC::setNumChannels( size_t newNumChannels )
	{
	buffer.resize( newNumChannels );
	}

void PVOC::setNumFrames( size_t newNumFrames )
	{
	for( auto & channel : buffer )
		channel.resize( newNumFrames );
	}

void PVOC::setNumBins( size_t newNumBins )
	{
	for( auto & channel : buffer )
		for( auto & frame : channel )
			frame.resize( newNumBins );
	}

void PVOC::setBufferSize( size_t newNumChannels, size_t newNumFrames, size_t newNumBins )
	{
	setNumChannels	( newNumChannels	);
	setNumFrames	( newNumFrames		);
	setNumBins		( newNumBins		);
	}

void PVOC::setBufferSize( const PVOC & other )
	{
	setNumChannels	( other.getNumChannels()	);
	setNumFrames	( other.getNumFrames()		);
	setNumBins		( other.getNumBins()		);
	}

void PVOC::setFrameSize( size_t newFrameSize )
	{
	setNumBins( newFrameSize );
	}

void PVOC::setSampleRate( size_t newSampleRate )
	{
	sampleRate = newSampleRate;
	}

void PVOC::copyFormat( const PVOC & other )
	{
	sampleRate = other.sampleRate;
	overlaps = other.overlaps;
	}

void PVOC::clearBuffer()
	{
	for( auto & channel : buffer )
		for( auto & frame : channel )
			for( auto & bin : frame )
				bin = {0,0};
	}

//======================================================
//	Conversions
//======================================================

//Phase Vocoder resynthesis. Comments are light, Audio::getPVOC has most of the info you'd need.
Audio PVOC::getAudio() const
	{
	std::cout << "Getting audio ... ";
	const size_t hopSize = getNumBins() / overlaps;
	const size_t frameSize = getNumBins();

	Audio out;
	out.setSampleRate( getSampleRate() );
	out.setBufferSize( getNumChannels(), getNumFrames() * hopSize );

	const std::vector<double> window = window::Hann( frameSize );

	std::complex<double> *fftIns[2];
	std::complex<double> *fftOut = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * frameSize );
	fftIns[0] = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * frameSize );
	fftIns[1] = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * frameSize );
	fftw_plan plans[2]; 
	plans[0] = fftw_plan_dft_1d( int(frameSize), (fftw_complex*) fftIns[0], (fftw_complex*) fftOut, FFTW_BACKWARD, FFTW_MEASURE );
	plans[1] = fftw_plan_dft_1d( int(frameSize), (fftw_complex*) fftIns[1], (fftw_complex*) fftOut, FFTW_BACKWARD, FFTW_MEASURE );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		std::fill( fftIns[1], fftIns[1] + frameSize, std::complex<double>(0,0) );

		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			const size_t activePlan = frame%2;
			const size_t prevPlan = (frame+1)%2;

			//reconstruct magnitude, phase, frequency info from magnitude, true frequency stored in buffer into ifft input
			//It's the same shit as in analyze but backwards so I'm not writing up every step
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const double magnitude = buffer[channel][frame][bin].magnitude;

				const double trueFreq = buffer[channel][frame][bin].frequency;
				const double binDeviation = trueFreq / double(getBinWidth()) - double(bin);
				const double deltaPhase = 2.0 * pi * binDeviation / double(overlaps);
				const double expectedPhaseDiff = double(bin) * ( 2.0 * pi / double(overlaps) );
				const double phaseDiff = deltaPhase + expectedPhaseDiff;
				const double phase = phaseDiff + arg( fftIns[prevPlan][bin] );
				
				fftIns[activePlan][bin] = std::polar( magnitude, phase );
				}

			//ifft 
			fftw_execute( plans[activePlan] );

			//Accumulate ifft output into audio buffer
			for( size_t relativeTimeFrame = 0; relativeTimeFrame < frameSize; ++relativeTimeFrame )
				{
				const size_t actualTimeFrame = frame * hopSize + relativeTimeFrame;
				if( actualTimeFrame < out.getNumFrames() )
					{
					//Don't forget to window, to normalize by 1/N due to fft scaling, and to divide by overlaps for oversampling
					///why the fuck are we windowing this, and I don't think the windows sum to a constant overlaps...
					const double normalizedOutputReal = fftOut[relativeTimeFrame].real() / ( getNumBins() * overlaps );
					out.buffer[channel][actualTimeFrame] += window[relativeTimeFrame] * normalizedOutputReal;
					}
				}
			}
		}

	fftw_destroy_plan( plans[1] );
	fftw_destroy_plan( plans[0] );
	fftw_free( fftIns[1] );
	fftw_free( fftIns[0] );
	fftw_free( fftOut );

	std::cout << "Done\n";
	return out;
	}


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

	return { char((Rs + m) * 255),  char((Gs + m) * 255),  char((Bs + m) * 255) };
	}
const PVOC &  PVOC::getSpectrograph( const std::string & fileName ) const
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

	for( size_t bin = 0; bin < getFrameSize()/2; ++bin )
		for( size_t frame = 0; frame < numFrames; ++frame )	
			{
			const double value = sqrt(sqrt(buffer[0][frame][bin].magnitude / maxMag ) );
			//bin deviation from bin center on range [-.5,.5]
			const double deviation = (buffer[0][frame][bin].frequency / double(getBinWidth()) - double(bin)) / double( overlaps );
			const int hue = int( deviation * 720.0 + 360 )  % 360;
			const std::array<char,3> rgb = HSVtoRGB( hue, 0.8, value );
			file.write( rgb.data(), 3 );
			}

	file.close();

	std::cout << "Done\n";
	return *this;
	}



} // End namespace xcdp