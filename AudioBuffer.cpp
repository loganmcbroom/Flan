#include "Audio.h"

#include <iostream>
#include <algorithm>
#include <complex>

#include <sndfile.h>
#include <fftw3.h>

#include "PVOC.h"
#include "WindowFunctions.h"

using namespace xcdp;

Audio::Audio() 
	: sampleRate( 48000 )
	, buffer(1)
	{}
Audio::Audio( const std::string & filePath ) 
	: sampleRate( 48000 )
	{
	load( filePath );
	}
//Audio::Audio( const PVOC & spectra )
//	: Audio( spectra.getAudio() )
//	{}

//const Audio::ChannelProxy Audio::operator[]( size_t channel ) const
//	{
//	return Audio::ChannelProxy( *this, channel );
//	}

//======================================================
//	I/O
//======================================================
bool Audio::load( const std::string & filePath ) 
	{
	//Open file and check validity, save the sample rate
	SF_INFO info;
	SNDFILE * file = sf_open( filePath.data(), SFM_READ, &info ); 
	if( file == nullptr )
		{
		std::cout << filePath << " could not be opened.\n";
		return false;
		}
	sampleRate = info.samplerate;

	//Create temporary buffer for interleaved data in file
	std::vector< double > interleavedBuffer( info.frames * info.channels );

	//Read the entire file into that buffer, then close the file
	sf_readf_double( file, interleavedBuffer.data(), info.frames );
	if( sf_close( file ) != 0 )
		{
		std::cout << "Error closing " << filePath << ".\n";
		return false;
		}

	//Convert interleaved data into the buffer
	buffer.resize( info.channels );
	for( size_t channel = 0; channel < info.channels; ++channel )
		{
		buffer[channel].resize( info.frames );
		for( size_t frame = 0; frame < size_t(info.frames); ++frame )
			buffer[channel][frame] = interleavedBuffer[ frame * info.channels + channel ];
		}
	
	return true;
	}

bool Audio::save( const std::string & filePath ) 
	{
	//Check that nothing silly is going on with the file formatting
	SF_INFO info = {};
	info.channels	= int( getNumChannels() );
	info.frames		= int( getNumFrames()	);
	info.samplerate = int( getSampleRate() );
	info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;
	if( !sf_format_check( &info ) )
		{
		std::cout << "Sound file formatting invalid while attempting to save to " << filePath << ",\n";
		printSummary();
		return false;
		}

	//Create a temporary buffer for interleaved data and copy the buffer in
	std::vector< double > interleavedBuffer( info.frames * info.channels );
	for( size_t channel = 0; channel < size_t(info.channels); ++channel )
		for( size_t frame = 0; frame < size_t(info.frames); ++frame )
			interleavedBuffer[ frame * info.channels + channel ] = buffer[channel][frame];

	//Clip all samples in the interleaved buffer
	std::for_each( interleavedBuffer.begin(), interleavedBuffer.end(), []( double & s )
		{
		s = std::clamp( s, -1.0, 1.0 );
		});

	//Open the file and write in the interleaved buffer
	SF_INFO outInfo = info;
	SNDFILE * file = sf_open( filePath.data(), SFM_WRITE, &outInfo );
	if( file == nullptr )
		{
		std::cout << filePath << " could not be opened for saving.\n";
		return false;
		}
	if( sf_writef_double( file, interleavedBuffer.data(), info.frames ) != info.frames )	
		{
		std::cout << "Error writing data into " << filePath << ".\n";
		return false;
		}
	sf_close( file );
	
	return true;
	}

void Audio::printSummary() const
	{
	std::cout << "\n===================================================\n";
	std::cout << "Channels: " << getNumChannels() << " Frames: " << getNumFrames() << " Sample Rate: " << getSampleRate() << "\n";
	std::cout << "===================================================\n\n";
	}

//======================================================
//	Getters
//======================================================
double Audio::getSample( size_t channel, size_t frame ) const 
	{
	return buffer[channel][frame];
	}

size_t Audio::getNumChannels() const 
	{
	return buffer.size();
	}

size_t Audio::getNumFrames() const 
	{
	return buffer[0].size();
	}

size_t Audio::getSampleRate() const 
	{
	return sampleRate;
	}

double Audio::getTimeOfFrame( size_t sample ) const
	{
	return double( sample ) / double( getSampleRate() );
	}

// TODO?:
//		Zero stuff input buffer???
// pvoc_anal - serious shit for serious sound playas. This is a phase vocoder, which is basically a short time fft with a special sauce as
//    can be read about in the big text block below, or online here http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/
PVOC Audio::getPVOC( const size_t frameSize, const size_t overlaps ) const
	{
	std::cout << "Getting PVOC ... ";
	//Figure out the hop size and number of hops
	const size_t hopSize = frameSize / overlaps;
	const size_t numHops = size_t( ceil( getNumFrames() / hopSize ) );

	PVOC out;
	out.sampleRate = getSampleRate();
	out.overlaps = overlaps;
	out.setBufferSize( getNumChannels(), numHops, frameSize );

	//Ggenerate window
	const std::vector<double> window = window::Hann( frameSize );

	//Allocate fftw input buffer, output buffers, and plans. These can be reused a bunch of times, hence FFTW_MEASURE.
	//	Why the heck do we need two outs, two plans? I'm glad you asked! The phase vocoder needs information from the previous
	//	FFT frame to correctly predict the true frequency of partials. We will alternate output buffers (hence plans)
	//	for each fft, using the output of the previously used buffer to get that info.
	//  fftOuts[0] will be the first active plan, so we need to fill fftOuts[1] with zeros
	std::complex<double> *fftIn  = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * frameSize );
	std::complex<double> *fftOuts[2];
	fftOuts[0] = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * frameSize );
	fftOuts[1] = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * frameSize );
	fftw_plan plans[2]; 
	plans[0] = fftw_plan_dft_1d( int( frameSize ), (fftw_complex*) fftIn, (fftw_complex*) fftOuts[0], FFTW_FORWARD, FFTW_MEASURE );
	plans[1] = fftw_plan_dft_1d( int(frameSize ) , (fftw_complex*) fftIn, (fftw_complex*) fftOuts[1], FFTW_FORWARD, FFTW_MEASURE );

	//Explicit conversions outside of loops for optimization
	const size_t numFrames = getNumFrames();
	const double binWidthD = double(getSampleRate()) / double( frameSize );

	//For each channel, do the whole dang thing
	for( size_t channel = 0; channel < getNumChannels(); ++ channel )
		{
		//Zero filling fftOuts[1] as mentioned above
		std::fill( fftOuts[1], fftOuts[1] + frameSize, std::complex<double>(0,0) );

		//For each hop, fft and save into buffer
		for( size_t hop = 0; hop < numHops; ++hop )
			{
			const size_t activePlan = hop%2;
			const size_t prevPlan = (hop+1)%2;

			//Fill fft input buffer with windowed signal, weird for condition protects from bad in buffer access in the last few hops
			for( size_t relativeFrame = 0; relativeFrame < frameSize; ++relativeFrame )
				{
				size_t actualFrame = relativeFrame + hopSize * hop;
				if( actualFrame < numFrames )
					fftIn[relativeFrame] = { buffer[channel][actualFrame] * window[relativeFrame], 0 };
				else
					fftIn[relativeFrame] = {0,0};
				}
	
			//FFT that shit
			fftw_execute( plans[activePlan] );

			//Get the data!
			for( size_t bin = 0; bin < frameSize; ++bin )
				{
				//Magnitude is just magnitude
				out.buffer[channel][hop][bin].magnitude = abs( fftOuts[activePlan][bin] );

				//Phase is not just phase. 
				///Arg of each output is being calculated twice, any way to avoid that?
				///This method assumes the true frequency is locally constant, it would be pretty badass
				///  to model it with something else.
				//What's going on here?
				//	First, get the bin phase for the current and previous frame. Find the difference.
				//	We are expecting the phase of the ideal bin wave to move around due to using overlapping windows, so
				//    figure out what that expected movement is (or at least an equivalent angle).
				//  Next we get deltaPhase, this tells us how far from the expected phase shift our partial strayed.
				//	Up to this point we have been using angles potentially outside [-pi,pi], e.g. expectedPhaseDiff
				//    might be quite large when we really mean to talk about the equivalent angle in [-pi,pi]. We now
				//	  fix all of this by wrapping deltaPhase into [-pi,pi]
				//  Next we need how far our true frequency deviates from the ideal bin center frequency.
				//    Dividing by 2pi gives a deviation in [-1/2,1/2], and we multiply by overlaps as our actual bin
				//    deviation is overlaps times larger than was accounted for (recall we measured the phase
				//    difference between two overlapped frames).
				//  Finally the actual computed frequence is the bin center plus the deviation from that center
				//    times the width of a single bin.
				//  And it's as simple as that!
				
				const double phaseDiff = arg( fftOuts[activePlan][bin] ) - arg( fftOuts[prevPlan][bin] );
				const double expectedPhaseDiff = double(bin) * ( 2.0 * pi / double(overlaps) );
				const double deltaPhase = phaseDiff - expectedPhaseDiff;
				const double wrappedDeltaPhase = deltaPhase-(2.0*pi)*floor((deltaPhase+pi) / ( 2.0 * pi ) );
				const double binDeviation = double(overlaps) * wrappedDeltaPhase / ( 2.0 * pi );
				out.buffer[channel][hop][bin].frequency = ( double(bin) + binDeviation ) * binWidthD;

				//std::cout << wrappedDeltaPhase << "\n";
				}
			}
		}

	fftw_destroy_plan( plans[1] );
	fftw_destroy_plan( plans[0] );
	fftw_free( fftOuts[1] );
	fftw_free( fftOuts[0] );
	fftw_free( fftIn );

	std::cout << "Done\n";
	return out;
	}

const std::vector<std::vector<double>>& xcdp::Audio::getBuffer() const
	{
	return buffer;
	}

//======================================================
//	Setters
//======================================================
void Audio::setSample( size_t channel, size_t frame, double sample ) 
	{
	buffer[channel][frame] = sample;
	}

void Audio::setNumChannels( size_t newNumChannels ) 
	{
	const size_t prevNumChannels = getNumChannels();
	buffer.resize( newNumChannels );
	if( newNumChannels > prevNumChannels )
		for( size_t channel = prevNumChannels; channel < newNumChannels; ++channel )
			std::fill( buffer[channel].begin(), buffer[channel].end(), 0 );
	}

void Audio::setNumFrames( size_t newNumFrames ) 
	{
	//Expand or contract, if expanding fill expansion with 0
	const size_t prevNumFrames = getNumFrames();
	for( auto & channel : buffer )
		{
		channel.resize( newNumFrames );
		if( newNumFrames > prevNumFrames )
			std::fill( channel.begin() + prevNumFrames, channel.end(), 0.0f );
		}
	}

void Audio::setBufferSize( size_t newNumChannels, size_t newNumFrames )
	{
	setNumChannels( newNumChannels );
	setNumFrames( newNumFrames );
	}

void xcdp::Audio::setBufferSize( const Audio& other )
	{
	setBufferSize( other.getNumChannels(), other.getNumFrames() );
	}

void Audio::setSampleRate( size_t newSampleRate ) 
	{
	sampleRate = newSampleRate;
	}

void Audio::copyFormat( const Audio & other )
	{
	setSampleRate( other.getSampleRate() );
	}