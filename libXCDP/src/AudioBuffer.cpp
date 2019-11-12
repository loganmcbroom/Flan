#include "AudioBuffer.h"

#include <iostream>
#include <algorithm>

#include <sndfile.h>

using namespace xcdp;

AudioBuffer::AudioBuffer( const std::string & filePath ) 
	{
	load( filePath );
	}
xcdp::AudioBuffer::AudioBuffer( const Format & other )
	: format( other )
	, buffer( getNumChannels() * getNumSamples() )
	{
	}

//======================================================
//	I/O
//======================================================
bool AudioBuffer::load( const std::string & filePath ) 
	{
	//Open file and check validity, save the sample rate
	SF_INFO info;
	SNDFILE * file = sf_open( filePath.data(), SFM_READ, &info ); 
	if( file == nullptr )
		{
		std::cout << filePath << " could not be opened.\n";
		return false;
		}

	//Copy file info into format
	format.sampleRate = info.samplerate;
	format.numChannels = info.channels;
	format.numSamples = info.frames;

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

	buffer.resize( info.channels * info.frames );
	for( size_t channel = 0; channel < info.channels; ++channel )
		{
		for( size_t frame = 0; frame < size_t(info.frames); ++frame )
			setSample( channel, frame, interleavedBuffer[ frame * info.channels + channel ] );
		}
	
	return true;
	}

bool AudioBuffer::save( const std::string & filePath ) const 
	{
	//Check that nothing silly is going on with the file formatting
	SF_INFO info = {};
	info.channels	= int( getNumChannels() );
	info.frames		= int( getNumSamples()	);
	info.samplerate = int( getSampleRate() );
	info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;
	if( !sf_format_check( &info ) )
		{
		std::cout << "Sound file formatting invalid while attempting to save to " << filePath << ",\n";
		printSummary();
		return false;
		}

	//Create a temporary buffer for interleaved data and copy the buffer in
	std::vector< double > interleavedBuffer( getNumSamples() * getNumChannels() );
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumSamples(); ++frame )
			interleavedBuffer[ frame * info.channels + channel ] = getSample( channel, frame );

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

void AudioBuffer::printSummary() const
	{
	std::cout << "\n===================================================\n";
	std::cout << "Channels: " << getNumChannels() << " Frames: " << getNumSamples() << " Sample Rate: " << getSampleRate() << "\n";
	std::cout << "===================================================\n\n";
	}

//======================================================
//	Getters
//======================================================
double AudioBuffer::getSample( size_t channel, size_t frame ) const 
	{
	//return buffer[channel][frame];
	return buffer[getPos( channel, frame )];
	}

AudioBuffer::Format xcdp::AudioBuffer::getFormat() const
	{
	return format;
	}

size_t AudioBuffer::getNumChannels() const
	{
	//return buffer.size();
	return format.numChannels;
	}

size_t AudioBuffer::getNumSamples() const 
	{
	//return buffer[0].size();
	return format.numSamples;
	}

size_t AudioBuffer::getSampleRate() const 
	{
	return format.sampleRate;
	}

double AudioBuffer::getTimeOfFrame( size_t sample ) const
	{
	return double( sample ) / double( getSampleRate() );
	}

double AudioBuffer::getMaxSampleMagnitude() const
	{
	return std::abs(*std::max_element( buffer.begin(), buffer.end(), []( double a, double b)
		{
		return std::abs( a ) < std::abs( b );
		}));
	}

//======================================================
//	Setters
//======================================================
void AudioBuffer::setSample( size_t channel, size_t frame, double sample ) 
	{
	//buffer[channel][frame] = sample;
	buffer[getPos( channel, frame )] = sample;
	}

/*
void AudioBuffer::setNumChannels( size_t newNumChannels ) 
	{
	//const size_t prevNumChannels = getNumChannels();
	//buffer.resize( newNumChannels );
	//if( newNumChannels > prevNumChannels )
	//	for( size_t channel = prevNumChannels; channel < newNumChannels; ++channel )
	//		std::fill( buffer[channel].begin(), buffer[channel].end(), 0 );
	setBufferSize( newNumChannels, numFrames );
	}

void AudioBuffer::setNumFrames( size_t newNumFrames ) 
	{
	//Expand or contract, if expanding fill expansion with 0
	//const size_t prevNumFrames = getNumSamples();
	//for( auto & channel : buffer )
	//	{
	//	channel.resize( newNumFrames );
	//	if( newNumFrames > prevNumFrames )
	//		std::fill( channel.begin() + prevNumFrames, channel.end(), 0.0f );
	//	}
	setBufferSize( numChannels, newNumFrames );
	}

void AudioBuffer::setBufferSize( size_t newNumChannels, size_t newNumFrames )
	{
	//setNumChannels( newNumChannels );
	//setNumFrames( newNumFrames );
	numChannels = newNumChannels;
	numFrames = newNumFrames;
	buffer.resize( numChannels * numFrames );
	}

void xcdp::AudioBuffer::setBufferSize( const AudioBuffer& other )
	{
	setBufferSize( other.getNumChannels(), other.getNumSamples() );
	}

void AudioBuffer::setSampleRate( size_t newSampleRate ) 
	{
	sampleRate = newSampleRate;
	}

void AudioBuffer::copyFormat( const AudioBuffer & other )
	{
	setSampleRate( other.getSampleRate() );
	}
*/

void xcdp::AudioBuffer::clearBuffer()
	{
	std::fill( buffer.begin(), buffer.end(), 0 );
	}

size_t xcdp::AudioBuffer::getPos( size_t c, size_t f ) const
	{
	return c * getNumSamples() + f;
	}
