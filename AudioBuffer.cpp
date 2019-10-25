#include "AudioBuffer.h"

#include <iostream>
#include <algorithm>

using namespace xcdp;

AudioBuffer::AudioBuffer() 
	: info( {} )
	, buffer()
	{
	}
xcdp::AudioBuffer::AudioBuffer( const std::string & filePath ) 
	: info( {} )
	{
	load( filePath );
	}
AudioBuffer::~AudioBuffer() 
	{
	}
AudioBuffer::AudioBuffer( const AudioBuffer & other ) 
	: info( other.info )
	, buffer( other.buffer )
	{
	}
AudioBuffer::AudioBuffer( AudioBuffer && other ) noexcept 
	: info( other.info )
	, buffer( std::move( other.buffer ) )
	{
	}
AudioBuffer & AudioBuffer::operator=( const AudioBuffer & other ) 
	{
	info = other.info;
	buffer = other.buffer;
	return *this;
	}
AudioBuffer & AudioBuffer::operator=( AudioBuffer && other ) noexcept 
	{
	if( &other == this ) 
		return *this;

	info = other.info;
	buffer = std::move( other.buffer );
	return *this;
	}

//======================================================
//	I/O
//======================================================
bool AudioBuffer::load( const std::string & filePath ) 
	{
	SNDFILE * file = sf_open( filePath.data(), SFM_READ, &info ); 
	if( file == nullptr )
		{
		std::cout << filePath << " could not be opened.\n";
		return false;
		}
	buffer.resize( info.frames * info.channels );
	sf_readf_float( file, buffer.data(), info.frames );
	if( sf_close( file ) != 0 )
		{
		std::cout << "Error closing " << filePath << ".\n";
		return false;
		}
	sf_close( file );
	return true;
	}

bool AudioBuffer::save( const std::string & filePath ) 
	{
	if( !sf_format_check( &info ) )
		{
		std::cout << "Sound file formatting invalid while attempting to save to " << filePath << ",\n";
		printSummary();
		return false;
		}
	for( float & s : buffer )
		if( std::abs( s ) > 1 )
			s = std::clamp( s, -1.0f, 1.0f );
	SF_INFO outInfo = info;
	SNDFILE * file = sf_open( filePath.data(), SFM_WRITE, &outInfo );
	if( file == nullptr )
		{
		std::cout << filePath << " could not be opened for saving.\n";
		return false;
		}
	if( sf_writef_float( file, buffer.data(), info.frames ) != info.frames )	
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
	std::cout << "Channels: " << info.channels << " Frames: " << info.frames << "\n";
	std::cout << "Sample Rate: " << info.samplerate << " Format: " << info.format << "\n";
	std::cout << "===================================================\n\n";
	}

//======================================================
//	Getters
//======================================================
float AudioBuffer::getSample( int channel, int frame ) const 
	{
	return buffer[ frame * info.channels + channel ];
	}

int AudioBuffer::getNumChannels() const 
	{
	return info.channels;
	}

int AudioBuffer::getNumFrames() const 
	{
	return info.frames;
	}

int AudioBuffer::getSampleRate() const 
	{
	return info.samplerate;
	}
double AudioBuffer::getTimeOfFrame( int sample ) const
	{
	return double( sample ) / double( getSampleRate() );
	}

//======================================================
//	Setters
//======================================================
void AudioBuffer::setSample( int channel, int frame, float sample ) 
	{
	buffer[ frame * info.channels + channel ] = sample;
	}

void AudioBuffer::setNumChannels( int newNumChannels ) 
	{
	//Because data is held in a single buffer in frame major order, this has to do serious data manipulation
	///TODO: implement
	info.channels = newNumChannels;
	}

void AudioBuffer::setNumFrames( int newNumFrames ) 
	{
	//Expand or contract, if expanding fill expansion, set metadata
	buffer.resize( getNumChannels() * newNumFrames );
	if( newNumFrames > getNumFrames() )
		std::fill( buffer.begin() + getNumFrames(), buffer.end(), 0.0f );
	info.frames = newNumFrames;
	}

void AudioBuffer::setBufferSize( int newNumChannels, int newNumFrames )
	{
	setNumChannels( newNumChannels );
	setNumFrames( newNumFrames );
	}

void AudioBuffer::setSampleRate( int newSampleRate ) 
	{
	info.samplerate = newSampleRate;
	}

void AudioBuffer::copyFormat( const AudioBuffer & other )
	{
	setSampleRate( other.getSampleRate() );
	info.format = other.info.format;
	}