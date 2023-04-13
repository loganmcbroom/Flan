#include "flan/Audio/AudioBuffer.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <ranges>

using namespace std::ranges;

namespace flan {

AudioBuffer::AudioBuffer()
	: format()
	, buffer()
	{}

AudioBuffer::AudioBuffer( std::vector<float> && tempBuffer, Channel numChannels, SampleRate sr )
	: format()
	, buffer( tempBuffer )
	{
	format.numChannels = numChannels;
	format.numFrames = tempBuffer.size() / numChannels;
	format.sampleRate = sr;
	}

AudioBuffer::AudioBuffer( const Format & other )
	: format( other )
	, buffer( getNumChannels() * getNumFrames() )
	{}

AudioBuffer::AudioBuffer( const std::string & filename )
	: format()
	, buffer()
	{
	load( filename );
	}

AudioBuffer AudioBuffer::copy() const
	{
	AudioBuffer out;
	out.format = format;
	out.buffer = buffer;
	return out;
	}

bool AudioBuffer::isNull() const
	{
	return buffer.empty() || getSampleRate() == 0;
	}

//======================================================
//	I/O
//======================================================

#ifdef USE_SNDFILE // IO using libsndfile  

#include <sndfile.h>

bool AudioBuffer::load( const std::string & filePath ) 
	{
	//Open file and check validity, save the sample rate
	SF_INFO info;
	SNDFILE * file = sf_open( filePath.data(), SFM_READ, &info ); 
	if( file == nullptr )
		{
		std::cout << filePath + " could not be opened.\n";
		return false;
		}

	//Copy file info into format
	format.sampleRate = info.samplerate;
	format.numChannels = info.channels;
	format.numFrames = info.frames;
	*this = AudioBuffer( format );

	//Create temporary buffer for interleaved data in file, read data in, close the file
	std::vector<float> interleavedBuffer( info.frames * info.channels );
	sf_readf_float( file, interleavedBuffer.data(), info.frames );

	if( sf_close( file ) != 0 )
		{
		std::cout << std::string( "Error closing " ) + filePath + ".\n";
		return false;
		}

	//Convert interleaved data in
	for( Channel channel = 0; channel < info.channels; ++channel )
		for( Frame frame = 0; frame < Frame(info.frames); ++frame )
			setSample( channel, frame, interleavedBuffer[ frame * info.channels + channel ] );

	return true;
	}

bool AudioBuffer::save( const std::string & filePath, int format ) const 
	{
	if( format == -1 ) format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;

	//Check that nothing silly is going on with the file formatting
	SF_INFO info = {};
	info.channels	= int( getNumChannels() );
	info.frames		= int( getNumFrames()	);
	info.samplerate = int( getSampleRate() );
	info.format = format;
	if( !sf_format_check( &info ) )
		{
		std::cout << std::string( "Sound file formatting invalid while attempting to save to " ) + filePath + ",\n";
		printSummary();
		return false;
		}

	//Create a temporary buffer for interleaved data and copy the buffer in
	std::vector<float> interleavedBuffer( getNumFrames() * getNumChannels() );
	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame : iota_view( 0, getNumFrames() ) )
			interleavedBuffer[ frame * info.channels + channel ] = getSample( channel, frame );

	//Clip all samples in the interleaved buffer
	std::for_each( interleavedBuffer.begin(), interleavedBuffer.end(), []( float & s )
		{
		s = std::clamp( s, -1.0f, 1.0f );
		});

	//Open the file and write in the interleaved buffer
	SF_INFO outInfo = info;
	SNDFILE * file = sf_open( filePath.data(), SFM_WRITE, &outInfo );
	if( file == nullptr )
		{
		std::cout << filePath + " could not be opened for saving.\n";
		return false;
		}
	if( sf_writef_float( file, interleavedBuffer.data(), info.frames ) != info.frames )	
		{
		std::cout << std::string( "Error writing data into " ) + filePath + ".\n";
		return false;
		}
	sf_close( file );
	
	return true;
	}

#else // io using custom wave handler

bool AudioBuffer::load( const std::string & filePath ) 
	{
	auto bail = [filePath]( const std::string & s )
		{
		std::cout << "Couldn't load " << filePath << ": " << s << std::endl;
		return false;
		};

	std::ifstream file( filePath, std::ios::binary );
	if( !file ) return bail( "Error opening " + filePath + "." );

	uint16_t int16Buffer;
	uint32_t int32Buffer;
	char strBuffer[4];

	//Read RIFF chunk
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "RIFF", 4 ) != 0 ) return bail( "The input was not a WAVE file." );
	file.read( strBuffer, 4 ); //Don't care about file size
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "WAVE", 4 ) != 0 ) return bail( "The input was not a WAVE file." );

	//Read fmt subchunk
	AudioBuffer::Format format;
	uint16_t wavFormat;
	uint16_t blockAlign;
	uint16_t bitsPerSample;
	file.read( strBuffer, 4 );
	file.read( (char * ) &int32Buffer, 4 ); //Chunk size
	file.read( (char * ) &int16Buffer, 2 ); wavFormat = int16Buffer;
	file.read( (char * ) &int16Buffer, 2 ); format.numChannels = int16Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.sampleRate = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); // byte rate
	file.read( (char * ) &int16Buffer, 2 ); blockAlign = int16Buffer;
	file.read( (char * ) &int16Buffer, 2 ); bitsPerSample = int16Buffer;

	if( bitsPerSample > 32 ) return bail( "Can't load a wave with more than 32 bits per sample." );

	while( true ) //Read chunks until data chunk shows up
		{
		file.read( (char * ) &strBuffer, 4 );
		file.read( (char * ) &int32Buffer, 4 );
		if( strncmp( strBuffer, "data", 4 ) == 0 ) break;
		else //Go to next chunk and check for eof
			{
			file.ignore( int32Buffer - 8 );
			if( file.eof() )
				return bail( "File had no data chunk" );
			}
		} 

	format.numFrames = int32Buffer / blockAlign;

	*this = AudioBuffer( format );

	if( wavFormat == 1 )
		{
		if( bitsPerSample == 8 )
			{
			const double limit = std::pow( 2.0, 8.0 - 1.0 );

			//Read data buffer
			for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
				for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
					{
					uint8_t intSample = 0;
					file.read( ( (char*) &intSample ), 1 );
					
					setSample( channel, frame, double( intSample - 128 ) / limit );
					}
			}
        else if( bitsPerSample == 16 )
			{
			const double limit = std::pow( 2.0, 16.0 - 1.0 );

			//Read data buffer
			for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
				for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
					{
					int16_t intSample = 0;
					file.read( ( (char*) &intSample ), 2 );
			
					setSample( channel, frame, double( intSample ) / limit );
					}
			}
        else if( bitsPerSample == 24 )
			{
			const double limit = std::pow( 2.0, 24.0 - 1.0 );

           //Read data buffer
			for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
				for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
					{
					int32_t intSample = 0;
					file.read( ( (char*) &intSample ), 3 );

					//If sign bit is on, this is negative and its highest byte should be filled via 2s compliment
					if( intSample & 0x800000 ) intSample |= 0xFF000000;
			
					setSample( channel, frame, double( intSample ) / limit );
					}
			}
		else if( bitsPerSample == 32 )
			{
			const double limit = std::pow( 2.0, 32.0 - 1.0 );

           //Read data buffer
			for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
				for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
					{
					int32_t intSample = 0;
					file.read( ( (char*) &intSample ), 4 );
			
					setSample( channel, frame, double( intSample ) / limit );
					}
			}
		else
			{
			return bail( "Only 8, 16, 24, and 32 bit wave files are supported without libsndfile enabled." );
			}
		}
	else if( wavFormat == 3 ) // float data
		{
		for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
			for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
				{
				float f;
				file.read( (char *) &f, 4 );
				setSample( channel, frame, littleEndianToCurrent( f ) );
				}
		}
	else
		return bail( std::string( "Can't load a wave with format " ) + std::to_string( wavFormat ) + "." );
	
	return true;
	}

bool AudioBuffer::save( const std::string & filePath, int format ) const 
	{
	//Make sure path ends in wav
	if( filePath.substr( filePath.size() - 4, 4 ) != ".wav" )
		{
		std::cout << "Libsndfile is required for saving audio in formats other than wave"
			<< " (attempted to save a " << filePath.substr( filePath.size() - 4, 4 ) << " file ).\n";
		return false;
		}

	const uint16_t byteDepth = 3;
	const double limit = std::pow( 2.0, 8 * byteDepth - 1 );

	//Convert buffer to 24bit ints
	std::vector<uint8_t> byteBuffer( buffer->size() * byteDepth );
	for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
		for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
			{
			double clamped = std::clamp( getSample( channel, frame ), -1.0f, 1.0f );
			clamped *= limit;
			int32_t intSample = clamped;

			byteBuffer[ ( frame * getNumChannels() + channel ) * byteDepth + 0 ] = (uint8_t) ( intSample >> 0  ) & 0xFF;
			byteBuffer[ ( frame * getNumChannels() + channel ) * byteDepth + 1 ] = (uint8_t) ( intSample >> 8  ) & 0xFF;
			byteBuffer[ ( frame * getNumChannels() + channel ) * byteDepth + 2 ] = (uint8_t) ( intSample >> 16 ) & 0xFF;
			}

	//Write to file
	return writeRIFF( filePath, "WAVE", byteBuffer.data(), byteBuffer.size(),
		{
		uint16_t( 1 ),													// Formatting
		uint16_t( getNumChannels() ),									// Channel Count
		uint32_t( getSampleRate() ),									// Sample rate
		uint32_t( getSampleRate() * getNumChannels() * byteDepth ),		// Byte rate
		uint16_t( getNumChannels() * byteDepth ),						// Blockalign
		uint16_t( byteDepth * 8 )										// Bits per sample
		});
	}

#endif

void AudioBuffer::printSummary() const
	{
	std::cout << *this;
	}

//======================================================
//	Getters
//======================================================
auto AudioBuffer::getSample( Channel channel, Frame frame ) const -> Sample
	{
	return buffer[getBufferPos( channel, frame )];
	}

AudioBuffer::Format flan::AudioBuffer::getFormat() const
	{
	return format;
	}

auto AudioBuffer::getNumChannels() const -> Channel
	{
	return format.numChannels;
	}

auto AudioBuffer::getNumFrames() const -> Frame
	{
	return format.numFrames;
	}

SampleRate AudioBuffer::getSampleRate() const 
	{
	return format.sampleRate;
	}

float AudioBuffer::frameToTime() const
	{
	return 1.0f / timeToFrame();
	}

float AudioBuffer::timeToFrame() const
	{
	return float( getSampleRate() );
	}

auto AudioBuffer::getLength() const -> Time
	{
	return getNumFrames() * frameToTime();
	}

float AudioBuffer::getMaxSampleMagnitude( Time startTime, Time endTime ) const
	{
	if( endTime == 0 ) endTime = getLength();
	auto startIter = buffer.begin() + startTime * timeToFrame();
	auto endIter   = buffer.begin() + endTime   * timeToFrame();
	return std::abs(*std::max_element( startIter, endIter, []( float a, float b )
		{
		return std::abs( a ) < std::abs( b );
		}));
	}

//======================================================
//	Setters
//======================================================
void AudioBuffer::setSample( Channel channel, Frame frame, Sample sample ) 
	{
	buffer[getBufferPos( channel, frame )] = sample;
	}

auto AudioBuffer::getSample( Channel channel, Frame frame ) -> Sample &
	{
	return buffer[getBufferPos( channel, frame )];
	}

void flan::AudioBuffer::clearBuffer()
	{
	std::fill( buffer.begin(), buffer.end(), 0 );
	}

auto AudioBuffer::getSamplePointer( Channel channel, Frame frame ) -> Sample *
	{ 
	return buffer.data() + getBufferPos( channel, frame );
	}

auto AudioBuffer::getSamplePointer( Channel channel, Frame frame ) const -> const Sample *
	{ 
	return buffer.data() + getBufferPos( channel, frame );
	}

std::vector<Sample> & AudioBuffer::getBuffer() 
	{ 
	return buffer; 
	}
const std::vector<Sample> & AudioBuffer::getBuffer() const
	{
	return buffer; 
	}

std::vector<Sample>::const_iterator AudioBuffer::channelBegin( Channel channel ) const
	{
	return buffer.begin() + channel * getNumFrames();
	}

std::vector<Sample>::const_iterator AudioBuffer::channelEnd( Channel channel ) const
	{
	return buffer.begin() + ( channel + 1 ) * getNumFrames();
	}

size_t AudioBuffer::getBufferPos( Channel channel, Frame sample ) const
	{
	return channel * getNumFrames() + sample;
	}

//======================================================
//	Global
//======================================================
std::ostream & operator<<( std::ostream & os, const AudioBuffer & audio )
	{
	os << "\n=========================== Audio Info ==========================="
	   << "\nChannels:\t"		<< audio.getNumChannels() 
	   << "\nSamples:\t"		<< audio.getNumFrames() 
	   << "\nSample Rate:\t"	<< audio.getSampleRate()
	   << "\n==================================================================" 
	   << "\n\n";
	return os;
	}

} // End namespace flan
