#include "xcdp/PVOCBuffer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

#include "xcdp/Utility.h"

namespace xcdp {

PVOCBuffer::PVOCBuffer( const Format & other )
	: format( other )
	, buffer( std::make_shared<std::vector<MF>>( getNumChannels() * getNumFrames() * getNumBins() ))
	{}

PVOCBuffer::PVOCBuffer( const std::string & filename )
	: format()
	, buffer()
	{
	load( filename );
	}

bool PVOCBuffer::save( const std::string & filename ) const
	{
	const int byteDepth = 3;
	const double limit = std::pow( 2, 8 * byteDepth - 1 );
	const float windowSize_f = getWindowSize();
	const float maxFreq_f = getSampleRate();

	//Convert buffer to 24bit signed int representation
	std::vector<uint8_t> bytes( buffer->size() * 2 * byteDepth );
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				int32_t m_32 = double( std::clamp( getMF( channel, frame, bin ).m / windowSize_f, -1.0f, 1.0f ) ) * limit;
				int32_t f_32 = double( std::clamp( getMF( channel, frame, bin ).f / maxFreq_f,	  -1.0f, 1.0f ) ) * limit;

				size_t pos = getBufferPos( channel, frame, bin ) * 2 * byteDepth;

				bytes[ pos + 0 ] = (uint8_t) ( m_32 >> 0  ) & 0xFF;
				bytes[ pos + 1 ] = (uint8_t) ( m_32 >> 8  ) & 0xFF;
				bytes[ pos + 2 ] = (uint8_t) ( m_32 >> 16 ) & 0xFF;

				bytes[ pos + 3 ] = (uint8_t) ( f_32 >> 0 )  & 0xFF;
				bytes[ pos + 4 ] = (uint8_t) ( f_32 >> 8 )  & 0xFF;
				bytes[ pos + 5 ] = (uint8_t) ( f_32 >> 16 ) & 0xFF;
				}

	//Write entire buffer to file
	writeRIFF( filename, "PVOC", bytes.data(), bytes.size(), 
		{
		(uint16_t) 1,					//Formatting
		(uint16_t) getNumChannels(),	//Channel Count
		(uint32_t) getNumFrames(),		//Number of frames
		(uint32_t) getNumBins(),		//Number of bins per frame
		(uint32_t) getSampleRate(),		//Sample Rate, note this is the sample rate of the audio
		(uint32_t) getOverlaps(),		//Number of windows overlapping at any frame
		(uint32_t) 24,					//Bit depth, note that each bin contains two of this
		(uint16_t) 1					//Window type indicator. 1 = hann.
		});

	return true;
	}

bool PVOCBuffer::load( const std::string & filename )
	{
	auto bail = []( const std::string & s )
		{
		std::cout << s << std::endl;
		return false;
		};

	std::ifstream file( filename, std::ios::binary );
	if( !file ) return bail( "Error opening " + filename + " to load PVOC." );

	uint16_t int16Buffer;
	uint32_t int32Buffer;
	char strBuffer[4];

	//Read RIFF chunk
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "RIFF", 4 ) != 0 ) return bail( filename + " isn't a correctly formatted RIFF file.\n" );
	file.read( strBuffer, 4 ); //Don't care about file size
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "PVOC", 4 ) != 0 ) return bail( filename + " isn't a PVOC file.\n"  );

	//Read fmt subchunk
	PVOCBuffer::Format format;
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "fmt ", 4 ) != 0 ) return bail( filename + " isn't formatted correctly (\"fmt \" wasn't at the start of the format chunk).\n" );
	file.read( (char * ) &int32Buffer, 4 ); //Chunk size
	file.read( (char * ) &int16Buffer, 2 );  if( int16Buffer != 1 ) return bail( "Formatting must be 1 (signed int)." );
	file.read( (char * ) &int16Buffer, 2 ); format.numChannels = int16Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.numFrames = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.numBins = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.sampleRate = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.overlaps = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); if( int32Buffer != 24 ) return bail( "Bit depth must be 24." );
	file.read( (char * ) &int16Buffer, 2 ); if( int16Buffer != 1 ) return bail( "PVOC window must be 1 (Hann)." );
	*this = PVOCBuffer( format );

	//Read data subchunk
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "data", 4 ) != 0 ) return bail( filename + " isn't a correctly formatted PVOC file (\"data\" wasn't at the start of the data chunk).\n" );
	file.read( (char * ) &int32Buffer, 4 );

	
	//Read buffer
	const double limit = std::pow( 2, 23 );
	const float windowSize_f = getWindowSize();
	const float maxFreq_f = getSampleRate();

	auto getFloat = [limit, &file, windowSize_f]( float div )
		{
		int32_t i = 0;
		file.read( ( (char*) &i ), 3 );
		if( i & 0x800000 ) i |= 0xFF000000;
		return float( double( i ) / limit ) * div;
		};

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				setMF( channel, frame, bin, { getFloat( windowSize_f ), getFloat( maxFreq_f ) } );

	return true;		
	}

void PVOCBuffer::printSummary() const
	{
	std::cout << *this;
	}

//======================================================
//	Getters
//======================================================

PVOCBuffer::MF PVOCBuffer::getMF( size_t channel, size_t frame, size_t bin ) const
	{
	return (*buffer)[getBufferPos( channel, frame, bin )];
	}

size_t PVOCBuffer::getNumChannels() const
	{
	return format.numChannels;
	}
size_t PVOCBuffer::getNumFrames() const
	{
	return format.numFrames;
	}
size_t PVOCBuffer::getNumBins() const
	{
	return format.numBins;
	}
size_t PVOCBuffer::getWindowSize() const
	{
	return ( getNumBins() - 1 ) * 2;
	}
PVOCBuffer::Format PVOCBuffer::getFormat() const
	{
	return format;
	}
size_t PVOCBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}
size_t PVOCBuffer::getOverlaps() const
	{
	return format.overlaps;
	}

float PVOCBuffer::getLength() const
	{
	return getNumFrames() * frameToTime();
	}

float PVOCBuffer::getMaxPartialMagnitude() const
	{
	float maxMagnitude = 0;
	for( auto i : *buffer )
		{
		const float trueMag = std::abs( i.m );
		if( trueMag > maxMagnitude )
			maxMagnitude = trueMag;
		}
	return maxMagnitude;
	}


float PVOCBuffer::timeToFrame() const
	{
	return float(getSampleRate()) * float(getOverlaps()) / float(getWindowSize());
	}
float PVOCBuffer::frameToTime() const
	{
	return 1.0f / timeToFrame();
	}

float PVOCBuffer::frequencyToBin() const
	{
	return 1.0f / binToFrequency();
	}
float PVOCBuffer::binToFrequency() const
	{
	return float( getSampleRate() ) / float( getWindowSize() );
	}

//======================================================
//	Setters
//======================================================

void PVOCBuffer::setMF( size_t channel, size_t frame, size_t bin, MF value )
	{
	(*buffer)[getBufferPos( channel, frame, bin )] = value;
	}
PVOCBuffer::MF & PVOCBuffer::getMF( size_t channel, size_t frame, size_t bin )
	{
	return (*buffer)[getBufferPos( channel, frame, bin )];
	}

void PVOCBuffer::clearBuffer()
	{
	std::fill( buffer->begin(), buffer->end(), MF({ 0,0 }) );
	}

PVOCBuffer::MF * PVOCBuffer::getMFPointer( size_t channel, size_t frame, size_t bin )
	{
	return buffer->data() + getBufferPos( channel, frame, bin );
	}

const PVOCBuffer::MF * PVOCBuffer::getMFPointer( size_t channel, size_t frame, size_t bin  ) const
	{
	return buffer->data() + getBufferPos( channel, frame, bin );
	}

size_t PVOCBuffer::getBufferPos( size_t c, size_t f, size_t b ) const
	{
	return c * getNumFrames() * getNumBins() + f * getNumBins() + b;
	}

//======================================================
//	Global
//======================================================

std::ostream & operator<<( std::ostream & os, const PVOCBuffer & pvoc )
	{
	os << "\n=========================== PVOCBuffer Info ==========================="
	   << "\nChannels:\t"				<< pvoc.getNumChannels() 
	   << "\nSamples:\t"				<< pvoc.getNumFrames() 
	   << "\nBins:\t"					<< pvoc.getNumBins() 
	   << "\nFrames/second:\t"			<< pvoc.timeToFrame() 
	   << "\nBins/Frequency:\t"			<< pvoc.frequencyToBin() 
	   << "\nOverlaps:\t"				<< pvoc.getOverlaps() 
	   << "\n=======================================================================" 
	   << "\n\n";
	return os;
	}

} // End namespace xcdp