#include "PVOCBuffer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Utility.h"

namespace xcdp {

PVOCBuffer::PVOCBuffer( const Format & other )
	: format( other )
	, buffer( std::make_shared<std::vector<MFPair>>( getNumChannels() * getNumFrames() * getNumBins() ))
	{}

PVOCBuffer::PVOCBuffer( const std::string & filename )
	: format()
	, buffer()
	{
	load( filename );
	}

//======================================================
//	IO - Custom RIFF format 
// We could use the CDP .ana format, but to quote their documentation:
//		- "The Spectral Domain processes require ana format inputs. This is all you really need to know about them."
//======================================================

bool PVOCBuffer::save( const std::string & filename ) const
	{
	auto bail = []( const std::string & s )
		{
		std::cout << s << std::endl;
		return false;
		};

	std::ofstream file( filename, std::ios::binary );
	if( !file ) return bail( "Error opening " + filename + " to write BMP." );

	const uint32_t RIFFChunkSize = 12;
	const uint32_t FMTChunkSize = 34;
	const uint32_t dataChunkSize = 8 + buffer->size() * 2 * sizeof(float);

	//Write RIFF chunk
	unsigned char RIFFChunk[ RIFFChunkSize ] = { 0 };
		writeBytes( RIFFChunk + 0, "RIFF"										); //RIFF ID 
		writeBytes( RIFFChunk + 4, (uint32_t) 4 + FMTChunkSize + dataChunkSize ); //RIFF Chunk Size
		writeBytes( RIFFChunk + 8, "PVOC"										); //File Type
	file.write( (const char *) RIFFChunk, RIFFChunkSize );

	//Write format information chunk
	unsigned char FMTChunk[ FMTChunkSize ] = { 0 };
		writeBytes( FMTChunk +   0, "fmt "						); //RIFF ID
		writeBytes( FMTChunk +   4, (uint32_t) FMTChunkSize - 8	); //RIFF Chunk Size
		writeBytes( FMTChunk +   8, (uint16_t) 1				); //Formatting
		writeBytes( FMTChunk +  10, (uint16_t) getNumChannels()	); //Channel Count
		writeBytes( FMTChunk +  12, (uint32_t) getNumFrames()	); //Number of frames
		writeBytes( FMTChunk +  16, (uint32_t) getNumBins()		); //Number of bins per frame
		writeBytes( FMTChunk +  20, (uint32_t) getSampleRate()	); //Sample Rate, note this is the sample rate of the audio, not the analysis frame rate
		writeBytes( FMTChunk +  24, (uint32_t) getOverlaps()	); //Number of windows overlapping at any frame
		writeBytes( FMTChunk +  28, (uint32_t) 64				); //Bit depth, note that each bin contains two of this in magnitude frequency order
		writeBytes( FMTChunk +  32, (uint16_t) 1				); //Window type indicator. 1 = hann.
	file.write( (const char *) FMTChunk, FMTChunkSize );

	unsigned char dataChunkInfo[ 8 ] = { 0 };
		writeBytes( dataChunkInfo + 0, "data"						); //RIFF ID
		writeBytes( dataChunkInfo + 4, (uint32_t) dataChunkSize - 8	); //RIFF Chunk Size
	file.write( (const char *) dataChunkInfo, 8 );

	for( int channel = 0; channel < getNumChannels(); ++channel )
		for( int frame = 0; frame < getNumFrames(); ++frame )
			for( int bin = 0; bin < getNumBins(); ++bin )
				{
				const auto mf = getBin( channel, frame, bin );

				const float m = makeLittleEndian( mf.magnitude );
				const float f = makeLittleEndian( mf.frequency );

				file.write( (char *) &m, sizeof( float ) ); 
				file.write( (char *) &f, sizeof( float ) ); 
				}

	file.close();
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
	if( !file ) return bail( "Error opening " + filename + " to write BMP." );

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
	file.read( (char * ) &int16Buffer, 2 ); //formatting
	file.read( (char * ) &int16Buffer, 2 ); format.numChannels = int16Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.numFrames = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.numBins = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.sampleRate = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.overlaps = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); if( int32Buffer != 32 ) return bail( "Using floats only here" ); //bit depth
	file.read( (char * ) &int16Buffer, 2 ); if( int16Buffer != 1 ) return bail( "Just use Hann window, nerd" ); //window type
	*this = PVOCBuffer( format );

	//Read data subchunk
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "data", 4 ) != 0 ) return bail( filename + " isn't a correctly formatted PVOC file (\"data\" wasn't at the start of the data chunk).\n" );
	file.read( (char * ) &int32Buffer, 4 );
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				float m, f;
				file.read( (char *) &m, sizeof( float ) );
				file.read( (char *) &f, sizeof( float ) );
				setBin( channel, frame, bin, { littleEndianToCurrent(m), littleEndianToCurrent(f) } );
				}
	}

//======================================================
//	Getters
//======================================================

PVOCBuffer::MFPair PVOCBuffer::getBin( size_t channel, size_t frame, size_t bin ) const
	{
	return (*buffer)[getPos( channel, frame, bin )];
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
size_t PVOCBuffer::getFrameSize() const
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
float PVOCBuffer::getBinWidth() const
	{
	return float( getSampleRate() ) / float( getFrameSize() );
	}
float PVOCBuffer::getMaxPartialMagnitude() const
	{
	float maxMagnitude = 0;
	for( auto i : *buffer )
		{
		const float trueMag = std::abs( i.magnitude );
		if( trueMag > maxMagnitude )
			maxMagnitude = trueMag;
		}
	return maxMagnitude;
	}

size_t PVOCBuffer::timeToFrame( float t ) const
	{
	return (size_t) timeToFrame_d( t );
	}
float PVOCBuffer::timeToFrame_d( float t ) const
	{
	return t * float(getSampleRate()) * float(getOverlaps()) / float(getFrameSize());
	}
float PVOCBuffer::frameToTime( size_t frame ) const
	{
	return frame * float(getFrameSize()) / float(getSampleRate()) / float(getOverlaps());
	}

size_t PVOCBuffer::frequencyToBin( float freq ) const
	{
	return (size_t) frequencyToBin_d( freq );
	}
float PVOCBuffer::frequencyToBin_d( float freq ) const
	{
	return freq / getBinWidth();
	}
float PVOCBuffer::binToFrequency( size_t bin ) const
	{
	return float( bin ) * getBinWidth();
	}

//======================================================
//	Setters
//======================================================

void PVOCBuffer::setBin( size_t channel, size_t frame, size_t bin, MFPair value )
	{
	(*buffer)[getPos( channel, frame, bin )] = value;
	}
PVOCBuffer::MFPair & PVOCBuffer::getBin( size_t channel, size_t frame, size_t bin )
	{
	return (*buffer)[getPos( channel, frame, bin )];
	}

void PVOCBuffer::clearBuffer()
	{
	std::fill( buffer->begin(), buffer->end(), MFPair({ 0,0 }) );
	}

size_t PVOCBuffer::getPos( size_t c, size_t f, size_t b ) const
	{
	return c * getNumFrames() * getNumBins() + f * getNumBins() + b;
	}

} // End namespace xcdp