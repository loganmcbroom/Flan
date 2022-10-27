#include "flan/PVOCBuffer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

#include "flan/Utility/Bytes.h"

namespace flan {

PVOCBuffer::PVOCBuffer()
	: format()
	, buffer( std::make_shared< std::vector<MF> >() )
	{}

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

PVOCBuffer PVOCBuffer::deepCopy() const
	{
	PVOCBuffer out;
	out.format = format;
	*out.buffer = *buffer; // Deep copy
	return out;
	}

bool PVOCBuffer::isNull() const
	{
	return !buffer || buffer->size() == 0 || getSampleRate() == 0;
	}

// PVOC-EX structure
//struct WAVEFORMATPVOCEX					// 80 bytes
//	{
//	struct WAVEFORMATEXTENSIBLE				// 40 bytes
//		{
//		struct WAVEFORMATEX						// 18 bytes, info for renderer as well as for flan
//			{ 
//			uint16_t formatTag = 0xFFFE;			// WAVE_FORMAT_EXTENSIBLE macro
//			uint16_t numChannels;					// Number of channels
//			uint32_t sampleRate;					// Sample rate
//			uint32_t byteRate;						// Bytes per second
//			uint16_t blockAlign;					// Block align
//			uint16_t bitsPerSample;					// Bits per sample
//			uint16_t cbSize;						// Bytes after WAVEFORMATEX data
//			} waveEx; 
//		union									// 2 bytes
//			{								 
//			uint16_t bitsPerSample;					// It's bitsPerSample again
//			uint16_t samplesPerBlock;				// Unused here
//			uint16_t reserved;						// Unused here
//			} samples;
//		uint32_t channelMask;					// 4 bytes, describes which channel positions are present                             
//		struct GUID								// 16 bytes
//			{ 
//			uint32_t data1	= 0x8312b9c2;						
//			uint16_t data2	= 0x2e6e;						
//			uint16_t data3	= 0x11d4;						
//			uint8_t  data4[8] = { 0xa8, 0x24, 0xde, 0x5b, 0x96, 0xc3, 0xab, 0x21 } ;					
//			} guid;							
//		} waveExt;
//	uint32_t version  = 1;				// 4 bytes
//	uint32_t dataSize = 32;				// 4 bytes, size of PVOCDATA data block
//	struct PVOCDATA							// 32 bytes
//		{				
//		uint16_t dataType;						// 0 (float) or 1 (double)
//		uint16_t analFormat;					// 0 = Amp/Freq, 1 = Amp/Phase, 2 = Complex 
//		uint16_t sourceFormat;					// 1 = WAVE_FORMAT_PCM, 3 = WAVE_FORMAT_IEEE_FLOAT
//		uint16_t windowType;					// window type, another fucking macro
//		uint32_t analysisBins;					// implicit FFT size = (nAnalysisBins-1) * 2
//		uint32_t windowlen;						// analysis window length, samples, NB may be <> FFT size 
//		uint32_t analStep;						// audio frames per flan frame
//		uint32_t frameAlign;					// usually nAnalysisBins * 2 * sizeof(float) 
//		float analysisRate;						// Sample rate / overlaps
//		float windowParam;						// default 0.0f unless needed 
//		} flanData;
//	};

bool PVOCBuffer::save( const std::string & filename ) const
	{
	const int byteDepth = 3;
	const double limit = std::pow( 2, 8 * byteDepth - 1 );
	const float windowSize_f = getDFTSize();
	const float maxFreq_f = getSampleRate();

	//Convert buffer to 24bit signed int representation
	std::vector<uint8_t> bytes( buffer->size() * 2 * byteDepth );
	for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
		for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
			for( uint32_t bin = 0; bin < getNumBins(); ++bin )
				{
				int32_t m_32 = double( std::clamp( getMF( channel, frame, bin ).m / windowSize_f, -1.0f, 1.0f ) ) * limit;
				int32_t f_32 = double( std::clamp( getMF( channel, frame, bin ).f / maxFreq_f,	  -1.0f, 1.0f ) ) * limit;

				uint32_t pos = getBufferPos( channel, frame, bin ) * 2 * byteDepth;

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
		(uint16_t) 1,					// Formatting
		(uint16_t) getNumChannels(),	// Channel Count
		(uint32_t) getNumFrames(),		// Number of frames
		(uint32_t) getNumBins(),		// Number of bins per frame
		(uint32_t) getSampleRate(),		// Sample Rate, note this is the sample rate of the audio
		(uint32_t) getHopSize(),		// Number of Audio frames jumped per dft
		(uint32_t) getWindowSize(),		// The number of audio frames used per fft. Used when audio data is zero padded.
		(uint32_t) 24,					// Bit depth, note that each bin contains two of this
		(uint16_t) 1					// Window type indicator. 1 = hann.
		});

	return true;

	//const int bytesPerSample = sizeof( float );

	//WAVEFORMATPVOCEX flanFormat;

	////Fill format structure	
	//flanFormat.waveExt.waveEx.numChannels		= getNumChannels(); 
	//flanFormat.waveExt.waveEx.stimampleRate		= getSampleRate(); // Frames/Sec
	//flanFormat.waveExt.waveEx.byteRate			= getNumChannels() * bytesPerSample * timeToFrame(); 
	//flanFormat.waveExt.waveEx.blockAlign		= getNumChannels() * bytesPerSample;
	//flanFormat.waveExt.waveEx.bitsPerSample		= 8 * ifbytesPerSample; 
	//flanFormat.waveExt.waveEx.cbSize			= 62;
	//flanFormat.waveExt.samples.bitsPerSample	= 8 * bytesPerSample;
	//flanFormat.waveExt.channelMask				= 0; // I want nothing to do with channel masking
	//flanFormat.flanData.dataType				= 0; // Float			
	//flanFormat.flanData.analFormat				= 0; // Amp/Freq			
	//flanFormat.flanData.sourceFormat			= 3; // WAVE_FORMAT_IEEE_FLOAT
	//flanFormat.flanData.windowType				= 2; // Hanning window			
	//flanFormat.flanData.analysisBins			= getNumBins();				
	//flanFormat.flanData.windowlen				= getWindowSize();				
	//flanFormat.flanData.analStep				= getWindowSize() / getOverlaps();				
	//flanFormat.flanData.frameAlign				= getNumBins() * 2 * bytesPerSample;				
	//flanFormat.flanData.analysisRate			= float( getSampleRate() ) / float( getWindowSize() / getOverlaps() );
	//flanFormat.flanData.windowParam				= 0; //Other values can be used for other windows		
	
	//Handle window chunk?

	//Write format structure to file
	//std::ofstream file( filename, std::ios::binary );
	//if( !file )
	//	{
	//	std::cout << "Error opening " + filename + " to write RIFF.\n";
	//	return false;
	//	}

	//unsigned char * write;

	////Write RIFF chunk
	//unsigned char RIFFChunk[ RIFFChunkSize ] = { 0 };
	//write = RIFFChunk;
	//writeBytes( write, "RIFF"		); write += 4; //RIFF ID 
	//	writeBytes( write, (uint32_t) 4 ); write += 4; //RIFF Chunk Size
	//	writeBytes( write, type		    ); write += 4;  //File Type
	//file.write( (const char *) RIFFChunk, RIFFChunkSize );

	////Write format chunk
	//std::vector<unsigned char> FMTChunk( FMTChunkSize, 0 );
	//write = FMTChunk.data();
	//	writeBytes( write, "fmt "						); write += 4;
	//	writeBytes( write, (uint32_t) FMTChunkSize - 8	); write += 4;
	//	for( int i = 0; i < format.size(); ++i )
	//		{
	//		if( format[i].numBytes == 2 ) writeBytes( write, uint16_t( format[i].value ) ); 
	//		else						  writeBytes( write, uint32_t( format[i].value ) ); 
	//		write += format[i].numBytes;
	//		}
	//file.write( (const char *) FMTChunk.data(), FMTChunkSize );

	////Write data chunk
	//unsigned char dataChunkInfo[ 8 ] = { 0 };
	//write = dataChunkInfo;
	//	writeBytes( write, "data"						); write += 4;
	//	writeBytes( write, (uint32_t) dataChunkSize - 8	); write += 4;
	//file.write( (const char *) dataChunkInfo, 8 );

	//file.write( (const char *) data, dataSize );

	//file.close();

	//Write buffer to file

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
	file.read( strBuffer, 4 ); if( std::strncmp( strBuffer, "RIFF", 4 ) != 0 ) return bail( filename + " isn't a correctly formatted RIFF file.\n" );
	file.read( strBuffer, 4 ); //Don't care about file size
	file.read( strBuffer, 4 ); if( std::strncmp( strBuffer, "PVOC", 4 ) != 0 ) return bail( filename + " isn't a PVOC file.\n"  );

	//Read fmt subchunk
	PVOCBuffer::Format format;
	file.read( strBuffer, 4 ); if( std::strncmp( strBuffer, "fmt ", 4 ) != 0 ) return bail( filename + " isn't formatted correctly (\"fmt \" wasn't at the start of the format chunk).\n" );
	file.read( (char * ) &int32Buffer, 4 ); //Chunk size
	file.read( (char * ) &int16Buffer, 2 );  if( int16Buffer != 1 ) return bail( "Formatting must be 1 (signed int)." );
	file.read( (char * ) &int16Buffer, 2 ); format.numChannels = int16Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.numFrames = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.numBins = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.sampleRate = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.hopSize = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.windowSize = int32Buffer;
	file.read( (char * ) &int32Buffer, 4 ); if( int32Buffer != 24 ) return bail( "Bit depth must be 24." );
	file.read( (char * ) &int16Buffer, 2 ); if( int16Buffer != 1 ) return bail( "PVOC window must be 1 (Hann)." );
	*this = PVOCBuffer( format );

	//Read data subchunk
	file.read( strBuffer, 4 ); if( strncmp( strBuffer, "data", 4 ) != 0 ) return bail( filename + " isn't a correctly formatted PVOC file (\"data\" wasn't at the start of the data chunk).\n" );
	file.read( (char * ) &int32Buffer, 4 );

	
	//Read buffer
	const double limit = std::pow( 2, 23 );
	const float windowSize_f = getDFTSize();
	const float maxFreq_f = getSampleRate();

	auto getFloat = [limit, &file]( float div )
		{
		int32_t i = 0;
		file.read( ( (char*) &i ), 3 );
		if( i & 0x800000 ) i |= 0xFF000000;
		return float( double( i ) / limit ) * div;
		};

	for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
		for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
			for( uint32_t bin = 0; bin < getNumBins(); ++bin )
				setMF( channel, frame, bin, { getFloat( windowSize_f ), getFloat( maxFreq_f ) } );

	return true;		

	//auto bail = [&filename]( const std::string & s )
	//	{
	//	std::cout << "Error loading PVOC data from " << filename << ": " << s << std::endl;
	//	return false;
	//	};

	//std::ifstream file( filename, std::ios::binary );
	//if( !file ) return bail(  "Couldn't open file." );

	////Read entire format chunk
	//WAVEFORMATPVOCEX fileFormat;
	//file.read( (char *) &fileFormat, 80 ); 

	//for( int i = 0; i < 80; ++i )
	//	std::cout << (uint16_t)((uint8_t *) &fileFormat)[i] << " ";

	//// Check format chunk and create PVOCBuffer::Format
	//PVOCBuffer::Format flanFormat;
	//if( fileFormat.waveExt.waveEx.formatTag != 0xFFFE ) return bail( "WAVE-EX tag wasn't 0xFFFE" );
	//flanFormat.numChannels = fileFormat.waveExt.waveEx.numChannels;
	//flanFormat.sampleRate = fileFormat.waveExt.waveEx.sampleRate;
	//if( fileFormat.waveExt.waveEx.bitsPerSample != 32 ) return bail( "flan only supports 32-bit floats (bit rate wasn't 32)." );
	//if( fileFormat.waveExt.waveEx.cbSize != 62 ) return bail( "WAVE-EX cbSize wasn't 62." );
	//if( fileFormat.waveExt.waveEx.cbSize != 62 ) std::cout << "Discarding channel mask when loading PVOC data (unsupported)." << std::endl;
	//WAVEFORMATPVOCEX::WAVEFORMATEXTENSIBLE::GUID correctGUID;
	//if( !memcmp( &fileFormat.waveExt.guid, &correctGUID, sizeof( WAVEFORMATPVOCEX::WAVEFORMATEXTENSIBLE::GUID ) ) )
	//	return bail( "Incorrect GUID." );
	//if( fileFormat.version != 1 ) return bail( "PVOC-EX versions above 1 aren't supported." );
	//if( fileFormat.dataSize != 32 ) return bail( "PVOC-EX versions above 1 aren't supported." );

	//if( fileFormat.flanData.dataType != 0 ) return bail( "PVOC data type must be float." );
	//if( fileFormat.flanData.analFormat != 0 ) return bail( "PVOC analysis format must be Amp/Freq." );
	//if( fileFormat.flanData.sourceFormat != 3 ) return bail( "PVOC source format must be float." );
	//if( fileFormat.flanData.windowType != 2 ) return bail( "Only Hann window is currently supported." );
	//flanFormat.numBins = fileFormat.flanData.analysisBins;
	//flanFormat.overlaps = ( float( flanFormat.numBins ) - 1.0f ) * 2.0f / float( fileFormat.flanData.analStep );		
	////float windowParam; //Unneeded until other windows are implemented


	////flanFormat.numFrames = ;
	//
	////Read buffer
	//for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
	//	for( uint32_t frame = 0; frame < getNumFrames(); ++frame )
	//		for( uint32_t bin = 0; bin < getNumBins(); ++bin )
	//			{
	//			//setMF( channel, frame, bin, { getFloat( windowSize_f ), getFloat( maxFreq_f ) } );
	//			}

	//return true;		
	}

void PVOCBuffer::printSummary() const
	{
	std::cout << *this;
	}

//======================================================
//	Getters
//======================================================

PVOCBuffer::MF PVOCBuffer::getMF( Channel channel, Frame frame, Bin bin ) const
	{
	return (*buffer)[getBufferPos( channel, frame, bin )];
	}

auto PVOCBuffer::getNumChannels() const -> Channel
	{
	return format.numChannels;
	}

auto PVOCBuffer::getNumFrames() const -> Frame 
	{
	return format.numFrames;
	}

auto PVOCBuffer::getNumBins() const -> Bin
	{
	return format.numBins;
	}

Frame PVOCBuffer::getDFTSize() const
	{
	return ( getNumBins() - 1 ) * 2;
	}

Frame PVOCBuffer::getWindowSize() const
	{
	return format.windowSize;
	}

PVOCBuffer::Format PVOCBuffer::getFormat() const
	{
	return format;
	}

uint32_t PVOCBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}

Frame PVOCBuffer::getHopSize() const
	{
	return format.hopSize;
	}

Time PVOCBuffer::getLength() const
	{
	return getNumFrames() * frameToTime();
	}

Frequency PVOCBuffer::getHeight() const
	{
	return float( getNumBins() ) * binToFrequency();
	}

Magnitude PVOCBuffer::getMaxPartialMagnitude() const
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

Magnitude PVOCBuffer::getMaxPartialMagnitude( uint32_t startFrame, uint32_t endFrame, uint32_t startBin, uint32_t endBin ) const
	{
	if( endFrame == 0 ) endFrame = getNumFrames();
	if( endBin == 0 ) endBin = getNumBins();

	float maxMagnitude = 0;

	for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
		for( uint32_t frame = startFrame; frame < endFrame; ++frame )
			for( uint32_t bin = startBin; bin < endBin; ++bin )
			    {
				const float trueMag = std::abs( getMF( channel, frame, bin ).m );
				if( trueMag > maxMagnitude )
					maxMagnitude = trueMag;

				}

	return maxMagnitude;
	}

float PVOCBuffer::timeToFrame() const
	{
	return float( getSampleRate() ) / float( getHopSize() );
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
	return float( getSampleRate() ) / float( getDFTSize() );
	}

float PVOCBuffer::getFrequencyOffset( Channel c, Frame f, Bin b ) const
	{
	return getMF( c, f, b ).f - b * binToFrequency();
	}

Channel PVOCBuffer::boundChannel( Channel c ) const
	{
	return std::clamp( c, 0, getNumChannels() - 1 );
	}

Frame PVOCBuffer::boundFrame( Frame f ) const
	{
	return std::clamp( f, 0, getNumFrames() - 1 );
	}

Bin PVOCBuffer::boundBin( Bin b ) const
	{
	return std::clamp( b, 0, getNumBins() - 1 );
	}

//======================================================
//	Setters
//======================================================

void PVOCBuffer::setMF( Channel channel, Frame frame, Bin bin, MF value )
	{
	(*buffer)[getBufferPos( channel, frame, bin )] = value;
	}
PVOCBuffer::MF & PVOCBuffer::getMF( Channel channel, Frame frame, Bin bin )
	{
	return (*buffer)[getBufferPos( channel, frame, bin )];
	}

void PVOCBuffer::clearBuffer()
	{
	std::fill( buffer->begin(), buffer->end(), MF({ 0,0 }) );
	}

PVOCBuffer::MF * PVOCBuffer::getMFPointer( Channel channel, Frame frame, Bin bin )
	{
	return buffer->data() + getBufferPos( channel, frame, bin );
	}

const PVOCBuffer::MF * PVOCBuffer::getMFPointer( Channel channel, Frame frame, Bin bin  ) const
	{
	return buffer->data() + getBufferPos( channel, frame, bin );
	}

std::vector<PVOCBuffer::MF>::iterator PVOCBuffer::begin()
	{
	return buffer->begin();
	}

std::vector<PVOCBuffer::MF>::iterator PVOCBuffer::end()
	{
	return buffer->end();
	}

std::vector<PVOCBuffer::MF>::const_iterator PVOCBuffer::begin() const
	{
	return buffer->begin();
	}

std::vector<PVOCBuffer::MF>::const_iterator PVOCBuffer::end() const
	{
	return buffer->end();
	}

size_t PVOCBuffer::getBufferPos( Channel c, Frame f, Bin b ) const
	{
	return c * getNumFrames() * getNumBins() + f * getNumBins() + b;
	}

//======================================================
//	Global
//======================================================

std::ostream & operator<<( std::ostream & os, const PVOCBuffer & flan )
	{
	os << "\n=========================== PVOCBuffer Info ==========================="
	   << "\nChannels:\t"				<< flan.getNumChannels() 
	   << "\nSamples:\t"				<< flan.getNumFrames() 
	   << "\nBins:\t"					<< flan.getNumBins() 
	   << "\nFrames/second:\t"			<< flan.timeToFrame() 
	   << "\nBins/Frequency:\t"			<< flan.frequencyToBin() 
	   << "\nHop size:\t"				<< flan.getHopSize() 
	   << "\nDFT size:\t"				<< flan.getDFTSize() 
	   << "\n=======================================================================" 
	   << "\n\n";
	return os;
	}

} // End namespace flan
