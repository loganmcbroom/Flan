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

AudioBuffer::AudioBuffer( std::vector<float> && temp_buffer, Channel num_channels, FrameRate sr )
	: format()
	, buffer( temp_buffer )
	{
	format.num_channels = num_channels;
	format.num_frames = temp_buffer.size() / num_channels;
	format.sample_rate = sr;
	}

AudioBuffer::AudioBuffer( const Format & other )
	: format( other )
	, buffer( get_num_channels() * get_num_frames() )
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

bool AudioBuffer::is_null() const
	{
	return buffer.empty() || get_sample_rate() == 0;
	}

//======================================================
//	I/O
//======================================================

#ifdef USE_SNDFILE // IO using libsndfile  

#include <sndfile.h>

bool AudioBuffer::load( const std::string & filepath ) 
	{
	//Open file and check validity, save the sample rate
	SF_INFO info;
	SNDFILE * file = sf_open( filepath.data(), SFM_READ, &info ); 
	if( file == nullptr )
		{
		std::cout << filepath + " could not be opened.\n";
		return false;
		}

	//Copy file info into format
	format.sample_rate = info.samplerate;
	format.num_channels = info.channels;
	format.num_frames = info.frames;
	*this = AudioBuffer( format );

	//Create temporary buffer for interleaved data in file, read data in, close the file
	std::vector<float> interleaved_buffer( info.frames * info.channels );
	sf_readf_float( file, interleaved_buffer.data(), info.frames );

	if( sf_close( file ) != 0 )
		{
		std::cout << std::string( "Error closing " ) + filepath + ".\n";
		return false;
		}

	//Convert interleaved data in
	for( Channel channel = 0; channel < info.channels; ++channel )
		for( Frame frame = 0; frame < Frame(info.frames); ++frame )
			set_sample( channel, frame, interleaved_buffer[ frame * info.channels + channel ] );

	return true;
	}

bool AudioBuffer::save( const std::string & filepath, int format ) const 
	{
	if( format == -1 ) format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;

	//Check that nothing silly is going on with the file formatting
	SF_INFO info = {};
	info.channels	= int( get_num_channels() );
	info.frames		= int( get_num_frames()	);
	info.samplerate = int( get_sample_rate() );
	info.format = format;
	if( !sf_format_check( &info ) )
		{
		std::cout << std::string( "Sound file formatting invalid while attempting to save to " ) + filepath + ",\n";
		print_summary();
		return false;
		}

	//Create a temporary buffer for interleaved data and copy the buffer in
	std::vector<float> interleaved_buffer( get_num_frames() * get_num_channels() );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		for( Frame frame : iota_view( 0, get_num_frames() ) )
			interleaved_buffer[ frame * info.channels + channel ] = get_sample( channel, frame );

	//Clip all samples in the interleaved buffer
	std::for_each( interleaved_buffer.begin(), interleaved_buffer.end(), []( float & s )
		{
		s = std::clamp( s, -1.0f, 1.0f );
		});

	//Open the file and write in the interleaved buffer
	SF_INFO out_info = info;
	SNDFILE * file = sf_open( filepath.data(), SFM_WRITE, &out_info );
	if( file == nullptr )
		{
		std::cout << filepath + " could not be opened for saving.\n";
		return false;
		}
	if( sf_writef_float( file, interleaved_buffer.data(), info.frames ) != info.frames )	
		{
		std::cout << std::string( "Error writing data into " ) + filepath + ".\n";
		return false;
		}
	sf_close( file );
	
	return true;
	}

#else // io using custom wave handler

bool AudioBuffer::load( const std::string & filepath ) 
	{
	auto bail = [filepath]( const std::string & s )
		{
		std::cout << "Couldn't load " << filepath << ": " << s << std::endl;
		return false;
		};

	std::ifstream file( filepath, std::ios::binary );
	if( !file ) return bail( "Error opening " + filepath + "." );

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
	file.read( (char * ) &int16Buffer, 2 ); format.num_channels = int16Buffer;
	file.read( (char * ) &int32Buffer, 4 ); format.sample_rate = int32Buffer;
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

	format.num_frames = int32Buffer / blockAlign;

	*this = AudioBuffer( format );

	if( wavFormat == 1 )
		{
		if( bitsPerSample == 8 )
			{
			const double limit = std::pow( 2.0, 8.0 - 1.0 );

			//Read data buffer
			for( uint32_t frame = 0; frame < get_num_frames(); ++frame )
				for( uint32_t channel = 0; channel < get_num_channels(); ++channel )
					{
					uint8_t intSample = 0;
					file.read( ( (char*) &intSample ), 1 );
					
					set_sample( channel, frame, double( intSample - 128 ) / limit );
					}
			}
        else if( bitsPerSample == 16 )
			{
			const double limit = std::pow( 2.0, 16.0 - 1.0 );

			//Read data buffer
			for( uint32_t frame = 0; frame < get_num_frames(); ++frame )
				for( uint32_t channel = 0; channel < get_num_channels(); ++channel )
					{
					int16_t intSample = 0;
					file.read( ( (char*) &intSample ), 2 );
			
					set_sample( channel, frame, double( intSample ) / limit );
					}
			}
        else if( bitsPerSample == 24 )
			{
			const double limit = std::pow( 2.0, 24.0 - 1.0 );

           //Read data buffer
			for( uint32_t frame = 0; frame < get_num_frames(); ++frame )
				for( uint32_t channel = 0; channel < get_num_channels(); ++channel )
					{
					int32_t intSample = 0;
					file.read( ( (char*) &intSample ), 3 );

					//If sign bit is on, this is negative and its highest byte should be filled via 2s compliment
					if( intSample & 0x800000 ) intSample |= 0xFF000000;
			
					set_sample( channel, frame, double( intSample ) / limit );
					}
			}
		else if( bitsPerSample == 32 )
			{
			const double limit = std::pow( 2.0, 32.0 - 1.0 );

           //Read data buffer
			for( uint32_t frame = 0; frame < get_num_frames(); ++frame )
				for( uint32_t channel = 0; channel < get_num_channels(); ++channel )
					{
					int32_t intSample = 0;
					file.read( ( (char*) &intSample ), 4 );
			
					set_sample( channel, frame, double( intSample ) / limit );
					}
			}
		else
			{
			return bail( "Only 8, 16, 24, and 32 bit wave files are supported without libsndfile enabled." );
			}
		}
	else if( wavFormat == 3 ) // float data
		{
		for( uint32_t frame = 0; frame < get_num_frames(); ++frame )
			for( uint32_t channel = 0; channel < get_num_channels(); ++channel )
				{
				float f;
				file.read( (char *) &f, 4 );
				set_sample( channel, frame, littleEndianToCurrent( f ) );
				}
		}
	else
		return bail( std::string( "Can't load a wave with format " ) + std::to_string( wavFormat ) + "." );
	
	return true;
	}

bool AudioBuffer::save( const std::string & filepath, int format ) const 
	{
	//Make sure path ends in wav
	if( filepath.substr( filepath.size() - 4, 4 ) != ".wav" )
		{
		std::cout << "Libsndfile is required for saving audio in formats other than wave"
			<< " (attempted to save a " << filepath.substr( filepath.size() - 4, 4 ) << " file ).\n";
		return false;
		}

	const uint16_t byteDepth = 3;
	const double limit = std::pow( 2.0, 8 * byteDepth - 1 );

	//Convert buffer to 24bit ints
	std::vector<uint8_t> byteBuffer( buffer->size() * byteDepth );
	for( uint32_t channel = 0; channel < get_num_channels(); ++channel )
		for( uint32_t frame = 0; frame < get_num_frames(); ++frame )
			{
			double clamped = std::clamp( get_sample( channel, frame ), -1.0f, 1.0f );
			clamped *= limit;
			int32_t intSample = clamped;

			byteBuffer[ ( frame * get_num_channels() + channel ) * byteDepth + 0 ] = (uint8_t) ( intSample >> 0  ) & 0xFF;
			byteBuffer[ ( frame * get_num_channels() + channel ) * byteDepth + 1 ] = (uint8_t) ( intSample >> 8  ) & 0xFF;
			byteBuffer[ ( frame * get_num_channels() + channel ) * byteDepth + 2 ] = (uint8_t) ( intSample >> 16 ) & 0xFF;
			}

	//Write to file
	return writeRIFF( filepath, "WAVE", byteBuffer.data(), byteBuffer.size(),
		{
		uint16_t( 1 ),													// Formatting
		uint16_t( get_num_channels() ),									// Channel Count
		uint32_t( get_sample_rate() ),									// Sample rate
		uint32_t( get_sample_rate() * get_num_channels() * byteDepth ),		// Byte rate
		uint16_t( get_num_channels() * byteDepth ),						// Blockalign
		uint16_t( byteDepth * 8 )										// Bits per sample
		});
	}

#endif

void AudioBuffer::print_summary() const
	{
	std::cout << *this;
	}

//======================================================
//	Getters
//======================================================
auto AudioBuffer::get_sample( Channel channel, Frame frame ) const -> Sample
	{
	return buffer[get_buffer_pos( channel, frame )];
	}

AudioBuffer::Format flan::AudioBuffer::get_format() const
	{
	return format;
	}

auto AudioBuffer::get_num_channels() const -> Channel
	{
	return format.num_channels;
	}

auto AudioBuffer::get_num_frames() const -> Frame
	{
	return format.num_frames;
	}

FrameRate AudioBuffer::get_sample_rate() const 
	{
	return format.sample_rate;
	}

Second AudioBuffer::frame_to_time( fFrame f ) const
	{
	return f / get_sample_rate();
	}

fFrame AudioBuffer::time_to_frame( Second t ) const
	{
	return t * float( get_sample_rate() );
	}

auto AudioBuffer::get_length() const -> Second
	{
	return frame_to_time( get_num_frames() );
	}

float AudioBuffer::get_max_sample_magnitude( Second start_time, Second end_time ) const
	{
	if( end_time == 0 ) end_time = get_length();
	auto start_frame = std::clamp( (Frame) time_to_frame( start_time ), 0, get_num_frames() - 1 );
	auto end_frame   = std::clamp( (Frame) time_to_frame( end_time   ), 0, get_num_frames() - 1 );
	Magnitude m = 0;
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		for( Frame frame = start_frame; frame < end_frame; ++frame )
			{
			const Sample & s = get_sample( channel, frame );
			if( std::abs( s ) > m )
				m = std::abs( s );
			}
	return m;
	}

//======================================================
//	Setters
//======================================================
void AudioBuffer::set_sample( Channel channel, Frame frame, Sample sample ) 
	{
	buffer[get_buffer_pos( channel, frame )] = sample;
	}

auto AudioBuffer::get_sample( Channel channel, Frame frame ) -> Sample &
	{
	return buffer[get_buffer_pos( channel, frame )];
	}

void flan::AudioBuffer::clear_buffer()
	{
	std::fill( buffer.begin(), buffer.end(), 0 );
	}

auto AudioBuffer::get_sample_pointer( Channel channel, Frame frame ) -> Sample *
	{ 
	return buffer.data() + get_buffer_pos( channel, frame );
	}

auto AudioBuffer::get_sample_pointer( Channel channel, Frame frame ) const -> const Sample *
	{ 
	return buffer.data() + get_buffer_pos( channel, frame );
	}

std::vector<Sample> & AudioBuffer::get_buffer() 
	{ 
	return buffer; 
	}
const std::vector<Sample> & AudioBuffer::get_buffer() const
	{
	return buffer; 
	}

std::vector<Sample>::const_iterator AudioBuffer::channel_begin( Channel channel ) const
	{
	return buffer.begin() + channel * get_num_frames();
	}

std::vector<Sample>::const_iterator AudioBuffer::channel_end( Channel channel ) const
	{
	return buffer.begin() + ( channel + 1 ) * get_num_frames();
	}

size_t AudioBuffer::get_buffer_pos( Channel channel, Frame sample ) const
	{
	return channel * get_num_frames() + sample;
	}

//======================================================
//	Global
//======================================================
std::ostream & operator<<( std::ostream & os, const AudioBuffer & audio )
	{
	os << "\n=========================== Audio Info ==========================="
	   << "\nChannels:\t"		<< audio.get_num_channels() 
	   << "\nSamples:\t"		<< audio.get_num_frames() 
	   << "\nSample Rate:\t"	<< audio.get_sample_rate()
	   << "\n==================================================================" 
	   << "\n\n";
	return os;
	}

} // End namespace flan
