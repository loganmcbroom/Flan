#include "flan/Audio/Audio.h"

using namespace flan;

Audio::Audio() 
	: AudioBuffer() 
	{}

Audio::Audio( AudioBuffer && other )
	: AudioBuffer( std::move( other ) )
	{
	}

Audio Audio::copy() const 
	{ 
	return AudioBuffer::copy(); 
	}

Audio Audio::create_null()
	{
	std::cout << "Null Audio created\n";
	return Audio();
	}

Audio Audio::create_from_buffer( std::vector<float> && buffer, Channel num_channels, FrameRate sr )
	{
	return AudioBuffer( std::move( buffer ), num_channels, sr );
	}

Audio Audio::create_from_format( const AudioBuffer::Format & other ) 
	{
	return AudioBuffer( other );
	}

Audio Audio::load_from_file( const std::string & filename ) 
	{
	return AudioBuffer( filename );
	}

Audio Audio::load_from_file( 
	const std::string & filename,
	SndfileStrings & strings
	)
	{
	return AudioBuffer( filename, strings );
	}

Audio Audio::create_empty_with_length( 
	Second length, 
	Channel num_channels, 
	FrameRate sample_rate 
	)
	{
	return create_empty_with_frames( length * sample_rate, num_channels, sample_rate );
	}

Audio Audio::create_empty_with_frames( 
	Frame num_frames, 
	Channel num_channels, 
	FrameRate sample_rate 
	)
	{
	Audio::Format format;
	format.num_channels = num_channels;
	format.num_frames = num_frames;
	format.sample_rate = sample_rate;
	Audio out( format );
	out.clear_buffer();
	return out;
	}