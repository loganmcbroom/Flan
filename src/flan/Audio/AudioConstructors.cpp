#include "flan/Audio/Audio.h"

using namespace flan;

Audio Audio::copy() const { return AudioBuffer::copy(); }

Audio::Audio() 
	: AudioBuffer() 
	{}

Audio::Audio( std::vector<float> && buffer, Channel num_channels, FrameRate sr )
	: AudioBuffer( std::move( buffer ), num_channels, sr )
	{}

Audio::Audio( const AudioBuffer::Format & other ) 
	: AudioBuffer( other ) 
	{}

Audio::Audio( const std::string & filename ) 
	: AudioBuffer( filename ) 
	{}

Audio::Audio( AudioBuffer && other )
	: AudioBuffer( std::move( other ) )
	{
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