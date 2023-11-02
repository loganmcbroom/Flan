#include "flan/Audio/Audio.h"
 
using namespace flan;

template<typename T>
static std::vector<const T *> get_pointers( const std::vector<T> & ts )
	{
	std::vector<const T *> ptrs( ts.size() );
	std::transform( ts.begin(), ts.end(), ptrs.begin(), []( const T & t ){ return &t; } );
	return std::move( ptrs );
	}

std::vector<Audio> Audio::split_channels(
	) const
	{
	Audio::Format format;
	format.num_channels = 1;
	format.num_frames = get_num_frames();
	format.sample_rate = get_sample_rate();

	std::vector<Audio> channels;
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		channels.emplace_back( Audio( format ) );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			channels[channel].get_sample( 0, frame ) = get_sample( channel, frame );
		}
	return std::move( channels );
	}

Audio Audio::combine_channels( 
	const std::vector<const Audio *> & channels 
	)
	{
	if( channels.empty() ) return Audio::create_null();

	// Use the total number of input channels, and the maximum number of input frames
	Audio::Format format;
	format.sample_rate = channels[0]->get_sample_rate();
	format.num_channels = std::accumulate( channels.begin(), channels.end(), 0, 
		[]( int c, const Audio * p ){ return c + p->get_num_channels(); } );
	format.num_frames = (*std::max_element( channels.begin(), channels.end(), []( const Audio * a, const Audio * b )
		{ return a->get_num_frames() < b->get_num_frames(); } ))->get_num_frames();
	Audio out( format ); 
	out.clear_buffer();

	Channel current_channel = 0;
	for( auto & a : channels )
		{
		for( Channel channel = 0; channel < a->get_num_channels(); ++channel )
			{
			for( Frame frame = 0; frame < a->get_num_frames(); ++frame )
				out.set_sample( current_channel, frame, a->get_sample( channel, frame ) );
			++current_channel;
			}
		}

	return out;
	}

 Audio Audio::combine_channels( const std::vector<Audio> & channels )
	{
	return combine_channels( get_pointers( channels ) );
	}