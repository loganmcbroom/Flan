#include "flan/Audio/Audio.h"

#include "WDL/resample.h"

using namespace flan;

Audio Audio::modify_boundaries( 
	Second start, 
	Second end 
	) const
	{
	flan_PROCESS_START( Audio() );

	const Frame num_out_frames = time_to_frame( -start ) + get_num_frames() + time_to_frame( end );
	if( num_out_frames <= 0 ) return Audio();

	Audio::Format format = get_format();
	format.num_frames = num_out_frames;
	Audio out( format );
	out.clear_buffer();

	out.mix_in_place( *this, -start );

	return out;
	}

Audio Audio::remove_edge_silence( 
	Amplitude non_silent_level 
	) const
	{
	Frame start = 0;
	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		{
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			if( get_sample( channel, frame ) > non_silent_level )
				goto exit_start_loop;
		++start;	
		}
	exit_start_loop:

	Frame end = get_num_frames();
	for( Frame frame = get_num_frames() - 1; frame >= 0; --frame )
		{
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			if( get_sample( channel, frame ) > non_silent_level )
				goto exit_end_loop;
		--end;	
		}
	exit_end_loop:

	return cut_frames( start, end );
	}

Audio Audio::remove_silence(
	Amplitude non_silent_level,
	Second minimum_gap
	) const
	{
	// Find all silent frames
	std::vector<Frame> silent_frames;
	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		{
		// Get max volume over channels
		Amplitude max_frame_amplitude = 0;
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			if( std::abs( get_sample( channel, frame ) ) > max_frame_amplitude )
				max_frame_amplitude = std::abs( get_sample( channel, frame ) );
		
		if( max_frame_amplitude <= non_silent_level )
			silent_frames.push_back( frame );
		}

	const Frame minimum_gap_frames = time_to_frame( minimum_gap );

	// Find all groups of silent frames longer than the gap
	std::vector<std::array<Frame, 2>> silent_sectors;
	for( Index i = 0; i < silent_frames.size(); ++i )
		{
		const Frame gap_start = silent_frames[i];
		while( i+1 < silent_frames.size() && silent_frames[i+1] - silent_frames[i] == 1 )
			++i;
		const Frame gap_end = silent_frames[i] + 1; // The gap end is, in c++ style, the frame after the last frame to be removed

		if( gap_end - gap_start > minimum_gap_frames )
			silent_sectors.push_back( {gap_start, gap_end} );
		}
	if( silent_sectors.empty() )
		return copy();

	const Frame num_removed_frames = std::accumulate( silent_sectors.begin(), silent_sectors.end(), 0, 
		[]( Frame f, std::array<Frame, 2> s ){ return f + s[1] - s[0]; } );

	Audio::Format format = get_format();
	format.num_frames -= num_removed_frames;
	Audio out( format );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		Index gap_index = 0;
		Frame in_frame = 0;
		Frame out_frame = 0;
		while( true )
			{
			// If we hit a silent gap, skip ahead in the input
			if( in_frame == silent_sectors[gap_index][0] )
				{
				in_frame += silent_sectors[gap_index][1] - silent_sectors[gap_index][0];
				++gap_index;
				}

			if( in_frame >= get_num_frames() || out_frame >= out.get_num_frames() )
				break;

			out.get_sample( channel, out_frame ) = get_sample( channel, in_frame );

			++in_frame;
			++out_frame;
			}
		}

	return out;
	}

Audio Audio::reverse(
	) const
	{
	flan_PROCESS_START( Audio() );

	Audio out( get_format() );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		const auto iter = get_buffer().begin() + channel * get_num_frames();
		std::copy( std::execution::par_unseq, iter, iter + get_num_frames(), out.get_buffer().rbegin() + channel * get_num_frames() );
		}

	return out;
	}

Audio Audio::cut( 
	Second start_time, 
	Second end_time, 
	Second start_fade, 
	Second end_fade 
	) const
	{
	flan_PROCESS_START( Audio() );

	return cut_frames( 
		time_to_frame( start_time	), 
		time_to_frame( end_time 	), 
		time_to_frame( start_fade 	), 
		time_to_frame( end_fade 	)
		);
	}

Audio Audio::cut_frames( 
	Frame start, 
	Frame end, 
	Frame start_fade, 
	Frame end_fade 
	) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	if( end <= start ) return Audio();
	const Frame start_frame = std::clamp( start, 0, get_num_frames() - 1 );
	const Frame end_frame   = std::clamp( end,   0, get_num_frames() - 1 );
	// Any fade errors are validated in Audio::fades
	
	auto format = get_format();
	format.num_frames = end - start;
	Audio out( format );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.get_num_frames() ), [&]( Frame frame )
			{
			out.set_sample( channel, frame, get_sample( channel, start_frame + frame ) );
			} );

	out.fade_frames_in_place( start_fade, end_fade, Interpolators::sqrt );
	return out;
	}

Audio Audio::repitch( const Function<Second, float> & factor, Second granularity_in_seconds, WDLResampleType quality ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	Frame granularity_in_frames = time_to_frame( granularity_in_seconds );
	if( granularity_in_frames < 1 ) granularity_in_frames = 1;

	// This will be inverted momentarily, hence the name
	auto factor_sampled_inverted = factor.sample( 0, std::ceil( get_num_frames() / float( granularity_in_frames ) ), granularity_in_seconds );
	std::for_each( std::execution::par_unseq, factor_sampled_inverted.begin(), factor_sampled_inverted.end(), []( float & sample ){
		constexpr float bound = 1000.0f; // Arbitrary
		sample = std::clamp( 1.0f / sample, 1.0f / bound, bound ); 
		} );

	// Estimate output length
	const Frame num_out_frames = std::ceil( std::accumulate( factor_sampled_inverted.begin(), factor_sampled_inverted.end(), 0.0f ) * granularity_in_frames );

	auto format = get_format();
	format.num_frames = num_out_frames;
	Audio out( format );
	
	WDL_Resampler rs;
		 if( quality == WDLResampleType::Sinc ) 			rs.SetMode( true,  0, true, 64 ); // reasonable sinc resampling
	else if( quality == WDLResampleType::Linear ) 			rs.SetMode( true,  1, false    ); // linear interpolation with simple filter
	else if( quality == WDLResampleType::Uninterpolated ) 	rs.SetMode( false, 0, false    ); // no interpolation, useful for dirty sound

	const Channel num_channels = get_num_channels();
	const Frame in_frames = get_num_frames();
	std::vector<Sample> rsoutbuf( num_channels * granularity_in_frames );

	Frame in_frame = 0;
	Frame out_frame = 0;
	while( in_frame < get_num_frames() )
		{
		const double factor_to_use = factor_sampled_inverted[ std::floor( in_frame / float( granularity_in_frames ) ) ];
		rs.SetRates( get_sample_rate(), double( get_sample_rate() ) * factor_to_use );
		WDL_ResampleSample* rsinbuf = nullptr;
		const Frame wanted = rs.ResamplePrepare( granularity_in_frames, num_channels, &rsinbuf );

		for( Channel channel = 0; channel < num_channels; ++channel )
			for( Frame frame = 0; frame < wanted; ++frame )
				{
				if( in_frame + frame < in_frames )
					rsinbuf[ frame * num_channels + channel ] = get_sample( channel, in_frame + frame );
				else
					rsinbuf[ frame * num_channels + channel ] = 0.0;
				}

		// Process
		rs.ResampleOut( rsoutbuf.data(), wanted, granularity_in_frames, num_channels );

		// Copy resampled data to output buffer
		for( Channel channel = 0; channel < num_channels; ++channel )
			for( Frame frame = 0; frame < granularity_in_frames; ++frame )
				if( out_frame + frame < num_out_frames )
					out.set_sample( channel, out_frame + frame, rsoutbuf[ frame * num_channels + channel ] );

		out_frame += granularity_in_frames;
		in_frame += wanted;
		}

	return out;
	}

Audio Audio::iterate( 
	uint32_t n, 
	const AudioMod & mod, 
	bool feedback 
	) const
	{
	flan_PROCESS_START( Audio() );
	// Half a length is removed to avoid length rounding changing the number of events generated
	return synthesize_grains_with_feedback_mod( get_length() * ( float( n ) - 0.5f ), 1.0f / get_length(), 0, mod, feedback, get_sample_rate() );
	//return synthesize_grains_from_entire_audio( get_length() * ( float( n ) - 0.5f ), 1.0f / get_length(), 0, mod, feedback );
	}

Audio Audio::delay( 
	Second length, 
	const Function<Second, Second> & delay_time, 
	const Function<Second, float> & decay, 
	const AudioMod & mod 
	) const
	{
	flan_PROCESS_START( Audio() );

	if( length <= 0 ) return Audio();

	auto delay_time_sampled = delay_time.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );
	auto decay_sampled = decay.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );

	auto events_per_second = [&]( float t )
		{ 
		const Second delay_time_c = delay_time_sampled[ time_to_frame( t ) ];
		if( delay_time_c <= 0 ) return 1.0f / float( get_sample_rate() );
		return 1.0f / delay_time_c;
		};

	const AudioMod delay_mod = [&]( Audio & in, Second t )
		{ 
		if( t == 0 ) return;

		if( ! mod.is_null() )
			mod( in, t );

		// Modify gain
		const float decay_t = decay_sampled[ time_to_frame( t ) ];
		std::for_each( in.get_buffer().begin(), in.get_buffer().end(), [decay_t]( Sample & s ){ s *= decay_t; } );
		};

	return synthesize_grains_with_feedback_mod( length, events_per_second, 0, delay_mod, true, get_sample_rate() );
	}

// Audio Audio::stereo_delay(
// 	Second length,
// 	const Function<Second, Second> & l_time,
// 	const Function<Second, Second> & r_time,
// 	const Function<Second, Amplitude> & decay
// 	) const
// 	{
// 	auto safe_get_sample = [&]( Channel c, Frame f )
// 		{ 
// 		if( f < 0 || get_num_frames() <= f ) return 0.0f;
// 		else return get_sample( c, f );
// 		};

// 	if( get_num_channels() != 2 ) return Audio();

// 	Audio out = Audio::create_empty_with_frames( time_to_frame( length ), 2, get_sample_rate() );

// 	auto l_time_sampled = out.sample_function_over_domain( l_time );
// 	auto r_time_sampled = out.sample_function_over_domain( r_time );
// 	auto decay_sampled = out.sample_function_over_domain( decay );

// 	const Frame l_buffer_size = time_to_frame( *std::max_element( l_time_sampled.begin(), l_time_sampled.end() ) );
// 	const Frame r_buffer_size = time_to_frame( *std::max_element( r_time_sampled.begin(), r_time_sampled.end() ) );
// 	if( l_buffer_size <= 0 || r_buffer_size <= 0 ) return Audio();
// 	std::vector<Sample> l_buffer( l_buffer_size, 0 );
// 	std::vector<Sample> r_buffer( r_buffer_size, 0 );
	
// 	for( Frame out_frame = 0; out_frame < out.get_num_frames(); ++out_frame )
// 		{
// 		const Frame l_follow_frame = ( out_frame + l_buffer_size - Frame( time_to_frame( l_time_sampled[out_frame] ) ) ) % l_buffer_size;
// 		const Frame r_follow_frame = ( out_frame + r_buffer_size - Frame( time_to_frame( r_time_sampled[out_frame] ) ) ) % r_buffer_size;
// 		const Frame l_lead_frame = out_frame % l_buffer_size;
// 		const Frame r_lead_frame = out_frame % r_buffer_size;

// 		// Read from buffers into output
// 		out.get_sample( 0, out_frame ) = l_buffer[l_lead_frame];
// 		out.get_sample( 1, out_frame ) = r_buffer[r_lead_frame];

// 		// Write lead input and follow cross buffer to buffer
// 		l_buffer[l_lead_frame] = safe_get_sample( 0, out_frame ) + r_buffer[r_follow_frame] * decay_sampled[out_frame];
// 		r_buffer[r_lead_frame] = safe_get_sample( 1, out_frame ) + l_buffer[l_follow_frame] * decay_sampled[out_frame];
// 		}

// 	return out;
// 	}

std::vector<Audio> Audio::chop( 
	Second slice_length, 
	Second fade 
	) const
	{
	flan_PROCESS_START( std::vector<Audio>() );

	// Input validation
	if( slice_length <= 0 ) return std::vector<Audio>();

	// Get all but last slice
	std::vector<Audio> outs;
	Second slice_position = 0;
	while( slice_position + slice_length < get_length() )
		{
		outs.push_back( cut( slice_position, slice_position + slice_length, fade, fade ) );
		slice_position += slice_length;
		}

	return outs;
	}

Audio Audio::rearrange( 
	Second slice_length, 
	Second fade 
	) const
	{
	flan_PROCESS_START( Audio() );

	// + fade here accounts for + fade / 2 at both ends
	std::vector<Audio> chops = chop( slice_length + fade, fade );
	if( chops.size() < 2 )
		return Audio();

	chops.pop_back(); // Final slice usually isn't the correct length

	std::random_device rd;
	std::mt19937 g( rd() );
	std::shuffle( chops.begin(), chops.end(), g );

	return Audio::join( std::move( chops ), -fade );
	}
