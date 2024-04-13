#include "flan/Audio/Audio.h"

#include <ranges>

#include "WDL/resample.h"

using namespace flan;
using namespace std::ranges;

template<typename T>
static std::vector<const T *> get_pointers( const std::vector<T> & ts )
	{
	std::vector<const T *> ptrs( ts.size() );
	std::transform( ts.begin(), ts.end(), ptrs.begin(), []( const T & t ){ return &t; } );
	return std::move( ptrs );
	}

Audio Audio::modify_boundaries( 
	Second start, 
	Second end 
	) const
	{
	if( is_null() ) return Audio::create_null();

	const Frame num_out_frames = time_to_frame( -start ) + get_num_frames() + time_to_frame( end );
	if( num_out_frames <= 0 ) return Audio::create_null();

	Audio::Format format = get_format();
	format.num_frames = num_out_frames;
	Audio out( format );
	out.clear_buffer();

	out.mix_in_place( *this, -start );

	return out;
	}

Audio Audio::remove_edge_silence( 
	Amplitude non_silent_level,
	Second fade_in_time
	) const
	{
	Frame start_frame = 0;
	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		{
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			if( get_sample( channel, frame ) > non_silent_level )
				goto exit_start_loop;
		++start_frame;	
		}
	exit_start_loop:

	Frame end_frame = get_num_frames();
	for( Frame frame = get_num_frames() - 1; frame >= 0; --frame )
		{
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			if( get_sample( channel, frame ) > non_silent_level )
				goto exit_end_loop;
		--end_frame;	
		}
	exit_end_loop:

	const Frame fade_frames = time_to_frame( fade_in_time );
	const Frame start_fade_frames = start_frame - fade_frames < 0? start_frame : fade_frames;
	const Frame end_fade_frames = end_frame + fade_frames >= get_num_frames()? get_num_frames() - end_frame : fade_frames;
	return cut_frames( start_frame - fade_frames, end_frame + fade_frames, start_fade_frames, end_fade_frames );
	}

Audio Audio::remove_silence(
	Amplitude non_silent_level,
	Second minimum_gap,
	Second fade_in_time
	) const
	{
	const Frame gap_frames = time_to_frame( minimum_gap );
	const Frame crossfade_frames_2 = time_to_frame( fade_in_time/2 );

	Frame noisy_start = std::numeric_limits<Frame>::min();
	Frame last_noisy_frame = 0;
	bool in_quiet_sector = true;
	std::vector<std::array<Frame, 2>> noisy_chunk_frames;
	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		{
		// Check for noisy frame
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			{
			if( get_sample( channel, frame ) > non_silent_level )
				{
				last_noisy_frame = frame;
				if( in_quiet_sector )
					{
					noisy_start = frame;
					in_quiet_sector = false;
					}
				break;
				}
			}

		if( !in_quiet_sector && frame - last_noisy_frame > gap_frames ) // If we hit a gap
			{
			in_quiet_sector = true;
			noisy_chunk_frames.push_back( { noisy_start, last_noisy_frame } );
			}
		}
	if( !in_quiet_sector )
		noisy_chunk_frames.push_back( { noisy_start, get_num_frames() } );

	if( noisy_chunk_frames.empty() ) return Audio::create_null();

	std::vector<Second> fade_ins( noisy_chunk_frames.size() + 1, fade_in_time );

	// When asking for fade-ins on chunks, we need to make sure none of them are trying to reach beyond the barrier
	for( int i = 0; i < noisy_chunk_frames.size(); ++i )
		{
		const Frame fade_frames_i = time_to_frame( fade_ins[i] );
		fade_ins[i] = frame_to_time( noisy_chunk_frames[i][0] - fade_frames_i < 0? 
			noisy_chunk_frames[i][0]: 
			fade_frames_i );
		fade_ins[i+1] = frame_to_time( noisy_chunk_frames[i][1] + fade_frames_i >= get_num_frames()? 
			get_num_frames() - noisy_chunk_frames[i][1]: 
			fade_frames_i );
		}

	// Do the actual cutting
	std::vector<Audio> noisy_chunks;
	for( int i = 0; i < noisy_chunk_frames.size(); ++i )
		{
		const Frame left_fade_frames  = time_to_frame( fade_ins[i  ] );
		const Frame right_fade_frames = time_to_frame( fade_ins[i+1] );
		noisy_chunks.push_back( cut_frames( 
			noisy_chunk_frames[i][0] - left_fade_frames, 
			noisy_chunk_frames[i][1] + right_fade_frames,
			left_fade_frames,
			right_fade_frames ) );
		}

	std::vector<Second> offsets( fade_ins.size() );
	for( int i = 0; i < fade_ins.size(); ++i )
		offsets[i] = -2*fade_ins[i];
	return Audio::join( get_pointers( noisy_chunks ), offsets );
	}

Audio Audio::reverse(
	) const
	{
	if( is_null() ) return Audio::create_null();

	Audio out( get_format() );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		const auto iter = get_buffer().begin() + channel * get_num_frames();
		std::copy( FLAN_PAR_UNSEQ iter, iter + get_num_frames(), out.get_buffer().rbegin() + channel * get_num_frames() );
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
	if( is_null() ) return Audio::create_null();

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
	if( is_null() ) return Audio::create_null();

	// Input validation
	if( end <= start ) return Audio::create_null();
	start = std::clamp( start, 0, get_num_frames() - 1 );
	end   = std::clamp( end,   0, get_num_frames() - 1 );
	// Any fade errors are validated in Audio::fades
	
	auto format = get_format();
	format.num_frames = end - start;
	Audio out( format );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( out.get_num_frames() ), [&]( Frame frame )
			{
			out.set_sample( channel, frame, get_sample( channel, start + frame ) );
			} );

	out.fade_frames_in_place( start_fade, end_fade, Interpolators::sqrt );
	return out;
	}

Audio Audio::repitch( const Function<Second, float> & factor, Second granularity_in_seconds, WDLResampleType quality ) const
	{
	if( is_null() ) return Audio::create_null();

	// Input validation
	Frame granularity_in_frames = time_to_frame( granularity_in_seconds );
	if( granularity_in_frames < 1 ) granularity_in_frames = 1;

	// This will be inverted momentarily, hence the name
	auto factor_sampled_inverted = factor.sample( 0, std::ceil( get_num_frames() / float( granularity_in_frames ) ), granularity_in_seconds );
	factor_sampled_inverted.for_each( []( float & sample ){
		constexpr float bound = 1000.0f; // Arbitrary
		sample = std::clamp( 1.0f / sample, 1.0f / bound, bound ); 
		} );

	// Estimate output length
	const Frame num_out_frames = std::ceil( factor_sampled_inverted.accumulate() * granularity_in_frames );

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
	int32_t n, 
	Second crossfade_time,
	const AudioMod & mod, 
	bool feedback 
	) const
	{
	if( is_null() ) return Audio::create_null();
	if( n < 1 ) return Audio::create_null();

	// Without a mod, we are just repeating the input, which join can do already
	if( mod.is_null() )
		return Audio::join( std::vector<const Audio *>( n, this ), -crossfade_time );

	// If mod changes input lengths, the time passed to mod will be wrong. The alternative is always executing in linear_sequence.
	std::vector<Audio> this_modifieds( n );
	flan::for_each_i( n, feedback? ExecutionPolicy::Linear_Sequenced : mod.get_execution_policy(), [&]( int i )
		{
		this_modifieds[i] = i > 0 && feedback ? this_modifieds[i-1].copy() : copy();
		mod( this_modifieds[i], i*get_length() );
		} );

	return Audio::join( this_modifieds, -crossfade_time );
	}

Audio Audio::delay( 
	Second length, 
	const Function<Second, Second> & delay_time, 
	const Function<Second, float> & decay, 
	const AudioMod & mod 
	) const
	{
	if( is_null() ) return Audio::create_null();

	if( length <= 0 ) return Audio::create_null();

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

	return texture( length, events_per_second, 0, delay_mod, true );
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

// 	if( get_num_channels() != 2 ) return Audio::create_null();

// 	Audio out = Audio::create_empty_with_frames( time_to_frame( length ), 2, get_sample_rate() );

// 	auto l_time_sampled = out.sample_function_over_domain( l_time );
// 	auto r_time_sampled = out.sample_function_over_domain( r_time );
// 	auto decay_sampled = out.sample_function_over_domain( decay );

// 	const Frame l_buffer_size = time_to_frame( *std::max_element( l_time_sampled.begin(), l_time_sampled.end() ) );
// 	const Frame r_buffer_size = time_to_frame( *std::max_element( r_time_sampled.begin(), r_time_sampled.end() ) );
// 	if( l_buffer_size <= 0 || r_buffer_size <= 0 ) return Audio::create_null();
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

std::vector<Audio> Audio::split_at_times(
	std::vector<Second> split_times,
	Second fade
	) const
	{
	if( is_null() ) return std::vector<Audio>();

	const Frame fade_frames = time_to_frame( fade );

	std::sort( split_times.begin(), split_times.end() );

	std::vector<Frame> split_frames;
	split_frames.push_back( 0 );
	for( Second t : split_times )
		{
		const Frame f = time_to_frame( t );
		if( f <= 0 ) continue;
		if( get_num_frames() <= f ) break;
		split_frames.push_back( f );
		}
	split_frames.push_back( get_num_frames() );

	std::vector<Audio> outs( split_frames.size() - 1 );
	flan::for_each_i( split_frames.size() - 1, ExecutionPolicy::Parallel_Unsequenced, [&]( int i )
		{
		outs[i] = cut_frames( split_frames[i], split_frames[i+1], fade_frames, fade_frames );
		} );

	return outs;
	}

std::vector<Audio> Audio::split_with_lengths(
	std::vector<Second> split_lengths,
	Second fade
	) const
	{
	for( Second & t : split_lengths ) if( t < 0 ) t = 0; 

	// Create split times from split lengths
	std::vector<Second> split_times;
	std::partial_sum( split_lengths.begin(), split_lengths.end(), std::back_inserter( split_times ) );
	return split_at_times( split_times, fade );
	}

std::vector<Audio> Audio::split_with_equal_lengths( 
	Second slice_length, 
	Second fade 
	) const
	{
	if( slice_length <= 0 ) return std::vector<Audio>();
	std::vector<Second> slice_lengths( std::ceil( get_length() / slice_length ), slice_length );
	return split_with_lengths( slice_lengths, fade );
	}

Audio Audio::random_chunks(
	Second length,
	const Function<Second, Second> & chunk_length, // Seconds per chunk
	const Function<Second, Second> & fade,
	const AudioMod & mod
	) const
	{
	if( is_null() || length <= 0 ) return Audio::create_null();

	// Integrate chunk length to find chunk frames
	std::vector<Frame> chunk_start_frames;
	double chunk_accumulator = 1;
	for( Frame out_frame = 0; out_frame < time_to_frame( length ); ++out_frame )
		{
		if( chunk_accumulator >= 1 )
			{
			chunk_start_frames.push_back( out_frame );
			chunk_accumulator = std::fmod( chunk_accumulator, 1.0 );
			}
		
		// The lower chunk limit is arbitrary, but without a limit a small chunk size and large fade size will use massive amounts of memory
		const double seconds_per_chunk_c = std::clamp( chunk_length( frame_to_time( out_frame ) ), frame_to_time( 32 ), get_length() );
		const double chunks_per_frame_c = 1.0 / seconds_per_chunk_c / get_sample_rate();
		chunk_accumulator += chunks_per_frame_c; // Integrate over frames to get chunks
		}
	chunk_start_frames.push_back( time_to_frame( length ) );

	// Convert from frame positions to frame lengths
	std::vector<Frame> chunk_frame_sizes( chunk_start_frames.size() - 1 );
	for( int i = 0; i < chunk_frame_sizes.size(); ++i ) 
		chunk_frame_sizes[i] = chunk_start_frames[i+1] - chunk_start_frames[i];

	// Get crossfade times
	std::vector<Second> crossfade_times;
	for( int i = 0; i < chunk_start_frames.size(); ++i ) 
		crossfade_times.push_back( fade( frame_to_time( chunk_start_frames[i] ) ) );

	// Get random chunks and mod them if needed
	std::random_device rd;
	std::mt19937 rng( rd() );
	std::vector<Audio> chunks;
	for( int i = 0; i < chunk_frame_sizes.size(); ++i ) 
		{
		const Frame desired_frames = chunk_frame_sizes[i] + time_to_frame( ( crossfade_times[i] + crossfade_times[i+1] ) / 2 );
		const Frame start_frame = desired_frames >= get_num_frames() ? 0 : rng() % ( get_num_frames() - desired_frames );
		Audio chunk = cut_frames( 
			start_frame, 
			start_frame + desired_frames,
			time_to_frame( crossfade_times[i] ),
			time_to_frame( crossfade_times[i+1] ) );
		if( !mod.is_null() ) 
			{
			const Second chunk_start_time = frame_to_time( chunk_start_frames[i] );
			mod( chunk, chunk_start_time );
			}
		chunks.push_back( std::move( chunk ) );
		}	

	std::vector<Second> offsets( crossfade_times.size() );
	for( int i = 0; i < crossfade_times.size(); ++i )
		offsets[i] = -crossfade_times[i];
	return Audio::join( get_pointers( chunks ), offsets );
	}
