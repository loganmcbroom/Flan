#include "Audio.h"

#include "WDL/resample.h"

using namespace flan; 

static const float sound_mps = 343; // Sound speed in air at 68 degrees F

Audio Audio::pan( const Function<Second, float> & pan_amount ) const
	{
	if( is_null() ) return Audio::create_null();

	if( get_num_channels() != 1 && get_num_channels() != 2 )
		return Audio::create_null();

	Audio out = get_num_channels() == 1 ? convert_to_stereo() : copy();
	out.pan_in_place( pan_amount );
	return out;
	}

Audio& Audio::pan_in_place( const Function<Second, float> & pan_amount )
	{
	if( is_null() ) { std::cout << "Null input" << std::endl; return *this; }

	if( get_num_channels() != 2 ) return *this; 

	const auto pan_amount_sampled = sample_function_over_domain( pan_amount );

	for( Channel channel = 0; channel < 2; ++channel )
		{
		flan::for_each_i( get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame frame )
			{
			const float pan = pan_amount_sampled[frame] / 2.0f + 0.5f; // Convert [-1,1] to [0,1]
			const float pan_func_input = channel == 0 ? pan : ( 1.0f - pan );
			const float scale = Interpolator::sine2()( pan_func_input );
			get_sample( channel, frame ) *= scale;
			} );
		}
	return *this;
	}

Audio Audio::widen( const Function<Second, float> & widen_amount ) const
	{
	return convert_to_mid_side().pan( widen_amount ).convert_to_left_right();
	}

/*
The spacialization functions here are largely inaccurate to measured expirimental data. While Flan usually aims to generate the highest quality 
sound, spatialization functions have nearly unlimited room for improvements. For a reasonably accurate continuous model of (1), even computing 
the filter curve requires ~20 normal distribution calls, which then would need to be then converted into a filter. This could of course be simplified 
and optimized, but there are a number of additional issues with using the sampled data, such as it not stretching past 12kH. In (2), a more accurate 
model of ITD is given and tested, but using that model would make the ITD computations for moving sound sources significantly more complex, especially 
when taking the doppler effect into account.

The strategy used here for ITD is simply to compute the linear distance from the sound source to each ear, find the time it would take for sound to 
travel that distance through air, and delay the monophonic input by that amount at each frame. Using repitch to manage the delays handles timing, the 
doppler effect, and any concerns of artifacts caused by rounding to the nearest frame. An additional 16-pole low pass iir filter is applied to handle 
high frequencies falling off faster over large distances. That filter is an approximation of the response given in (3).

The strategy for ILD is a combination of data taken from (1) and my own expirements with simplifying that data to a simple monotonic filter. Sound 
sources pointing directly into an ear (~75 degrees azimuthal angle) are unaffected, and sounds opposite those are maximally filtered, with intermediate 
angles mixing these signals	sinusoidally as a function of angle. The maximal filter is a 1 pole lowpass iir filter with a cutoff of 500hZ.

(1) Transformation of sound-pressure level from the free field to the eardrum presented in numerical form, Shaw and Vaillancourt, 1984
(2) Model for the interaural time differences in the azimuthal plane, G. Kuhn, 1976
(3) https://en.wikibooks.org/wiki/Engineering_Acoustics/Outdoor_Sound_Propagation#Attenuation_by_atmospheric_absorption.5B5.5D_.5B6.5D_.5B7.5D
*/

// Audio Audio::filter_pinna( const Function<Second, Meter> & height ) const
// 	{
// 	auto elevation_angles = sample_function_over_domain( height )
// 		.transform( []( Meter m ){ return std::atan( m ); } );

// 	auto main_gain  = [&]( Second t ){ return -5.0f + elevation_angles[time_to_frame(t)]/(pi/2) * 10.0f; };
// 	auto thin_gain  = [&]( Second t ){ return main_gain(t) * 0.8f; };
// 	auto broad_gain = [&]( Second t ){ return main_gain(t) * 0.1f; };

// 	// Damping factors are from qr=1/2
// 	return
// 		 filter_2pole_bandshelf( 8000, 0.25f, main_gain )  // Main shelf
// 		.filter_2pole_bandshelf( 10000, 0.03f, thin_gain )  // Thin shelf
// 		.filter_2pole_bandshelf( 3500, 0.7f, broad_gain ); // Broad shelf
// 	}

// Frequency dependent level falloff over distance is approximated with a cascade of 16 lowpass filters
// See https://en.wikibooks.org/wiki/Engineering_Acoustics/Outdoor_Sound_Propagation#Attenuation_by_atmospheric_absorption.5B5.5D_.5B6.5D_.5B7.5D
// Through a great deal of experimentation I've found the following filter and cutoff function approximate the correct falloff to an acceptable degree.
static Audio atmospheric_dispersion( const Audio & me, const FunctionSample<vec2d> & position_samples )
	{
	auto distance_sqrts_sampled = position_samples.transform(
		[]( vec2d v ) -> Meter { return std::sqrt( v.mag() ); } );

	auto falloff_cutoff = [&]( Second t )
		{ 
		const Meter dr = distance_sqrts_sampled[me.time_to_frame( t )];
		return 207365 / dr; 
		};
	return me.filter_1pole_repeat_low( falloff_cutoff, 16 );
	}

// Level falloff over distance
// Signals in a computer have no level unit in terms of the physical world, so I make the arbitrary choice here that a signal
// with a level of 1 played 1 meter from an ear will remain at level 1, with an epsilon to avoid infinite gain.
static void falloff( Audio & me, const FunctionSample<vec2d> & ps )
	{
	const float epsilon = 0.00001f;
	me.modify_volume_in_place( [&]( Second t )
		{ 
		const Frame frame_t = std::floor( me.time_to_frame(t) );
		const vec2d position_t = ps[frame_t];
		const Meter dist = position_t.mag();
		return 1.0f / ( dist + epsilon ); 
		} );
	}

static Audio head_ild( const Audio & me, const FunctionSample<vec2d> & ps, Radian direct_angle )
	{
	auto mix_sample = ps.transform( [&]( vec2d p )
		{
		const vec2 planar_pos = vec2( p.x(), p.y() );
		const Radian angle_away_from_direct = std::arg( static_cast<std::complex<float>>( planar_pos ) ) - direct_angle;
		return 0.5f + 0.5f * std::cos( angle_away_from_direct );
		} );

	auto mix = mix_sample.to_time_function( me.get_sample_rate() );

	return std::move( me
		.filter_1pole_lowpass( 500, 1 )
		.modify_volume_in_place( [&]( Second t ){ return 1.0f - mix(t); } )
		.mix_in_place( me, 0, mix ) );
	};

// ITD
// Note, this process will break if ps moves at or above the speed of sound
static Audio head_itd( const Audio & me, const FunctionSample<vec2d> & positions_sampled )
	{
	const Frame granularity_frames = 32;

	if( positions_sampled.is_constant() )
		{
		auto p = positions_sampled.get_constant();
		const Meter dist = p.mag();
		const Second delay = me.time_to_frame( dist / sound_mps );
		
		Audio::Format format;
		format.sample_rate = me.get_sample_rate();
		format.num_channels = 1;
		format.num_frames = me.get_num_frames() + delay;
		Audio out( format );
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame ) 
			out.set_sample( 0, frame + delay, me.get_sample( 0, frame ) );
		return out;
		}

	auto & ps = positions_sampled.get_vector();

	// Get the change in distance between the source and reciever for each granularity_frames jump
	// While we're at it, find out how long the output buffer will need to be
	Second max_needed_length = 0;
	std::vector<double> relative_distance_change_per_chunk( 1, 0.0 );
	double previous_distance = ps[0].mag(); // Assume no movement is happening on frame 0
	for( Frame frame = granularity_frames; frame < me.get_num_frames(); frame += granularity_frames )
		{
		const double current_distance = ps[frame].mag();
		relative_distance_change_per_chunk.push_back( current_distance - previous_distance );
		previous_distance = current_distance;
		
		const Second max_frame_recieve_time = me.frame_to_time(frame + granularity_frames) + current_distance / sound_mps;
		if( max_needed_length < max_frame_recieve_time ) max_needed_length = max_frame_recieve_time;
		}

	Audio::Format format;
	format.num_channels = 1;
	format.num_frames = std::ceil( me.time_to_frame( max_needed_length ) );
	format.sample_rate = me.get_sample_rate();
	Audio out( format );
	out.clear_buffer();

	const Second initial_delay = ps[0].mag() / sound_mps;
	const Frame initial_delay_frames = me.time_to_frame( initial_delay );

	// Get the time stretch required at each frame
	std::vector<double> stretches( relative_distance_change_per_chunk.size() );
	for( Frame frame = 0; frame < stretches.size(); ++frame )
		stretches[frame] = 1.0 / ( 1.0 - relative_distance_change_per_chunk[frame] / granularity_frames / double( sound_mps ) * me.get_sample_rate() );
	const double max_stretch = *std::max_element( stretches.begin(), stretches.end() );
	
	WDL_Resampler rs;
	rs.SetMode( true,  0, true, 32 ); // Reasonable sinc resampling
	rs.SetFeedMode( true );
	std::vector<Sample> rs_out_buffer( std::ceil( granularity_frames * max_stretch ) );

	Frame out_frame = 0;
	for( Frame in_frame = 0; in_frame < me.get_num_frames(); in_frame += granularity_frames )
		{
		const double stretch_c = stretches[in_frame/granularity_frames];
		rs.SetRates( me.get_sample_rate(), double( me.get_sample_rate() ) * stretch_c );
		WDL_ResampleSample* rsinbuf = nullptr;
		rs.ResamplePrepare( granularity_frames, 1, &rsinbuf );

		for( Frame frame = 0; frame < granularity_frames; ++frame )
			{
			if( in_frame + frame < me.get_num_frames() )
				rsinbuf[ frame ] = me.get_sample( 0, in_frame + frame );
			else
				rsinbuf[ frame ] = 0.0;
			}

		// Process
		const Frame num_out_samples = rs.ResampleOut( rs_out_buffer.data(), granularity_frames, std::ceil( granularity_frames*stretch_c ), 1 );

		// Copy resampled data to output buffer
		for( Frame frame = 0; frame < num_out_samples; ++frame )
			if( out_frame + frame < out.get_num_frames() )
				out.get_sample( 0, out_frame + frame ) = rs_out_buffer[frame];

		out_frame += num_out_samples;
		}

	return out;
	}

Audio Audio::stereo_spatialize(
	const Function<Second, vec2> & position,
	Meter head_width,
	const Function<Second, float> & speed_limit
	) const
	{
	if( get_num_channels() != 1 )
		{
		std::cout << "Audio::spacialize only operates on mono inputs." << std::endl;
		return Audio::create_null();
		}

	auto position_sampled = sample_function_over_domain( position ).convert_to<vec2d>();

	// Position needs to be speed limited to just under the speed of sound
	if( !position.is_constant() )
		{
		// This is an unusually large epsilon. When a sound moves directly away from an ear, the closer it is to the speed
		// of sound, the larger the itd resample output buffer needs to be. At the time of writing this the resample chunk size
		// is 64 frames, around one millisecond. Having too small of an epsilon would allow this output buffer to become overly
		// large, and pointlessly so. Resampling audio beyond even 100 times makes it essentially inaudible, and even with
		// an epsilon of 1, stretches with expansions of above 300 are possible.
		const float epsilon = 1;
		auto speed_limit_sampled = sample_function_over_domain( speed_limit );
		auto & ps = position_sampled.get_vector();
		for( Frame frame = 1; frame < get_num_frames(); ++frame )
			{
			const auto movement = ps[frame] - ps[frame-1];
			const Meter mag = movement.mag();
			// Division by sr converts from meters per second to meters per frame
			const float speed_limit_c = std::clamp( speed_limit_sampled[frame], 0.0f, sound_mps - epsilon ) / get_sample_rate();
			if( mag > speed_limit_c )
				ps[frame] = ps[frame-1] + movement / mag * speed_limit_c;
			}
		}

	auto do_the_whole_thing_for_one_ear = [&]( bool is_left_ear, Radian ear_direction )
		{
		const vec2d ear_position = vec2d( 0, ( is_left_ear? 1 : -1 ) * head_width/2.0f );
	
		const auto relative_positions = position_sampled.transform( [&]( vec2d v ){ return v - ear_position; } );

		//Audio buffer = atmospheric_dispersion( *this, relative_positions );
		Audio buffer = head_ild( *this, relative_positions, ear_direction );
		falloff( buffer, relative_positions );
		buffer = head_itd( buffer, relative_positions ); // This handles doppler

		return buffer;
		};

	Audio l, r;
	flan::for_each_i( 2, ExecutionPolicy::Linear_Sequenced, [&]( int i )
		{
		if( i == 0 ) l = do_the_whole_thing_for_one_ear( true,   75.0f * pi2 / 360.0f );
		else 		 r = do_the_whole_thing_for_one_ear( false, -75.0f * pi2 / 360.0f );
		});

	return Audio::combine_channels( l, r );
	}

