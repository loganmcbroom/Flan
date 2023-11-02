#include "Audio.h"

#include "WDL/resample.h"

using namespace flan; 

Audio Audio::pan( const Function<Second, float> & pan_amount ) const
	{
	if( is_null() ) return Audio();

	if( get_num_channels() != 1 && get_num_channels() != 2 )
		return Audio();

	Audio out = get_num_channels() == 1 ? convert_to_stereo() : copy();
	out.pan_in_place( pan_amount );
	return out;
	}

Audio& Audio::pan_in_place( const Function<Second, float> & pan_amount )
	{
	std::cout << __FUNCTION__ << std::endl; 
	if( is_null() ) { std::cout << "Null input" << std::endl; return *this; }

	if( get_num_channels() != 2 ) return *this; 

	auto pan_amount_sampled = sample_function_over_domain( pan_amount );

	for( Channel channel = 0; channel < 2; ++channel )
		{
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
			{
			const float pan = pan_amount_sampled[frame] / 2.0f + 1.0f; // Convert [-1,1] to [0,1]
			const float pan_func_input = channel == 0 ? pan : ( 1.0f - pan );
			const float scale = Interpolators::sine2( pan_func_input );
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

// I'm leaving this comment block here as a reminder that I am out of my mind and need to check my priorities sometimes
// There will be no atmospheric dispersion filters, why would we need that?
	// Frequency dependent level falloff over distance is approximated with a cascade of 16 lowpass filters
	// See https://en.wikibooks.org/wiki/Engineering_Acoustics/Outdoor_Sound_Propagation#Attenuation_by_atmospheric_absorption.5B5.5D_.5B6.5D_.5B7.5D
	// Through a great deal of experimentation I've found the following IIR filter and cutoff function approximate the correct falloff to an acceptable degree.
	// This has nearly no effect on close sounds, and for far sounds ear distance is negligable, so we only compute on head center.
	// auto falloff_cutoff = [&]( Second t )
	// 	{ 
	// 	const vec2 position_t = position_sampled[time_to_frame(t)];
	// 	return std::min( get_sample_rate(), 207365 / std::sqrt( position_t.mag() ) ); 
	// 	};
	// auto falloff_filtered = filter_1pole_repeat( falloff_cutoff, 16 );

Audio Audio::stereo_spatialize_variable( 
	const Function<Second, vec2> & position 
	) const
	{
	if( get_num_channels() != 1 )
		{
		std::cout << "Audio::spacialize only operates on mono inputs." << std::endl;
		return Audio::create_null();
		}

	const float sound_mps = 331.29f;

	// Based on a human head being 18cm from ear to ear
	const vec2 l_ear = { 0.0f,  0.09f };
	const vec2 r_ear = { 0.0f, -0.09f };

	auto position_sampled = sample_function_over_domain( position );

	Audio l_buffer = copy();
	Audio r_buffer = copy();

	// Level falloff over distance
	// Signals in a computer have no level unit in terms of the physical world, so I make the arbitrary choice here that a signal
	// with a level of 1 played 1 meter from an ear will remain at level 1, with an epsilon to avoid infinite gain.
	l_buffer.modify_volume_in_place( [&]( Second t )
		{ 
		const vec2 position_t = position_sampled[time_to_frame(t)];
		const Meter l_dist = ( l_ear - position_t ).mag();
		return 1.0f / ( l_dist + 0.00001 ); 
		} );
	r_buffer.modify_volume_in_place( [&]( Second t )
		{ 
		const vec2 position_t = position_sampled[time_to_frame(t)];
		const Meter r_dist = ( r_ear - position_t ).mag();
		return 1.0f / ( r_dist + 0.00001 ); 
		} );


	// Head ILD
	auto head_ild_func = [&]( Audio & buffer, Radian direct_angle, vec2 ear )
		{
		auto mix = [&]( Second t )
			{ 
			const vec2 position_t = position_sampled[time_to_frame(t)];
			const Radian buffer_angle_away_from_direct = std::arg( static_cast<std::complex<float>>( position_t - ear ) ) - direct_angle;
			return 0.5f + 0.5f * std::cos( buffer_angle_away_from_direct );
			};
		buffer = std::move( buffer
			.filter_1pole_lowpass( 500, 1 )
			.modify_volume_in_place( [&]( Second t ){ return 1.0f - mix(t); } )
			.mix_in_place( buffer, 0, mix ) );
		};
	head_ild_func( l_buffer,  75.0f * pi2 / 360.0f, l_ear );
	head_ild_func( r_buffer, -75.0f * pi2 / 360.0f, r_ear );


	// ITD
	// Assuming sound sources stay under the speed of sound, this process could be handled entirely with one repitch call per channel.
	// The trouble is that a source moving faster than sound towards a reciever will cause sound to arrive backwards, potentially overlapping
	// later audio. Repitch is just a resampler call and isn't set up to handle overlapping outputs in that way. Additionally, if a "position
	// flickering" effect is desired, the source may be moving discontinuously around space, in which case the doppler effect shouldn't be 
	// applied. Some speed, then, must be chosen as the limit between continous movement and a "teleport", and we can conveniently choose the
	// speed of sound as the boundary to handle both issues. Note, the actual limit used is an epsilon below the speed of sound, as moving exactly the
	// speed of sound towards a reciever would require an infinite repitch input.
	// Because moving faster than sound is considered teleportation, a source moving towards a reciever supersonically over multiple frames
	// will cause distortion due to frame rounding.

	// Find any frames on which that source moves faster than sound
	std::vector<Frame> teleport_frames;
	teleport_frames.push_back( 0 ); // The source initially is "teleported" to its starting location
	for( Frame frame = 1; frame < get_num_frames(); ++frame )
		{
		const Meter position_movement_distance = ( position_sampled[frame] - position_sampled[frame-1] ).mag();
		if( position_movement_distance > ( sound_mps - 0.01 ) * frame_to_time( 1 )  ) // If the movement over one frame was faster than sound - epsilon
			teleport_frames.push_back( frame );
		}

	auto head_itd_func = [&]( const Audio & buffer, vec2 ear ) -> Audio
		{
		// Get the change in distance between the source and reciever on each frame
		// While we're at it, find out how long the output buffer will need to be
		Second max_needed_length = 0;
		std::vector<Meter> relative_distance_change( buffer.get_num_frames() );
		Meter previous_distance = ( position_sampled[0] - ear ).mag(); // Assume no movement is happening on frame 0
		for( Frame frame = 0; frame < buffer.get_num_frames(); ++frame )
			{
			const Meter current_distance = ( position_sampled[frame] - ear ).mag();
			relative_distance_change[frame] = current_distance - previous_distance;
			previous_distance = current_distance;
			
			const Second frame_recieve_time = buffer.frame_to_time(frame) + current_distance / sound_mps;
			if( max_needed_length < frame_recieve_time ) max_needed_length = frame_recieve_time;
			}

		Audio::Format format;
		format.num_channels = 1;
		format.num_frames = std::ceil( buffer.time_to_frame( max_needed_length ) );
		format.sample_rate = buffer.get_sample_rate();
		Audio out( format );
		out.clear_buffer();

		// For each teleport
		for( Index teleport_index = 0; teleport_index < teleport_frames.size(); ++teleport_index )
			{
			const Frame teleport_start_frame = teleport_frames[teleport_index];
			const Frame teleport_end_frame = teleport_index == teleport_frames.size() - 1 ? buffer.get_num_frames() : teleport_frames[teleport_index + 1];
			const Second initial_delay = ( position_sampled[teleport_start_frame] - ear ).mag() / sound_mps;
			const Frame initial_delay_frames = buffer.time_to_frame( initial_delay );

			// Get the time stretch required at each frame
			std::vector<double> stretches( teleport_end_frame - teleport_start_frame );
			for( Frame frame = 0; frame < stretches.size(); ++frame )
				stretches[frame] = 1.0 - relative_distance_change[frame] / double( sound_mps ) * buffer.get_sample_rate();

			// The normal repitch would require repeated buffer copies for each section, so we will have to roll our own for this algorithm

			// Estimate output length
			const Frame num_out_frames = std::ceil( std::accumulate( stretches.begin(), stretches.end(), 0.0f ) );
			
			WDL_Resampler rs;
			rs.SetMode( true,  0, true, 32 ); // reasonable sinc resampling

			const Frame in_frames = teleport_end_frame - teleport_start_frame;
			WDL_ResampleSample rs_out_buffer;

			Frame in_frame = teleport_start_frame;
			Frame out_frame = teleport_start_frame + initial_delay_frames;
			while( in_frame < teleport_end_frame )
				{
				const double stretch_t = 1.0f / stretches[in_frame - teleport_start_frame];
				rs.SetRates( get_sample_rate(), get_sample_rate() * stretch_t );
				WDL_ResampleSample* rs_in_buffer = nullptr;
				const Frame wanted = rs.ResamplePrepare( 1, 1, &rs_in_buffer );

				for( Frame frame = 0; frame < wanted; ++frame )
					{
					if( in_frame + frame < in_frames )
						rs_in_buffer[frame] = buffer.get_sample( 0, in_frame + frame );
					else
						rs_in_buffer[frame] = 0.0;
					}

				// Process
				rs.ResampleOut( &rs_out_buffer, wanted, 1, 1 );

				// Copy resampled data to output buffer
				if( out_frame < out.get_num_frames() )
					out.set_sample( 0, out_frame, rs_out_buffer );

				out_frame += 1;
				in_frame += wanted;
				}
			}

		return out;
		};
	const Audio buffer_l_repitched = head_itd_func( l_buffer, l_ear );
	const Audio buffer_r_repitched = head_itd_func( r_buffer, r_ear );

	return Audio::combine_channels( buffer_l_repitched, buffer_r_repitched );
	}

Audio Audio::stereo_spatialize( 
	vec2 position 
	) const
	{
	if( get_num_channels() != 1 )
		{
		std::cout << "Audio::stereo_spatialize only operates on mono inputs." << std::endl;
		return Audio();
		}

	const float sound_mps = 331.29f;

	// Based on a human head being 18cm from ear to ear
	const vec2 l_ear = { 0.0f,  0.09f };
	const vec2 r_ear = { 0.0f, -0.09f };

	const Meter l_dist = ( l_ear - position ).mag();
	const Meter r_dist = ( r_ear - position ).mag();
	if( l_dist < 0.00001 || r_dist < 0.00001 )
		return Audio(); 

	// Frame rounding here can cause around a quarter inch positional descrepency, less than we are concerned with.
	const Frame l_delay = time_to_frame( l_dist / sound_mps );
	const Frame r_delay = time_to_frame( r_dist / sound_mps );

	// ITD
	Audio::Format buffer_format;
	buffer_format.sample_rate = get_sample_rate();
	buffer_format.num_channels = 1;
	buffer_format.num_frames = get_num_frames() + std::max( l_delay, r_delay );
	Audio buffer_l( buffer_format );
	Audio buffer_r( buffer_format );
	for( Frame frame = 0; frame < get_num_frames(); ++frame ) buffer_l.set_sample( 0, frame + l_delay, get_sample( 0, frame ) );
	for( Frame frame = 0; frame < get_num_frames(); ++frame ) buffer_r.set_sample( 0, frame + r_delay, get_sample( 0, frame ) );
	
	// Level falloff over distance
	// Signals in a computer have no level unit in terms of the physical world, so I make the arbitrary choice here that a signal
	// with a level of 1 played 1 meter from an ear will remain at level 1, with an epsilon to avoid infinite gain.
	buffer_l.modify_volume_in_place( 1.0f / ( l_dist + 0.00001 ) );
	buffer_r.modify_volume_in_place( 1.0f / ( r_dist + 0.00001 ) );

	// Head ILD
	auto head_ild_func = [&]( Audio & buffer, Radian direct_angle )
		{
		const Radian buffer_angle_away_from_direct = std::arg( static_cast<std::complex<float>>( position - l_ear ) ) - direct_angle;
		const float mix = 0.5f + 0.5f * std::cos( buffer_angle_away_from_direct );
		buffer = std::move( buffer
			.filter_1pole_lowpass( 500, 1 )
			.modify_volume_in_place( 1.0f - mix )
			.mix_in_place( buffer, 0, mix ) );
		};
	head_ild_func( buffer_l,  75 );
	head_ild_func( buffer_r, -75 );

	return Audio::combine_channels( buffer_l, buffer_r );
	}