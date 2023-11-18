#include "flan/Audio/Audio.h"

#include <flan/FFTHelper.h>
#include <algorithm>

#include "flan/Graph.h"
#include "flan/WindowFunctions.h"
#include "bmp/bitmap_image.hpp"

#undef max
#undef min

using namespace flan;

Graph Audio::convert_to_graph( 
	Interval I, 
	Pixel width, 
	Pixel height, 
	Graph::WaveformMode mode, 
	float timeline_scale 
	) const
	{
	if( is_null() ) return Graph( width, height );

	if( I.x2 == -1 ) I.x2 = get_length();
	
	const Frame start_frame = time_to_frame( I.x1 );
	const Frame end_frame = time_to_frame( I.x2 );

	Graph g( width, height );
	g.fill_image( Color::from_hsv( 0, 0, .04 ) );
	g.add_full_split_view_y( I * Interval( -1, 1 ), get_num_channels() );
		
	std::vector<const float *> channel_starts( get_num_channels() );
	for( int i = 0; i < get_num_channels(); ++i )
		channel_starts[i] = get_sample_pointer( i, 0 );
	g.draw_waveforms( channel_starts, get_num_frames(), { 0, -1, get_length(), 1 }, 0, mode );

	if( timeline_scale > 0 )
		{
		const float bigJump = std::pow( 4.0f, std::floor( std::log2( I.w() ) / 2 - 0.5f ) );
		g.draw_x_ticks( bigJump / 4.0f, 1.0, 1, timeline_scale / 2, 0, -1, Color::from_hsv( 0, 0, .6 ), 0 );
		g.draw_x_ticks( bigJump, 1, 1.0, timeline_scale, 0, -1, Color::from_hsv( 0, 0, 1 ), timeline_scale );
		}

	return g;
	}

void Audio::save_to_bmp( const std::string & filename, Interval I, int width, int height ) const
	{
		auto b = convert_to_graph( I, width, height );
	b.save_image( filename );
	}

Graph Audio::convert_to_spectrum_graph(
	Pixel width, 
	Pixel height,
	Frame smoothing_frames
	) const
	{
	// Graph grid
	Graph g( width, height );
	const float spectrum_length_log = std::log2( get_sample_rate() / 2.0f );
	g.add_full_split_view_y( { 4, -0.1, spectrum_length_log, 1.1f }, get_num_channels() );
	g.fill_image( Color::from_hsv( 0, 0, .05 ) );
	g.draw_log_grid_x( 1, 2, Graph::Plane::All, Color::from_hsv( 0, 0, .1 ) );	
	g.draw_linear_grid_y( 0.1, Graph::Plane::All, Color::from_hsv( 0, 0, .1 ) );	
	g.draw_linear_grid_x( 1.0f, Graph::Plane::All, Color::from_hsv( 0, 0, .25 ) );

	FFTHelper fft( power_of_2_container( get_num_frames() ), true, false, false );
	Audio spectrum = Audio::create_empty_with_frames( fft.complex_buffer_size(), get_num_channels(), get_sample_rate() ); // Hijacking Audio as spectrum
	
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// FFT input
		std::fill( fft.real_begin(), fft.real_end(), 0 );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			fft.get_real_buffer()[frame] = get_sample( channel, frame );
		fft.r2c_execute();

		// Convert fft out into magnitudes
		for( Bin bin = 0; bin < fft.complex_buffer_size(); ++bin )
			spectrum.get_sample( channel, bin ) = std::abs( fft.get_complex_buffer()[bin] );
		}

	// Smooth spectrum
	Audio hann_window = Audio::create_empty_with_frames( smoothing_frames, 1, get_sample_rate() );
	for( Frame i = 0; i < hann_window.get_num_frames(); ++i )
		hann_window.get_sample( 0, i ) = Windows::hann( float( i ) / ( smoothing_frames - 1 ) );
	spectrum = spectrum.convolve( hann_window );

	// Scale magnitudes to 1 and bring up low data
	const Magnitude max_mag = spectrum.get_max_sample_magnitude();
	for( auto & m : spectrum.get_buffer() ) 
		m = std::pow( m / max_mag, .5f );

	std::vector<Function<Second, Sample>> channel_funcs;
	for( Channel channel = 0; channel < spectrum.get_num_channels(); ++channel )
		channel_funcs.emplace_back( [&, channel]( Second x ) -> Sample
			{
			const Frequency freq = std::pow( 2.0f, x );
			const Bin bin = freq / ( float( get_sample_rate() ) / float( fft.real_buffer_size() ) );
			if( bin < 0 || spectrum.get_num_frames() <= bin ) return 0.0f;
			return spectrum.get_sample( channel, bin );
			} );
	g.draw_waveforms( channel_funcs, { 4, -1, spectrum_length_log, 1.0f }, 0, Graph::WaveformMode::Direct );

	// for( Channel channel = 0; channel < get_num_channels(); ++channel )
	// 	{
	// 	// Graph spectrum
	// 	g.draw_waveform( [&]( float x )
	// 		{ 
	// 		const Frequency freq = std::pow( 10.0f, x );
	// 		const Bin bin = freq / ( float( get_sample_rate() ) / float( fft.real_buffer_size() ) );
	// 		if( bin < 0 || fft.complex_buffer_size() <= bin ) return 0.0f;
	// 		return magnitudes[channel * fft.complex_buffer_size() + bin];
	// 		}, Interval::R, channel, Color::from_hsv( 360.0f * channel / get_num_channels(), .8, .65  ) );
	// 	}

	g.draw_axes( Graph::Plane::All, Color::White );
	g.draw_x_ticks( 1.0f, 0, 2.0f, 6, 6, Graph::Plane::All, Color::White, 12 );

	return g;
	}

void Audio::save_spectrum_to_bmp( 
	const std::string & filename, 
	Pixel width, 
	Pixel height,
	Frame smoothing_frames
	) const	
	{
		convert_to_spectrum_graph( width, height, smoothing_frames ).save_image( filename );
	}
