#include <iostream>
#include <algorithm>
#include <numeric>
#include <ranges>

#include "flan/Audio/Audio.h"
#include "flan/PV/PV.h" 
#include "flan/SPV/SPV.h" 
#include "flan/SQPV/SQPV.h"
#include "flan/Function.h"
#include "flan/Wavetable.h"
#include "flan/Synthesis.h"
#include "flan/Graph.h"
#include "flan/WindowFunctions.h"
#include "flan/Utility/Timer.h"
#include "flan/DSPUtility.h"

using namespace flan;

void play( const Audio & toPlay );
void graph( const std::string & f );
void graph( const bitmap_image & b );
template<typename T, typename Callable>
void test( const T & a, int n, Callable c );
void frequency_response_1d( Second length, std::function< std::vector<Audio> ( const Audio & )> filter );
void frequency_response_1d( Second length, std::function< Audio ( const Audio & )> filter );
void frequency_response_2d( Second length, std::function< Audio ( const Audio & )> filter );
void graph_spectrum( const Audio & a );

const float pi = std::acos( -1.0f );

int main()
	{
	auto synth = Audio::synthesize_waveform( waveforms::sine, .5, []( Second t ){ return 111; } ).set_volume( 1 );
	//auto bah = Audio( "Bah.wav" ).set_volume( .9 );

	// play( bah
	// 	.convert_to_PV()
	// 	.select( 5, []( TF tf ){ return TF{ tf.t / 4, std::sin(tf.t*50)*100*tf.t + tf.f }; } )
	// 	.convert_to_audio()
	// 	.set_volume( .7 ) );

	//play( bah.stereo_delay( 5, .5, 1, 0.5 ).set_volume( .5 ) );

	// play( Audio::join(
	// 	//synth,
	// 	synth.add_moisture( 1, []( Second t ){ return 20; }, []( Second t ){ return 2; }, []( Second t ){ return std::sin(pi2*t*t); } )
	// 	).set_volume( .7 )

	// 	);

	//const Second length = 35;
	// play( PV::synthesize( length, 
	// 	[]( Second t ){ return 40 + std::pow( 2.0f, t/4 ); }, 
	// 	[]( std::pair<Second, Harmonic> sh ){ return 1.0f / ( 1.0f + sh.second ); },
	// 	[]( Second t ){ return 20; },
	// 	[]( TF tf ){ return 90 * std::pow( 2.0f, -tf.t ); }
	// 	).convert_to_audio()
	// 	//.modify_volume_in_place( []( Second t ){ return std::sin( std::pow( 2.0f, t / 3 ) ); } )
	// 	.stereo_spatialize( []( Second t ){ return vec2{ std::sin( std::pow( 2.0f, t / 3 ) ), std::cos( std::pow( 2.0f, t / 3 ) ) }; } )
	// 	.compress( -6 )
	// 	.set_volume( .9 ) );

	// graph( synth.set_volume( 1 ).get_amplitude_envelope().convert_to_graph() );

	// auto env = bah.get_amplitude_envelope();

	// play( bah.waveshape( [&]( Second t, Sample s )
	// 	{ 
	// 	const Sample a = s;
	// 	const Sample b = std::abs( s );
	// 	const float mix = env( t );
	// 	return b * mix + a * (1.0f - mix);
	// 	} ) );

	// auto shaper = []( Second t, Sample s ){ return std::clamp( 4.0f*s, -1.0f, 1.0f ); };
	// graph( Audio::combine_channels( 
	// 	synth,
	// 	synth.waveshape( shaper, 1 ),
	// 	synth.waveshape( shaper, 4 ),
	// 	synth.waveshape( shaper, 16 )
	// 	).convert_to_PV().convert_to_graph() );


	return 0;
	}


// #include <Windows.h>
// #include <Mmsystem.h>
// void play( const Audio & toPlay )
// 	{
// 	if( !toPlay.save( "temp_file_save.wav" ) )
//         return;
// 	std::cout << "Playing sound ... " << std::endl;
// 	if( ! PlaySound("temp_file_save.wav", nullptr, SND_FILENAME) )
// 		std::cout << "Error playing sound\n";
// 	}

void graph( const std::string & f )
	{
	system( (std::string("start ") + f).c_str() );
	}

void graph( const bitmap_image & b )
	{
	b.save_image( "Temp.bmp" );
	system( (std::string("start ") + "Temp.bmp").c_str() );
	}

template<typename T, typename Callable>
void test( const T & a, int n, Callable c )
	{
	Timer t;
	t.start();
	for( const int i : std::views::iota( 0, n ) )
		c( a );
	t.stop();
	std::cout << t.elapsed_ms() / n << " ms per call" << std::endl;
	}

void frequency_response_1d( Second length, std::function< std::vector<Audio> ( const Audio & )> filter )
	{
	const FrameRate sample_rate = 48000;
	
	const Audio sine = Audio::synthesize_waveform( waveforms::sine, length, 
		[&]( Second t ){ return std::pow( sample_rate / 2, t / length ); },  
		sample_rate );

	const auto filtered = filter( sine );	

	Graph g;
	g.set_view( { -1, -12, length, 6 } );
	g.fill_image( Color::White );
	g.draw_linear_grid( 1, 1, 0, Color( 200, 200, 200 ) );
	g.draw_axes( 0, Color::Black );
	for( auto & a : filtered )
		{
		const auto a_f = a.get_amplitude_envelope();
		const auto a_f_dB = [&]( Frequency freq ) -> Decibel { return 20.0f * log10( a_f( freq ) ); };
		g.draw_function( a_f_dB, { 0, length }, -1, Color::Black );
		}
	graph( g );
	}

void frequency_response_1d( Second length, std::function< Audio ( const Audio & )> filter )
	{
	return frequency_response_1d( length, [&]( const Audio & a )
		{ 
		std::vector<Audio> out;
		out.emplace_back( filter( a ) );
		return out;
		} );
	}

void frequency_response_2d( Second length, std::function< Audio ( const Audio & )> filter ) 
	{
	Audio::Format format;
	format.num_channels = 1;
	format.num_frames = length * 48000;
	format.sample_rate = 48000;
	Audio white_noise( format );
	for( Frame frame = 0; frame < white_noise.get_num_frames(); ++frame )
		{
		white_noise.get_sample( 0, frame ) = ( std::rand() / float( RAND_MAX ) ) * 2.0f - 1.0f;
		}

	auto filtered_pv = filter( white_noise ).convert_to_PV();
	const float max_mag = filtered_pv.get_max_partial_magnitude();

	std::vector<Function<vec2, Magnitude>> fs;
	fs.push_back( [&filtered_pv, max_mag]( vec2 xy )
		{
		const float m = filtered_pv.get_MF( 0, filtered_pv.time_to_frame( xy.x() ), filtered_pv.frequency_to_bin( xy.y() ) ).m / max_mag;
		return std::abs( m );
		} );

	Graph g;
	g.set_view( { 0, 0, filtered_pv.get_length(), filtered_pv.get_height() } ); 
	g.draw_spectrograms( fs, { 0, 0, filtered_pv.get_length(), filtered_pv.get_height() }, 0 );

	graph( g );
	}

