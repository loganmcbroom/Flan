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
void frequency_response_1d( std::function< std::vector<Audio> ( const Audio & )> filter );
void frequency_response_1d( std::function< Audio ( const Audio & )> filter );
void frequency_response_2d( Second length, std::function< Audio ( const Audio & )> filter );

const float pi = std::acos( -1.0f );

int main()
	{
	auto synth = Audio::synthesize( Func1x1::saw, 1, []( Second t ){ return 111; } ).setVolume( .5 );
	auto bah = Audio( "Bah.wav" ).convertToMono().setVolume( .9 );

	frequency_response_1d(
		[]( const Audio & a )
			{
			std::vector<Audio> out;
			out.emplace_back( a.filter_butterworth_lowshelf( 2, 256, -12 ) );
			return out;
			}
		);

	return 0;
	}


#include <Windows.h>
#include <Mmsystem.h>
void play( const Audio & toPlay )
	{
	if( !toPlay.save( "temp_file_save.wav" ) )
        return;
	std::cout << "Playing sound ... " << std::endl;
	if( ! PlaySound("temp_file_save.wav", nullptr, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	}

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
	std::cout << t.elapsedMilliseconds() / n << " ms per call" << std::endl;
	}

void frequency_response_1d( std::function< std::vector<Audio> ( const Audio & )> filter )
	{
	const FrameRate sample_rate = 48000;
	const Second length = std::log2( sample_rate / 2 );
	const Audio sine = Audio::synthesize( Func1x1::sine, length, 
		[&]( Second t ){ return std::pow( 2.0f, t ); }, 
		sample_rate );

	const auto filtered = filter( sine );	

	Graph g;
	g.setView( { -1, -18*3, length, 12 } );
	g.fillImage( Color::White );
	g.drawLinearGrid( 1, 1, 0, Color( 200, 200, 200 ) );
	g.drawAxes( 0, Color::Black );
	for( auto & a : filtered )
		{
		const Func1x1 a_f = a.convertToFunction( .2 );
		const Func1x1 a_f_dB = [&]( Frequency freq ) -> Decibel { return 20.0f * log10( a_f( freq ) ); };
		g.drawFunction( a_f_dB, { 0, length }, -1, Color::Black );
		}
	graph( g );
	}

void frequency_response_1d( std::function< Audio ( const Audio & )> filter )
	{
	return frequency_response_1d( [&]( const Audio & a )
		{ 
		std::vector<Audio> out;
		out.emplace_back( filter( a ) );
		return out;
		} );
	}

void frequency_response_2d( Second length, std::function< Audio ( const Audio & )> filter ) 
	{
	Audio::Format format;
	format.numChannels = 1;
	format.numFrames = length * 48000;
	format.sampleRate = 48000;
	Audio white_noise( format );
	for( Frame frame = 0; frame < white_noise.getNumFrames(); ++frame )
		{
		white_noise.getSample( 0, frame ) = ( std::rand() / float( RAND_MAX ) ) * 2.0f - 1.0f;
		}

	auto filtered_pv = filter( white_noise ).convertToPV();
	const float maxMag = filtered_pv.getMaxPartialMagnitude();

	std::vector<Func2x1> fs;
	fs.push_back( [&filtered_pv, maxMag]( Second x, Frequency y )
		{
		const float m = filtered_pv.getMF( 0, filtered_pv.timeToFrame( x ), filtered_pv.frequencyToBin( y ) ).m / maxMag;
		return std::abs( m );
		} );

	Graph g;
	g.setView( { 0, 0, filtered_pv.getLength(), filtered_pv.getHeight() } ); 
	g.drawSpectrograms( fs, { 0, 0, filtered_pv.getLength(), filtered_pv.getHeight() }, 0 );

	graph( g );
	}

