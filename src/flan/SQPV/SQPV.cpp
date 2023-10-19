// #include "SQPV.h"

// #include <iostream>

// #include "flan/Graph.h"

// using namespace flan;

// SQPV::SQPV()
// 	: SQPVBuffer()
// 	{}

// SQPV::SQPV( const Format & f )
// 	: SQPVBuffer( f )
// 	{}

// SQPV::SQPV( SQPVBuffer && other )
// 	: SQPVBuffer( std::move( other ) )
// 	{}

// Graph SQPV::convert_to_graph() const
// 	{
// 	// Temp
// 	Rect D;
// 	Pixel width = 1920;
// 	Pixel height = 1080; 
// 	float timeline_scale = 0;

// 	flan_PROCESS_START( Graph( width, height ) );

// 	D.d.x1 = 0;
// 	D.d.x2 = get_length();
// 	D.r.x1 = getPitchBandwidth().first;
// 	D.r.x2 = getPitchBandwidth().second;

// 	// Validation, Conversion
// 	const Frame start_frame	= std::clamp( int( time_to_frame( D.x1() ) ), 0, get_num_frames() - 1 );
// 	const Frame end_frame	= std::clamp( int( time_to_frame( D.x2() ) ), 0, get_num_frames() - 1 );
// 	const Bin start_bin		= std::clamp( int( pitchToBin( D.y1() ) ), 0, get_num_bins() - 1 );
// 	const Bin end_bin		= std::clamp( int( pitchToBin( D.y2() ) ), 0, get_num_bins() - 1 );

// 	// Convert PV data to Value/Hue data
// 	const Magnitude max_mag = get_max_partial_magnitude( start_frame, end_frame, start_bin, end_bin );
// 	std::vector<Function<vec2, float>> fs;
// 	for( Channel channel = 0; channel < get_num_channels(); ++channel )
// 		fs.emplace_back( [channel, max_mag, this]( Second x, UnsignedPitch y )
// 			{
// 			const float m = getMP( channel, time_to_frame( x ), pitchToBin( y ) ).m;
// 			return std::sqrt( std::abs( m ) / max_mag ); // sqrt brings up dark areas
// 			} );

// 	Graph g( width, height );
// 	g.add_full_split_view_y( D, get_num_channels() ); 
// 	if( max_mag != 0 )
// 		g.draw_spectrograms( fs, D, 0 );
	
// 	// Tick drawing
// 	if( timeline_scale > 0 )
// 		{
// 		const float bigTimeJump = std::pow( 4.0f, std::floor( std::log2( D.w() ) / 2 - 0.5f ) );
// 		g.draw_x_ticks( bigTimeJump / 4.0f, D.y2(), 1.0f, timeline_scale / 2, 0, -1, Color::from_hsv( 0, 0, .6 ), 0 );
// 		g.draw_x_ticks( bigTimeJump, D.y2(), 1.0f, timeline_scale, 0, -1, Color::from_hsv( 0, 0, 1 ), timeline_scale );
// 		//g.draw_y_ticks( 2000, 0, 0, 14, -1, Color::from_hsv( 0, 0, 1 ), true );
// 		}

// 	return g;
// 	}

// SQPV SQPV::modify_pitch( const Function<std::pair<Second, Pitch>, Pitch> & mod ) const
// 	{
// 	flan_PROCESS_START( SQPV() );

// 	SQPV out = copy();

// 	for( Channel channel = 0; channel < get_num_channels(); ++channel )
// 		{
// 		runtime_execution_policy_handler( mod.get_execution_policy(), [&]( auto policy ){
// 		std::for_each( policy, iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
// 			{
// 			for( Bin bin = 0; bin < get_num_bins(); ++bin )
// 				{
// 				//out.getMP( channel, frame, bin ).p = mod( frame_to_time( frame ), out.getMP( channel, frame, bin ).p );
// 				}
// 			} ); } );
// 		}

// 	return out;
// 	}

// SQPV SQPV::select( Second length, const Function<std::pair<Second, Pitch>, std::pair<Second, Pitch>> & selector ) const
// 	{
// 	flan_PROCESS_START( SQPV() );

// 	// Input validation
// 	if( length <= 0 ) return SQPV();

// 	auto format = get_format();
// 	format.num_frames = time_to_frame( length );
// 	SQPV out( format );

// 	// selector_sampled contains Second and Pitch
// 	const auto selector_sampled = selector.sample( 0, out.get_num_frames(), frame_to_time( 1 ), 
// 		getPitchBandwidth().first * getBinsPerOctave(), 
// 		getPitchBandwidth().second * getBinsPerOctave(), 
// 		1.0f / getBinsPerOctave() );

// 	for( Channel channel = 0; channel < out.get_num_channels(); ++channel )
// 		{
// 		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.get_num_frames() ), [&]( Frame frame )
// 			{
// 			for( Bin bin = 0; bin < out.get_num_bins(); ++bin )
// 				{
// 				const vec2 s = selector_sampled[ out.get_buffer_pos( 0, frame, bin ) ];
// 				const fFrame selectedFrame = time_to_frame( s.x() );
// 				const UnsignedPitch selectedPitch = s.y();
// 				const fBin selectedBin = pitchToBin( selectedPitch );

// 				if( selectedFrame < 0 || get_num_frames() - 1 <= selectedFrame || bin < 0 || get_num_bins() - 1 <= bin )
// 				 	out.getMP( channel, frame, bin ) = { 0, 0 };
// 				else
// 					{	
// 					const fFrame mix = selectedFrame - std::floor( selectedFrame );
// 					const MP left_mp = getMP( channel, std::floor( selectedFrame ), bin );
// 					const MP right_mp = getMP( channel, std::floor( selectedFrame ) + 1, bin );

// 					const Magnitude w1 = ( 1.0f - mix ) * left_mp.m;
// 					const Magnitude w2 = mix * right_mp.m;

// 					MP assign_mp = left_mp;
// 					assign_mp.m = w1 + w2;
// 					assign_mp.p = w1 > w2 ? left_mp.p : right_mp.p;
					
// 					// const Pitch pitch_shift = binToPitch( bin ) - selectedPitch;
// 					// assign_mp.p += pitch_shift;
// 					out.getMP( channel, frame, bin ) = assign_mp;
// 					}
// 				}
// 			} );
// 		}

// 	return out;
// 	}
