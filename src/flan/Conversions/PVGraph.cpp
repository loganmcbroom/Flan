#include "flan/PV/PV.h"
#include "flan/Graph.h"

using namespace flan;

Graph PV::convert_to_graph( Rect D, Pixel width, Pixel height, float timeline_scale ) const
	{
	flan_PROCESS_START( Graph( width, height ) );

	if( D.x2() == -1 ) D.d.x2 = get_length();
	if( D.y2() == -1 ) D.r.x2 = get_height();

	// Validation, Conversion
	const Frame start_frame	= std::clamp( int( time_to_frame( D.x1() ) ), 0, get_num_frames() - 1 );
	const Frame end_frame	= std::clamp( int( time_to_frame( D.x2() ) ), 0, get_num_frames() - 1 );
	const Bin start_bin		= std::clamp( int( frequency_to_bin( D.y1() ) ), 0, get_num_bins() - 1 );
	const Bin end_bin		= std::clamp( int( frequency_to_bin( D.y2() ) ), 0, get_num_bins() - 1 );

	// Convert PV data to Value/Hue data
	const float max_mag = get_max_partial_magnitude( start_frame, end_frame, start_bin, end_bin );
	std::vector<Function<vec2, float>> fs;
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		fs.emplace_back( [channel, max_mag, this]( Second x, Frequency y )
			{
			const Magnitude m = get_MF( channel, time_to_frame( x ), frequency_to_bin( y ) ).m;
			return std::sqrt( std::abs( m ) / max_mag ) * log2( 2.0f + y ) / 4.0f; // sqrt brings up dark areas, log scaling brings up high frequencies
			} );

	Graph g( width, height );
	g.add_full_split_view_y( D, get_num_channels() ); 
	if( max_mag != 0 )
		g.draw_spectrograms( fs, { 0, 0, get_length(), get_height() }, 0 );
	
	// Tick drawing
	if( timeline_scale > 0 )
		{
		const float bigTimeJump = std::pow( 4.0f, std::floor( std::log2( D.w() ) / 2 - 0.5f ) );
		g.draw_x_ticks( bigTimeJump / 4.0f, D.y2(), 1.0f, timeline_scale / 2, 0, -1, Color::from_hsv( 0, 0, .6 ), 0 );
		g.draw_x_ticks( bigTimeJump, D.y2(), 1.0f, timeline_scale, 0, -1, Color::from_hsv( 0, 0, 1 ), timeline_scale );
		//g.draw_y_ticks( 2000, 0, 0, 14, -1, Color::from_hsv( 0, 0, 1 ), true );
		}

	return g;
	}

const PV & PV::save_to_bmp( const std::string & fileName, Rect D, Pixel width, Pixel height ) const
	{
	flan_FUNCTION_LOG;
	convert_to_graph( D, width, height ).save_image( fileName );
	return *this;
	}