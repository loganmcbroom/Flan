#include "flan/Graph.h"

#include <algorithm>

#include "flan/Function.h"
#include "flan/Utility/Interpolator.h"
#include "flan/Utility/iota_iter.h"
#include "flan/Utility/execution.h"

using namespace flan;

const Graph::Plane Graph::Plane::All = Graph::Plane( -1 );

Graph::Graph( Pixel _width, Pixel _height )
	: bitmap_image( _width == -1? Graph::default_width : _width, _height == -1 ? Graph::default_height : _height )
	, views()
	{
	}

void Graph::add_view( View view, Plane plane )
	{
	views.emplace_back( plane, view );
	}

void Graph::set_view( Rect rect )
	{
	views.clear();
	views.emplace_back( -1, View( rect, Rect( 0.0f, 0.0f, (float) width(), (float) height() ) ) );
	}

void Graph::add_split_view_y( View view, int num_views, Plane start_plane )
	{
	float y_c = view.V.y1();
	const float viewHeight = view.V.h() / float( num_views );
	for( int v = 0; v < num_views; ++v )
		{
		const View newView( view.U, { view.V.x1(), y_c, view.V.x2(), y_c + viewHeight } );
		add_view( newView, start_plane + v );
		y_c += viewHeight;
		}
	}

void Graph::add_full_split_view_y( Rect rect, int num_views, Plane start_plane )
	{
	add_split_view_y( View( rect, Rect( 0, 0, width(), height() ) ), num_views, start_plane );
	}

bool Graph::do_planes_match( Plane p1, Plane p2 ) const
	{
	if( p1 == Plane::All || p2 == Plane::All ) return true;
	else return p1 == p2;
	}

std::vector<std::pair<Graph::Plane,View>> Graph::get_intersecting_views( Rect U, Plane plane ) const
	{
	std::vector<std::pair<Plane,View>> active_views;
	for( auto & v : views )
		{
		if( do_planes_match( v.first, plane ) && v.second.U.intersect( U ).valid() )
			active_views.emplace_back( v );
		}
	return active_views;
	}

//======================================================================================================================================================
// Waveforms
//======================================================================================================================================================

void Graph::draw_waveform( const Function<float, float> & data, Rect rect, Plane plane, Color c, WaveformMode mode )
	{
	
	const auto active_views = get_intersecting_views( rect, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;
		const Rect drawRect = rect.intersect( view.U );

		const Pixel startPixel = std::ceil(  view.xUToV( drawRect.x1() ) );
		const Pixel endPixel   = std::floor( view.xUToV( drawRect.x2() ) );
		const Pixel yMidPixel = view.yUToV( rect.r.midpoint() );

		// Oversample data into buffer
		//for( Pixel x = startPixel; x < endPixel; ++x )
		runtime_execution_policy_handler( data.get_execution_policy(), [&]( auto policy ){
		std::for_each( FLAN_POLICY iota_iter( startPixel ), iota_iter( endPixel ), [&]( Pixel x )
			{
			//Get average of subsamples
			float sum = 0;
			for( uint32_t sample = 0; sample < oversample; ++sample )
				{
				const float subsample = view.xVToU( float( x ) + float( sample ) / oversample );
				const float data_c = data( subsample );
				sum += mode == WaveformMode::Direct ? data_c : std::abs( data_c );
				}
			sum /= oversample;

			const Pixel yOffset_pixels = view.hUToV( std::clamp( sum, -1.0f, 1.0f ) * rect.h() / 2.0f );

			auto set_pixelWithColor = [this, yMidPixel, &rect, &view, c]( Pixel x, Pixel y )
				{ 
				if( y < view.V.y1() ||  view.V.y2() <= y ) return;
				//const float r = std::abs( y - yMidPixel ) / view.hUToV( rect.h() / 2.0f );
				set_pixel( x, y, c );
				};

			if( mode == WaveformMode::Direct )
				{
				if( yOffset_pixels < 0 )
					for( Pixel y = yMidPixel; y >= yMidPixel + yOffset_pixels; --y ) 
						set_pixelWithColor( x, y );
				else
					for( Pixel y = yMidPixel; y <= yMidPixel + yOffset_pixels; ++y ) 
						set_pixelWithColor( x, y );
				}
			else
				for( Pixel y = yMidPixel - yOffset_pixels; y <= yMidPixel + yOffset_pixels; ++y ) 
					set_pixelWithColor( x, y );
			} ); } );
		}
	}

void Graph::draw_waveform( const float * data, int n, Rect rect, Plane plane, Color c, WaveformMode mode )
	{
	draw_waveform( [data, n, &rect]( float x )
		{ 
		const int i = std::floor( ( x - rect.x1() ) / rect.w() * n );
		//if( i < 0 || i >= n ) return 0.0f;
		return data[i]; 
		}, 
		rect, plane, c, mode );
	}

void Graph::draw_waveforms( const std::vector<Function<float, float>> & fs, Rect rect, Plane start_plane, WaveformMode mode )
	{
	for( int f = 0; f < fs.size(); ++f )
		{
		const Color c = Color::from_hsv( 360.0f * f / fs.size(), .8, .65 );
		draw_waveform( fs[f], rect, start_plane + f, c, mode );
		}
	}

void Graph::draw_waveforms( const std::vector<const float *> & fs, int n, Rect rect, Plane start_plane, WaveformMode mode )
	{
	for( int f = 0; f < fs.size(); ++f )
		{
		const Color c = Color::from_hsv( 360.0f * f / fs.size(), .8, .65 );
		draw_waveform( fs[f], n, rect, start_plane + f, c, mode );
		}
	}


//======================================================================================================================================================
// Spectrograms
//======================================================================================================================================================

void Graph::draw_spectrogram( const Function<vec2, float> & data, Rect rect, Plane plane, float hue )
	{
	
	const int oversample_c = std::sqrt( oversample );

	const auto active_views = get_intersecting_views( rect, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;
		const Rect drawRect = rect.intersect( view.U );

		const Pixel startXPixel = std::ceil(  view.xUToV( drawRect.x1() ) );
		const Pixel endXPixel	= std::floor( view.xUToV( drawRect.x2() ) );
		const Pixel startYPixel = std::ceil(  view.yUToV( drawRect.y1() ) );
		const Pixel endYPixel	= std::floor( view.yUToV( drawRect.y2() ) );

		runtime_execution_policy_handler( data.get_execution_policy(), [&]( auto policy ){
		std::for_each( FLAN_POLICY iota_iter( startXPixel ), iota_iter( endXPixel ), [&]( Pixel x )
			{
			for( Pixel y = startYPixel; y < endYPixel; ++y )
				{
				// Average data over pixel
				float mag = 0;
				for( uint32_t ySample = 0; ySample < oversample_c; ++ySample )
					for( uint32_t xSample = 0; xSample < oversample_c; ++xSample )
						{
						const float xSubsample = view.xVToU( x + float( xSample ) / oversample_c );
						const float ySubsample = view.yVToU( y + float( ySample ) / oversample_c );
						mag += data( xSubsample, ySubsample );
						}
				mag /= ( oversample_c * oversample_c ); // Mag is in range [0, 1] after division

				const Color color = Color::from_hsv( hue, 1.0f, std::clamp( mag, 0.0f, 1.0f ) );
				set_pixel( x, y, color );
				}
			} ); } );
		}
	}

void Graph::draw_spectrogram( const float * data, int n, int m, Rect rect, Plane plane, float hue )
	{
	Function<vec2, float> func = [&data, n, m, &rect]( vec2 xy ) -> float
		{ 
		const int i = std::floor( ( xy.x() - rect.x1() ) / rect.w() * n );
		const int j = std::floor( ( xy.y() - rect.y1() ) / rect.h() * m );
		//if( i < 0 || i >= n || j < 0 || j >= m ) return 0.0f;
		return data[ i * m + j ]; 
		};

	draw_spectrogram( func, rect, plane, hue );
	}

void Graph::draw_spectrograms( const std::vector<Function<vec2, float>> & fs, Rect rect, Plane start_plane )
	{
	for( int f = 0; f < fs.size(); ++f )
		{
		const float hue = 360.0f * f / fs.size();
		draw_spectrogram( fs[f], rect, start_plane + f, hue );
		}
	}

void Graph::draw_spectrograms( const std::vector<const float *> & fs, int n, int m, Rect rect, Plane start_plane )
	{
	for( int f = 0; f < fs.size(); ++f )
		{
		const float hue = 360.0f * f / fs.size();
		draw_spectrogram( fs[f], n, m, rect, start_plane + f, hue );
		}
	}


//======================================================================================================================================================
// Functions
//======================================================================================================================================================

void Graph::draw_function( const Function<float, float> & f, Interval domain, Plane plane, Color c )
	{
	
	//Draw function
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto active_views = get_intersecting_views( domain * Interval::R, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;
		const Interval drawDomain = domain.intersect( view.U.d );

		const float pixelAdvance = view.wVToU( 1 );

		float x1 = drawDomain.x1;
		float y1 = f( x1 );
		for( float x2 = drawDomain.x1 + pixelAdvance; x2 < drawDomain.x2; x2 += pixelAdvance )	
			{
			const float y2 = f( x2 );
			if( view.U.contains( x1, y1 ) && view.U.contains( x2, y2 ) ) 
				draw_line_segment( view, draw, x1, y1, x2, y2 );
			x1 = x2;
			y1 = y2;
			}
		}
	}

void Graph::draw_function( const std::vector<std::pair<float,float>> & data, Plane plane, Color c )
	{
	auto f = interpolate_points( data );

	auto xComp = []( const std::pair<float,float> & p,  const std::pair<float,float> & s ){ return p.first < s.first; };
	const float left	= std::min_element( data.begin(), data.end(), xComp )->first;
	const float right	= std::max_element( data.begin(), data.end(), xComp )->first;
	
	draw_function( f, Interval( left, right ), plane, c );
	}

void Graph::draw_functions( const std::vector<Function<float, float>> & fs, const std::vector<Interval> & domains, Plane plane )
	{

	// Draw waveforms
	for( int f = 0; f < fs.size(); ++f )
		{
		const float hue = 360.0f * f / fs.size();
		draw_function( fs[f], f < domains.size() ? domains[f] : Interval::R, plane, Color::from_hsv( hue, 1, 1 ) );
		}
	}


//======================================================================================================================================================
// Primitives
//======================================================================================================================================================

void Graph::set_pixel( Pixel x, Pixel y, Color c )
	{
	bitmap_image::set_pixel( x, height() - 1 - y, c );
	}

void Graph::set_point( const View & v, float x, float y, Color c )
	{
	set_pixel( v.xUToV( x ), v.yUToV( y ), c );
	}

void Graph::draw_horizontal_line( const View & v, image_drawer & draw, float x1, float x2, float y )
	{
	draw.horiztonal_line_segment( v.xUToV( x1 ), v.xUToV( x2 ), height() - 1 - v.yUToV( y ) );
	}

void Graph::draw_vertical_line( const View & v, image_drawer & draw, float y1, float y2, float x )
	{
	draw.vertical_line_segment( height() - 1 - v.yUToV( y2 ), height() - 1 - v.yUToV( y1 ), v.xUToV( x ) );
	}

void Graph::draw_line_segment( const View & v, image_drawer & draw, float x1, float y1, float x2, float y2 )	
	{
	draw.line_segment( v.xUToV( x1 ), height() - 1 - v.yUToV( y1 ), v.xUToV( x2 ), height() - 1 - v.yUToV( y2 ) );
	}

void Graph::set_rect( const View & v, Rect rect, Color c )
	{
	const Rect drawRect = rect.intersect( v.U );
	set_region( 
		v.xUToV( drawRect.x1() ), 
		height() - v.yUToV( drawRect.y2() ), // +1 handles an edge case created by bmp plane flip
		v.wUToV( drawRect.w() ), 
		v.wUToV( drawRect.h() ), 
		c.red, c.green, c.blue );
	}

void Graph::fill_image( Color c )
	{
	set_region( 0, 0, width(), height(), c.red, c.green, c.blue );
	}


//======================================================================================================================================================
// Additional Drawing Routines
//======================================================================================================================================================

void Graph::draw_axes( Plane plane, Color c )
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;
		if( view.U.x1() <= 0 && view.U.x2() > 0 ) 
			draw_vertical_line( view, draw, view.U.y1(), view.U.y2(), 0 );
		if( view.U.y1() <= 0 && view.U.y2() > 0 ) 
			draw_horizontal_line( view, draw, view.U.x1(), view.U.x2(), 0 );
		}
	}

void Graph::draw_linear_grid_x( 
	float x_jump_size, 
	Plane plane, 
	Color c 
	)
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	if( x_jump_size <= 0 ) return;

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;		
		const float x_start = std::ceil(  view.U.x1() / x_jump_size ) * x_jump_size;
		const float x_end   = std::floor( view.U.x2() / x_jump_size ) * x_jump_size;
		for( float x = x_start; x <= x_end; x += x_jump_size ) 
			draw_vertical_line( view, draw, view.U.y1(), view.U.y2(), x );
		}
	}

void Graph::draw_linear_grid_y( 
	float y_jump_size, 
	Plane plane, 
	Color c 
	)
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	if( y_jump_size <= 0 ) return;

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;		
		const float y_start = std::ceil(  view.U.y1() / y_jump_size ) * y_jump_size;
		const float y_end   = std::floor( view.U.y2() / y_jump_size ) * y_jump_size;
		for( float y = y_start; y <= y_end; y += y_jump_size ) 
			draw_horizontal_line( view, draw, view.U.x1(), view.U.x2(), y );
		}
	}

void Graph::draw_linear_grid( float x_jump_size, float y_jump_size, Plane plane, Color c )
	{
	draw_linear_grid_x( x_jump_size, plane, c );
	draw_linear_grid_y( y_jump_size, plane, c );
	}

void Graph::draw_log_grid_x( 
	float x_jump_size, 
	uint32_t lines_per_step,
	Plane plane, 
	Color c 
	)
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	if( x_jump_size <= 0 ) return;

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;		
		const float x_start = std::floor( view.U.x1() / x_jump_size ) * x_jump_size;
		const float x_end   = std::ceil( view.U.x2() / x_jump_size ) * x_jump_size;
		for( float x_linear = x_start; x_linear <= x_end; x_linear += x_jump_size )
			{
			for( int step = 0; step < lines_per_step; ++step )
				{
				const float x = x_linear + std::log( 1.0f + step ) / std::log( lines_per_step );
				if( view.U.x1() <= x && x < view.U.x2() )
					draw_vertical_line( view, draw, view.U.y1(), view.U.y2(), x );
				}
			}
		}
	}

void Graph::draw_log_grid_y( 
	float y_jump_size, 
	uint32_t lines_per_step,
	Plane plane, 
	Color c 
	)
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	if( y_jump_size <= 0 ) return;

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;		
		const float y_start = std::floor( view.U.y1() / y_jump_size ) * y_jump_size;
		const float y_end   = std::ceil( view.U.y2() / y_jump_size ) * y_jump_size;
		for( float y_linear = y_start; y_linear <= y_end; y_linear += y_jump_size )
			{
			for( int step = 0; step < lines_per_step; ++step )
				{
				const float y = y_linear + std::log( 1.0f + step ) / std::log( lines_per_step );
				if( view.U.y1() <= y && y < view.U.y2() )
					draw_horizontal_line( view, draw, view.U.x1(), view.U.x2(), y );
				}
			}
		}
	}

void Graph::draw_x_ticks( 
	float jump, 
	float y, 
	float scale_base,
	Pixel offset_down, 
	Pixel offset_up, 
	Plane plane, 
	Color c, 
	float number_scale 
	)
	{
	if( jump <= 0 ) return;

	image_drawer draw( *this );
	draw.pen_color( c );

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;

		const float y_start = std::clamp( y - view.hVToU( offset_down ), view.U.y1(), view.U.y2() );
		const float y_end   = std::clamp( y + view.hVToU( offset_up   ), view.U.y1(), view.U.y2() );

		const float x_start = std::ceil(  view.U.x1() / jump ) * jump;
		const float x_end   = std::floor( view.U.x2() / jump ) * jump;
		for( float x = x_start; x <= x_end; x += jump ) 
			{
			draw_vertical_line( view, draw, y_start, y_end, x );
			if( number_scale > 0 ) 
				{
				const float number = std::pow( scale_base, x );
				draw_float( { x, y_start - view.hVToU( 12 ) }, number_scale * 4 / 5, number_scale, number, plane, c );
				}
			}
		}
	}

void Graph::draw_y_ticks( 
	float jump, 
	float x, 
	float scale_base,
	Pixel offset_left, 
	Pixel offset_right, 
	Plane plane, 
	Color c, 
	float number_scale 
	)
	{
	if( jump <= 0 ) return;

	image_drawer draw( *this );
	draw.pen_color( c );

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;

		const float x_start = std::clamp( x - view.wVToU( offset_left  ), view.U.x1(), view.U.x2() );
		const float x_end   = std::clamp( x + view.wVToU( offset_right ), view.U.x1(), view.U.x2() );

		const float y_start = std::ceil(  view.U.y1() / jump ) * jump;
		const float y_end   = std::floor( view.U.y2() / jump ) * jump;
		for( float y = y_start; y <= y_end; y += jump ) 
			{
			draw_horizontal_line( view, draw, x_start, x_end, y );
			if( number_scale > 0 ) 
				{
				const float number = std::pow( scale_base, y );
				draw_float( { x_end, y - .5f * view.hVToU( 10 ) }, number_scale * 4 / 5, number_scale, number, plane, c );
				}
			}
		}
	}

void Graph::draw_point( const vec2 & p, Pixel r, Plane plane, Color c )
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;

		const Pixel xMid = view.xUToV( p.x() );
		const Pixel yMid = view.yUToV( p.y() );
		const Pixel x_start  = std::clamp( xMid - r, (int) view.V.x1(), int( view.V.x2() ) - 1 );
		const Pixel x_end	= std::clamp( xMid + r, (int) view.V.x1(), int( view.V.x2() ) - 1 );

		for( Pixel x = x_start; x <= x_end; ++x )
			{
			const Pixel dis = std::abs( x - xMid );
			const Pixel offset = std::floor( std::sqrt( r * r - dis * dis ) );
			const Pixel y_start = std::clamp( yMid - offset, (int) view.V.y1(), (int) view.V.y2() - 1 );
			const Pixel y_end   = std::clamp( yMid + offset, (int) view.V.y1(), (int) view.V.y2() - 1 );
			draw.vertical_line_segment( height() - 1 - y_end, height() - 1 - y_start, x );
			}
		}
	}

void Graph::draw_points( const std::vector<vec2> & ps, Pixel r, Plane plane, Color c )
	{
	for( auto & p : ps )
		draw_point( p, r, plane, c );
	}

static int get_num_digits( float x )  
	{  
    return 
		(x < 1e1 ? 1 :   
        (x < 1e2 ? 2 :   
        (x < 1e3 ? 3 :   
        (x < 1e4 ? 4 :   
        (x < 1e5 ? 5 :   
        (x < 1e6 ? 6 :   
        (x < 1e7 ? 7 :  
        (x < 1e8 ? 8 :  
        (x < 1e9 ? 9 :  
        (x < 1e10 ? 10 :  
        (x < 1e11 ? 11 :  
        (x < 1e12 ? 12 :  
        (x < 1e13 ? 13 :  
        (x < 1e14 ? 14 :  
        (x < 1e15 ? 15 :  
        (x < 1e16 ? 16 :  
        (x < 1e17 ? 17 :  
        18 )))))))))))))))));  
	} 

static std::vector<int8_t> get_digits( int n, int num_digits )
	{
	std::vector<int8_t> digits( num_digits );
	size_t n_copy = n; //Need a copy to work with
	for( size_t i = 0; i < num_digits; ++i )
		{
		digits[digits.size() - 1 - i] = n_copy % 10;
		n_copy = std::floor( n_copy / 10 );
		}
	return digits;
	}

void Graph::draw_float( vec2 pos, Pixel digit_width, Pixel digit_height, float number, Plane plane, Color c )
	{
	// Digit seperation
	const bool negative = number < 0;
	if( number < 0 ) number = -number;
	const int q = std::floor( number );
	const int r = std::round( ( number - q ) * 1000.0f );
	std::vector<int8_t> wholeDigits = get_digits( q, get_num_digits( q ) );
	std::vector<int8_t> remDigits = get_digits( r, 3 );
	std::vector<int8_t> digits;
	if( negative ) digits.push_back( -1 );
	digits.insert( digits.end(), wholeDigits.begin(), wholeDigits.end() );
	digits.push_back( 10 );
	digits.insert( digits.end(), remDigits.begin(), remDigits.end() );

	// Manual drawing routines for each digit
	image_drawer draw( *this );
	draw.pen_width( 1 );
	draw.pen_color( c );

	const auto active_views = get_intersecting_views( Rect::R2, plane );
	for( auto & plane_view : active_views )
		{
		const View & view = plane_view.second;

		const float w = view.wVToU( digit_width );
		const float y2 = pos.y() + view.hVToU( digit_height );

		auto draw_path = [this, &draw, &view]( const Rect & r, const std::vector<vec2> & ps )
			{
			for( int i = 1; i < ps.size(); ++i )
				draw_line_segment( view, draw, 
					r.x1() + r.w() * ps[i-1].x(), 
					r.y1() + r.h() * ps[i-1].y(), 
					r.x1() + r.w() * ps[i].x(), 
					r.y1() + r.h() * ps[i].y() );
			};

		float xPos = pos.x();
		for( auto d : digits )
			{
			const Rect fullDigitRect = { xPos, pos.y(), xPos + w, y2 };
			const Rect digitRect = fullDigitRect.intersect( view.U );
			if( fullDigitRect != digitRect ) return;

			constexpr float x1 = .15, x2 = .85;
			switch( d )
				{
				case -1: // Used for minus
					draw_path( digitRect, { {x1, .5}, {x2, .5} } );
					break;
				case 0: 
					draw_path( digitRect, { {x1, 0}, {x2, 0}, {x2, 1}, {x1, 1}, {x1, 0}, {x1, 1} } );
					break;
				case 1: 
					draw_path( digitRect, { {.5, 0}, {.5, 1} } );
					break;
				case 2: 
					draw_path( digitRect, { {x1, 1}, {x2, 1}, {x2, .5}, {x1, .5}, {x1, 0}, {x2, 0} } );
					break;
				case 3: 
					draw_path( digitRect, { {x1, 1}, {x2, 1}, {x2, .5}, {x1, .5}, {x2, .5}, {x2, 0}, {x1, 0} } );
					break;
				case 4: 
					draw_path( digitRect, { {x1, 1}, {x1, .5}, {x2, .5}, {x2, 1}, {x2, 0} } );
					break;
				case 5: 
					draw_path( digitRect, { {x2, 1}, {x1, 1}, {x1, .5}, {x2, .5}, {x2, 0}, {x1, 0} } );
					break;
				case 6: 
					draw_path( digitRect, { {x2, 1}, {x1, 1}, {x1, 0}, {x2, 0}, {x2, .5}, {x1, .5} } );
					break;
				case 7: 
					draw_path( digitRect, { {x1, 1}, {x2, 1}, {.5, 0} } );
					break;
				case 8: 
					draw_path( digitRect, { {x2, .5}, {x2, 1}, {x1, 1}, {x1, 0}, {x2, 0}, {x2, .5}, {x1, .5} } );
					break;
				case 9: 
					draw_path( digitRect, { {x2, .5}, {x1, .5}, {x1, 1}, {x2, 1}, {x2, 0} } );
					break;
				case 10: // Dot
					draw_path( digitRect, { {.4, 0}, {.6, 0}, {.6, .2}, {.4, .2}, {.4, 0} } );
					break;
				default:
					draw_path( digitRect, { {x1, 0}, {x2, 0}, {x1, 1}, {x2, 1}, {x1, 0} } );
					break;
				}

			xPos += w;
			};
		}
	}