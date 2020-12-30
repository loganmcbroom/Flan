#include "..\include\xcdp\Graph.h"

#include <algorithm>

#include "xcdp/Function.h"

using namespace xcdp;

Graph::Graph( Pixel _width, Pixel _height, Rect _view )
	: bitmap_image( _width == -1? Graph::DefaultWidth : _width, _height == -1 ? Graph::DefaultHeight : _height )
	, views()
	{
	}

void Graph::addView( View view, int plane )
	{
	views.emplace_back( plane, view );
	}

void Graph::setView( Rect rect )
	{
	views.clear();
	views.emplace_back( -1, View( rect, Rect( 0.0f, 0.0f, (float) width(), (float) height() ) ) );
	}

void Graph::addSplitViewY( View view, int numViews, int startPlane )
	{
	float y_c = view.V.y1();
	const float viewHeight = view.V.h() / float( numViews );
	for( int v = 0; v < numViews; ++v )
		{
		const View newView( view.U, { view.V.x1(), y_c, view.V.x2(), y_c + viewHeight } );
		addView( newView, startPlane + v );
		y_c += viewHeight;
		}
	}

void Graph::addFullSplitViewY( Rect rect, int numViews, int startPlane )
	{
	addSplitViewY( View( rect, Rect( 0, 0, width(), height() ) ), numViews, startPlane );
	}

bool Graph::doPlanesMatch( int p1, int p2 ) const
	{
	if( p1 == -1 || p2 == -1 ) return true;
	else return p1 == p2;
	}

std::vector<std::pair<int,View>> Graph::getIntersectingViews( Rect U, int plane ) const
	{
	std::vector<std::pair<int,View>> activeViews;
	for( auto & v : views )
		{
		if( doPlanesMatch( v.first, plane ) && v.second.U.intersect( U ).valid() )
			activeViews.emplace_back( v );
		}
	return activeViews;
	}

//======================================================================================================================================================
// Waveforms
//======================================================================================================================================================

void Graph::drawWaveform( Func1x1 data, Rect rect, int plane, WaveformMode mode, XCDP_CANCEL_ARG_CPP )
	{
	XCDP_FUNCTION_LOG;

	const auto activeViews = getIntersectingViews( rect, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;
		const Rect drawRect = rect.intersect( view.U );

		const Pixel startPixel = std::ceil(  view.xUToV( drawRect.x1() ) );
		const Pixel endPixel   = std::floor( view.xUToV( drawRect.x2() ) );
		const Pixel yMidPixel = view.yUToV( rect.r.midpoint() );

		// Oversample data into buffer
		for( Pixel x = startPixel; x < endPixel; ++x )
			{
			XCDP_CANCEL_POINT();

			//Get average of subsamples
			float sum = 0;
			for( uint32_t sample = 0; sample < oversample; ++sample )
				{
				const float subsample = view.xVToU( float( x ) + float( sample ) / oversample );
				const float data_c = data( subsample );
				sum += mode == WaveformMode::Direct? data_c : std::abs( data_c );
				}
			sum /= oversample;

			const Pixel yOffsetPixels = view.hUToV( std::clamp( sum, -1.0f, 1.0f ) * rect.h() / 2.0f );

			auto setPixelWithColor = [this, yMidPixel, &rect, &view]( Pixel x, Pixel y )
				{ 
				if( y < view.V.y1() ||  view.V.y2() <= y ) return;
				const float r = std::abs( y - yMidPixel ) / view.hUToV( rect.h() / 2.0f );
				const Color c = Color::fromHSV( 90.0f * r + hue, .8f, .65 ); 
				setPixel( x, y, c );
				};

			if( mode == WaveformMode::Direct )
				{
				if( yOffsetPixels < 0 )
					for( Pixel y = yMidPixel; y >= yMidPixel + yOffsetPixels; --y ) setPixelWithColor( x, y );
				else
					for( Pixel y = yMidPixel; y <= yMidPixel + yOffsetPixels; ++y ) setPixelWithColor( x, y );
				}
			else
				for( Pixel y = yMidPixel - yOffsetPixels; y <= yMidPixel + yOffsetPixels; ++y ) 
					setPixelWithColor( x, y );
			}
		}
	}

void Graph::drawWaveform( const float * data, int n, Rect area, int plane, WaveformMode mode, XCDP_CANCEL_ARG_CPP )
	{
	drawWaveforms( { data }, n, area, plane, mode, canceller );
	}

void Graph::drawWaveforms( std::vector<Func1x1> fs, Rect rect, int startPlane, WaveformMode mode, XCDP_CANCEL_ARG_CPP )
	{
	for( int f = 0; f < fs.size(); ++f )
		{
		hue = 360.0f * f / fs.size();
		drawWaveform( fs[f], rect, startPlane + f, mode, canceller );
		}
	}

void Graph::drawWaveforms( std::vector<const float *> datas, int n, Rect rect, int startPlane, WaveformMode mode, XCDP_CANCEL_ARG_CPP )
	{
	std::vector<Func1x1> fs;
	for( auto & d : datas ) 
		fs.push_back( [&d, n, &rect]( float x )
			{ 
			const int i = std::floor( ( x - rect.x1() ) / rect.w() * n );
			//if( i < 0 || i >= n ) return 0.0f;
			return d[i]; 
			} );
	drawWaveforms( fs, rect, startPlane, mode, canceller );
	}


//======================================================================================================================================================
// Spectrograms
//======================================================================================================================================================

void Graph::drawSpectrogram( Func2x1 data, Rect rect, int plane, XCDP_CANCEL_ARG_CPP )
	{
	XCDP_FUNCTION_LOG;

	const int oversample_c = std::sqrt( oversample );

	const auto activeViews = getIntersectingViews( rect, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;
		const Rect drawRect = rect.intersect( view.U );

		const Pixel startXPixel = std::ceil(  view.xUToV( drawRect.x1() ) );
		const Pixel endXPixel	= std::floor( view.xUToV( drawRect.x2() ) );
		const Pixel startYPixel = std::ceil(  view.yUToV( drawRect.y1() ) );
		const Pixel endYPixel	= std::floor( view.yUToV( drawRect.y2() ) );

		for( Pixel y = startYPixel; y < endYPixel; ++y )
			{
			for( Pixel x = startXPixel; x < endXPixel; ++x )
				{
				XCDP_CANCEL_POINT();
				
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

				const Color color = Color::fromHSV( hue, 1.0f, std::clamp( mag, 0.0f, 1.0f ) );
				setPixel( x, y, color );
				}
			}
		}
	}

void Graph::drawSpectrogram( const float * data, int n, int m, Rect rect, int plane, XCDP_CANCEL_ARG_CPP )
	{
	drawSpectrograms( { data }, n, m, rect, plane, canceller );
	}

void Graph::drawSpectrograms( std::vector<Func2x1> fs, Rect rect, int startPlane, XCDP_CANCEL_ARG_CPP )
	{
	for( int f = 0; f < fs.size(); ++f )
		{
		hue = 360.0f * f / fs.size();
		drawSpectrogram( fs[f], rect, startPlane + f, canceller );
		}
	}

void Graph::drawSpectrograms( std::vector<const float *> datas, int n, int m, Rect rect, int startPlane, XCDP_CANCEL_ARG_CPP )
	{
	std::vector<Func2x1> fs;
	for( auto & d : datas ) 
		fs.push_back( [&d, n, m, &rect]( float x, float y )
			{ 
			const int i = std::floor( ( x - rect.x1() ) / rect.w() * n );
			const int j = std::floor( ( y - rect.y1() ) / rect.h() * m );
			//if( i < 0 || i >= n || j < 0 || j >= m ) return 0.0f;
			return d[ i * m + j ]; 
			} );
	drawSpectrograms( fs, rect, canceller );
	}


//======================================================================================================================================================
// Functions
//======================================================================================================================================================

void Graph::drawFunction( const Func1x1 & f, Interval domain, int plane, Color c, XCDP_CANCEL_ARG_CPP )
	{
	XCDP_FUNCTION_LOG;

	//Draw function
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( domain * Interval::R, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;
		const Interval drawDomain = domain.intersect( view.U.d );

		const float pixelAdvance = view.wVToU( 1 );

		float x1 = drawDomain.x1;
		float y1 = f( x1 );
		for( float x2 = drawDomain.x1 + pixelAdvance; x2 < drawDomain.x2; x2 += pixelAdvance )	
			{
			XCDP_CANCEL_POINT();
			const float y2 = f( x2 );
			if( view.U.contains( x1, y1 ) && view.U.contains( x2, y2 ) ) 
				drawLineSegment( view, draw, x1, y1, x2, y2 );
			x1 = x2;
			y1 = y2;
			}
		}
	}

void Graph::drawFunction( const std::vector<std::pair<float,float>> & data, int plane, Color c, XCDP_CANCEL_ARG_CPP )
	{
	auto f = Func1x1::interpolatePoints( data );

	auto xComp = []( const std::pair<float,float> & p,  const std::pair<float,float> & s ){ return p.first < s.first; };
	const float left	= std::min_element( data.begin(), data.end(), xComp )->first;
	const float right	= std::max_element( data.begin(), data.end(), xComp )->first;
	
	drawFunction( f, Interval( left, right ), plane, c, canceller );
	}

void Graph::drawFunctions( const std::vector<Func1x1> & fs, std::vector<Interval> domains, int plane, XCDP_CANCEL_ARG_CPP )
	{
	if( domains.size() < fs.size() )
		domains.resize( fs.size(), Interval::R );

	// Draw waveforms
	for( int f = 0; f < fs.size(); ++f )
		{
		hue = 360.0f * f / fs.size();
		drawFunction( fs[f], domains[f], plane, Color::fromHSV( hue, 1, 1 ), canceller );
		}
	}


//======================================================================================================================================================
// Primitives
//======================================================================================================================================================

void Graph::setPixel( Pixel x, Pixel y, Color c )
	{
	set_pixel( x, height() - 1 - y, c );
	}

void Graph::setPoint( const View & v, float x, float y, Color c )
	{
	setPixel( v.xUToV( x ), v.yUToV( y ), c );
	}

void Graph::drawHorizontalLine( const View & v, image_drawer & draw, float x1, float x2, float y )
	{
	draw.horiztonal_line_segment( v.xUToV( x1 ), v.xUToV( x2 ), height() - 1 - v.yUToV( y ) );
	}

void Graph::drawVerticalLine( const View & v, image_drawer & draw, float y1, float y2, float x )
	{
	draw.vertical_line_segment( height() - 1 - v.yUToV( y2 ), height() - 1 - v.yUToV( y1 ), v.xUToV( x ) );
	}

void Graph::drawLineSegment( const View & v, image_drawer & draw, float x1, float y1, float x2, float y2 )	
	{
	draw.line_segment( v.xUToV( x1 ), height() - 1 - v.yUToV( y1 ), v.xUToV( x2 ), height() - 1 - v.yUToV( y2 ) );
	}

void Graph::setRect( const View & v, Rect rect, Color c )
	{
	const Rect drawRect = rect.intersect( v.U );
	set_region( 
		v.xUToV( drawRect.x1() ), 
		height() - v.yUToV( drawRect.y2() ), // +1 handles an edge case created by bmp plane flip
		v.wUToV( drawRect.w() ), 
		v.wUToV( drawRect.h() ), 
		c.red, c.green, c.blue );
	}

void Graph::fillImage( Color c )
	{
	set_region( 0, 0, width(), height(), c.red, c.green, c.blue );
	}


//======================================================================================================================================================
// Additional Drawing Routines
//======================================================================================================================================================

void Graph::drawAxes( int plane, Color c )
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( Rect::R2, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;
		if( view.U.y1() <= 0 && view.U.y2() > 0 ) drawVerticalLine(   view, draw, view.U.y1(), view.U.y2(), 0 );
		if( view.U.x1() <= 0 && view.U.x2() > 0 ) drawHorizontalLine( view, draw, view.U.x1(), view.U.x2(), 0 );
		}
	}

void Graph::drawLinearGrid( float xJumpSize, float yJumpSize, int plane, Color c )
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( Rect::R2, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;

		if( xJumpSize > 0 )
			{
			const float xStart = std::ceil(  view.U.x1() / xJumpSize ) * xJumpSize;
			const float xEnd   = std::floor( view.U.x2() / xJumpSize ) * xJumpSize;
			for( float x = xStart; x <= xEnd; x += xJumpSize ) 
				drawVerticalLine( view, draw, view.U.y1(), view.U.y2(), x );
			}
		
		if( xJumpSize > 0 )
			{
			const float yStart = std::ceil(  view.U.y1() / yJumpSize ) * yJumpSize;
			const float yEnd   = std::floor( view.U.y2() / yJumpSize ) * yJumpSize;
			for( float y = yStart; y <= yEnd; y += yJumpSize ) 
				drawHorizontalLine( view, draw, view.U.x1(), view.U.x2(), y );
			}
		}
	}

void Graph::drawXTicks( float jump, float y, Pixel offsetDown, Pixel offsetUp, int plane, Color c, bool showNumbers )
	{
	if( jump <= 0 ) return;

	image_drawer draw( *this );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( Rect::R2, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;

		const float yStart = std::clamp( y - view.hVToU( offsetDown ), view.U.y1(), view.U.y2() );
		const float yEnd   = std::clamp( y + view.hVToU( offsetUp   ), view.U.y1(), view.U.y2() );

		const float xStart = std::ceil(  view.U.x1() / jump ) * jump;
		const float xEnd   = std::floor( view.U.x2() / jump ) * jump;
		for( float x = xStart; x <= xEnd; x += jump ) 
			{
			drawVerticalLine( view, draw, yStart, yEnd, x );
			if( showNumbers ) drawFloat( { x, yStart - view.hVToU( 12 ) }, 8, 10, x, plane, c );
			}
		}
	}

void Graph::drawYTicks( float jump, float x, Pixel offsetLeft, Pixel offsetRight, int plane, Color c, bool showNumbers )
	{
	if( jump <= 0 ) return;

	image_drawer draw( *this );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( Rect::R2, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;

		const float xStart = std::clamp( x - view.wVToU( offsetLeft  ), view.U.x1(), view.U.x2() );
		const float xEnd   = std::clamp( x + view.wVToU( offsetRight ), view.U.x1(), view.U.x2() );

		const float yStart = std::ceil(  view.U.y1() / jump ) * jump;
		const float yEnd   = std::floor( view.U.y2() / jump ) * jump;
		for( float y = yStart; y <= yEnd; y += jump ) 
			{
			drawHorizontalLine( view, draw, xStart, xEnd, y );
			if( showNumbers ) drawFloat( { xEnd, y - .5f * view.hVToU( 10 ) }, 8, 10, y, plane, c );
			}
		}
	}

void Graph::drawPoint( const vec2 & p, Pixel r, int plane, Color c )
	{
	image_drawer draw( *this );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( Rect::R2, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;

		const Pixel xMid = view.xUToV( p.x() );
		const Pixel yMid = view.yUToV( p.y() );
		const Pixel xStart  = std::clamp( xMid - r, (int) view.V.x1(), int( view.V.x2() ) - 1 );
		const Pixel xEnd	= std::clamp( xMid + r, (int) view.V.x1(), int( view.V.x2() ) - 1 );

		for( Pixel x = xStart; x <= xEnd; ++x )
			{
			const Pixel dis = std::abs( x - xMid );
			const Pixel offset = std::floor( std::sqrt( r * r - dis * dis ) );
			const Pixel yStart = std::clamp( yMid - offset, (int) view.V.y1(), (int) view.V.y2() - 1 );
			const Pixel yEnd   = std::clamp( yMid + offset, (int) view.V.y1(), (int) view.V.y2() - 1 );
			draw.vertical_line_segment( height() - 1 - yEnd, height() - 1 - yStart, x );
			}
		}
	}

void Graph::drawPoints( const std::vector<vec2> & ps, Pixel r, int plane, Color c )
	{
	for( auto & p : ps )
		drawPoint( p, r, plane, c );
	}

static int getNumDigits( float x )  
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

static std::vector<int8_t> getDigits( int n, int numDigits )
	{
	std::vector<int8_t> digits( numDigits );
	size_t n_copy = n; //Need a copy to work with
	for( size_t i = 0; i < numDigits; ++i )
		{
		digits[digits.size() - 1 - i] = n_copy % 10;
		n_copy = std::floor( n_copy / 10 );
		}
	return digits;
	}

void Graph::drawFloat( vec2 pos, Pixel digitWidth, Pixel digitHeight, float number, int plane, Color c, XCDP_CANCEL_ARG_CPP )
	{
	// Digit seperation
	const bool negative = number < 0;
	if( number < 0 ) number = -number;
	const int q = std::floor( number );
	const int r = std::round( ( number - q ) * 1000.0f );
	std::vector<int8_t> wholeDigits = getDigits( q, getNumDigits( q ) );
	std::vector<int8_t> remDigits = getDigits( r, 3 );
	std::vector<int8_t> digits;
	if( negative ) digits.push_back( -1 );
	digits.insert( digits.end(), wholeDigits.begin(), wholeDigits.end() );
	digits.push_back( 10 );
	digits.insert( digits.end(), remDigits.begin(), remDigits.end() );

	// Manual drawing routines for each digit
	image_drawer draw( *this );
	draw.pen_width( 1 );
	draw.pen_color( c );

	const auto activeViews = getIntersectingViews( Rect::R2, plane );
	for( auto & planeView : activeViews )
		{
		const View & view = planeView.second;

		const float w = view.wVToU( digitWidth );
		const float y2 = pos.y() + view.hVToU( digitHeight );

		auto drawPath = [this, &draw, &view]( const Rect & r, const std::vector<vec2> & ps )
			{
			for( int i = 1; i < ps.size(); ++i )
				drawLineSegment( view, draw, 
					r.x1() + r.w() * ps[i-1].x(), 
					r.y1() + r.h() * ps[i-1].y(), 
					r.x1() + r.w() * ps[i].x(), 
					r.y1() + r.h() * ps[i].y() );
			};

		float xPos = pos.x();
		for( auto d : digits )
			{
			XCDP_CANCEL_POINT();

			const Rect fullDigitRect = { xPos, pos.y(), xPos + w, y2 };
			const Rect digitRect = fullDigitRect.intersect( view.U );
			if( fullDigitRect != digitRect ) return;

			constexpr float x1 = .15, x2 = .85;
			switch( d )
				{
				case -1: // Used for minus
					drawPath( digitRect, { {x1, .5}, {x2, .5} } );
					break;
				case 0: 
					drawPath( digitRect, { {x1, 0}, {x2, 0}, {x2, 1}, {x1, 1}, {x1, 0}, {x1, 1} } );
					break;
				case 1: 
					drawPath( digitRect, { {.5, 0}, {.5, 1} } );
					break;
				case 2: 
					drawPath( digitRect, { {x1, 1}, {x2, 1}, {x2, .5}, {x1, .5}, {x1, 0}, {x2, 0} } );
					break;
				case 3: 
					drawPath( digitRect, { {x1, 1}, {x2, 1}, {x2, .5}, {x1, .5}, {x2, .5}, {x2, 0}, {x1, 0} } );
					break;
				case 4: 
					drawPath( digitRect, { {x1, 1}, {x1, .5}, {x2, .5}, {x2, 1}, {x2, 0} } );
					break;
				case 5: 
					drawPath( digitRect, { {x2, 1}, {x1, 1}, {x1, .5}, {x2, .5}, {x2, 0}, {x1, 0} } );
					break;
				case 6: 
					drawPath( digitRect, { {x2, 1}, {x1, 1}, {x1, 0}, {x2, 0}, {x2, .5}, {x1, .5} } );
					break;
				case 7: 
					drawPath( digitRect, { {x1, 1}, {x2, 1}, {.5, 0} } );
					break;
				case 8: 
					drawPath( digitRect, { {x2, .5}, {x2, 1}, {x1, 1}, {x1, 0}, {x2, 0}, {x2, .5}, {x1, .5} } );
					break;
				case 9: 
					drawPath( digitRect, { {x2, .5}, {x1, .5}, {x1, 1}, {x2, 1}, {x2, 0} } );
					break;
				case 10: // Dot
					drawPath( digitRect, { {.4, 0}, {.6, 0}, {.6, .2}, {.4, .2}, {.4, 0} } );
					break;
				default:
					drawPath( digitRect, { {x1, 0}, {x2, 0}, {x1, 1}, {x2, 1}, {x1, 0} } );
					break;
				}

			xPos += w;
			};
		}
	}