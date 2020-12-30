#include "xcdp/Function.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "bmp/bitmap_image.hpp"
#include "xcdp/spline.h"
#include "xcdp/Graph.h"

static const float pi = acos( -1.0f );

namespace xcdp {

Graph Func1x1::convertToGraph( Rect view, Interval domain, Pixel width, Pixel height, XCDP_CANCEL_ARG_CPP ) const
	{
	Graph g( width, height );
	g.setView( view );
	g.fillImage( Color::White );
	g.drawAxes( 0, Color::Black );
	g.drawLinearGrid( 1, 1, 0, Color( 200, 200, 200 ) );
	g.drawFunction( *this, domain, -1, Color::Black, canceller );
	return g;
	}

Func1x1 Func1x1::saveAsBMP( const std::string & filename, Rect view, Interval domain, Pixel width, Pixel height, XCDP_CANCEL_ARG_CPP ) const
	{
	auto bmp = convertToGraph( view, domain, width, height, canceller );
	bmp.save_image( filename );
	return *this;
	}

Func1x1 Func1x1::ADSR( float a, float d, float s, float r, float sLvl, float aExp, float dExp, float rExp )
	{
	return [a, d, s, r, sLvl, aExp, dExp, rExp]( float t )
		{
			 if( t < 0 || t > a+d+s+r ) return (float) 0.0;
		else if( t < a			   ) return std::pow( t / a, aExp );
		else if( t < a + d		   ) return std::pow( (float) 1.0 - ( t - a ) / d, dExp ) * ( (float) 1.0 - sLvl ) + sLvl;
		else if( t < a + d + s	   ) return sLvl;
		else if( t < a + d + s + r ) return std::pow( (float) 1.0 - ( t - a - d - s ) / r, rExp ) * sLvl;
		else return 0.0f;
		};
	}

Func1x1 Func1x1::periodize( Func1x1 period )
	{
	return [=]( float t )
		{
		const float p = period( t );
		return f( std::fmod( t, p ) + ( t < 0.0f ? p : 0.0f ) );
		};
	}

Func1x1 Func1x1::interpolatePoints( const std::vector< std::pair< float, float > > & ps, Interpolator interp )
	{
	return [ps, interp]( float t )
		{
		if( ps.size() == 0 ) return 0.0f;
		if( t <= ps.front().first ) return ps.front().second;
		if( ps.back().first <= t ) return ps.back().second;

		const auto p2 = std::lower_bound( ps.begin(), ps.end(), t, []( const std::pair< float, float > & p, float t ){ return p.first < t; } );
		const auto p1 = p2 - 1;
		const float mix = interp( ( t - p1->first ) / ( p2->first - p1->first ) );
		return ( 1.0f - mix ) * p1->second + mix * p2->second;
		};
	}

Func1x1 Func1x1::interpolatePoints( float deltaX, const std::vector< float > ys, Interpolator interp )
	{
	std::vector<std::pair< float, float >> points( ys.size() );

	float xAccum = 0;
	std::transform( ys.begin(), ys.end(), points.begin(), [deltaX, &xAccum]( float y )
		{ 
		auto point = std::pair<float,float>( xAccum, y );
		xAccum += deltaX;
		return point;
		} );

	return interpolatePoints( points, interp );
	}

Func1x1 Func1x1::spline( const std::vector< std::pair< float, float > > points )
	{
	// Split points into xs and ys
	std::vector<double> xs( points.size() );
	std::vector<double> ys( points.size() );
	std::transform( points.begin(), points.end(), xs.begin(), []( std::pair<float,float> p ){ return p.first; } );
	std::transform( points.begin(), points.end(), ys.begin(), []( std::pair<float,float> p ){ return p.second; } );

	// Generate tk::spline from xs and ys
	tk::spline spline;
	spline.set_points( xs, ys );

	return [spline]( float t )
		{
		return spline( t );
		};
	}

//Func1x1 Func1x1::uniformDistribution( Func1x1 lowerBound, Func1x1 upperBound )
//	{
//	static std::default_random_engine rng( std::time( nullptr ) );
//	return [lowerBound, upperBound]( float t )
//		{ 
//		return std::uniform_real_distribution<float>( lowerBound( t ), upperBound( t ) )( rng ); 
//		};
//	}
//
//Func1x1 Func1x1::normalDistribution( Func1x1 mean, Func1x1 sigma )
//	{
//	static std::default_random_engine rng( std::time( nullptr ) );
//	return [mean, sigma]( float t )
//		{ 
//		return std::normal_distribution<float>( mean( t ), sigma( t ) )( rng ); 
//		};
//	}

}