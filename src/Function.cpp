#include "flan/Function.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "bmp/bitmap_image.hpp"
#include "flan/Graph.h"

static const float pi = std::acos( -1.0f );
static const float pi2 = 2.0f * pi;

namespace flan {

Graph Func1x1::convertToGraph( Rect view, Interval domain, Pixel width, Pixel height, flan_CANCEL_ARG_CPP ) const
	{
	Graph g( width, height );
	g.setView( view );
	g.fillImage( Color::White );
	g.drawAxes( 0, Color::Black );
	g.drawLinearGrid( 1, 1, 0, Color( 200, 200, 200 ) );
	g.drawFunction( *this, domain, -1, Color::Black, canceller );
	return g;
	}

Func1x1 Func1x1::saveAsBMP( const std::string & filename, Rect view, Interval domain, Pixel width, Pixel height, flan_CANCEL_ARG_CPP ) const
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


Func1x1 Func1x1::sine 		= []( float t ){ return std::sin( pi2 * t ); };
Func1x1 Func1x1::square		= []( float t ){ return t < 0.5f ? -1.0f : 1.0f; };
Func1x1 Func1x1::saw 		= []( float t ){ return -1.0f + 2.0f * t; };
Func1x1 Func1x1::triangle 	= []( float t ){ return t < 0.5f ? -1.0f + 4.0f * t : 3.0f - 4.0f * t; };

}