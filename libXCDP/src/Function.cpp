#include "xcdp/Function.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

static const float pi = acos( -1.0f );

namespace xcdp {

void Func1x1::graph( const std::string & filename,
		float left, float right, float bottom, float top, size_t resolution ) const
	{
	std::cout << "Generating function graph ... ";

	if( left > right ) std::swap( left, right );
	if( bottom > top ) std::swap( bottom, top );
	const float res_f = float( resolution );
	const size_t width  = (right - left) * float( resolution );
	const size_t height = (top - bottom) * float( resolution );

	std::vector<float> outputs( width, 0 );
	for( size_t x = 0; x < outputs.size(); ++x )
		outputs[x] = std::round( ( f( float(x) / res_f + left ) - bottom ) * res_f );

	//y -> x to match write order
	std::vector<std::vector<std::array<uint8_t,3>>> data( width, std::vector( height, std::array<uint8_t,3>( {0,0,0} ) ) );
	for( uint32_t y = 0; y < height; ++y )
		{
		for( uint32_t x = 0; x < width; ++x )	
			{
			float value = ( y == outputs[x] ? 0.0f : 1.0f);
			if( x + left   * resolution == 0 ) value = 0.0f;
			if( y + bottom * resolution == 0 ) value = 0.0f;
			data[x][y] = HSVtoRGB( 0, 0.0, value );
			}
		}

	writeBMP( filename, data );

	std::cout << "Done\n";
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

Func1x1 Func1x1::interpolatePoints( const std::vector< std::pair< float, float > > ps, Interpolator interp )
	{
	return [ps, interp]( float t )
		{
		if( ps.size() == 0 ) return 0.0f;
		if( t < ps[0].first || ps[ps.size()-1].first < t ) return 0.0f;
		for( size_t i = 1; i < ps.size(); ++i )
			if( t < ps[i].first )
				{
				const float mix = interp( ( t - ps[i-1].first ) / ( ps[i].first - ps[i-1].first ) );
				return ( 1.0f - mix ) * ps[i-1].second + mix * ps[i].second;
				}
		return 0.0f;
		};
	}


}