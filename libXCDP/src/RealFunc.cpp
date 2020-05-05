#include "RealFunc.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

static const double pi = acos( -1 );

namespace xcdp {

void RealFunc::graph( const std::string & filename,
		float left, float right, float bottom, float top, size_t resolution ) const
	{
	std::cout << "Generating function graph ... ";

	if( left > right ) std::swap( left, right );
	if( bottom > top ) std::swap( bottom, top );
	const double resD = double( resolution );
	const size_t width  = (right - left) * resolution;
	const size_t height = (top - bottom) * resolution;

	
	std::vector<double> outputs( width, 0 );
	for( size_t x = 0; x < outputs.size(); ++x )
		outputs[x] = int( ( f( double(x) / resD + left ) - bottom ) * resD );

	//y -> x to match write order
	std::vector<std::vector<std::array<uint8_t,3>>> data( width, std::vector( height, std::array<uint8_t,3>( {0,0,0} ) ) );
	for( int y = 0; y < height; ++y )
		{
		for( int x = 0; x < width; ++x )	
			{
			double value = ( y == outputs[x] ? 0.0 : 1.0);
			if( x + left   * resolution == 0 ) value = 0.0;
			if( y + bottom * resolution == 0 ) value = 0.0;
			data[x][y] = HSVtoRGB( 0, 0.0, value );
			}
		}

	writeBMP( filename, data );

	std::cout << "Done\n";
	}

RealFunc RealFunc::ADSR( float a, float d, float s, float r, float sLvl, float aExp, float dExp, float rExp )
	{
	return [a, d, s, r, sLvl, aExp, dExp, rExp]( float t )
		{
			 if( t < 0 || t > a+d+s+r ) return (float) 0.0;
		else if( t < a			   ) return std::pow( t / a, aExp );
		else if( t < a + d		   ) return std::pow( (float) 1.0 - ( t - a ) / d, dExp ) * ( (float) 1.0 - sLvl ) + sLvl;
		else if( t < a + d + s	   ) return sLvl;
		else if( t < a + d + s + r ) return std::pow( (float) 1.0 - ( t - a - d - s ) / r, rExp ) * sLvl;
		};
	}

RealFunc RealFunc::oscillate( RealFunc min, RealFunc max, RealFunc period, RealFunc wave )
	{
	return [&]( float t )
		{
		const double center = ( max(t) + min(t) ) / 2.0;
		const double amp    = ( max(t) - min(t) ) / 2.0;
		return center + amp * wave( 2.0 * pi * t / period(t) );
		};
	}

RealFunc RealFunc::interpolatePoints( const std::vector< std::pair< float, float > > ps, Interpolator interp )
	{
	return [ps, interp]( float t )
		{
		if( ps.size() == 0 ) return 0.0f;
		if( t < ps[0].first || ps[ps.size()-1].first < t ) return 0.0f;
		for( size_t i = 1; i < ps.size(); ++i )
			if( t < ps[i].first )
				return interp( ( t - ps[i-1].first ) / ( ps[i].first - ps[i-1].first ), ps[i-1].second, ps[i].second );
		};
	}


}