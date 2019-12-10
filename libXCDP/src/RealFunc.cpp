#include "RealFunc.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

namespace xcdp {

void RealFunc::graph( const std::string & filename,
		double left, double right, double bottom, double top, size_t resolution ) const
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
	std::vector<std::vector<std::array<char,3>>> data( width, std::vector( height, std::array<char,3>( {0,0,0} ) ) );
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

RealFunc RealFunc::ADSR( double a, double d, double s, double r, double sLvl, double aExp, double dExp, double rExp )
	{
	return [a, d, s, r, sLvl, aExp, dExp, rExp]( double t )
		{
			 if( t < 0 || t > a+d+s+r ) return 0.0;
		else if( t < a			   ) return std::pow( t / a, aExp );
		else if( t < a + d		   ) return std::pow( 1.0 - ( t - a ) / d, dExp ) * ( 1.0 - sLvl ) + sLvl;
		else if( t < a + d + s	   ) return sLvl;
		else if( t < a + d + s + r ) return std::pow( 1.0 - ( t - a - d - s ) / r, rExp ) * sLvl;
		};
	}

RealFunc RealFunc::interpolatePoints( const std::vector< std::pair< double, double > > ps, Interpolator interp )
	{
	return [ps, interp]( double t )
		{
		if( ps.size() == 0 ) return 0.0;
		if( t < ps[0].first || ps[ps.size()-1].first < t ) return 0.0;
		for( size_t i = 1; i < ps.size(); ++i )
			if( t < ps[i].first )
				return interp( ( t - ps[i-1].first ) / ( ps[i].first - ps[i-1].first ), ps[i-1].second, ps[i].second );
		};
	}


}