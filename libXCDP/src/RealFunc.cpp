#include "RealFunc.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>

using namespace xcdp;

std::array<char,3> HSVtoRGB( int H, double S, double V ) 
	{
	double C = S * V;
	double X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	double m = V - C;
	double Rs, Gs, Bs;

		 if( H >= 0   && H < 60  ) { Rs = C; Gs = X; Bs = 0; }
	else if( H >= 60  && H < 120 ) { Rs = X; Gs = C; Bs = 0; }
	else if( H >= 120 && H < 180 ) { Rs = 0; Gs = C; Bs = X; }
	else if( H >= 180 && H < 240 ) { Rs = 0; Gs = X; Bs = C; }
	else if( H >= 240 && H < 300 ) { Rs = X; Gs = 0; Bs = C; }
	else						   { Rs = C; Gs = 0; Bs = X; }

	return std::array<char,3>({ char((Rs + m) * 255),  char((Gs + m) * 255),  char((Bs + m) * 255) });
	}

void RealFunc::graph( const std::string & filename,
		double left, double right, double bottom, double top, size_t resolution ) const
	{
	std::cout << "Generating function graph ... ";
	std::ofstream file( filename, std::ios::binary );
	if( !file )
		{
		std::cout << "Error opening " << filename << " to save RealFunc graph.\n";
		return;
		}

	if( left > right ) std::swap( left, right );
	if( bottom > top ) std::swap( bottom, top );

	const double resD = double( resolution );
	const size_t width  = (right - left) * resolution;
	const size_t height = (top - bottom) * resolution;

	//The image header
	unsigned char header[ 18 ] = { 0 };
	header[  2 ] = 1;  //Uncompressed, RGB image
	header[ 12 ] = ( width		     )	& 0xFF;
	header[ 13 ] = ( width		>> 8 )	& 0xFF;
	header[ 14 ] = ( height		     )	& 0xFF;
	header[ 15 ] = ( height		>> 8 )	& 0xFF;
	header[ 16 ] = 24;  //Bits per pixel

	file.write( (const char *) header, 18 );

	std::vector<int> outputs( width, 0 );
	for( size_t x = 0; x < outputs.size(); ++x )
		outputs[x] = int( ( f( double(x) / resD + left ) - bottom ) * resD );

	//y -> x to match tga write order
	for( int y = 0; y < height; ++y )
		{
		for( int x = 0; x < width; ++x )	
			{
			double value = ( y == outputs[x] ? 0.0 : 1.0);
			if( x + left   * resolution == 0 ) value = 0.0;
			if( y + bottom * resolution == 0 ) value = 0.0;
			const std::array<char,3> rgb = HSVtoRGB( 0, 0.0, value );
			file.write( rgb.data(), 3 );
			}
		}

	file.close();

	std::cout << "Done\n";
	}

RealFunc xcdp::ADSR( double a, double d, double s, double r, double sLvl, double aExp, double dExp, double rExp )
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
