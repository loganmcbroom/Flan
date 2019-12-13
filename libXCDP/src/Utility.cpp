#include "Utility.h"

#include <random>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstring>

#include <bitset> //temp include

static double pi = acos(-1);

namespace xcdp {

const Interpolator Interpolators::constant = []( double mix, double a, double b ) 
	{ 
	return a;
	};

const Interpolator Interpolators::linear = []( double mix, double a, double b ) 
	{ 
	return ( 1.0 - mix ) * a + mix * b;
	};

const Interpolator Interpolators::sine = []( double mix, double a, double b ) 
	{ 
	return ( 1.0 - cos( pi * mix ) ) / 2.0 * ( b - a ) + a; 
	};


///Code duplication?
const Perturber Perturbers::normalDist = []( double x, double amount )
	{
	static std::default_random_engine rng( time(nullptr) );
	static std::normal_distribution<double> dist( 0, 1.0 );
	return x + dist(rng) * amount;
	};

const Perturber Perturbers::normalDistUp = []( double x, double amount )
	{
	static std::default_random_engine rng( time(nullptr) );
	static std::normal_distribution<double> dist( 0, 1.0 );
	return x + abs(dist(rng) * amount);
	};

const Perturber Perturbers::normalDistDown = []( double x, double amount )
	{
	static std::default_random_engine rng( time(nullptr) );
	static std::normal_distribution<double> dist( 0, 1.0 );
	return x - abs(dist(rng) * amount);
	};

const Perturber Perturbers::identity = []( double x, double amount ){ return x; };

std::array<uint8_t,3> HSVtoRGB( int H, double S, double V ) 
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

	return std::array<uint8_t,3>({ uint8_t((Rs + m) * 255),  uint8_t((Gs + m) * 255),  uint8_t((Bs + m) * 255) });
	}

template <typename T>
void write_bytes( unsigned char * target, T u )
	{
    static_assert( CHAR_BIT == 8, "CHAR_BIT != 8" );

    for( size_t k = 0; k < sizeof(T); k++ )
        target[k] = ( u >> 8*k ) & 0xFF;

    return;
	}

void write_bytes( unsigned char * target, const char * source )
	{
    static_assert( CHAR_BIT == 8, "CHAR_BIT != 8" );

    for( size_t k = 0; k < strlen(source); k++ )
        target[k] = source[k];

    return;
	}

bool writeBMP( const std::string & filename, const std::vector<std::vector<std::array<uint8_t,3>>> & data )
	{
	auto bail = []( const std::string & s )
		{
		std::cout << s << std::endl;
		return false;
		};

	std::ofstream file( filename, std::ios::binary );
	if( !file ) return bail( "Error opening " + filename + " to write BMP." );

	const uint32_t width  = data.size();
	if( width == 0 ) return bail( "Couldn't write BMP, data had no width" );
	const uint32_t height = data[0].size();
	const uint32_t dataSize = width * height * 4;
	const uint32_t headerSize = 14;
	const uint32_t DIBheaderSize = 40;

	unsigned char header[ headerSize ] = { 0 };
		write_bytes( header +  0, "BM" );
		write_bytes( header +  2, headerSize + DIBheaderSize + dataSize );
		write_bytes( header + 10, headerSize + DIBheaderSize );
	file.write( (const char *) header, headerSize );

	unsigned char DIBheader[ DIBheaderSize ] = { 0 };
		write_bytes( DIBheader +  0, DIBheaderSize	);  //size
		write_bytes( DIBheader +  4, width			);  //width
		write_bytes( DIBheader +  8, height			);  //height
		write_bytes( DIBheader + 12, (uint16_t) 1	);  //number of planes
		write_bytes( DIBheader + 14, (uint16_t) 32	);  //bits per pixel
		write_bytes( DIBheader + 16, (uint32_t) 0	);  //Compression type
		write_bytes( DIBheader + 20, (uint32_t) 0	);  //compressed size (0 for 0 type compression)
		write_bytes( DIBheader + 24, (uint32_t) 1	);  //x resolution ???
		write_bytes( DIBheader + 28, (uint32_t) 1	);  //y resolution ???
		write_bytes( DIBheader + 32, (uint32_t) 0	);  //number of used colors (0 for 2^n)
		write_bytes( DIBheader + 36, (uint32_t) 0	);  //number of important colors
	file.write( (const char *) DIBheader, DIBheaderSize );

	//y -> x to match write order
	for( int y = 0; y < height; ++y )
		for( int x = 0; x < width; ++x )
			{
			const uint32_t alpha = 0b00000000;
			file.write( (char *) &data[x][y][2], 1 ); 
			file.write( (char *) &data[x][y][1], 1 ); 
			file.write( (char *) &data[x][y][0], 1 ); 
			file.write( (char *) &alpha, 1 );
			}

	file.close();
	return true;
	}

}