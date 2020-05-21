#include "xcdp/Utility.h"

#include <random>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstring>

#include <bitset> //temp include

static float pi = std::acos( -1.0f );

namespace xcdp {

const Interpolator Interpolators::constant = []( float x ) 
	{ 
	return 0.0f;
	};

const Interpolator Interpolators::linear = []( float x ) 
	{ 
	return x;
	};

const Interpolator Interpolators::sine = []( float x ) 
	{ 
	return ( 1.0f - cos( pi * x ) ) / 2.0f; 
	};

const Distribution Distributions::normal = []()
	{
	static std::default_random_engine rng( time( nullptr ) );
	static std::normal_distribution<float> dist( 0, 1.0 );
	return dist(rng);
	};

const Distribution Distributions::identity = []()
	{ 
	return 0; 
	};

std::array<uint8_t,3> HSVtoRGB( int H, float S, float V ) 
	{
	float C = S * V;
	float X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	float m = V - C;
	float Rs, Gs, Bs;

		 if( H >= 0   && H < 60  ) { Rs = C; Gs = X; Bs = 0; }
	else if( H >= 60  && H < 120 ) { Rs = X; Gs = C; Bs = 0; }
	else if( H >= 120 && H < 180 ) { Rs = 0; Gs = C; Bs = X; }
	else if( H >= 180 && H < 240 ) { Rs = 0; Gs = X; Bs = C; }
	else if( H >= 240 && H < 300 ) { Rs = X; Gs = 0; Bs = C; }
	else						   { Rs = C; Gs = 0; Bs = X; }

	return std::array<uint8_t,3>({ uint8_t((Rs + m) * 255),  uint8_t((Gs + m) * 255),  uint8_t((Bs + m) * 255) });
	}

bool isLittleEndian()
	{
    short int number = 0x1;
    char *numPtr = (char*)&number;
    return (numPtr[0] == 1);
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
		writeBytes( header +  0, "BM" );
		writeBytes( header +  2, headerSize + DIBheaderSize + dataSize );
		writeBytes( header + 10, headerSize + DIBheaderSize );
	file.write( (const char *) header, headerSize );

	unsigned char DIBheader[ DIBheaderSize ] = { 0 };
		writeBytes( DIBheader +  0, DIBheaderSize	);  //size
		writeBytes( DIBheader +  4, width			);  //width
		writeBytes( DIBheader +  8, height			);  //height
		writeBytes( DIBheader + 12, (uint16_t) 1	);  //number of planes
		writeBytes( DIBheader + 14, (uint16_t) 32	);  //bits per pixel
		writeBytes( DIBheader + 16, (uint32_t) 0	);  //Compression type
		writeBytes( DIBheader + 20, (uint32_t) 0	);  //compressed size (0 for 0 type compression)
		writeBytes( DIBheader + 24, (uint32_t) 1	);  //x resolution ???
		writeBytes( DIBheader + 28, (uint32_t) 1	);  //y resolution ???
		writeBytes( DIBheader + 32, (uint32_t) 0	);  //number of used colors (0 for 2^n)
		writeBytes( DIBheader + 36, (uint32_t) 0	);  //number of important colors
	file.write( (const char *) DIBheader, DIBheaderSize );

	//y -> x to match write order
	for( int y = 0; y < height; ++y )
		for( int x = 0; x < width; ++x )
			{
			const uint8_t alpha = 0b00;
			file.write( (char *) &data[x][y][2], 1 ); 
			file.write( (char *) &data[x][y][1], 1 ); 
			file.write( (char *) &data[x][y][0], 1 ); 
			file.write( (char *) &alpha, 1 );
			}

	file.close();
	return true;
	}

}