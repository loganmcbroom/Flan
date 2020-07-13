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

const Interpolator Interpolators::midpoint = []( float x ) 
	{ 
	return 0.5f;
	};

const Interpolator Interpolators::nearest = []( float x ) 
	{ 
	return std::round( x );
	};

const Interpolator Interpolators::floor = []( float x ) 
	{ 
	return 0.0f;
	};

const Interpolator Interpolators::ceil = []( float x ) 
	{ 
	return 1.0f;
	};

const Interpolator Interpolators::linear = []( float x ) 
	{ 
	return x;
	};

const Interpolator Interpolators::smoothstep = []( float x ) 
	{ 
	return x * x * ( 3.0f - 2.0f * x );
	};

const Interpolator Interpolators::smootherstep = []( float x ) 
	{ 
	return x * x * x * ( x * ( x * 6.0f - 15.0f ) + 10.0f );
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

bool writeBMP( const std::string & filename, size_t width, const std::vector<std::array<uint8_t,3>> & data )
	{
	auto bail = []( const std::string & s )
		{
		std::cout << s << std::endl;
		return false;
		};

	if( width == 0 ) return bail( "Couldn't write BMP with width 0." );
	if( data.size() == 0 ) return bail( "Couldn't write BMP, no data was sent." );

	std::ofstream file( filename, std::ios::binary );
	if( !file ) return bail( "Error opening " + filename + " to write BMP." );

	const uint32_t height = data.size() / width;
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
			file.write( (char *) &data[x * height + y][2], 1 ); 
			file.write( (char *) &data[x * height + y][1], 1 ); 
			file.write( (char *) &data[x * height + y][0], 1 ); 
			file.write( (char *) &alpha, 1 );
			}

	file.close();
	return true;
	}

bool writeRIFF( const std::string & filename, const char type[4], const void * data, size_t dataSize, std::vector<RIFFData> format )
	{
	std::ofstream file( filename, std::ios::binary );
	if( !file )
		{
		std::cout << "Error opening " + filename + " to write RIFF.\n";
		return false;
		}

	const uint32_t RIFFChunkSize = 12;
	uint32_t FMTChunkSize = 8;
	for( auto & i : format ) FMTChunkSize += i.numBytes;
	const uint32_t dataChunkSize = uint32_t( 8 + dataSize );

	unsigned char * write;

	//Write RIFF chunk
	unsigned char RIFFChunk[ RIFFChunkSize ] = { 0 };
	write = RIFFChunk;
		writeBytes( write, "RIFF"		); write += 4; //RIFF ID 
		writeBytes( write, (uint32_t) 4 ); write += 4; //RIFF Chunk Size
		writeBytes( write, type		    ); write += 4;  //File Type
	file.write( (const char *) RIFFChunk, RIFFChunkSize );

	//Write format chunk
	std::vector<unsigned char> FMTChunk( FMTChunkSize, 0 );
	write = FMTChunk.data();
		writeBytes( write, "fmt "						); write += 4;
		writeBytes( write, (uint32_t) FMTChunkSize - 8	); write += 4;
		for( int i = 0; i < format.size(); ++i )
			{
			if( format[i].numBytes == 2 ) writeBytes( write, uint16_t( format[i].value ) ); 
			else						  writeBytes( write, uint32_t( format[i].value ) ); 
			write += format[i].numBytes;
			}
	file.write( (const char *) FMTChunk.data(), FMTChunkSize );

	//Write data chunk
	unsigned char dataChunkInfo[ 8 ] = { 0 };
	write = dataChunkInfo;
		writeBytes( write, "data"						); write += 4;
		writeBytes( write, (uint32_t) dataChunkSize - 8	); write += 4;
	file.write( (const char *) dataChunkInfo, 8 );

	file.write( (const char *) data, dataSize );

	file.close();
	return true;
	}

}