#include "flan/Utility/Bytes.h"

#include <cstdint>
#include <iostream>
#include <fstream>

namespace flan {

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
	std::vector<uint8_t> FMTChunk( FMTChunkSize, 0 );
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