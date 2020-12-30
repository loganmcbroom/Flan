#pragma once

#include <vector>
#include <array>
#include <string>

namespace xcdp {

/** Returns true if the system is little-endian. */
bool isLittleEndian();

/** Reverses the endianness of u. */
template <typename T>
T swapEndian( T u )
	{
    union
    {
        T u;
        unsigned char u8[ sizeof(T) ];
    } source, dest;

    source.u = u;

    for( size_t k = 0; k < sizeof(T); k++ )
        dest.u8[k] = source.u8[ sizeof(T) - k - 1 ];

    return dest.u;
	}

/** Forces u to be little-endian. */
template <typename T>
T makeLittleEndian( T u )
	{
	return isLittleEndian() ? u : swapEndian( u );
	}

/** Forces a little-endian u to match the endianness the enviornment uses. */
template <typename T>
T littleEndianToCurrent( T u )
	{
	return isLittleEndian() ? u : swapEndian( u );
	}

/** Write bytes (in little-endian) to target. */
template <typename T>
inline void writeBytes( unsigned char * target, T u )
	{
	u = makeLittleEndian( u );

	union
    {
        T u;
        unsigned char u8[ sizeof(T) ];
    } u2;

	u2.u = u;

    for( size_t k = 0; k < sizeof(T); k++ )
        target[k] = u2.u8[k];
	}

/** Write bytes from a c-string to target. */
inline void writeBytes( unsigned char * target, const char * source )
	{
    for( size_t k = 0; k < strlen(source); k++ )
        target[k] = source[k];
	}

/** DEPRICATED - bmps are now handled by bitmap_image
 *  Write a BMP stored in memory to disk.
 *  \param filename File path to write to. Should use the extension bmp.
 *  \param width The width of the output. Height will be determined by the size of data.
 *  \param data The data to write. The data should be in x-major order, meaning each vertical line should
 *      be written before the next. Each data point should be in RGB order.
 */
bool writeBMP( const std::string & filename, size_t width, const std::vector<std::array<uint8_t,3>> & data );

/** \cond */
struct RIFFData
    {
    RIFFData( uint32_t v ) : value( v ), numBytes( 4 ) {}
    RIFFData( uint16_t v ) : value( v ), numBytes( 2 ) {}
    uint32_t value;
    uint16_t numBytes; //Must be 2 or 4
    };
/** \endcond */
bool writeRIFF( const std::string & filename, const char type[4], const void * data, size_t dataSize, std::vector<RIFFData> format );

}