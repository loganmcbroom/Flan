#pragma once

//Any small misc code can live here 

#include <functional>
#include <vector>
#include <string>
#include <array>
#include <chrono>

namespace xcdp {

/** Interpolators describe how values between known data points should approximate those points near it.
 * These functions are passed values on [0,1] and are expected to return values on that same range.
 * OpenCL based algorithms sample Iterators passed to them, so rapidly changing user-defined
 * Iterators may not work as expected with those algorithms.
 */
using Interpolator = std::function< float ( float ) >;

namespace Interpolators 
{
/** Zero returning function. */
extern const Interpolator constant;

/** Input returning function. */
extern const Interpolator linear;

/** Input returning function. */
extern const Interpolator linear;

/** Smoothstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
extern const Interpolator smoothstep;

/** Smootherstep, as seen here: https://en.wikipedia.org/wiki/Smoothstep */
extern const Interpolator smootherstep;

/** Sine based interpolation. This moves from a valley to the following peak of a sinusoid. */
extern const Interpolator sine;
}

/** Distributions generate random data around 0 according to a distribution function.
 * These are used 
 */
using Distribution = std::function< float () >;

namespace Distributions 
{
/** Zero returning function. */
extern const Distribution identity;
	
/** Normal distribution with a standard deviation of 1. */
extern const Distribution normal;
}

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

/** Convert HSV (Hue-Saturation-Value) data to RGB (Red-Green-Blue) data. */
std::array<uint8_t,3> HSVtoRGB( int H, float S, float V );

/** Write a BMP stored in memory to disk.
 *  \param filename File path to write to. Should use the extension bmp.
 *  \param width The width of the output. Height will be determined by the size of data.
 *  \param data The data to write. The data should be in x-major order, meaning each vertical line should
 *      be written before the next. Each data point should be in RGB order.
 */
bool writeBMP( const std::string & filename, size_t width, const std::vector<std::array<uint8_t,3>> & data );

struct RIFFData
    {
    RIFFData( uint32_t v ) : value( v ), numBytes( 4 ) {}
    RIFFData( uint16_t v ) : value( v ), numBytes( 2 ) {}
    uint32_t value;
    uint16_t numBytes; //Must be 2 or 4
    };
bool writeRIFF( const std::string & filename, const char type[4], const void * data, size_t dataSize, std::vector<RIFFData> format );

class Timer
{
public:
    void start()
        {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
        }
    
    void stop()
        {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
        }
    
    float elapsedMilliseconds()
        {
        std::chrono::time_point<std::chrono::system_clock> endTime;
        
        if(m_bRunning)
            endTime = std::chrono::system_clock::now();
        else
            endTime = m_EndTime;
        
        return std::chrono::duration_cast<std::chrono::milliseconds>( endTime - m_StartTime ).count();
     }
    
    float elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0f;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};

}
