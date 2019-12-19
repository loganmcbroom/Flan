#pragma once

//Any small misc code can live here 

#include <functional>
#include <vector>
#include <string>
#include <array>
#include <chrono>

namespace xcdp {

typedef std::function< double ( double mix, double a, double b ) > Interpolator;
namespace Interpolators 
	{
	extern const Interpolator 
		constant, 
		linear, 
		sine;
	}

typedef std::function< double ( double, double ) > Perturber;
namespace Perturbers {
extern const Perturber
	identity, 
	normalDist, 
	normalDistUp, 
	normalDistDown;
}

bool isLittleEndian();

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

template <typename T>
T makeLittleEndian( T u )
	{
	return isLittleEndian() ? u : swapEndian( u );
	}

template <typename T>
T littleEndianToCurrent( T u )
	{
	return isLittleEndian() ? u : swapEndian( u );
	}

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

inline void writeBytes( unsigned char * target, const char * source )
	{
    for( size_t k = 0; k < strlen(source); k++ )
        target[k] = source[k];
	}

std::array<uint8_t,3> HSVtoRGB( int H, double S, double V );
bool writeBMP( const std::string & filename, const std::vector<std::vector<std::array<uint8_t,3>>> & data );

class elapsedtimer
{
public:
	elapsedtimer()
	{
		m_begin = std::chrono::steady_clock::now();
	}
	long long elapsed()
	{
		auto now = std::chrono::steady_clock::now();
		return std::chrono::duration_cast<std::chrono::milliseconds>(now - m_begin).count();
	}
private:
	std::chrono::steady_clock::time_point m_begin;
};

}
