#include "xcdp/DSPUtility.h"

#include <algorithm>

#include "xcdp/FFTHelper.h"
#include "xcdp/Utility/vec2.h"

namespace xcdp {

size_t autocorrelationFFTSize( size_t windowSize )
	{
	return std::pow( 2, 1 + (int) std::ceil( std::log2( windowSize ) ) );
	}

std::vector<float> autocorrelation( const float * signal, Frame n, std::shared_ptr<FFTHelper> fft ) 
	{
	// Forced power of 2 for faster fft, at least twice as big to avoid time-aliasing
	const Frame realBufferSize = autocorrelationFFTSize( n );

	//Allocate fft
	if( ! fft ) fft = std::make_shared<FFTHelper>( realBufferSize, true, true, false );

	// Read buffer data ( and zero fill padding ) into fftIn and execute fft
	std::copy( signal, signal + n, fft->realBuffer );
	std::fill( fft->realBegin() + n, fft->realEnd(), 0 );
	fft->r2cExecute();

	// Transform out buffer to squared magnitudes
	std::for_each( fft->complexBegin(), fft->complexEnd(), []( std::complex<float> & c ){ c = std::norm( c ); } );

	// FFT back and copy to output, normalizing
	fft->c2rExecute();
	std::vector<float> out( realBufferSize / 2 );
	std::transform( fft->realBegin(), fft->realBegin() + out.size(), out.begin(), [realBufferSize]( float f ){ return f / realBufferSize; } );

	return out;
	}

std::pair<float,float> parabolicInterpolation( float y0, float y1, float y2, int x1 )  
    {
    const float delta_x = 0.5f * ( y0 - y2 ) / ( y0 - 2 * y1 + y2 );
    const float X = x1 + delta_x;
    const float Y = y1 - 0.25f * ( y0 - y2 ) * delta_x;
    return { X, Y };
    }

std::pair<float,float> parabolicInterpolation( const std::vector<float> & d, int tau ) 
	{
    return parabolicInterpolation( d[tau-1], d[tau], d[tau+1], tau );
	}

std::vector<vec2> findPeaks( std::function< float ( int ) > data, int size, int maxPeaks, bool ampOrder, bool interpolate ) 
    {
    if( maxPeaks == -1 ) maxPeaks = size / 2;

    std::vector<vec2> peaks;
    peaks.reserve( size );

    if( size < 2 ) return peaks;

    const float scale = 1.0f;//range / float( size );

    int i = 0;

    // first check the boundaries:
    if (i+1 < size && data( i + 1 ) > data( i + 1 ) ) 
            peaks.push_back( vec2( i*scale, data( i ) ) );        

    while(true) 
        {
        while( i+1 < size-1 && data( i ) >= data( i + 1 ) ) ++i; // going down
        while( i+1 < size-1 && data( i ) <  data( i + 1 ) ) ++i; // now we're climbing

        int j = i;
        while( j+1 < size-1 && ( data( j ) == data( j + 1 ) ) ) ++j; // not anymore, go through the plateau

        // end of plateau, do we go up or down?
        if( j+1 < size-1 && data( j + 1 ) < data( j ) ) 
            { // going down again
            float resultBin = 0.0;
            float resultVal = 0.0;

            if( j != i ) 
                { // plateau peak between i and j
                if( interpolate ) resultBin = (i + j) * 0.5;
                else resultBin = i;
                resultVal = data( i );
                }
            else 
                { // interpolate peak at i-1, i and i+1
                if( interpolate ) 
                    {
                    auto binVal = parabolicInterpolation( data( j - 1 ), data( j ), data( j + 1 ), j );
                    resultBin = binVal.first;
                    resultVal = binVal.second;
                    }
                else 
                    {
                    resultBin = j;
                    resultVal = data( j );
                    }
                }

            const float resultPos = resultBin * scale;
            if( resultPos >= size )  break;
            peaks.push_back( vec2( resultPos, resultVal ) );
            }

        // nothing found, start loop again
        i = j;

        if( i+1 >= size-1 ) 
            { // check the one just before the last position
            if( i == size-2 && data( i - 1 ) < data( i ) && data( i + 1 ) < data( i ) ) 
                {
                float resultBin = 0.0;
                float resultVal = 0.0;
                if( interpolate ) 
                    {
                    auto binVal = parabolicInterpolation( data( j - 1 ), data( j ), data( j + 1 ), j );
                    resultBin = binVal.first;
                    resultVal = binVal.second;
                    }
                else 
                    {
                    resultBin = i;
                    resultVal = data( i );
                    }
                peaks.push_back( vec2( resultBin*scale, resultVal ) );
                }
            break;
            }
        }

    // check upper boundary here, so peaks are already sorted by position
    const float pos = size / scale;
    if( size - 2 < pos && pos <= size - 1 && data( size - 1 ) > data( size - 2 ) ) 
        peaks.push_back( vec2( ( size - 1 ) * scale, data( size - 1 ) ) );

    if( ampOrder ) std::sort( peaks.begin(), peaks.end(), []( auto & l, auto & r ){ return l.y() > r.y(); } );

    // we only want this many peaks
    const size_t nWantedPeaks = std::min( (size_t) maxPeaks, peaks.size() );
    peaks.resize( nWantedPeaks );

    return peaks;
    }

std::vector<vec2> findPeaks( const std::vector<float> & data, int maxPeaks, bool ampOrder, bool interpolate )
    {
    return findPeaks( [&data]( int i ){ return data[i]; }, maxPeaks, ampOrder, interpolate );
    }

std::vector<vec2> findValleys( std::function< float ( int ) > data, int size, int maxPeaks, bool ampOrder, bool interpolate )
    {
    std::vector<vec2> flippedPeaks = findPeaks( [&data]( int i ){ return -data( i ); }, size, maxPeaks, ampOrder, interpolate );
    std::for_each( flippedPeaks.begin(), flippedPeaks.end(), []( vec2 & v ){ v.y() *= -1; } );
    return flippedPeaks;
    }
   
std::vector<vec2> findValleys( const std::vector<float> & data, int maxPeaks, bool ampOrder, bool interpolate )
    {
    return findValleys( [&data]( int i ){ return data[i]; }, data.size(), maxPeaks, ampOrder, interpolate );
    }

float mean( std::function< float ( int ) > data, int n )
    {
    if( n == 0 ) return 0;
    float sum = 0;
    for( int i = 0; i < n; ++i )
        sum += data( i );

    return sum / n;
    }

vec2 meanAndStandardDeviation( std::function< float ( int ) > data, int n )
    {
    if( n == 0 ) return { 0, 0 };
    const float mean_c = mean( data, n );

    float diffSum = 0;
    for( int i = 0; i < n; ++i )
        {
        const float d = data( i ) - mean_c;
        diffSum += d * d;
        }

    const float variance = diffSum / n;
        
    return { mean_c, std::sqrt( variance ) };
    }

float mean( std::vector<float> data )
    {
    return mean( [&data]( int i ){ return data[i]; }, data.size() );
    }

vec2 meanAndStandardDeviation( std::vector<float> data )
    { 
    return meanAndStandardDeviation( [&data]( int i ){ return data[i]; }, data.size() );
    }

};