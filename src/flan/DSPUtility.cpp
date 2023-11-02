#include "flan/DSPUtility.h"

#include <algorithm>
#include <ranges>

#include "flan/FFTHelper.h"
#include "flan/Utility/vec2.h"
#include "flan/Utility/iota_iter.h"
#include "flan/Utility/execution.h"

namespace flan {

// std::vector<float> autocorrelation( const float * signal, Frame n, std::shared_ptr<FFTHelper> fft ) 
// 	{
// 	// Forced power of 2 for faster fft, at least twice as big to avoid time-aliasing
// 	const Frame real_buffer_size = autocorrelation_fft_size( n );

// 	//Allocate fft
// 	if( ! fft ) fft = std::make_shared<FFTHelper>( real_buffer_size, true, true, false );

// 	// Read buffer data ( and zero fill padding ) into fftIn and execute fft
// 	std::copy( std::execution::par_unseq, signal, signal + n, fft->real_begin() );
// 	std::fill( std::execution::par_unseq, fft->real_begin() + n, fft->real_end(), 0 );
// 	fft->r2c_execute();

// 	// Transform out buffer to squared magnitudes
// 	std::for_each( std::execution::par_unseq, fft->complex_begin(), fft->complex_end(), []( std::complex<float> & c ){ c = std::norm( c ); } );

// 	// FFT back and copy to output, normalizing
// 	fft->c2r_execute();
// 	std::vector<float> out( real_buffer_size / 2 );
// 	std::transform( std::execution::par_unseq, fft->real_begin(), fft->real_begin() + out.size(), out.begin(), [real_buffer_size]( float f ){ return f / real_buffer_size; } );

// 	return out;
// 	}

std::pair<float,float> parabolic_interpolation( float y0, float y1, float y2, int x1 )  
    {
    const float delta_x = 0.5f * ( y0 - y2 ) / ( y0 - 2 * y1 + y2 );
    const float X = x1 + delta_x;
    const float Y = y1 - 0.25f * ( y0 - y2 ) * delta_x;
    return { X, Y };
    }

std::pair<float,float> parabolic_interpolation( const std::vector<float> & d, int tau ) 
	{
    return parabolic_interpolation( d[tau-1], d[tau], d[tau+1], tau );
	}

std::pair<float,float> parabolic_interpolation( std::function< float ( int ) > f, int tau ) 
	{
    return parabolic_interpolation( f(tau-1), f(tau), f(tau+1), tau );
	}

std::vector<vec2> find_peaks( std::function< float ( int ) > data, int size,  int maxPeaks, bool ampOrder, bool interpolate ) 
    {
    if( maxPeaks == -1 ) maxPeaks = size / 2;

    std::vector<vec2> peaks;
    if( size < 2 ) return peaks;
    peaks.reserve( size );

    // For each data point excluding first and last, check if that point is a peak
    std::mutex mutex;
    std::for_each( FLAN_PAR_SEQ iota_iter( 1 ), iota_iter( size - 1 ), [&]( const Frame frame )
        {
        const float frameVal = data(frame);

        // Find the first data in one direction that is less than the data for frame, OR return -1 if the data ever goes up (aka frame is not a peak)
        auto finder = [&]( const bool goRight )
            {
            Frame newFrame = frame;
            while( true )
                {
                goRight? newFrame++ : newFrame--;
                if( newFrame < 0 || size <= newFrame ) return -1; // We got to the edge of the data and it just stayed constant, return -1

                const float newVal = data(newFrame);
                if( newVal > frameVal ) return -1; // Frame was not a peak, return -1
                if( newVal < frameVal ) break; // No need to go farther, we found it
                }
            return newFrame;
            };

        // Find first left data point that is less than frame height, or return 
        const Frame leftFrame = finder( false );
        if( leftFrame == -1 ) return;

        // Find first right data point that less than frame height, or return
        const Frame rightFrame = finder( true );
        if( rightFrame == -1 ) return;

        if( ( rightFrame - leftFrame ) > 2 ) // If we were in a plateau
            {
            // If frame isn't in the middle of right and left we can just leave
            // If we didn't do this check, all the frames in a plateau would be marked as peaks when we only want one
            // It isn't super important that the middle of the plateau is the frame used, but it makes the most sense
            const float plateauMean = ( rightFrame + leftFrame ) * 0.5;
            if( frame != std::floor( plateauMean ) ) return;

            std::lock_guard<std::mutex> lock(mutex);
            peaks.emplace_back( interpolate ? plateauMean : frame, data(frame) );
            }
        else // Not in a plateau
            {
            if( interpolate )
                {
                std::lock_guard<std::mutex> lock(mutex);
                const auto interpolatedData = parabolic_interpolation( data, frame );
                peaks.emplace_back( interpolatedData.first, interpolatedData.second );
                }
            else
                {
                std::lock_guard<std::mutex> lock(mutex);
                peaks.emplace_back( frame, data(frame) );
                }
            }
        } );

    if( ampOrder ) // Sort peaks by descending y value if requested
        std::sort( FLAN_PAR_UNSEQ peaks.begin(), peaks.end(), []( auto & l, auto & r ){ return l.y() > r.y(); } );
    else // Otherwise sort by ascending x value
        std::sort( FLAN_PAR_UNSEQ peaks.begin(), peaks.end(), []( auto & l, auto & r ){ return l.x() < r.x(); } );

    // We only want this many peaks
    const size_t nWantedPeaks = std::min( (size_t) maxPeaks, peaks.size() );
    peaks.resize( nWantedPeaks );

    return peaks;
    }

std::vector<vec2> find_peaks( const std::vector<float> & data, int maxPeaks, bool ampOrder, bool interpolate )
    {
    return find_peaks( [&data]( int i ){ return data[i]; }, data.size(), maxPeaks, ampOrder, interpolate );
    }

std::vector<vec2> find_valleys( std::function< float ( int ) > data, int size, int maxPeaks, bool ampOrder, bool interpolate )
    {
    std::vector<vec2> flippedPeaks = find_peaks( [&data]( int i ){ return -data( i ); }, size, maxPeaks, ampOrder, interpolate );
    std::for_each( FLAN_PAR_UNSEQ flippedPeaks.begin(), flippedPeaks.end(), []( vec2 & v ){ v.y() *= -1; } );
    return flippedPeaks;
    }
   
std::vector<vec2> find_valleys( const std::vector<float> & data, int maxPeaks, bool ampOrder, bool interpolate )
    {
    return find_valleys( [&data]( int i ){ return data[i]; }, data.size(), maxPeaks, ampOrder, interpolate );
    }

float mean( const std::vector<float> & data )
    {
    if( data.size() == 0 ) return 0;
    float sum = 0;
    for( float x : data )
        sum += x;

    return sum / data.size();
    }

float mean( std::function< float ( int ) > data, int n )
    {
    // Sample function data
    if( n <= 0 ) return 0;
    float sum = 0;
    for( const int i : std::views::iota( 0, n ) )
        sum += data(i);

    return sum / n;
    }

vec2 mean_and_sd( std::function< float ( int ) > data, int n )
    { 
    if( n <= 0 ) return { 0, 0 };
    const float mean_c = mean( data, n );

    float diffSum = 0;
    for( const int i : std::views::iota( 0, n ) )
        {
        const float d = data(i) - mean_c;
        diffSum += d * d;
        }

    const float variance = diffSum / n;
        
    return { mean_c, std::sqrt( variance ) };
    }

vec2 mean_and_sd( const std::vector<float> & data  )
    {
    return mean_and_sd( [&data]( int i ){ return data[i]; }, data.size() );
    }

};