#include "flan/DSPUtility.h"

#include <algorithm>
#include <ranges>

#include "flan/FFTHelper.h"
#include "flan/Utility/vec2.h"
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

std::vector<vec2> find_peaks( std::function< float ( int ) > data, int size, int max_peaks, bool amp_order, bool interpolate ) 
    {
    if( max_peaks == -1 ) max_peaks = size / 2;

    std::vector<vec2> peaks;
    if( size < 2 ) return peaks;
    peaks.reserve( size );

    for( Frame run_start = 1; run_start < size - 1; )
        {
        const float run_value = data( run_start );
        Frame run_end = run_start;
        while( run_end + 1 < size && data( run_end + 1 ) == run_value )
            ++run_end;

        const Frame left_frame = run_start - 1;
        const Frame right_frame = run_end + 1;
        if( data( left_frame ) < run_value && right_frame < size && data( right_frame ) < run_value )
            {
            if( run_start == run_end )
                {
                if( interpolate )
                    {
                    const auto interpolatedData = parabolic_interpolation( data, run_start );
                    peaks.emplace_back( interpolatedData.first, interpolatedData.second );
                    }
                else
                    peaks.emplace_back( static_cast<float>( run_start ), run_value );
                }
            else
                {
                const float plateau_mean = ( right_frame + left_frame ) * 0.5f;
                peaks.emplace_back( interpolate ? plateau_mean : std::floor( plateau_mean ), run_value );
                }
            }

        run_start = run_end + 1;
        }

    if( amp_order ) // Sort peaks by descending y value if requested
        std::sort( FLAN_PAR_UNSEQ peaks.begin(), peaks.end(), []( auto & l, auto & r ){ return l.y() > r.y(); } );
    else // Otherwise sort by ascending x value
        std::sort( FLAN_PAR_UNSEQ peaks.begin(), peaks.end(), []( auto & l, auto & r ){ return l.x() < r.x(); } );

    // We only want this many peaks
    const size_t n_wanted_peaks = std::min( (size_t) max_peaks, peaks.size() );
    peaks.resize( n_wanted_peaks );

    return peaks;
    }

std::vector<vec2> find_peaks( const std::vector<float> & data, int max_peaks, bool amp_order, bool interpolate )
    {
    return find_peaks( [&data]( int i ){ return data[i]; }, data.size(), max_peaks, amp_order, interpolate );
    }

std::vector<vec2> find_valleys( std::function< float ( int ) > data, int size, int max_peaks, bool amp_order, bool interpolate )
    {
    std::vector<vec2> flippedPeaks = find_peaks( [&data]( int i ){ return -data( i ); }, size, max_peaks, amp_order, interpolate );
    std::for_each( FLAN_PAR_UNSEQ flippedPeaks.begin(), flippedPeaks.end(), []( vec2 & v ){ v.y() *= -1; } );
    return flippedPeaks;
    }
   
std::vector<vec2> find_valleys( const std::vector<float> & data, int max_peaks, bool amp_order, bool interpolate )
    {
    return find_valleys( [&data]( int i ){ return data[i]; }, data.size(), max_peaks, amp_order, interpolate );
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