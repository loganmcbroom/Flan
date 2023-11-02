#include "flan/Audio/Audio.h"

#include <iostream>
#include <set>
#include <limits>
#include <numeric>
#include <execution>

#include "flan/WindowFunctions.h"
#include "flan/FFTHelper.h"
#include "flan/DSPUtility.h"

#include "flan/Graph.h"

using namespace flan;

// Computes sum of squared differences for set number of lag values
static std::vector<float> compute_d( const float * in, Frame n )
	{
    const size_t fft_size = n; //std::pow( 2, 1 + (int) std::ceil( std::log2( n ) ) );
	auto fft = std::make_shared<FFTHelper>( fft_size, true, true, false );

    // Calculate power terms
	std::vector<float> power_terms( n / 2 );
    power_terms[0] = 0.0;
    for( int j = 0; j < n / 2; ++j ) 
        power_terms[0] += in[j] * in[j];
    for( int tau = 1; tau < power_terms.size(); ++tau ) 
        power_terms[tau] = power_terms[tau-1] - in[tau-1] * in[tau-1] + in[tau-1 + n / 2] * in[tau-1 + n / 2];  
		
    // Get modified autocorrelation. Note this needs two ffts, one of a full window and one of a zero padded half window
    
    // Get half fft
    std::vector<std::complex<float>> half_fft( fft->complex_buffer_size() );
    std::copy( in, in + n/2, fft->real_begin() );
    std::fill( fft->real_begin() + n / 2, fft->real_end(), 0 );
    fft->r2c_execute();
    std::copy( fft->complex_begin(), fft->complex_end(), half_fft.begin() );

    // Get full fft
    std::copy( in, in + n, fft->real_begin() );
    fft->r2c_execute();

    // Multiply
    for( int i = 0; i < fft->complex_buffer_size(); ++i ) 
		fft->get_complex_buffer()[i] *= std::conj( half_fft[i] );

    // IFFT
    fft->c2r_execute();
    
    // Find difference function d.
	std::vector<float> d( n / 2 );
    for( int j = 0; j < n / 2; ++j ) 
        d[j] = power_terms[0] + power_terms[j] - 2 * fft->get_real_buffer()[j] / n;

	return d;
	}

static std::vector<float> compute_d_prime( const float * in, Frame n )
	{
	// Get d
	std::vector<float> d = compute_d( in, n );

	// Transform d into d' via cumulative mean normalization
    d[0] = 1;
    double sum = 0;
    for( int tau = 1; tau < d.size(); ++tau ) 
		{
        sum += d[tau];
        if( sum == 0 ) d[tau] = 1;
        else d[tau] *= tau / sum;
		} 

	return d;
	}

//static Frame detectAutocorrelationPeak( float * data, Frame start, Frame end, float confidenceRation = 0.5 )
//	{
//	const float confidenceRequired = data[0] * confidenceRation;
//
//	// Discard initial peak. If ac data is never negative, we don't really care about it
//	Frame searchFrame = 0;
//	while( searchFrame < end && data[searchFrame] > 0 )  ++searchFrame; 
//	if( searchFrame < start ) searchFrame = start;
//
//	// Find potential peak frames using local maxima
//	std::set<Frame> peakFrames;
//	const Frame window_size = 128;
//	const Frame hopSize = window_size / 2;
//	for( ; searchFrame < end; searchFrame += hopSize )
//		{
//		// If local maxima is close to the center, it's a fine guess
//		auto localMaxFrameIter = std::max_element( data + searchFrame, data + searchFrame + window_size );
//		if( std::abs( std::distance( data + searchFrame + window_size / 2, localMaxFrameIter ) ) <= window_size / 4 )
//			{
//			const Frame localMaxFrame = std::distance( data, localMaxFrameIter );
//			peakFrames.insert( localMaxFrame );
//			}
//		}
//
//	// Find max peak among scaled peaks
//	auto distribution = [end]( Frame x )
//		{
//		return 1.0f - float( x ) / float( end );
//		};
//	Frame bestFrame = -1;
//	float peakAmp = 0;
//	for( auto f : peakFrames )
//		{
//		const float scaledAmp_c = data[f] * distribution( f );
//		if( scaledAmp_c > peakAmp && scaledAmp_c > confidenceRequired )
//			{
//			peakAmp = scaledAmp_c;
//			bestFrame = f;
//			}
//		}
//
//	return bestFrame;
//	}

std::vector<float> Audio::get_total_energy() const
	{
	std::vector<float> energies( get_num_channels() );
	flan::for_each_i( get_num_channels(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel ) {
		energies[channel] = std::accumulate( channel_begin( channel ), channel_end( channel ), 0.0f, 
		[]( float x, Sample s ){ return x + s * s; } );
		} );
	return energies;
	}

std::vector<float> Audio::get_energy_difference( const Audio & other ) const
	{	
	const Audio resampled =	get_sample_rate() == other.get_sample_rate() ? Audio::create_null() : other.resample( get_sample_rate() );
	const Audio * sr_correct_source = get_sample_rate() == other.get_sample_rate() ? &other : &resampled;
	return Audio::mix( std::vector<const Audio *>{ this, sr_correct_source }, {0, 0}, { 1, -1 } ).get_total_energy();
	}

float Audio::get_local_wavelength( Channel channel, Frame start, Frame window_size ) const
	{
	const float absolute_cutoff = 0.2f;

	if( is_null() ) return 0;

	// Get d' - Cumulative mean normalized difference function (search for YIN algorithm for details)
	const std::vector<float> d_prime = compute_d_prime( get_sample_pointer( channel, start ), window_size );
	
	// Get local minima of d_prime, sorted from largest to smallest
	const std::vector<vec2> minima = find_valleys( d_prime );
	
	// For wavelength finding, we want absolute minimum, as long as it's low enough
	vec2 best_minima = {0,-1};
	for( const vec2 & m : minima )
		if( m.y() < absolute_cutoff ) 
			best_minima = m;
    
    return best_minima.x();
	}

std::vector<float> Audio::get_local_wavelengths( Channel channel, Frame start, Frame end, Frame window_size, Frame hop, flan_CANCEL_ARG_CPP ) const
	{
	if( is_null() ) return std::vector<float>();


    if( end == -1 ) end = get_num_frames();

    std::vector<Frequency> out;
    for( Frame frame = start; frame + window_size < end; frame += hop )
		{
		flan_CANCEL_POINT( std::vector<float>() );
        out.push_back( get_local_wavelength( channel, frame, window_size ) );
		}
        
    return out;
	}

float Audio::get_average_wavelength( Channel channel, float min_active_ratio, float max_length_sigma, 
	Frame start, Frame end, Frame window_size, Frame hop, flan_CANCEL_ARG_CPP ) const
	{
	if( is_null() ) return 0;
	return get_average_wavelength( get_local_wavelengths( channel, start, end, window_size, hop, canceller ), min_active_ratio, max_length_sigma );
	}

float Audio::get_average_wavelength( const std::vector<float> & locals, float min_active_ratio, float max_length_sigma ) const
	{
	if( is_null() ) return 0;

	// Copy valid wavelengths
	std::vector<float> valid_lengths;

	const int num_valids = locals.size() - std::count( locals.begin(), locals.end(), -1 );
	if( num_valids <= min_active_ratio * locals.size() ) return -1; // Check that a sufficient amount of data had a detectable wavelength

	std::copy_if( locals.begin(), locals.end(), std::back_inserter( valid_lengths ), []( float x ){ return x != -1; } );

	const vec2 msd = mean_and_sd( valid_lengths );
	if( max_length_sigma != -1 && msd.y() > max_length_sigma ) return -1; // Check that the lengths didn't vary too much

	return msd.x();
	}

Frequency Audio::get_local_frequency( Channel channel, Frame start, Frame window_size, flan_CANCEL_ARG_CPP ) const
	{
	const float minimum_proximity_cutoff = 0.2f;
	const float absolute_cutoff = 0.15f;

    if( is_null() ) return 0;

	// Get d' - Cumulative mean normalized difference function (search for YIN algorithm for details)
	const std::vector<float> d_prime = compute_d_prime( get_sample_pointer( channel, start ), window_size );
	
	// Get local minima of d_prime, sorted from largest to smallest
	const std::vector<vec2> minima = find_valleys( d_prime, -1, true );

	// Filter out anything that isn't at least within 20% of the minimum minima.
	// Also filter anything that isn't below a global cutoff ( d_prime is normalized to some degree ).
	const float cutoff = std::min( minima.back().y() * ( 1.0f + minimum_proximity_cutoff ), absolute_cutoff );
	
	// For frequency finding, we want small wavelength to minimize octave errors
	vec2 best_minima = { float( window_size ),-1 };
	for( const vec2 & m : minima )
		if( m.y() < cutoff && m.x() < best_minima.x() ) 
			best_minima = m;
    
    return float( get_sample_rate() ) / best_minima.x();
    }

std::vector<Frequency> Audio::get_local_frequencies( Channel channel, Frame start, Frame end, Frame window_size, Frame hop, flan_CANCEL_ARG_CPP ) const
    {
	if( is_null() ) return std::vector<Frequency>();

    if( end == -1 ) end = get_num_frames();

    std::vector<Frequency> out( (end - window_size - start) / hop ); // This will round toward zero
	flan::for_each_i( out.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int i )
		{
		const Frame frame = start + i * hop;
        out[i] = get_local_frequency( channel, frame, window_size, canceller );
		} );
        
    return out;
    }

Function<Second, Amplitude> Audio::get_amplitude_envelope(
	Second window_width
	) const
	{
	if( is_null() ) return 0;

	if( window_width <= 0 ) return 0;

	Audio mono = convert_to_mono();

	// Rectify
	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		if( mono.get_sample( 0, frame ) < 0 )
			mono.get_sample( 0, frame ) *= -1;

	const fFrame window_width_frames = time_to_frame( window_width );

	Audio hann_window = Audio::create_empty_with_frames( window_width_frames, 1, get_sample_rate() );
	for( Frame i = 0; i < hann_window.get_num_frames(); ++i )
		hann_window.get_sample( 0, i ) = Windows::hann( float( i ) / ( window_width_frames - 1 ) );

	auto convolved = mono.convolve( hann_window );
	
	std::vector<float> ys = std::move( convolved.get_buffer() );
	return [ys = std::move( ys ), sr = get_sample_rate() ]( float t )
		{
		// Get the image sample index as a float
		const float x = t * sr;
		if( 0 <= x && x < ys.size() - 1 )
			{
			Frame x1 = std::floor( x );
			float y1 = ys[x1];
			float y2 = ys[x1 + 1];
			return y1 + ( y2 - y1 ) * ( x - x1 );
			}
		else return 0.0f;
		};

	// if( is_null() ) return 0;

	// Audio out = convert_to_mono();

	// // Rectify
	// for( Frame frame = 0; frame < get_num_frames(); ++frame )
	// 	if( out.get_sample( 0, frame ) < 0 )
	// 		out.get_sample( 0, frame ) *= -1;

	// out = out.filter_1pole_lowpass( 10, 4 );

	// // Needing this is a c++ issue, not stealing the vector makes the lamba capture initialization get flagged as trying to call Audio copy ctor
	// std::vector<Sample> buffer = std::move( out.get_buffer() );

	// return [buffer = std::move( buffer ), sample_rate = get_sample_rate()]( Second t )
	// 	{
	// 	const Frame frame = t * sample_rate;
	// 	if( 0 <= frame && frame < buffer.size() )
	// 		return buffer[frame];
	// 	return 0.0f;
	// 	};
	}
	
Function<Second, Amplitude> Audio::get_frequency_envelope( 
	) const
	{
	const int hop_size = 128;
	const std::vector<Frequency> frequencies = convert_to_mono().get_local_frequencies( 0, 0, -1, 2048, hop_size );

	return [frequencies = std::move( frequencies ), sample_rate = get_sample_rate(), hop_size]( Second t )
		{
		const float x = t * sample_rate / hop_size;
		if( 0 <= x && x < frequencies.size() - 1 )
			{
			const Frame x1 = std::floor( x );
			const float y1 = frequencies[x1];
			const float y2 = frequencies[x1 + 1];
			const float m = ( y2 - y1 ) / 1;
			return y1 + m * ( x - x1 ); // Point-slope form for lerp
			}
		else return 0.0f;
		};
	}
	