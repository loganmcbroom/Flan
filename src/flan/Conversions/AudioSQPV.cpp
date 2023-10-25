// #include "flan/Audio/Audio.h"
// #include "flan/SQPV/SQPV.h"

// #include <algorithm>
// #include <cassert>
// #include <cmath>
// #include <complex>
// #include <cstdlib>
// #include <deque>
// #include <map>
// #include <iostream>

// #include "flan/phase_vocoder.h"

// using namespace flan;

// // Container for three complex numbers with component-wise arithmetic overloads
// template<typename T>
// struct Twiddle 
// 	{
// 	using cpx = std::complex<T>;
// 	std::array<cpx, 3> data;
// 	cpx & operator[]( size_t i ) { return data[i]; }
// 	const cpx & operator[]( size_t i ) const { return data[i]; }

// 	Twiddle operator+( const Twiddle & other ) const
// 		{
// 		return {
// 			data[0] + other.data[0],
// 			data[1] + other.data[1],
// 			data[2] + other.data[2],
// 			};
// 		};
	
// 	Twiddle operator*( const Twiddle & other ) const
// 		{
// 		return {
// 			data[0] * other.data[0],
// 			data[1] * other.data[1],
// 			data[2] * other.data[2],
// 			};
// 		};

// 	Twiddle operator+( cpx x ) const
// 		{
// 		return {
// 			data[0] + x,
// 			data[1] + x,
// 			data[2] + x,
// 			};
// 		}
// 	};

// template<typename T>
// static inline Twiddle<T> make_twiddle( T quality, T period )
// 	{
// 	return {
// 		std::polar( T( 1 ), pi2 * ( quality - 1 ) / period ),
// 		std::polar( T( 1 ), pi2 * ( quality + 0 ) / period ),
// 		std::polar( T( 1 ), pi2 * ( quality + 1 ) / period ),
// 		};
// 	}

// SQPV Audio::convert_to_SQPV( std::pair<Frequency, Frequency> bandwidth, Bin bins_per_octave ) const
// 	{
// 	/*
// 	See "Sliding With A Constant-Q" - https://www.dafx.de/paper-archive/2008/papers/dafx08_63.pdf - for details on this algorithm.
// 	*/

// 	
// 	SQPV::Format format;
// 	format.num_channels = get_num_channels();
// 	format.num_frames = get_num_frames();
// 	format.bins_per_octave = bins_per_octave;
// 	format.bandwidth = bandwidth;
// 	format.sample_rate = get_sample_rate();
// 	SQPV out( format );

// 	std::pair<double, double> window = std::make_pair( 0.5, -0.5 ); // hann window
// 	const Frame N = out.getPeriod( 0 );	
// 	const std::complex<double> fiddle = std::polar( 1.0, double( -pi2 * out.getQ() ) );

// 	for( Channel channel = 0; channel < get_num_channels(); ++channel )
// 		{
// 		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.get_num_bins() ), [&]( Bin bin )
// 			{
// 			const Frame N_k = out.getPeriod( bin );
// 			const Twiddle twiddle = make_twiddle<double>( out.getQ(), N_k );

// 			const Frequency binFrequency = out.getBinFrequency( bin );
// 			Twiddle<double> F_tk = {0,0,0};
// 			double phase_buffer = 0;
// 			if( bin == out.get_num_bins() - 1 )
// 				std::cout << "f";
			
// 			for( Frame frame = std::floor(-N_k/2-1); frame < out.get_num_frames(); ++frame )
// 				{
// 				auto get_sample_s = [&]( Frame frame ){ return frame < 0 || out.get_num_frames() <= frame ? 0.0f : get_sample( channel, frame ); };
// 				const double new_sample = get_sample_s( frame + N_k/2.0f );
// 				const double old_sample = get_sample_s( frame - N_k/2.0f );

// 				const auto fiddled = ( fiddle * new_sample - old_sample )/ double( N_k );
// 				F_tk = twiddle * ( F_tk + fiddled  ); // (6)

// 				if( frame >= 0 )
// 					{
// 					// If you are reading the paper on this step note that a factor of j is missing for several steps of the derivation, but it is correct.
// 					const std::complex<float> F_tk_windowed = static_cast<std::complex<float>>( 
// 						window.first * F_tk[1] + window.second / 2.0 * ( F_tk[0] + F_tk[2] ) ); // (10)
// 					const MF mf = phase_vocoder( phase_buffer, F_tk_windowed, binFrequency, out.get_analysis_rate(), out.get_sample_rate() );
// 					out.getMP( channel, frame, bin ) = { mf.m, out.frequencyToPitch( mf.f ) };
// 					}
// 				}
// 			} );
// 		}

// 	return out;
// 	}

// SQPV Audio::convert_to_ms_SQPV( std::pair<Frequency, Frequency> bandwidth, Bin bins_per_octave ) const
// 	{
// 	return convert_to_mid_side().convert_to_SQPV( bandwidth, bins_per_octave );
// 	}

// Audio SQPV::convert_to_audio() const
// 	{
// 	
// 	// Compute twiddles
// 	std::vector<std::complex<float>> twiddles( get_num_bins() );
// 	for( Bin bin = 0; bin < get_num_bins(); ++bin )
// 		twiddles[bin] = std::polar( 1.0f, pi2 * getQ() / getPeriod( bin ) ); // Rotation in (11)

// 	// Setup done
// 	Audio::Format format;
// 	format.num_channels = get_num_channels();
// 	format.num_frames = get_num_frames();
// 	format.sample_rate = get_sample_rate();
// 	Audio out( format );

// 	for( Channel channel = 0; channel < get_num_channels(); ++channel )
// 		{
// 		// Invert phase vocoding
// 		std::vector<std::complex<float>> cpxBuffer( get_num_frames() * get_num_bins(), 0 );
// 		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )		
// 			{
// 			double phase_buffer = 0;
// 			const Radian phase_diff_0 = pitchToFrequency( getMP( channel, 0, bin ).p ) / get_analysis_rate() * pi2;
// 			for( Frame frame = 0; frame < get_num_frames(); ++frame )
// 				{
// 				const MP mp = getMP( channel, frame, bin );
// 				cpxBuffer[get_buffer_pos( 0, frame, bin )] = inverse_phase_vocoder( phase_buffer, { mp.m, pitchToFrequency( mp.p ) }, get_analysis_rate() );
// 				}
// 			} );

// 		// Invert scqt
// 		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
// 			{
// 			Sample sample = 0;
// 			for( Bin bin = 0; bin < get_num_bins(); ++bin )
// 				sample += ( cpxBuffer[ get_buffer_pos( 0, frame, bin ) ] * twiddles[bin] ).real();
// 			out.get_sample( channel, frame ) = sample;
// 			} );
// 		}

// 	return out;
// 	}

// Audio SQPV::convert_to_lr_audio() const
// 	{
// 	return convert_to_audio().convert_to_left_right();
// 	}