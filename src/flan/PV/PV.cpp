#include "flan/PV/PV.h"

#include <iostream>
#include <complex>
#include <array>
#include <numeric>
#include <ranges>

#include "flan/FFTHelper.h"
#include "flan/Utility/buffer_access.h"
#include "flan/WindowFunctions.h"

using namespace std::ranges;

namespace flan {

PV PV::get_frame( Second time ) const
	{
	if( is_null() ) return PV();

	const float selectedFrame = std::clamp( time_to_frame( time ), 0.0f, float( get_num_frames() - 1 ) );

	auto format = get_format();
	format.num_frames = 1;
	PV out( format );

	for( Channel channel = 0; channel < out.get_num_channels(); ++channel )
		for( Bin bin = 0; bin < out.get_num_bins(); ++bin )
			out.get_MF( channel, 0, bin ) = getBinInterpolated( channel, selectedFrame, bin );

	return out;
	}

MF PV::getBinInterpolated( Channel channel, float frame, float bin, Interpolator i ) const
	{
	const std::array< MF, 4 > p = 
		{
		get_MF( channel, floor( frame ), floor( bin ) ),
		get_MF( channel, ceil ( frame ), floor( bin ) ),
		get_MF( channel, ceil ( frame ), ceil ( bin ) ),
		get_MF( channel, floor( frame ), ceil ( bin ) )
		};

	const float l = i( frame - floor( frame ) );
	const float m = i( bin   - floor( bin   ) );
	const float mI = 1.0f - m;
	const float lI = 1.0f - l;

	return {
		   mI * ( lI * p[0].m + l * p[1].m ) + m * ( lI * p[3].m + l * p[2].m ),
		   mI * ( lI * p[0].f + l * p[1].f ) + m * ( lI * p[3].f + l * p[2].f )
		   };
	}

MF PV::getBinInterpolated( Channel channel, float frame, Bin bin, Interpolator i ) const
	{
	const MF & l = get_MF( channel, floor( frame ), bin );
	const MF & h = get_MF( channel, ceil ( frame ), bin );

	const float mix = i( frame - floor( frame ) );

	return { 
		   ( 1.0f - mix ) * l.m + mix * h.m,
		   ( 1.0f - mix ) * l.f + mix * h.f
		   };
	}

MF PV::getBinInterpolated( Channel channel, Frame frame, float bin, Interpolator i ) const
	{
	const auto l = get_MF( channel, frame, floor( bin ) );
	const auto h = get_MF( channel, frame, ceil ( bin ) );

	const float mix = i( bin - floor( bin ) );

	return  { 
			( 1.0f - mix ) * l.m + mix * h.m,
			( 1.0f - mix ) * l.f + mix * h.f
			};
	}

//========================================================================
//	Selection
//========================================================================

PV PV::select( 
	Second length, 
	const Function<TF, TF> & selector
	) const
	{
	if( is_null() ) return PV();
	if( length <= 0.0f ) return PV();

	auto format = get_format();
	format.num_frames = time_to_frame( length );
	PV out( format );

	const auto selector_sampled = out.sample_function_over_domain( selector );

	for( Channel channel = 0; channel < out.get_num_channels(); ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( out.get_num_frames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < out.get_num_bins(); ++bin )
				{
				const TF s = selector_sampled[ get_buffer_pos( 0, frame, bin ) ];
				const Frame selection_frame = time_to_frame( s.t );
				const Bin selection_bin = frequency_to_bin( s.f );

				if( selection_frame  < 0 || get_num_frames() - 1 <= selection_frame ||
				    selection_bin < 0 || get_num_bins()   - 1 <= selection_bin )
					continue;

				MF selected_mf = get_MF( channel, selection_frame, selection_bin );
				if( s.f > 1 )
					selected_mf.f *= bin_to_frequency( bin ) / s.f;
				out.get_MF( channel, frame, bin ) = selected_mf;
				}
			} );

	return out;
	}

PV PV::freeze( 		
	const std::vector<Second> & pause_times,
	const std::vector<Second> & pause_lengths ) const
	{
	if( is_null() ) return PV();

	if( pause_lengths.size() != pause_times.size() ) 
		{
		std::cerr << "Error in flan::PV::freeze: pause_times and pause_lengths were not the same size.";
		return PV();
		}

	// Combine timing info into single vector for easy sorting later
	std::vector<std::array<Second, 2>> timing( pause_times.size() );
	std::transform( pause_times.begin(), pause_times.end(), pause_lengths.begin(), timing.begin(), []( Second a, Second b ){ return std::array<Second, 2>{ a, b }; } );

	using FramePair = std::array<Frame,2>;

	// Convert times into frames
	std::vector<FramePair> timing_frames( timing.size() );
	std::ranges::transform( timing, timing_frames.begin(), [this]( const std::array<Second,2> & ts )
		{ 
		return FramePair{ 
			std::clamp( Frame( time_to_frame( ts[0] ) ), 0, Frame( get_num_frames() - 1 ) ), 
			std::max  ( Frame( time_to_frame( ts[1] ) ), 0 ) };
		});

	// Sort by occurance frame
	std::sort( FLAN_PAR_UNSEQ timing_frames.begin(), timing_frames.end(), []( FramePair & a, FramePair & b ){ return a[0] < b[0]; } );

	// Remove simultaneous events
	auto timingEnd = std::unique( timing_frames.begin(), timing_frames.end(), []( const FramePair & a, const FramePair & b ) {
		return a[0] == b[0];
		} );
	timing_frames.erase( timingEnd, timing_frames.end() ); 

	float total_freeze_frames = 0;
	for( auto & i : timing_frames ) total_freeze_frames += i[1];

	auto format = get_format();
	format.num_frames += total_freeze_frames;
	PV out( format );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		uint32_t timing_index = 0;
		for( Frame in_frame = 0, out_frame = 0; in_frame < get_num_frames(); ++in_frame )
			{
			// If it's time to freeze, freeze
			if( timing_index < timing_frames.size() && in_frame == timing_frames[timing_index][0] )
				{
				for( Frame freeze_frame = 0; freeze_frame < timing_frames[timing_index][1]; ++freeze_frame )
					{
					for( Bin bin = 0; bin < get_num_bins(); ++bin )
						out.set_MF( channel, out_frame, bin, get_MF( channel, in_frame, bin ) );
					++out_frame;
					}
				++timing_index;
				}
			else // It's not time to freeze
				{
				for( Bin bin = 0; bin < get_num_bins(); ++bin )
					out.set_MF( channel, out_frame, bin, get_MF( channel, in_frame, bin ) );
				++out_frame;
				}
			}
		}

	return out;
	}


//========================================================================
// Combinations
//========================================================================

PV PV::replace_amplitudes( const PV & amp_source, const Function<TF, float> & amount ) const
	{
	if( is_null() ) return PV();

	// Input validation
	if( amp_source.is_null() ) return PV();
	auto amount_sampled = sample_function_over_domain( amount );
	std::for_each( amount_sampled.begin(), amount_sampled.end(), []( float & amount ){ amount = std::clamp( amount, 0.0f, 1.0f ); } );

	PV out( get_format() );
	out.clear_buffer();

	const uint32_t num_channels = std::min( amp_source.get_num_channels(), get_num_channels() );
	const uint32_t num_frames   = std::min( amp_source.get_num_frames(),   get_num_frames()   );
	const uint32_t num_bins     = std::min( amp_source.get_num_bins(),     get_num_bins()     );

	for( Channel channel = 0; channel < num_channels; ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( num_frames ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < num_bins; ++bin )
				{
				const auto currentBin = get_MF( channel, frame, bin );
				const float amount_c = amount_sampled[buffer_access( bin, frame, get_num_bins() )];
				out.set_MF( channel, frame, bin,
					{
					amp_source.get_MF( channel, frame, bin ).m * amount_c + currentBin.m * ( 1.0f - amount_c ),
					currentBin.f
					} );
				}
			} );
	return out;
	}

PV PV::subtract_amplitudes( const PV & amp_source, const Function<TF, float> & amount ) const
	{
	if( is_null() ) return PV();

	// Input validation
	if( amp_source.is_null() ) return PV();
	auto amount_sampled = sample_function_over_domain( amount );

	PV out = copy();

	const uint32_t num_channels = std::min( amp_source.get_num_channels(), get_num_channels() );
	const uint32_t num_frames   = std::min( amp_source.get_num_frames(),   get_num_frames()   );
	const uint32_t num_bins     = std::min( amp_source.get_num_bins(),     get_num_bins()     );

	for( Channel channel = 0; channel < num_channels; ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( num_frames ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < num_bins; ++bin )
				{
				auto & outBin = out.get_MF( channel, frame, bin );
				const float amount_c = amount_sampled[buffer_access( bin, frame, get_num_bins() )];
				outBin.m = std::abs( outBin.m - amp_source.get_MF( channel, frame, bin ).m * amount_c );
				} 
			} );

	return out;
	}


//========================================================================
// Generation
//========================================================================

PV PV::synthesize( 
	Second length, 
	const Function<Second, Frequency> & freq, 
	const Function<std::pair<Second, Harmonic>, Magnitude> & harmonic_weights,
	const Function<Second, Frequency> & harmonic_bandwidth,
	const Function<TF, Frequency> & harmonic_frequency_std_dev
	)
	{
	PV::Format format;
	format.num_bins = 2049;
	format.num_channels = 1;
	format.sample_rate = 48000;
	format.analysis_rate = format.sample_rate / 128;
	format.num_frames = length * format.analysis_rate;
	format.window_size = 2048;

	PV out( format );
	out.clear_buffer();
	const float scale = std::sqrt( out.get_dft_size() ); // Don't think too much about this constant
	const Frequency min_frequency = out.get_height() / out.get_num_bins() / 2.0f;

	// Sample frequency and only allow frequencies above min_frequency
	auto frequency_sampled = freq.sample( 0, out.time_to_frame( length ), out.frame_to_time( 1 ) );
	std::for_each( frequency_sampled.begin(), frequency_sampled.end(), [&]( Frequency & f ){ f = std::max( f, min_frequency ); } );

	// Get number of harmonics per frame
	std::vector<Harmonic> num_harmonics_per_frame( out.get_num_frames() );
	for( Frame frame = 0; frame < num_harmonics_per_frame.size(); ++frame ) 
		num_harmonics_per_frame[frame] = std::floor( out.get_height() / frequency_sampled[frame] );

	// Partial sum the number of harmonics per frame to get index offsets into a sampled copy of harmonic_weights
	// Using partial_sum instead of exclusive_scan here, the offset by one index has to happen somewhere
	std::vector<std::size_t> harmonic_weight_sample_partial_sums( out.get_num_frames() );
	std::partial_sum( num_harmonics_per_frame.begin(), num_harmonics_per_frame.end(), harmonic_weight_sample_partial_sums.begin() );
	auto harmonic_buffer_position = [&]( Frame f, Harmonic h ) -> std::size_t
		{
		if( f == 0 ) return h;
		return harmonic_weight_sample_partial_sums[f-1] + h; 
		};

	// Sample harmonic_weights
	std::vector<Magnitude> harmonic_weight_sampled( harmonic_weight_sample_partial_sums.back() );
	for( Frame frame = 0; frame < out.get_num_frames(); ++frame )
		{
		const Harmonic num_harmonics_on_frame = num_harmonics_per_frame[frame];
		for( Harmonic harmonic = 0; harmonic < num_harmonics_on_frame; ++harmonic )
			harmonic_weight_sampled[harmonic_buffer_position( frame, harmonic )] = harmonic_weights( { out.frame_to_time( frame ), harmonic+1 } );
		}

	auto harmonic_bandwidth_sampled = out.sample_function_over_time_domain( harmonic_bandwidth );
	auto harmonic_frequency_std_dev_sampled = out.sample_function_over_domain( harmonic_frequency_std_dev );

	std::default_random_engine rng( std::time( nullptr ) );

	std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( out.get_num_frames() ), [&]( Frame frame )
		{
		const Frequency base_freq_c = frequency_sampled[frame];
		const Harmonic num_harmonics = num_harmonics_per_frame[frame];
		for( Harmonic harmonic = 0; harmonic < num_harmonics; ++harmonic )
			{
			const Magnitude peak_magnitude = harmonic_weight_sampled[harmonic_buffer_position(frame, harmonic)] * scale;
			const Frequency central_harmonic_frequency = base_freq_c * ( harmonic + 1 );
			const Frequency bandwidth = harmonic_bandwidth_sampled[frame] / 2.0f;

			const Frequency low_frequency = central_harmonic_frequency - bandwidth;
			const Frequency high_frequency = central_harmonic_frequency + bandwidth;
			const Bin low_bin  = std::max( 0, 						(Bin) std::ceil ( out.frequency_to_bin( low_frequency ) ) );
			const Bin high_bin = std::min( out.get_num_bins() - 1, 	(Bin) std::floor( out.frequency_to_bin( high_frequency ) ) );

			for( Bin bin = low_bin; bin <= high_bin; ++bin )
				{
				const float window_position = ( out.bin_to_frequency( bin ) - low_frequency ) / ( high_frequency - low_frequency );
				const Magnitude bin_magnitude = peak_magnitude * Windows::hann( window_position );

				const Frequency harmonic_frequency_std_dev_c = harmonic_frequency_std_dev_sampled[out.get_buffer_pos( 0, frame, bin )];
				const Frequency scattered_frequency = harmonic_frequency_std_dev_c <= 0 ?
					central_harmonic_frequency :
					std::normal_distribution<Frequency>( central_harmonic_frequency, harmonic_frequency_std_dev_c )( rng );

				out.get_MF( 0, frame, bin ) = { bin_magnitude, scattered_frequency };
				}
			}
		} );
	
	return out;
	}

//========================================================================
// Uncategorized
//========================================================================

static PV harmonic_scaler( 
	const PV & me, 
	const Function<std::pair<Second, Harmonic>, float> & series, 
	std::function< Frequency ( Frequency, Harmonic )> harmonic_func, 
	Harmonic num_harmonics 
	)
	{
	PV out( me.get_format() );
	out.clear_buffer();

	std::vector<Frequency> series_sampled( num_harmonics * me.get_num_frames() );
	runtime_execution_policy_handler( series.get_execution_policy(), [&]( auto policy )
		{
		std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( me.get_num_frames() ), [&]( Frame frame )
			{
			for( Harmonic harmonic = 0; harmonic < num_harmonics; ++harmonic )
				series_sampled[harmonic + frame*num_harmonics] = series( std::pair( me.frame_to_time( frame ), harmonic ) );
			} );
		} );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( me.get_num_frames() ), [&]( Frame frame )
			{
			auto series_sampled_for_frame = series_sampled.begin() + frame * num_harmonics;

			for( Bin bin = 0; bin < me.get_num_bins(); ++bin )
				{
				const MF & source = me.get_MF( channel, frame, bin );
				if( source.f <= 1.0f ) continue;
				
				for( Harmonic harmonic = 0; harmonic < num_harmonics; ++harmonic )
					{
					const Frequency harmonic_frequency = harmonic_func( source.f, harmonic+1 );
					const Bin harmonic_bin = me.frequency_to_bin( harmonic_frequency );
					if( harmonic_bin >= me.get_num_bins() ) break;

					MF & dest = out.get_MF( channel, frame, harmonic_bin );
					const Magnitude overwrite_mag = source.m * series_sampled_for_frame[harmonic];
					if( dest.m < overwrite_mag )
						dest = { overwrite_mag, harmonic_frequency };
					}
				}
			} );

	return out;
	}

PV PV::add_octaves( const Function<std::pair<Second, Harmonic>, float> & series ) const
	{
	if( is_null() ) return PV();
	return harmonic_scaler( *this, series, []( Frequency f, Harmonic h ){ return f * std::pow( 2, h ); }, std::ceil( std::log2( get_height() ) ) );
	}

PV PV::add_harmonics( const Function<std::pair<Second, Harmonic>, float> & series ) const
	{
	if( is_null() ) return PV();
	return harmonic_scaler( *this, series, []( Frequency f, Harmonic h ){ return f * ( h + 1 ); }, get_num_bins() );
	}

PV PV::shape( const Function<MF,MF> & shaper, bool use_shift_alignment ) const
	{
	if( is_null() ) return PV();

	PV out( get_format() );
	out.clear_buffer();

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		runtime_execution_policy_handler( shaper.get_execution_policy(), [&]( auto policy ){
		std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < get_num_bins(); ++bin )
				{
				const MF inMF = get_MF( channel, frame, bin );
				const MF shapedMF = shaper( inMF );
				
				if( use_shift_alignment )
					{
					const Bin binShift = bin - frequency_to_bin( inMF.f );
					const Bin shapedMFBin = frequency_to_bin( shapedMF.f ) + binShift;
					if( shapedMFBin < 0 || get_num_bins() <= shapedMFBin ) 
						continue;

					MF & outMF = out.get_MF( channel, frame, shapedMFBin );
					if( shapedMF.m > outMF.m )
						outMF = shapedMF;
					}
				else 
					{
					out.get_MF( channel, frame, bin ) = shapedMF;
					}
				}
			} ); } );
		}
			
	return out;
	}

// PV PV::perturb( 
// 	const Function<TF, MF> & mf_std_dev,
// 	float damping 
// 	) const
// 	{
// 	/*
// 	Velocity base with no const axis is most promising so far.
// 	Currently being computed with a summation over both x and y of the random sampling for smoothing.
// 	Other smoothing should give other interesting results.
// 	Magnitude is currently unimplimented, test other overall process options like simplex and then add mag handling.
// 	*/

// 	if( is_null() ) return PV();

// 	const float epsilon = 0.00001;

// 	// Input validation
// 	auto mf_std_dev_sampled = sample_function_over_domain( mf_std_dev );
// 	std::ranges::for_each( mf_std_dev_sampled, []( MF & mf )
// 		{ 
// 		mf.m = std::max( mf.m, 0.0f ); 
// 		mf.f = std::max( mf.f, 0.0f ); 
// 		} );

// 	// Initialize random engine
// 	std::default_random_engine rng( std::time( nullptr ) );

// 	// Sample distributions using random engine. This is done outside the processing loop because random engines aren't thread safe.
// 	std::vector<float> frqAccels( get_num_bins() * get_num_frames() );
// 	for( Bin i : iota_view( 0u, frqAccels.size() ) )
// 		{
// 		const float frequency_std_dev_c = frequency_std_dev_sampled[i];
// 		frqAccels[i] = frequency_std_dev_c < epsilon ? 0 : std::normal_distribution<float>( 0, frequency_std_dev_c / 20 )( rng );
// 		}

// 	std::vector<float> frqVelocs( frqAccels.size() );
// 	for( Bin bin = 0; bin < get_num_bins(); ++bin )
// 		{
// 		int prevPos = buffer_access( bin, 0, get_num_bins() );
// 		frqVelocs[prevPos] = frqAccels[prevPos];

// 		for( Frame frame = 0; frame < get_num_frames(); ++frame )
// 			{
// 			const int buffer_pos = buffer_access( bin, frame, get_num_bins() );
// 			frqVelocs[buffer_pos] = ( frqAccels[buffer_pos] + frqVelocs[prevPos] ) * damping;
// 			prevPos = buffer_pos;
// 			}
// 		}

// 	std::vector<float> frqOffsets( frqAccels.size() );
// 	for( Frame frame = 0; frame < get_num_frames(); ++frame )
// 		{
// 		int prevPos = buffer_access( 0, frame, get_num_bins() );
// 		frqOffsets[prevPos] = frqVelocs[prevPos];

// 		for( Bin bin = 0; bin < get_num_bins(); ++bin )
// 			{
// 			const int buffer_pos = buffer_access( bin, frame, get_num_bins() );
// 			frqOffsets[buffer_pos] = ( frqVelocs[buffer_pos] + frqOffsets[prevPos] ) * damping;
// 			prevPos = buffer_pos;
// 			}
// 		}

// 	PV out( get_format() );

// 	for( Channel channel : iota_view( 0, get_num_channels() ) ) {
// 		float magOffset  = 0;
// 		float freqOffset = 0;
// 		for( Frame frame : iota_view( 0, get_num_frames() ) )	{
// 			const auto buffer_pos = buffer_access( 0, frame, get_num_bins() ); 

// 			const float magnitude_std_dev_c = magnitude_std_dev_sampled[buffer_pos];
// 			const float frequency_std_dev_c = frequency_std_dev_sampled[buffer_pos];

// 			magOffset += magnitude_std_dev_c < epsilon ? 0 : std::normal_distribution<float>( 0, magnitude_std_dev_c / 20 )( rng );
// 			//frqOffset += frequency_std_dev_c < epsilon ? 0 : std::normal_distribution<float>( 0, frequency_std_dev_c / 20 )( rng );

// 			for( Bin bin = 0; bin < get_num_bins(); ++bin ) {
// 				const MF & currentBin = get_MF( channel, frame, bin );
				
				
// 				const Magnitude mag = currentBin.m + magOffset;
// 				const Frequency freq = currentBin.f + frqOffsets[bin] * 200;
// 				//const Frequency freq = currentBin.f * exp2( frqOffsets[frame] );
// 				out.get_MF( channel, frame, bin ) = { mag, freq };
// 				}
// 			}
// 		}

// 	return out;
// 	}

PV predicateNLoudestPartials( const PV & me, const Function<Second, Bin> & num_bins, std::function< bool (Bin, Bin) > predicate )
	{
	// Input validation
	auto num_binsSamples = num_bins.sample( 0, me.get_num_frames(), me.frame_to_time( 1 ) );
	std::ranges::for_each( num_binsSamples, [&]( Bin & b ){ b = std::clamp( b, 0, me.get_num_frames() ); } );

	PV out( me.get_format() );
	out.clear_buffer();

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( me.get_num_frames() ), [&]( Frame frame )
			{
			// Create a buffer storing the frame's magnitudes and which bin that magnitude came from
			using BinMag = std::pair<Bin, Magnitude>;
			std::vector<BinMag> indexAndVolumes( me.get_num_bins() );
			for( Bin bin = 0; bin < me.get_num_bins(); ++bin )
				indexAndVolumes[bin] = { bin, me.get_MF( channel, frame, bin ).m };

			// Sort that buffer by magnitudes
			std::sort( FLAN_PAR_UNSEQ indexAndVolumes.begin(), indexAndVolumes.end(), []( BinMag & a, BinMag & b )
				{ return abs(a.second) > abs(b.second); } );
				
			// The sorted buffer now contains bin indices in loudness order, so go back through and use the predicate 
			// to decide which to copy to the output.
			for( Bin bin = 0; bin < me.get_num_bins(); ++bin )
				{
				const Bin actualBin = indexAndVolumes[bin].first;
				const MF currentMF = me.get_MF( channel, frame, actualBin );
				if( predicate( bin, num_binsSamples[frame] ) )
					out.set_MF( channel, frame, actualBin, currentMF );
				else
					out.set_MF( channel, frame, actualBin, { 0.0, currentMF.f } );
				}
			} );

	return out;
	}

PV PV::retain_n_loudest_partials( const Function<Second, Bin> & num_bins ) const
	{
	if( is_null() ) return PV();
	return predicateNLoudestPartials( *this, num_bins, []( Bin a, Bin b ){ return a < b; } );
	}

PV PV::remove_n_loudest_partials( const Function<Second, Bin> & num_bins ) const
	{
	if( is_null() ) return PV();
	return predicateNLoudestPartials( *this, num_bins, []( Bin a, Bin b ){ return a >= b; } );
	}

PV PV::resonate( Second length, const Function<TF, float> & decay ) const
	{
	if( is_null() ) return PV();

	// Input validation
	if( length < 0 ) 
		length = 0;

	auto format = get_format();
	format.num_frames = get_num_frames() + ceil( time_to_frame( length )  );
	PV out( format );

	auto decay_sampled = out.sample_function_over_domain( decay );
	std::ranges::for_each( decay_sampled, []( float & x ){ x = std::clamp( x, 0.0f, 1.0f ); } );

	// Copy first frame into out
	for( Channel channel = 0; channel < out.get_num_channels(); ++channel )
		for( Bin bin = 0; bin < out.get_num_bins(); ++bin )
			{
			out.set_MF( channel, 0, bin, get_MF( channel, 0, bin ) );
			}
	
	const float secondsPerFrame_c = frame_to_time( 1 );

	for( Channel channel = 0; channel < out.get_num_channels(); ++channel )
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( out.get_num_bins() ), [&]( Bin bin )
			{
			for( Frame frame = 1; frame < out.get_num_frames(); ++frame )
				{
				const float decay_t = std::pow( decay_sampled[buffer_access( bin, frame, get_num_bins() )], secondsPerFrame_c );
				const float decayed_amp = out.get_MF( channel, frame - 1, bin ).m * decay_t;
				if( frame < get_num_frames() && get_MF( channel, frame, bin ).m > decayed_amp )
					out.set_MF( channel, frame, bin, get_MF( channel, frame, bin ) );
				else
					out.set_MF( channel, frame, bin, { decayed_amp, out.get_MF( channel, frame - 1, bin ).f } );
				}
			} );

	return out;
	}

} //End namespace flan
