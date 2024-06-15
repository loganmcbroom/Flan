#include "flan/PV/PV.h"

#include <iostream>
#include <numeric>

#include "flan/DSPUtility.h"
#include "flan/WindowFunctions.h"
#include "flan/Utility/vv_iterator.h"
#include "flan/Utility/iota_iter.h"
#include "flan/Utility/execution.h"

static const float pi = std::acos( -1.0f );

namespace flan {

// Check if a and b are within a half note of one another
static bool notesClose( Frequency a, Frequency b )
	{
	static const float lo = std::pow( 2.0f, -1.0f / 24.0f );
	static const float hi = std::pow( 2.0f,  1.0f / 24.0f );

	if( b < 0.01 ) return false; // We have no interest in notes that low
	const float r = a / b;
	return lo < r && r < hi; 
	}

// Returns normalized pitch salience
PV::Salience PV::get_salience( Channel channel, Frequency min_frequency, Frequency max_frequency ) const
	{
	auto hann_dft2 = []( float f ) -> float
		{
		if( f == 0 ) return 1.0f;
		if( std::abs( f ) == 1.0f ) return 0.5f;
		return std::sin( pi * f ) / ( pi * f * ( 1.0f - f * f ) );
		};

	if( is_null() )
		return Salience();

	// Suggested optimal parameters from Salamon & Gomez
	const int binEffectDist = 10;
	const int Nh = 2 * 10;
	const float alpha = 0.8f;
	//const float beta = 1.0f; - Unimplemented in code
	const float gamma = 40.0f;
	const float eTestFactor = std::pow( 10, gamma / 20.0f );

	// Precomute powers of alpha
	std::array<float,Nh> alphaPowers{};
	alphaPowers[0] = 1.0f;
	for( int i = 1; i < Nh; ++i ) alphaPowers[i] = alphaPowers[i-1] * alpha;

	// Precompute cosine outputs
	std::array<float, 2 * binEffectDist + 1> gOutputs;
	std::for_each( FLAN_PAR_UNSEQ iota_iter(0), iota_iter( binEffectDist+1 ), [&]( int i )
		{ gOutputs[i] = .5f * ( 1.0f + std::cos( float( i ) / binEffectDist * pi / 2.0f ) ); } );

	auto B = [log2MinFreq = std::log2( min_frequency )]( Frequency f ) -> Bin
		{
		return std::round( 10.0f * 12.0f * ( std::log2( f ) - log2MinFreq ) );
		};

	Salience salience;
	salience.num_bins = B( max_frequency );
	salience.num_frames = get_num_frames();
	salience.buffer.resize( salience.num_bins * salience.num_frames, 0 );
		
	std::for_each( FLAN_PAR_UNSEQ iota_iter(0), iota_iter(salience.num_frames), [&]( Frame frame )
		{
		const Magnitude a_M = get_max_partial_magnitude( frame, frame + 1 ); // Get mazimum magnitude for this frame
		const float aiLimit = a_M / eTestFactor; // Checking a_i > limit is equivalent to checking 10log_20(aM/ai) < gamma

		// Get frame-wise peaks
		const MF * framePtr = get_MF_pointer( channel, frame, 0 );
		auto magPeakFunc = [framePtr]( int i ) -> float { return framePtr[i].m; };
		const std::vector<vec2> peaks = find_peaks( magPeakFunc, get_num_bins(), -1, false, false );

		for( const vec2 & i : peaks )
			{
			if( i.y() < aiLimit ) continue;

			// Instantaneous amplitude correction
			const Bin bin = i.x(); // PV bin
			const float iF = framePtr[bin].f;
			const float binOffset = frequency_to_bin( iF ) - i.x();
			const float kernelFactor = hann_dft2( binOffset * get_window_size() / get_dft_size() );
			const float iM = kernelFactor >= 0.5f? i.y() / kernelFactor : 0;
				
			for( int h = 0; h < Nh; ++h ) // For each subharmonic the current peak could be a harmonic of:
				{
				const Bin B_c = B( iF / ( h + 1 ) ); // Get subharmonic pitch bin
				if( B_c < 0 ) break;

				// Contribute to the salience of the subharmonic for all nearby bins
				// Delta checking is handled in loop bounds
				const Bin loopStart = std::max( 0, B_c - binEffectDist );
				const Bin loopEnd = std::min( salience.num_bins - 1, B_c + binEffectDist );
				for( Bin b = loopStart; b <= loopEnd; ++b )
					salience.get( frame, b ) += gOutputs[std::abs(B_c - b)] * alphaPowers[h] * iM;
				}
			}
		} );

	// Normalize
	const float max = *std::max_element( salience.buffer.begin(), salience.buffer.end() );
	std::for_each( salience.buffer.begin(), salience.buffer.end(), [max]( float & x ){ x /= max; } );

	return salience;
	}
	
std::vector<PV::Contour> PV::get_contours( Channel channel, Frequency min_frequency, Frequency max_frequency, 
	Frame filter_short, float filter_quiet, flan_CANCEL_ARG_CPP ) const
	{
	const float TPlus = 0.9f;
	const float TSigma = 0.9f;
	const float pitchBinInCents = 10;
	const float maxdeltaPitch = 80; // cents
	const Frame maxGapLength = time_to_frame( .1f );

	const Salience salience = get_salience( channel, min_frequency, max_frequency );
	if( salience.buffer.empty() ) return std::vector<PV::Contour>();

	// Get S+ and S-. SPlus is first used to store all peaks, then peaks are moved to S- as needed.
	// Each vec2 contains float representations of pitch bin and amplitude. Intepolation during peak finding means these can take non-integer values.
	std::vector<std::vector<vec2>> SPlus( get_num_frames() );
	std::vector<std::vector<vec2>> SMinus( get_num_frames() );
	
	// Frame-wise peak finding and thresholding
	for( Frame frame = 0; frame < salience.num_frames; ++frame )
		{
		flan_CANCEL_POINT( std::vector<PV::Contour>() );

		const auto salienceBegin = salience.buffer.begin() + frame * salience.num_bins;
		SPlus[frame] = find_peaks( [salienceBegin]( int i ){ return salienceBegin[i]; }, salience.num_bins, -1, true, true );
		const float threshold = TPlus * *std::max_element( salienceBegin, salienceBegin + salience.num_bins );

		// Filter peaks below threshhold into S-
		while( ! SPlus[frame].empty() && SPlus[frame].back().y() < threshold )
			{
			SMinus[frame].push_back( SPlus[frame].back() );
			SPlus[frame].pop_back();
			}
		}

	// Now all peaks above frame-wise threshhold are in SPlus, find the mean and SD of those
	const int numPeaks = std::accumulate( SPlus.begin(), SPlus.end(), 0, []( float s, const std::vector<vec2> & v ){ return s + v.size(); } );
	const float magSum = std::accumulate( vv_iterator<vec2>::begin( SPlus ), vv_iterator<vec2>::end( SPlus ), 0.0f, []( float s, const vec2 & v ){ return s + v.y(); } );
	if( numPeaks == 0 ) return std::vector<PV::Contour>();
	const float mean = magSum / numPeaks;
	const float diffSum = std::accumulate( vv_iterator<vec2>::begin( SPlus ), vv_iterator<vec2>::end( SPlus ), 0.0f, [mean]( float s, const vec2 & v )
		{ 
		const float d = v.y() - mean;
		return s + d * d;
		} );
    const float sigma = std::sqrt( diffSum / numPeaks );

	// Filter peaks into S- that are too far below global mean
	const float globalThreshold = mean - TSigma * sigma;
	for( int f = 0; f < get_num_frames(); ++f ) 
		{
		for( int p = SPlus[f].size() - 1; p >= 0 && SPlus[f][p].y() < globalThreshold; --p )
			{
			SMinus[f].push_back( SPlus[f].back() );
			SPlus[f].pop_back();
			}
		}

	// Countour generation
	std::vector<Contour> contours;
	while( true )
		{
		flan_CANCEL_POINT( std::vector<PV::Contour>() );

		// Find current SPlus max. We only need to check the first element of each frame, as each frame is sorted by magnitude.
		const int maxSPlusFrame = std::distance( SPlus.begin(), std::max_element( SPlus.begin(), SPlus.end(), []( const std::vector<vec2> & v, const std::vector<vec2> & w )
			{ 
			return ( v.empty()? 0 : v[0].y() )
				 < ( w.empty()? 0 : w[0].y() ); 
			} ) );
		if( maxSPlusFrame >= SPlus.size() || SPlus[maxSPlusFrame].empty() ) break;

		// Push back new contour
		contours.emplace_back();
		Contour & contour = contours.back();
		contour.bins.push_back( SPlus[maxSPlusFrame][0] );
		SPlus[maxSPlusFrame].erase( SPlus[maxSPlusFrame].begin() );

		auto continuityExtender = [&]( Frame start, Frame end )
			{
			const bool forward = end > start;

			float currentPitchBin = contour.bins.back().x();
			Frame currentGap = 0;
			auto continuityFinder = [&currentPitchBin, maxdeltaPitch, pitchBinInCents]( const vec2 & v )
				{ 
				return std::abs( v.x() - currentPitchBin ) < maxdeltaPitch / pitchBinInCents; 
				};

			for( Frame frame = start; frame != end && currentGap < maxGapLength; forward? ++frame : --frame )
				{
				// Look for peak near current peak pitch on next frame in SPlus
				const auto newSPlusPitchBin = std::find_if( SPlus[frame].begin(), SPlus[frame].end(), continuityFinder );
				if( newSPlusPitchBin != SPlus[frame].end() ) // Found in SPlus, sweet. Update and continue.
					{
					contour.bins.push_back( *newSPlusPitchBin );
					currentPitchBin = contour.bins.back().x();
					SPlus[frame].erase( newSPlusPitchBin );
					currentGap = 0;
					}
				else // None found! Search S-.
					{
					const auto newSMinusPitchBin = std::find_if( SMinus[frame].begin(), SMinus[frame].end(), continuityFinder );
					if( newSMinusPitchBin != SMinus[frame].end() ) // Backup found, but update time limit until we find in SPlus.
						{
						contour.bins.push_back( *newSMinusPitchBin );
						currentPitchBin = contour.bins.back().x();
						SMinus[frame].erase( newSMinusPitchBin );
						++currentGap;
						}
					else break; // Aint nothin here for us, bail.
					}
				}
			};

		// Read right to left until the salience peaks run out
		continuityExtender( maxSPlusFrame - 1, -1 );

		// Set contour start bin
		contour.start_frame = maxSPlusFrame + 1 - contour.bins.size();

		// Things are loaded backwards into contour, flip.
		std::reverse( contour.bins.begin(), contour.bins.end() );

		// Repeat search but left to right
		continuityExtender( maxSPlusFrame + 1, SPlus.size() );

		// Filter by length
		if( contours.back().bins.size() < filter_short )
			contours.pop_back();

		// Contour info generation
		auto gainAccess = [&contour]( int i ){ return contour.bins[i].y(); };
		const vec2 gainMSD = mean_and_sd( gainAccess, contour.bins.size() );
		contour.salience_mean = gainMSD.x();
		contour.salience_std_dev = gainMSD.y();

		auto pitchAccess = [&contour]( int i ){ return contour.bins[i].x(); };
		const vec2 pitchMSD = mean_and_sd( pitchAccess, contour.bins.size() );
		contour.pitch_mean = pitchMSD.x();
		contour.pitch_std_dev = pitchMSD.y();
		};

	// Filter by salience
	std::vector<Contour> contoursFiltered;
	const float maxMeanSalience = std::max_element( contours.begin(), contours.end(), []( const Contour & a, const Contour & b )
		{ return a.salience_mean < b.salience_mean; } )->salience_mean;
	const float minSalience = maxMeanSalience / filter_quiet;
	for( auto & c : contours )
		if( c.salience_mean >= minSalience )
			contoursFiltered.emplace_back( std::move( c ) );
	
	return contoursFiltered;
	}

PV PV::prism( const PrismFunc & prism_func, bool use_local_contour_time, flan_CANCEL_ARG_CPP ) const
	{
	if( is_null() ) return PV();

	const Frequency min_frequency = 55.0;
	const Frequency max_frequency = 1760.0;

	PV out = copy();

	auto pitch_bin_to_freq = [min_frequency]( float bin ){ return min_frequency * std::pow( 2.0f, bin / 120.0f ); };

	// The pitch bins from the contour aren't exact enough, so we need to more accurately estimate the base frequency with a weighted mean of nearby bins
	auto get_accurate_freq = [this]( Channel channel, Frame frame, Frequency approx_freq )
		{
		float total_weighted_freq = 0;
		float total_mag = 0;
		for( Bin bin = 0; bin < get_num_bins(); ++bin )
			{
			const MF & source_MF = get_MF( channel, frame, bin );
			if( source_MF.f <= 0 ) continue;
			if( notesClose( source_MF.f, approx_freq ) ) // Assert the bin has a frequency within a half a note of the target harmonic
				{
				const Magnitude abs_mag = std::abs( source_MF.m );
				total_weighted_freq += source_MF.f * abs_mag;
				total_mag += abs_mag;
				}
			}
		if( total_mag == 0 ) return 0.0f;
		return total_weighted_freq / total_mag;
		};

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		std::vector<Contour> contours = get_contours( channel, min_frequency, max_frequency, 60, 20, canceller );
		if( contours.empty() ) return PV();
		std::sort( FLAN_PAR_UNSEQ contours.begin(), contours.end(), []( const Contour & a, const Contour & b ){ return a.start_frame < b.start_frame; } );
		for( Index contour_index = 0; contour_index < contours.size(); ++contour_index )
			{
			const Contour & contour = contours[contour_index];

			flan::for_each_i( contour.bins.size(), prism_func.get_execution_policy(), [&]( Frame contour_frame ) {
				const Frame frame = contour_frame + contour.start_frame;
				const MF * const source_frame_ptr = get_MF_pointer( channel, frame, 0 );

				// The pitch bin from the contour isn't exact enough, so we need to more accurately estimate the base frequency with a weighted mean of nearby bins
				const Frequency approx_base_freq = pitch_bin_to_freq( contour.bins[contour_frame].x() );
				const Frequency base_freq = get_accurate_freq( channel, frame, approx_base_freq );
				if( base_freq < 1.0f )
					return;
				const size_t max_num_harmonics = std::floor( get_height() / base_freq );

				// First we will locate everything that needs to be changed and clear it from out
				// We need to find all bins before changing them because we don't want to overwrite bins that need to be changed
				std::vector<std::vector<Bin>> bins_to_change_vec( max_num_harmonics ); // Harmonic -> Bin vector
				for( auto & v : bins_to_change_vec )
					v.reserve( 20 );

				for( Harmonic harmonic = 0; harmonic < bins_to_change_vec.size(); ++harmonic ) // For each harmonic
					{
					const Frequency freq = base_freq * ( harmonic + 1 );

					// For bins near the harmonic, we find all bins with a freq close to the harmonic
					std::vector<Bin> & bins_to_change = bins_to_change_vec[harmonic];
					const Bin start_bin = bound_bin( frequency_to_bin( freq ) - 10 );
					const Bin end_bin   = bound_bin( frequency_to_bin( freq ) + 10 );
					const MF * source_ptr = get_MF_pointer( channel, frame, start_bin );
					for( Bin bin = start_bin; bin <= end_bin; ++bin, ++source_ptr )
						{
						const MF & source_MF = *source_ptr;
						if( source_MF.f <= 0 ) continue;
						// Assert the bin has a frequency within a half a note of the target harmonic
						if( notesClose( source_MF.f, freq ) ) 
							bins_to_change.push_back( bin );
						}

					// Clear source bins from out
					MF * const clear_frame_ptr = out.get_MF_pointer( channel, frame, 0 );
					for( auto bin : bins_to_change )
						( clear_frame_ptr + bin )->m = 0;
					}

				// For each harmonic, find the dominant bin
				std::vector<Bin> harmonic_max_mag_bins( max_num_harmonics );
				std::vector<Magnitude> harmonic_max_mags( max_num_harmonics, 0 );
				for( Harmonic harmonic = 0; harmonic < harmonic_max_mag_bins.size(); ++harmonic )
					{
					const std::vector<Bin> & bins_to_change = bins_to_change_vec[harmonic];
					if( bins_to_change.empty() ) continue;
					const Bin harmonic_max_mag_bin = *std::max_element( bins_to_change.begin(), bins_to_change.end(), [source_frame_ptr]( Bin a, Bin b )
						{ return ( source_frame_ptr + a )->m < ( source_frame_ptr + b )->m; } );
					harmonic_max_mag_bins[harmonic] = harmonic_max_mag_bin;
					harmonic_max_mags[harmonic] = ( source_frame_ptr + harmonic_max_mag_bin )->m;
					if( harmonic_max_mags[harmonic] < 0.01 ) // This is an arbitrary delta that needs improvement
						{
						harmonic_max_mags[harmonic] = 0;
						}
					}

				// Now we go back through everything in the frame that needs writing and write it
				for( Harmonic harmonic = 0; harmonic < harmonic_max_mag_bins.size(); ++harmonic )
					{
					const Frequency freq = base_freq * ( harmonic + 1 );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

					const MF modified_harmonic = prism_func( 
						contour_index, 
						frame_to_time( use_local_contour_time ? contour_frame : frame ), 
						harmonic + 1, 
						base_freq, 
						harmonic_max_mags );
						
					if( modified_harmonic.f < 0 ) 
						continue;
					
					if( harmonic_max_mags[harmonic] != 0 )
						{
						const std::vector<Bin> & bins_to_change = bins_to_change_vec[harmonic];

						const Bin new_max_mag_bin = modified_harmonic.f / freq * harmonic_max_mag_bins[harmonic];
						const Bin bin_shift = new_max_mag_bin - harmonic_max_mag_bins[harmonic];
						const Frequency f_scale = modified_harmonic.f / freq;
						const Frequency m_scale = modified_harmonic.m / harmonic_max_mags[harmonic];
						for( auto bin : bins_to_change )
							{
							const Bin new_bin = bin + bin_shift;
							if( new_bin < 0 || new_bin >= get_num_bins() )	
								continue;

							// Set new bin to updated value
							//std::lock_guard<std::mutex> lock( mutexBuffer[frame] );
							const MF & source_MF = *( source_frame_ptr + bin );
							MF & dest_MF = out.get_MF( channel, frame, new_bin );
							if( dest_MF.m < source_MF.m * m_scale )
								dest_MF = { source_MF.m * m_scale, source_MF.f * f_scale };
							}
						}
					else // There is no harmonic to scale, so create one
						{
						const Frequency frequency_bandwidth = 10;
						const Frequency low_frequency = modified_harmonic.f - frequency_bandwidth / 2;
						const Frequency high_frequency = modified_harmonic.f + frequency_bandwidth / 2;
						const Bin low_bin  = std::max( 0, 						(Bin) std::ceil ( out.frequency_to_bin( low_frequency ) ) );
						const Bin high_bin = std::min( out.get_num_bins() - 1, 	(Bin) std::floor( out.frequency_to_bin( high_frequency ) ) );

						for( Bin bin = low_bin; bin <= high_bin; ++bin )
							{
							const float window_position = ( out.bin_to_frequency( bin ) - low_frequency ) / frequency_bandwidth;
							const Magnitude bin_magnitude = modified_harmonic.m * Windows::hann( window_position );
							out.get_MF( channel, frame, bin ) = { bin_magnitude, modified_harmonic.f };
							}
						}
					}
				});
			}
		}

	return out;
	}

}