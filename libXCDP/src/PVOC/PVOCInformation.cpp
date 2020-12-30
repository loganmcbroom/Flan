#include "xcdp/PVOC.h"

#include <iostream>
#include <numeric>

#include "xcdp/DSPUtility.h"
#include "xcdp/WindowFunctions.h"
#include "xcdp/vv_iterator.h"


static const float pi = std::acos( -1.0f );

namespace xcdp {

struct Salience
	{
	float & get( Frame f, Bin b ) { return buffer[ f * numBins + b ]; }
	Frame numFrames;
	Bin numBins;
	std::vector<float> buffer;
	};

// Returns normalized pitch salience
Salience getSalience( const PVOC & me, Channel channel, Frequency minFreq, Frequency maxFreq )
	{
	// Suggested optimal parameters from Salamon & Gomez
	const int Nh = 20;
	const float alpha = 0.8f;
	//const float beta = 1.0f; - Unimplemented in code
	const float gamma = 40.0f;
	const float eTestFactor = std::pow( 10, gamma / 20.0f );

	// Precomute powers of alpha
	std::array<float,Nh> alphaPowers{};
	alphaPowers[0] = 1.0f;
	for( int i = 1; i < Nh; ++i ) alphaPowers[i] = alphaPowers[i-1] * alpha;

	// Precompute cosine outputs
	std::array<float, 2 * 10 + 1> gOutputs;
	for( int i = 0; i <= 10; ++i )
		gOutputs[i] = .5f * ( 1.0f + std::cos( float( i ) / 10.0f * pi / 2.0f ) );

	auto B = [minFreq]( Frequency f ) -> Bin
		{
		return std::floor( 10.0f * 12.0f * std::log2( f / minFreq ) );
		};

	Salience salience;
	salience.numBins = B( maxFreq );
	salience.numFrames = me.getNumFrames();
	salience.buffer.resize( salience.numBins * salience.numFrames, 0 );
		
	for( Frame frame = 0; frame < salience.numFrames; ++frame )
		{
		const float a_M = me.getMaxPartialMagnitude( frame, frame + 1 );
		const float aiLimit = a_M / eTestFactor; // Checking a_i > limit is equivalent to checking 10log_20(aM/ai) < gamma

		// Get frame-wise peaks.
		const PVOC::MF * framePtr = me.getMFPointer( channel, frame, 0 );
		auto magPeakFunc = [framePtr]( int i ) -> float { return framePtr[i].m; };
		const std::vector<vec2> peaks = findPeaks( magPeakFunc, me.getNumBins(), 128, true, false ); // The peak limit here is arbitrary

		for( const vec2 & i : peaks )
			{
			if( i.y() < aiLimit ) continue;

			// Instantaneous amplitude correction
			const Bin bin = i.x();
			const float iF = framePtr[bin].f;
			const float binOffset = iF * me.frequencyToBin() - i.x();
			const float kernelFactor = Windows::HannDFT2( binOffset * me.getWindowSize() / me.getDFTSize() );
			const float iM = kernelFactor >= 0.5f? i.y() / kernelFactor : 0;
				
			for( int h = 0; h < Nh; ++h )
				{
				const Bin B_c = B( iF / ( h + 1 ) );
				if( B_c < 0 ) break;
				// Delta checking is handled in loop bounds
				for( Bin b = std::max( 0, B_c - 10 ); b <= std::min( salience.numBins - 1, B_c + 10 ); ++b )
					salience.get( frame, b ) += gOutputs[std::abs(B_c - b)] * alphaPowers[h] * iM;
				}
			}
		}

	// Normalize
	const float max = *std::max_element( salience.buffer.begin(), salience.buffer.end() );
	std::for_each( salience.buffer.begin(), salience.buffer.end(), [max]( float & x ){ x /= max; } );

	return salience;
	}

std::vector<PVOC::Contour> PVOC::getContours( Channel channel, Frequency minFreq, Frequency maxFreq ) const
	{
	const float TPlus = 0.9f;
	const float TSigma = 0.9f;
	const float pitchBinInCents = 10;
	const float maxdeltaPitch = 80; // cents
	const Frame maxGapLength = .1f * timeToFrame();

	Salience salience = getSalience( *this, channel, minFreq, maxFreq );

	// Get S+ and S-. SPlus is first used to store all peaks, then peaks are moved to S- as needed.
	// Each vec2 contains float representations of pitch bin and amplitude. Intepolation during peak finding means these can take non-integer values.
	std::vector<std::vector<vec2>> SPlus( getNumFrames() );
	std::vector<std::vector<vec2>> SMinus( getNumFrames() );
	
	// Frame-wise peak finding and thresholding
	for( Frame frame = 0; frame < salience.numFrames; ++frame )
		{
		const auto salienceBegin = salience.buffer.begin() + frame * salience.numBins;
		SPlus[frame] = findPeaks( [salienceBegin]( int i ){ return salienceBegin[i]; }, salience.numBins, -1, true, true );
		const float threshold = TPlus * *std::max_element( salienceBegin, salienceBegin + salience.numBins );

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
	const float mean = magSum / numPeaks;
	const float diffSum = std::accumulate( vv_iterator<vec2>::begin( SPlus ), vv_iterator<vec2>::end( SPlus ), 0.0f, [mean]( float s, const vec2 & v )
		{ 
		const float d = v.y() - mean;
		return s + d * d;
		} );
    const float sigma = std::sqrt( diffSum / numPeaks );

	// Filter peaks into S- that are too far below global mean
	const float globalThreshold = mean - TSigma * sigma;
	for( int f = 0; f < getNumFrames(); ++f ) 
		for( int p = SPlus[f].size() - 1; p >= 0 && SPlus[f][p].y() < globalThreshold; --p )
			{
			SMinus[f].push_back( SPlus[f].back() );
			SPlus[f].pop_back();
			}

	// Countour generation
	std::vector<Contour> contours;
	while( true )
		{
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

		auto continuityExtender = [&SPlus, &SMinus, maxSPlusFrame, maxdeltaPitch, pitchBinInCents, maxGapLength, &contour]( Frame start, Frame end )
			{
			const bool forward = end > start;

			float currentPitchBin = contour.bins.back().x();
			Frame currentGap = 0;
			auto continuityFinder = [&currentPitchBin, maxdeltaPitch, pitchBinInCents]( const vec2 & v ){ return std::abs( v.x() - currentPitchBin ) < maxdeltaPitch / pitchBinInCents; };
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
		contour.startFrame = maxSPlusFrame + 1 - contour.bins.size();

		// Things are loaded backwards into contour, flip.
		std::reverse( contour.bins.begin(), contour.bins.end() );

		// Repeat search but left to right
		continuityExtender( maxSPlusFrame + 1, SPlus.size() );
		};

	// Contour info generation
	for( auto & c : contours )
		{
		auto pitchAccess = [&c]( int i ){ return c.bins[i].x(); };
		const vec2 pitchMSD = meanAndStandardDeviation( pitchAccess, c.bins.size() );
		c.pitchMean = pitchMSD.x();
		c.pitchSD = pitchMSD.y();

		auto gainAccess = [&c]( int i ){ return c.bins[i].y(); };
		const vec2 gainMSD = meanAndStandardDeviation( gainAccess, c.bins.size() );
		c.salienceMean = gainMSD.x();
		c.salienceSD = gainMSD.y();
		}
	
	return contours;
	}

}