#include "flan/PVOC.h"

#include <iostream>
#include <numeric>

#include "flan/DSPUtility.h"
#include "flan/WindowFunctions.h"
#include "flan/vv_iterator.h"
#include "flan/CLContext.h"
#include "flan/CLProgs.h"

static const float pi = std::acos( -1.0f );

namespace flan {

// Returns normalized pitch salience
PVOC::Salience PVOC::getSalience( Channel channel, Frequency minFreq, Frequency maxFreq, flan_CANCEL_ARG_CPP ) const
	{
	if( isNull() )
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
	for( int i = 0; i <= binEffectDist; ++i )
		gOutputs[i] = .5f * ( 1.0f + std::cos( float( i ) / binEffectDist * pi / 2.0f ) );

	auto B = [minFreq]( Frequency f ) -> Bin
		{
		return std::round( 10.0f * 12.0f * std::log2( f / minFreq ) );
		};

	Salience salience;
	salience.numBins = B( maxFreq );
	salience.numFrames = getNumFrames();
	salience.buffer.resize( salience.numBins * salience.numFrames, 0 );
		
	for( Frame frame = 0; frame < salience.numFrames; ++frame )
		{
		const float a_M = getMaxPartialMagnitude( frame, frame + 1 );
		const float aiLimit = a_M / eTestFactor; // Checking a_i > limit is equivalent to checking 10log_20(aM/ai) < gamma

		// Get frame-wise peaks.
		const PVOC::MF * framePtr = getMFPointer( channel, frame, 0 );
		auto magPeakFunc = [framePtr]( int i ) -> float { return framePtr[i].m; };
		const std::vector<vec2> peaks = findPeaks( magPeakFunc, getNumBins(), -1, false, false ); // The peak limit here is arbitrary

		for( const vec2 & i : peaks )
			{
			if( i.y() < aiLimit ) continue;

			// Instantaneous amplitude correction
			const Bin bin = i.x(); // PVOC bin
			const float iF = framePtr[bin].f;
			const float binOffset = iF * frequencyToBin() - i.x();
			const float kernelFactor = Windows::HannDFT2( binOffset * getWindowSize() / getDFTSize() );
			const float iM = kernelFactor >= 0.5f? i.y() / kernelFactor : 0;
				
			for( int h = 0; h < Nh; ++h ) // For each subharmonic the current peak could be a harmonic of:
				{
				flan_CANCEL_POINT( Salience() );

				const Bin B_c = B( iF / ( h + 1 ) ); // Get subharmonic pitch bin
				if( B_c < 0 ) break;

				// Contribute to the salience of the subharmonic for all nearby bins
				// Delta checking is handled in loop bounds
				for( Bin b = std::max( 0, B_c - binEffectDist ); b <= std::min( salience.numBins - 1, B_c + binEffectDist ); ++b )
					salience.get( frame, b ) += gOutputs[std::abs(B_c - b)] * alphaPowers[h] * iM;
				}
			}
		}

	// Normalize
	const float max = *std::max_element( salience.buffer.begin(), salience.buffer.end() );
	std::for_each( salience.buffer.begin(), salience.buffer.end(), [max]( float & x ){ x /= max; } );

	return salience;
	}
	
PVOC::Salience PVOC::getSalience_gpu( Channel channel, Frequency minFreq, Frequency maxFreq, flan_CANCEL_ARG_CPP ) const
	{
	if( isNull() )
		return Salience();

	// Suggested optimal parameters from Salamon & Gomez
	const int binEffectDist = 10;
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
	std::array<float, 2 * binEffectDist + 1> gOutputs;
	for( int i = 0; i <= binEffectDist; ++i )
		gOutputs[i] = .5f * ( 1.0f + std::cos( float( i ) / binEffectDist * pi / 2.0f ) );

	auto B = [minFreq]( Frequency f ) -> Bin
		{
		return std::round( 10.0f * 12.0f * std::log2( f / minFreq ) );
		};

	Salience salience;
	salience.numBins = B( maxFreq );
	salience.numFrames = getNumFrames();
	salience.buffer.resize( salience.numBins * salience.numFrames, 0 );

	// Prepare OpenCL
	CLContext cl;
	ProgramHelper programHelper( cl, CLProgs::PVOC_getSalience );
	
	// Create cl buffers
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * getNumBins() * getNumFrames() ); // Input for entire channel
	cl::Buffer clOut( cl.context, CL_MEM_READ_WRITE, sizeof( float ) * salience.buffer.size() ); // Get back Salience
	cl::Buffer clAlphas( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * alphaPowers.size() );
	cl::Buffer clCosines( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * gOutputs.size() );

	cl::Kernel cl_getSalience( programHelper.program, "getSalience" );
	cl_getSalience.setArg( 0, clIn );
	cl_getSalience.setArg( 1, clOut );
	cl_getSalience.setArg( 2, clAlphas );
	cl_getSalience.setArg( 3, clCosines );
	cl_getSalience.setArg( 4, (int) getNumBins() );
	cl_getSalience.setArg( 5, (int) salience.numBins );
	cl_getSalience.setArg( 6, (float) minFreq );
	cl_getSalience.setArg( 7, (float) frequencyToBin() );
	cl_getSalience.setArg( 8, (float) getWindowSize() / getDFTSize() );
	flan_CANCEL_POINT( Salience() );
		
	// Copy data to cl buffer
	cl::copy( cl.queue, getMFPointer( channel, 0, 0 ), getMFPointer( channel, 0, 0 ) + getNumBins() * getNumFrames(), clIn );
	cl::copy( cl.queue, alphaPowers.begin(), alphaPowers.end(), clAlphas );
	cl::copy( cl.queue, gOutputs.begin(), gOutputs.end(), clCosines );
	flan_CANCEL_POINT( Salience() );
	cl.queue.enqueueBarrierWithWaitList();

	// Call gpu function and wait for it to finish
	// Range args are global offset, global work size, and local work size
	cl.queue.enqueueNDRangeKernel( cl_getSalience, cl::NullRange, cl::NDRange( getNumFrames() ), cl::NullRange );
	flan_CANCEL_POINT( Salience() );
	cl.queue.enqueueBarrierWithWaitList();

	// Copy data back to host and wait for copy
	cl::copy( cl.queue, clOut, salience.buffer.begin(), salience.buffer.end() );
	cl.queue.enqueueBarrierWithWaitList();

	// Normalize
	const float max = *std::max_element( salience.buffer.begin(), salience.buffer.end() );
	std::for_each( salience.buffer.begin(), salience.buffer.end(), [max]( float & x ){ x /= max; } );

	return salience;
	}

std::vector<PVOC::Contour> PVOC::getContours( Channel channel, Frequency minFreq, Frequency maxFreq, 
	Frame filterShort, float filterQuiet, flan_CANCEL_ARG_CPP ) const
	{
	const float TPlus = 0.9f;
	const float TSigma = 0.9f;
	const float pitchBinInCents = 10;
	const float maxdeltaPitch = 80; // cents
	const Frame maxGapLength = .1f * timeToFrame();

	Salience salience = getSalience_gpu( channel, minFreq, maxFreq, canceller );

	// Get S+ and S-. SPlus is first used to store all peaks, then peaks are moved to S- as needed.
	// Each vec2 contains float representations of pitch bin and amplitude. Intepolation during peak finding means these can take non-integer values.
	std::vector<std::vector<vec2>> SPlus( getNumFrames() );
	std::vector<std::vector<vec2>> SMinus( getNumFrames() );
	
	// Frame-wise peak finding and thresholding
	for( Frame frame = 0; frame < salience.numFrames; ++frame )
		{
		flan_CANCEL_POINT( std::vector<PVOC::Contour>() );
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
		{
		flan_CANCEL_POINT( std::vector<PVOC::Contour>() );
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
		flan_CANCEL_POINT( std::vector<PVOC::Contour>() );

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

		// Filter by length
		if( contours.back().bins.size() < filterShort )
			contours.pop_back();

		// Contour info generation
		auto gainAccess = [&contour]( int i ){ return contour.bins[i].y(); };
		const vec2 gainMSD = meanAndStandardDeviation( gainAccess, contour.bins.size() );
		contour.salienceMean = gainMSD.x();
		contour.salienceSD = gainMSD.y();

		auto pitchAccess = [&contour]( int i ){ return contour.bins[i].x(); };
		const vec2 pitchMSD = meanAndStandardDeviation( pitchAccess, contour.bins.size() );
		contour.pitchMean = pitchMSD.x();
		contour.pitchSD = pitchMSD.y();
		};

	// Filter by salience
	std::vector<Contour> contoursFiltered;
	const float maxMeanSalience = std::max_element( contours.begin(), contours.end(), []( const Contour & a, const Contour & b )
		{ return a.salienceMean < b.salienceMean; } )->salienceMean;
	const float minSalience = maxMeanSalience / filterQuiet;
	for( auto & c : contours )
		if( c.salienceMean >= minSalience )
			contoursFiltered.emplace_back( std::move( c ) );
	
	return contoursFiltered;
	}

PVOC PVOC::prism( PrismFunc f, bool perNote, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	const Frequency minFreq = 55.0;
	const Frequency maxFreq = 1760.0;

	PVOC out = deepCopy();

	auto notesClose = [lo = std::pow( 2.0f, -1.0f / 24.0f ), hi = std::pow( 2.0f, 1.0f / 24.0f )]( Frequency a, Frequency b ) 
		{ 
		const float r = a / b;
		return lo < r && r < hi; 
		};

	auto pitchBinToFreq = [minFreq]( float bin ){ return minFreq * std::pow( 2.0f, bin / 120.0f); };

	auto getAccurateFreq = [this, notesClose]( Channel channel, Frame frame, Frequency approxFreq )
		{
		std::vector<Frequency> approxBaseFreqs;
		float totalMag = 0;
		for( Bin bin = 0; bin < getNumBins(); ++bin )
			{
			const MF & sourceMF = getMF( channel, frame, bin );
			if( sourceMF.f <= 0 ) continue;
			if( notesClose( sourceMF.f, approxFreq ) ) // Assert the bin has a frequency within a half a note of the target harmonic
				{
				const Magnitude absMag = std::abs( sourceMF.m );
				approxBaseFreqs.push_back( sourceMF.f * absMag );
				totalMag += absMag;
				}
			}
		if( totalMag == 0 ) return 0.0f;
		return std::accumulate( approxBaseFreqs.begin(), approxBaseFreqs.end(), 0.0f ) / totalMag;
		};

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		std::vector<Contour> contours = getContours( channel, minFreq, maxFreq, 60, 20, canceller );
		for( auto & c : contours )
			{
			for( Frame contourIndex = 0; contourIndex < c.bins.size(); ++contourIndex )
				{
				const Frame frame = contourIndex + c.startFrame;
				const MF * const sourceFramePtr = getMFPointer( channel, frame, 0 );

				// The pitch bin from the contour isn't exact enough, so we need to more accurately estimate the base frequency with a weighted mean of nearby bins
				const Frequency approxBaseFreq = pitchBinToFreq( c.bins[contourIndex].x() );
				const Frequency baseFreq = getAccurateFreq( channel, frame, approxBaseFreq );
				if( baseFreq < 1.0f )
					continue;

				// First we will locate everything that needs to be changed and clear it from out
				// We need to find all bins before changing them because we could overwrite bins that need to be changed
				std::vector<std::vector<Bin>> binsToChangeVec( std::ceil( getHeight() / baseFreq ) ); // Harmonic -> Bin vector
				for( auto & v : binsToChangeVec )
					v.reserve( 20 );

				for( int harmonic = 0; harmonic < std::floor( getHeight() / baseFreq ); ++harmonic ) // For each harmonic
					{
					Frequency freq = baseFreq * ( harmonic + 1 );
					flan_CANCEL_POINT( PVOC() );

					// For bins near the harmonic, we find all bins with a freq close to the harmonic
					std::vector<Bin> & binsToChange = binsToChangeVec[harmonic];
					const Bin startBin = boundBin( freq * frequencyToBin() - 10 );
					const Bin endBin   = boundBin( freq * frequencyToBin() + 10 );
					const MF * sourcePtr = getMFPointer( channel, frame, startBin );
					for( Bin bin = startBin; bin <= endBin; ++bin, ++sourcePtr )
						{
						const MF & sourceMF = *sourcePtr;
						if( sourceMF.f <= 0 ) continue;
						// Assert the bin has a frequency within a half a note of the target harmonic
						if( notesClose( sourceMF.f, freq ) ) 
							binsToChange.push_back( bin );
						}

					// Clear source bins from out
					MF * const clearFramePtr = out.getMFPointer( channel, frame, 0 );
					for( auto bin : binsToChange )
						( clearFramePtr + bin )->m = 0;
					}

				// For each harmonic, find the dominant bin
				std::vector<Bin> harmonicMaxMagBins( getHeight() / baseFreq );
				std::vector<Magnitude> harmonicMaxMags( harmonicMaxMagBins.size(), 0 );
				for( int harmonic = 0; harmonic < std::floor( getHeight() / baseFreq ); ++harmonic ) // For each harmonic
					{
					Frequency freq = baseFreq * ( harmonic + 1 );
					const std::vector<Bin> & binsToChange = binsToChangeVec[harmonic];
					if( binsToChange.empty() ) continue;
					const Bin harmonicMaxMagBin = *std::max_element( binsToChange.begin(), binsToChange.end(), [sourceFramePtr]( Bin a, Bin b )
						{ return ( sourceFramePtr + a )->m < ( sourceFramePtr + b )->m; } );
					harmonicMaxMagBins[harmonic] = harmonicMaxMagBin;
					harmonicMaxMags[harmonic] = ( sourceFramePtr + harmonicMaxMagBin )->m;
					if( harmonicMaxMags[harmonic] < 0.01 ) // This is an arbitrary delta that needs improvement
						{
						harmonicMaxMags[harmonic] = 0;
						}
					}

				// Now we go back through everything in the frame that needs writing and write it
				for( int harmonic = 0; harmonic < std::floor( getHeight() / baseFreq ); ++harmonic ) // For each harmonic
					{
					Frequency freq = baseFreq * ( harmonic + 1 );
					flan_CANCEL_POINT( PVOC() );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

					const MF modifiedHarmonic = f( ( perNote? contourIndex : frame ) * frameToTime(), harmonic, baseFreq, harmonicMaxMags );
					if( modifiedHarmonic.f < 0 ) 
						continue;
					
					if( harmonicMaxMags[harmonic] != 0 )
						{
						const std::vector<Bin> & binsToChange = binsToChangeVec[harmonic];

						const Bin newMaxMagBin = modifiedHarmonic.f / freq * harmonicMaxMagBins[harmonic];
						const Bin binShift = newMaxMagBin - harmonicMaxMagBins[harmonic];
						const Frequency fScale = modifiedHarmonic.f / freq;
						const Frequency mScale = modifiedHarmonic.m / harmonicMaxMags[harmonic];
						for( auto bin : binsToChange )
							{
							const Bin newBin = bin + binShift;
							if( newBin < 0 || newBin >= getNumBins() )	
								continue;

							// Set new bin to updated value
							const MF & sourceMF = *( sourceFramePtr + bin );
							MF & destMF = out.getMF( channel, frame, newBin );
							if( destMF.m < sourceMF.m * mScale )
								destMF = { sourceMF.m * mScale, sourceMF.f * fScale };
							}
						}
					else // There is no harmonic to scale, so create one
						{
						MF & destMF = out.getMF( channel, frame, modifiedHarmonic.f * frequencyToBin() );
						destMF = modifiedHarmonic;
						}
					}
				}
			}
		}

	return out;
	}

}