#include "flan/PV/PV.h"

#include <iostream>
#include <algorithm>
#include <ranges>

#include "spline/spline.h"
#include "flan/Utility/iota_iter.h"

using namespace std::ranges;

namespace flan {

PV PV::modify( const Func2x2 & mod, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	const Bin numBins = getNumBins();
	const Frame numFrames = getNumFrames();
	const uint32_t inChannelDataCount = numFrames * numBins;

	//Sample mod function and convert output to frame/bin
	auto modSamples = mod.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
	std::for_each( std::execution::par_unseq, modSamples.begin(), modSamples.end(), [&]( vec2 & v ){ 
		v.x() *= timeToFrame();
		v.y() *= frequencyToBin(); 
		} );

	// Get the largest frame value among all mod samples, rounded up
	const float lastOutputFrame = std::ceil( std::ranges::max_element( modSamples, std::ranges::less(), []( vec2 v ){ return v.x(); } )->x() );

	if( lastOutputFrame * frameToTime() > 60.0f * 10.0f ) // Outfile longer than 10 minutes?
		{
		std::cout << "PV::modify tried to make a file longer than 10 minutes, which is currently disabled";
		return PV();
		}

	auto format = getFormat();
	format.numFrames = ceil( lastOutputFrame );
	PV out( format );
	out.clearBuffer();

	// Using a mutex per frame has proven to be faster in practice than other lock sets
	// A C array is used because mutex cannot be copied or moved, ruling out std::vector, and the size isn't known at compile-time, ruling out std::array
	// C++ should probably cover this use case but no std container does
	auto mutexBuffer = std::unique_ptr<std::mutex[]>( new std::mutex[format.numFrames] );

	std::vector<MF> inModified( inChannelDataCount );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		// Calculate and copy mapped frequencies and input magnitudes
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		auto modDataSamplesWriteHead = inModified.begin();
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				{
				*modDataSamplesWriteHead = 
					{
					inBufferReadHead->m,
					mod( frame * frameToTime(), inBufferReadHead->f ).y()
					};
				}
		
		std::for_each( std::execution::par, iota_iter(1), iota_iter(numFrames), [&]( Frame frame )
			{
			//std::for_each( std::execution::par, iota_iter(1), iota_iter(numBins), [&]( Bin bin )
			for( Bin bin = 1; bin < numBins; ++bin )
				{
				// For each square in the input, the mod maps that square to this quad
				const std::array<vec2, 4> p = {
					modSamples[ ( frame - 1 ) * numBins + bin - 1 ],
					modSamples[ ( frame - 0 ) * numBins + bin - 1 ],
					modSamples[ ( frame - 0 ) * numBins + bin - 0 ],
					modSamples[ ( frame - 1 ) * numBins + bin - 0 ] };
		
				const std::array<MF, 4> pMF = {
					inModified[ ( frame - 1 ) * numBins + bin - 1 ],
					inModified[ ( frame - 0 ) * numBins + bin - 1 ],
					inModified[ ( frame - 0 ) * numBins + bin - 0 ],
					inModified[ ( frame - 1 ) * numBins + bin - 0 ] };
		
				const vec2 D12 = p[1] - p[0];
				const vec2 D23 = p[2] - p[1];
				const vec2 D34 = p[3] - p[2];
				const vec2 D41 = p[0] - p[3];
			
				// Find bounding box containing quad to iterate over
				const Frame minx = fmax( floor( fmin( fmin( p[0].x(), p[1].x() ), fmin( p[2].x(), p[3].x() ) ) ), 0 );
				const Bin 	miny = fmax( floor( fmin( fmin( p[0].y(), p[1].y() ), fmin( p[2].y(), p[3].y() ) ) ), 0 );
				const Frame maxx = fmin( ceil(  fmax( fmax( p[0].x(), p[1].x() ), fmax( p[2].x(), p[3].x() ) ) ), format.numFrames - 1 );
				const Bin 	maxy = fmin( ceil(  fmax( fmax( p[0].y(), p[1].y() ), fmax( p[2].y(), p[3].y() ) ) ), numBins - 1 );
			
				// Iterate over bounding box
				for( Frame x = minx; x <= maxx; ++x )
					{
					for( Bin y = miny; y <= maxy; ++y )
						{	
						// Test if (x,y) is inside the quad
						// Ref: http://paulbourke.net/geometry/polygonmesh/#insidepoly
						bool c = false;
						if( ( ( p[0].y() <= y && y < p[3].y() ) || ( p[3].y() <= y && y < p[0].y() ) ) && ( x < D41.x() / D41.y() * ( y - p[0].y() ) + p[0].x() ) ) c = !c;
						if( ( ( p[1].y() <= y && y < p[0].y() ) || ( p[0].y() <= y && y < p[1].y() ) ) && ( x < D12.x() / D12.y() * ( y - p[1].y() ) + p[1].x() ) ) c = !c;
						if( ( ( p[2].y() <= y && y < p[1].y() ) || ( p[1].y() <= y && y < p[2].y() ) ) && ( x < D23.x() / D23.y() * ( y - p[2].y() ) + p[2].x() ) ) c = !c;
						if( ( ( p[3].y() <= y && y < p[2].y() ) || ( p[2].y() <= y && y < p[3].y() ) ) && ( x < D34.x() / D34.y() * ( y - p[3].y() ) + p[3].x() ) ) c = !c;
				
						if( c )
							{
							// Quadrilateral interpolation
							// Ref: https://www.particleincell.com/2012/quad-interpolation/
					
							const std::array<float, 4> alpha = { p[0].x(), p[1].x() - p[0].x(), p[3].x() - p[0].x(), p[0].x() - p[1].x() + p[2].x() - p[3].x() };
							const std::array<float, 4> beta  = { p[0].y(), p[1].y() - p[0].y(), p[3].y() - p[0].y(), p[0].y() - p[1].y() + p[2].y() - p[3].y() };
					
							const float quadA = alpha[3] * beta[2] - alpha[2] * beta[3];
							const float quadB = alpha[3] * beta[0] - alpha[0] * beta[3] 
											  + alpha[1] * beta[2] - alpha[2] * beta[1]
											  + x 		 * beta[3] - alpha[3] * y;
							const float quadC = alpha[1] * beta[0] - alpha[0] * beta[1]
											  + x		 * beta[1] - alpha[1] * y;
					
							float m;
							if( quadA == 0.0f )	
								{
								if( quadB == 0.0f )
									break;
								m = -quadC / quadB;
								}
							else
								{
								const float descriminant = quadB * quadB - 4.0f * quadA * quadC;
								if( descriminant < 0 ) break;
								m = ( -quadB + sqrt( descriminant ) ) / ( 2.0f * quadA ); // Quadratic formula
								}
							const float lDenominator = alpha[1] + alpha[3] * m;
							if( lDenominator == 0 ) break;
							const float l = ( x - alpha[0] - alpha[2] * m ) / lDenominator;
					
							// Make sure (l,m) is within a unit square + epsilon
							const float epsilon = 0.0001f;
							if( fabs( l - 0.5f ) > 0.5f + epsilon || fabs( m - 0.5f ) > 0.5f + epsilon ) break;

							const float interpL = interp( l );
							const float interpM = interp( m );

							const std::array<Magnitude, 4> w = {
								( 1.0f - interpL ) * ( 1.0f - interpM ) * pMF[0].m,
								(        interpL ) * ( 1.0f - interpM ) * pMF[1].m,
								(        interpL ) * (        interpM ) * pMF[2].m,
								( 1.0f - interpL ) * (        interpM ) * pMF[3].m };
							const float totalWeight = w[0] + w[1] + w[2] + w[3];
							if( totalWeight <= 0.0f ) break;

							// Choosing the MF assigned in this algorithm has many less than perfect options.
							// WFS: Weighted frequency sum gives inharmonic frequencies
							// MM: Max-mag requires tracking additional info
							// FIMM: In freq-innacurate max-mag the current quad's max-mag is compared to the write location. If it is louder,
							// overwrite frequency and then add to the magnitude. The problem with this strategy is that if a bin has quiet
							// data written to it more than once followed by loud data, a quiet frequency might be used. 
							// MIMM: Mag-innacurate max-mag is the same idea but with overwriting the magnitude instead of adding to it. 
							// Similar to freq-innacurate, this causes magnitude inaccuracies with repeated writes.
							// MIMM is being used currently.

							const auto maxWeightIter = std::max_element( w.begin(), w.end() );
							const int maxWeightedQuadMagIndex = std::distance( w.begin(), maxWeightIter );
							const Frequency loudestQuadFrequency = pMF[maxWeightedQuadMagIndex].f;
							
							// Lock mutex for frame x
							std::lock_guard<std::mutex> lock( mutexBuffer.get()[x] );
							MF & outMF = out.getMF( channel, x, y );
							if( *maxWeightIter > outMF.m )
								outMF = { *maxWeightIter, loudestQuadFrequency };
								

							// Weighted Frequency Sum: kind of a muddy sound
							// const float weightedFreqSum = 
							// 			  pMF[0].f * w[0]
							// 			+ pMF[1].f * w[1]
							// 			+ pMF[2].f * w[2]
							// 			+ pMF[3].f * w[3];
							// // Lock mutex for frame x
							// std::lock_guard<std::mutex> lock( mutexBuffer.get()[x] );
							// MF & outMF = out.getMF( channel, x, y );
							// outMF.f = ( outMF.f * outMF.m + weightedFreqSum ) / ( outMF.m + totalWeight );
							// outMF.m += totalWeight;
							}
						}
					}
				}
			} );
		}

	return out;
	}

PV PV::modifyFrequency( const Func2x1 & mod, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	PV out( getFormat() );
	out.clearBuffer();

	const Channel numChannels = getNumChannels();
	const Frame numFrames = getNumFrames();
	const Bin numBins = getNumBins();
	const uint32_t inChannelDataCount = numFrames * numBins;

	//Sample mod function and convert output to frame/bin
	auto modSamples = mod.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
	std::for_each( std::execution::par_unseq, modSamples.begin(), modSamples.end(), [&]( Frequency & f ){ 
		f *= frequencyToBin(); 
		} );

	std::vector<MF> inModified( inChannelDataCount ); // Buffer for sampling mod on input frequencies

	for( Channel channel = 0; channel < numChannels; ++channel )
		{
		// Sample mod over input frequencies into inModified
		// This is computed early because the mod function isn't guaranteed to be thread safe
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		auto modDataSamplesWriteHead = inModified.begin();
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				*modDataSamplesWriteHead = {
					inBufferReadHead->m,
					mod( frame * frameToTime(), inBufferReadHead->f ) 
					};

		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numFrames ), [&]( Frame frame )
			{
			// For each adjacent pair of bins
			for( Bin bin = 1; bin < numBins; ++bin )
				{
				// Get where those bins were mapped to by the mod function
				const int hiBinIndex = frame * numBins + bin;
				const float loBin = modSamples[hiBinIndex - 1];
				const float hiBin = modSamples[hiBinIndex    ];
				const bool forward = hiBin > loBin; // Check if the bins were mapped upside down

				const Bin loBinRound = forward ? std::ceil( loBin ) : std::floor( loBin );
				const Bin hiBinRound = forward ? std::ceil( hiBin ) : std::floor( hiBin );
				const Bin startBin = std::clamp( loBinRound, 0, out.getNumBins() - 1 );
				const Bin endBin   = std::clamp( hiBinRound, 0, out.getNumBins() - 1 );

				const MF & loMF = inModified[ buffer_access( bin - 1, frame, getNumBins() ) ];
				const MF & hiMF = inModified[ buffer_access( bin    , frame, getNumBins() ) ];

				for( Bin y = startBin; y != endBin; forward? ++y : --y )
					{
					const float mix = interp( ( float( y ) - loBin ) / ( hiBin - loBin ) );

					const Magnitude w0 = ( 1.0f - mix ) * loMF.m;
					const Magnitude w1 = (        mix ) * hiMF.m;
					
					const MF & maxMagMF = w0 < w1 ? loMF : hiMF;
					MF & outMF = out.getMF( channel, frame, y );
					if( maxMagMF.m > outMF.m )
						{
						outMF.m += maxMagMF.m;
						outMF.f = maxMagMF.f;
						}

					// WFS:
					// const Magnitude totalWeight = w0 + w1;
					// const float weightedFreqSum = w0 * loMF.f + w1 * hiMF.f;
					// auto & outMF = out.getMF( channel, frame, y );
					// outMF.f = ( outMF.m * outMF.f + weightedFreqSum ) / ( outMF.m + totalWeight );
					// outMF.m += totalWeight;
					}
				}
			} );
		}

	return out;
	}

PV PV::modifyTime( const Func2x1 & mod, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	const uint32_t numChannels = getNumChannels();
	const uint32_t numFrames = getNumFrames();
	const uint32_t numBins = getNumBins();

	// Sample mod function and convert output to frame/bin
	auto modSamples = mod.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
	std::for_each( std::execution::par_unseq, modSamples.begin(), modSamples.end(), [&]( Time & t ){ 
		t *= timeToFrame();
		} );

	// Get the largest frame value among all mod samples, rounded up
	const float lastOutputFrame = std::ceil( *std::ranges::max_element( modSamples ) );

	auto format = getFormat();
	format.numFrames = lastOutputFrame;
	PV out( format );
	out.clearBuffer();

	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numChannels ), [&]( Channel channel )
		{
		// Note the bin loop happens first here. There are a lot of considerations involved in the ordering, but
		//  it comes down to time displacement happening along frames. Paralellizing the frames requires synchronization
		//  that isn't worth the benifits, and the parallel loop being first is faster. The algorithm isn't cache friendly 
		//  either way because we store PV data in frame-major order.
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numBins ), [&]( Bin bin )
			{
			// For each adjacent pair of frames
			std::for_each( std::execution::seq, iota_iter( 1 ), iota_iter( numFrames ), [&]( Frame frame )
				{
				const float lFrame = modSamples[ buffer_access( bin, frame - 1, numBins ) ];
				const float rFrame = modSamples[ buffer_access( bin, frame    , numBins ) ];
				const bool forward = rFrame > lFrame;

				const Frame startFrame = forward ? std::ceil( lFrame ) : std::floor( lFrame );
				const Frame endFrame   = forward ? std::ceil( rFrame ) : std::floor( rFrame );

				const MF & lMF = getMF( channel, frame - 1, bin );
				const MF & rMF = getMF( channel, frame,     bin );

				for( Frame x = startFrame; x != endFrame; forward? ++x : --x )
					{
					if( x < 0 || out.getNumFrames() <= x ) continue;	

					const float mix = interp( ( x - lFrame ) / ( rFrame - lFrame ) );
					const Magnitude w0 = ( 1.0f - mix ) * lMF.m;
					const Magnitude w1 = (        mix ) * rMF.m;
					const Magnitude totalWeight = w0 + w1;
					const float weightedFreqSum = w0 * lMF.f + w1 * rMF.f;

					auto & outMF = out.getMF( channel, x, bin );
					outMF.f = ( outMF.f * outMF.m + weightedFreqSum ) / ( outMF.m + totalWeight );
					outMF.m += totalWeight;
					}
				} );
			} );
		} );

	return out;
	}

PV PV::repitch( const Func2x1 & factor, Interpolator interp ) const
	{
	flan_FUNCTION_LOG;
	return modifyFrequency( [&]( vec2 tf ){ return factor( tf ) * tf.f(); }, interp );
	}

PV PV::stretch( const Func2x1 & factor, Interpolator interp ) const
	{
	flan_FUNCTION_LOG;
	return modifyTime( [&factor]( vec2 tf ){ return factor( tf.x(), tf.y() ) * tf.x(); }, interp );
	}

PV PV::stretch_spline( const Func1x1 & interpolation ) const
	{
	flan_PROCESS_START( PV() );

	const auto safeInterpolation = [&interpolation, this]( Frame frame )
		{
		return std::max( uint32_t( interpolation( frame * frameToTime() ) ), 1u );
		};
	
	//Set up output and spline x coordinates
	auto format = getFormat();
	format.numFrames = 0;

	std::vector<double> Xs( getNumFrames() );
	for( Frame frame = 0; frame < getNumFrames()-1; ++frame )
		{
		Xs[frame] = format.numFrames;
		format.numFrames += safeInterpolation( frame );
		}
	Xs[getNumFrames()-1] = format.numFrames;
	PV out( format );

	//Allocate Y coordinate vectors
	std::vector<double> magnitudeYs( Xs.size() );
	std::vector<double> frequencyYs( Xs.size() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Channel -> bin is correct here, despite non-sequential access to flan buffer
		for( Bin bin = 0; bin < getNumBins(); ++bin )
			{
			//For each x coordinate get corresponding frequency and magnitude
			for( Frame frame : iota_view( 0, getNumFrames() ) )
				{
				magnitudeYs[frame] = getMF( channel, frame, bin ).m;
				frequencyYs[frame] = getMF( channel, frame, bin ).f;
				}

			//Do the splines
			tk::spline magnitudeSpline, frequencySpline;
			magnitudeSpline.set_points(Xs,magnitudeYs);
			frequencySpline.set_points(Xs,frequencyYs);

			//Evaluate splines into out
			for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
				{
				out.setMF( channel, frame, bin, 
					{ 
					(float) magnitudeSpline(frame),
					(float) frequencySpline(frame)
					} );
				}
			}
		}
			
	return out;
	}

PV PV::desample( const Func2x1 & factor, Interpolator interp ) const
	{
	const Channel numChannels = getNumChannels();
	const Frame numFrames = getNumFrames();
	const Bin numBins = getNumBins();

	// Sample factor over frames and bins
	auto factorSamples = factor.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );

	flan_PROCESS_START( PV() );

	PV out( getFormat() );
	out.clearBuffer();

	for( Channel channel = 0; channel < numChannels; ++channel )
		{
		// This stores which frames should be interpolated between
		std::vector<std::vector<Frame>> selectedFrames( numBins );

		// This stores an accumulator for each bin which is integrated into
		std::vector<float> accums( numBins, 1 ); // Start at one to trigger interp from starting frame

		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numBins ), [&]( Bin bin )
			{
			// Factor is integrated, when the accumulator passes 1 the currest frame is selected as an interpolation endpoint
			for( Frame frame = 0; frame < numFrames; ++frame )
				{
				float & accum = accums[bin];
				if( factorSamples[ frame * numBins + bin ] < 1.0f ) continue;
				const float factor_c = factorSamples[ frame * numBins + bin ];
				accum += factor_c;
				if( accum >= 1.0f )
					{
					selectedFrames[bin].push_back( frame );
					accum -= 1.0f; // No need to worry about accum being bigger than 2, factor_c can't be bigger than 1
					}
				}
			} );

		//for( Bin bin = 0; bin < numBins; ++bin )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numBins ), [&]( Bin bin )
			{
			// For each pair of selected interpolation endpoints in selectedFrames, intepolate the frames between them
			std::vector<Frame> & selectedFrames_c = selectedFrames[bin];
			if( selectedFrames_c.size() < 2 ) return;
			for( const int selectedFrameIndex : std::views::iota( 0ull, selectedFrames_c.size() - 1 ) )
				{
				const Frame lFrame = selectedFrames_c[selectedFrameIndex    ];
				const Frame rFrame = selectedFrames_c[selectedFrameIndex + 1];
				const MF lMF = getMF( channel, lFrame, bin );
				const MF rMF = getMF( channel, rFrame, bin );
				for( Frame frame = lFrame; frame < rFrame; ++frame )
					{
					const float mix = interp( float( frame - lFrame ) / ( rFrame - lFrame ) ); 
					const float w0 = ( 1.0f - mix ) * lMF.m;
					const float w1 = (        mix ) * rMF.m;
					const MF mf = 
						{
						w0 + w1,
						w0 > w1 ? lMF.f : rMF.f
						};
					// Weighting by frequency is another option, it doesn't noticably change the output for inputs I've tested though
					// const MF mf = 
					// 	{
					// 	w0 + w1,
					// 	( w0 * lMF.f + w1 * rMF.f ) / ( w0 + w1 )
					// 	};
					out.setMF( channel, frame, bin, mf );
					}
				}
			} );
		}

	return out;
	}

PV PV::timeExtrapolate( Time startTime, Time endTime, Time extrapolationTime, Interpolator interpolator ) const
	{
	flan_PROCESS_START( PV() );

	// Input validation
	startTime = std::clamp( startTime, 0.0f, getLength() );
	if( endTime == -1 ) endTime = getLength();
	endTime	  = std::clamp( endTime,   0.0f, getLength() );
	if( startTime >= endTime ) return PV();
	if( extrapolationTime <= 0 ) return PV();

	const Frame startFrame	= timeToFrame() * startTime;			//x_1
	const Frame endFrame	= timeToFrame() * endTime;				//x_2
	const Frame extFrames	= timeToFrame() * extrapolationTime;	//Extrapolated frames to generate

	auto format = getFormat();
	format.numFrames = endFrame + extFrames;
	PV out( format );
	out.clearBuffer();

	// Sample interpolator
	std::vector<float> interpSamples( out.getNumFrames() - startFrame );
	for( Frame frame = 0; frame < interpSamples.size(); ++frame )
		interpSamples[frame] = interpolator( float( frame - startFrame ) / float( endFrame - startFrame ) );

	for( Channel channel = 0; channel < getNumChannels(); ++channel ) 
		{
		// Copy up to the start frame from input into output
		std::copy( std::execution::par_unseq, channelBegin( channel ), channelBegin( channel ) + startFrame * getNumBins(), out.channelBegin( channel ) );
			
		//Interpolate and continue beyond endFrame to extrapolate
		std::for_each( std::execution::par_unseq, iota_iter( startFrame ), iota_iter( out.getNumFrames() ), [&]( Frame frame ) 
			{
			const float mix = interpSamples[ frame - startFrame ];

			for( Bin bin = 0; bin < getNumBins(); ++bin ) 
				{
				const auto leftMF  =  getMF( channel, startFrame, bin );
				const auto rightMF =  getMF( channel, endFrame, bin );

				// This algorithm decides bin based on frequency. Phase vocoder outputs can have multiple bins in a region
				// containing the same frequency, and squishing all those into one bin will cause the output to sound less full.
				// To fix that, we note for each bin at endTime how far off the frequency in that bin is from where it "should" be. That shift
				// is later added back in when deciding what bin should be used for the current frequency, maintaining fidelity.
				const Bin rightBinShift = bin - rightMF.f * frequencyToBin();

				const MF extrapMF = { std::abs( ( 1.0f - mix ) * leftMF.m + mix * rightMF.m ),
					          					( 1.0f - mix ) * leftMF.f + mix * rightMF.f };
				const Bin extrapBin = extrapMF.f * frequencyToBin() + rightBinShift;
				if( extrapBin < 0 || extrapBin >= getNumBins() ) continue;

				MF & outMF = out.getMF( channel, frame, extrapBin );

				if( extrapMF.m > outMF.m ) 
					outMF = extrapMF;
				}
			} );
		}
	return out;
	}

}