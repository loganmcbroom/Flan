#include "flan/Audio/Audio.h"

#include <iostream>

#include "WDL/resample.h"
#include "r8brain/CDSPResampler.h"

#include "flan/WindowFunctions.h"
#include "flan/Utility/iota_iter.h"

using namespace flan;

Audio Audio::resample( FrameRate newSampleRate ) const
	{
	flan_PROCESS_START( Audio() );

	if( newSampleRate == getSampleRate() ) 
		return copy();

	auto format = getFormat();
	format.numFrames *= newSampleRate / float( getSampleRate() );
	format.sampleRate = newSampleRate;
	Audio out( format );

	r8b::CDSPResampler resampler( getSampleRate(), newSampleRate, getNumFrames() );
	resampler.oneshot<float, float>( &getBuffer()[0], getBuffer().size(), &out.getBuffer()[0], out.getBuffer().size() );

	return out;
	}

Audio Audio::convertToMidSide() const
	{
	flan_PROCESS_START( Audio() );

	if( getNumChannels() != 2 )
		{
		std::cout << "Can't transform non-stereo Audio between Mid-Side and Left-Right formats." << std::endl;
		return copy();
		}
	Audio out( getFormat() );
	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
		{
		out.setSample( 0, frame, getSample( 0, frame ) + getSample( 1, frame ) );
		out.setSample( 1, frame, getSample( 0, frame ) - getSample( 1, frame ) );
		} );
	return out;
	}

Audio Audio::convertToLeftRight() const
	{
	flan_FUNCTION_LOG;
	return convertToMidSide();
	}

Audio Audio::convertToStereo() const
	{
	flan_PROCESS_START( Audio() );

	auto format = getFormat();
	format.numChannels = 2;
	Audio out( format );

	switch( getNumChannels() )
		{
		case 1:
			{
			// Copy input channel one into the output channel one while scaling down by sqrt2
			std::transform( std::execution::par_unseq, getBuffer().begin(), getBuffer().end(), out.getBuffer().begin(), 
				[sqrt2 = std::sqrt( 2.0f )]( Sample s ){ return s / sqrt2; } );

			// // Copy output channel one to output channel two
			std::copy( std::execution::par_unseq, out.getBuffer().begin(), out.getBuffer().begin() + getNumFrames(), out.getBuffer().begin() + getNumFrames() );
			break;
			}
		case 2:
			return copy();
		default:
			std::cout << "I don't know how to convert that number of channels to stereo." << std::endl;
			
		}
	return out;
	}

Audio Audio::convertToMono() const
	{
	flan_PROCESS_START( Audio() );

	auto format = getFormat();
	format.numChannels = 1;
	Audio out( format );

	const float sqrtChannels = std::sqrt( getNumChannels() );

	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
		{
		Sample sampleAccumulator = 0;
		for( Channel channel = 0; channel < getNumChannels(); ++channel )
			sampleAccumulator += getSample( channel, frame );
		out.setSample( 0, frame, sampleAccumulator / sqrtChannels );
		} );

	return out;
	}

Func1x1 Audio::convertToFunction( Second granularity ) const
	{
	flan_PROCESS_START( 0 );

	if( granularity <= 0 ) return 0;

	const fFrame granularityFrames = timeToFrame( granularity );

	auto safeGetSample = [this]( Frame frame )
		{
		if( frame < 0 || frame >= getNumFrames() )
			return 0.0f;
		else
			return getSample( 0, frame );
		};

	const Frame numSamples = std::ceil( float( getNumFrames() ) / granularityFrames );
	std::vector<float> ys( numSamples );

	for( Frame x = 0; x < numSamples; ++x )
		{
		const fFrame sampleFrame = x * granularityFrames;
		float sum = 0;
		for( Frame i = - granularityFrames; i <= granularityFrames; ++i )
			{
			const Sample sample = std::abs( safeGetSample( sampleFrame + i ) );
			const Sample windowSample = Windows::Hann( ( float( i ) / granularityFrames + 1 ) * .5f );
			sum += sample * windowSample;
			}
		
		// Division by granularityFrames is to account for the weighting via the hann window
		// The integral of that window is 1/2, but we are stretching it by granularity frames in both directions
		// The pi/2 factor handles that the average absolute value of a sinusoid with amplitude 1 is 2/pi
		ys[x] = sum / granularityFrames * pi / 2.0f;
		}
	
	const float timeToFrame_f = timeToFrame( 1 );
	return [ys = std::move( ys ), timeToFrame_f, granularityFrames ]( float t )
		{
		// Get the image sample index as a float
		const float x = t * timeToFrame_f / granularityFrames;
		if( 0 <= x && x < ys.size() - 1 )
			{
			Frame x1 = std::floor( x );
			float y1 = ys[x1];
			float y2 = ys[x1 + 1];
			return y1 + ( y2 - y1 ) * ( x - x1 );
			}
		else return 0.0f;
		};
	}

// class QVOCBuffer {
// 	public:
// 		/*
// 		There is no natural ordering of QVOC data. The data is stored here in a single buffer so that ascending channel -> frame -> bin 
// 		iteration gives sequential memory access.
// 		*/

// 		// Takes a single channel of data to collect size info, doesn't actually get the data
// 		QVOCBuffer( const std::vector<std::vector<std::complex<float>>> & cqOut )
// 			{
// 			// Fill index buffer with the indices of the starts of each frame of data
// 			index_buffer.reserve( cqOut.size() );
// 			index_buffer.push_back( 0 );
// 			std::transform( cqOut.begin(), cqOut.end(), std::back_inserter( index_buffer ), 
// 				[]( const std::vector<std::complex<float>> & v ){ return v.size(); } );

// 			std::partial_sum( index_buffer.begin(), index_buffer.end(), index_buffer.begin() );

// 			buffer.resize( getChannelSize() );
// 			}

// 		MF & getMF( Channel channel, Frame frame, Bin bin )
// 			{
// 			return buffer[ getChannelSize() * channel + index_buffer[frame] + bin ];
// 			}

// 		size_t getChannelSize() const { return index_buffer.back(); }
// 		Frame getNumFrames() const { return index_buffer.size() - 1; }
// 		Bin getNumBins( Frame f ) const { return index_buffer[f+1] - index_buffer[f]; }

// 	private:
// 		std::vector<size_t> index_buffer;
// 		std::vector<MF> buffer;
// };

// #include "constantQ/src/ConstantQ.h"
// #include "constantQ/src/CQInverse.h"
// Audio Audio::convertToQVOC() const
// 	{
// 	CQParameters params( getSampleRate(), 20, 20000, 128 );

// 	// Forward constant-q transform
// 	ConstantQ cq( params );
// 	auto cqOut = cq.process( getBuffer() );
// 	auto cqOutR = cq.getRemainingOutput();
// 	cqOut.insert( cqOut.end(), cqOutR.begin(), cqOutR.end() );
// 	auto cqOutProperties = cq.getKernalProperties();
// 	const float dftSize = cqOutProperties.fftSize; 

// 	// Transform from cq output format to QVOCBuffer (MF format)
// 	QVOCBuffer QVOC( cqOut );
// 	std::vector<float> phaseBuffer( cq.getTotalBins(), 0 );
// 	for( Frame frame = 0; frame < cqOut.size(); ++frame )
// 		{
// 		for( Bin bin = 0; bin < cqOut[frame].size(); ++bin )
// 			{
// 			const Bin flippedBin = cqOut[frame].size() - 1 - bin;
// 			const float hopSize = cq.getColumnHop( bin );

// 			const float magnitude = std::abs( cqOut[frame][flippedBin] );

// 			// We can't use frequency estimation with too low of a temporal resolution
// 			if( hopSize >= dftSize )
// 				QVOC.getMF( 0, frame, bin ) = { magnitude, cq.binToFrequency( bin ) };

// 			const float phase = std::arg( cqOut[frame][flippedBin] );
// 			const float phaseDiff = phase - phaseBuffer[bin]; // Actual change in radians for the current hop
// 			phaseBuffer[bin] = phase; // Update phase buffer
// 			const float expectedPhaseDiff = bin * pi2 * hopSize / dftSize; // Expected change in radians for the current hop.
// 			const float deltaPhase = phaseDiff - expectedPhaseDiff; // The difference between the actual and the expected phase changes, aka how far off was our sinusoid from the expected.
// 			const float wrappedDeltaPhase = deltaPhase - pi2 * std::round( deltaPhase / pi2 ); // Wrap delta phase into [-pi,pi]
// 			const float binDeviation = wrappedDeltaPhase / pi2 * dftSize / hopSize;
// 			const Frequency frequency = cq.binToFrequency( bin + binDeviation );

// 			QVOC.getMF( 0, frame, bin ) = { magnitude, frequency };
// 			}
// 		}

// 	// Modify qvoc data
// 	for( Frame frame = 0; frame < QVOC.getNumFrames(); ++frame )
// 		{
// 		for( Bin bin = 0; bin < QVOC.getNumBins( frame ) - 1; ++bin )
// 			{
// 			// QVOC.getMF( 0, frame, bin ) = QVOC.getMF( 0, frame, bin + 1 );
// 			// QVOC.getMF( 0, frame, bin ).f /= std::pow( 2.0f, 1.0f / 12 );
// 			}
// 		}

// 	// Transform QVOCBuffer back to cq format
// 	std::for_each( phaseBuffer.begin(), phaseBuffer.end(), []( float & x ){ x = 0; } );

// 	// We need to construct a vector of vectors with the correct sizes
// 	std::vector<std::vector<std::complex<float>>> cqIn( QVOC.getNumFrames(), std::vector<std::complex<float>>{} );
// 	for( Frame frame = 0; frame < cqIn.size(); ++frame )
// 		cqIn[frame].resize( QVOC.getNumBins( frame ), 0 );

// 	// Apply the MF->complex transform
// 	for( Frame frame = 0; frame < cqIn.size(); ++frame )
// 		{
// 		for( Bin bin = 0; bin < cqIn[frame].size(); ++bin )
// 			{
// 			const Bin flippedBin = cqOut[frame].size() - 1 - bin;
// 			const float hopSize = cq.getColumnHop( bin );

// 			const MF mf = QVOC.getMF( 0, frame, bin );

// 			phaseBuffer[bin] += cq.frequencyToBin( mf.f ) * pi2 * hopSize / dftSize;
			
// 			cqIn[frame][flippedBin] = std::polar( mf.m, phaseBuffer[bin] );
// 			}
// 		}

// 	// Inverse constant-q transform
// 	CQInverse cqi( params );
// 	auto cqiOut = cqi.process( cqIn );
// 	auto cqiOutR = cqi.getRemainingOutput();
// 	cqiOut.insert( cqiOut.end(), cqiOutR.begin(), cqiOutR.end() );

// 	// Create output Audio
// 	auto format = getFormat();
// 	format.numFrames = cqiOut.size();
// 	Audio out( format );

// 	// Copy to output
// 	for( Frame frame = 0; frame < cqiOut.size(); ++frame )
// 		out.getSample( 0, frame ) = cqiOut[frame];

// 	return out;
// 	}

