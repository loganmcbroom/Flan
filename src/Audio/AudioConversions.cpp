#include "flan/Audio.h"

#include <iostream>

#include "bmp/bitmap_image.hpp"
#include "flan/WindowFunctions.h"
#include "flan/PVOC.h"
#include "flan/CLContext.h"
#include "flan/CLProgs.h"
#include "flan/FFTHelper.h"
#include "flan/Graph.h"

static const float pi = std::acos( -1.0f );
static const float pi2 = pi * 2.0f;

using namespace flan;

Audio Audio::convertToMidSide( flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	if( getNumChannels() != 2 )
		{
		std::cout << "Can't transform non-stereo Audio between Mid-Side and Left-Right formats." << std::endl;
		return *this;
		}
	Audio out( getFormat() );
	for( Frame frame = 0; frame < getNumFrames(); ++frame )
		{
		flan_CANCEL_POINT( Audio() );
		out.setSample( 0, frame, getSample( 0, frame ) + getSample( 1, frame ) );
		out.setSample( 1, frame, getSample( 0, frame ) - getSample( 1, frame ) );
		}
	return out;
	}

Audio Audio::convertToLeftRight( flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;
	return convertToMidSide( canceller );
	}

Audio Audio::convertToStereo( flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	auto format = getFormat();
	format.numChannels = 2;
	Audio out( format );

	switch( getNumChannels() )
		{
		case 1:
			{
			for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
				{
				flan_CANCEL_POINT( Audio() );
				out.setSample( 0, frame, getSample( 0, frame ) / sqrt(2.0f) );
				out.setSample( 1, frame, getSample( 0, frame ) / sqrt(2.0f) );
				}
			break;
			}
		case 2:
			return *this;
		default:
			std::cout << "I don't know how to convert that number of channels to stereo." << std::endl;
			
		}
	return out;
	}

Audio Audio::convertToMono( flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	auto format = getFormat();
	format.numChannels = 1;
	Audio out( format );

	const float sqrtChannels = std::sqrt( getNumChannels() );

	for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
		{
		flan_CANCEL_POINT( Audio() );
		Sample sampleAccumulator = 0;
		for( Channel channel = 0; channel < getNumChannels(); ++channel )
			sampleAccumulator += getSample( channel, frame );
		out.setSample( 0, frame, sampleAccumulator / sqrtChannels );
		}

	return out;
	}

Graph Audio::convertToGraph( Interval I, Pixel width, Pixel height, float timelineScale, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Graph( width, height ) );

	if( I.x2 == -1 ) I.x2 = getLength();
	
	const Frame startFrame = I.x1 * timeToFrame();
	const Frame endFrame = I.x2 * timeToFrame();

	Graph g( width, height );
	g.fillImage( Color::fromHSV( 0, 0, .04 ) );
	g.addFullSplitViewY( I * Interval( -1, 1 ), getNumChannels() );
		
	std::vector<const float *> channelStarts( getNumChannels() );
	for( int i = 0; i < getNumChannels(); ++i )
		channelStarts[i] = getSamplePointer( i, 0 );
	g.drawWaveforms( channelStarts, getNumFrames(), { 0, -1, getLength(), 1 }, 0, Graph::WaveformMode::Symmetric, canceller );

	if( timelineScale > 0 )
		{
		const float bigJump = std::pow( 4.0f, std::floor( std::log2( I.w() ) / 2 - 0.5f ) );
		g.drawXTicks( bigJump / 4.0f, 1, timelineScale / 2, 0, -1, Color::fromHSV( 0, 0, .6 ), 0 );
		g.drawXTicks( bigJump, 1, timelineScale, 0, -1, Color::fromHSV( 0, 0, 1 ), timelineScale );
		}

	return g;
	}

Audio Audio::saveToBMP( const std::string & filename, Interval I, int width, int height, flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;
	auto b = convertToGraph( I, width, height, canceller );
	b.save_image( filename );
	return *this;
	}

PVOC Audio::convertToPVOC_cpu( Frame windowSize, Frame hopSize, Frame dftSize, std::shared_ptr<FFTHelper> fft, flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;

	const Bin numBins = dftSize / 2 + 1;
	const Frame numFrames = getNumFrames();
	//+1 since we analyze at start and end times
	const uint32_t numHops = uint32_t( std::ceil( numFrames / hopSize ) ) + 1; 
	const float binWidth = float( getSampleRate() ) / float( dftSize );
	
	// Set up output
	PVOCBuffer::Format PVOCFormat;
	PVOCFormat.numChannels = getNumChannels();
	PVOCFormat.numFrames = numHops;
	PVOCFormat.numBins = numBins;
	PVOCFormat.sampleRate = getSampleRate();
	PVOCFormat.hopSize = hopSize;
	PVOCFormat.windowSize = windowSize;
	PVOC out( PVOCFormat );

	//Sample Hann window
	std::vector<float> hannWindow( windowSize );
	for( uint32_t i = 0; i < windowSize; ++i )
		hannWindow[i] = Windows::Hann( float( i ) / float( windowSize - 1 ) );

	//Allocate fft buffers and phase buffer
	std::vector<float> phaseBuffer( numBins );
	if( ! fft ) fft = std::make_shared<FFTHelper>( dftSize, true, false, false );

	//For each channel, do the whole thing
	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Set initial phase to 0
		for( Bin bin = 0; bin < numBins; ++bin )
			phaseBuffer[bin] = 0;

		//For each hop, fft and save into buffer
		for( int hop = 0; hop < numHops; ++hop )
			{
			flan_CANCEL_POINT( out );

			//Fill fft input buffer with windowed signal
			const Frame inFrameStart = hopSize * hop - windowSize / 2;
			const Frame inFrameEnd = inFrameStart + windowSize;
			const Frame inFrameStart_s = std::max( inFrameStart, 0 );
			const Frame inFrameEnd_s   = std::min( inFrameEnd, int( numFrames ) );

			const int32_t fftStart = inFrameStart_s - inFrameStart;
			const int32_t fftEnd = inFrameEnd_s - inFrameStart;

			// Clear any fft buffer frames that won't be written to
			if( inFrameStart != inFrameStart_s )
				std::fill( fft->realBegin(), fft->realBegin() + fftStart, 0 );
			if( inFrameEnd != inFrameEnd_s )
				std::fill( fft->realBegin() + fftEnd, fft->realEnd(), 0 );

			// Write input frames to fft buffer
			const Sample * inFrameRead = getSamplePointer( channel, inFrameStart_s );
			for( Frame fftFrame = fftStart; fftFrame < fftEnd; ++fftFrame )
				{
				fft->realBuffer[fftFrame] = *inFrameRead * hannWindow[fftFrame];
				++inFrameRead;
				}
	
			fft->r2cExecute();

			//For each bin, compute frequency and magnitude
			for( Bin bin = 0; bin < numBins; ++bin )
				{
				//	First, get the bin phase for the current and previous frame. Find the difference.
				//	We are expecting the phase of the ideal bin wave to move around due to using overlapping windows, so
				//    figure out what that expected movement is (or at least an equivalent angle).
				//  Next we get deltaPhase, this tells us how far from the expected phase shift our partial strayed.
				//	Up to this point we have been using angles potentially outside [-pi,pi], e.g. expectedPhaseDiff
				//    might be quite large when we really mean to talk about the equivalent angle in [-pi,pi]. We
				//	  fix all of this by wrapping deltaPhase into [-pi,pi]
				//  Next we need how far our true frequency deviates from the ideal bin center frequency.
				//    Dividing by 2pi gives a deviation in [-1/2,1/2], and we multiply by overlaps as our actual bin
				//    deviation is overlaps times larger than was accounted for (recall we measured the phase
				//    difference between two overlapped frames).
				//  Finally, the actual computed frequence is the bin center plus the deviation from that center
				//    times the width of a single bin.
				const float bin_f = float( bin );

				const float magnitude = std::abs( fft->complexBuffer[bin] );

				const float phase = std::arg( fft->complexBuffer[bin] ); // Current phase
				const float phaseDiff = phase - phaseBuffer[bin]; // Actual change in radians per hop
				phaseBuffer[bin] = phase; // Update phase buffer
				const float expectedPhaseDiff = bin_f * pi2 * hopSize / dftSize; // Expected change in radians per hop for this bin.
				const float deltaPhase = phaseDiff - expectedPhaseDiff; // The difference between the actual and the expected phase changes, aka how far off was our sinusoid from the expected.
				const float wrappedDeltaPhase = deltaPhase - pi2 * std::round( deltaPhase / pi2 ); // Wrap delta phase into [-pi,pi]
				const float binDeviation = wrappedDeltaPhase / pi2 * dftSize / hopSize;
				const float frequency = ( bin_f + binDeviation ) * binWidth;

				out.setMF( channel, hop, bin, { magnitude, frequency } );
				}
			}
		}

	return out;
	}

PVOC Audio::convertToPVOC( Frame windowSize, Frame hopSize, Frame dftSize, std::shared_ptr<FFTHelper> fft, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	if( dftSize < windowSize )
		dftSize = windowSize;

#ifndef USE_OPENCL
	return convertToPVOC_cpu( windowSize, hopSize, dftSize, fft, canceller );
#else

	if( ! isOpenCLAvailable() )
		{
		std::cout << "OpenCL Unavailable, using cpu backup routine";
		return convertToPVOC_cpu( windowSize, hopSize, dftSize, fft, canceller );
		}

	const Bin numBins = dftSize / 2 + 1;
	const Frame numFrames = getNumFrames();
	//+1 since we analyze at start and end times
	const uint32_t numHops = uint32_t( ceil( numFrames / hopSize ) ) + 1; 
	const Frequency binWidth = float( getSampleRate() ) / float( dftSize );
	const float overlaps = float( dftSize ) / hopSize; // Used in cl code
	
	// Set up output
	PVOCBuffer::Format PVOCFormat;
	PVOCFormat.numChannels = getNumChannels();
	PVOCFormat.numFrames = numHops;
	PVOCFormat.numBins = numBins;
	PVOCFormat.sampleRate = getSampleRate();
	PVOCFormat.hopSize = hopSize;
	PVOCFormat.windowSize = windowSize;
	PVOC out( PVOCFormat );
	PVOCBuffer::MF * outBufferWriteHead = out.getMFPointer( 0, 0, 0 );

	//Sample Hann window
	std::vector<float> hannWindow( windowSize );
	for( uint32_t i = 0; i < windowSize; ++i )
		hannWindow[i] = Windows::Hann( float( i ) / float( windowSize - 1 ) );

	// Prepare fftw
	if( ! fft  ) fft = std::make_shared<FFTHelper>( dftSize, true, false, false );
	
	// Prepare OpenCL
	CLContext cl;
	ProgramHelper programHelper( cl, CLProgs::Audio_convertToPVOC );
	const uint32_t outChannelDataCount = numBins * numHops;
	cl::Buffer clIn( cl.context, CL_MEM_READ_WRITE, sizeof( std::complex<float> ) * outChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::KernelFunctor< cl::Buffer > cl_fftToPhase( programHelper.program, "fftToPhase" );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, float, float, int > cl_phaseToFreq( programHelper.program, "phaseToFreq" );
	
	for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
		{
		for( int hop = 0; hop < numHops; ++hop )
			{
			//Fill fft input buffer with windowed signal
			const Frame inFrameStart = hopSize * hop - windowSize / 2;
			const Frame inFrameEnd = inFrameStart + windowSize;
			const Frame inFrameStart_s = std::max( inFrameStart, 0 );
			const Frame inFrameEnd_s   = std::min( inFrameEnd, int( numFrames ) );

			const int32_t dataStart = inFrameStart_s - inFrameStart;
			const int32_t dataEnd = inFrameEnd_s - inFrameStart;

			// Clear any fft buffer frames that won't be written to
			if( inFrameStart != inFrameStart_s )
				std::fill( fft->realBegin(), fft->realBegin() + dataStart, 0 );
			if( inFrameEnd != inFrameEnd_s )
				std::fill( fft->realBegin() + dataEnd, fft->realEnd(), 0 );

			// Write input frames to fft buffer
			const Sample * inFrameRead = getSamplePointer( channel, inFrameStart_s );
			for( Frame fftFrame = dataStart; fftFrame < dataEnd; ++fftFrame )
				{
				fft->realBuffer[fftFrame] = *inFrameRead * hannWindow[fftFrame];
				++inFrameRead;
				}

			flan_CANCEL_POINT( PVOC() );
			cl.queue.enqueueBarrierWithWaitList(); // Can't start writing to fftOut until copying to device from previous iteration completes
			fft->r2cExecute();

			//Copy fft output to device memory
			cl.queue.enqueueWriteBuffer( clIn, false, 
				sizeof( std::complex<float> ) * hop * numBins, 
				sizeof( std::complex<float> ) * numBins, 
				fft->complexBuffer );
			}

		//transform spectral buffers into phase/mag buffers
		flan_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl_fftToPhase( cl::EnqueueArgs( cl.queue, cl::NDRange( outChannelDataCount ) ), clIn );

		//second call computes all frequencies ( and magnitudes to avoid 2 copies from device ) into second buffer
		flan_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl_phaseToFreq( cl::EnqueueArgs( cl.queue, cl::NDRange( outChannelDataCount ) ), 
			clIn, 
			clOut,
			overlaps,
			binWidth,
			numBins );

		flan_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut, outBufferWriteHead, outBufferWriteHead + outChannelDataCount );
		outBufferWriteHead += outChannelDataCount;
		};

	flan_CANCEL_POINT( PVOC() );
	cl.queue.enqueueBarrierWithWaitList();
	return out;
#endif
	}

Func1x1 Audio::convertToFunction( Time granularity, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( 0 );

	if( granularity <= 0 ) return 0;

	const float granularityFrames = granularity * timeToFrame();

	auto safeGetSample = [this]( uint32_t frame )
		{
		if( frame < 0 || frame >= getNumFrames() )
			return 0.0f;
		else
			return getSample( 0, frame );
		};

	std::vector<float> ys( std::ceil( float( getNumFrames() ) / granularityFrames ) );

	for( uint32_t x = 0; x < ys.size(); ++x )
		{
		flan_CANCEL_POINT( Func1x1() );
		const float sampleFrame = x * granularityFrames;
		float sum = 0;
		for( int32_t i = - granularityFrames; i <= granularityFrames; ++i )
			{
			const Sample sample = std::abs( safeGetSample( sampleFrame + i ) );
			const float windowSample = Windows::Hann( ( float( i ) / granularityFrames + 1 ) * .5f );
			// granularityFrames / 2 is the integral of the Hann window
			sum += sample * windowSample * 2.0f / granularityFrames;
			}
		ys[x] = sum;
		}
	
	const float timeToFrame_f = timeToFrame();
	return [ ys, timeToFrame_f, granularityFrames ]( float t )
		{
		float x = t * timeToFrame_f / granularityFrames;
		if( 0 <= x && x < ys.size() - 1 )
			{
			float x1 = std::floor( x );
			float x2 = x1 + 1;
			float y1 = ys[x1];
			float y2 = ys[x2];
			return y1 + ( y2 - y1 ) * ( x - x1 );
			}
		else return 0.0f;
		};
	}