#include "xcdp/PVOC.h"

#include "xcdp/WindowFunctions.h"
#include "xcdp/Audio.h"
#include "xcdp/CLContext.h"
#include "xcdp/CLProgs.h"
#include "xcdp/FFTHelper.h"
#include "xcdp/Graph.h"

using namespace xcdp;

static const float pi = std::acos( -1.0f );
static const float pi2 = pi * 2.0f;


Audio PVOC::convertToAudio_cpu( std::shared_ptr<FFTHelper> fft, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_FUNCTION_LOG;

	const Bin numBins = getNumBins();
	const Frame dftSize = getDFTSize();
	const Frame windowSize = getWindowSize();
	const Frame hopSize = getHopSize();
	const int numHops = getNumFrames();
	const Frequency binWidth = binToFrequency();
	
	AudioBuffer::Format audioFormat;
	audioFormat.numChannels = getNumChannels();
	audioFormat.numFrames = numHops * hopSize;
	audioFormat.sampleRate = getSampleRate();
	Audio out( audioFormat );

	//Sample Hann window
	std::vector<float> hannWindow( windowSize );
	const float windowScale = 2.67f / ( dftSize * windowSize / hopSize ); // I don't know why, converting audio->pvoc->audio reduces volume by an input specific amount, usually about 2.67
	for( Frame i = 0; i < windowSize; ++i )
		hannWindow[i] = Windows::Hann( float( i ) / float( windowSize - 1 ) ) * windowScale;

	std::vector<float> phaseBuffer( numBins );
	if( ! fft ) fft = std::make_shared<FFTHelper>( dftSize, false, true, false );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		// Initial phase? Using 0, but maybe something else would work better.
		std::fill( phaseBuffer.begin(), phaseBuffer.end(), 0 );

		for( int hop = 0; hop < numHops; ++hop )
			{
			XCDP_CANCEL_POINT( Audio() );
			//reconstruct phase and convert to rectangular, store in fftIn
			//It's the same as in the analysis but backwards so I'm not writing up every step
			for( Bin bin = 0; bin < numBins; ++bin )
				{
				const float bin_f = float( bin );

				const float magnitude = getMF( channel, hop, bin ).m;
				const float frequency = getMF( channel, hop, bin ).f;
				const float binDeviation = frequency / binWidth - bin_f;
				//const float deltaPhase = binDeviation * pi2 * hopSize / dftSize;
				//const float expectedPhaseDiff = bin_f * pi2 * hopSize / dftSize;
				const float phaseDiff = ( binDeviation + bin_f ) * pi2 * hopSize / dftSize;
				phaseBuffer[bin] += phaseDiff;
				
				fft->complexBuffer[bin] = std::polar( magnitude, phaseBuffer[bin] );
				}

			fft->c2rExecute();

			//Accumulate ifft output into audio buffer
			const Frame outFrameStart = hopSize * hop - windowSize / 2;
			const Frame outFrameEnd = outFrameStart + windowSize;
			const Frame outFrameStart_s = std::max( outFrameStart, 0 );
			const Frame outFrameEnd_s   = std::min( outFrameEnd, out.getNumFrames() );

			const Frame fftStart = outFrameStart_s - outFrameStart;
			const Frame fftEnd = outFrameEnd_s - outFrameStart;

			float * writeHead = out.getSamplePointer( channel, outFrameStart_s );
			for( Frame fftFrame = fftStart; fftFrame < fftEnd; ++fftFrame )
				{
				*writeHead += fft->realBuffer[fftFrame] * hannWindow[fftFrame];
				++writeHead;
				}
			}
		}

	return out;
	}

Audio PVOC::convertToAudio( std::shared_ptr<FFTHelper> fft, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( Audio() );

#ifndef USE_OPENCL
	return convertToAudio_cpu( *this, fft, canceller );
#else

	const Bin numBins = getNumBins();
	const Frame dftSize = getDFTSize();
	const Frame windowSize = getWindowSize();
	const Frame hopSize = getHopSize();
	const uint32_t numHops = getNumFrames();
	const Frequency binWidth = binToFrequency();
	const float overlaps = float( dftSize ) / hopSize;
	const uint32_t outChannelDataCount = numBins * numHops;
	const MF * inBufferReadHead = getMFPointer( 0, 0, 0 );
	
	AudioBuffer::Format audioFormat;
	audioFormat.numChannels = getNumChannels();
	audioFormat.numFrames = numHops * hopSize;
	audioFormat.sampleRate = getSampleRate();
	Audio out( audioFormat );

	//Sample Hann window
	std::vector<float> hannWindow( windowSize );
	const float windowScale = 2.67f / ( dftSize * windowSize / hopSize ); // I don't know why, converting audio->pvoc->audio reduces volume by an input specific amount, usually about 2.67
	for( uint32_t i = 0; i < windowSize; ++i )
		hannWindow[i] = Windows::Hann( float( i ) / float( windowSize - 1 ) ) * windowScale;

	if( ! fft ) fft = std::make_shared<FFTHelper>( dftSize, false, true, false );

	//prepare OpenCL
	CLContext cl;
	ProgramHelper programHelper( cl, CLProgs::PVOC_convertToAudio );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_READ_WRITE, sizeof( std::complex<float> ) * outChannelDataCount );
	cl::Buffer clPhase( cl.context, CL_MEM_READ_WRITE, sizeof( float ) * numBins );
	cl::KernelFunctor< cl::Buffer > cl_phaseToFFT( programHelper.program, "phaseToFFT" );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, float, float, int > 
		cl_freqToPhase( programHelper.program, "freqToPhase" );

	std::vector<std::complex<float>> clOutHostBuff( outChannelDataCount );
	
	for( uint32_t channel = 0; channel < getNumChannels(); ++channel )
		{
		cl.queue.enqueueFillBuffer( clPhase, 0, 0, sizeof( float ) * numBins );

		// Copy channel into OpenCL device
		cl::copy( cl.queue, inBufferReadHead, inBufferReadHead + outChannelDataCount, clIn );
		inBufferReadHead += outChannelDataCount;
		XCDP_CANCEL_POINT( Audio() );
		cl.queue.enqueueBarrierWithWaitList();

		// Convert frequencies into phases, copy polar into clOut
		//	This must be ordered due to the phase accumulator
		for( uint32_t hop = 0; hop < numHops; ++hop )
			{
			cl_freqToPhase( cl::EnqueueArgs( cl.queue, cl::NDRange( hop * numBins ), cl::NDRange( numBins ), cl::NDRange() ),
				clIn, 
				clOut, 
				clPhase,
				overlaps,
				binWidth,
				numBins );
			XCDP_CANCEL_POINT( Audio() );
			cl.queue.enqueueBarrierWithWaitList();
			}

		//Convert polar to rectangular in-place
		cl_phaseToFFT( cl::EnqueueArgs( cl.queue, cl::NDRange( outChannelDataCount ) ), clOut );
		XCDP_CANCEL_POINT( Audio() );
		cl.queue.enqueueBarrierWithWaitList();

		// Copy cl out to host mem (faster than copying out of device many times)
		cl::copy( cl.queue, clOut, clOutHostBuff.begin(), clOutHostBuff.end() );
		XCDP_CANCEL_POINT( Audio() );
		cl.queue.enqueueBarrierWithWaitList();

		auto clOutHostBuffRead = clOutHostBuff.begin();

		for( uint32_t hop = 0; hop < numHops; ++hop )
			{
			XCDP_CANCEL_POINT( Audio() );

			std::copy( clOutHostBuffRead, clOutHostBuffRead + numBins, fft->complexBuffer );
			clOutHostBuffRead += numBins;

			fft->c2rExecute();

			// Accumulate ifft output into audio buffer
			const Frame outFrameStart = hopSize * hop - windowSize / 2;
			const Frame outFrameEnd = outFrameStart + windowSize;
			const Frame outFrameStart_s = std::max( outFrameStart, 0 );
			const Frame outFrameEnd_s   = std::min( outFrameEnd, out.getNumFrames() );

			// Where in the fft buffer we are reading from
			const Frame fftStart = outFrameStart_s - outFrameStart;
			const Frame fftEnd = outFrameEnd_s - outFrameStart;

			float * writeHead = out.getSamplePointer( channel, outFrameStart_s );
			for( Frame fftFrame = fftStart; fftFrame < fftEnd; ++fftFrame )
				{
				*writeHead += fft->realBuffer[fftFrame] * hannWindow[fftFrame];
				++writeHead;
				}
			}
		}

	return out;
#endif
	}
	
Graph PVOC::convertToGraph( Rect D, Pixel width, Pixel height, float timelineScale, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( Graph( width, height ) );

	if( D.x2() == -1 ) D.d.x2 = getLength();
	if( D.y2() == -1 ) D.r.x2 = getHeight();

	// Validation, Conversion
	const float startFrame	= std::clamp( int( D.x1() * timeToFrame() ), 0, getNumFrames() - 1 );
	const float endFrame	= std::clamp( int( D.x2() * timeToFrame() ), 0, getNumFrames() - 1 );
	const float startBin	= std::clamp( int( D.y1() * frequencyToBin() ), 0, getNumBins() - 1 );
	const float endBin		= std::clamp( int( D.y2() * frequencyToBin() ), 0, getNumBins() - 1 );

	// Convert PVOC data to Value/Hue data
	const float maxMag = getMaxPartialMagnitude( startFrame, endFrame, startBin, endBin );
	std::vector<Func2x1> fs;
	for( int i = 0; i < getNumChannels(); ++i )
		fs.push_back( [this, maxMag, i]( float x, float y )
			{
			const float m = getMF( i, x * timeToFrame(), y * frequencyToBin() ).m;
			return std::sqrt( std::abs( m ) / maxMag ) * log2( 2.0f + y ) / 4.0f; // sqrt brings up dark areas, log scaling brings up high frequencies
			} );

	Graph g( width, height );
	g.addFullSplitViewY( D, getNumChannels() ); 
	if( maxMag != 0 )
		g.drawSpectrograms( fs, { 0, 0, getLength(), getHeight() }, canceller );
	
	// Tick drawing
	if( timelineScale > 0 )
		{
		const float bigTimeJump = std::pow( 4.0f, std::floor( std::log2( D.w() ) / 2 - 0.5f ) );
		g.drawXTicks( bigTimeJump / 4.0f, D.y2(), timelineScale / 2, 0, -1, Color::fromHSV( 0, 0, .6 ), 0 );
		g.drawXTicks( bigTimeJump, D.y2(), timelineScale, 0, -1, Color::fromHSV( 0, 0, 1 ), timelineScale );
		//g.drawYTicks( 2000, 0, 0, 14, -1, Color::fromHSV( 0, 0, 1 ), true );
		}

	return g;
	}

const PVOC & PVOC::saveToBMP( const std::string & fileName, Rect D, Pixel width, Pixel height, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_FUNCTION_LOG;
	convertToGraph( D, width, height, canceller ).save_image( fileName );
	return *this;
	}