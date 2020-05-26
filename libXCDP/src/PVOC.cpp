#include "xcdp/PVOC.h"

#include <iostream>
#include <complex>
#include <array>

#include "xcdp/spline.h"

#include "xcdp/WindowFunctions.h"
#include "xcdp/Audio.h"
#include "xcdp/CLContext.h"
#include "xcdp/CLProgs.h"

namespace xcdp {

const float pi = std::acos( -1.0f );

//========================================================================
//	Conversions
//========================================================================

/** \cond */
#ifdef USE_FFTW

#include <fftw3.h>
class convertToAudio_FFTHelper
	{
public:
	convertToAudio_FFTHelper( size_t windowSize, size_t numBins )
		: in_( (std::complex<float>*) fftwf_alloc_complex( numBins ) )
		, out_( fftwf_alloc_real( windowSize ) )
		, plan( fftwf_plan_dft_c2r_1d( int(windowSize), (fftwf_complex*) in_, out_, FFTW_MEASURE ) )
		{}
	~convertToAudio_FFTHelper()
		{
		fftwf_destroy_plan( plan );
		fftwf_free( out_ );
		fftwf_free( in_ );
		}

	void execute() { fftwf_execute( plan ); }

	std::complex<float> * inBuffer() { return in_; }

	std::complex<float> in( size_t n ) { return in_[n]; }
	void setIn( size_t n, std::complex<float> v ) { in_[n] = v; }

	float out( size_t n ) { return out_[n]; }
	float setOut( size_t n, float v ) { out_[n] = v; }

private:
	std::complex<float> * in_;
	float * out_;
	fftwf_plan plan;
	};

#else

#define DJ_FFT_IMPLEMENTATION // define this in exactly *one* .cpp file
#include "xcdp/dj_fft.h";

class convertToAudio_FFTHelper
	{
public:
	convertToAudio_FFTHelper( size_t windowSize, size_t numBins ) 
		: in_( windowSize )
		, out_( windowSize, 0 )
		{}

	void execute() 
		{ 
		const std::vector<std::complex<float>> out_c = dj::fft1d( in_, dj::fft_dir::DIR_BWD );
		std::transform( out_c.begin(), out_c.end(), out_.begin(), []( std::complex<float> c ){ return c.real(); } );
		}

	std::complex<float> * inBuffer() { return in_.data(); }

	std::complex<float> in( size_t n ) { return in_[n]; }
	void setIn( size_t n, std::complex<float> v ) 
		{ 
		in_[n] = v; 
		in_[ in_.size() - 1 - n ] = std::conj( v ); 
		}

	float out( size_t n ) { return out_[n]; }
	void setOut( size_t n, float v ) { out_[n] = v; }

private:
	std::vector<std::complex<float>> in_;
	std::vector<float> out_;
	};

#endif

/** \endcond */

Audio convertToAudio_cpu( const PVOC & me )
	{
	std::cout << "Getting audio ... \n";
	const size_t numBins = me.getNumBins();
	const size_t windowSize = me.getWindowSize();
	const size_t hopSize = windowSize / me.getOverlaps();
	const size_t numHops = me.getNumFrames();
	const float binWidth_f = me.binToFrequency();
	
	AudioBuffer::Format audioFormat;
	audioFormat.numChannels = me.getNumChannels();
	audioFormat.numFrames = numHops * hopSize;
	audioFormat.sampleRate = me.getSampleRate();
	Audio out( audioFormat );

	const float windowScale = windowSize * me.getOverlaps();

	std::vector<float> phaseBuffer( numBins );
	convertToAudio_FFTHelper fft( windowSize, numBins );

	for( size_t channel = 0; channel < me.getNumChannels(); ++channel )
		{
		///initial phase?
		std::fill( phaseBuffer.begin(), phaseBuffer.end(), 0 );

		for( size_t hop = 0; hop < numHops; ++hop )
			{
			//reconstruct phase and convert to rectangular, store in fftIn
			//It's the same as in the analysis but backwards so I'm not writing up every step
			for( size_t bin = 0; bin < numBins; ++bin )
				{
				const float magnitude = me.getMF( channel, hop, bin ).m;

				const float trueFreq = me.getMF( channel, hop, bin ).f;
				const float binDeviation = trueFreq / binWidth_f - float(bin);
				const float deltaPhase = 2.0 * pi * binDeviation / float(me.getOverlaps());
				const float expectedPhaseDiff = float(bin) * ( 2.0f * pi / float(me.getOverlaps()) );
				const float phaseDiff = deltaPhase + expectedPhaseDiff;
				phaseBuffer[bin] += phaseDiff;
				
				fft.setIn( bin, std::polar( magnitude, phaseBuffer[bin] ) );
				}

			fft.execute();

			//Accumulate ifft output into audio buffer
			for( size_t relativeSample = 0; relativeSample < windowSize; ++relativeSample )
				{
				// - overlaps/2 is a consequency of marginal frames
				const int actualSample = (int(hop)-int(me.getOverlaps()/2)) * int(hopSize) + int(relativeSample);
				if( 0 <= actualSample && actualSample < out.getNumFrames() )
					{ 
					out.setSample( channel, actualSample, out.getSample( channel, actualSample ) +
						window::Hann( float(relativeSample) / float(windowSize) ) / windowScale * fft.out( relativeSample ) );
					}
				}
			}
		}

	return out;
	}

Audio PVOC::convertToAudio() const
	{

	std::cout << "PVOC::convertToAudio ... \n";

#ifndef USE_OPENCL
	return convertToAudio_cpu( *this );
#endif

	const size_t numBins = getNumBins();
	const int numBins_i = int(numBins);
	const size_t windowSize = getWindowSize();
	const size_t hopSize = windowSize / getOverlaps();
	const size_t numHops = getNumFrames();
	const float binWidth_f = binToFrequency();
	const size_t outChannelDataCount = numBins * numHops;
	const float overlaps_f = float( getOverlaps() );
	const MF * inBufferReadHead = getMFPointer( 0, 0, 0 );
	
	AudioBuffer::Format audioFormat;
	audioFormat.numChannels = getNumChannels();
	audioFormat.numFrames = numHops * hopSize;
	audioFormat.sampleRate = getSampleRate();
	Audio out( audioFormat );

	const float windowScale = float( windowSize * getOverlaps() );

	//Sample Hann window
	std::vector<float> hannWindow( windowSize );
	for( size_t i = 0; i < windowSize; ++i )
		hannWindow[i] = window::Hann( float( i ) / float( windowSize ) ) / windowScale ;

	convertToAudio_FFTHelper fft( windowSize, numBins );

	//prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_convertToAudio );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_READ_WRITE, sizeof( std::complex<float> ) * outChannelDataCount );
	cl::Buffer clPhase( cl.context, CL_MEM_READ_WRITE, sizeof( float ) * numBins );
	cl::KernelFunctor< cl::Buffer > cl_phaseToFFT( programHelper.program, "phaseToFFT" );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, float, float, int > 
		cl_freqToPhase( programHelper.program, "freqToPhase" );

	std::vector<std::complex<float>> clOutHostBuff( outChannelDataCount );
	
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		cl.queue.enqueueFillBuffer( clPhase, 0, 0, sizeof( float ) * numBins );

		//Copy buffer into OpenCL device
		cl::copy( cl.queue, inBufferReadHead, inBufferReadHead + outChannelDataCount, clIn );
		inBufferReadHead += outChannelDataCount;
		cl.queue.enqueueBarrierWithWaitList();

		//Convert frequencies into phases, copy polar into clOut
		//	This must be ordered due to the phase accumulator
		for( size_t hop = 0; hop < numHops; ++hop )
			{
			cl_freqToPhase( cl::EnqueueArgs( cl.queue, cl::NDRange( hop * numBins ), cl::NDRange( numBins ), cl::NDRange() ),
				clIn, 
				clOut, 
				clPhase,
				overlaps_f,
				binWidth_f,
				numBins_i );
			cl.queue.enqueueBarrierWithWaitList();
			}

		//Convert polar to rectangular in-place
		cl_phaseToFFT( cl::EnqueueArgs( cl.queue, cl::NDRange( outChannelDataCount ) ), clOut );
		cl.queue.enqueueBarrierWithWaitList();

		//copy cl out to host mem (faster than copying out of device many times)
		cl::copy( cl.queue, clOut, clOutHostBuff.begin(), clOutHostBuff.end() );
		cl.queue.enqueueBarrierWithWaitList();

		auto clOutHostBuffRead = clOutHostBuff.begin();

		for( size_t hop = 0; hop < numHops; ++hop )
			{
			std::copy( clOutHostBuffRead, clOutHostBuffRead + numBins, fft.inBuffer() );
//#ifndef USE_FFTW //dj_fft can only handle power of two ffts, so we need to recover hermitian symmetry here.
			for( int i = 0; i < numBins - 1; ++i )
				fft.setIn( i, fft.in(i) );
//#endif
			clOutHostBuffRead += numBins;

			fft.execute();

			//Accumulate ifft output into audio buffer
			int actualSample = ( int(hop)-int(getOverlaps()/2) ) * int(hopSize);
			float * inputSampleHead = out.getSamplePointer( channel, actualSample );
			for( size_t relativeSample = 0; relativeSample < windowSize; ++relativeSample, ++actualSample )
				if( 0 <= actualSample && actualSample < out.getNumFrames() )
					*(inputSampleHead + relativeSample) += fft.out( relativeSample ) * hannWindow[relativeSample];
			}
		}

	return out;
	}

const PVOC & PVOC::graph( const std::string & fileName ) const
	{
	std::cout << "PVOC::graph ... \n";

	size_t numFrames = getNumFrames();
	float maxMag = getMaxPartialMagnitude();

	//bin -> frame to match write order
	std::vector<std::array<uint8_t,3>> data( getNumFrames()*getNumBins(), std::array<uint8_t,3>{0,0,0} );
	if( maxMag > 0 )
		for( size_t bin = 0; bin < getNumBins(); ++bin )
			for( size_t frame = 0; frame < numFrames; ++frame )	
				{
				/* The choice of value and hue here is an aesthetic one
				 * I've tried to make both low and high frequency information relatively visible
				 * without oversaturating the highs
				 * Feel free to change the initial hue below if you don't like matrix green
				*/

				const float initialHueAngle = 110;

				const float normMag = std::abs(getMF( 0, frame, bin ).m) / maxMag;
				const float value = std::pow( std::clamp( log2( float(bin) + 3.0f )*normMag, 0.0f, 1.0f ), 0.4f );

				//bin deviation from bin center on range [-.5,.5]
				const float deviation = ( getMF(0, frame, bin).f / float(binToFrequency()) - float(bin)) / float( getOverlaps() );
				const int hueAngle = int( deviation * 1800.0 + initialHueAngle );
				const int wrappedHueAngle = hueAngle - int(360.0*floor( float(hueAngle) / 360.0 ));

				data[frame*getNumBins() + bin] = HSVtoRGB( wrappedHueAngle, 1.0, value );
				}

	writeBMP( fileName, getNumFrames(), data );
	return *this;
	}

//========================================================================
//	Selection
//========================================================================

PVOC PVOC::getFrame( float time ) const
	{
	const float selectedFrame = std::clamp( timeToFrame() * time, 0.0f, float( getNumFrames() - 1 ) );

	auto format = getFormat();
	format.numFrames = 1;
	PVOC out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t bin = 0; bin < out.getNumBins(); ++bin )
			out.setMF( channel, 0, bin, getBinInterpolated( channel, selectedFrame, bin ) );

	return out;
	}

PVOC PVOC::select( float length, Func2x2 selector, Interpolator interp ) const
	{
	std::cout << "PVOC::select ... \n";
	auto format = getFormat();
	format.numFrames = timeToFrame() * length;
	PVOC out( format );

	const float timeToFrameConst = timeToFrame();
	const float frameToTimeConst = frameToTime();
	const float binToFreqConst = binToFrequency();
	const float freqToBinConst = frequencyToBin();

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
			{
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				//No need to clamp, getBinInterpolated does bounds checking
				vec2 selectedPoint = selector( frameToTimeConst * frame, binToFreqConst * bin );
				selectedPoint.x() *= timeToFrameConst;
				selectedPoint.y() *= freqToBinConst;
				out.setMF( channel, frame, bin,  
					getBinInterpolated( channel, selectedPoint.x(), selectedPoint.y(), interp ) );
				}
			}

	return out;
	}

PVOC PVOC::freeze( const std::vector<std::array<float,2>> & timing ) const
	{
	std::cout << "PVOC::freeze ... \n";

	const float timeToFrameConst = timeToFrame();

	//Convert times into frames
	std::vector<std::array<size_t,2>> timingFrames( timing.size() );
	std::transform( timing.begin(), timing.end(), timingFrames.begin(), [timeToFrameConst]( const std::array<float,2> & i )
		{ 
		return std::array<size_t,2>{ 
			size_t( double( i[0] ) * timeToFrameConst ), 
			size_t( double( i[1] ) * timeToFrameConst ) };
		});

	//Sort by occurance frame
	std::sort( timingFrames.begin(), timingFrames.end(), []( std::array<size_t,2> & a, std::array<size_t,2> & b )
		{ return a[0] < b[0]; } );

	float totalFreezeTime = 0;
	for( auto & i : timing ) totalFreezeTime += i[1];

	auto format = getFormat();
	format.numFrames += ceil( totalFreezeTime * timeToFrameConst );
	PVOC out( format );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		bool holding = false;
		size_t numFramesHeld = 0;
		size_t timingIndex = 0;
		for( size_t inFrame = 0, outFrame = 0; inFrame < getNumFrames(); ++outFrame )
			{
			if( holding )
				{
				++numFramesHeld;
				if( numFramesHeld > timingFrames[timingIndex][1] )
					{
					holding = false;
					++timingIndex;
					numFramesHeld = 0;
					}
				}
			else
				{
				++inFrame;
				if( inFrame > timingFrames[timingIndex][0] )
					{
					holding = true;
					}
				}

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				out.setMF( channel, outFrame, bin, getMF( channel, inFrame, bin ) );
			}
		}

	return out;
	}

//========================================================================
//	Resampling
//========================================================================

PVOC PVOC::modify( Func2x2 mod, Interpolator interp ) const
	{
	std::cout << "PVOC::modify ... \n";

	const size_t numBins = getNumBins();
	const size_t numFrames = getNumFrames();
	const size_t inChannelDataCount = numFrames * numBins;

	const float frameToTimeConst = float( getWindowSize() ) / float( getSampleRate() ) / float( getOverlaps() );
	const float timeToFrameConst = 1.0f / frameToTimeConst;
	const float binToFrequencyConst = binToFrequency();
	const float freqToBinConst = 1.0f / binToFrequencyConst;

	//Sample mod functions
	std::vector<MF> inModified( inChannelDataCount );
	std::vector<vec2> modPointSamples( inChannelDataCount );
	float lastOutputFrame = 0;
	for( size_t frame = 0; frame < numFrames; ++frame )
		for( size_t bin = 0; bin < numBins; ++bin )
			{
			auto modOut = mod( frame * frameToTimeConst, bin * binToFrequencyConst );
			modOut.x() *= timeToFrameConst;
			modOut.y() *= freqToBinConst;
			
			modPointSamples[ frame * numBins + bin ] = modOut;
			if( lastOutputFrame < modOut.x() )
				lastOutputFrame = modOut.x();
				
			}
	++lastOutputFrame;
	if( lastOutputFrame * frameToTimeConst > 600.0f ) // outfile longer than 10 minutes?
		{
		std::cout << "PVOC::modify_cpu tried to make a file longer than 10 minutes, which is currently disabled";
		return PVOC();
		}

	auto format = getFormat();
	format.numFrames = ceil( lastOutputFrame );
	PVOC out( format );
	out.clearBuffer();
	const size_t outChannelDataCount = out.getNumFrames() * out.getNumBins();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Calculate and copy mapped frequencies and input magnitudes
		auto modDataSamplesWriteHead = inModified.begin();
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		for( size_t frame = 0; frame < numFrames; ++frame )
			for( size_t bin = 0; bin < numBins; ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				{
				*modDataSamplesWriteHead = 
					{
					inBufferReadHead->m,
					mod( frame * frameToTimeConst, inBufferReadHead->f ).y()
					};
				}
		
		for( size_t frame = 1; frame < numFrames; ++frame )
			for( size_t bin = 1; bin < numBins; ++bin )
				{
				const std::array<vec2, 4> p = {
					modPointSamples[ ( frame - 1 ) * numBins + bin - 1 ],
					modPointSamples[ ( frame - 0 ) * numBins + bin - 1 ],
					modPointSamples[ ( frame - 0 ) * numBins + bin - 0 ],
					modPointSamples[ ( frame - 1 ) * numBins + bin - 0 ] };
		
				const std::array<MF, 4> pMF = {
					inModified[ ( frame - 1 ) * numBins + bin - 1 ],
					inModified[ ( frame - 0 ) * numBins + bin - 1 ],
					inModified[ ( frame - 0 ) * numBins + bin - 0 ],
					inModified[ ( frame - 1 ) * numBins + bin - 0 ] };
		
				const vec2 D12 = p[1] - p[0];
				const vec2 D23 = p[2] - p[1];
				const vec2 D34 = p[3] - p[2];
				const vec2 D41 = p[0] - p[3];
			
				// Find box containing shape to iterate over
				const int minx = fmax( floor( fmin( fmin( p[0].x(), p[1].x() ), fmin( p[2].x(), p[3].x() ) ) ), 0 );
				const int miny = fmax( floor( fmin( fmin( p[0].y(), p[1].y() ), fmin( p[2].y(), p[3].y() ) ) ), 0 );
				const int maxx = fmin( ceil(  fmax( fmax( p[0].x(), p[1].x() ), fmax( p[2].x(), p[3].x() ) ) ), format.numFrames - 1 );
				const int maxy = fmin( ceil(  fmax( fmax( p[0].y(), p[1].y() ), fmax( p[2].y(), p[3].y() ) ) ), numBins - 1 );
			
				// Iterate over box
				for( int x = minx; x <= maxx; ++x )
					for( int y = miny; y <= maxy; ++y )
						{	
						// Ref: http://paulbourke.net/geometry/polygonmesh/#insidepoly
						bool c = false;
						if( ( ( p[0].y() <= y && y < p[3].y() ) || ( p[3].y() <= y && y < p[0].y() ) ) && ( x < D41.x() / D41.y() * ( y - p[0].y() ) + p[0].x() ) ) c = !c;
						if( ( ( p[1].y() <= y && y < p[0].y() ) || ( p[0].y() <= y && y < p[1].y() ) ) && ( x < D12.x() / D12.y() * ( y - p[1].y() ) + p[1].x() ) ) c = !c;
						if( ( ( p[2].y() <= y && y < p[1].y() ) || ( p[1].y() <= y && y < p[2].y() ) ) && ( x < D23.x() / D23.y() * ( y - p[2].y() ) + p[2].x() ) ) c = !c;
						if( ( ( p[3].y() <= y && y < p[2].y() ) || ( p[2].y() <= y && y < p[3].y() ) ) && ( x < D34.x() / D34.y() * ( y - p[3].y() ) + p[3].x() ) ) c = !c;
				
						if( c )
							{
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
									continue;
								m = -quadC / quadB;
								}
							else
								{
								const float descriminant = quadB * quadB - 4.0f * quadA * quadC;
								if( descriminant < 0 ) continue;
								m = ( -quadB + sqrt( descriminant ) ) / ( 2.0f * quadA );
								}
							if( alpha[1] + alpha[3] * m == 0 ) continue;
							const float l = ( x - alpha[0] - alpha[2] * m ) / ( alpha[1] + alpha[3] * m );
					
							const float epsilon = 0.0001f;
							if( fabs( l - 0.5f ) > 0.5f + epsilon || fabs( m - 0.5f ) > 0.5f + epsilon ) continue;

							const std::array<float, 4> w = {
								( 1.0f - interp( l ) ) * ( 1.0f - interp( m ) ) * pMF[0].m,
								(        interp( l ) ) * ( 1.0f - interp( m ) ) * pMF[1].m,
								(        interp( l ) ) * (        interp( m ) ) * pMF[2].m,
								( 1.0f - interp( l ) ) * (        interp( m ) ) * pMF[3].m };
							const float totalWeight = w[0] + w[1] + w[2] + w[3];
							if( totalWeight <= 0 ) continue;
					
							auto & outMF = out.getMF( channel, x, y );

							// Frequency is found by taking a weighted sum of the contributing frequencies of the current step 
							//	and what is already in the output MF
							outMF.f = ( outMF.f * outMF.m 
											  + pMF[0].f * w[0]
											  + pMF[1].f * w[1]
											  + pMF[2].f * w[2]
											  + pMF[3].f * w[3] ) 
											  / ( outMF.m + totalWeight );
							outMF.m += totalWeight;
							}
						}
				}
		}

	return out;
	}

PVOC PVOC::modifyFrequency( Func2x1 outFreqFunc, Interpolator interp ) const
	{
	std::cout << "PVOC::modifyFrequency ... \n";
	PVOC out( getFormat() );
	out.clearBuffer();

	const size_t numChannels = getNumChannels();
	const size_t numFrames = getNumFrames();
	const size_t numBins = getNumBins();

	const float freqToBinConst = frequencyToBin();
	const float binToFreqConst = binToFrequency();
	const float frameToTimeConst = frameToTime();

	for( size_t channel = 0; channel < numChannels; ++channel )
		for( size_t frame = 0; frame < numFrames; ++frame )
			{
			float prevOutBin = freqToBinConst * outFreqFunc( frameToTimeConst * frame, 0.0f );
			for( size_t bin = 1; bin < numBins; ++bin )
				{
				const float outBin = freqToBinConst * outFreqFunc( frameToTimeConst * frame, binToFreqConst * bin );
				const bool forward = outBin > prevOutBin;

				const long long prevOutBinRound = forward? std::ceil( prevOutBin ) : std::floor( prevOutBin );
				const long long outBinRound     = forward? std::ceil( outBin     ) : std::floor( outBin     );
				const size_t startBin = std::clamp( prevOutBinRound, 0ll, long long(out.getNumBins() - 1) );
				const size_t endBin   = std::clamp( outBinRound    , 0ll, long long(out.getNumBins() - 1) );

				for( size_t interpBin = startBin; interpBin != endBin; forward? ++interpBin : --interpBin )
					{
					if( interpBin < 0 || getNumBins() <= interpBin ) continue;

					const MF & currentMF = getMF( channel, frame, bin );
					const MF & prevMF = getMF( channel, frame, bin-1 );

					const float mix = interp( ( float( interpBin ) - prevOutBin ) / ( outBin - prevOutBin ) );
					const float interpMagnitude = ( 1.0f - mix ) * prevMF.m + mix * currentMF.m;

					auto & outMF = out.getMF( channel, frame, interpBin );

					//if( prevMF.m * ( 1.0 - mix ) > currentMF.m * ( 0.0 + mix ) )
					//	outMF.f = outFreqFunc( frameToTime( frame ), prevMF.f );
					//else
					//	outMF.f = outFreqFunc( frameToTime( frame ), currentMF.f );

					const float w0 = ( 1.0f - mix ) * prevMF.m;
					const float w1 = (        mix ) * currentMF.m;

					outMF.f = ( outMF.m * outMF.f + w0 * prevMF.f + w1 * currentMF.f ) / ( outMF.m + w0 + w1 );
					outMF.m += interpMagnitude;
					}
				prevOutBin = outBin;
				}
			}

	return out;
	}

PVOC PVOC::modifyTime( Func2x1 outPosFunc, Interpolator interp ) const
	{
	std::cout << "PVOC::modifyTime ... \n";

	const size_t numChannels = getNumChannels();
	const size_t numFrames = getNumFrames();
	const size_t numBins = getNumBins();

	//const float freqToBinConst = frequencyToBin();
	const float binToFreqConst = binToFrequency();
	const float frameToTimeConst = frameToTime();
	const float timeToFrameConst = timeToFrame();

	//Find farthest output frame and sample mod function
	std::vector<float> modFrames( numFrames * numBins );
	float lastOutputFrame = 0;
	for( size_t frame = 0; frame < getNumFrames(); ++frame )
		for( size_t bin = 0; bin < getNumBins(); ++bin )
			{
			const float modFrame = outPosFunc( frame * frameToTimeConst, bin * binToFreqConst ) * timeToFrameConst;
			if( lastOutputFrame < modFrame )
				lastOutputFrame = modFrame;
			modFrames[ frame * numBins + bin ] = modFrame;
			}

	auto format = getFormat();
	format.numFrames = lastOutputFrame;
	PVOC out( format );
	out.clearBuffer();

	for( size_t channel = 0; channel < numChannels; ++channel )
		{
		for( size_t frame = 1; frame < numFrames; ++frame ) // Start at 1 because for each frame we will access previous frame
			{
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const float prevModFrame = modFrames[ ( frame - 1 ) * numBins + bin ];
				const float modFrame = modFrames[ frame * numBins + bin ];
				const bool forward = modFrame > prevModFrame;

				const long long startFrame = forward? std::ceil( prevModFrame ) : std::floor( prevModFrame );
				const long long endFrame   = forward? std::ceil( modFrame     ) : std::floor( modFrame     );

				for( long long outFrame = startFrame; outFrame != endFrame; forward? ++outFrame : --outFrame )
					{
					if( outFrame < 0 || out.getNumFrames() <= outFrame ) continue;

					const MF & currentMF = getMF( channel, frame, bin );
					const MF & prevMF = getMF( channel, frame - 1, bin );

					const float mix = interp( ( modFrame - prevModFrame ) / ( outFrame - prevModFrame ) );
					const float w0 = ( 1.0f -  mix ) * prevMF.m;
					const float w1 = mix * currentMF.m;

					//change to new freq technique
					auto & outMF = out.getMF( channel, outFrame, bin );
					outMF.f = ( outMF.f * outMF.m + w0 * prevMF.f + w1 * currentMF.f ) / ( outMF.m + w0 + w1 );
					outMF.m += w0 + w1;
					}
				}
			}
		}

	return out;
	}

#ifdef USE_OPENCL

PVOC PVOC::modify_cl( Func2x2 mod, Interpolator interp  ) const
	{
	std::cout << "PVOC::modify_cl ... \n";

	using float2 = std::array<float,2>;

	const size_t numBins = getNumBins();
	const int numBins_i = numBins;
	const size_t numFrames = getNumFrames();
	const size_t inChannelDataCount = numFrames * numBins;

	const float frameToTimeConst = float( getWindowSize() ) / float( getSampleRate() ) / float( getOverlaps() );
	const float timeToFrameConst = 1.0f / frameToTimeConst;
	const float binToFrequencyConst = binToFrequency();
	const float freqToBinConst = 1.0f / binToFrequencyConst;

	//Sample mod functions
	std::vector<float2> inModified( inChannelDataCount );
	std::vector<vec2> modPointSamples( inChannelDataCount );
	float lastOutputFrame = 0;
	for( size_t frame = 0; frame < numFrames; ++frame )
		for( size_t bin = 0; bin < numBins; ++bin )
			{
			auto modOut = mod( frame * frameToTimeConst, bin * binToFrequencyConst );
			modOut.x() *= timeToFrameConst;
			modOut.y() *= freqToBinConst;
			
			modPointSamples[ frame * numBins + bin ] = modOut;
			if( lastOutputFrame < modOut.x() )
				{
				//std::cout << modOut.x() << std::endl;
				lastOutputFrame = modOut.x();
				}
			}
	++lastOutputFrame;
	if( lastOutputFrame * frameToTimeConst > 600.0f ) // outfile longer than 10 minutes?
		{
		std::cout << "PVOC::modify tried to make a file longer than 10 minutes, which is currently disabled";
		return PVOC();
		}

	auto format = getFormat();
	format.numFrames = ceil( lastOutputFrame );
	PVOC out( format );
	//out.clearBuffer();
	const size_t outChannelDataCount = out.getNumFrames() * out.getNumBins();
	const int outNumFrames_i = format.numFrames;

	//Prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_modify );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::Buffer clModPointSamples( cl.context, CL_MEM_READ_ONLY, sizeof( float2 ) * inChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, int, int > 
		cl_modify( programHelper.program, "modify" );

	//Copy stretch samples to device
	cl::copy( cl.queue, modPointSamples.begin(), modPointSamples.end(), clModPointSamples );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Clear cl out buffer ( the entire buffer might not be written to so it must be scrubbed )
		cl.queue.enqueueBarrierWithWaitList(); //Can't clear until out copy from previous iteration is complete
		cl.queue.enqueueFillBuffer( clOut, 0, 0, sizeof( PVOCBuffer::MF ) * outChannelDataCount );

		//Calculate and copy mapped frequencies and input magnitudes
		auto modDataSamplesWriteHead = inModified.begin();
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		for( size_t frame = 0; frame < numFrames; ++frame )
			for( size_t bin = 0; bin < numBins; ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				{
				*modDataSamplesWriteHead = 
					{
					inBufferReadHead->m,
					mod( frame * frameToTimeConst, inBufferReadHead->f ).y()
					};
				}
		cl::copy( cl.queue, inModified.begin(), inModified.end(), clIn );
		
		//Compute on device
		cl.queue.enqueueBarrierWithWaitList();
		cl_modify( cl::EnqueueArgs( cl.queue, cl::NDRange( inChannelDataCount ) ), 
			clIn,
			clModPointSamples, 
			clOut, 
			numBins_i,
			outNumFrames_i );
			
		//copy cl output buffer to out
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut,
			out.getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * outChannelDataCount, 
			out.getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * outChannelDataCount );
		}

	return out;
	}

PVOC PVOC::modifyFrequency_cl( Func2x1 outFreqFunc, Interpolator interp  ) const
	{
	std::cout << "PVOC::modifyFrequency_cl ... \n";
	PVOC out( getFormat() );
	out.clearBuffer();

	const size_t numBins = getNumBins();
	const int numBins_i = numBins;
	const size_t numFrames = getNumFrames();
	const size_t inChannelDataCount = numFrames * numBins;

	const float frameToTimeConst = float( getWindowSize() ) / float( getSampleRate() ) / float( getOverlaps() );
	const float binToFrequencyConst = binToFrequency();
	const float freqToBinConst = 1.0f / binToFrequency();

	//Sample stretch function
	std::vector<float> outBinBuffer( inChannelDataCount );
	std::vector<float> mappedFrequencies( inChannelDataCount );
	for( size_t frame = 0; frame < numFrames; ++frame )
		for( size_t bin = 0; bin < numBins; ++bin )
			outBinBuffer[ frame * numBins + bin ] = outFreqFunc( frame * frameToTimeConst, bin * binToFrequencyConst ) * freqToBinConst;

	//Prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_modifyFrequency );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::Buffer clMappedFrequencies( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * inChannelDataCount );
	cl::Buffer clF( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * inChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int > 
		cl_modifyFrequency( programHelper.program, "modifyFrequency" );

	//Copy stretch samples to device
	cl::copy( cl.queue, outBinBuffer.begin(), outBinBuffer.end(), clF );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Clear cl out buffer ( the entire buffer might not be written to so it must be scrubbed )
		cl.queue.enqueueBarrierWithWaitList(); //Can't clear until out copy from previous iteration is complete
		cl.queue.enqueueFillBuffer( clOut, 0, 0, sizeof( PVOCBuffer::MF ) * inChannelDataCount );

		//copy input buffer to device
		cl::copy( cl.queue, 
			getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * inChannelDataCount, 
			getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * inChannelDataCount, 
			clIn );

		//Calculate mapped frequencies
		auto mappedFrequenciesWriteHead = mappedFrequencies.begin();
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		for( size_t frame = 0; frame < numFrames; ++frame )
			for( size_t bin = 0; bin < numBins; ++bin, ++mappedFrequenciesWriteHead, ++inBufferReadHead )
				*mappedFrequenciesWriteHead = outFreqFunc( frame * frameToTimeConst, inBufferReadHead->f );

		//Copy mapped frequencies
		cl::copy( cl.queue, mappedFrequencies.begin(), mappedFrequencies.end(), clMappedFrequencies );
		
		//Compute on device
		cl.queue.enqueueBarrierWithWaitList();
		cl_modifyFrequency( cl::EnqueueArgs( cl.queue, cl::NDRange( inChannelDataCount ) ), 
			clIn, 
			clMappedFrequencies,
			clF, 
			clOut, 
			numBins_i );
			
		//copy cl output buffer to out
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut,
			out.getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * inChannelDataCount, 
			out.getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * inChannelDataCount );
		}

	return out;
	}

PVOC PVOC::modifyTime_cl( Func2x1 outPosFunc, Interpolator interp  ) const
	{
	std::cout << "PVOC::modifyTime_cl ... \n";

	const size_t numBins = getNumBins();
	const size_t numFrames = getNumFrames();
	const size_t inChannelDataCount = numFrames * numBins;
	const float frameToTimeConst = float( getWindowSize() ) / float( getSampleRate() ) / float( getOverlaps() );
	const float binToFrequencyConst = binToFrequency();
	const float timeToFrameConst = 1.0f / frameToTimeConst;

	//Sample stretch function
	std::vector<float> outPosBuffer( inChannelDataCount );
	float lastOutputTime = 0;
	for( size_t frame = 0; frame < numFrames; ++frame )
		for( size_t bin = 0; bin < numBins; ++bin )
			{
			const float outPos = outPosFunc( frame * frameToTimeConst, bin * binToFrequencyConst ) * timeToFrameConst;
			outPosBuffer[ frame * numBins + bin ] = outPos;
			if( outPos > lastOutputTime )
				lastOutputTime = outPos;
			}
	++lastOutputTime;

	auto format = getFormat();
	format.numFrames = ceil( lastOutputTime );
	PVOC out( format );
	const size_t outChannelDataCount = out.getNumFrames() * out.getNumBins();

	//prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_modifyTime );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::Buffer clF( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * inChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, int, int > 
		cl_modifyTime( programHelper.program, "modifyTime" );

	//stretch function samples only need one copy, they're the same for all channels
	cl::copy( cl.queue, outPosBuffer.begin(), outPosBuffer.end(), clF );

	const int numBins_i = numBins;
	const int numPrevFrames_i = out.getNumFrames();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Clear cl out buffer ( the entire buffer might not be written to )
		cl.queue.enqueueBarrierWithWaitList(); //Can't clear until out copy from previous iteration is complete
		cl.queue.enqueueFillBuffer( clOut, 0, 0, sizeof( PVOCBuffer::MF ) * outChannelDataCount );

		//copy input buffer to device
		cl::copy( cl.queue, 
			getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * inChannelDataCount, 
			getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * inChannelDataCount, 
			clIn );
		
		//Compute on device
		cl.queue.enqueueBarrierWithWaitList();
		cl_modifyTime( cl::EnqueueArgs( cl.queue, cl::NDRange( inChannelDataCount ) ), 
			clIn, 
			clF, 
			clOut, 
			numBins_i,
			numPrevFrames_i );
			
		//copy cl output buffer to out
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut,
			out.getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * outChannelDataCount, 
			out.getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * outChannelDataCount );
		}

	cl.queue.enqueueBarrierWithWaitList();
	return out;
	}

#else

//Dummy functions for disabled opencl
PVOC PVOC::modify_cl( Func2x2 mod, Interpolator interp  ) const { return *this; }
PVOC PVOC::modifyFrequency_cl( Func2x1 outFreqFunc, Interpolator interp  ) const { return *this; }
PVOC PVOC::modifyTime_cl( Func2x1 outPosFunc, Interpolator interp  ) const { return *this; }

#endif

PVOC PVOC::repitch( Func2x1 factor, Interpolator interp ) const
	{
	std::cout << "PVOC::repitch ... \n\t";
	return modifyFrequency_cl( [&factor]( vec2 tf ){ return factor(tf.x(),tf.y())*tf.y(); } );
	}

PVOC PVOC::stretch( Func2x1 factor, Interpolator interp ) const
	{
	std::cout << "PVOC::stretch ... \n\t";
	return modifyTime_cl( [&factor]( vec2 tf ){ return factor(tf.x(),tf.y())*tf.x(); }, interp );
	}

PVOC PVOC::stretch_spline( Func1x1 interpolation ) const
	{
	std::cout << "PVOC::stretch_spline ... \n";

	const float frameToTimeConst = frameToTime();
	const float binToFreqConst = binToFrequency();

	const auto safeInterpolation = [&interpolation, this, frameToTimeConst, binToFreqConst]( size_t frame )
		{
		return std::max( size_t( interpolation( frame * frameToTimeConst ) ), 1ull );
		};
	
	//Set up output and spline x coordinates
	auto format = getFormat();
	format.numFrames = 0;
	std::vector<double> Xs( getNumFrames() );
	for( size_t frame = 0; frame < getNumFrames()-1; ++frame )
		{
		Xs[frame] = format.numFrames;
		format.numFrames += safeInterpolation( frame );
		}
	Xs[getNumFrames()-1] = format.numFrames;
	PVOC out( format );

	//Allocate Y coordinate vectors
	std::vector<double> magnitudeYs( Xs.size() );
	std::vector<double> frequencyYs( Xs.size() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Channel -> bin is correct here, despite non-sequential access to pvoc buffer
		for( size_t bin = 0; bin < getNumBins(); ++bin )
			{
			//For each x coordinate get corresponding frequency and magnitude
			for( size_t frame = 0; frame < getNumFrames(); ++frame )
				{
				magnitudeYs[frame] = getMF( channel, frame, bin ).m;
				frequencyYs[frame] = getMF( channel, frame, bin ).f;
				}

			//Do the splines
			tk::spline magnitudeSpline, frequencySpline;
			magnitudeSpline.set_points(Xs,magnitudeYs);
			frequencySpline.set_points(Xs,frequencyYs);

			//Evaluate splines into out
			for( size_t frame = 0; frame < out.getNumFrames(); ++frame )
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

PVOC::MF PVOC::getBinInterpolated( size_t channel, float frame, float bin, Interpolator i ) const
	{
	frame = std::clamp( frame, 0.0f, float( getNumFrames() - 1 ) );
	bin   = std::clamp( bin,   0.0f, float( getNumBins()   - 1 ) );

	const std::array< MF, 4 > p = 
		{
		getMF( channel, floor( frame ), floor( bin ) ),
		getMF( channel, ceil ( frame ), floor( bin ) ),
		getMF( channel, ceil ( frame ), ceil ( bin ) ),
		getMF( channel, floor( frame ), ceil ( bin ) )
		};

	const float l = i( frame - floor( frame ) );
	const float m = i( bin   - floor( bin   ) );

	return {
		   ( 1.0f - m ) * ( ( 1.0f - l ) * p[0].m + l * p[1].m ) + m * ( ( 1.0f - l ) * p[3].m + l * p[2].m ),
		   ( 1.0f - m ) * ( ( 1.0f - l ) * p[0].f + l * p[1].f ) + m * ( ( 1.0f - l ) * p[3].f + l * p[2].f )
		   };
	}

PVOC::MF PVOC::getBinInterpolated( size_t channel, float frame, size_t bin, Interpolator i ) const
	{
	frame = std::clamp( frame, 0.0f, float( getNumFrames() - 1 ) );

	const MF & l = getMF( channel, floor( frame ), bin );
	const MF & h = getMF( channel, ceil ( frame ), bin );

	const float mix = i( frame - floor( frame ) );

	return { 
		   ( 1.0f - mix ) * l.m + mix * h.m,
		   ( 1.0f - mix ) * l.f + mix * h.f
		   };
	}

PVOC::MF PVOC::getBinInterpolated( size_t channel, size_t frame, float bin, Interpolator i ) const
	{
	bin = std::clamp( bin, 0.0f, float( getNumBins() - 1 ) );

	const auto l = getMF( channel, frame, floor( bin ) );
	const auto h = getMF( channel, frame, ceil ( bin ) );

	const float mix = i( bin - floor( bin ) );

	return  { 
			( 1.0f - mix ) * l.m + mix * h.m,
			( 1.0f - mix ) * l.f + mix * h.f
			};
	}

PVOC PVOC::desample( Func2x1 factor, Interpolator interpolator ) const
	{
	std::cout << "PVOC::desample ... \n\t";
	auto downsample = stretch( [&factor]( vec2 tf ){ return 1.0f / factor(tf.x(),tf.y()); }, interpolator );
	std::cout << "\t";
	return downsample.stretch( factor,                                                      interpolator );
	}

PVOC PVOC::timeExtrapolate( float startTime, float endTime, float extrapolationTime, Interpolator interpolator ) const
	{
	std::cout << "PVOC::timeExtrapolate ... \n";

	//Linear extrapolation via lines matching original PVOC at x_1 and x_2
	const size_t startFrame = timeToFrame() * double( startTime );			//x_1
	const size_t endFrame	= timeToFrame() * double( endTime );			//x_2
	const size_t extFrames	= timeToFrame() * double( extrapolationTime );	//Extrapolated frames to generate

	auto format = getFormat();
	format.numFrames = endFrame + extFrames;
	PVOC out( format );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Copy up to start frame
		for( size_t frame = 0; frame < startFrame; ++frame )
			for( size_t bin = 0; bin < getNumBins(); ++bin )
				out.setMF( channel, frame, bin, getMF( channel, frame, bin ) );
			
		//Interpolate and continue beyond endFrame to extrapolate
		for( size_t frame = startFrame; frame < out.getNumFrames(); ++frame )
			{
			const float mix = interpolator (float( frame - startFrame ) / float( endFrame - startFrame ) ); 

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const auto leftBin  =  getMF( channel, startFrame, bin );
				const auto rightBin =  getMF( channel, endFrame, bin );

				out.setMF( channel, frame, bin, 
					{ 
					std::abs( ( 1.0f - mix ) * leftBin.m + mix * rightBin.m ), 
					          ( 1.0f - mix ) * leftBin.f + mix * rightBin.f 
					} );
				}
			}
		}
	return out;
	}

//========================================================================
// Combinations
//========================================================================

PVOC PVOC::replaceAmplitudes( const PVOC & ampSource, Func2x1 amount ) const
	{
	std::cout << "PVOC::replaceAmplitudes ... \n";

	const float frameToTimeConst = frameToTime();
	const float binToFreqConst = binToFrequency();

	PVOC out( getFormat() );

	if( ampSource.getNumChannels() < getNumChannels()
	 || ampSource.getNumFrames()   < getNumFrames() 
	 || ampSource.getNumBins()	   < getNumBins() )
		out.clearBuffer();

	const size_t numChannels = std::min( ampSource.getNumChannels(), getNumChannels() );
	const size_t numFrames	 = std::min( ampSource.getNumFrames(),   getNumFrames()   );
	const size_t numBins     = std::min( ampSource.getNumBins(),     getNumBins()     );

	for( size_t channel = 0; channel < numChannels; ++channel )
		for( size_t frame = 0; frame < numFrames; ++frame )
			for( size_t bin = 0; bin < numBins; ++bin )
				{
				const auto currentBin = getMF( channel, frame, bin );
				const float currentAmount = amount( frame * frameToTimeConst, bin * binToFreqConst );
				out.setMF( channel, frame, bin,
					{
					float( ampSource.getMF( channel, frame, bin ).m * currentAmount + currentBin.m * ( 1.0f - currentAmount ) ),
					float( currentBin.f )
					} );
				}

	return out;
	}

PVOC PVOC::subtractAmplitudes( const PVOC & other, Func2x1 amount ) const
	{
	std::cout << "PVOC::subtractAmplitudes ... \n";

	const float frameToTimeConst = frameToTime();
	const float binToFreqConst = binToFrequency();

	PVOC out( *this );

	const size_t numChannels = std::min( other.getNumChannels(), getNumChannels() );
	const size_t numFrames	 = std::min( other.getNumFrames(),   getNumFrames()   );
	const size_t numBins     = std::min( other.getNumBins(),     getNumBins()     );

	for( size_t channel = 0; channel < numChannels; ++channel )
		for( size_t frame = 0; frame < numFrames; ++frame )
			for( size_t bin = 0; bin < numBins; ++bin )
				{
				auto & outBin = out.getMF( channel, frame, bin );
				const float currentAmount = amount( frame * frameToTimeConst, bin * binToFreqConst );
				outBin.m = std::abs( outBin.m - other.getMF( channel, frame, bin ).m * currentAmount );
				}

	return out;
	}

//========================================================================
// Uncategorized
//========================================================================

PVOC PVOC::perturb( Func1x1 magAmount, Func1x1 frqAmount, 
					Distribution magPerturber, Distribution frqPerturber ) const
	{
	std::cout << "PVOC::perturb ... \n";

	const float frameToTimeConst = frameToTime();

	PVOC out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			const float currentMagAmount = magAmount( frameToTimeConst * frame );
			const float currentFrqAmount = frqAmount( frameToTimeConst * frame );

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				const auto currentBin = getMF( channel, frame, bin );
				out.setMF( channel, frame, bin,
					{
					currentBin.m + magPerturber() * currentMagAmount / log2( 1.0f + bin ),
					currentBin.f + frqPerturber() * currentFrqAmount
					} );
				}
			}

	return out;
	}

PVOC predicateNLoudestPartials( const PVOC & me, Func1x1 numBins, std::function< bool (size_t, size_t) > predicate )
	{
	const float frameToTimeConst = me.frameToTime();

	PVOC out( me.getFormat() );
	out.clearBuffer();

	auto safeNumBins = [&numBins, me, frameToTimeConst]( size_t frame )
		{
		return size_t( std::clamp( numBins( frameToTimeConst * frame ), 0.0f, float( me.getNumBins() ) ) );
		};

	typedef std::pair<size_t, float> indexVolume;
	std::vector<indexVolume> indexAndVolumes( me.getNumBins() );

	for( size_t channel = 0; channel < me.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < me.getNumFrames(); ++frame )
			{
			for( size_t bin = 0; bin < me.getNumBins(); ++bin )
				{
				indexAndVolumes[bin].first = bin;
				indexAndVolumes[bin].second = me.getMF( channel, frame, bin ).m;
				}
			std::sort( indexAndVolumes.begin(), indexAndVolumes.end(), []( indexVolume& a, indexVolume& b ){ return abs(a.second) > abs(b.second); } );
			for( size_t bin = 0; bin < me.getNumBins(); ++bin )
				{
				const size_t actualBin = indexAndVolumes[bin].first;
				if( predicate( bin, safeNumBins( frame ) ) )
					out.setMF( channel, frame, actualBin, me.getMF( channel, frame, actualBin ) );
				else
					out.setMF( channel, frame, actualBin, { 0.0, me.getMF( channel, frame, actualBin ).f } );
				}
			}

	return out;
	}

PVOC PVOC::retainNLoudestPartials( Func1x1 numBins ) const
	{
	std::cout << "PVOC::retainNLoudestPartials ... \n";
	return predicateNLoudestPartials( *this, numBins, []( size_t a, size_t b ){ return a < b; } );
	}

PVOC PVOC::removeNLoudestPartials( Func1x1 numBins ) const
	{
	std::cout << "PVOC::removeNLoudestPartials ... \n";
	return predicateNLoudestPartials( *this, numBins, []( size_t a, size_t b ){ return a >= b; } );
	}

PVOC PVOC::resonate( Func2x1 decay, float length ) const
	{
	std::cout << "PVOC::resonate ... \n";

	length = std::max( length, float(0.0) );
	auto format = getFormat();
	format.numFrames = std::max( size_t( ceil( timeToFrame() * length ) ), 1ull );
	PVOC out( format );

	//Copy first frame into out
	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t bin = 0; bin < out.getNumBins(); ++bin )
			out.setMF( channel, 0, bin, getMF( channel, 0, bin ) );
	
	const float frameToTimeConst = frameToTime();
	const float secondsPerFrame = frameToTime();
	const float binToFreqConst = binToFrequency();

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 1; frame < out.getNumFrames(); ++frame )
			for( size_t bin = 0; bin < out.getNumBins(); ++bin )
				{
				const float currentDecay = pow( decay( frameToTimeConst * frame, binToFreqConst * bin  ), secondsPerFrame );
				const float decayedAmp = out.getMF( channel, frame - 1, bin ).m * currentDecay;
				if( frame < getNumFrames() && getMF( channel, frame, bin ).m > decayedAmp )
					out.setMF( channel, frame, bin, getMF( channel, frame, bin ) );
				else
					out.setMF( channel, frame, bin, { decayedAmp, out.getMF( channel, frame - 1, bin ).f } );
				}

	return out;
	}

} //End namespace xcdp