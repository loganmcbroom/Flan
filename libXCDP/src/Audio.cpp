#include "xcdp/Audio.h"

#include <iostream>
#include <algorithm>
#include <complex>

#include <fftw3.h>

#include "xcdp/Function.h"
#include "xcdp/PVOC.h"
#include "xcdp/Spectrum.h"
#include "xcdp/WindowFunctions.h"
#include "WDL/resample.h"
#include "xcdp/CLContext.h"
#include "xcdp/CLProgs.h"

using namespace xcdp;

const float pi = std::acos( -1.0f );

Audio::Audio() 
	: AudioBuffer() 
	{}

Audio::Audio( const AudioBuffer::Format & other ) 
	: AudioBuffer( other ) 
	{}

Audio::Audio( const std::string & filename ) 
	: AudioBuffer( filename ) 
	{}

//========================================================
// Helper Functions
//========================================================

void printToLog( std::string s )
	{
	std::cout << s << std::endl;
	}
bool doSampleRatesMatch( Audio::Vec ins )
	{
	//Check if all sample rates match the first file
	size_t sampleRate = ins[0].getSampleRate();
	for( auto & in : ins ) 
		if( in.getSampleRate() != sampleRate )
			{
			printToLog( "Mismatched sample rates" );
			return false;
			}
	return true;
	}
bool doChannelCountsMatch( Audio::Vec ins )
	{
	if( ins.size() == 0 ) return true;

	//Check if all channel counts match the first one
	size_t numchannels = ins[0].getNumChannels();
	for( auto & in : ins ) 
		if( in.getNumChannels() != numchannels )
			{
			printToLog( "Mismatched channel count" );
			return false;
			}
	return true;
	}

size_t getMaxNumChannels( Audio::Vec ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumChannels() < b.getNumChannels();
		} )->getNumChannels();
	}
size_t getMaxNumSamples( Audio::Vec ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumFrames() < b.getNumFrames();
		} )->getNumFrames();
	}

//========================================================
// Overloads
//========================================================

Audio Audio::operator+( const Audio & other ) const
	{
	return mix( {*this, other} );
	}

Audio Audio::operator-() const	
	{
	return invertPhase();
	}

//========================================================
// Information
//========================================================

float Audio::getTotalEnergy() const
    {
    double sum = 0;
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumFrames(); ++sample )
            sum += pow( getSample( channel, sample ), 2 );
    
    return sum;
    }

Audio Audio::graph( const std::string & filename, size_t width, size_t height ) const
	{
	std::cout << "Audio::graph ... " << std::endl;
	const auto audioColor = HSVtoRGB( 0, .8f, .65f ); // redish
	const auto backgroundColor = HSVtoRGB( 180, .8f, .1f ); // blueish

	if( getNumFrames() == 0 ) return *this;
	width = std::min( width, getNumFrames() );

	std::vector<std::vector<std::array<uint8_t,3>>> data( width, std::vector( height*getNumChannels(), backgroundColor ) );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		std::vector<float> energies( width, 0.0 );
		for( size_t sample = 0; sample < getNumFrames(); ++sample )
			energies[ float( sample ) / float( getNumFrames() * width ) ] += getSample( channel, sample );
		float maxEnergy = 0;
		for( size_t x = 0; x < width; ++x )
			maxEnergy = std::max( maxEnergy, std::abs( energies[x] ) );

		for( size_t x = 0; x < width; ++x )
			{
			const int height2 = height / 2;
			data[x][height2 + height * channel ] = {0,0,0};
			for( long long int y = 0; std::abs( y ) < std::abs( energies[x] / maxEnergy * height2 ) * 0.9; energies[x] > 0 ? ++y : --y )
				{
				data[x][ height2 + y + height * channel ] = audioColor;
				}
			}
		}

	writeBMP( filename, data );

	return *this;
	}

//========================================================
// Conversions
//========================================================

PVOC convertToPVOC_cpu( const Audio & me, const size_t frameSize, const size_t overlaps )
	{
	//Figure out some helpful quantities
	const size_t numBins = frameSize / 2 + 1;
	const size_t hopSize = frameSize / overlaps;
	const size_t numFrames = me.getNumFrames();
	//+1 since we analyze at start and end times
	const size_t numHops = size_t( ceil( numFrames / hopSize ) ) + 1; 
	const float binWidthD = float( me.getSampleRate() ) / float( frameSize );
	
	PVOCBuffer::Format PVOCFormat;
	PVOCFormat.numChannels = me.getNumChannels();
	PVOCFormat.numFrames = numHops;
	PVOCFormat.numBins = numBins;
	PVOCFormat.sampleRate = me.getSampleRate();
	PVOCFormat.overlaps = overlaps;
	PVOC out( PVOCFormat );

	//Allocate fftw input buffer, output buffer, phase buffer, and plan
	std::vector<float> phaseBuffer( numBins );
	float * fftIn = fftwf_alloc_real( frameSize );
	std::complex<float> *fftOut = (std::complex<float>*) fftwf_alloc_complex( numBins );
	fftwf_plan plan = fftwf_plan_dft_r2c_1d( int( frameSize ), fftIn, (fftwf_complex*) fftOut, FFTW_MEASURE );

	//For each channel, do the whole thing
	for( size_t channel = 0; channel < me.getNumChannels(); ++channel )
		{
		//Set initial phase to 0
		for( size_t bin = 0; bin < numBins; ++bin )
			phaseBuffer[bin] = 0;

		//For each hop, fft and save into buffer
		for( size_t hop = 0; hop < numHops; ++hop )
			{
			//Fill fft input buffer with windowed signal, condition protects from bad buffer access in the last few hops
			for( size_t relativeSample = 0; relativeSample < frameSize; ++relativeSample )
				{
				// - overlaps/2 to center window around analysis point
				const int actualSample = int(relativeSample) + int(hopSize) * (int(hop) - int(overlaps/2));
				if( 0 <= actualSample && actualSample < numFrames )
					fftIn[relativeSample] = me.getSample( channel, actualSample ) * window::Hann( float(relativeSample) / float(frameSize) );
				else
					fftIn[relativeSample] = 0;
				}
	
			fftwf_execute( plan );

			//For each bin, compute frequency and magnitude
			for( size_t bin = 0; bin < numBins; ++bin )
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
				const float phase = arg( fftOut[bin] );
				const float phaseDiff = phase - phaseBuffer[bin];
				phaseBuffer[bin] = phase;
				const float expectedPhaseDiff = float(bin) * ( 2.0 * pi / float(overlaps) );
				const float deltaPhase = phaseDiff - expectedPhaseDiff;
				const float wrappedDeltaPhase = deltaPhase-(2.0*pi)*floor((deltaPhase+pi) / ( 2.0 * pi ) );
				const float binDeviation = float(overlaps) * wrappedDeltaPhase / ( 2.0 * pi );
				const float frequency = ( float(bin) + binDeviation ) * binWidthD;

				out.setMF( channel, hop, bin, { abs( fftOut[bin] ), frequency } );
				}
			}
		}

	fftwf_destroy_plan( plan );
	fftwf_free( fftOut );
	fftwf_free( fftIn );

	return out;
	}

PVOC Audio::convertToPVOC( const size_t frameSize, const size_t overlaps ) const
	{
	struct convertToPVOC_FFTWHelper 
		{
		convertToPVOC_FFTWHelper( size_t frameSize, size_t numBins )
			: in( fftwf_alloc_real( frameSize ) )
			, out( (std::complex<float>*) fftwf_alloc_complex( numBins ) )
			, plan( fftwf_plan_dft_r2c_1d( int( frameSize ), in, (fftwf_complex*) out, FFTW_MEASURE ) )
			{}
		~convertToPVOC_FFTWHelper()
			{
			fftwf_destroy_plan( plan );
			fftwf_free( out );
			fftwf_free( in );
			}

		void execute() { fftwf_execute( plan ); }

		float * in;
		std::complex<float> * out;
		fftwf_plan plan;
		};

	std::cout << "Audio::convertToPVOC ... " << std::endl;

#ifndef USE_OPENCL
	return convertToPVOC_cpu( *this, frameSize, overlaps );
#else

	//Figure out some helpful quantities
	const size_t numBins = frameSize / 2 + 1;
	const int numBins_i = numBins;
	const float overlaps_f = overlaps;
	const size_t hopSize = frameSize / overlaps;
	const size_t numFrames = getNumFrames();
	//+1 since we analyze at start and end times
	const size_t numHops = size_t( ceil( numFrames / hopSize ) ) + 1; 
	const float binWidth_f = float( getSampleRate() ) / float( frameSize );
	const size_t outChannelDataCount = numBins * numHops;
	
	PVOCBuffer::Format PVOCFormat;
	PVOCFormat.numChannels = getNumChannels();
	PVOCFormat.numFrames = numHops;
	PVOCFormat.numBins = numBins;
	PVOCFormat.sampleRate = getSampleRate();
	PVOCFormat.overlaps = overlaps;
	PVOC out( PVOCFormat );
	PVOCBuffer::MF * outBufferWriteHead = out.getMFPointer( 0, 0, 0 );

	//Sample Hann window
	std::vector<float> hannWindow( frameSize );
	for( size_t i = 0; i < frameSize; ++i )
		hannWindow[i] = window::Hann( float( i ) / float( frameSize ) );

	//Prepare fftw
	convertToPVOC_FFTWHelper fftw( frameSize, numBins );
	
	//prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::Audio_convertToPVOC );
	cl::Buffer clIn( cl.context, CL_MEM_READ_WRITE, sizeof( std::complex<float> ) * outChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::KernelFunctor< cl::Buffer > cl_fftToPhase( programHelper.program, "fftToPhase" );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, float, float, int > 
		cl_phaseToFreq( programHelper.program, "phaseToFreq" );
	
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		for( size_t hop = 0; hop < numHops; ++hop )
			{
			//Fill fft input buffer with windowed signal
			int actualSample = int(hopSize) * (int(hop) - int(overlaps/2));
			const float * inputSampleHead = getSamplePointer( channel, actualSample );
			for( size_t relativeSample = 0; relativeSample < frameSize; ++relativeSample, ++actualSample )
				{
				if( 0 <= actualSample && actualSample < numFrames )
					fftw.in[relativeSample] = *(inputSampleHead + relativeSample) * hannWindow[relativeSample];
				else
					fftw.in[relativeSample] = 0;
				}

			cl.queue.enqueueBarrierWithWaitList(); //Can't start writing to fftOut until copying to device from previous iteration completes
			fftw.execute();

			//Copy fft output to device memory
			
			cl.queue.enqueueWriteBuffer( clIn, false, 
				sizeof( std::complex<float> ) * hop * numBins, 
				sizeof( std::complex<float> ) * numBins, 
				fftw.out );
			}

		//transform spectral buffers into phase/mag buffers
		cl.queue.enqueueBarrierWithWaitList();
		cl_fftToPhase( cl::EnqueueArgs( cl.queue, cl::NDRange( outChannelDataCount ) ), clIn );

		//second call computes all frequencies ( and magnitudes to avoid 2 copies from device ) into second buffer
		cl.queue.enqueueBarrierWithWaitList();
		cl_phaseToFreq( cl::EnqueueArgs( cl.queue, cl::NDRange( outChannelDataCount ) ), 
			clIn, 
			clOut,
			overlaps_f,
			binWidth_f,
			numBins_i );

		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut, outBufferWriteHead, outBufferWriteHead + outChannelDataCount );
		outBufferWriteHead += outChannelDataCount;
		};

	cl.queue.enqueueBarrierWithWaitList();
	return out;
#endif
	}

//Spectrum Audio::convertToSpectrum() const
//	{
//	const size_t outnumBins = std::floor( float( getNumFrames() ) / 2.0f ) + 1;
//
//	Spectrum::Format format;
//	format.numChannels = getNumChannels();
//	format.numBins = outnumBins;
//	format.sampleRate = getSampleRate();
//	Spectrum out( format );

	//Allocate fftw buffers and plan
	//real_t * fftIn = fftw_alloc_real( getNumFrames() );
	//std::complex<real_t> * fftOut = (std::complex<real_t>*) fftw_alloc_complex( outnumBins );
	//fftw_plan plan = fftw_plan_dft_r2c_1d( int( getNumFrames() ), fftIn, (fftw_complex*) fftOut, FFTW_ESTIMATE );

	//for( size_t channel = 0; channel < getNumChannels(); ++channel )
	//	{
	//	//Read buffer data into fftIn
	//	for( size_t sample = 0; sample < getNumFrames(); ++sample )
	//		fftIn[sample] = getSample( channel, sample );

	//	//fft
	//	fftw_execute( plan );

	//	//Read fftOut into out spectrum
	//	for( size_t bin = 0; bin < outnumBins; ++bin )
	//		out.setSpectra( channel, bin, fftOut[bin] );
	//	}

	//fftw_destroy_plan( plan );
	//fftw_free( fftOut );
	//fftw_free( fftIn );

	//return out;
	//}

Audio Audio::convertToMidSide() const
	{
	std::cout << "Audio::convertToMidSide ... " << std::endl;
	if( getNumChannels() != 2 )
		{
		printToLog( "Can't transform non-stereo Audio between Mid-Side and Left-Right formats. \n" );
		return *this;
		}
	Audio out( getFormat() );
	for( size_t sample = 0; sample < getNumFrames(); ++sample )
		{
		out.setSample( 0, sample, getSample( 0, sample ) + getSample( 1, sample ) );
		out.setSample( 1, sample, getSample( 0, sample ) - getSample( 1, sample ) );
		}
	return out;
	}

Audio Audio::convertToLeftRight() const
	{
	std::cout << "Audio::convertToLeftRight ... " << std::endl << "\t";
	return convertToMidSide();
	}

Audio Audio::convertToStereo() const
	{
	std::cout << "Audio::convertToStereo ... " << std::endl;
	auto format = getFormat();
	format.numChannels = 2;
	Audio out( format );

	switch( getNumChannels() )
		{
		case 1:
			{
			for( size_t sample = 0; sample < out.getNumFrames(); ++sample )
				{
				out.setSample( 0, sample, getSample( 0, sample ) / sqrt(2.0f) );
				out.setSample( 1, sample, getSample( 0, sample ) / sqrt(2.0f) );
				}
			break;
			}
		case 2:
			{
			out = *this;
			break;
			}
		default:
			{
			printToLog( "I don't know how to convert that number of channels to stereo." );
			}
		}
	return out;
	}

Audio Audio::convertToMono() const
	{
	std::cout << "Audio::convertToMono ... " << std::endl;
	auto format = getFormat();
	format.numChannels = 1;
	Audio out( format );

	for( size_t sample = 0; sample < out.getNumFrames(); ++sample )
		{
		float sampleAccumulator = 0;
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			sampleAccumulator += getSample( channel, sample );
		out.setSample( 0, sample, sampleAccumulator / sqrt( getNumChannels() ) );
		}

	return out;
	}

//========================================================
// Procs
//========================================================

Audio Audio::invertPhase() const
	{
	std::cout << "Audio::invertPhase ... " << std::endl;
	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumFrames(); ++sample )
			out.setSample( channel, sample, -getSample( channel, sample ) );

	return out;
	}

Audio Audio::modifyVolume( Func1x1 volumeLevel ) const
	{
	std::cout << "Audio::modifyVolume ... " << std::endl;

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumFrames(); ++sample )
			{
			float calculatedSample = getSample( channel, sample ) * volumeLevel( frameToTime() * sample );
			out.setSample( channel, sample, calculatedSample );
			}

	return out;
	}

Audio Audio::setVolume( Func1x1 level ) const
	{
	std::cout << "Audio::setVolume ... " << std::endl << "\t";
	// Divide by getMaxSampleMagnitude to normalize, multiply by level to set
	const float maxMag = getMaxSampleMagnitude();
	if( maxMag == 0 ) return *this;
	return modifyVolume( [ &level, maxMag ]( float t ){ return level(t) / maxMag; } );
	}

Audio Audio::waveshape( Func1x1 shaper ) const
	{
	std::cout << "Audio::waveshape ... " << std::endl;

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumFrames(); ++sample )
			{
			out.setSample( channel, sample, shaper( getSample( channel, sample ) ) );
			}

	return out;
	}

Audio Audio::pan( Func1x1 panAmount ) const
	{
	std::cout << "Audio::pan ... " << std::endl;

	//Stereo panning algorithm
	auto stereoPan = [this]( Func1x1 panAmount )
		{
		Audio out( getFormat() );

		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			for( size_t frame = 0; frame < getNumFrames(); ++frame )
				{
				float sample = getSample( channel, frame )
					* sin( pi / 4.0f * ( std::clamp( panAmount( frameToTime() * frame ), -1.0f, 1.0f ) + 3.0f - float(channel)*2.0f ) )
					* sqrt( 2.0f );
				out.setSample( channel, frame, sample );
				}
		return out;
		};
	
	switch( getNumChannels() )
		{
		case 1: //Panning mono: cast to stereo and do stereo pan
			return convertToStereo().pan( panAmount );
			break;
		case 2: //Panning stereo: Use stereo pan as defined above
			return stereoPan( panAmount );
			break;
		default:
			printToLog( "I don't know how to pan that number of channels" );
			return *this;
		}
	}

Audio Audio::widen( Func1x1 widenAmount ) const
	{
	return convertToMidSide().pan( widenAmount ).convertToLeftRight();
	}

Audio Audio::iterate( size_t n, Audio::Mod mod, bool fbIterate ) const
	{
	std::cout << "Audio::iterate ... " << std::endl;

	if( mod == nullptr )
		{
		auto format = getFormat();
		format.numFrames = getNumFrames() * n;
		Audio out( format );
			
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			size_t outSample = 0;
			for( size_t i = 0; i < n; ++i )
				{
				for( size_t inSample = 0; inSample < getNumFrames(); ++inSample )
					{
					out.setSample( channel, outSample, getSample( channel, inSample ) );
					++outSample;
					}
				}
			}
		return out;
		}
	else
		{
		//For each iteration, apply mod, either to input, or previous mod output based on fbIterate
		std::vector<Audio> modOutputs;
		size_t totalSamplesGenerated = 0;
		modOutputs.push_back( mod( *this, 0 ) );
		totalSamplesGenerated += modOutputs[0].getNumFrames();
		for( size_t i = 1; i < n; ++i )
			{
			modOutputs.push_back( mod( fbIterate? modOutputs[i-1] : *this, i ) );
			totalSamplesGenerated += modOutputs[i].getNumFrames();
			}
	
		auto format = getFormat();
		format.numFrames = totalSamplesGenerated;
		Audio out( format );
			
		//Append mod outputs to out
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			size_t outSample = 0;
			for( size_t i = 0; i < n; ++i )
				{
				for( size_t modSample = 0; modSample < modOutputs[i].getNumFrames(); ++modSample )
					{
					out.setSample( channel, outSample, modOutputs[i].getSample( channel, modSample ) );
					++outSample;
					}
				}
			}
		return out;
		}
	}

Audio Audio::reverse() const
	{
	std::cout << "Audio::reverse ... " << std::endl;

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumFrames(); ++sample )
			{
			out.setSample( channel, sample, getSample( channel, getNumFrames() - 1 - sample ) );
			}

	return out;
	}

Audio Audio::cut( float startTime, float endTime ) const
	{
	std::cout << "Audio::cut ... " << std::endl;
	//input validity checking
	if( endTime < startTime ) endTime = startTime;
	const size_t startSample = std::max( size_t(0),      size_t( ceil ( timeToFrame() * startTime ) ) );
	const size_t endSample   = std::min( getNumFrames(), size_t( ceil ( timeToFrame() * endTime   ) ) );

	auto format = getFormat();
	format.numFrames =  endSample - startSample;
	Audio out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t sample = 0; sample < out.getNumFrames(); ++sample )
			out.setSample( channel, sample, getSample( channel, startSample + sample ) );

	return out;
	}

Audio Audio::repitch( Func1x1 factor, size_t granul, size_t qual ) const
	{
	std::cout << "Audio::repitch ... " << std::endl;

	auto safeFactor = [this, &factor]( float count )
		{ 
		return std::max( 1.0f / factor( float(count) / getSampleRate() ), 0.001f ); 
		};

	// estimate output length
	float acclen = 0.0;
	for( int count = 0; count < getNumFrames(); count += granul )
		acclen += granul * safeFactor( count );

	auto format = getFormat();
	format.numFrames = size_t(ceil(acclen));
	Audio out( format );
	
	WDL_Resampler rs;
		 if( qual == 0 ) rs.SetMode( true,  0, true, 64 ); // reasonable sinc resampling
	else if( qual == 1 ) rs.SetMode( true,  1, false    ); // linear interpolation with simple filter
	else if( qual == 2 ) rs.SetMode( false, 0, false    ); // no interpolation, useful for dirty sound
	else 
		{
		std::cout << "Unknown quality in Audio::repitch: " << qual << std::endl;
		return Audio();
		}
	const size_t chs = getNumChannels();
	const size_t insamples = getNumFrames();
	std::vector<float> rsoutbuf( chs * granul );

	size_t count = 0;
	size_t outcount = 0;
	while( count < getNumFrames() )
		{
		const double factor_to_use = safeFactor( count );
		rs.SetRates( getSampleRate(), double( getSampleRate() ) * factor_to_use );
		WDL_ResampleSample* rsinbuf = nullptr;
		const int wanted = rs.ResamplePrepare( granul, chs, &rsinbuf );

		for( size_t channel = 0; channel < chs; ++channel )
			for( size_t j = 0; j < wanted; ++j )
				{
				if( count + j < insamples )
					rsinbuf[ j * chs + channel ] = getSample( channel, count + j );
				else
					rsinbuf[ j * chs + channel ] = 0.0;
				}

		rs.ResampleOut( rsoutbuf.data(), wanted, granul, chs );
		for( size_t i = 0; i < chs; ++i )
			for( size_t j = 0; j < granul; ++j )
				if( outcount + j < acclen )
					out.setSample( i, outcount + j, rsoutbuf[ j  *chs + i ] );

		outcount += granul;
		count += wanted;
		}

	return out;
	}

Audio Audio::convolve( const std::vector<Func1x1> & ir ) const
	{
	std::cout << "Audio::convolve ... " << std::endl;
	auto format = getFormat();
	format.numFrames = getNumFrames() + ir.size();
	Audio out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t outSample = 0; outSample < out.getNumFrames(); ++outSample )
			{
			float convolutionSum = 0;
			for( int irSample = 0; irSample < ir.size(); ++irSample )
				{
				const int inSample = outSample + 1 - irSample;
				if( 0 <= inSample && inSample < getNumFrames() )
					convolutionSum += ir[irSample]( out.frameToTime() * outSample ) * getSample( channel, inSample );
				}
			out.setSample( channel, outSample, convolutionSum );
			}

	return out;
	}

Audio Audio::delay( float delayTime, size_t numDelays, float decay, Audio::Mod mod, bool fbIterate ) const
	{
	std::cout << "Audio::delay ... " << std::endl;

	if( mod == nullptr )
		{
		auto format = getFormat();
		format.numFrames = getNumFrames() + timeToFrame() * delayTime * numDelays;
		Audio out( format );

		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			for( size_t delay = 0; delay < numDelays+1; ++delay )
				{
				const float fbDecay = pow( decay, delay );
				const size_t outSampleOffset = delay * timeToFrame() * delayTime;
				for( size_t inSample = 0; inSample < getNumFrames(); ++inSample )
					{
					out.getSample( channel, inSample + outSampleOffset ) += getSample( channel, inSample ) * fbDecay;
					}
				}
			}

		return out;
		}
	else
		{
		std::vector<Audio> streams;
		streams.push_back( *this );
		if( fbIterate )
			for( size_t stream = 1; stream <= numDelays; ++stream )
				streams.push_back( mod( streams[stream-1], stream ).modifyVolume( decay ) );
		else
			for( size_t stream = 1; stream <= numDelays; ++stream )
				streams.push_back( mod( *this, stream ).modifyVolume( pow( decay, stream ) ) );

		std::vector<float> delayTimes( streams.size() );
		for( size_t delay = 0; delay < delayTimes.size(); ++delay )
			delayTimes[delay] = float( delay + 1 ) * delayTime;

		return mix( streams, {}, delayTimes);
		}
	}

Audio Audio::fades( float fadeTime, bool start, bool end ) const
	{
	std::cout << "Audio::fades ... " << std::endl;
	fadeTime = std::min( fadeTime, float( getLength() ) / 2.0f );

	Audio out( *this );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < timeToFrame() * fadeTime; ++sample )
			{
			const float fadeAmt = 0.5f - cos( pi * frameToTime() * sample / fadeTime ) / 2.0f;
			if( start ) out.getSample( channel, sample						) *= fadeAmt;
			if( end	  ) out.getSample( channel, getNumFrames() - 1 - sample ) *= fadeAmt;
			}

	return out;
	}

Audio Audio::lowPass( Func1x1 cutoff, size_t taps ) const
	{
	std::cout << "Audio::lowPass ... " << std::endl;
	std::vector<Func1x1> ir;
	for( size_t tap = 0; tap <= taps; ++tap ) 
		ir.push_back( [&cutoff, this, taps, tap]( float t )
			{
			const float n = float(tap) - int(taps) / 2;
			const float N = getSampleRate();
			const float K = 2.0f * cutoff(t);
			if( n == 0 ) return K / N;
			else return std::sin( pi * n * K / N ) / ( N * sin( pi * n / N ) ) * window::Hann( float(tap) / float(taps) ); 
			} );

	return convolve( ir );
	}

//========================================================
// Multi-In Procs
//========================================================

Audio Audio::mix( Audio::Vec ins, std::vector<Func1x1> balances, std::vector<float> startTimes  )
	{
	std::cout << "Audio::mix ... " << std::endl;
	if( ins.size() == 0 ) return Audio();
	const float oneOverN = 1.0f / float( ins.size() );

	if( balances.size() < ins.size() ) 
		balances.insert( balances.end(), ins.size() - balances.size(), 
			[oneOverN]( float ){ return oneOverN; } );
	if( startTimes.size() < ins.size() )
		startTimes.insert( startTimes.end(), ins.size() - startTimes.size(), 0.0 );

	auto format = ins[0].getFormat();
	format.numChannels = getMaxNumChannels( ins );
	//Find how long the longest thing being mixed will last
	for( size_t in = 0; in < ins.size(); ++in )
		{
		format.numFrames = std::max( 
			format.numFrames, 
			ins[in].getNumFrames() + size_t( ins[in].timeToFrame() * startTimes[in] ) );
		}
	Audio out( format );
	out.clearBuffer();

	for( size_t in = 0; in < ins.size(); ++in )
		for( size_t channel = 0; channel < ins[in].getNumChannels(); ++channel )
			for( size_t sample = 0; sample < ins[in].getNumFrames(); ++sample )
				{
				const size_t outSample = sample + ins[in].timeToFrame() * startTimes[in];
				out.getSample( channel, outSample ) += 
					ins[in].getSample( channel, sample ) * balances[in]( ins[0].frameToTime() * outSample );
				}

	return out;
	}

Audio Audio::join( Audio::Vec ins )
	{
	std::cout << "Audio::join ... " << std::endl;
	if( ins.size() == 0 ) return Audio();

	auto format = ins[0].getFormat();
	format.numChannels = getMaxNumChannels( ins );
	int numSamples = 0;
	for( auto & in : ins ) numSamples += in.getNumFrames();
	format.numFrames  = numSamples;
	Audio out( format );

	int currentOutStartSample = 0;

	//For each in, copy in to out
	for( size_t in = 0; in < ins.size(); ++in )
		{
		for( size_t channel = 0; channel < ins[in].getNumChannels(); ++channel )
			{
			for( size_t sample = 0; sample < ins[in].getNumFrames(); ++sample )
				{
				out.setSample( channel, currentOutStartSample + sample, ins[in].getSample( channel, sample ) );
				}
			}
		currentOutStartSample += ins[in].getNumFrames();
		}

	return out;
	}
