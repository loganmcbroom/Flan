#include "Audio.h"

#include <iostream>
#include <algorithm>
#include <complex>

#include <samplerate.h>
#include <fftw3.h>

#include "RealFunc.h"
#include "PVOC.h"
#include "WindowFunctions.h"

using namespace xcdp;

const double pi = std::acos( -1.0 );

//========================================================
// Utility Functions
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
		return a.getNumSamples() < b.getNumSamples();
		} )->getNumSamples();
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

Audio Audio::operator-( const Audio other ) const	
	{
	return *this + -other;
	}

//========================================================
// Information
//========================================================

double Audio::getTotalEnergy() const
    {
    double sum = 0;
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
            sum += pow( getSample( channel, sample ), 2 );
    
    return sum;
    }

//========================================================
// Conversions
//========================================================

PVOC Audio::convertToPVOC( const size_t frameSize, const size_t overlaps ) const
	{
	std::cout << "Getting PVOC ... \n";
	//Figure out some helpfull quantities
	const size_t numBins = frameSize / 2 + 1;
	const size_t hopSize = frameSize / overlaps;
	const size_t numFrames = getNumSamples();
	//+1 since we analyze at start and end times
	const size_t numHops = size_t( ceil( numFrames / hopSize ) ) + 1; 
	const double binWidthD = double( getSampleRate() ) / double( frameSize );
	
	PVOCBuffer::Format PVOCFormat;
	PVOCFormat.numChannels = getNumChannels();
	PVOCFormat.numFrames = numHops;
	PVOCFormat.numBins = numBins;
	PVOCFormat.sampleRate = getSampleRate();
	PVOCFormat.overlaps = overlaps;
	PVOC out( PVOCFormat );

	//Generate window
	const std::vector<double> window = window::Hann( frameSize );

	//Allocate fftw input buffer, output buffer, phase buffer, and plan
	std::vector<double> phaseBuffer( numBins );
	double * fftIn = fftw_alloc_real( frameSize );
	std::complex<double> *fftOut = (std::complex<double>*) fftw_alloc_complex( numBins );
	fftw_plan plan = fftw_plan_dft_r2c_1d( int( frameSize ), fftIn, (fftw_complex*) fftOut, FFTW_MEASURE );

	//For each channel, do the whole thing
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
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
					fftIn[relativeSample] = getSample( channel, actualSample ) * window[relativeSample];
				else
					fftIn[relativeSample] = 0;
				}
	
			fftw_execute( plan );

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
				const double phase = arg( fftOut[bin] );
				const double phaseDiff = phase - phaseBuffer[bin];
				phaseBuffer[bin] = phase;
				const double expectedPhaseDiff = double(bin) * ( 2.0 * pi / double(overlaps) );
				const double deltaPhase = phaseDiff - expectedPhaseDiff;
				const double wrappedDeltaPhase = deltaPhase-(2.0*pi)*floor((deltaPhase+pi) / ( 2.0 * pi ) );
				const double binDeviation = double(overlaps) * wrappedDeltaPhase / ( 2.0 * pi );
				const double frequency = ( double(bin) + binDeviation ) * binWidthD;

				out.setBin( channel, hop, bin, { abs( fftOut[bin] ), frequency } );
				}
			}
			
		}

	fftw_destroy_plan( plan );
	fftw_free( fftOut );
	fftw_free( fftIn );

	return out;
	}

Audio Audio::convertToMidSide() const
	{
	if( getNumChannels() != 2 )
		{
		printToLog( "Can't transform non-stereo Audio between Mid-Side and Left-Right formats. \n" );
		return *this;
		}
	Audio out( getFormat() );
	for( size_t sample = 0; sample < getNumSamples(); ++sample )
		{
		out.setSample( 0, sample, getSample( 0, sample ) + getSample( 1, sample ) );
		out.setSample( 1, sample, getSample( 0, sample ) - getSample( 1, sample ) );
		}
	return out;
	}
Audio Audio::convertToLeftRight() const
	{
	return convertToMidSide();
	}

Audio Audio::convertToStereo() const
	{
	auto format = getFormat();
	format.numChannels = 2;
	Audio out( format );

	switch( getNumChannels() )
		{
		case 1:
			{
			for( size_t sample = 0; sample < out.getNumSamples(); ++sample )
				{
				out.setSample( 0, sample, getSample( 0, sample ) / sqrt(2) );
				out.setSample( 1, sample, getSample( 0, sample ) / sqrt(2) );
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
	auto format = getFormat();
	format.numChannels = 1;
	Audio out( format );

	for( size_t sample = 0; sample < out.getNumSamples(); ++sample )
		{
		double sampleAccumulator = 0;
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
	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
			out.setSample( channel, sample, -getSample( channel, sample ) );

	return out;
	}

Audio Audio::modifyVolume( RealFunc volumeLevel ) const
	{
	std::cout << "Modifying volume ... \n";

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
			{
			double calculatedSample = getSample( channel, sample )*volumeLevel(sampleToTime(sample));
			out.setSample( channel, sample, calculatedSample );
			}

	return out;
	}

Audio Audio::setVolume( RealFunc level ) const
	{
	// Divide by getMaxSampleMagnitude to normalize, multiply by level to set
	const double maxMag = getMaxSampleMagnitude();
	return modifyVolume( [ &level, maxMag ]( double t ){ return level(t) / maxMag; } );
	}

Audio Audio::waveshape( RealFunc shaper ) const
	{
	std::cout << "Waveshaping ... \n";

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
			{
			out.setSample( channel, sample, shaper( getSample( channel, sample ) ) );
			}

	return out;
	}

Audio Audio::pan( RealFunc panAmount ) const
	{
	std::cout << "Panning ... \n";

	//Stereo panning algorithm
	auto stereoPan = [this]( RealFunc panAmount )
		{
		Audio out( getFormat() );

		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			for( size_t frame = 0; frame < getNumSamples(); ++frame )
				{
				double sample = getSample( channel, frame )
					* sin( pi / 4.0 * ( std::clamp( panAmount( sampleToTime(frame)), -1.0, 1.0 ) + 3.0 - double(channel)*2.0 ) )
					* sqrt( 2.0 );
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

Audio Audio::widen( RealFunc widenAmount ) const
	{
	return convertToMidSide().pan( widenAmount ).convertToLeftRight();
	}

Audio Audio::iterate( size_t n, Audio::Mod mod, bool fbIterate ) const
	{
	std::cout << "Iterating ... \n";

	if( mod == nullptr )
		{
		auto format = getFormat();
		format.numSamples = getNumSamples() * n;
		Audio out( format );
			
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			size_t outSample = 0;
			for( size_t i = 0; i < n; ++i )
				{
				for( size_t inSample = 0; inSample < getNumSamples(); ++inSample )
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
		totalSamplesGenerated += modOutputs[0].getNumSamples();
		for( size_t i = 1; i < n; ++i )
			{
			modOutputs.push_back( mod( fbIterate? modOutputs[i-1] : *this, i ) );
			totalSamplesGenerated += modOutputs[i].getNumSamples();
			}
	
		auto format = getFormat();
		format.numSamples = totalSamplesGenerated;
		Audio out( format );
			
		//Append mod outputs to out
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			size_t outSample = 0;
			for( size_t i = 0; i < n; ++i )
				{
				for( size_t modSample = 0; modSample < modOutputs[i].getNumSamples(); ++modSample )
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
	std::cout << "Reversing ... \n";

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
			{
			out.setSample( channel, sample, getSample( channel, getNumSamples() - 1 - sample ) );
			}

	return out;
	}

Audio Audio::cut( double startTime, double endTime ) const
	{
	std::cout << "Cutting ... \n";
	//input validity checking
	if( endTime < startTime ) endTime = startTime;
	const size_t startSample = std::max( size_t(0),       timeToSample( startTime ) );
	const size_t endSample   = std::min( getNumSamples(), timeToSample( endTime   ) );

	auto format = getFormat();
	format.numSamples =  endSample - startSample;
	Audio out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t sample = 0; sample < out.getNumSamples(); ++sample )
			out.setSample( channel, sample, getSample( channel, startSample + sample ) );

	return out;
	}

Audio Audio::repitch( RealFunc factor ) const
	{
	std::cout << "Repitching ... \n";
	//We can't have true continuous pitch shifting, so we approximate the scaling factor
	// piecewise linearly. Blocksize decides how often the continuous factor is sampled
	const size_t blockSize = getSampleRate() / 100; //100 times per second
	const size_t numBlocks = size_t( floor( getNumSamples() / blockSize ) );

	//Protection from bad factor outputs
	auto safeFactor = [factor]( double t ){ return std::clamp( factor(t), .001, 1000.0 ); };

	//Estimate output buffer length via integration (trapezoidal approximation)
	//libsamplerate uses linear pitch shifting as well so this should be within a sample
	//of the actual needed buffer length. It is rounded up, and there is an additional 
	//protection against bad array access later.
	double outLength = 0.5 / safeFactor( 0 ) 
					 + 0.5 / safeFactor( sampleToTime( blockSize * numBlocks ) ); // first and last terms
	for( size_t i = 1; i < numBlocks; ++i ) //Middle terms
		outLength += 1.0 / safeFactor( sampleToTime( i * blockSize ) );
	outLength *= blockSize; //Multiply by trapezoid base size
	//Handle remaining ( in length modulo blockSize ) samples
	outLength += ( getNumSamples() % blockSize ) 
					* ( 0.5 / safeFactor( sampleToTime( numBlocks * blockSize ) ) 
					  + 0.5 / safeFactor( sampleToTime( getNumSamples() )   ) );

	auto format = getFormat();
	format.numSamples = size_t( ceil( outLength ) );
	Audio out( format );


	//Converter initialization
	int error = 0;
	SRC_STATE * state = src_new( 1, int( getNumChannels() ), &error ); //type 0 resampling, best quality
	if( state == nullptr )
		{
		printToLog( "Error creating sample rate converter in repitch.\n" );
		return *this;
		}

	//Conversion settings
	SRC_DATA data = {};

	//What the heck are these two? These, and the block+1 and numBlocks+1 later, are 
	// there to correctly lerp between pitches. If we don't start at 0 and set the initial
	// pitch, the first block will have constant pitch
	data.input_frames = 0; 
	src_process( state, &data );

	//Interleaved float buffers
	std::vector<float> inF ( getNumSamples() * getNumChannels() );
	std::vector<float> outF( int( ceil(outLength) ) * out.getNumChannels() );

	//Load data into float buffer
	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumSamples(); ++frame )
			inF[getNumChannels()*frame+channel] = float( getSample( channel, frame ) );

	//Set up data struct used by libsamplerate
	data.data_in = inF.data(); //We will increment these as needed in the main loop
	data.data_out = outF.data();
	data.input_frames  = long( blockSize ); //Process one block at a time
	data.output_frames = long( out.getNumSamples() );
	data.end_of_input = 0;

	//Process blocks
	for( size_t block = 0; block < numBlocks; ++block )
		{
		//Note the 1/factor, we are increasing pitch by decreasing sample rate and then playing back at
		//the original sample rate
		data.src_ratio = 1.0 / safeFactor( sampleToTime( (block + 1) * blockSize ) );
		src_process( state, &data );
		data.data_in  += blockSize * getNumChannels();
		data.data_out += data.output_frames_gen * out.getNumChannels();
		data.output_frames -= data.output_frames_gen; //This protects from bad array access
		}

	//Process remainder of in samples that didn't fill a block
	data.end_of_input = 1;
	data.input_frames  = long( getNumSamples() % blockSize );
	data.src_ratio = 1.0 / safeFactor( sampleToTime( (numBlocks + 1) * blockSize ) );
	src_process( state, &data );

	src_delete( state );

	//Transform back, again if we sould use src on doubles this would be waaay faster
	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t frame = 0; frame < out.getNumSamples(); ++frame )
			{
			//if( channel == 0 and frame < 100 )
			//	std::cout << " " << inF[getNumChannels()*frame+channel];
			out.setSample( channel, frame, (double) outF[getNumChannels()*frame+channel] );
			}

	return out;                   
	}

Audio Audio::convolve( const std::vector<double> & ir ) const
	{
	std::cout << "Convolving ... \n";
	auto format = getFormat();
	format.numSamples = getNumSamples() + ir.size();
	Audio out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t sample = 0; sample < out.getNumSamples(); ++sample )
			{
			double convolutionSum = 0;
			for( int i = 0; i < ir.size(); ++i )
				{
				const int inSample = sample + 1 - i;
				if( 0 <= inSample && inSample < getNumSamples() )
					convolutionSum += ir[i] * getSample( channel, inSample );
				}
			out.setSample( channel, sample, convolutionSum );
			}

	return out;
	}

Audio Audio::delay( double delayTime, size_t numDelays, double decay, Audio::Mod mod, bool fbIterate ) const
	{
	std::cout << "Delaying ... \n";

	if( mod == nullptr )
		{
		auto format = getFormat();
		format.numSamples = getNumSamples() + timeToSample( delayTime ) * numDelays;
		Audio out( format );

		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			for( size_t delay = 0; delay < numDelays+1; ++delay )
				{
				const double fbDecay = pow( decay, delay );
				const size_t outSampleOffset = delay * timeToSample(delayTime);
				for( size_t inSample = 0; inSample < getNumSamples(); ++inSample )
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

		std::vector<double> delayTimes( streams.size() );
		for( size_t delay = 0; delay < delayTimes.size(); ++delay )
			delayTimes[delay] = double( delay + 1 ) * delayTime;

		return mix( streams, {}, delayTimes);
		}
	}

Audio Audio::fades( double fadeTime ) const
	{
	fadeTime = std::min( fadeTime, getLength() / 2.0 );

	Audio out( *this );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t sample = 0; sample < timeToSample( fadeTime ); ++sample )
			{
			const double fadeAmt = 0.5 - cos( pi * sampleToTime( sample ) / fadeTime ) / 2.0;
			out.getSample( channel, sample						 ) *= fadeAmt;
			out.getSample( channel, getNumSamples() - 1 - sample ) *= fadeAmt;
			}

	return out;
	}
	

//========================================================
// Multi-In Procs
//========================================================

Audio Audio::mix( Audio::Vec ins, std::vector< RealFunc > balances, std::vector< double > startTimes  )
	{
	if( ins.size() == 0 ) return Audio();
	const double oneOverN = 1.0 / double( ins.size() );

	if( balances.size() < ins.size() ) 
		balances.insert( balances.end(), ins.size() - balances.size(), 
			[oneOverN]( double ){ return oneOverN; } );
	if( startTimes.size() < ins.size() )
		startTimes.insert( startTimes.end(), ins.size() - startTimes.size(), 0.0 );

	auto format = ins[0].getFormat();
	format.numChannels = getMaxNumChannels( ins );
	//Find how long the longest thing being mixed will last
	for( size_t in = 0; in < ins.size(); ++in )
		{
		format.numSamples = std::max( 
			format.numSamples, 
			ins[in].getNumSamples() + ins[in].timeToSample(startTimes[in]) );
		}
	Audio out( format );
	out.clearBuffer();

	for( size_t in = 0; in < ins.size(); ++in )
		for( size_t channel = 0; channel < ins[in].getNumChannels(); ++channel )
			for( size_t sample = 0; sample < ins[in].getNumSamples(); ++sample )
				{
				const size_t outSample = sample + ins[in].timeToSample(startTimes[in]);
				out.getSample( channel, outSample ) += 
					ins[in].getSample( channel, sample ) * balances[in]( ins[0].sampleToTime( outSample ) );
				}

	return out;
	}

Audio Audio::join( Audio::Vec ins )
	{
	if( ins.size() == 0 ) return Audio();

	auto format = ins[0].getFormat();
	format.numChannels = getMaxNumChannels( ins );
	int numSamples = 0;
	for( auto & in : ins ) numSamples += in.getNumSamples();
	format.numSamples  = numSamples;
	Audio out( format );

	int currentOutStartSample = 0;

	//For each in, copy in to out
	for( size_t in = 0; in < ins.size(); ++in )
		{
		for( size_t channel = 0; channel < ins[in].getNumChannels(); ++channel )
			{
			for( size_t sample = 0; sample < ins[in].getNumSamples(); ++sample )
				{
				out.setSample( channel, currentOutStartSample + sample, ins[in].getSample( channel, sample ) );
				}
			}
		currentOutStartSample += ins[in].getNumSamples();
		}

	return out;
	}