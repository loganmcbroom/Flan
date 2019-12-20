#include "Audio.h"

#include <iostream>
#include <algorithm>
#include <complex>

#include <samplerate.h>
#include <fftw3.h>

#include "RealFunc.h"
#include "PVOC.h"
#include "Spectrum.h"
#include "WindowFunctions.h"
#include "WDL/resample.h"


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

Audio Audio::graph( const std::string & filename, size_t width, size_t height ) const
	{
	const auto audioColor = HSVtoRGB( 0, .8, .65 ); // redish
	const auto backgroundColor = HSVtoRGB( 180, .8, .1 ); // blueish

	if( getNumSamples() == 0 ) return *this;
	width = std::min( width, getNumSamples() );

	std::vector<std::vector<std::array<uint8_t,3>>> data( width, std::vector( height*getNumChannels(), backgroundColor ) );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		std::vector<double> energies( width, 0.0 );
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
			energies[ double(sample) / double(getNumSamples()) * width ] += getSample( channel, sample );
		double maxEnergy = 0;
		for( size_t x = 0; x < width; ++x )
			maxEnergy = std::max( maxEnergy, std::abs( energies[x] ) );

		for( size_t x = 0; x < width; ++x )
			{
			const int height2 = height / 2;
			data[x][height2 + height * channel ] = {0,0,0};
			for( int y = 0; std::abs( y ) < std::abs( energies[x] / maxEnergy * height2 ) * 0.9; energies[x] > 0 ? ++y : --y )
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
					fftIn[relativeSample] = getSample( channel, actualSample ) * window::Hann( double(relativeSample) / double(frameSize) );
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
Spectrum Audio::convertToSpectrum() const
	{
	const size_t outnumBins = std::floor( double( getNumSamples() ) / 2.0 ) + 1;

	Spectrum::Format format;
	format.numChannels = getNumChannels();
	format.numBins = outnumBins;
	format.sampleRate = getSampleRate();
	Spectrum out( format );

	//Allocate fftw buffers and plan
	double * fftIn = fftw_alloc_real( getNumSamples() );
	std::complex<double> * fftOut = (std::complex<double>*) fftw_alloc_complex( outnumBins );
	fftw_plan plan = fftw_plan_dft_r2c_1d( int( getNumSamples() ), fftIn, (fftw_complex*) fftOut, FFTW_ESTIMATE );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		{
		//Read buffer data into fftIn
		for( size_t sample = 0; sample < getNumSamples(); ++sample )
			fftIn[sample] = getSample( channel, sample );

		//fft
		fftw_execute( plan );

		//Read fftOut into out spectrum
		for( size_t bin = 0; bin < outnumBins; ++bin )
			out.setSpectra( channel, bin, fftOut[bin] );
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
	if( maxMag == 0 ) return *this;
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

Audio Audio::repitch( RealFunc factor, size_t granul, size_t qual ) const
	{
	std::cout << "Repitching...\n";

	auto safeFactor = [this, &factor]( double count )
		{ 
		return std::max( 1.0 / factor( double(count) / getSampleRate() ), 0.001 ); 
		};

	// estimate output length
	double acclen = 0.0;
	for( int count = 0; count < getNumSamples(); count += granul )
		acclen += granul * safeFactor( count );

	auto format = getFormat();
	format.numSamples = size_t(ceil(acclen));
	Audio out( format );
	
	WDL_Resampler rs;
		 if( qual == 0 ) rs.SetMode( true,  0, true, 64 ); // reasonable sinc resampling
	else if( qual == 1 ) rs.SetMode( true,  1, false    ); // linear interpolation with simple filter
	else if( qual == 2 ) rs.SetMode( false, 0, false    ); // no interpolation, useful for dirty sound
	const size_t chs = getNumChannels();
	const size_t insamples = getNumSamples();
	std::vector<double> rsoutbuf( chs * granul );

	int count = 0;
	int outcount = 0;
	while( count < getNumSamples() )
		{
		const double factor_to_use = safeFactor( count );
		rs.SetRates( getSampleRate(), getSampleRate() * factor_to_use );
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
		for( int i = 0; i < chs; ++i )
			for( int j = 0; j < granul; ++j )
				if( outcount + j < acclen )
					out.setSample( i, outcount + j, rsoutbuf[ j  *chs + i ] );

		outcount += granul;
		count += wanted;
		}

	return out;
	}

Audio Audio::convolve( const std::vector<RealFunc> & ir ) const
	{
	std::cout << "Convolving ... \n";
	auto format = getFormat();
	format.numSamples = getNumSamples() + ir.size();
	Audio out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t outSample = 0; outSample < out.getNumSamples(); ++outSample )
			{
			double convolutionSum = 0;
			for( int irSample = 0; irSample < ir.size(); ++irSample )
				{
				const int inSample = outSample + 1 - irSample;
				if( 0 <= inSample && inSample < getNumSamples() )
					convolutionSum += ir[irSample]( out.sampleToTime( outSample ) ) * getSample( channel, inSample );
				}
			out.setSample( channel, outSample, convolutionSum );
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

Audio Audio::lowPass( RealFunc cutoff, size_t taps ) const
	{
	std::vector<RealFunc> ir;
	for( size_t tap = 0; tap <= taps; ++tap ) 
		ir.push_back( [&cutoff, this, taps, tap]( double t )
			{
			const double n = double(tap) - int(taps) / 2;
			const double N = getSampleRate();
			const double K = 2.0 * cutoff(t);
			if( n == 0 ) return K / N;
			else return sin( pi * n * K / N ) / ( N * sin( pi * n / N ) ) * window::Hann( double(tap) / double(taps) ); 
			} );

	return convolve( ir );
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