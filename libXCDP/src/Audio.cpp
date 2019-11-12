/** TODO:
	housekeep chans
	extend loop
	filter userbank
	filter varibank(2)
	RealFunc ADSR( double A, double D, double S, double R, double Aexp = 0, double Dexp = 0, double Rexp = 0 );
	distort average
	distort interpolate
	distort multiply
*/

#include "Audio.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <complex>

#include <samplerate.h>
#include <fftw3.h>

#include "PVOC.h"
#include "WindowFunctions.h"

namespace xcdp {

const double pi = std::acos( -1.0 );

//========================================================
// Utility Functions
//========================================================

// You want print to log? $5 print to log best deal in town.
void printToLog( std::string s )
	{
	std::cout << s << std::endl;
	}
bool doSampleRatesMatch( const std::vector<Audio> & ins )
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
bool doChannelCountsMatch( const std::vector<Audio> & ins )
	{
	//check if all channel counts match the first file
	size_t numchannels = ins[0].getNumChannels();
	for( auto & in : ins ) 
		if( in.getNumChannels() != numchannels )
			{
			printToLog( "Mismatched channel count" );
			return false;
			}
	return true;
	}

size_t getMaxNumChannels( const std::vector<Audio> & ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumChannels() < b.getNumChannels();
		} )->getNumChannels();
	}
size_t getMaxNumSamples( const std::vector<Audio> &  ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumSamples() < b.getNumSamples();
		} )->getNumSamples();
	}

Audio Audio::mono2Stereo( ) const
	{
	if( getNumChannels() != 1 )
		printToLog( "That wasn't a mono file you tried to make stereo just now.\n" );

	auto format = getFormat();
	format.numChannels = 2;
	Audio stereo( format );

	for( size_t frame = 0; frame < stereo.getNumSamples(); ++frame )
		{
		stereo.setSample( 0, frame, getSample( 0, frame ) );
		stereo.setSample( 1, frame, getSample( 0, frame ) );
		}

		return stereo;
	}

//========================================================
// Conversions
//========================================================

PVOC Audio::convertToPVOC( const size_t frameSize, const size_t overlaps ) const
	{
	std::cout << "Getting PVOC ... ";
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
		
		for( size_t relativeSample = 0; relativeSample < frameSize; ++relativeSample )
			{
			// - overlaps/2 to center window around analysis point
			const int actualSample = int(relativeSample) + int(hopSize) * (-1 - int(overlaps/2));
			if( 0 <= actualSample && actualSample < numFrames )
				fftIn[relativeSample] = getSample( channel, actualSample ) * window[relativeSample];
			else
				fftIn[relativeSample] = 0;
			}

		fftw_execute( plan );
		for( size_t bin = 0; bin < numBins; ++bin )
			phaseBuffer[bin] = arg( fftOut[bin] );

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

	std::cout << "Done\n";
	return out;
	}

//========================================================
// Procs
//========================================================

Audio Audio::modifyVolume( RealFunc volumeLevel ) const
	{
	std::cout << "Modifying volume ... ";

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumSamples(); ++frame )
			{
			double calculatedSample = getSample( channel, frame )*volumeLevel(getTimeOfFrame(frame));
			out.setSample( channel, frame, calculatedSample );
			}

	std::cout << "Done\n";
	return out;
	}
Audio Audio::modifyVolume( double volumeLevel ) const
	{
	return modifyVolume( [ volumeLevel ]( double ){ return volumeLevel; } );
	}

Audio Audio::setVolume( double level ) const
	{
	// Divide by maxVolume to normalize, multiply by level to set
	double volume = level / getMaxSampleMagnitude();
	return modifyVolume( [ volume ]( double t ){ return volume; } );
	}

Audio Audio::waveshape( RealFunc shaper ) const
	{
	std::cout << "Waveshaping ... ";

	Audio out( getFormat() );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumSamples(); ++frame )
			{
			out.setSample( channel, frame, shaper( getSample( channel, frame ) ) );
			}

	std::cout << "Done\n";
	return out;
	}

Audio Audio::pan( RealFunc panAmount ) const
	{
	std::cout << "Panning ... ";

	//Stereo panning algorithm
	auto stereoPan = [this]( RealFunc panAmount )
		{
		Audio out( getFormat() );

		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			for( size_t frame = 0; frame < getNumSamples(); ++frame )
				{
				double sample = getSample( channel, frame )
					* sin( pi / 4.0 * ( std::clamp( panAmount( getTimeOfFrame(frame)), -1.0, 1.0 ) + 3.0 - double(channel)*2.0 ) )
					* sqrt( 2.0 );
				out.setSample( channel, frame, sample );
				}
		std::cout << "Done\n";
		return out;
		};
	
	switch( getNumChannels() )
		{
		case 1: //Panning mono: cast to stereo and do stereo pan
			return mono2Stereo().pan( panAmount );
			break;
		case 2: //Panning stereo: Use stereo pan as defined above
			return stereoPan( panAmount );
			break;
		default:
			printToLog( "I don't know how to pan that number of channels" );
			return *this;
		}
	}
Audio Audio::pan( double panAmount ) const
	{
	return pan( [panAmount](double){ return panAmount; } );
	}

Audio Audio::iterate( size_t n, std::function< Audio ( const Audio &, size_t n) > mod ) const
	{
	std::cout << "Iterating ... ";
	
	if( mod == 0 ) //No mod
		{
		auto format = getFormat();
		format.numSamples = getNumSamples() * n;
		Audio out( format );

		//Append this to out n times
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			size_t outFrame = 0;
			for( size_t i = 0; i < n; ++i )
				{
				for( size_t inFrame = 0; inFrame < getNumSamples(); ++inFrame )
					{
					out.setSample( channel, outFrame, getSample( channel, inFrame ) );
					++outFrame;
					}
				}
			}
		std::cout << "Done\n";
		return out;
		}
	else //Mod present
		{
		//For each iteration, apply mod to input
		std::vector<Audio> modOutputs;
		size_t totalFramesGenerated = 0;
		for( size_t i = 0; i < n; ++i )
			{
			modOutputs.push_back( mod( *this, i ) );
			totalFramesGenerated += modOutputs[i].getNumSamples();
			}
		auto format = getFormat();
		format.numSamples = totalFramesGenerated;
		Audio out( format );
			
		//Append mod outputs to out
		for( size_t channel = 0; channel < getNumChannels(); ++channel )
			{
			size_t outSample = 0;
			for( size_t i = 0; i < n; ++i )
				{
				for( size_t modSample = 0; modSample < modOutputs[i].getNumSamples(); ++modSample )
					{
					//std::cout << outSample << " ";
					out.setSample( channel, outSample, modOutputs[i].getSample( channel, modSample ) );
					++outSample;
					}
				}
			}
		std::cout << "Done\n";
		return out;
		}
	}

Audio Audio::cutAtSamples( size_t startSample, size_t endSample ) const
	{
	std::cout << "Cutting at samples ... ";
	//input validity checking
	startSample = std::min( size_t(0), startSample );
	endSample = std::max( getNumSamples(), endSample );
	if( endSample < startSample ) endSample = startSample;

	auto format = getFormat();
	format.numSamples =  endSample - startSample;
	Audio out( format );

	for( size_t frame = 0; frame < out.getNumSamples(); ++frame )
		for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
			out.setSample( channel, frame, getSample( channel, startSample + frame ) );

	std::cout << "Done\n";
	return out;
	}
Audio Audio::cutAtTimes( double startTime, double endTime ) const
	{
	//Input time validity is checked in cutAtSamples
	return cutAtSamples( size_t(startTime / getSampleRate()), size_t(endTime / getSampleRate()) );
	}

Audio Audio::repitch( RealFunc factor ) const
	{
	std::cout << "Repitching ... ";
	//We can't have true continuous pitch shifting, so we approximate the scaling factor
	// piecewise linearly. Blocksize decides how often the continuous factor is sampled
	const size_t blockSize = getSampleRate() / 32; //32 times per second
	const size_t numBlocks = size_t( floor( getNumSamples() / blockSize ) );

	//Protection from bad factor outputs
	auto safeFactor = [factor]( double t ){ return std::clamp( factor(t), .001, 1000.0 ); };

	//Estimate output buffer length via integration (trapezoidal approximation)
	//libsamplerate uses linear pitch shifting as well so this should be within a sample
	//of the actual needed buffer length. It is rounded up, and there is an additional 
	//protection against bad array access later.
	double outLength = 0.5 / safeFactor( 0 ) 
					 + 0.5 / safeFactor( getTimeOfFrame( blockSize * numBlocks ) ); // first and last terms
	for( size_t i = 1; i < numBlocks; ++i ) //Middle terms
		outLength += 1.0 / safeFactor( getTimeOfFrame( i * blockSize ) );
	outLength *= blockSize; //Multiply by trapezoid base size
	//Handle remaining ( in length modulo blockSize ) samples
	outLength += ( getNumSamples() % blockSize ) 
					* ( 0.5 / safeFactor( getTimeOfFrame( numBlocks * blockSize ) ) 
					  + 0.5 / safeFactor( getTimeOfFrame( getNumSamples() )   ) );

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
	data.input_frames  = 0; 
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
		data.src_ratio = 1.0 / safeFactor( getTimeOfFrame( (block + 1) * blockSize ) );
		src_process( state, &data );
		data.data_in  += blockSize * getNumChannels();
		data.data_out += data.output_frames_gen * out.getNumChannels();
		data.output_frames -= data.output_frames_gen; //This protects from bad array access
		}

	//Process remainder of in samples that didn't fill a block
	data.end_of_input = 1;
	data.input_frames  = long( getNumSamples() % blockSize );
	data.src_ratio = 1.0 / safeFactor( getTimeOfFrame( (numBlocks + 1) * blockSize ) );
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

	std::cout << "Done\n";
	return out;                   
	}
Audio Audio::repitch( double factor ) const
	{
	return repitch( [factor]( double t ){ return factor; } );
	}

Audio Audio::convolve( const std::vector<double> & ir ) const
	{
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
					{
					convolutionSum += ir[i] * getSample( channel, inSample );
					}
				}
			out.setSample( channel, sample, convolutionSum );
			}
	return out;
	}

//========================================================
// TODO
//========================================================

//RealFunc ADSR( double A, double D, double S, double R, double Aexp = 0, double Dexp = 0, double Rexp = 0 );
/*
Audio Audio::mix( ProcInputVec ins, std::vector< RealFunc > balances ) const
	{
	//If we didn't get a balances input, construct a default one, otherwise check input validity
	if( balances.size() == 0 )
		{
		double insSize = ins.size();
		auto defaultBalanceFunc = [insSize]( double )
			{ 
			return 1.0 / insSize; //1/n trades off natural sound for guaranteed no clip
			};
		balances.assign( ins.size(), defaultBalanceFunc );
		}
	else if( balances.size() < ins.size() )
		{
		printToLog( "Too few balances sent to mix, padding with 1/n's. \n" );
		double oneOverN = 1.0 / double( ins.size() ); 
		balances.insert( balances.end(), ins.size() - balances.size(), [oneOverN]( double )
			{ 
			return oneOverN;
			} );
		}
	else if( balances.size() > ins.size() )
		{
		printToLog( "Too Many balances sent to mix. Not really a problem, just thought you should know. \n" );
		}

	if( !doChannelCountsMatch( ins ) )
		{
		printToLog( "You can't mix files with differing channel counts at the moment (Is there a good way?).\n" );
		return Audio();
		}

	int numChannels = getMaxNumChannels( ins );
	int numSamples  = getMaxNumSamples ( ins );

	Audio out;
	out.copyFormat( ins[0] );
	out.setBufferSize( numChannels, numSamples );

	for( size_t channel = 0; channel < numChannels; ++channel )
		for( size_t frame = 0; frame < numSamples; ++frame )
			{
			double mixedSample = 0.0;
			for( size_t i = 0; i < ins.size(); ++i )
				mixedSample += balances[i]( ins[0].getTimeOfFrame( frame ) ) * ins[i].getSample( channel, frame );
			out.setSample( channel, frame, std::clamp( mixedSample, -1.0, 1.0 ) );
			}

	return out;
	}
Audio mix( ProcInputVec ins, const std::vector< double > & balances ) const
	{
	std::vector< RealFunc > funcs( balances.size() );
	for( unsigned size_t i = 0; i < balances.size(); ++i )
		funcs[i] = [i, &balances](double){ return balances[i]; };

	return mix( ins, funcs );
	}
Audio join( AudioVec ins ) 
	{
	int numChannels = getMaxNumChannels( ins );
	int numSamples = 0;
	for( auto & in : ins ) numSamples += in.getNumSamples();

	Audio out;
	out.copyFormat( ins[0] );
	out.setBufferSize( numChannels, numSamples );

	int currentOutStartSample = 0;

	//For each in, copy in into out
	for( size_t in = 0; in < ins.size(); ++in )
		{
		for( size_t frame = 0; frame < ins[in].getNumSamples(); ++frame )
			for( size_t channel = 0; channel < ins[in].getNumChannels(); ++channel )
				out.setSample( channel, currentOutStartSample + frame, ins[in].getSample( channel, frame ) );
		currentOutStartSample += ins[in].getNumSamples();
		}

	return out;
	}
*/

} //end namespace xcdp