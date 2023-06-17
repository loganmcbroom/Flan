#include "SQPVBuffer.h"

#include <algorithm>
#include <execution>

using namespace flan;

SQPVBuffer::SQPVBuffer()
	: format()
	, pitch_bandwidth()
	, numBins()
	, Q()
	, bin_frequencies()
	, buffer()
	{}

SQPVBuffer::SQPVBuffer( const Format & f )
	: format( f )
	, pitch_bandwidth( { frequencyToPitch( getFrequencyBandwidth().first  ).p, 
						 frequencyToPitch( getFrequencyBandwidth().second ).p } )
	, numBins( std::ceil( frequencyToBin( getFrequencyBandwidth().second ) ) )
	, Q( 1.0f / ( std::pow( 2.0f, 1.0f / f.binsPerOctave ) - 1.0f ) )
	, bin_frequencies( [&]()
		{ 
		std::vector<Frequency> out( getNumBins() );
		for( Bin bin = 0; bin < out.size(); ++bin )
			out[bin] = binToFrequency( bin ); 
		return out;
		}() )
	, buffer( getNumChannels() * getNumFrames() * getNumBins() )
	{}

MP SQPVBuffer::getMP( Channel c, Frame f, Bin b ) const
	{
	return buffer[getBufferPos( c, f, b )];
	}

MP & SQPVBuffer::getMP( Channel c, Frame f, Bin b )
	{
	return buffer[getBufferPos( c, f, b )];
	}

const SQPVBuffer::Format & SQPVBuffer::getFormat() const
	{
	return format;
	}

fFrame SQPVBuffer::timeToFrame( Second t ) const
	{
	return t * getSampleRate();
	}

Second SQPVBuffer::frameToTime( fFrame f ) const
	{
	return f / getSampleRate();
	}

fBin SQPVBuffer::frequencyToBin( Frequency f ) const
	{
	return pitchToBin( frequencyToPitch( f ).p );
	}

Frequency SQPVBuffer::binToFrequency( fBin bin ) const
	{
	return pitchToFrequency( Pitch( binToPitch( bin ), true ) );
	}

Pitch SQPVBuffer::frequencyToPitch( Frequency f ) const
	{
	if( f == 0 ) return std::numeric_limits<UnsignedPitch>::min();
	return Pitch( std::log2( std::abs( f ) ), f >= 0 );
	}

Frequency SQPVBuffer::pitchToFrequency( Pitch p ) const
	{
	return std::pow( 2.0f, p.p ) * ( p.positive_frequency ? 1 : -1 );
	}

UnsignedPitch SQPVBuffer::binToPitch( fBin b ) const
	{
	return b / getBinsPerOctave() + getPitchBandwidth().first;
	}

fBin SQPVBuffer::pitchToBin( UnsignedPitch p ) const
	{
	return ( p - getPitchBandwidth().first ) * getBinsPerOctave();
	}

size_t SQPVBuffer::getBufferPos( Channel c, Frame f , Bin b ) const
	{
	return ( c * getNumFrames() + f ) * getNumBins() + b;
	}

Channel SQPVBuffer::getNumChannels() const
	{
	return format.numChannels;
	}

Frame SQPVBuffer::getNumFrames() const
	{
	return format.numFrames;
	}

Bin	SQPVBuffer::getNumBins() const
	{
	return numBins;
	}

FrameRate SQPVBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}

FrameRate SQPVBuffer::getAnalysisRate() const
	{
	return format.sampleRate;
	}

std::pair<Frequency, Frequency> SQPVBuffer::getFrequencyBandwidth() const
	{
	return format.bandwidth;
	}

std::pair<UnsignedPitch, UnsignedPitch> SQPVBuffer::getPitchBandwidth() const
	{
	return pitch_bandwidth;
	}

float SQPVBuffer::getBinsPerOctave() const
	{
	return format.binsPerOctave;
	}

Cycle SQPVBuffer::getQ() const
	{
	return Q;
	}

Frequency SQPVBuffer::getBinFrequency( Bin bin ) const
	{
	return bin_frequencies[bin];
	}

bool SQPVBuffer::isNull() const
	{
	return getSampleRate() <= 0
		|| buffer.size() == 0
		|| getFrequencyBandwidth().first >= getFrequencyBandwidth().second;
	}

void SQPVBuffer::clearBuffer()
	{
	std::fill( std::execution::par_unseq, buffer.begin(), buffer.end(), MP{ 0, Pitch( 0, true ) } );
	}

SQPVBuffer SQPVBuffer::copy() const
	{
	SQPVBuffer out( getFormat() );
	out.buffer = buffer; // Deep copy
	return out;
	}

Magnitude SQPVBuffer::getMaxPartialMagnitude() const
	{
	Magnitude maxMagnitude = 0;
	for( const MP & mf : buffer )
		{
		const Magnitude trueMag = std::abs( mf.m );
		if( trueMag > maxMagnitude )
			maxMagnitude = trueMag;
		}
	return maxMagnitude;
	}

Magnitude SQPVBuffer::getMaxPartialMagnitude( Frame startFrame, Frame endFrame, Bin startBin, Bin endBin ) const
	{
	if( endFrame == 0 ) endFrame = getNumFrames();
	if( endBin == 0 ) endBin = getNumBins();

	Magnitude maxMagnitude = 0;

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = startFrame; frame < endFrame; ++frame )
			for( Bin bin = startBin; bin < endBin; ++bin )
			    {
				const Magnitude trueMag = std::abs( getMP( channel, frame, bin ).m );
				if( trueMag > maxMagnitude )
					maxMagnitude = trueMag;
				}

	return maxMagnitude;
	}

Frame SQPVBuffer::getPeriod( Bin bin ) const
	{
	return std::ceil( getQ() * getSampleRate() / getBinFrequency( bin ) );
	}
