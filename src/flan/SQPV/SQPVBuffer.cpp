#include "SQPVBuffer.h"

#include <algorithm>
#include <execution>

using namespace flan;

SQPVBuffer::SQPVBuffer()
	: format()
	, numBins()
	, quality()
	, buffer()
	{}

SQPVBuffer::SQPVBuffer( const Format & f )
	: format( f )
	, numBins( std::ceil( f.binsPerOctave * std::log2( f.bandwidth.second / f.bandwidth.first ) ) )
	, quality( 1.0f / ( std::pow(2.0f, 1.0f / f.binsPerOctave) - 1.0f ) )
	, buffer( getNumChannels() * getNumFrames() * getNumBins() )
	{}

MF SQPVBuffer::getMF( Channel c, Frame f, Bin b ) const
	{
	return buffer[buffer_access( c, f, b )];
	}

MF & SQPVBuffer::getMF( Channel c, Frame f, Bin b )
	{
	return buffer[buffer_access( c, f, b )];
	}

const SQPVBuffer::Format & SQPVBuffer::getFormat() const
	{
	return format;
	}

fFrame SQPVBuffer::timeToFrame( Time t ) const
	{
	return t * getSampleRate();
	}

Time SQPVBuffer::frameToTime( fFrame f ) const
	{
	return f / getSampleRate();
	}

fBin SQPVBuffer::frequencyToBin( Frequency f ) const
	{
	return getBinsPerOctave() * std::log2( f / getBandwidth().first );
	}

Frequency SQPVBuffer::binToFrequency( fBin bin ) const
	{
	return getBandwidth().first * std::pow( 2.0f, bin / getBinsPerOctave() );
	}

size_t SQPVBuffer::buffer_access( Channel c, Frame f , Bin b ) const
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

SampleRate SQPVBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}

std::pair<Frequency, Frequency> SQPVBuffer::getBandwidth() const
	{
	return format.bandwidth;
	}

float SQPVBuffer::getBinsPerOctave() const
	{
	return format.binsPerOctave;
	}

float SQPVBuffer::getQuality() const
	{
	return quality;
	}

bool SQPVBuffer::isNull() const
	{
	return getSampleRate() <= 0
		|| buffer.size() == 0
		|| getBandwidth().first >= getBandwidth().second;
	}

void SQPVBuffer::clearBuffer()
	{
	std::fill( std::execution::par_unseq, buffer.begin(), buffer.end(), MF{ 0, 0 } );
	}

SQPVBuffer SQPVBuffer::copy() const
	{
	SQPVBuffer out( getFormat() );
	out.buffer = buffer; // Deep copy
	return out;
	}

