#include "SPVBuffer.h"

#include <execution>

using namespace flan;

SPVBuffer::SPVBuffer()
	: format()
	, buffer()
	{}

SPVBuffer::SPVBuffer( Format _format )
	: format( _format )
	, buffer( getNumChannels() * getNumFrames() * getNumBins() )
	{}

const SPVBuffer::Format & SPVBuffer::getFormat() const 
	{
	return format;
	}

fFrame SPVBuffer::timeToFrame( Second t ) const 
	{
	return t * getSampleRate();
	}

Second SPVBuffer::frameToTime( fFrame f ) const 
	{
	return f / getSampleRate();
	}

fBin SPVBuffer::frequencyToBin( Frequency f ) const 
	{
	return f * getNumBins() / getSampleRate();
	}

Frequency SPVBuffer::binToFrequency( fBin b ) const 
	{
	return b * getSampleRate() / getNumBins();
	}

Channel SPVBuffer::getNumChannels() const
	{
	return format.numChannels;
	}

Frame SPVBuffer::getNumFrames() const
	{
	return format.numFrames;
	}

Bin	SPVBuffer::getNumBins() const
	{
	return format.numBins;
	}

FrameRate SPVBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}

FrameRate SPVBuffer::getAnalysisRate() const
	{
	return format.sampleRate;
	}

bool SPVBuffer::isNull() const 
	{ 
	return getSampleRate() <= 0
		|| buffer.size() == 0;
	}

MF SPVBuffer::getMF( Channel channel, Frame frame, Bin bin ) const 
	{
	return buffer[ getBufferPos( channel, frame, bin ) ];
	}

MF & SPVBuffer::getMF( Channel channel, Frame frame, Bin bin )
	{
	return buffer[ getBufferPos( channel, frame, bin ) ];
	}

void SPVBuffer::clearBuffer()
	{
	std::fill( std::execution::par_unseq, buffer.begin(), buffer.end(), MF{0,0} );
	}
	
SPVBuffer SPVBuffer::copy() const
	{
	SPVBuffer out;
	out.format = format;
	out.buffer = buffer; // Deep copy
	return out;
	}

size_t SPVBuffer::getBufferPos( Channel c, Frame f, Bin b ) const
	{
	return ( c * getNumFrames() + f ) * getNumBins() + b;
	}

std::vector<MF> & SPVBuffer::getBuffer()
	{
	return buffer;
	}


