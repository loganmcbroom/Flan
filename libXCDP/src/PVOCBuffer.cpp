#include "PVOCBuffer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

namespace xcdp {

//======================================================
//	Getters
//======================================================

PVOCBuffer::PVOCBuffer( const Format & other )
	: format( other )
	, buffer( getNumChannels() * getNumFrames() * getNumBins() )
	{}

PVOCBuffer::MFPair PVOCBuffer::getBin( size_t channel, size_t frame, size_t bin ) const
	{
	return buffer[getPos( channel, frame, bin )];
	}

size_t PVOCBuffer::getNumChannels() const
	{
	return format.numChannels;
	}
size_t PVOCBuffer::getNumFrames() const
	{
	return format.numFrames;
	}
size_t PVOCBuffer::getNumBins() const
	{
	return format.numBins;
	}
size_t PVOCBuffer::getFrameSize() const
	{
	return ( getNumBins() - 1 ) * 2;
	}
PVOCBuffer::Format PVOCBuffer::getFormat() const
	{
	return format;
	}
size_t PVOCBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}
size_t PVOCBuffer::getOverlaps() const
	{
	return format.overlaps;
	}
double PVOCBuffer::getBinWidth() const
	{
	return double( getSampleRate() ) / double( getFrameSize() );
	}
double PVOCBuffer::getMaxPartialMagnitude() const
	{
	double maxMagnitude = 0;
	for( auto i : buffer )
		{
		const double trueMag = std::abs( i.magnitude );
		if( trueMag > maxMagnitude )
			maxMagnitude = trueMag;
		}
	return maxMagnitude;
	}

size_t PVOCBuffer::timeToFrame( double t ) const
	{
	return size_t( t * double(getSampleRate()) * double(getOverlaps()) / double(getFrameSize()) );
	}

double PVOCBuffer::frameToTime( size_t frame ) const
	{
	return frame * double(getFrameSize()) / double(getSampleRate()) / double(getOverlaps());
	}

//======================================================
//	Setters
//======================================================

void PVOCBuffer::setBin( size_t channel, size_t frame, size_t bin, MFPair value )
	{
	buffer[getPos( channel, frame, bin )] = value;
	}
PVOCBuffer::MFPair & PVOCBuffer::getBin( size_t channel, size_t frame, size_t bin )
	{
	return buffer[getPos( channel, frame, bin )];
	}

void PVOCBuffer::clearBuffer()
	{
	std::fill( buffer.begin(), buffer.end(), MFPair({ 0,0 }) );
	}

size_t PVOCBuffer::getPos( size_t c, size_t f, size_t b ) const
	{
	return c * getNumFrames() * getNumBins() + f * getNumBins() + b;
	}

} // End namespace xcdp