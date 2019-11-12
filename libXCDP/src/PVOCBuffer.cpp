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
		if( i.magnitude > maxMagnitude )
			maxMagnitude = i.magnitude;
	return maxMagnitude;
	}

long long PVOCBuffer::timeToFrame( double t ) const
{
	// sample  /  hopSize 
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

/*
void PVOCBuffer::setNumChannels( size_t newNumChannels )
	{
	setBufferSize( newNumChannels, numFrames, numBins );
	}

void PVOCBuffer::setNumFrames( size_t newNumFrames )
	{
	setBufferSize( numChannels, newNumFrames, numBins );
	}

void PVOCBuffer::setNumBins( size_t newNumBins )
	{
	setBufferSize( numChannels, numFrames, newNumBins );
	}

void PVOCBuffer::setBufferSize( size_t newNumChannels, size_t newNumFrames, size_t newNumBins )
	{
	format.numChannels = newNumChannels;
	format.numFrames = newNumFrames;
	format.numBins = newNumBins;
	buffer.resize( numChannels * numFrames * numBins );
	}

void PVOCBuffer::setBufferSize( const PVOCBuffer & other )
	{
	//setNumChannels	( other.getNumChannels()	);
	//setNumFrames	( other.getNumFrames()		);
	//setNumBins		( other.getNumBins()		);
	setBufferSize( other.getNumChannels(), other.getNumFrames(),  other.getNumBins() );
	}

void PVOCBuffer::setSampleRate( size_t newSampleRate )
	{
	format.sampleRate = newSampleRate;
	}

void PVOCBuffer::setOverlaps( size_t newOverlaps )
	{
	overlaps = newOverlaps;
	}

void PVOCBuffer::copyFormat( const PVOCBuffer & other )
	{
	sampleRate = other.sampleRate;
	overlaps = other.overlaps;
	}
*/

void PVOCBuffer::clearBuffer()
	{
	//for( auto & channel : buffer )
	//	for( auto & frame : channel )
	//		for( auto & bin : frame )
	//			bin = {0,0};
	std::fill( buffer.begin(), buffer.end(), MFPair({ 0,0 }) );
	}

size_t PVOCBuffer::getPos( size_t c, size_t f, size_t b ) const
	{
	return c * getNumFrames() * getNumBins() + f * getNumBins() + b;
	}

} // End namespace xcdp