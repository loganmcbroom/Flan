#include "PVOC.h"

#include <algorithm>
#include <iostream>

/** TODO:

	blur blur
	blur chorus
	blur scatter
	stretch time
	formants get
	formants put
	hilite trace
	morph morph
	superaccu
	focus exag
	combine diff

*/

namespace xcdp {

PVOC PVOC::timeAverage( int factor ) const
	{
	std::cout << "Time averaging ... ";
	PVOC out;
	out.copyFormat( *this );
	out.setBufferSize( *this );

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			{
			size_t startFrame = floor( frame / factor ) * factor;
			size_t endFrame	  = std::min( getNumFrames()-1, startFrame + factor );
			const double interpolation = double(frame % factor) / double(factor);

			for( size_t bin = 0; bin < getNumBins(); ++bin )
				{
				out.buffer[channel][frame][bin] = 
					{ 
					  buffer[channel][startFrame][bin].magnitude*(1.0-interpolation) 
					+ buffer[channel][endFrame  ][bin].magnitude*interpolation,
					//  buffer[channel][startFrame][bin].frequency*(1.0-interpolation) 
					//+ buffer[channel][endFrame  ][bin].frequency*interpolation
					buffer[channel][startFrame][bin].frequency
					};
				}
			}

	std::cout << "Done\n";
	return out;
	}

PVOC PVOC::repitch( double factor )
	{
	std::cout << "Repitching ... ";
	PVOC out;
	out.copyFormat( *this );
	out.setBufferSize( *this );
	out.clearBuffer();

	for( size_t channel = 0; channel < getNumChannels(); ++channel )
		for( size_t frame = 0; frame < getNumFrames(); ++frame )
			for( size_t bin = 0; bin < getNumBins()/2; ++bin )
				{
				size_t index = size_t( bin * factor );
				if( index < getNumBins() )
					{
					out.buffer[channel][frame][index].magnitude += buffer[channel][frame][bin].magnitude;
					out.buffer[channel][frame][index].frequency =  buffer[channel][frame][bin].frequency * factor;
					}
				}

	std::cout << "Done\n";
	return out;
	}

//PVOC PVOC::gate( RealFunc cutoff ) const
//	{
//	SpectralBuffer out;
//	out.copyFormat( in );
//	out.setBufferSize( in );
//
//	for( int channel = 0; channel < in.getNumChannels(); ++channel )
//		for( int bin = 0; bin < out.getNumBins(); ++bin )
//			out[channel][bin] = abs(in[channel][bin]) > cutoff( in.getBinFreq( bin ) ) ? in[channel][bin] : 0;
//	
//
//	return out;
//	}
//PVOC PVOC::gate( double cutoff ) const
//	{
//	return gate( in, [cutoff]( double f ){ return cutoff; } );
//	}
//
//PVOC PVOC::invert( int lowerBound, int upperBound  ) const
//	{
//	lowerBound = std::clamp( lowerBound, 0, int( getNumBins() - 1 ) );
//	upperBound = std::clamp( upperBound, 0, int( getNumBins() - 1 ) );
//	upperBound = std::max( lowerBound, upperBound );
//
//	PVOC out( in );
//	for( int channel = 0; channel < out.getNumChannels(); ++channel )
//		std::reverse( out.buffer[channel].begin() + lowerBound, out.buffer[channel].begin() + upperBound );
//
//	return out;
//	}

} //End namespace xcdp