
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

#include <algorithm>
#include <iostream>

#include "SpectralProcs.h"

namespace xcdp::spec {

SpectralBuffer average( SpecInput in, int factor )
	{
	factor = std::max( 0, factor );
	SpectralBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		{
		for( int bin = 0; bin < out.getNumBins(); ++bin )
			{
			int start = std::max( 0,				bin - factor	 );
			int end	  = std::min( in.getNumBins(),	bin + factor + 1 );
			double averageNorm = 0;
			for( int i = start; i < end; ++i )
				averageNorm += abs( in.getBin( channel, i ) );
			averageNorm /= double( 1 + 2 * factor );
			out.setBin( channel, bin, 
				std::polar( averageNorm, arg( in.getBin( channel, bin ) ) ) );
			}
		}

	return out;
	}

SpectralBuffer gate( SpecInput in, RealFunc cutoff )
	{
	SpectralBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		for( int bin = 0; bin < out.getNumBins(); ++bin )
			out[channel][bin] = abs(in[channel][bin]) > cutoff( in.getBinFreq( bin ) ) ? in[channel][bin] : 0;
	

	return out;
	}
SpectralBuffer gate( SpecInput in, double cutoff )
	{
	return gate( in, [cutoff]( double f ){ return cutoff; } );
	}

SpectralBuffer invert( SpecInput in, int lowerBound, int upperBound  )
	{
	lowerBound = std::clamp( lowerBound, 0, in.getNumBins()-1 );
	upperBound = std::clamp( upperBound, 0, in.getNumBins()-1 );
	upperBound = std::max( lowerBound, upperBound );

	SpectralBuffer out( in );
	for( int channel = 0; channel < out.getNumChannels(); ++channel )
		std::reverse( out.buffer[channel].begin() + lowerBound, out.buffer[channel].begin() + upperBound );

	return out;
	}

} //End namespace xcdp