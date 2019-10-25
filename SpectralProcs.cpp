
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

#include "SpectralProcs.h"

namespace xcdp {

SpectralBuffer blur_average( SpecInput in, int factor )
	{
	SpectralBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in.getNumChannels(), in.getNumBins() );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		{
		for( int bin = 0; bin < out.getNumBins(); ++bin )
			{
			int start = std::max( 0,				bin - factor );
			int end	  = std::min( in.getNumBins(),	bin + factor );
			std::complex<float> average = 0;
			for( int i = start; i < end; ++i )
				average += in.getBin( channel, i );
			out.setBin( channel, bin, average / std::complex<float>( end - start, 0 ) );
			}
		}

	return out;
	}

} //End namespace xcdp