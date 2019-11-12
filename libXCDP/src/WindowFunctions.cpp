#include "WindowFunctions.h"

namespace xcdp::window {

const std::vector<double> Hann( size_t frameSize )
	{
	 std::vector<double> out( frameSize );
	for( size_t sample = 0; sample < frameSize; ++sample )
		 out[sample] = 0.5 * ( 1.0 - cos( 2.0 * pi * double(sample) / double(frameSize) ) );
	return out;
	}

} // End namespace xcdp::window