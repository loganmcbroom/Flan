#include "flan/Function.h"

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <string>

namespace flan {

Function<Second, Amplitude> ADSR( 
	Second a, 
	Second d, 
	Second s, 
	Second r, 
	Amplitude sLvl, 
	float aExp, 
	float dExp, 
	float rExp )
	{
	return [a, d, s, r, sLvl, aExp, dExp, rExp]( Second t ) -> Amplitude
		{
			 if( t < 0 || t > a+d+s+r ) return (float) 0.0;
		else if( t < a			   ) return std::pow( t / a, aExp );
		else if( t < a + d		   ) return std::pow( (float) 1.0 - ( t - a ) / d, dExp ) * ( (float) 1.0 - sLvl ) + sLvl;
		else if( t < a + d + s	   ) return sLvl;
		else if( t < a + d + s + r ) return std::pow( (float) 1.0 - ( t - a - d - s ) / r, rExp ) * sLvl;
		else return 0.0f;
		};
	}

namespace waveforms {

const std::function<Amplitude (Second)> sine 	 = []( Second t ){ const Second t0 = std::fmod( t, 1.0f ); return std::sin( pi2 * t0 ); };
const std::function<Amplitude (Second)> square	 = []( Second t ){ const Second t0 = std::fmod( t, 1.0f ); return t0 < 0.5f ? -1.0f : 1.0f; };
const std::function<Amplitude (Second)> saw 	 = []( Second t ){ const Second t0 = std::fmod( t, 1.0f ); return -1.0f + 2.0f * t0; };
const std::function<Amplitude (Second)> triangle = []( Second t ){ const Second t0 = std::fmod( t, 1.0f ); return t0 < 0.5f ? -1.0f + 4.0f * t0 : 3.0f - 4.0f * t0; };

}

}