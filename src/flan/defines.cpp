#include "defines.h"

using namespace flan;

Amplitude flan::decibel_to_amplitude( Decibel d ) { return std::pow( 10.0f, d / 20.0f ); }
Decibel flan::amplitude_to_decibel( Amplitude a ) { return 20.0f * std::log10( a ); }