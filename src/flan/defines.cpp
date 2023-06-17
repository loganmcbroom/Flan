#include "defines.h"

using namespace flan;

Amplitude flan::decibelToAmplitude( Decibel d ) { return std::pow( 10.0f, d / 20.0f ); }
Decibel flan::amplitudeToDecibel( Amplitude a ) { return 20.0f * std::log10( a ); }