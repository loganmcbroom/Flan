#pragma once

#include "SpectrumBuffer.h"

namespace xcdp {

class Audio;

class Spectrum : public SpectrumBuffer
{
public:
	Spectrum( const Format f ) : SpectrumBuffer( f ) {}

	Audio convertToAudio();

	//static Spectrum 

	//Spectrum multiply( const Spectrum & filter );
};

}