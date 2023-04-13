#pragma once

#include "SQPVBuffer.h"

#include "flan/Function.h"

namespace flan {

class Audio;

class SQPV : public SQPVBuffer
{
public:
	SQPV();
	SQPV( const Format & );
	SQPV( SQPVBuffer && );

	Audio convertToAudio() const;
	Audio convertToLeftRightAudio() const;

	SQPV modifyFrequency( const Func2x1 & mod ) const;
};

}