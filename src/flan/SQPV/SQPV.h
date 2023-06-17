#pragma once

#include "SQPVBuffer.h"

#include "flan/Function.h"

namespace flan {

class Audio;
class Graph;

class SQPV : public SQPVBuffer
{
public:
	SQPV();
	SQPV( const Format & );
	SQPV( SQPVBuffer && );

	Audio convertToAudio() const;
	Audio convertToLeftRightAudio() const;
	Graph convertToGraph() const;

	SQPV modifyPitch( const Func2x1 & mod ) const;
	SQPV repitch( const Func2x1 & mod ) const;

	SQPV select( Second length, const Func2x2 & selector ) const;
};

}