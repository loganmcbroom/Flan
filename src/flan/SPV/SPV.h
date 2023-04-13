#pragma once

#include "SPVBuffer.h"

#include "flan/Function.h"
#include "flan/Utility/Interpolator.h"

namespace flan {

class SPV : public SPVBuffer
{
public:
	SPV();
	SPV( SPVBuffer && );
	SPV( Format );

	Audio convertToAudio();
	Audio convertToLeftRightAudio();

	SPV modifyFrequency( const Func2x1 & mod ) const;
	SPV repitch( const Func2x1 & mod ) const;
};

}