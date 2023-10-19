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

	Audio convert_to_audio();
	Audio convert_to_lr_audio();

	SPV modify_frequency( const Function<TF, Frequency> & mod ) const;
	SPV repitch( const Function<TF, Frequency> & mod ) const;
};

}