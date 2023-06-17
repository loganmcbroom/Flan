#pragma once

#include <functional>

#include "flan/PV/PVBuffer.h"
#include "flan/Utility/execution.h"

namespace flan {

/** Function type exclusive to PV::prism. 
*  The identity PrismFunc is []( int n, Second t, int h, Frequency f, const std::vector<Magnitude> & hM ){ return MF{ hM[h], f }; }
*/
class PrismFunc {
public:
	using FuncType = std::function<MF ( Index note, Second, Index harmonic, Frequency contourFreq, const std::vector<Magnitude> & harmonicMagnitudes )>;

	template<typename C>
	requires std::convertible_to<C, FuncType>
	PrismFunc( C f_, ExecutionPolicy p_ = ExecutionPolicy::Parallel_Unsequenced ) 
		: func( f_ ) 
		, policy( p_ )
		{}

	MF operator()( int n, Second t, int h, Frequency f, const std::vector<Magnitude> & hM ) const { return func(n, t, h, f, hM); }
	ExecutionPolicy getExecutionPolicy() const { return policy; }

private:
	FuncType func;
	ExecutionPolicy policy;
};

}