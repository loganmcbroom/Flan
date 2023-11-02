#pragma once

#include <functional>

#include "flan/defines.h"
#include "flan/Utility/execution.h"

namespace flan {

/** Function type exclusive to PV::prism. 
*  The identity PrismFunc is []( int n, Second t, int h, Frequency f, const std::vector<Magnitude> & hM ){ return MF{ hM[h], f }; }
*/
class PrismFunc {
public:
	using FuncType = std::function<MF ( Index note, Second, Harmonic harmonic, Frequency contourFreq, const std::vector<Magnitude> & harmonicMagnitudes )>;

	PrismFunc( const PrismFunc & ) 				= delete;
	PrismFunc & operator=( const PrismFunc & ) 	= delete;
	PrismFunc( PrismFunc && ) 					= default;
	PrismFunc & operator=( PrismFunc && ) 		= default;
	~PrismFunc() 								= default;

	template<typename C>
	requires 
#ifndef __APPLE__
		std::convertible_to<C, FuncType> && 
#endif
		(!std::same_as<C, std::nullptr_t>)
	PrismFunc( C f_, ExecutionPolicy p_ = ExecutionPolicy::Parallel_Unsequenced ) 
		: func( f_ ) 
		, policy( p_ )
		, null( false )
		{}
	PrismFunc();

	MF operator()( 
		int n, 
		Second t, 
		int h, 
		Frequency f, 
		const std::vector<Magnitude> & hM 
		) const;

	bool is_null() const;

	ExecutionPolicy get_execution_policy() const;

private:
	FuncType func;
	ExecutionPolicy policy;
	const bool null;
};

}