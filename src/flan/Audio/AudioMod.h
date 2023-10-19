#pragma once 

#include <functional>

#include "flan/defines.h"
#include "flan/Utility/execution.h"

namespace flan {

class Audio;

class AudioMod
{
public:
	using FuncType = std::function<void ( Audio &, Second )>;

	template<typename C>
	requires std::convertible_to<C, FuncType> && (!std::same_as<C, std::nullptr_t>)
	AudioMod( C _f, ExecutionPolicy _p = ExecutionPolicy::Parallel_Unsequenced )
		: f( _f )
		, policy( _p )
		, null( false )
		{
		}	
	AudioMod();
	void operator()( Audio & in, Second t ) const;
	bool is_null() const;
	ExecutionPolicy get_execution_policy() const;

private:
	FuncType f;
	ExecutionPolicy policy;
	const bool null;
};

}