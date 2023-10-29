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

	AudioMod( const AudioMod & ) 				= delete;
	AudioMod & operator=( const AudioMod & ) 	= delete;
	AudioMod( AudioMod && ) 					= default;
	AudioMod & operator=( AudioMod && ) 		= default;
	~AudioMod() 								= default;

	template<typename C>
	requires 
#ifndef __APPLE__
		std::convertible_to<C, FuncType> && 
#endif
		(!std::same_as<C, std::nullptr_t>)
	AudioMod( C _f, ExecutionPolicy _p = ExecutionPolicy::Parallel_Unsequenced )
		: func( _f )
		, policy( _p )
		, null( false )
		{
		}	
	AudioMod();

	void operator()( Audio & in, Second t ) const;
	bool is_null() const;
	ExecutionPolicy get_execution_policy() const;

private:
	FuncType func;
	ExecutionPolicy policy;
	const bool null;
};

}