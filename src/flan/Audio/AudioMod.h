#pragma once 

#include <functional>

#include "flan/defines.h"
#include "flan/Utility/execution.h"

namespace flan {

template<typename T>
class SoundMod
{
public:
	using FuncType = std::function<void ( T &, Second )>;

	SoundMod( const SoundMod & ) 				= delete;
	SoundMod & operator=( const SoundMod & ) 	= delete;
	SoundMod( SoundMod && ) 					= default;
	SoundMod & operator=( SoundMod && ) 		= default;
	~SoundMod() 								= default;

	template<typename C>
	requires 
#ifndef __APPLE__
		std::convertible_to<C, FuncType> && 
#endif
		(!std::same_as<C, std::nullptr_t>)
	SoundMod( C _f, ExecutionPolicy _p = ExecutionPolicy::Parallel_Unsequenced )
		: func( _f )
		, policy( _p )
		, null( false )
		{
		}	
	SoundMod();

	void operator()( T & in, Second t ) const;
	bool is_null() const;
	ExecutionPolicy get_execution_policy() const;

private:
	FuncType func;
	ExecutionPolicy policy;
	const bool null;
};

class Audio;
using AudioMod = SoundMod<Audio>;
class PV;
using PVMod = SoundMod<PV>;

}