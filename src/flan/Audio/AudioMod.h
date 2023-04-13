#pragma once 

#include <functional>

#include "flan/Utility/execution.h"

namespace flan {

class Audio;

class AudioMod
{
public:
	using FuncType = std::function<Audio ( const Audio &, float )>;

	template<typename C>
	requires std::convertible_to<C, FuncType> && (!std::same_as<C, std::nullptr_t>)
	AudioMod( C _f, ExecutionPolicy _p = ExecutionPolicy::Parallel_Unsequenced )
		: f( _f )
		, policy( _p )
		, valid( true )
		{
		}

	AudioMod();
	Audio operator()( const Audio & in, float t ) const;
	AudioMod operator()( AudioMod other ) const;
	AudioMod operator|( AudioMod other ) const;
	operator bool() const;
	bool isValid() const;
	ExecutionPolicy getExecutionPolicy() const;

private:
	FuncType f;
	ExecutionPolicy policy;
	const bool valid;
};

}