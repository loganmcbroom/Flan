#include "AudioMod.h"

#include "flan/Audio/Audio.h"

using namespace flan;

AudioMod::AudioMod()
	: func( nullptr )
	, policy( ExecutionPolicy::Parallel_Unsequenced )
	, null( true )
	{
	}

void AudioMod::operator()( Audio & in, float t ) const
	{
	if( is_null() ) return;
	func( in, t );
	}

bool AudioMod::is_null() const 
	{ 
	return null; 
	}

ExecutionPolicy AudioMod::get_execution_policy() const 
	{ 
	return policy;
	}