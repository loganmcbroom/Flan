#include "AudioMod.h"

#include "flan/Audio/Audio.h"

using namespace flan;

template<typename T>
SoundMod<T>::SoundMod()
	: func( nullptr )
	, policy( ExecutionPolicy::Parallel_Unsequenced )
	, null( true )
	{
	}

template<typename T>
void SoundMod<T>::operator()( T & in, float t ) const
	{
	if( is_null() ) return;
	func( in, t );
	}

template<typename T>
bool SoundMod<T>::is_null() const 
	{ 
	return null; 
	}

template<typename T>
ExecutionPolicy SoundMod<T>::get_execution_policy() const 
	{ 
	return policy;
	}

template class SoundMod<Audio>;
template class SoundMod<PV>;