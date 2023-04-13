#include "AudioMod.h"

#include "flan/Audio/Audio.h"

using namespace flan;

AudioMod::AudioMod()
	: f( nullptr )
	, policy( ExecutionPolicy::Parallel_Unsequenced )
	, valid( false )
	{
	}

Audio AudioMod::operator()( const Audio & in, float t ) const
	{
	if( ! isValid() ) return Audio();

	return f( in, t );
	}

AudioMod AudioMod::operator()( AudioMod other ) const
	{
	if( ! isValid() || ! other.isValid() ) return AudioMod();

	return AudioMod( [this, other = std::move(other)]( const Audio & in, float t )
		{
		return operator()( other( in, t ), t );
		}, lowestExecution( getExecutionPolicy(), other.getExecutionPolicy() ) );
	}

AudioMod AudioMod::operator|( AudioMod other ) const
	{
	return operator()( std::move( other ) );
	}

AudioMod::operator bool() const { return valid; }
bool AudioMod::isValid() const { return valid; }

ExecutionPolicy AudioMod::getExecutionPolicy() const { return policy; }