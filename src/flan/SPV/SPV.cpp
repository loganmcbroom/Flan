#include "SPV.h"

#include <iostream>

#include "flan/Utility/execution.h"

using namespace flan;

SPV::SPV()
	: SPVBuffer()
	{}

SPV::SPV( SPVBuffer && other )
	: SPVBuffer( std::move( other ) )
	{}

SPV::SPV( Format _format )
	: SPVBuffer( _format )
	{}

SPV SPV::modifyFrequency( const Func2x1 & mod ) const
	{
	flan_PROCESS_START( SPV() );

	SPV out = copy();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		runtimeExecutionPolicyHandler( mod.getExecutionPolicy(), [&]( auto policy ){
		std::for_each( policy, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			const Second time = frameToTime( frame );
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				out.getMF( channel, frame, bin ).f = mod( time, out.getMF( channel, frame, bin ).f );
			} ); } );
		}

	return out;
	}

SPV SPV::repitch( const Func2x1 & mod ) const
	{
	return modifyFrequency( [&]( vec2 tf ){ return tf.f() * mod( tf ); } );
	}