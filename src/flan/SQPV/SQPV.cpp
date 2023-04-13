#include "SQPV.h"

#include <iostream>


using namespace flan;

SQPV::SQPV()
	: SQPVBuffer()
	{}

SQPV::SQPV( const Format & f )
	: SQPVBuffer( f )
	{}

SQPV::SQPV( SQPVBuffer && other )
	: SQPVBuffer( std::move( other ) )
	{}

SQPV SQPV::modifyFrequency( const Func2x1 & mod ) const
	{
	flan_PROCESS_START( SQPV() );

	SQPV out = copy();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		runtimeExecutionPolicyHandler( mod.getExecutionPolicy(), [&]( auto policy ){
		std::for_each( policy, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				out.getMF( channel, frame, bin ).f = mod( frameToTime( frame ), out.getMF( channel, frame, bin ).f );
			} ); } );
		}

	return out;
	}