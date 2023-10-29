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

SPV SPV::modify_frequency( const Function<TF, Second> & mod ) const
	{
	if( is_null() ) return SPV();

	SPV out = copy();

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		runtime_execution_policy_handler( mod.get_execution_policy(), [&]( auto policy ){
		std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
			{
			const Second time = frame_to_time( frame );
			for( Bin bin = 0; bin < get_num_bins(); ++bin )
				out.get_MF( channel, frame, bin ).f = mod( TF{time, out.get_MF( channel, frame, bin ).f} );
			} ); } );
		}

	return out;
	}

SPV SPV::repitch( const Function<TF, Frequency> & mod ) const
	{
	return modify_frequency( [&]( TF tf ){ return tf.f * mod( tf ); } );
	}