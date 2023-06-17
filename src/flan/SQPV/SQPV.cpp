#include "SQPV.h"

#include <iostream>

#include "flan/Graph.h"

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

Graph SQPV::convertToGraph() const
	{
	// Temp
	Rect D;
	Pixel width = 1920;
	Pixel height = 1080; 
	float timelineScale = 0;

	flan_PROCESS_START( Graph( width, height ) );

	D.d.x1 = 0;
	D.d.x2 = getLength();
	D.r.x1 = getPitchBandwidth().first;
	D.r.x2 = getPitchBandwidth().second;

	// Validation, Conversion
	const Frame startFrame	= std::clamp( int( timeToFrame( D.x1() ) ), 0, getNumFrames() - 1 );
	const Frame endFrame	= std::clamp( int( timeToFrame( D.x2() ) ), 0, getNumFrames() - 1 );
	const Bin startBin		= std::clamp( int( pitchToBin( D.y1() ) ), 0, getNumBins() - 1 );
	const Bin endBin		= std::clamp( int( pitchToBin( D.y2() ) ), 0, getNumBins() - 1 );

	// Convert PV data to Value/Hue data
	const Magnitude maxMag = getMaxPartialMagnitude( startFrame, endFrame, startBin, endBin );
	std::vector<Func2x1> fs;
	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		fs.push_back( [channel, maxMag, this]( Second x, UnsignedPitch y )
			{
			const float m = getMP( channel, timeToFrame( x ), pitchToBin( y ) ).m;
			return std::sqrt( std::abs( m ) / maxMag ); // sqrt brings up dark areas
			} );

	Graph g( width, height );
	g.addFullSplitViewY( D, getNumChannels() ); 
	if( maxMag != 0 )
		g.drawSpectrograms( fs, D, 0 );
	
	// Tick drawing
	if( timelineScale > 0 )
		{
		const float bigTimeJump = std::pow( 4.0f, std::floor( std::log2( D.w() ) / 2 - 0.5f ) );
		g.drawXTicks( bigTimeJump / 4.0f, D.y2(), timelineScale / 2, 0, -1, Color::fromHSV( 0, 0, .6 ), 0 );
		g.drawXTicks( bigTimeJump, D.y2(), timelineScale, 0, -1, Color::fromHSV( 0, 0, 1 ), timelineScale );
		//g.drawYTicks( 2000, 0, 0, 14, -1, Color::fromHSV( 0, 0, 1 ), true );
		}

	return g;
	}

SQPV SQPV::modifyPitch( const Func2x1 & mod ) const
	{
	flan_PROCESS_START( SQPV() );

	SQPV out = copy();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		runtimeExecutionPolicyHandler( mod.getExecutionPolicy(), [&]( auto policy ){
		std::for_each( policy, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				//out.getMP( channel, frame, bin ).p = mod( frameToTime( frame ), out.getMP( channel, frame, bin ).p );
				}
			} ); } );
		}

	return out;
	}

SQPV SQPV::select( Second length, const Func2x2 & selector ) const
	{
	flan_PROCESS_START( SQPV() );

	// Input validation
	if( length <= 0 ) return SQPV();

	auto format = getFormat();
	format.numFrames = timeToFrame( length );
	SQPV out( format );

	// selectorSamples contains Second and Pitch
	const auto selectorSamples = selector.sample( 0, out.getNumFrames(), frameToTime( 1 ), 
		getPitchBandwidth().first * getBinsPerOctave(), 
		getPitchBandwidth().second * getBinsPerOctave(), 
		1.0f / getBinsPerOctave() );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		{
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < out.getNumBins(); ++bin )
				{
				const vec2 s = selectorSamples[ out.getBufferPos( 0, frame, bin ) ];
				const fFrame selectedFrame = timeToFrame( s.x() );
				const UnsignedPitch selectedPitch = s.y();
				const fBin selectedBin = pitchToBin( selectedPitch );

				if( selectedFrame < 0 || getNumFrames() - 1 <= selectedFrame || bin < 0 || getNumBins() - 1 <= bin )
				 	out.getMP( channel, frame, bin ) = { 0, 0 };
				else
					{	
					const fFrame mix = selectedFrame - std::floor( selectedFrame );
					const MP left_mp = getMP( channel, std::floor( selectedFrame ), bin );
					const MP right_mp = getMP( channel, std::floor( selectedFrame ) + 1, bin );

					const Magnitude w1 = ( 1.0f - mix ) * left_mp.m;
					const Magnitude w2 = mix * right_mp.m;

					MP assign_mp = left_mp;
					assign_mp.m = w1 + w2;
					assign_mp.p = w1 > w2 ? left_mp.p : right_mp.p;
					
					// const Pitch pitch_shift = binToPitch( bin ) - selectedPitch;
					// assign_mp.p += pitch_shift;
					out.getMP( channel, frame, bin ) = assign_mp;
					}
				}
			} );
		}

	return out;
	}
