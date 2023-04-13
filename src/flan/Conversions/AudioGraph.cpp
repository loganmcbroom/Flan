#include "flan/Audio/Audio.h"

#include "flan/Graph.h"

#include "bmp/bitmap_image.hpp"

using namespace flan;

Graph Audio::convertToGraph( Interval I, Pixel width, Pixel height, float timelineScale ) const
	{
	flan_PROCESS_START( Graph( width, height ) );

	if( I.x2 == -1 ) I.x2 = getLength();
	
	const Frame startFrame = timeToFrame( I.x1 );
	const Frame endFrame = timeToFrame( I.x2 );

	Graph g( width, height );
	g.fillImage( Color::fromHSV( 0, 0, .04 ) );
	g.addFullSplitViewY( I * Interval( -1, 1 ), getNumChannels() );
		
	std::vector<const float *> channelStarts( getNumChannels() );
	for( int i = 0; i < getNumChannels(); ++i )
		channelStarts[i] = getSamplePointer( i, 0 );
	g.drawWaveforms( channelStarts, getNumFrames(), { 0, -1, getLength(), 1 }, 0, Graph::WaveformMode::Symmetric );

	if( timelineScale > 0 )
		{
		const float bigJump = std::pow( 4.0f, std::floor( std::log2( I.w() ) / 2 - 0.5f ) );
		g.drawXTicks( bigJump / 4.0f, 1, timelineScale / 2, 0, -1, Color::fromHSV( 0, 0, .6 ), 0 );
		g.drawXTicks( bigJump, 1, timelineScale, 0, -1, Color::fromHSV( 0, 0, 1 ), timelineScale );
		}

	return g;
	}

void Audio::saveToBMP( const std::string & filename, Interval I, int width, int height ) const
	{
	flan_FUNCTION_LOG;
	auto b = convertToGraph( I, width, height );
	b.save_image( filename );
	}