#include "flan/PV/PV.h"
#include "flan/Graph.h"

using namespace flan;

Graph PV::convertToGraph( Rect D, Pixel width, Pixel height, float timelineScale ) const
	{
	flan_PROCESS_START( Graph( width, height ) );

	if( D.x2() == -1 ) D.d.x2 = getLength();
	if( D.y2() == -1 ) D.r.x2 = getHeight();

	// Validation, Conversion
	const float startFrame	= std::clamp( int( D.x1() * timeToFrame() ), 0, getNumFrames() - 1 );
	const float endFrame	= std::clamp( int( D.x2() * timeToFrame() ), 0, getNumFrames() - 1 );
	const float startBin	= std::clamp( int( D.y1() * frequencyToBin() ), 0, getNumBins() - 1 );
	const float endBin		= std::clamp( int( D.y2() * frequencyToBin() ), 0, getNumBins() - 1 );

	// Convert PV data to Value/Hue data
	const float maxMag = getMaxPartialMagnitude( startFrame, endFrame, startBin, endBin );
	std::vector<Func2x1> fs;
	for( int i = 0; i < getNumChannels(); ++i )
		fs.push_back( [this, maxMag, i]( float x, float y )
			{
			const float m = getMF( i, x * timeToFrame(), y * frequencyToBin() ).m;
			return std::sqrt( std::abs( m ) / maxMag ) * log2( 2.0f + y ) / 4.0f; // sqrt brings up dark areas, log scaling brings up high frequencies
			} );

	Graph g( width, height );
	g.addFullSplitViewY( D, getNumChannels() ); 
	if( maxMag != 0 )
		g.drawSpectrograms( fs, { 0, 0, getLength(), getHeight() }, 0 );
	
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

const PV & PV::saveToBMP( const std::string & fileName, Rect D, Pixel width, Pixel height ) const
	{
	flan_FUNCTION_LOG;
	convertToGraph( D, width, height ).save_image( fileName );
	return *this;
	}