#include <iostream>
#include <fstream>
#include <string>
#include <Windows.h>
#include <Mmsystem.h>

#include <sndfile.h>

#include "Audio.h"
#include "PVOC.h"

using namespace xcdp;

void spectrographFinancial( const std::string & fileName )
	{
	//open data
	std::ifstream spyFile( fileName );
	if( ! spyFile.is_open() )
		{
		std::cout << "Error opening " << fileName <<"\n";
		return;
		}
	const size_t numLines = std::count(std::istreambuf_iterator<char>(spyFile), 
             std::istreambuf_iterator<char>(), '\n');
	spyFile.clear();                 // clear fail and eof bits
	spyFile.seekg(0, std::ios::beg); // back to the start!

	std::string lineString( 128, 'a' );
	std::getline( spyFile, lineString );

	//set up "Audio" buffer
	Audio financial;
	financial.setBufferSize( 1, numLines );
	size_t line = 0;
	while( std::getline( spyFile, lineString ) )
		{
		size_t highStartPos = 0;
		for( size_t i = 0; i < 3; ++i )
			highStartPos = lineString.find( '	', highStartPos+1 );
		size_t highEndPos = lineString.find( '	', highStartPos+1 );
		size_t lowEndPos = lineString.find( '	', highEndPos+1 );

		double high = std::stod( lineString.substr(highStartPos, highEndPos-highStartPos ) );
		double low  = std::stod( lineString.substr(highEndPos, lowEndPos-highEndPos ) );

		financial.buffer[0][line] = (high-low)/2.0;

		++line;
		}
	financial.getPVOC( 512, 128 ).getSpectrograph();

	system( "\"C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe\" spectrograph.tga" );
	}

int main( int argc, char *argv[] )
	{
	Audio( "bah.wav" )
	.iterate( 3, []( const Audio & in, size_t n ){ return in.repitch( n+1 ); } )
	.getPVOC( 1024, 4 )
	//.getSpectrograph();
	.repitch( .8 )
	.timeAverage( 50 )
	.getAudio()
	.setVolume( .8 )
	.save( "test.wav" );

	//system( "\"C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe\" spectrograph.tga" );

	std::cout << "Playing sound ... ";
	if( ! PlaySound("test.wav", NULL, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	std::cout << "Done\n";

	return 0;
	}