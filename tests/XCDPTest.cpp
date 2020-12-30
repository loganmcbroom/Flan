/* TODO:
	filter package - just do it in PVOC?
	continuity maintenance in wavelength/frequency finders
	octave error minimization in yin-based pitch finder (fast contour finder)
	audio to function of frequency
	convolve constant signals via fft
	Add \mainpage with examples somewhere
*/

#include <iostream>
#include <algorithm>
#include <numeric>

#include "xcdp/Audio.h"
#include "xcdp/PVOC.h" 
#include "xcdp/Function.h"
#include "xcdp/Wavetable.h"
#include "xcdp/Synthesis.h"
#include "xcdp/Graph.h"
#include "xcdp/FFTHelper.h"
#include "xcdp/WindowFunctions.h"

using namespace xcdp;

void play( const Audio & toPlay );
void graph( const std::string & f );
void graph( const bitmap_image & b );

//float normal_pdf( float x, float s )
//{
//    static const float inv_sqrt_2pi = 0.3989422804014327;
//    float a = x / s;
//
//    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
//}
//
//
//size_t getNumDigits( size_t x )  
//	{  
//    return 
//		(x < 1e1 ? 1 :   
//        (x < 1e2 ? 2 :   
//        (x < 1e3 ? 3 :   
//        (x < 1e4 ? 4 :   
//        (x < 1e5 ? 5 :   
//        (x < 1e6 ? 6 :   
//        (x < 1e7 ? 7 :  
//        (x < 1e8 ? 8 :  
//        (x < 1e9 ? 9 :  
//        (x < 1e10 ? 10 :  
//        (x < 1e11 ? 11 :  
//        (x < 1e12 ? 12 :  
//        (x < 1e13 ? 13 :  
//        (x < 1e14 ? 14 :  
//        (x < 1e15 ? 15 :  
//        (x < 1e16 ? 16 :  
//        (x < 1e17 ? 17 :  
//        18 )))))))))))))))));  
//	}  
//
//
//
//std::vector<int> getDigits( size_t n )
//	{
//	const size_t numDigits = getNumDigits( n );
//	std::vector<int> digits( numDigits );
//	size_t n_copy = n; //Need a copy to work with
//	for( size_t i = 0; i < numDigits; ++i )
//		{
//		digits[i] = n_copy % 10;
//		n_copy = std::floor( n_copy / 10 );
//		}
//	return digits;
//	}
//
//bool isPalindrome( size_t n )
//	{
//	auto digits = getDigits( n );
//	const int half = std::floor( digits.size() / 2 );
//	for( int i = 0; i < half; ++i )
//		{
//		if( digits[i] != digits[ digits.size() - 1 - i] )
//			return false;
//		}
//	return true;
//	}
//
//size_t kondo( const size_t n )
//	{
//	// Extract digits
//	auto digits = getDigits( n );
//
//	// Clean digits
//	std::sort( digits.begin(), digits.end() );
//
//	// Reassemble 
//	size_t m = 0;
//	for( size_t i = 0; i < digits.size(); ++i )
//		{
//		m *= 10;
//		m += digits[i];
//		}
//
//	return n + m;
//	}
//
//size_t kondoIter( size_t n, size_t limit )
//	{
//	while( n < limit ) 
//		{
//		n = kondo( n );
//		}
//	return n;
//	}
//
//void main()
//	{
//	const size_t numTests = 300;
//
//	for( size_t i = 1; i < numTests; ++i )
//		{
//		size_t n = i;
//		while( n < 1e17 ) 
//			{
//			n = kondo(n);
//			if( isPalindrome( n ) )
//				std::cout << i << " arrived at palindrome " << n << "\n";	
//			}
//		std::cout << i << "\n";
//		}
//	return;
//
//	//Second size_t to track number of times the end appears
//	std::vector<std::array< size_t, 2 >> ends;
//	//std::vector<float> ratios( numTests );
//
//	for( size_t i = 1; i < numTests; ++i )
//		{
//		size_t n = i;
//		while( n < 1e17 ) n = kondo( n );
//		auto end = std::find_if( ends.begin(), ends.end(), [n]( std::array< size_t, 2 > & a )
//			{
//			return a[0] == n;
//			} );
//		if( end == ends.end() )
//			ends.push_back( {n, 1} );
//		else
//			(*end)[1]++;
//			
//		//ratios[i] = float( ends.size() ) / float( i );
//		}
//
//	std::sort( ends.begin(), ends.end(), []( std::array< size_t, 2 > & a, std::array< size_t, 2 > & b )
//		{
//		return a[1] > b[1];
//		} );
//
//	std::cout << "\n\nTotal distinct endings: " << ends.size() << "\n\n";
//
//	for( int i = 0; i < ends.size(); ++i )
//		std::cout << ends[i][0] << " with " << ends[i][1] << " finishes\n";
//
//	//Func1x1 ratioFunc = [ratios]( float x )
//	//	{ 
//	//	if( x < 0 || ratios.size() <= x ) return 0.0f;
//	//	return ratios[ std::floor(x) ]; 
//	//	};
//
//	//ratioFunc.graph( "temp.bmp", float(numTests) * -0.1, numTests * 1.1, -.02, .1 );
//	//graph( "temp.bmp" );
//	}

#include "bmp/bitmap_image.hpp"

void main()
	{
	auto bah = Audio( "Audio/Bah.wav" );
	auto bai = Audio( "Audio/Bai.wav" );
	auto meow = Audio( "Audio/meow.wav" );
	auto slaw = Audio( "Audio/slaw.wav" );
	auto yeahh = Audio( "Audio/feel.wav" );
	auto sine = Synthesis::sine( 1, []( float t ){ return 1000; } );

	auto b = yeahh.convertToPVOC();//.desample( 200 );
	graph( b.convertToGraph() );
	play( b.convertToAudio() );

	}


#include <Windows.h>
#include <Mmsystem.h>
void play( const Audio & toPlay )
	{
	if( !toPlay.save( "TempFileSave.wav" ) )
        return;
	std::cout << "Playing sound ... \n";
	if( ! PlaySound("TempFileSave.wav", nullptr, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	}

void graph( const std::string & f )
	{
	system( (std::string("start ") + f).c_str() );
	}

void graph( const bitmap_image & b )
	{
	b.save_image( "Temp.bmp" );
	system( (std::string("start ") + "Temp.bmp").c_str() );
	}