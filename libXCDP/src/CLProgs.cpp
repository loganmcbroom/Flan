#include "xcdp/CLProgs.h"

namespace xcdp {
namespace CLProgs {

extern const char * Audio_convertToPVOC = R"(
__kernel void fftToPhase( global read_write float2 * fft )
	{
	const int i = get_global_id( 0 );
	const float phase = atan2( fft[i].y, fft[i].x );
	const float mag = hypot( fft[i].x, fft[i].y );
	fft[i].x = mag;
	fft[i].y = phase;
	}

__kernel void phaseToFreq( global read_write float2 * MagPhase, 
						   global write_only float2 * MagFreq,
						   float overlaps,
						   float binWidth,
						   int numBins )
	{
	const float pi = 3.141592653f;
	const float tau = 2.0f * pi;

	const int i = get_global_id( 0 );
	const float bin = i % numBins;
	
	float phaseDiff = MagPhase[i].y;
	if( i >= numBins )
		phaseDiff -= MagPhase[i-numBins].y;
	const float expectedPhaseDiff = bin * tau / overlaps;
	const float deltaPhase = phaseDiff - expectedPhaseDiff;
	const float wrappedDeltaPhase = deltaPhase - tau * floor( ( deltaPhase + pi ) / tau );
	const float binDeviation = overlaps * wrappedDeltaPhase / tau;
	const float frequency = ( bin + binDeviation ) * binWidth;
	
	MagFreq[i].x = MagPhase[i].x;
	MagFreq[i].y = frequency;
	}
	)";

extern const char * PVOC_convertToAudio = R"(
__kernel void phaseToFFT( global read_write float2 * MagPhase )
	{
	const int i = get_global_id( 0 );
	const float x = MagPhase[i].x * cos( MagPhase[i].y );
	const float y = MagPhase[i].x * sin( MagPhase[i].y );
	MagPhase[i].x = x;
	MagPhase[i].y = y;
	}

__kernel void freqToPhase( global read_write float2 * MagFreq,
						   global write_only float2 * MagPhase,
						   global read_write float * phaseAccum,
						   float overlaps,
						   float binWidth,
						   int numBins )
	{
	//load these from host
	//const float overlaps = 16.0f;
	//const float binWidth = 44100.0f / 2048.0f; //samplerate / framesize
	//const int numBins = 1025;
	
	const float pi = 3.141592653f;
	const float tau = 2.0f * pi;

	const int i = get_global_id( 0 );
	const int bin = i % numBins;
	const float binf = bin;
	
	const float binDeviation = MagFreq[i].y / binWidth - binf;
	const float deltaPhase = tau * binDeviation / overlaps;
	const float expectedPhaseDiff = binf * tau / overlaps;
	const float phaseDiff = deltaPhase + expectedPhaseDiff;
	phaseAccum[bin] += phaseDiff;
	
	MagPhase[i].x = MagFreq[i].x;
	MagPhase[i].y = phaseAccum[bin];
	}
	)";

extern const char * PVOC_modify = R"(
kernel void modify( global read_only float2 * in,
					global read_only float2 * modPointSamples,
					global write_only float2 * out,
					int numBins,
					int outNumFrames )
	{
	const int i = get_global_id( 0 );
	const int bin = i % numBins;
	if( bin != 0 && i >= numBins )
		{
		const int frame = floor( (float) i / numBins );
		
		/*
		p4 p3
		p1 p2
		*/
		
		const float2 p1 = modPointSamples[ ( frame - 1 ) * numBins + bin - 1 ];
		const float2 p2 = modPointSamples[ ( frame - 0 ) * numBins + bin - 1 ];
		const float2 p3 = modPointSamples[ ( frame - 0 ) * numBins + bin - 0 ];
		const float2 p4 = modPointSamples[ ( frame - 1 ) * numBins + bin - 0 ];
		
		const float2 p1MF = in[ ( frame - 1 ) * numBins + bin - 1 ];
		const float2 p2MF = in[ ( frame - 0 ) * numBins + bin - 1 ];
		const float2 p3MF = in[ ( frame - 0 ) * numBins + bin - 0 ];
		const float2 p4MF = in[ ( frame - 1 ) * numBins + bin - 0 ];
		
		const float2 D12 = p2 - p1;
		const float2 D23 = p3 - p2;
		const float2 D34 = p4 - p3;
		const float2 D41 = p1 - p4;
			
		// Find box containing shape to iterate over
		const int minx = fmax( floor( fmin( fmin( p1.x, p2.x ), fmin( p3.x, p4.x ) ) ), 0 );
		const int miny = fmax( floor( fmin( fmin( p1.y, p2.y ), fmin( p3.y, p4.y ) ) ), 0 );
		const int maxx = fmin( ceil(  fmax( fmax( p1.x, p2.x ), fmax( p3.x, p4.x ) ) ), outNumFrames - 1 );
		const int maxy = fmin( ceil(  fmax( fmax( p1.y, p2.y ), fmax( p3.y, p4.y ) ) ), numBins - 1 );
			
		// Iterate over box
		for( int x = minx; x <= maxx; ++x )
			for( int y = miny; y <= maxy; ++y )
				{	
				// Ref: http://paulbourke.net/geometry/polygonmesh/#insidepoly
				bool c = false;
				if( ( ( p1.y <= y && y < p4.y ) || ( p4.y <= y && y < p1.y ) ) && ( x < D41.x / D41.y * ( y - p1.y ) + p1.x ) ) c = !c;
				if( ( ( p2.y <= y && y < p1.y ) || ( p1.y <= y && y < p2.y ) ) && ( x < D12.x / D12.y * ( y - p2.y ) + p2.x ) ) c = !c;
				if( ( ( p3.y <= y && y < p2.y ) || ( p2.y <= y && y < p3.y ) ) && ( x < D23.x / D23.y * ( y - p3.y ) + p3.x ) ) c = !c;
				if( ( ( p4.y <= y && y < p3.y ) || ( p3.y <= y && y < p4.y ) ) && ( x < D34.x / D34.y * ( y - p4.y ) + p4.x ) ) c = !c;
				
				if( c )
					{
					// Ref: https://www.particleincell.com/2012/quad-interpolation/
					
					const float4 alpha = (float4)( p1.x, p2.x - p1.x, p4.x - p1.x, p1.x - p2.x + p3.x - p4.x );
					const float4 beta  = (float4)( p1.y, p2.y - p1.y, p4.y - p1.y, p1.y - p2.y + p3.y - p4.y );
					
					const float quadA = alpha.s3 * beta.s2 - alpha.s2 * beta.s3;
					const float quadB = alpha.s3 * beta.s0 - alpha.s0 * beta.s3 
									  + alpha.s1 * beta.s2 - alpha.s2 * beta.s1
									  + x 		 * beta.s3 - alpha.s3 * y;
					const float quadC = alpha.s1 * beta.s0 - alpha.s0 * beta.s1
									  + x		 * beta.s1 - alpha.s1 * y;
					
					float m;
					if( quadA == 0.0f )	
						{
						if( quadB == 0.0f )
							return;
						m = -quadC / quadB;
						}
					else
						{
						const float descriminant = quadB * quadB - 4.0f * quadA * quadC;
						if( descriminant < 0 ) return;
						m = ( -quadB + sqrt( descriminant ) ) / ( 2.0f * quadA );
						}
					if( alpha.s1 + alpha.s3 * m == 0 ) return;
					const float l = ( x - alpha.s0 - alpha.s2 * m ) / ( alpha.s1 + alpha.s3 * m );
					
					const float epsilon = 0.001f;
					if( fabs( l - 0.5f ) > 0.5f + epsilon || fabs( m - 0.5f ) > 0.5f + epsilon ) return;
					
					const float w1 = ( 1.0f - l ) * ( 1.0f - m ) * p1MF.x;
					const float w2 = ( 0.0f + l ) * ( 1.0f - m ) * p2MF.x;
					const float w3 = ( 0.0f + l ) * ( 0.0f + m ) * p3MF.x;
					const float w4 = ( 1.0f - l ) * ( 0.0f + m ) * p4MF.x;
					const float maxW = max( max(w1, w2), max(w3, w4) );
					
					out[ x * numBins + y ].x = w1 + w2 + w3 + w4;
						 if( maxW == w1 ) { out[ x * numBins + y ].y = p1MF.y; }
					else if( maxW == w2 ) { out[ x * numBins + y ].y = p2MF.y; }
					else if( maxW == w3 ) { out[ x * numBins + y ].y = p3MF.y; }
					else   				  { out[ x * numBins + y ].y = p4MF.y; }
					}
				}
		}
	}
	)";

extern const char * PVOC_modifyTime = R"(
kernel void modifyTime( global read_only float2 * in,
						  global read_only float * f,
						  global write_only float2 * out,
						  int numBins,
						  int numOutFrames )
	{
	const int i = get_global_id( 0 );
	if( i > numBins )
		{
		const int bin = i % numBins;
		
		const float2 MF = in[i];
		const float2 prevMF = in[ i - numBins ];
		
		const float outPos = f[i];
		const float prevOutPos = f[ i - numBins ];
		const bool forward = outPos > prevOutPos;

		const int startFrame = forward? ceil( prevOutPos ) : floor( prevOutPos );
		const int endFrame   = forward? ceil( outPos     ) : floor( outPos     );

		for( int outFrame = startFrame; outFrame != endFrame; forward? ++outFrame : --outFrame )
			{
			if( 0 <= outFrame ) // && outFrame < numOutFrames )
				{
				const float mix = ( (float) outFrame - prevOutPos ) / ( outPos - prevOutPos );
				out[outFrame * numBins + bin] = ( 1.0f - mix ) * prevMF + mix * MF;
				}
			}
		}
	}
	)";

extern const char * PVOC_modifyFrequency = R"(
kernel void modifyFrequency( global read_only float2 * in,
						  global read_only float * mappedFrequencies,
						  global read_only float * f,
						  global write_only float2 * out,
						  int numBins )
	{
	const int i = get_global_id( 0 );
	const int bin = i % numBins;
	if( bin != 0 )
		{
		const int frame = floor( (float) i / numBins );
			
		const float2 MF = in[i];
		const float2 prevMF = in[ i - 1 ];
		
		const float outBin = f[ i ];
		const float prevOutBin = f[ i - 1 ];
		const bool forward = outBin > prevOutBin;

		int startBin = forward? ceil( prevOutBin ) : floor( prevOutBin );
		int endBin   = forward? ceil( outBin     ) : floor( outBin     );
		startBin = clamp( startBin, 0, numBins-1 );
		endBin   = clamp( endBin,   0, numBins-1 );

		for( int writeBin = startBin; writeBin != endBin; forward? ++writeBin : --writeBin )
			{
			const float mix = ( (float) writeBin - prevOutBin ) / ( outBin - prevOutBin );
			const float interpMag = ( 1.0f - mix ) * prevMF.x + mix * MF.x;
			
			const size_t outIndex = frame * numBins + writeBin;

			out[outIndex].x = interpMag;
			out[outIndex].y = prevMF.x * ( 1.0 - mix ) > MF.x * mix ? 
				mappedFrequencies[i-1]:
				mappedFrequencies[i  ];
			}
		}
	}
	)";

};};