#include "xcdp/PVOC.h"

#include <iostream>

#include "xcdp/spline.h"
#include "xcdp/CLContext.h"
#include "xcdp/CLProgs.h"

namespace xcdp {

PVOC PVOC::modify( Func2x2 mod, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	const Bin numBins = getNumBins();
	const Frame numFrames = getNumFrames();
	const uint32_t inChannelDataCount = numFrames * numBins;

	//Sample mod functions
	std::vector<MF> inModified( inChannelDataCount );
	std::vector<vec2> modPointSamples( inChannelDataCount );
	float lastOutputFrame = 0;
	for( Frame frame = 0; frame < numFrames; ++frame )
		for( Bin bin = 0; bin < numBins; ++bin )
			{
			XCDP_CANCEL_POINT( PVOC() );
			vec2 modOut = mod( frame * frameToTime(), bin * binToFrequency() );
			modOut.x() *= timeToFrame();
			modOut.y() *= frequencyToBin();
			
			modPointSamples[ frame * numBins + bin ] = modOut;
			if( lastOutputFrame < modOut.x() )
				lastOutputFrame = modOut.x();
				
			}
	++lastOutputFrame;
	if( lastOutputFrame * frameToTime() > 600.0f ) // Outfile longer than 10 minutes?
		{
		std::cout << "PVOC::modify_cpu tried to make a file longer than 10 minutes, which is currently disabled";
		return PVOC();
		}

	auto format = getFormat();
	format.numFrames = ceil( lastOutputFrame );
	PVOC out( format );
	out.clearBuffer();
	const uint32_t outChannelDataCount = out.getNumFrames() * out.getNumBins();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Calculate and copy mapped frequencies and input magnitudes
		auto modDataSamplesWriteHead = inModified.begin();
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				{
				*modDataSamplesWriteHead = 
					{
					inBufferReadHead->m,
					mod( frame * frameToTime(), inBufferReadHead->f ).y()
					};
				}
		
		for( Frame frame = 1; frame < numFrames; ++frame )
			for( Bin bin = 1; bin < numBins; ++bin )
				{
				const std::array<vec2, 4> p = {
					modPointSamples[ ( frame - 1 ) * numBins + bin - 1 ],
					modPointSamples[ ( frame - 0 ) * numBins + bin - 1 ],
					modPointSamples[ ( frame - 0 ) * numBins + bin - 0 ],
					modPointSamples[ ( frame - 1 ) * numBins + bin - 0 ] };
		
				const std::array<MF, 4> pMF = {
					inModified[ ( frame - 1 ) * numBins + bin - 1 ],
					inModified[ ( frame - 0 ) * numBins + bin - 1 ],
					inModified[ ( frame - 0 ) * numBins + bin - 0 ],
					inModified[ ( frame - 1 ) * numBins + bin - 0 ] };
		
				const vec2 D12 = p[1] - p[0];
				const vec2 D23 = p[2] - p[1];
				const vec2 D34 = p[3] - p[2];
				const vec2 D41 = p[0] - p[3];
			
				// Find box containing quad to iterate over
				const int minx = fmax( floor( fmin( fmin( p[0].x(), p[1].x() ), fmin( p[2].x(), p[3].x() ) ) ), 0 );
				const int miny = fmax( floor( fmin( fmin( p[0].y(), p[1].y() ), fmin( p[2].y(), p[3].y() ) ) ), 0 );
				const int maxx = fmin( ceil(  fmax( fmax( p[0].x(), p[1].x() ), fmax( p[2].x(), p[3].x() ) ) ), format.numFrames - 1 );
				const int maxy = fmin( ceil(  fmax( fmax( p[0].y(), p[1].y() ), fmax( p[2].y(), p[3].y() ) ) ), numBins - 1 );
			
				// Iterate over box
				for( int x = minx; x <= maxx; ++x )
					for( int y = miny; y <= maxy; ++y )
						{	
						XCDP_CANCEL_POINT( PVOC() );
						// Test if (x,y) is inside the quad
						// Ref: http://paulbourke.net/geometry/polygonmesh/#insidepoly
						bool c = false;
						if( ( ( p[0].y() <= y && y < p[3].y() ) || ( p[3].y() <= y && y < p[0].y() ) ) && ( x < D41.x() / D41.y() * ( y - p[0].y() ) + p[0].x() ) ) c = !c;
						if( ( ( p[1].y() <= y && y < p[0].y() ) || ( p[0].y() <= y && y < p[1].y() ) ) && ( x < D12.x() / D12.y() * ( y - p[1].y() ) + p[1].x() ) ) c = !c;
						if( ( ( p[2].y() <= y && y < p[1].y() ) || ( p[1].y() <= y && y < p[2].y() ) ) && ( x < D23.x() / D23.y() * ( y - p[2].y() ) + p[2].x() ) ) c = !c;
						if( ( ( p[3].y() <= y && y < p[2].y() ) || ( p[2].y() <= y && y < p[3].y() ) ) && ( x < D34.x() / D34.y() * ( y - p[3].y() ) + p[3].x() ) ) c = !c;
				
						if( c )
							{
							// Quadrilateral interpolation
							// Ref: https://www.particleincell.com/2012/quad-interpolation/
					
							const std::array<float, 4> alpha = { p[0].x(), p[1].x() - p[0].x(), p[3].x() - p[0].x(), p[0].x() - p[1].x() + p[2].x() - p[3].x() };
							const std::array<float, 4> beta  = { p[0].y(), p[1].y() - p[0].y(), p[3].y() - p[0].y(), p[0].y() - p[1].y() + p[2].y() - p[3].y() };
					
							const float quadA = alpha[3] * beta[2] - alpha[2] * beta[3];
							const float quadB = alpha[3] * beta[0] - alpha[0] * beta[3] 
											  + alpha[1] * beta[2] - alpha[2] * beta[1]
											  + x 		 * beta[3] - alpha[3] * y;
							const float quadC = alpha[1] * beta[0] - alpha[0] * beta[1]
												+ x		 * beta[1] - alpha[1] * y;
					
							float m;
							if( quadA == 0.0f )	
								{
								if( quadB == 0.0f )
									continue;
								m = -quadC / quadB;
								}
							else
								{
								const float descriminant = quadB * quadB - 4.0f * quadA * quadC;
								if( descriminant < 0 ) continue;
								m = ( -quadB + sqrt( descriminant ) ) / ( 2.0f * quadA );
								}
							if( alpha[1] + alpha[3] * m == 0 ) continue;
							const float l = ( x - alpha[0] - alpha[2] * m ) / ( alpha[1] + alpha[3] * m );
					
							const float epsilon = 0.0001f;
							if( fabs( l - 0.5f ) > 0.5f + epsilon || fabs( m - 0.5f ) > 0.5f + epsilon ) continue;

							const float interpL = interp( l );
							const float interpM = interp( m );

							const std::array<float, 4> w = {
								( 1.0f - interpL ) * ( 1.0f - interpM ) * pMF[0].m,
								(        interpL ) * ( 1.0f - interpM ) * pMF[1].m,
								(        interpL ) * (        interpM ) * pMF[2].m,
								( 1.0f - interpL ) * (        interpM ) * pMF[3].m };
							const float totalWeight = w[0] + w[1] + w[2] + w[3];
							if( totalWeight <= 0 ) continue;
					
							auto & outMF = out.getMF( channel, x, y );

							// Frequency is found by taking a weighted sum of the contributing frequencies of the current step 
							//	and what is already in the output MF
							outMF.f = ( outMF.f * outMF.m 
											  + pMF[0].f * w[0]
											  + pMF[1].f * w[1]
											  + pMF[2].f * w[2]
											  + pMF[3].f * w[3] ) 
											  / ( outMF.m + totalWeight );
							outMF.m += totalWeight;
							}
						}
				}
		}

	return out;
	}

PVOC PVOC::modifyFrequency( Func2x1 outFreqFunc, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	PVOC out( getFormat() );
	out.clearBuffer();

	const uint32_t numChannels = getNumChannels();
	const uint32_t numFrames = getNumFrames();
	const uint32_t numBins = getNumBins();

	for( Channel channel = 0; channel < numChannels; ++channel )
		for( Frame frame = 0; frame < numFrames; ++frame )
			{
			float prevOutBin = frequencyToBin() * outFreqFunc( frameToTime() * frame, 0.0f );
			for( Bin bin = 1; bin < numBins; ++bin )
				{
				const float outBin = frequencyToBin() * outFreqFunc( frameToTime() * frame, binToFrequency() * bin );
				const bool forward = outBin > prevOutBin;

				const long long prevOutBinRound = forward? std::ceil( prevOutBin ) : std::floor( prevOutBin );
				const long long outBinRound     = forward? std::ceil( outBin     ) : std::floor( outBin     );
				const uint32_t startBin = std::clamp( prevOutBinRound, 0ll, long long(out.getNumBins() - 1) );
				const uint32_t endBin   = std::clamp( outBinRound    , 0ll, long long(out.getNumBins() - 1) );

				for( uint32_t interpBin = startBin; interpBin != endBin; forward? ++interpBin : --interpBin )
					{
					XCDP_CANCEL_POINT( PVOC() );
					if( interpBin < 0 || getNumBins() <= interpBin ) continue;

					const MF & currentMF = getMF( channel, frame, bin );
					const MF & prevMF = getMF( channel, frame, bin-1 );

					const float mix = interp( ( float( interpBin ) - prevOutBin ) / ( outBin - prevOutBin ) );
					const float interpMagnitude = ( 1.0f - mix ) * prevMF.m + mix * currentMF.m;

					auto & outMF = out.getMF( channel, frame, interpBin );

					//if( prevMF.m * ( 1.0 - mix ) > currentMF.m * ( 0.0 + mix ) )
					//	outMF.f = outFreqFunc( frameToTime( frame ), prevMF.f );
					//else
					//	outMF.f = outFreqFunc( frameToTime( frame ), currentMF.f );

					const float w0 = ( 1.0f - mix ) * prevMF.m;
					const float w1 = (        mix ) * currentMF.m;

					outMF.f = ( outMF.m * outMF.f + w0 * prevMF.f + w1 * currentMF.f ) / ( outMF.m + w0 + w1 );
					outMF.m += interpMagnitude;
					}
				prevOutBin = outBin;
				}
			}

	return out;
	}

PVOC PVOC::modifyTime( Func2x1 outPosFunc, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	const uint32_t numChannels = getNumChannels();
	const uint32_t numFrames = getNumFrames();
	const uint32_t numBins = getNumBins();

	//Find farthest output frame and sample mod function
	std::vector<float> modFrames( numFrames * numBins );
	float lastOutputFrame = 0;
	for( Frame frame = 0; frame < getNumFrames(); ++frame )
		for( Bin bin = 0; bin < getNumBins(); ++bin )
			{
			XCDP_CANCEL_POINT( PVOC() );
			const float modFrame = outPosFunc( frame * frameToTime(), bin * binToFrequency() ) * timeToFrame();
			if( lastOutputFrame < modFrame )
				lastOutputFrame = modFrame;
			modFrames[ frame * numBins + bin ] = modFrame;
			}

	auto format = getFormat();
	format.numFrames = lastOutputFrame;
	PVOC out( format );
	out.clearBuffer();

	for( Channel channel = 0; channel < numChannels; ++channel )
		{
		for( Frame frame = 1; frame < numFrames; ++frame ) // Start at 1 because for each frame we will access previous frame
			{
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				const float prevModFrame = modFrames[ ( frame - 1 ) * numBins + bin ];
				const float modFrame = modFrames[ frame * numBins + bin ];
				const bool forward = modFrame > prevModFrame;

				const long long startFrame = forward? std::ceil( prevModFrame ) : std::floor( prevModFrame );
				const long long endFrame   = forward? std::ceil( modFrame     ) : std::floor( modFrame     );

				for( long long outFrame = startFrame; outFrame != endFrame; forward? ++outFrame : --outFrame )
					{
					XCDP_CANCEL_POINT( PVOC() );
					if( outFrame < 0 || out.getNumFrames() <= outFrame ) continue;

					const MF & currentMF = getMF( channel, frame, bin );
					const MF & prevMF = getMF( channel, frame - 1, bin );

					const float mix = interp( ( outFrame - prevModFrame ) / ( modFrame - prevModFrame ) );
					const float w0 = ( 1.0f -  mix ) * prevMF.m;
					const float w1 = mix * currentMF.m;

					//change to new freq technique
					auto & outMF = out.getMF( channel, outFrame, bin );
					outMF.f = ( outMF.f * outMF.m + w0 * prevMF.f + w1 * currentMF.f ) / ( outMF.m + w0 + w1 );
					outMF.m += w0 + w1;
					}
				}
			}
		}

	return out;
	}

#ifdef USE_OPENCL

PVOC PVOC::modify_cl( Func2x2 mod, Interpolator interp, XCDP_CANCEL_ARG_CPP  ) const
	{
	XCDP_PROCESS_START( PVOC() );

	using float2 = std::array<float,2>;

	const uint32_t numBins = getNumBins();
	const int numBins_i = numBins;
	const uint32_t numFrames = getNumFrames();
	const uint32_t inChannelDataCount = numFrames * numBins;

	//Sample mod functions
	std::vector<float2> inModified( inChannelDataCount );
	std::vector<vec2> modPointSamples( inChannelDataCount );
	float lastOutputFrame = 0;
	for( Frame frame = 0; frame < numFrames; ++frame )
		for( Bin bin = 0; bin < numBins; ++bin )
			{
			XCDP_CANCEL_POINT( PVOC() );
			auto modOut = mod( frame * frameToTime(), bin * binToFrequency() );
			modOut.x() *= timeToFrame();
			modOut.y() *= frequencyToBin();
			
			modPointSamples[ frame * numBins + bin ] = modOut;
			if( lastOutputFrame < modOut.x() )
				{
				//std::cout << modOut.x() << std::endl;
				lastOutputFrame = modOut.x();
				}
			}
	++lastOutputFrame;
	if( lastOutputFrame * frameToTime() > 600.0f ) // outfile longer than 10 minutes?
		{
		std::cout << "PVOC::modify tried to make a file longer than 10 minutes, which is currently disabled";
		return PVOC();
		}

	auto format = getFormat();
	format.numFrames = ceil( lastOutputFrame );
	PVOC out( format );
	//out.clearBuffer();
	const uint32_t outChannelDataCount = out.getNumFrames() * out.getNumBins();
	const int outNumFrames_i = format.numFrames;

	//Prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_modify );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::Buffer clModPointSamples( cl.context, CL_MEM_READ_ONLY, sizeof( float2 ) * inChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, int, int > 
		cl_modify( programHelper.program, "modify" );

	//Copy stretch samples to device
	cl::copy( cl.queue, modPointSamples.begin(), modPointSamples.end(), clModPointSamples );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Clear cl out buffer ( the entire buffer might not be written to so it must be scrubbed )
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList(); //Can't clear until out copy from previous iteration is complete
		cl.queue.enqueueFillBuffer( clOut, 0, 0, sizeof( PVOCBuffer::MF ) * outChannelDataCount );

		//Calculate and copy mapped frequencies and input magnitudes
		auto modDataSamplesWriteHead = inModified.begin();
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				{
				XCDP_CANCEL_POINT( PVOC() );
				*modDataSamplesWriteHead = 
					{
					inBufferReadHead->m,
					mod( frame * frameToTime(), inBufferReadHead->f ).y()
					};
				}
		cl::copy( cl.queue, inModified.begin(), inModified.end(), clIn );
		
		//Compute on device
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl_modify( cl::EnqueueArgs( cl.queue, cl::NDRange( inChannelDataCount ) ), 
			clIn,
			clModPointSamples, 
			clOut, 
			numBins_i,
			outNumFrames_i );
			
		//copy cl output buffer to out
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut,
			out.getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * outChannelDataCount, 
			out.getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * outChannelDataCount );
		}

	return out;
	}

PVOC PVOC::modifyFrequency_cl( Func2x1 outFreqFunc, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	PVOC out( getFormat() );
	out.clearBuffer();

	const uint32_t numBins = getNumBins();
	const int numBins_i = numBins;
	const uint32_t numFrames = getNumFrames();
	const uint32_t inChannelDataCount = numFrames * numBins;

	//Sample stretch function
	std::vector<float> outBinBuffer( inChannelDataCount );
	std::vector<float> mappedFrequencies( inChannelDataCount );
	for( Frame frame = 0; frame < numFrames; ++frame )
		for( Bin bin = 0; bin < numBins; ++bin )
			outBinBuffer[ frame * numBins + bin ] = outFreqFunc( frame * frameToTime(), bin * binToFrequency() ) * frequencyToBin();

	//Prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_modifyFrequency );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::Buffer clMappedFrequencies( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * inChannelDataCount );
	cl::Buffer clF( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * inChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, int > 
		cl_modifyFrequency( programHelper.program, "modifyFrequency" );

	//Copy stretch samples to device
	cl::copy( cl.queue, outBinBuffer.begin(), outBinBuffer.end(), clF );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Clear cl out buffer ( the entire buffer might not be written to so it must be scrubbed )
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList(); //Can't clear until out copy from previous iteration is complete
		cl.queue.enqueueFillBuffer( clOut, 0, 0, sizeof( PVOCBuffer::MF ) * inChannelDataCount );

		//copy input buffer to device
		cl::copy( cl.queue, 
			getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * inChannelDataCount, 
			getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * inChannelDataCount, 
			clIn );

		//Calculate mapped frequencies
		auto mappedFrequenciesWriteHead = mappedFrequencies.begin();
		auto inBufferReadHead = getMFPointer( 0, 0, 0 ) + channel * inChannelDataCount;
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin, ++mappedFrequenciesWriteHead, ++inBufferReadHead )
				{
				XCDP_CANCEL_POINT( PVOC() );
				*mappedFrequenciesWriteHead = outFreqFunc( frame * frameToTime(), inBufferReadHead->f );
				}

		//Copy mapped frequencies
		cl::copy( cl.queue, mappedFrequencies.begin(), mappedFrequencies.end(), clMappedFrequencies );
		
		//Compute on device
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl_modifyFrequency( cl::EnqueueArgs( cl.queue, cl::NDRange( inChannelDataCount ) ), 
			clIn, 
			clMappedFrequencies,
			clF, 
			clOut, 
			numBins_i );
			
		//copy cl output buffer to out
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut,
			out.getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * inChannelDataCount, 
			out.getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * inChannelDataCount );
		}

	return out;
	}

PVOC PVOC::modifyTime_cl( Func2x1 outPosFunc, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	const uint32_t numBins = getNumBins();
	const uint32_t numFrames = getNumFrames();
	const uint32_t inChannelDataCount = numFrames * numBins;

	//Sample stretch function
	std::vector<float> outPosBuffer( inChannelDataCount );
	float lastOutputTime = 0;
	for( Frame frame = 0; frame < numFrames; ++frame )
		for( Bin bin = 0; bin < numBins; ++bin )
			{
			XCDP_CANCEL_POINT( PVOC() );
			const float outPos = outPosFunc( frame * frameToTime(), bin * binToFrequency() ) * timeToFrame();
			outPosBuffer[ frame * numBins + bin ] = outPos;
			if( outPos > lastOutputTime )
				lastOutputTime = outPos;
			}
	++lastOutputTime;

	auto format = getFormat();
	format.numFrames = ceil( lastOutputTime );
	PVOC out( format );
	const uint32_t outChannelDataCount = out.getNumFrames() * out.getNumBins();

	//prepare OpenCL
	auto cl = CLContext::get();
	static ProgramHelper programHelper( CLProgs::PVOC_modifyTime );
	cl::Buffer clIn( cl.context, CL_MEM_READ_ONLY, sizeof( PVOCBuffer::MF ) * inChannelDataCount );
	cl::Buffer clF( cl.context, CL_MEM_READ_ONLY, sizeof( float ) * inChannelDataCount );
	cl::Buffer clOut( cl.context, CL_MEM_WRITE_ONLY, sizeof( PVOCBuffer::MF ) * outChannelDataCount );
	cl::KernelFunctor< cl::Buffer, cl::Buffer, cl::Buffer, int, int > 
		cl_modifyTime( programHelper.program, "modifyTime" );

	//stretch function samples only need one copy, they're the same for all channels
	cl::copy( cl.queue, outPosBuffer.begin(), outPosBuffer.end(), clF );

	const int numBins_i = numBins;
	const int numPrevFrames_i = out.getNumFrames();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Clear cl out buffer ( the entire buffer might not be written to )
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList(); //Can't clear until out copy from previous iteration is complete
		cl.queue.enqueueFillBuffer( clOut, 0, 0, sizeof( PVOCBuffer::MF ) * outChannelDataCount );

		//copy input buffer to device
		cl::copy( cl.queue, 
			getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * inChannelDataCount, 
			getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * inChannelDataCount, 
			clIn );
		
		//Compute on device
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl_modifyTime( cl::EnqueueArgs( cl.queue, cl::NDRange( inChannelDataCount ) ), 
			clIn, 
			clF, 
			clOut, 
			numBins_i,
			numPrevFrames_i );
			
		//copy cl output buffer to out
		XCDP_CANCEL_POINT( PVOC() );
		cl.queue.enqueueBarrierWithWaitList();
		cl::copy( cl.queue, clOut,
			out.getMFPointer( 0, 0, 0 ) + ( channel + 0 ) * outChannelDataCount, 
			out.getMFPointer( 0, 0, 0 ) + ( channel + 1 ) * outChannelDataCount );
		}

	XCDP_CANCEL_POINT( PVOC() );
	cl.queue.enqueueBarrierWithWaitList();
	return out;
	}

#else

//Dummy functions for disabled opencl
PVOC PVOC::modify_cl( Func2x2 mod, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const 
	{ std::cout << "Warning: modify_cl called but OpenCL is disabled.\n"; return *this; }
PVOC PVOC::modifyFrequency_cl( Func2x1 outFreqFunc, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const 
	{ std::cout << "Warning: modifyFrequency_cl called but OpenCL is disabled.\n"; return *this; }
PVOC PVOC::modifyTime_cl( Func2x1 outPosFunc, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const 
	{ std::cout << "Warning: modifyTime_cl called but OpenCL is disabled.\n"; return *this; }

#endif

PVOC PVOC::repitch( Func2x1 factor, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_FUNCTION_LOG;
	return modifyFrequency_cl( [&factor]( vec2 tf ){ return factor(tf.x(),tf.y())*tf.y(); },
		interp, canceller );
	}

PVOC PVOC::stretch( Func2x1 factor, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_FUNCTION_LOG;
	return modifyTime_cl( [&factor]( vec2 tf ){ return factor(tf.x(),tf.y())*tf.x(); }, 
		interp, canceller );
	}

PVOC PVOC::stretch_spline( Func1x1 interpolation, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	const auto safeInterpolation = [&interpolation, this]( Frame frame )
		{
		return std::max( uint32_t( interpolation( frame * frameToTime() ) ), 1u );
		};
	
	//Set up output and spline x coordinates
	auto format = getFormat();
	format.numFrames = 0;
	std::vector<double> Xs( getNumFrames() );
	for( Frame frame = 0; frame < getNumFrames()-1; ++frame )
		{
		Xs[frame] = format.numFrames;
		format.numFrames += safeInterpolation( frame );
		}
	Xs[getNumFrames()-1] = format.numFrames;
	PVOC out( format );

	//Allocate Y coordinate vectors
	std::vector<double> magnitudeYs( Xs.size() );
	std::vector<double> frequencyYs( Xs.size() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Channel -> bin is correct here, despite non-sequential access to pvoc buffer
		for( Bin bin = 0; bin < getNumBins(); ++bin )
			{
			//For each x coordinate get corresponding frequency and magnitude
			for( Frame frame = 0; frame < getNumFrames(); ++frame )
				{
				XCDP_CANCEL_POINT( PVOC() );
				magnitudeYs[frame] = getMF( channel, frame, bin ).m;
				frequencyYs[frame] = getMF( channel, frame, bin ).f;
				}

			//Do the splines
			tk::spline magnitudeSpline, frequencySpline;
			magnitudeSpline.set_points(Xs,magnitudeYs);
			frequencySpline.set_points(Xs,frequencyYs);

			//Evaluate splines into out
			for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
				{
				XCDP_CANCEL_POINT( PVOC() );
				out.setMF( channel, frame, bin, 
					{ 
					(float) magnitudeSpline(frame),
					(float) frequencySpline(frame)
					} );
				}
			}
		}
			
	return out;
	}

PVOC PVOC::desample( Func2x1 factor, Interpolator interp, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	auto safeFactor = [&factor]( vec2 tf )
		{
		return std::max( 1.0f, factor( tf.x(),tf.y() ) );
		};
	auto downsample = stretch( [&safeFactor]( vec2 tf )
		{
		return 1.0f / safeFactor( tf ); 
		}, interp, canceller );
	return downsample.stretch( safeFactor, interp, canceller );
	}

PVOC PVOC::timeExtrapolate( float startTime, float endTime, float extrapolationTime, 
	Interpolator interpolator, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( PVOC() );

	// Input validation
	startTime = std::clamp( startTime, 0.0f, getLength() );
	endTime	  = std::clamp( endTime,   0.0f, getLength() );
	if( startTime >= endTime ) return PVOC();
	if( extrapolationTime <= 0 ) return PVOC();

	//Linear extrapolation via lines matching original PVOC at x_1 and x_2
	const uint32_t startFrame	= timeToFrame() * startTime;			//x_1
	const uint32_t endFrame		= timeToFrame() * endTime;				//x_2
	const uint32_t extFrames	= timeToFrame() * extrapolationTime;	//Extrapolated frames to generate

	auto format = getFormat();
	format.numFrames = endFrame + extFrames;
	PVOC out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Copy up to start frame
		for( Frame frame = 0; frame < startFrame; ++frame )
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				XCDP_CANCEL_POINT( PVOC() );
				out.setMF( channel, frame, bin, getMF( channel, frame, bin ) );
				}
			
		//Interpolate and continue beyond endFrame to extrapolate
		for( Frame frame = startFrame; frame < out.getNumFrames(); ++frame )
			{
			const float mix = interpolator (float( frame - startFrame ) / float( endFrame - startFrame ) ); 

			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				XCDP_CANCEL_POINT( PVOC() );
				const auto leftBin  =  getMF( channel, startFrame, bin );
				const auto rightBin =  getMF( channel, endFrame, bin );

				out.setMF( channel, frame, bin, 
					{ 
					std::abs( ( 1.0f - mix ) * leftBin.m + mix * rightBin.m ), 
					          ( 1.0f - mix ) * leftBin.f + mix * rightBin.f 
					} );
				}
			}
		}
	return out;
	}


}