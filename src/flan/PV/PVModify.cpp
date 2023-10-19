#include "flan/PV/PV.h"

#include <iostream>
#include <algorithm>
#include <ranges>

#include "spline/spline.h"
#include "flan/Utility/iota_iter.h"

using namespace std::ranges;

namespace flan {

PV PV::modify( const Function<TF, TF> & mod, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	const uint32_t in_channel_data_count = get_num_frames() * get_num_bins();

	//Sample mod function and convert output to frame/bin
	auto mod_sampled = sample_function_over_domain( mod );
	std::for_each( std::execution::par_unseq, mod_sampled.begin(), mod_sampled.end(), [&]( TF & v ){ 
		v.t = time_to_frame( v.t );
		v.f = frequency_to_bin( v.f ); 
		} );

	// Get the largest frame value among all mod samples, rounded up
	const float last_output_frame = std::ceil( std::ranges::max_element( mod_sampled, std::ranges::less(), []( TF v ){ return v.t; } )->t );

	if( frame_to_time( last_output_frame ) > 60.0f * 10.0f ) // Outfile longer than 10 minutes?
		{
		std::cout << "PV::modify tried to make a file longer than 10 minutes, which is currently disabled";
		return PV();
		}

	auto format = get_format();
	format.num_frames = ceil( last_output_frame );
	PV out( format );
	out.clear_buffer();

	// Using a mutex per frame has proven to be faster in practice than other lock sets
	// A C array is used because mutex cannot be copied or moved, ruling out std::vector, and the size isn't known at compile-time, ruling out std::array
	// C++ should probably cover this use case but no std container does
	auto mutexBuffer = std::unique_ptr<std::mutex[]>( new std::mutex[format.num_frames] );

	std::vector<MF> inModified( in_channel_data_count );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// Calculate and copy mapped frequencies and input magnitudes
		auto inBufferReadHead = get_MFPointer( 0, 0, 0 ) + channel * in_channel_data_count;
		auto modDataSamplesWriteHead = inModified.begin();
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			for( Bin bin = 0; bin < get_num_bins(); ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				{
				*modDataSamplesWriteHead = 
					{
					inBufferReadHead->m,
					mod( TF{ frame_to_time( frame ), inBufferReadHead->f } ).f
					};
				}
		
		std::for_each( std::execution::par, iota_iter( 1 ), iota_iter( get_num_frames() ), [&]( Frame frame )
			{
			//std::for_each( std::execution::par, iota_iter(1), iota_iter(get_num_bins()), [&]( Bin bin )
			for( Bin bin = 1; bin < get_num_bins(); ++bin )
				{
				// For each square in the input, the mod maps that square to this quad
				auto to_vec2 = []( TF tf ){ return vec2{ tf.t, tf.f }; };
				const std::array<vec2, 4> p = {
					to_vec2( mod_sampled[ ( frame - 1 ) * get_num_bins() + bin - 1 ] ),
					to_vec2( mod_sampled[ ( frame - 0 ) * get_num_bins() + bin - 1 ] ),
					to_vec2( mod_sampled[ ( frame - 0 ) * get_num_bins() + bin - 0 ] ),
					to_vec2( mod_sampled[ ( frame - 1 ) * get_num_bins() + bin - 0 ] ) };
		
				const std::array<MF, 4> pMF = {
					inModified[ ( frame - 1 ) * get_num_bins() + bin - 1 ],
					inModified[ ( frame - 0 ) * get_num_bins() + bin - 1 ],
					inModified[ ( frame - 0 ) * get_num_bins() + bin - 0 ],
					inModified[ ( frame - 1 ) * get_num_bins() + bin - 0 ] };
		
				const vec2 D12 = p[1] - p[0];
				const vec2 D23 = p[2] - p[1];
				const vec2 D34 = p[3] - p[2];
				const vec2 D41 = p[0] - p[3];
			
				// Find bounding box containing quad to iterate over
				const Frame minx = fmax( floor( fmin( fmin( p[0].x(), p[1].x() ), fmin( p[2].x(), p[3].x() ) ) ), 0 );
				const Bin 	miny = fmax( floor( fmin( fmin( p[0].y(), p[1].y() ), fmin( p[2].y(), p[3].y() ) ) ), 0 );
				const Frame maxx = fmin( ceil(  fmax( fmax( p[0].x(), p[1].x() ), fmax( p[2].x(), p[3].x() ) ) ), out.get_num_frames() - 1 );
				const Bin 	maxy = fmin( ceil(  fmax( fmax( p[0].y(), p[1].y() ), fmax( p[2].y(), p[3].y() ) ) ), out.get_num_bins() - 1 );
			
				// Iterate over bounding box
				for( Frame x = minx; x <= maxx; ++x )
					{
					for( Bin y = miny; y <= maxy; ++y )
						{	
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
									break;
								m = -quadC / quadB;
								}
							else
								{
								const float descriminant = quadB * quadB - 4.0f * quadA * quadC;
								if( descriminant < 0 ) break;
								m = ( -quadB + sqrt( descriminant ) ) / ( 2.0f * quadA ); // Quadratic formula
								}
							const float lDenominator = alpha[1] + alpha[3] * m;
							if( lDenominator == 0 ) break;
							const float l = ( x - alpha[0] - alpha[2] * m ) / lDenominator;
					
							// Make sure (l,m) is within a unit square + epsilon
							const float epsilon = 0.0001f;
							if( fabs( l - 0.5f ) > 0.5f + epsilon || fabs( m - 0.5f ) > 0.5f + epsilon ) break;

							const float interpL = interp( l );
							const float interpM = interp( m );

							const std::array<Magnitude, 4> w = {
								( 1.0f - interpL ) * ( 1.0f - interpM ) * pMF[0].m,
								(        interpL ) * ( 1.0f - interpM ) * pMF[1].m,
								(        interpL ) * (        interpM ) * pMF[2].m,
								( 1.0f - interpL ) * (        interpM ) * pMF[3].m };
							const float totalWeight = w[0] + w[1] + w[2] + w[3];
							if( totalWeight <= 0.0f ) break;

							// Choosing the MF assigned in this algorithm has many less than perfect options.
							// WFS: Weighted frequency sum gives inharmonic frequencies
							// MM: Max-mag requires tracking additional info
							// FIMM: In freq-innacurate max-mag the current quad's max-mag is compared to the write location. If it is louder,
							// overwrite frequency and then add to the magnitude. The problem with this strategy is that if a bin has quiet
							// data written to it more than once followed by loud data, a quiet frequency might be used. 
							// MIMM: Mag-innacurate max-mag is the same idea but with overwriting the magnitude instead of adding to it. 
							// Similar to freq-innacurate, this causes magnitude inaccuracies with repeated writes.
							// MIMM is being used currently.

							const auto maxWeightIter = std::max_element( w.begin(), w.end() );
							const int maxWeightedQuadMagIndex = std::distance( w.begin(), maxWeightIter );
							const Frequency loudestQuadFrequency = pMF[maxWeightedQuadMagIndex].f;
							
							// Lock mutex for frame x
							std::lock_guard<std::mutex> lock( mutexBuffer.get()[x] );
							MF & outMF = out.get_MF( channel, x, y );
							if( *maxWeightIter > outMF.m )
								outMF = { *maxWeightIter, loudestQuadFrequency };
								

							// Weighted Frequency Sum: kind of a muddy sound
							// const float weightedFreqSum = 
							// 			  pMF[0].f * w[0]
							// 			+ pMF[1].f * w[1]
							// 			+ pMF[2].f * w[2]
							// 			+ pMF[3].f * w[3];
							// // Lock mutex for frame x
							// std::lock_guard<std::mutex> lock( mutexBuffer.get()[x] );
							// MF & outMF = out.get_MF( channel, x, y );
							// outMF.f = ( outMF.f * outMF.m + weightedFreqSum ) / ( outMF.m + totalWeight );
							// outMF.m += totalWeight;
							}
						}
					}
				}
			} );
		}

	return out;
	}

PV PV::modify_frequency( const Function<TF, Frequency> & mod, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	PV out( get_format() );
	out.clear_buffer();

	const size_t in_channel_data_count = get_num_frames() * get_num_bins();

	//Sample mod function and convert output to frame/bin
	auto mod_sampled = sample_function_over_domain( mod );
	std::for_each( std::execution::par_unseq, mod_sampled.begin(), mod_sampled.end(), [&]( Frequency & f ){ 
		f = frequency_to_bin( f ); 
		} );

	std::vector<MF> inModified( in_channel_data_count ); // Buffer for sampling mod on input frequencies

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// Sample mod over input frequencies into inModified
		// This is computed early because the mod function isn't guaranteed to be thread safe
		auto inBufferReadHead = get_MFPointer( 0, 0, 0 ) + channel * in_channel_data_count;
		auto modDataSamplesWriteHead = inModified.begin();
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			for( Bin bin = 0; bin < get_num_bins(); ++bin, ++modDataSamplesWriteHead, ++inBufferReadHead )
				*modDataSamplesWriteHead = {
					inBufferReadHead->m,
					mod( TF( frame_to_time( frame ), inBufferReadHead->f ) ) 
					};

		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
			{
			// For each adjacent pair of bins
			for( Bin bin = 1; bin < get_num_bins(); ++bin )
				{
				// Get where those bins were mapped to by the mod function
				const int hiBinIndex = frame * get_num_bins() + bin;
				const float loBin = mod_sampled[hiBinIndex - 1];
				const float hiBin = mod_sampled[hiBinIndex    ];
				const bool forward = hiBin > loBin; // Check if the bins were mapped upside down

				const Bin loBinRound = forward ? std::ceil( loBin ) : std::floor( loBin );
				const Bin hiBinRound = forward ? std::ceil( hiBin ) : std::floor( hiBin );
				const Bin start_bin = std::clamp( loBinRound, 0, out.get_num_bins() - 1 );
				const Bin end_bin   = std::clamp( hiBinRound, 0, out.get_num_bins() - 1 );

				const MF & loMF = inModified[ buffer_access( bin - 1, frame, get_num_bins() ) ];
				const MF & hiMF = inModified[ buffer_access( bin    , frame, get_num_bins() ) ];

				for( Bin y = start_bin; y != end_bin; forward? ++y : --y )
					{
					const float mix = interp( ( float( y ) - loBin ) / ( hiBin - loBin ) );

					const Magnitude w0 = ( 1.0f - mix ) * loMF.m;
					const Magnitude w1 = (        mix ) * hiMF.m;
					
					const MF & max_magMF = w0 < w1 ? loMF : hiMF;
					MF & outMF = out.get_MF( channel, frame, y );
					if( max_magMF.m > outMF.m )
						{
						outMF.m += max_magMF.m;
						outMF.f = max_magMF.f;
						}

					// WFS:
					// const Magnitude totalWeight = w0 + w1;
					// const float weightedFreqSum = w0 * loMF.f + w1 * hiMF.f;
					// auto & outMF = out.get_MF( channel, frame, y );
					// outMF.f = ( outMF.m * outMF.f + weightedFreqSum ) / ( outMF.m + totalWeight );
					// outMF.m += totalWeight;
					}
				}
			} );
		}

	return out;
	}

PV PV::modify_time( const Function<TF, Second> & mod, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	// Sample mod function and convert output to frame/bin
	auto mod_sampled = sample_function_over_domain( mod );
	std::for_each( std::execution::par_unseq, mod_sampled.begin(), mod_sampled.end(), [&]( Second & t ){ 
		t = time_to_frame( t );
		} );

	// Get the largest frame value among all mod samples, rounded up
	const float last_output_frame = std::ceil( *std::ranges::max_element( mod_sampled ) );

	auto format = get_format();
	format.num_frames = last_output_frame;
	PV out( format );
	out.clear_buffer();

	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_channels() ), [&]( Channel channel )
		{
		// Note the bin loop happens first here. There are a lot of considerations involved in the ordering, but
		//  it comes down to time displacement happening along frames. Paralellizing the frames requires synchronization
		//  that isn't worth the benifits, and the parallel loop being first is faster. The algorithm isn't cache friendly 
		//  either way because we store PV data in frame-major order.
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )
			{
			// For each adjacent pair of frames
			std::for_each( std::execution::seq, iota_iter( 1 ), iota_iter( get_num_frames() ), [&]( Frame frame )
				{
				const float lFrame = mod_sampled[ buffer_access( bin, frame - 1, get_num_bins() ) ];
				const float rFrame = mod_sampled[ buffer_access( bin, frame    , get_num_bins() ) ];
				const bool forward = rFrame > lFrame;

				const Frame start_frame = forward ? std::ceil( lFrame ) : std::floor( lFrame );
				const Frame end_frame   = forward ? std::ceil( rFrame ) : std::floor( rFrame );

				const MF & lMF = get_MF( channel, frame - 1, bin );
				const MF & rMF = get_MF( channel, frame,     bin );

				for( Frame x = start_frame; x != end_frame; forward? ++x : --x )
					{
					if( x < 0 || out.get_num_frames() <= x ) continue;	

					const float mix = interp( ( x - lFrame ) / ( rFrame - lFrame ) );
					const Magnitude w0 = ( 1.0f - mix ) * lMF.m;
					const Magnitude w1 = (        mix ) * rMF.m;
					const Magnitude totalWeight = w0 + w1;
					const float weightedFreqSum = w0 * lMF.f + w1 * rMF.f;

					auto & outMF = out.get_MF( channel, x, bin );
					outMF.f = ( outMF.f * outMF.m + weightedFreqSum ) / ( outMF.m + totalWeight );
					outMF.m += totalWeight;
					}
				} );
			} );
		} );

	return out;
	}

PV PV::repitch( const Function<TF, float> & factor, Interpolator interp ) const
	{
	flan_FUNCTION_LOG;
	return modify_frequency( [&]( TF tf ){ return factor( tf ) * tf.f; }, interp );
	}

PV PV::stretch( const Function<TF, float> & factor, Interpolator interp ) const
	{
	flan_FUNCTION_LOG;
	return modify_time( [&factor]( TF tf ){ return factor( tf ) * tf.t; }, interp );
	}

PV PV::stretch_spline( const Function<Second, float> & interpolation ) const
	{
	flan_PROCESS_START( PV() );

	const auto safeInterpolation = [&interpolation, this]( Frame frame )
		{
		return std::max( uint32_t( interpolation( frame * frame_to_time( 1 ) ) ), 1u );
		};
	
	//Set up output and spline x coordinates
	auto format = get_format();
	format.num_frames = 0;

	std::vector<double> Xs( get_num_frames() );
	for( Frame frame = 0; frame < get_num_frames()-1; ++frame )
		{
		Xs[frame] = format.num_frames;
		format.num_frames += safeInterpolation( frame );
		}
	Xs[get_num_frames()-1] = format.num_frames;
	PV out( format );

	//Allocate Y coordinate vectors
	std::vector<double> magnitudeYs( Xs.size() );
	std::vector<double> frequencyYs( Xs.size() );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		//Channel -> bin is correct here, despite non-sequential access to flan buffer
		for( Bin bin = 0; bin < get_num_bins(); ++bin )
			{
			//For each x coordinate get corresponding frequency and magnitude
			for( Frame frame : iota_view( 0, get_num_frames() ) )
				{
				magnitudeYs[frame] = get_MF( channel, frame, bin ).m;
				frequencyYs[frame] = get_MF( channel, frame, bin ).f;
				}

			//Do the splines
			tk::spline magnitudeSpline, frequencySpline;
			magnitudeSpline.set_points(Xs,magnitudeYs);
			frequencySpline.set_points(Xs,frequencyYs);

			//Evaluate splines into out
			for( Frame frame = 0; frame < out.get_num_frames(); ++frame )
				{
				out.set_MF( channel, frame, bin, 
					{ 
					(float) magnitudeSpline(frame),
					(float) frequencySpline(frame)
					} );
				}
			}
		}
			
	return out;
	}

PV PV::desample( const Function<TF, float> & events_per_second, Interpolator interp ) const
	{
	// Sample factor over frames and bins
	auto events_per_second_samples = sample_function_over_domain( events_per_second );

	flan_PROCESS_START( PV() );

	PV out( get_format() );
	out.clear_buffer();

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// This stores which frames should be interpolated between
		std::vector<std::vector<Frame>> selectedFrames( get_num_bins() );

		// Factor is integrated, when the accumulator passes 1 the currest frame is selected as an interpolation endpoint
		std::vector<float> accums( get_num_bins(), 1 ); // Start at one to trigger interp from starting frame
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )
			{
			for( Frame frame = 0; frame < get_num_frames(); ++frame )
				{
				float & accum = accums[bin];
				if( events_per_second_samples[buffer_access( bin, frame, get_num_bins() )] < 1.0f ) continue;
				const float factor_c = events_per_second_samples[buffer_access( bin, frame, get_num_bins() )];
				accum += factor_c * frame_to_time( 1 );
				if( accum >= 1.0f )
					{
					selectedFrames[bin].push_back( frame );
					accum -= 1.0f; // No need to worry about accum being bigger than 2, factor_c can't be bigger than 1
					}
				}
			} );

		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )
			{
			// For each pair of selected interpolation endpoints in selectedFrames, intepolate the frames between them
			std::vector<Frame> & selectedFrames_c = selectedFrames[bin];
			if( selectedFrames_c.size() < 2 ) return;
			for( const int selectedFrameIndex : std::views::iota( 0ull, selectedFrames_c.size() - 1 ) )
				{
				const Frame lFrame = selectedFrames_c[selectedFrameIndex    ];
				const Frame rFrame = selectedFrames_c[selectedFrameIndex + 1];
				const MF lMF = get_MF( channel, lFrame, bin );
				const MF rMF = get_MF( channel, rFrame, bin );
				for( Frame frame = lFrame; frame < rFrame; ++frame )
					{
					const float mix = interp( float( frame - lFrame ) / ( rFrame - lFrame ) ); 
					const float w0 = ( 1.0f - mix ) * lMF.m;
					const float w1 = (        mix ) * rMF.m;
					const MF mf = 
						{
						w0 + w1,
						w0 > w1 ? lMF.f : rMF.f
						};
					// Weighting by frequency is another option, it doesn't noticably change the output for inputs I've tested though
					// const MF mf = 
					// 	{
					// 	w0 + w1,
					// 	( w0 * lMF.f + w1 * rMF.f ) / ( w0 + w1 )
					// 	};
					out.set_MF( channel, frame, bin, mf );
					}
				}
			} );
		}

	return out;
	}

PV PV::time_extrapolate( Second start_time, Second end_time, Second extrapolationTime, Interpolator interpolator ) const
	{
	flan_PROCESS_START( PV() );

	// Input validation
	start_time = std::clamp( start_time, 0.0f, get_length() );
	if( end_time == -1 ) end_time = get_length();
	end_time	  = std::clamp( end_time,   0.0f, get_length() );
	if( start_time >= end_time ) return PV();
	if( extrapolationTime <= 0 ) return PV();

	const Frame start_frame	= time_to_frame( start_time );			//x_1
	const Frame end_frame	= time_to_frame( end_time );			//x_2
	const Frame extFrames	= time_to_frame( extrapolationTime );	//Extrapolated frames to generate

	auto format = get_format();
	format.num_frames = end_frame + extFrames;
	PV out( format );
	out.clear_buffer();

	// Sample interpolator
	std::vector<float> interpSamples( out.get_num_frames() - start_frame );
	for( Frame frame = 0; frame < interpSamples.size(); ++frame )
		interpSamples[frame] = interpolator( float( frame - start_frame ) / float( end_frame - start_frame ) );

	for( Channel channel = 0; channel < get_num_channels(); ++channel ) 
		{
		// Copy up to the start frame from input into output
		std::copy( std::execution::par_unseq, channel_begin( channel ), channel_begin( channel ) + start_frame * get_num_bins(), out.channel_begin( channel ) );
			
		//Interpolate and continue beyond end_frame to extrapolate
		std::for_each( std::execution::par_unseq, iota_iter( start_frame ), iota_iter( out.get_num_frames() ), [&]( Frame frame ) 
			{
			const float mix = interpSamples[ frame - start_frame ];

			for( Bin bin = 0; bin < get_num_bins(); ++bin ) 
				{
				const auto leftMF  =  get_MF( channel, start_frame, bin );
				const auto rightMF =  get_MF( channel, end_frame, bin );

				// This algorithm decides bin based on frequency. Phase vocoder outputs can have multiple bins in a region
				// containing the same frequency, and squishing all those into one bin will cause the output to sound less full.
				// To fix that, we note for each bin at end_time how far off the frequency in that bin is from where it "should" be. That shift
				// is later added back in when deciding what bin should be used for the current frequency, maintaining fidelity.
				const Bin rightBinShift = bin - frequency_to_bin( rightMF.f );

				const MF extrapMF = { std::abs( ( 1.0f - mix ) * leftMF.m + mix * rightMF.m ),
					          					( 1.0f - mix ) * leftMF.f + mix * rightMF.f };
				const Bin extrapBin = frequency_to_bin( extrapMF.f ) + rightBinShift;
				if( extrapBin < 0 || extrapBin >= get_num_bins() ) continue;

				MF & outMF = out.get_MF( channel, frame, extrapBin );

				if( extrapMF.m > outMF.m ) 
					outMF = extrapMF;
				}
			} );
		}
	return out;
	}

}