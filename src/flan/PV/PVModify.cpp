#include "flan/PV/PV.h"

#include <iostream>
#include <algorithm>
#include <ranges>

#include "spline/spline.h"
#include "flan/Utility/iota_iter.h"
#include "flan/Utility/buffer_access.h"

using namespace std::ranges;

namespace flan {

PV PV::modify( const Function<TF, TF> & mod, const Interpolator & interp ) const
	{
	if( is_null() ) return PV();

	const uint32_t in_channel_data_count = get_num_frames() * get_num_bins();

	//Sample mod function and convert output to frame/bin
	auto mod_sampled = sample_function_over_domain( mod );
	mod_sampled.for_each( [&]( TF & v ){ 
		v.t = time_to_frame( v.t );
		v.f = frequency_to_bin( v.f ); 
		} );

	// Get the largest frame value among all mod samples, rounded up
	const float last_output_frame = std::ceil( mod_sampled.maximum( []( TF v ){ return v.t; } ) );

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

	std::vector<MF> in_modified( in_channel_data_count );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// Calculate and copy mapped frequencies and input magnitudes
		auto in_buffer_read_head = get_MF_pointer( 0, 0, 0 ) + channel * in_channel_data_count;
		auto in_modified_write_head = in_modified.begin();
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			for( Bin bin = 0; bin < get_num_bins(); ++bin, ++in_modified_write_head, ++in_buffer_read_head )
				{
				*in_modified_write_head = 
					{
					in_buffer_read_head->m,
					mod( TF{ frame_to_time( frame ), in_buffer_read_head->f } ).f
					};
				}
		
		std::for_each( FLAN_PAR_SEQ iota_iter( 1 ), iota_iter( get_num_frames() ), [&]( Frame frame )
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
					in_modified[ ( frame - 1 ) * get_num_bins() + bin - 1 ],
					in_modified[ ( frame - 0 ) * get_num_bins() + bin - 1 ],
					in_modified[ ( frame - 0 ) * get_num_bins() + bin - 0 ],
					in_modified[ ( frame - 1 ) * get_num_bins() + bin - 0 ] };
		
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

PV modify_frequency_base( 
	const PV & me, 
	const FunctionSample2d<Frequency> & mod, 
	const std::vector<Frequency> input_frequencies_modded, 
	const Interpolator & interp )
	{
	if( me.is_null() ) return PV();

	PV out( me.get_format() );
	out.clear_buffer();

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		const auto in_modified = input_frequencies_modded.begin() + channel * me.get_num_frames() * me.get_num_bins();

		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( me.get_num_frames() ), [&]( Frame frame )
			{
			// For each adjacent pair of bins
			for( Bin bin = 1; bin < me.get_num_bins(); ++bin )
				{
				// Get where those bins were mapped to by the mod function
				const int hiBinIndex = frame * me.get_num_bins() + bin;
				const float loBin = me.frequency_to_bin( mod[hiBinIndex - 1] );
				const float hiBin = me.frequency_to_bin( mod[hiBinIndex    ] );
				const bool forward = hiBin > loBin; // Check if the bins were mapped upside down

				const Bin loBinRound = forward ? std::ceil( loBin ) : std::floor( loBin );
				const Bin hiBinRound = forward ? std::ceil( hiBin ) : std::floor( hiBin );
				const Bin start_bin = std::clamp( loBinRound, 0, out.get_num_bins() - 1 );
				const Bin end_bin   = std::clamp( hiBinRound, 0, out.get_num_bins() - 1 );

				const MF & loMF = MF{ me.get_MF( channel, frame, bin - 1 ).m, in_modified[buffer_access( bin - 1, frame, me.get_num_bins() )] };
				const MF & hiMF = MF{ me.get_MF( channel, frame, bin     ).m, in_modified[buffer_access( bin, frame, me.get_num_bins())] };

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

PV PV::modify_frequency( const Function<TF, Frequency> & mod, const Interpolator & interp ) const
	{
	const auto mod_sampled = sample_function_over_domain( mod );

	std::vector<Frequency> in_modified; // Buffer for sampling mod on input frequencies
	in_modified.reserve( get_buffer().size() );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			for( Bin bin = 0; bin < get_num_bins(); ++bin )
				in_modified.push_back( mod( TF{ frame_to_time( frame ), get_MF( channel, frame, bin ).f } ) );

	return modify_frequency_base( *this, mod_sampled, in_modified, interp );
	}

PV PV::repitch( const Function<TF, float> & local_expansion_factor, const Interpolator & interp ) const
	{
	auto factor_sampled = sample_function_over_domain( local_expansion_factor );

	// Partial integral with respect to time. Factor_sampled becomes TF -> Frame.
	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		for( Bin bin = 1; bin < get_num_bins(); ++bin )
			factor_sampled.at( frame, bin ) += factor_sampled.at( frame, bin - 1 );

	// Bin to frequency conversion. Factor_sampled becomes TF->Frequency.
	for( int i = 0; i < factor_sampled.size(); ++i )
		factor_sampled[i] = bin_to_frequency( factor_sampled[i] );
 
	// Lerp the freq integrated function to map the frequency values stored in the MFs
	std::vector<Frequency> in_modified; // Buffer for sampling mod on input frequencies
	in_modified.reserve( get_buffer().size() );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			for( Bin bin = 0; bin < get_num_bins(); ++bin )
				{
				const fBin fbin = std::clamp( frequency_to_bin( get_MF( channel, frame, bin ).f ), 0.0f, fBin( get_num_bins() - 1 ) - 0.0001f );
				const Bin lo = std::floor( fbin );
				const Bin hi = lo + 1;
				const Frequency lo_freq = factor_sampled.at( frame, lo );
				const Frequency hi_freq = factor_sampled.at( frame, hi );
				const float r = fbin - lo;
				const Frequency lerp = lo_freq * ( 1.0f - r ) + hi_freq * r;
				
				in_modified.push_back( lerp );
				}

	return modify_frequency_base( *this, factor_sampled, in_modified, interp );
	}

PV modify_time_base( const PV & me, const FunctionSample2d<Second> & mod, const Interpolator & interp )
	{
	if( me.is_null() ) return PV();

	// Get the largest frame value among all mod samples, rounded up
	const float last_output_frame = std::ceil( me.time_to_frame( mod.maximum() ) );

	auto format = me.get_format();
	format.num_frames = last_output_frame;
	PV out( format );
	out.clear_buffer();

	std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( me.get_num_channels() ), [&]( Channel channel )
		{
		// Note the bin loop happens first here. There are a lot of considerations involved in the ordering, but
		//  it comes down to time displacement happening along frames. Paralellizing the frames requires synchronization
		//  that isn't worth the benifits, and the parallel loop being first is faster. The algorithm isn't cache friendly 
		//  either way because we store PV data in frame-major order.
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( me.get_num_bins() ), [&]( Bin bin )
			{
			// For each adjacent pair of frames
			std::for_each( iota_iter( 1 ), iota_iter( me.get_num_frames() ), [&]( Frame frame )
				{
				const float lFrame = me.time_to_frame( mod.at( frame - 1, bin ) );
				const float rFrame = me.time_to_frame( mod.at( frame, bin ) );
				const bool forward = rFrame > lFrame;

				const Frame start_frame = forward ? std::ceil( lFrame ) : std::floor( lFrame );
				const Frame end_frame   = forward ? std::ceil( rFrame ) : std::floor( rFrame );

				const MF & lMF = me.get_MF( channel, frame - 1, bin );
				const MF & rMF = me.get_MF( channel, frame,     bin );

				for( Frame x = start_frame; x != end_frame; forward? ++x : --x )
					{
					if( x < 0 || out.get_num_frames() <= x ) continue;	

					const float mix = interp( ( x - lFrame ) / ( rFrame - lFrame ) );
					const Magnitude w0 = ( 1.0f - mix ) * lMF.m;
					const Magnitude w1 = (        mix ) * rMF.m;
					const Magnitude totalWeight = w0 + w1;
					const float weightedFreqSum = w0 * lMF.f + w1 * rMF.f;

					if( totalWeight == 0.0f )
						return;

					auto & outMF = out.get_MF( channel, x, bin );
					outMF.f = ( outMF.f * outMF.m + weightedFreqSum ) / ( outMF.m + totalWeight );
					outMF.m += totalWeight;
					}
				} );
			} );
		} );

	return out;
	}

PV PV::modify_time( const Function<TF, Second> & mod, const Interpolator & interp ) const
	{
	// Sample mod function and convert output to frame/bin
	const auto mod_sampled = sample_function_over_domain( mod );
	return modify_time_base( *this, mod_sampled, interp );
	}

PV PV::stretch( const Function<TF, float> & local_expansion_factor, const Interpolator & interp ) const
	{
	auto factor_sampled = sample_function_over_domain( local_expansion_factor );

	// Partial integral with respect to time. Factor_sampled becomes TF -> Frame.
	for( Bin bin = 0; bin < get_num_bins(); ++bin )
		for( Frame frame = 1; frame < get_num_frames(); ++frame )
			factor_sampled.at( frame, bin ) += factor_sampled.at( frame - 1, bin );

	// Frame to second conversion. Factor_sampled becomes TF->Second.
	for( int i = 0; i < factor_sampled.size(); ++i )
		factor_sampled[i] = frame_to_time( factor_sampled[i] );

	return modify_time_base( *this, factor_sampled, interp );
	}

PV PV::stretch_spline( const Function<Second, float> & interpolation ) const
	{
	if( is_null() ) return PV();

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

PV PV::desample( const Function<TF, float> & decimation_ratio, const Interpolator & interp ) const
	{
	if( is_null() ) return PV();

	// Sample factor over frames and bins
	auto decimation_ratio_samples = sample_function_over_domain( decimation_ratio );

	PV out( get_format() );
	out.clear_buffer();

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// This stores which frames should be interpolated between
		std::vector<std::vector<Frame>> selectedFrames( get_num_bins() );

		// Factor is integrated, when the accumulator passes 1 the currest frame is selected as an interpolation endpoint
		std::vector<float> accums( get_num_bins(), 1 ); // Start at one to trigger interp from starting frame
		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )
			{
			for( Frame frame = 0; frame < get_num_frames(); ++frame )
				{
				float & accum = accums[bin];
				const float factor_c = std::clamp( decimation_ratio_samples.at( frame, bin ), 0.0f, 1.0f );
				accum += factor_c;
				if( accum >= 1.0f )
					{
					selectedFrames[bin].push_back( frame );
					accum -= 1.0f; // No need to worry about accum being bigger than 2, factor_c can't be bigger than 1
					}
				}
			} );

		std::for_each( FLAN_PAR_UNSEQ iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )
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

PV PV::smear_time( 
	const Function<TF, Second> & smear_size, 
	const Function<TF, int> & granularity,
	const Function<Second, float> & distribution
	) const
	{
	if( is_null() ) return PV();

	auto granularity_sampled = sample_function_over_domain( granularity );
	granularity_sampled.for_each( [&]( int & i ){ i = std::max( i, 1 ); } ); 

	auto smear_size_sampled = sample_function_over_domain( smear_size );
	smear_size_sampled.for_each( [&]( Second & s ){ s = std::max( s, 0.0f ); } );

	/* 
	We could use the input size for the output but I'd rather capture the fade-in
	But because the smear size could change over time we need to figure out just how big the output needs to beeeeee
	And oh my god we need to do it for every bin and use the max of those, fantastic

	Something else odd about this algorithm is that without an additional parameter, we can't know if far off in either time direction the smear_size gets huge 
	and contains the input again. If it does, the user may have intended to include the sound that would generate, but we can't just sample the smear_size
	out to infinity. The compromise is to sample the smear size on the input domain, and assume that the size remains constant outside the input domain.
	*/
	std::vector<Frame> leftmost_frame( get_num_bins(), 0 );
	std::vector<Frame> rightmost_frame( get_num_bins(), get_num_frames() - 1 );
	flan::for_each_i( get_num_bins(), ExecutionPolicy::Parallel_Unsequenced, [&]( Bin bin )
		{
		Frame & l = leftmost_frame[bin];
		Frame & r = rightmost_frame[bin];
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			const Frame expansion_frames = time_to_frame( smear_size_sampled.at( frame, bin ) );
			const Frame l_c = frame - expansion_frames;
			const Frame r_c = frame + expansion_frames;
			l = std::min( l, l_c );
			r = std::max( r, r_c );
			}
		} );
	const Frame true_leftmost_frame = *std::min_element( leftmost_frame.begin(), leftmost_frame.end() );
	const Frame true_rightmost_frame = *std::max_element( rightmost_frame.begin(), rightmost_frame.end() );

	// The distribution function is assumed to be relatively smooth, and so it will be sampled now rather than directly accessed later
	const Frame max_smear_frames = time_to_frame( smear_size_sampled.maximum() );
	const Frame dist_samples_2 = max_smear_frames * 2;
	auto distribution_sampled = distribution.sample( -dist_samples_2, dist_samples_2, 1.0f / dist_samples_2 );

	Format format = get_format();
	format.num_frames = true_rightmost_frame - true_leftmost_frame;
	PV out( format );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		flan::for_each_i( out.get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame out_frame )
			{
			const Frame in_frame = std::clamp( out_frame + true_leftmost_frame, 0, get_num_frames() - 1 );

			for( Bin bin = 0; bin < get_num_bins(); ++bin )
				{
				// Here we take a time average of the surrounding MF data
				const Frame expansion_frames = time_to_frame( smear_size_sampled.at( in_frame, bin ) );
				double mag_sum = 0;
				double freq_sum = 0;
				double total_dist_weight = 0;
				double dist_weight_used = 0;
				const int granularity_c = granularity_sampled.at( in_frame, bin );
			 	for( Frame frame_offset = -expansion_frames; frame_offset < expansion_frames; frame_offset += granularity_c )
			 		{
					const float smear_size_c = smear_size_sampled.at( in_frame, bin ); // Forced >= 0 earlier
					const float dist_input = frame_to_time( frame_offset ) / smear_size_c;
					int dist_sampled_access = distribution_sampled.size() * 0.5f * ( 1 + dist_input );
					dist_sampled_access = std::clamp( dist_sampled_access, 0, int( distribution_sampled.size() - 1 ) );
					const float dist_c = distribution_sampled[dist_sampled_access];
					total_dist_weight += dist_c;

					const Frame source_frame = out_frame + true_leftmost_frame + frame_offset;
					if( source_frame < 0 || source_frame >= get_num_frames() ) continue;

					const MF mf_c = get_MF( channel, source_frame, bin );
					dist_weight_used += dist_c;
					mag_sum += mf_c.m * dist_c;
					freq_sum += mf_c.f * dist_c;
			 		}
				// In this averaging we treat magnitude and frequency differently.
				// If we are drawing data from outside the edge of the input, we count the magnitude as 0
				// The frequency we don't want to count at all, as using any number would be wrong, so nothing is added.
				if( total_dist_weight > 0.0 ) mag_sum /= total_dist_weight;
				if( dist_weight_used > 0.0 ) freq_sum /= dist_weight_used;

				out.get_MF( channel, out_frame, bin ) = MF( mag_sum, freq_sum );
			 	}
			} );

	return out;
	}

PV PV::time_extrapolate( Second start_time, Second end_time, Second extrapolationTime, const Interpolator & interpolator ) const
	{
	if( is_null() ) return PV();

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
		std::copy( FLAN_PAR_UNSEQ channel_begin( channel ), channel_begin( channel ) + start_frame * get_num_bins(), out.channel_begin( channel ) );
			
		//Interpolate and continue beyond end_frame to extrapolate
		std::for_each( FLAN_PAR_UNSEQ iota_iter( start_frame ), iota_iter( out.get_num_frames() ), [&]( Frame frame ) 
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