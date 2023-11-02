#include "flan/Audio/Audio.h"

#include <iostream>
#include <execution>

#include "WDL/resample.h"
#include "r8brain/CDSPResampler.h"

#include "flan/WindowFunctions.h"
#include "flan/Utility/iota_iter.h"

using namespace flan;

Audio Audio::resample( FrameRate new_sample_rate ) const
	{
	if( is_null() ) return Audio(); 

	if( new_sample_rate == get_sample_rate() ) 
		return copy();

	auto format = get_format();
	format.num_frames *= new_sample_rate / get_sample_rate();
	format.sample_rate = new_sample_rate;
	Audio out( format );

	r8b::CDSPResampler resampler( get_sample_rate(), new_sample_rate, get_num_frames() );
	resampler.oneshot<float, float>( &get_buffer()[0], get_buffer().size(), &out.get_buffer()[0], out.get_buffer().size() );

	return out;
	}

Audio Audio::convert_to_mid_side() const
	{
	if( is_null() ) return Audio();

	if( get_num_channels() != 2 )
		{
		std::cout << "Can't transform non-stereo Audio between Mid-Side and Left-Right formats." << std::endl;
		return copy();
		}
	Audio out( get_format() );
	flan::for_each_i( get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame frame )
		{
		out.set_sample( 0, frame, get_sample( 0, frame ) + get_sample( 1, frame ) );
		out.set_sample( 1, frame, get_sample( 0, frame ) - get_sample( 1, frame ) );
		} );
	return out;
	}

Audio Audio::convert_to_left_right() const
	{
		return convert_to_mid_side();
	}

Audio Audio::convert_to_stereo() const
	{
	if( is_null() ) return Audio();

	auto format = get_format();
	format.num_channels = 2;
	Audio out( format );

	switch( get_num_channels() )
		{
		case 1:
			{
			// Copy input channel one into the output channel one while scaling down by sqrt2
			std::transform( FLAN_PAR_UNSEQ get_buffer().begin(), get_buffer().end(), out.get_buffer().begin(), 
				[sqrt2 = std::sqrt( 2.0f )]( Sample s ){ return s / sqrt2; } );

			// Copy output channel one to output channel two
			std::copy( FLAN_PAR_UNSEQ out.get_buffer().begin(), out.get_buffer().begin() + get_num_frames(), out.get_buffer().begin() + get_num_frames() );
			break;
			}
		case 2:
			return copy();
		default:
			std::cout << "I don't know how to convert that number of channels to stereo." << std::endl;
			
		}
	return out;
	}

Audio Audio::convert_to_mono() const
	{
	if( is_null() ) return Audio();

	auto format = get_format();
	format.num_channels = 1;
	Audio out( format );

	flan::for_each_i( get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame frame )
		{
		Sample sampleAccumulator = 0;
		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			sampleAccumulator += get_sample( channel, frame );
		out.set_sample( 0, frame, sampleAccumulator / get_num_channels() );
		} );

	return out;
	}

Function<Second, Amplitude> Audio::convert_to_function(
	) const
	{
	if( is_null() ) return 0;

	Audio out = convert_to_mono();

	// Needing this is a c++ issue, not stealing the vector makes the lamba capture initialization get flagged as trying to call Audio copy ctor
	std::vector<Sample> buffer = std::move( out.get_buffer() );

	return [buffer = std::move( buffer ), sample_rate = get_sample_rate()]( Second t )
		{
		const Frame frame = t * sample_rate;
		if( 0 <= frame && frame < buffer.size() )
			return buffer[frame];
		return 0.0f;
		};
	}

// class QVOCBuffer {
// 	public:
// 		/*
// 		There is no natural ordering of QVOC data. The data is stored here in a single buffer so that ascending channel -> frame -> bin 
// 		iteration gives sequential memory access.
// 		*/

// 		// Takes a single channel of data to collect size info, doesn't actually get the data
// 		QVOCBuffer( const std::vector<std::vector<std::complex<float>>> & cqOut )
// 			{
// 			// Fill index buffer with the indices of the starts of each frame of data
// 			index_buffer.reserve( cqOut.size() );
// 			index_buffer.push_back( 0 );
// 			std::transform( cqOut.begin(), cqOut.end(), std::back_inserter( index_buffer ), 
// 				[]( const std::vector<std::complex<float>> & v ){ return v.size(); } );

// 			std::partial_sum( index_buffer.begin(), index_buffer.end(), index_buffer.begin() );

// 			buffer.resize( getChannelSize() );
// 			}

// 		MF & get_MF( Channel channel, Frame frame, Bin bin )
// 			{
// 			return buffer[ getChannelSize() * channel + index_buffer[frame] + bin ];
// 			}

// 		size_t getChannelSize() const { return index_buffer.back(); }
// 		Frame get_num_frames() const { return index_buffer.size() - 1; }
// 		Bin get_num_bins( Frame f ) const { return index_buffer[f+1] - index_buffer[f]; }

// 	private:
// 		std::vector<size_t> index_buffer;
// 		std::vector<MF> buffer;
// };

// #include "constantQ/src/ConstantQ.h"
// #include "constantQ/src/CQInverse.h"
// Audio Audio::convertToQVOC() const
// 	{
// 	CQParameters params( get_sample_rate(), 20, 20000, 128 );

// 	// Forward constant-q transform
// 	ConstantQ cq( params );
// 	auto cqOut = cq.process( get_buffer() );
// 	auto cqOutR = cq.getRemainingOutput();
// 	cqOut.insert( cqOut.end(), cqOutR.begin(), cqOutR.end() );
// 	auto cqOutProperties = cq.getKernalProperties();
// 	const float dft_size = cqOutProperties.fft_size; 

// 	// Transform from cq output format to QVOCBuffer (MF format)
// 	QVOCBuffer QVOC( cqOut );
// 	std::vector<float> phase_buffer( cq.getTotalBins(), 0 );
// 	for( Frame frame = 0; frame < cqOut.size(); ++frame )
// 		{
// 		for( Bin bin = 0; bin < cqOut[frame].size(); ++bin )
// 			{
// 			const Bin flippedBin = cqOut[frame].size() - 1 - bin;
// 			const float hopSize = cq.getColumnHop( bin );

// 			const float magnitude = std::abs( cqOut[frame][flippedBin] );

// 			// We can't use frequency estimation with too low of a temporal resolution
// 			if( hopSize >= dft_size )
// 				QVOC.get_MF( 0, frame, bin ) = { magnitude, cq.bin_to_frequency( bin ) };

// 			const float phase = std::arg( cqOut[frame][flippedBin] );
// 			const float phaseDiff = phase - phase_buffer[bin]; // Actual change in radians for the current hop
// 			phase_buffer[bin] = phase; // Update phase buffer
// 			const float expectedPhaseDiff = bin * pi2 * hopSize / dft_size; // Expected change in radians for the current hop.
// 			const float deltaPhase = phaseDiff - expectedPhaseDiff; // The difference between the actual and the expected phase changes, aka how far off was our sinusoid from the expected.
// 			const float wrappedDeltaPhase = deltaPhase - pi2 * std::round( deltaPhase / pi2 ); // Wrap delta phase into [-pi,pi]
// 			const float binDeviation = wrappedDeltaPhase / pi2 * dft_size / hopSize;
// 			const Frequency frequency = cq.bin_to_frequency( bin + binDeviation );

// 			QVOC.get_MF( 0, frame, bin ) = { magnitude, frequency };
// 			}
// 		}

// 	// Modify qvoc data
// 	for( Frame frame = 0; frame < QVOC.get_num_frames(); ++frame )
// 		{
// 		for( Bin bin = 0; bin < QVOC.get_num_bins( frame ) - 1; ++bin )
// 			{
// 			// QVOC.get_MF( 0, frame, bin ) = QVOC.get_MF( 0, frame, bin + 1 );
// 			// QVOC.get_MF( 0, frame, bin ).f /= std::pow( 2.0f, 1.0f / 12 );
// 			}
// 		}

// 	// Transform QVOCBuffer back to cq format
// 	std::for_each( phase_buffer.begin(), phase_buffer.end(), []( float & x ){ x = 0; } );

// 	// We need to construct a vector of vectors with the correct sizes
// 	std::vector<std::vector<std::complex<float>>> cqIn( QVOC.get_num_frames(), std::vector<std::complex<float>>{} );
// 	for( Frame frame = 0; frame < cqIn.size(); ++frame )
// 		cqIn[frame].resize( QVOC.get_num_bins( frame ), 0 );

// 	// Apply the MF->complex transform
// 	for( Frame frame = 0; frame < cqIn.size(); ++frame )
// 		{
// 		for( Bin bin = 0; bin < cqIn[frame].size(); ++bin )
// 			{
// 			const Bin flippedBin = cqOut[frame].size() - 1 - bin;
// 			const float hopSize = cq.getColumnHop( bin );

// 			const MF mf = QVOC.get_MF( 0, frame, bin );

// 			phase_buffer[bin] += cq.frequency_to_bin( mf.f ) * pi2 * hopSize / dft_size;
			
// 			cqIn[frame][flippedBin] = std::polar( mf.m, phase_buffer[bin] );
// 			}
// 		}

// 	// Inverse constant-q transform
// 	CQInverse cqi( params );
// 	auto cqiOut = cqi.process( cqIn );
// 	auto cqiOutR = cqi.getRemainingOutput();
// 	cqiOut.insert( cqiOut.end(), cqiOutR.begin(), cqiOutR.end() );

// 	// Create output Audio
// 	auto format = get_format();
// 	format.num_frames = cqiOut.size();
// 	Audio out( format );

// 	// Copy to output
// 	for( Frame frame = 0; frame < cqiOut.size(); ++frame )
// 		out.get_sample( 0, frame ) = cqiOut[frame];

// 	return out;
// 	}

