/*
    flan note: This is a transcribed and sanatized version of the pYIN 
    algorithm available here https://code.soundsoftware.ac.uk/projects/pyin

    pYIN - A fundamental frequency estimator for monophonic audio
    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "flan/Audio/Audio.h"

#include <iostream>
#include <numeric>
#include <random>

#include "flan/defines.h"
#include "flan/FFTHelper.h"
#include "flan/DSPUtility.h"

namespace flan {

//Audio::FreqProbs Audio::PYIN( Channel channel, Frame start, Frame window_size, std::shared_ptr<FFTHelper> fft ) const
//	{
//    flan_PROCESS_START( Audio::FreqProbs() );
//
//    if( !fft ) fft = std::make_shared<FFTHelper>( autocorrelation_fft_size( window_size / 2 ), true, true, false );
//
//	// Get d'
//	const std::vector<float> d_prime = compute_d_prime( get_sample_pointer( channel, start ), window_size, fft );
//    const std::vector<std::pair<int,float>> probs = PYIN_Probability( d_prime );
//
//    FreqProbs yin;
//
//    for( auto tauProb : probs )
//		{
//        const Frequency currentF0 = float( get_sample_rate() ) / parabolic_interpolation( d_prime, tauProb.first );
//        yin.push_back( { currentF0, tauProb.second } );
//		}
//    
//    return yin;
//    }

//std::vector<Audio::FreqProbs> Audio::PYINs( Channel channel, Frame start, Frame end, Frame window_size, Frame jump ) const
//    {
//    flan_PROCESS_START( std::vector<Audio::FreqProbs>() );
//
//    if( end == -1 ) end = get_num_frames();
//
//    auto fft = std::make_shared<FFTHelper>( autocorrelation_fft_size( window_size / 2 ), true, true, false );
//
//    std::vector<FreqProbs> fps;
//
//    for( Frame frame = start; frame + window_size < end; frame += jump )
//        fps.push_back( PYIN( channel, frame, window_size, fft ) );
//        
//    return fps;
//    }

//static std::vector<std::pair<int,float>> PYIN_Probability( const std::vector<float> & d_prime ) 
//	{
//    static float betaCDF[100] = {
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-06
//        , 2e-06, 3e-06, 5e-06, 8e-06, 1.2e-05, 1.7e-05, 2.4e-05, 3.4e-05, 4.7e-05, 6.5e-05, 
//        8.8e-05, 0.000118, 0.000157, 0.000208, 0.000273, 0.000355, 0.000459, 0.000589, 0.000751, 
//        0.000952, 0.001199, 0.001502, 0.001872, 0.002321, 0.002863, 0.003515, 0.004296, 0.005227, 
//        0.006332, 0.007639, 0.009178, 0.010984, 0.013096, 0.015556, 0.018411, 0.021713, 0.025519, 
//        0.029891, 0.034896, 0.040608, 0.047105, 0.054471, 0.062795, 0.072172, 0.082702, 0.094487, 
//        0.107635, 0.122256, 0.138461, 0.156361, 0.176065, 0.197679, 0.221303, 0.247025, 0.274923, 
//        0.305056, 0.337462, 0.372152, 0.409102, 0.448249, 0.48948, 0.532624, 0.577441, 0.623612, 
//        0.670722, 0.71825, 0.765548, 0.811825, 0.856126, 0.89731, 0.934022, 0.964668, 0.987383, 0.999997 };
//
//    std::gamma_distribution<float> gamma( 1.0f, 1.0f );
//
//    // For each local minima, assign peak probablity weight based on precomputed cumulative beta distribution
//    // Also, find sum and max of probabilities
//    std::vector<std::pair<int,float>> peakProb;
//    float sumProb = 0;
//    float maxProb = 0;
//    for( int tau = 2; tau + 1 < d_prime.size(); ++tau )
//        {
//        // If we are at a min
//        if( d_prime[tau+1] > d_prime[tau] && d_prime[tau-1] > d_prime[tau] ) 
//            {
//            // Get betaPDF from precomputed 100 entry table
//            const int betaIndex = 99 - std::floor( d_prime[tau] * 100 );
//            if( betaIndex < 16 ) continue; // Everything less has betaCDF too close to 0
//            peakProb.push_back( { tau, betaCDF[ betaIndex ] } );
//            sumProb += betaCDF[ betaIndex ];
//            if( betaCDF[ betaIndex ] > maxProb ) maxProb = betaCDF[ betaIndex ];  
//            }
//        }
//      
//    // Normalize so that sum of probabilities = maxProb
//    if( sumProb > 0 ) 
//        {
//        const float probScale = maxProb / sumProb;
//        std::for_each( peakProb.begin(), peakProb.end(), [probScale]( std::pair<int,float> & p ){ p.second *= probScale; } );
//        }
//        
//    return peakProb;
//	}

}