// #pragma once

// #include "SQPVBuffer.h"

// #include "flan/Function.h"

// namespace flan {

// class Audio;
// class Graph;

// class SQPV : public SQPVBuffer
// {
// public:
// 	SQPV();
// 	SQPV( const Format & );
// 	SQPV( SQPVBuffer && );

// 	Audio convert_to_audio() const;
// 	Audio convert_to_lr_audio() const;
// 	Graph convert_to_graph() const;

// 	SQPV modify_pitch( const Function<std::pair<Second, Pitch>, Pitch> & mod ) const;
// 	SQPV repitch( const Function<std::pair<Second, Pitch>, float> & mod ) const;

// 	SQPV select( Second length, const Function<std::pair<Second, Pitch>, std::pair<Second, Pitch>> & selector ) const;
// };

// }