#pragma once

#include <cstdint>

namespace flan {

/** Basic struct for color information */
struct Color 
	{ 
	Color( char _red, char _green, char _blue ) : red( _red ), green( _green ), blue( _blue ) {}

    template< typename T >
    Color( T c ) : red( c[0] ), green( c[1] ), blue( c[2] ) {}

	uint8_t red, green, blue; 

    /** Convert HSV (Hue-Saturation-Value) data to RGB (Red-Green-Blue) data. */
    static Color fromHSV( int H, float S, float V );

    static const Color White;
    static const Color Black;
    static const Color Red;
    static const Color Green;
    static const Color Blue;
	};

inline const Color Color::White = Color::fromHSV( 0, 0, 1 );
inline const Color Color::Black = Color::fromHSV( 0, 0, 0 );
inline const Color Color::Red   = Color::fromHSV( 0, 1, 1 );
inline const Color Color::Green = Color::fromHSV( 120, 1, 1 );
inline const Color Color::Blue  = Color::fromHSV( 240, 1, 1 );

}