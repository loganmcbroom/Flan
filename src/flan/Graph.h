#pragma once

#include <vector>
#include <functional>

#include "flan/defines.h"
#include "bmp/bitmap_image.hpp"

#include "flan/Utility/View.h"
#include "flan/Utility/vec2.h"
#include "flan/Utility/Color.h"
#include "flan/Function.h"

namespace flan {

/** Graph is a general graphing object. 
	This class stores image data in a std::shared_ptr.
*/
class Graph : public bitmap_image {
	public:
		Graph( const Graph & ) = delete;
		Graph & operator=( const Graph & ) = delete;
		Graph( Graph && ) = default;
		Graph & operator=( Graph && ) = default;
		~Graph() = default;

		struct Plane {
			int p;
			Plane( int i ) : p( i ) {}
			operator int() { return p; }
			bool operator==( Plane o ) const { return p == o.p; }
			static const Plane All;
			};

		/** This is used by waveform graphing functions.
		 *	Direct displays any sampled function data as-is.
		 *	Symmetric uses sample data magnitudes and displays the sampled data symmetrically around the waveform center line.
		 */
		enum class WaveformMode
			{
			Direct,
			Symmetric,
			};

		Graph( Pixel width = -1, Pixel height = -1 );

		//Graph & operator=( const Graph & other )
		//	{
		//	std::cout << "fuck you";
		//	return *this;
		//	}

		/** Add a view to the view list.
		 * \param view The view to add. This maps plane space to pixel space.
		 * \param plane The plane index of the view. -1 indicates all planes.
		 */
		void addView( View view, Plane plane = Plane::All );

		/** Set the view. All other views will be removed. The view will be active for all plane indices.
		 * \param view The plane area to map to all of pixel space.
		 */
		void setView( Rect view ); // note this will remove all other views

		/** Utility for splitting a view evenly along the y-axis into views on consecutive planes.
		 * \param view The view to add. This maps plane space to pixel space.
		 * \param numViews The number of views to split the provided view between.
		 * \param startPlane The plane index to start the split.
		 */
		void addSplitViewY( View view, int numViews, Plane startPlane = 0 );

		/** Utility for splitting the image evenly along the y-axis into views on consecutive planes.
		 * \param view The plane area to map to all of pixel space.
		 * \param numViews The number of views to split the provided view between.
		 * \param startPlane The plane index to start the split.
		 */
		void addFullSplitViewY( Rect view, int numViews, Plane startPlane = 0 );

		/** If p1 and p2 indicate the same plane, return true. Handles -1 indicating all planes.
		 */
		bool doPlanesMatch( Plane p1, Plane p2 ) const;

		/** Find all views which intersect a provided plane and Rect.
		 */
		std::vector<std::pair<Plane,View>> getIntersectingViews( Rect U, Plane plane ) const;


		//======================================================================================================================================================
		// Waveforms
		//======================================================================================================================================================

		/** Functional waveform graph.
		 * \param f The function to sample over [0,1]. Outputs will be clamped to [-1,1].
		 * \param rect The plane space to graph within.
		 */
		void drawWaveform( const Func1x1 & f, Rect rect = Rect(), Plane plane = -1, Color c = Color::White, WaveformMode mode = WaveformMode::Direct );

		/** Buffer waveform graph.
		 * \param data The buffer to sample. Outputs will be clamped to [-1,1].
		 * \param n The buffer size.
		 * \param rect The plane space to graph within.
		 */
		void drawWaveform( const float * data, int n, Rect rect = Rect(), Plane plane = -1, Color c = Color::White, WaveformMode mode = WaveformMode::Direct );

		/** Functional waveform graph. Functions will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param fs The functions to sample over [0,1]. Outputs will be clamped to [-1,1].
		 * \param rect The plane space to graph within.
		 */
		void drawWaveforms( const std::vector<Func1x1> & fs, Rect rect = Rect(), Plane startPlane = 0, WaveformMode mode = WaveformMode::Direct );

		/** Buffer waveform graph. Buffers will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param datas The buffers to sample. Outputs will be clamped to [-1,1].
		 * \param n The buffer size.
		 * \param rect The plane space to graph within.
		 */
		void drawWaveforms( const std::vector<const float *> & datas, int n, Rect rect = Rect(), Plane startPlane = 0, WaveformMode mode = WaveformMode::Direct );


		//======================================================================================================================================================
		// Spectrograms
		//======================================================================================================================================================

		/** Functional spectrum over time graph.
		 * \param f The function to sample over [0,1]X[0,1]. Outputs will be clamped to [0,1].
		 * \param rect The plane space to graph within.
		 * \param plane The plane to graph within.
		 */
		void drawSpectrogram( const Func2x1 & f, Rect rect = Rect(), Plane plane = Plane::All, float hue = 0 );

		/** Buffer spectrum over time graph.
		 * \param data The buffer to sample. This should have a time-major ordering. Outputs will be clamped to [-1,1].
		 * \param numFrames The number of time steps in the buffer.
		 * \param numBins The number of frequency steps in the buffer.
		 * \param rect The plane space to graph within.
		 */
		void drawSpectrogram( const float * data, int numFrames, int numBins, Rect rect = Rect(), Plane plane = Plane::All, float hue = 0 );

		/** Functional spectrum over time graph. Functions will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param fs The functions to sample over [0,1]X[0,1]. Outputs will be clamped to [0,1].
		 * \param rect The plane space to graph within.
		 */
		void drawSpectrograms( const std::vector<Func2x1> & fs, Rect rect = Rect(), Plane startPlane = 0 );

		/** Buffer spectrum over time graph. Buffers will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param datas The buffers to sample. These should have a time-major ordering. Outputs will be clamped to [-1,1].
		 * \param numFrames The number of time steps in the buffers.
		 * \param numBins The number of frequency steps in the buffers.
		 * \param rect The plane space to graph within.
		 */
		void drawSpectrograms( const std::vector<const float *> & datas, int numFrames, int numBins, Rect rect = Rect(), Plane startPlane = 0 );


		//======================================================================================================================================================
		// Functions
		//======================================================================================================================================================

		/** Function graphing.
		 * \param f The function to graph.
		 * \param domain The domain of the function.
		 */
		void drawFunction( const Func1x1 & f, Interval domain = Interval::R, Plane plane = Plane::All, Color c = Color::Black );

		/** Point interpolation graphing
		 * \param points The provided points are linearly interpolated and graphed with a domain exactly fitting the points.
		 */
		void drawFunction( const std::vector<std::pair<float,float>> & points, Plane plane = Plane::All, Color c = Color::Black );

		/** Function graphing. Functions will be graphed with maximally spaced hues.
		 * \param fs The functions to graph.
		 * \param domains The domains of the functions.
		 */
		void drawFunctions( const std::vector<Func1x1> & fs, const std::vector<Interval> & domains = {}, Plane plane = Plane::All );

		//======================================================================================================================================================
		// Primitives
		//======================================================================================================================================================

		// These all should in theory be acting on every view.
		// For effeciency, they don't bounds check. Appropriate views and bounds should be found using getIntersectingViews and iterated instead.

		void setPixel( Pixel x, Pixel y, Color c );
		void setPoint( const View & v, float x, float y, Color c );
		void drawHorizontalLine( const View & v, image_drawer & d, float x1, float x2, float y );
		void drawVerticalLine( const View & v, image_drawer & d, float y1, float y2, float x );
		void drawLineSegment( const View & v, image_drawer & d, float x1, float y1, float x2, float y2 );
		void setRect( const View & v, Rect rect, Color c );
		void fillImage( Color c );

		//======================================================================================================================================================
		// Additional Drawing Routines
		//======================================================================================================================================================

		void drawAxes( Plane plane = Plane::All, Color c = Color::Black );
		void drawLinearGrid( float xJumpSize = 1, float yJumpSize = 1, Plane plane = Plane::All, Color c = Color::fromHSV( 0, 0.0f, .7f ) );
		void drawXTicks( float jump, float y, Pixel offsetDown = 4, Pixel offsetUp = 4, Plane plane = Plane::All, Color c = Color::Black, float numberScale = 0 );
		void drawYTicks( float jump, float x, Pixel offsetLeft = 4, Pixel offsetRight = 4, Plane plane = Plane::All, Color c = Color::Black, float numberScale = 0 );

		/** Draw a single filled circular point, which has a size and shape independant of the view.
		 * \param p The point.
		 * \param radius The point radius in pixels.
		 */
		void drawPoint( const vec2 & p, Pixel radius = 6, Plane plane = Plane::All, Color c = Color::Black );

		/** Draw a set of points using drawPoint.
		 * \param ps The points.
		 * \param radius The point radius in pixels.
		 */
		void drawPoints( const std::vector<vec2> & ps, Pixel radius = 6, Plane plane = Plane::All, Color c = Color::Black );

		/** Very primitive text rendering. Seriously, I manually pathed out 0-9, dash, and dot, and draw it using lines. Looks cool though.
		 */
		void drawFloat( vec2 pos, Pixel digitWidth, Pixel digitHeight, float number,  Plane plane = Plane::All, Color c = Color::Black );

	private:
		std::vector<std::pair<Plane, View>> views;
		//float hue = 0;

		static uint32_t oversample;
		static Pixel DefaultWidth, DefaultHeight;
	};


inline uint32_t Graph::oversample = 4;
inline Pixel Graph::DefaultWidth = 1920;
inline Pixel Graph::DefaultHeight = 1080;

}