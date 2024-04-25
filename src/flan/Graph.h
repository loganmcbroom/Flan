#pragma once

#include <vector>
#include <functional>

#include "flan/defines.h"
#include "bmp/bitmap_image.hpp"

#include "flan/Utility/View.h"
#include "flan/Utility/vec2.h"
#include "flan/Utility/Color.h"

namespace flan {

template< typename I, typename O >
struct Function;

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
		void add_view( View view, Plane plane = Plane::All );

		/** Set the view. All other views will be removed. The view will be active for all plane indices.
		 * \param view The plane area to map to all of pixel space.
		 */
		void set_view( Rect view );

		/** Utility for splitting a view evenly along the y-axis into views on consecutive planes.
		 * \param view The view to add. This maps plane space to pixel space.
		 * \param num_views The number of views to split the provided view between.
		 * \param start_plane The plane index to start the split.
		 */
		void add_split_view_y( View view, int num_views, Plane start_plane = 0 );

		/** Utility for splitting the image evenly along the y-axis into views on consecutive planes.
		 * \param view The plane area to map to all of pixel space.
		 * \param num_views The number of views to split the provided view between.
		 * \param start_plane The plane index to start the split.
		 */
		void add_full_split_view_y( Rect view, int num_views, Plane start_plane = 0 );

		/** If p1 and p2 indicate the same plane, return true. Handles -1 indicating all planes.
		 */
		bool do_planes_match( Plane p1, Plane p2 ) const;

		/** Find all views which intersect a provided plane and Rect.
		 */
		std::vector<std::pair<Plane,View>> get_intersecting_views( Rect U, Plane plane ) const;


		//======================================================================================================================================================
		// Waveforms
		//======================================================================================================================================================

		/** Functional waveform graph.
		 * \param f The function to sample over [0,1]. Outputs will be clamped to [-1,1].
		 * \param rect The plane space to graph within.
		 */
		void draw_waveform( 
			const Function<Second, Sample> & f, 
			Rect rect = Rect(), 
			Plane plane = -1, 
			Color c = Color::White, 
			WaveformMode mode = WaveformMode::Direct, 
			uint32_t oversample = 4
			);

		/** Buffer waveform graph.
		 * \param data The buffer to sample. Outputs will be clamped to [-1,1].
		 * \param n The buffer size.
		 * \param rect The plane space to graph within.
		 */
		void draw_waveform( 
			const float * data, 
			int n, 
			Rect rect = Rect(), 
			Plane plane = -1, 
			Color c = Color::White, 
			WaveformMode mode = WaveformMode::Direct, 
			uint32_t oversample = 4 
			);

		/** Functional waveform graph. Functions will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param fs The functions to sample over [0,1]. Outputs will be clamped to [-1,1].
		 * \param rect The plane space to graph within.
		 */
		void draw_waveforms( 
			const std::vector<Function<float, float>> & fs, 
			Rect rect = Rect(), 
			Plane start_plane = 0, 
			WaveformMode mode = WaveformMode::Direct, 
			uint32_t oversample = 4 
			);

		/** Buffer waveform graph. Buffers will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param datas The buffers to sample. Outputs will be clamped to [-1,1].
		 * \param n The buffer size.
		 * \param rect The plane space to graph within.
		 */
		void draw_waveforms( 
			const std::vector<const float *> & datas, 
			int n, Rect rect = Rect(), 
			Plane start_plane = 0, 
			WaveformMode mode = WaveformMode::Direct, 
			uint32_t oversample = 4 
			);


		//======================================================================================================================================================
		// Spectrograms
		//======================================================================================================================================================

		/** Functional spectrum over time graph.
		 * \param f The function to sample over [0,1]X[0,1]. Outputs will be clamped to [0,1].
		 * \param rect The plane space to graph within.
		 * \param plane The plane to graph within.
		 */
		void draw_spectrogram( const Function<vec2, float> & f, Rect rect = Rect(), Plane plane = Plane::All, float hue = 0, uint32_t oversample = 4 );

		/** Buffer spectrum over time graph.
		 * \param data The buffer to sample. This should have a time-major ordering. Outputs will be clamped to [-1,1].
		 * \param num_frames The number of time steps in the buffer.
		 * \param num_bins The number of frequency steps in the buffer.
		 * \param rect The plane space to graph within.
		 */
		void draw_spectrogram( const float * data, int num_frames, int num_bins, Rect rect = Rect(), Plane plane = Plane::All, float hue = 0, uint32_t oversample = 4 );

		/** Functional spectrum over time graph. Functions will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param fs The functions to sample over [0,1]X[0,1]. Outputs will be clamped to [0,1].
		 * \param rect The plane space to graph within.
		 */
		void draw_spectrograms( const std::vector<Function<vec2, float>> & fs, Rect rect = Rect(), Plane start_plane = 0, uint32_t oversample = 4 );

		/** Buffer spectrum over time graph. Buffers will be graphed from bottom to top, splitting the rect evenly, and with maximally distant hues.
		 * \param datas The buffers to sample. These should have a time-major ordering. Outputs will be clamped to [-1,1].
		 * \param num_frames The number of time steps in the buffers.
		 * \param num_bins The number of frequency steps in the buffers.
		 * \param rect The plane space to graph within.
		 */
		void draw_spectrograms( 
			const std::vector<const float *> & datas, 
			int num_frames, 
			int num_bins, 
			Rect rect = Rect(), 
			Plane start_plane = 0, 
			uint32_t oversample = 4 
			);


		//======================================================================================================================================================
		// Functions
		//======================================================================================================================================================

		/** Function graphing.
		 * \param f The function to graph.
		 * \param domain The domain of the function.
		 */
		void draw_function( const Function<float, float> & f, Interval domain = Interval::R, Plane plane = Plane::All, Color c = Color::Black );

		/** Point interpolation graphing
		 * \param points The provided points are linearly interpolated and graphed with a domain exactly fitting the points.
		 */
		void draw_function( const std::vector<vec2> & points, Plane plane = Plane::All, Color c = Color::Black );

		/** Function graphing. Functions will be graphed with maximally spaced hues.
		 * \param fs The functions to graph.
		 * \param domains The domains of the functions.
		 */
		void draw_functions( const std::vector<Function<float, float>> & fs, const std::vector<Interval> & domains = {}, Plane plane = Plane::All );

		//======================================================================================================================================================
		// Primitives
		//======================================================================================================================================================

		// These all should in theory be acting on every view.
		// For effeciency, they don't bounds check. Appropriate views and bounds should be found using get_intersecting_views and iterated instead.

		void set_pixel( Pixel x, Pixel y, Color c );
		void set_point( const View & v, float x, float y, Color c );
		void draw_horizontal_line( const View & v, image_drawer & d, float x1, float x2, float y );
		void draw_vertical_line( const View & v, image_drawer & d, float y1, float y2, float x );
		void draw_line_segment( const View & v, image_drawer & d, float x1, float y1, float x2, float y2 );
		void set_rect( const View & v, Rect rect, Color c );
		void fill_image( Color c );

		//======================================================================================================================================================
		// Additional Drawing Routines
		//======================================================================================================================================================

		void draw_axes( Plane plane = Plane::All, Color c = Color::Black );

		void draw_linear_grid_x( float x_jump_size = 1, Plane plane = Plane::All, Color c = Color::from_hsv( 0, 0.0f, .7f ) );
		void draw_linear_grid_y( float y_jump_size = 1, Plane plane = Plane::All, Color c = Color::from_hsv( 0, 0.0f, .7f ) );
		void draw_linear_grid( float x_jump_size = 1, float y_jump_size = 1, Plane plane = Plane::All, Color c = Color::from_hsv( 0, 0.0f, .7f ) );

		void draw_x_ticks( 
			float jump, 
			float y, 
			float scale_base = 1.0f, 
			Pixel offset_down = 4, 
			Pixel offset_up = 4, 
			Plane plane = Plane::All, 
			Color c = Color::Black, 
			float number_scale = 0 
			);
		void draw_y_ticks( 
			float jump, 
			float x, 
			float scale_base = 1.0f,
			Pixel offset_left = 4, 
			Pixel offset_right = 4, 
			Plane plane = Plane::All, 
			Color c = Color::Black, 
			float number_scale = 0 
			);
		
		void draw_log_grid_x( float x_jump_size = 1, uint32_t lines_per_step = 10, Plane plane = Plane::All, Color c = Color::from_hsv( 0, 0.0f, .7f ) );
		void draw_log_grid_y( float y_jump_size = 1, uint32_t lines_per_step = 10, Plane plane = Plane::All, Color c = Color::from_hsv( 0, 0.0f, .7f ) );

		/** Draw a single filled circular point, which has a size and shape independant of the view.
		 * \param p The point.
		 * \param radius The point radius in pixels.
		 */
		void draw_point( const vec2 & p, Pixel radius = 6, Plane plane = Plane::All, Color c = Color::Black );

		/** Draw a set of points using draw_point.
		 * \param ps The points.
		 * \param radius The point radius in pixels.
		 */
		void draw_points( const std::vector<vec2> & ps, Pixel radius = 6, Plane plane = Plane::All, Color c = Color::Black );

		/** Very primitive text rendering. Seriously, I manually pathed out 0-9, dash, and dot, and draw it using lines. Looks cool though.
		 */
		void draw_float( vec2 pos, Pixel digit_width, Pixel digit_height, float number,  Plane plane = Plane::All, Color c = Color::Black );

	private:
		std::vector<std::pair<Plane, View>> views;
		//float hue = 0;

		static Pixel default_width, default_height;
	};

inline Pixel Graph::default_width = 1920;
inline Pixel Graph::default_height = 1080;

}