#pragma once

// #include <vector>
// #include <algorithm>
// #include <execution>

// namespace flan
// {

// template<typename T>
// class Buffer
// {
// public:
// 	Buffer( const Buffer & ) = delete;
// 	Buffer( Buffer && ) = default;
// 	Buffer& operator=( const Buffer & ) = delete;
// 	Buffer& operator=( Buffer && ) = default;
// 	~Buffer() = default;

// 	std::vector<std::views::take> channels()
// 		{
// 		std::vector<std::views::take> channel_vec( num_channels );
// 		for( const int i : std::views::iota( 0, num_channels ) )
// 			{
// 			channel_vec[i] = std::views::take
// 			}
// 		}

// 	std::vector<T>::iterator begin() { return buffer.begin(); }
// 	std::vector<T>::iterator end()   { return buffer.end();   }

// private:
// 	const size_t num_channels;
// 	const size_t channel_size;
// 	std::vector<T> buffer;
// };

// }