#pragma once

namespace flan {

/*
While working on algorithms there were a lot of times I noticed that everything would be much faster if I could just modify the input sound.
Methods don't allow you to check if the input object is a temporary.
I could use normal functions with an overload for rvalue refs, but then I can't use the much better method syntax.
The solution to getting both behaviours is pipes. A pipe is an object that has some algorithm inputs bound to it.
You can then pass a single type into the pipe and get a transformed object of the same type out.
You can also compose pipes to get new pipes.
*/

template<typename T>
class Pipe {
public:
	using FuncType = std::function<T ( T && )>;

	template<typename C>
	requires std::convertible_to<C, FuncType>
	Pipe( C _f )
		: f( _f )
		{
		}

	T operator()( T && a ) const { return f( std::move( a ) ); }
	T operator()( const T & a ) const { return f( a.copy() ); }
	

	Pipe<T> operator()( const Pipe<T> & p ) const 
		{
		return [&]( T && a ){ return operator()( p( a ) ); };
		}

	Pipe<T> operator>>( const Pipe<T> & p ) const { return operator()( p ); }
	Pipe<T> operator<<( const Pipe<T> & p ) const { return p( *this ); }

private:
	FuncType f;
};
template<typename T>
T operator>>( const T & a, const Pipe<T> & p ) { return p(a); }
template<typename T>
T operator<<( const Pipe<T> & p, const T & a ) { return a >> p; }

}