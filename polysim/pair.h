#ifndef PAIR_H
#define PAIR_H

#include <iostream>
#include <math.h>
#include "stdint.h"

typedef unsigned uint;

template <class T>	class Pair
{
public:
	union
	{
		struct
		{
			union
			{
				T X;
				T first;
			};
			union
			{
				T Y;
				T second;
			};			
		};
		
		T v[2];	// possibility to access same data array-like;
	};

	Pair():X(0),Y(0){}

	Pair(const T& a, const T& b):X(a),Y(b){ }
	
	Pair(const T& a):X(a),Y(a){ }
// 	Pair(T& a, T& b):X(a),X(b){ }

	template <class src_T>
	Pair(const Pair<src_T> & src)
	:X((T)src.X)
	,Y((T)src.Y)
	{}

	template <class src_T>
	Pair(src_T x,src_T y)
	:X((T)x)
	,Y((T)y)
	{}

	Pair(std::istream& stream)
	{
		stream >> *this;
	}

	template <class src_T>
	void Assign(const src_T& x, const src_T& y)
	{
		X = (T)x; Y = (T)y;
	}

	template <class src_T>
	Pair<T>& 	operator  = (const Pair < src_T > & src)
	{
		X = (T)src.X; Y = (T)src.Y; return *this;
	}

	template <class src_T>
	Pair<T>&	operator  = (const src_T	 	& src)
	{
		T copy((T)src);
		X = Y = copy;
		return *this;
	}

	template <class src_T>
	T	operator *  (const Pair<src_T>	& src) const
	{
		T sum = 0;
		sum  = X*((T)src.X);
		sum += Y*((T)src.Y);
		return sum;
	}

	template < class src_T >
	Pair < T > operator *
		 (const src_T  & src) const
	{
		Pair < T > product;
		product.X = X*((T)src);
		product.Y = Y*((T)src);
		return product;
	}

	template <class src_T>
	Pair<T> 	operator +  (const Pair < src_T > & src)
	{
		Pair<T> sum;
		sum.X = X + (T)src.X;
		sum.Y = Y + (T)src.Y;
		return sum;
	}

	template <class src_T>
	Pair<T> 	operator -  (const Pair < src_T > & src) const
	{
		Pair<T> sum;
		sum.X = X - (T)src.X;
		sum.Y = Y - (T)src.Y;
		return sum;
	}

	template <class src_T>
	void operator += (const Pair <src_T>  	& src)
	{

		X += (T)src.X;	Y += (T)src.Y;
	}

	template <class src_T>
	void operator %= (const Pair <src_T>  	& src)
	{

		X %= (T)src.X;	Y %= (T)src.Y;
	}
	
	template <class src_T>
	void operator /= (const Pair <src_T>  	& src)
	{

		X /= (T)src.X;	Y /= (T)src.Y;
	}
	
	template <class scalarT>
	void operator /= (const scalarT & src)
	{

		X /= (T)src;	Y /= (T)src;
	}


	template <class src_T>
	Pair<T>& 	operator -= (const Pair<src_T>  	& src)
	{
		X -= (T)src.X;	Y -= (T)src.Y; return *this;
	}

	template <class factor_T>
	Pair<T>&	operator *= (const factor_T  	& factor)
	{
		X = (T) (factor * (factor_T) X );
		Y = (T) (factor * (factor_T) Y );
		return *this;
	}

	template <class src_T>
	bool 	operator == (const Pair <src_T>  	& src) const
	{
		if (X != (T) src.X) return false;
		if (Y != (T) src.Y) return false;
		return true;
	}

	template <class src_T>
	bool 	operator != (const Pair <src_T>  	& src) const
	{
		if (X != (T) src.X) return true;
		if (Y != (T) src.Y) return true;
		return false;
	}

	bool CheckComponentLength(unsigned max)
	{
		if ((X<=(T)max)&&(X>=(T)(-max)) && (Y<=(T)max)&&(Y>=(T)(-max)))
			return true;
		else
			return false;
	}

	void Normalize()
	{
		T norm((*this)*(*this));
		double inverse_sqrt_norm(1.0/sqrt((double)norm));

		if(!isnan(inverse_sqrt_norm))
		{
			X *= inverse_sqrt_norm;
			Y *= inverse_sqrt_norm;
		}
	}
};

template < class T >
std::ostream& operator << (std::ostream& stream, const Pair<T> & p)
{
	stream << p.X << " " <<  p.Y;
	return stream;
}

template < class T >
bool operator < (const Pair <T> & _X, const Pair <T> & _Y)
{
	return  (_X.first < _Y.first || 
					(!( _Y.first < _X.first) && _X.second < _Y.second));
}

template < class T >
std::istream& operator >> (std::istream& stream, Pair<T> & p)
{
	stream >> p.X >> p.Y;
	return stream;
}

typedef Pair < int  >     IPair;
typedef Pair < unsigned >     UPair;
typedef Pair < int64_t >  ILPair;
typedef Pair < uint64_t > ULPair;
typedef Pair < float >    FPair;
typedef Pair < double >   DPair;


#endif
