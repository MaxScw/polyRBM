#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <math.h>
// #include "moves.h"
// #include "bondvector.h"
#include <stdint.h>
#include <sys/types.h>


// typedef hash < int > H;

template < class T >
class Abstract_Point
{
public:
	virtual T GetX() const = 0;
	virtual T GetY() const = 0;
	virtual T GetZ() const = 0;
};

/**
 * @brief A 3D point class of arbitrary primive type.
 * 
 * The point contains three coordinates x, y and z and one 
 * component, w, as metadata. Point<T> bundles simple arithmetic operations
 * affecting three coordinates x,y and z. The fourth component w can be used for
 * meta data and leads to an improved layout in memory especially for vectors
 * of Point<T>s .
 * 
 * @param T builtin-type forming the x, y an z component as well as the
 * metadata w.
**/
template <class T>
class Point
{
private:
// 	static const Moves < Point < T > >	moves;
    static std::string coordSep;

public:
    
	typedef T value_type;
    
    static void setSeparator(const std::string& sep){coordSep=sep;} 
    static const std::string& getSeparator(){return coordSep;} 
    
	union
	{
		struct
		{
 			/// @var X,Y,Z,W
			T X,Y,Z,W;
		};
		T v[4];	// possibility to access data array-like;
	};

	/// Get the ith component i=0 --> X, ..., i=3 --> W .
	const T& Get(unsigned i)	const {return v[i];}
	
	const T& GetX() 			const {return X;}
	const T& GetY() 			const {return Y;}
	const T& GetZ() 			const {return Z;}
	const T& GetW() 			const {return W;}
	
	/// Analog to Get(i).
	void Set(unsigned i, const T& val){v[i] = val;}
	void SetX(const T& x) {X = x;}
	void SetY(const T& y) {Y = y;}
	void SetZ(const T& z) {Z = z;}
	void SetW(const T& z) {W = z;}

	/// Components can be accessed directly by using Point < T > as an array.
	const T& operator [](unsigned i)const	{return  v[i];}

	/// The non-const version of operator[], probably it will/should be removed again.
	T& operator [](unsigned i)			{return  v[i];}

	/// Default constructor sets everything zero, including w.
	Point():X(0),Y(0),Z(0),W(0){}

// 	Point():X(0),Y(0),Z(0){}

// 	Point(const src_T& x, const src_T& y, const src_T& z)
// 	{
// 		SetX( (T) x );
// 		SetY( (T) y );
// 		SetZ( (T) z );
// 	}

// 	template <class src_T>
// 	Point(const src_T& val)
// 	{
// 		SetX( (T) val );
// 		SetY( (T) val );
// 		SetZ( (T) val );
// 	}

	/// Single T value constructor : every component, including w, is set to val.
	Point(const T& val):X(val),Y(val),Z(val),W(val){}
// 	Point(const T& val)
// 	:X(val),Y(val),Z(val),dummy(val){}
	
// 	Point(unsigned val):X(val),Y(val),Z(val),W(val){}
// 	Point(long unsigned val):X(val),Y(val),Z(val),W(val){}

	/** 
	* 	Single-parameter constructor for arbitrary Point-like classes. 
	*  The results of src.GetX(), ..., src.GetW() are used to initialize the own components.
	* 
	*  @param PT This can be an arbitrary class providing
	*  GetX(), ..., GetW() getters, such as Point<T>, BitPointDynamic or others written by the user.
	*  The return types of the getters are explicitly casted to T.
	* 
	**/
	template < class PT >
	Point(const PT & src):X(T ( src.GetX() )),Y ( T(src.GetY()) ), Z(T(src.GetZ())), W(T(src.GetW())){}
	
// 	Point(const store_& val)
// 	template < class T2 >
// <<<<<<< HEAD
	
// =======
// >>>>>>> 5c23fe699658a6829d195e047e2f442c1ce10337
	
	/// 3D constructor, w is set to T(0) .
	template <class T1, class T2, class T3 >
	Point(const T1& x, const T2& y, const T3& z)
	{
		SetX( T( x ));
		SetY( T( y ));
		SetZ( T( z ));
		SetW( T( 0 ));
	}

	/// 4D constructor
	template <class T1, class T2, class T3, class T4 >
	Point(const T1& x, const T2& y, const T3& z, const T4& d)
	{
		SetX( T( x ));
		SetY( T( y ));
		SetZ( T( z ));
		SetW( T( d ));
	}

	
// 	template <class src_T>
// 	Point(const Point < src_T > & src)
// 	{
// 		SetX( (T) src.GetX() );
// 		SetY( (T) src.GetY() );
// 		SetZ( (T) src.GetZ() );
// 		SetW( (T) src.GetW() );
// 	}

// 	template <class src_T>
// 	Point(const Abstract_Point < src_T > & src)
// 	{
// 		SetX( (T) src.GetX() );
// 		SetY( (T) src.GetY() );
// 		SetZ( (T) src.GetZ() );
// 	}

// 	template < class PT >
// 	Point& operator = (const PT& rhs)
// 	:{return *this;}

	
	
	template <class src_T>
	Point<T> operator = (const Point<src_T>	& src)
	{
		SetX(src.GetX());
		SetY(src.GetY());
		SetZ(src.GetZ());
		SetW(src.GetW());
		return *this;
	}

	template <class src_T>
	void Assign(const src_T& x, const src_T& y, const src_T& z, const src_T& w = 0)
	{
		SetX( (T) x );
		SetY( (T) y );
		SetZ( (T) z );
		SetW( (T) w );
	}

	/// Scalar product returning a sum of type T. 
	/// Components from src are converted to T before multiplying.
	/// @param src_T Underlying data type of src.
	/// @param src   object of class Point < src_T >
	template <class src_T>
	T operator * (const Point < src_T > & src) const
	{
		T sum = 0;
		sum  = X*( T (src.GetX()));
		sum += Y*( T (src.GetY()));
		sum += Z*( T (src.GetZ()));
		return sum;
	}
// 	template <class src_T>

	/// returns a newly constructed Point object as copy of the left hand side, 
	/// which is multiplied by @param factor. 
	Point operator * (const T& factor) const
	{
		Point result(*this);
		result *= factor;
		return result;
	}


// 	template <class src_T>
// 	template <class src_T>
// 	Point<T> 	operator +  (const Point<src_T>	& src) const
// 	{
// 		return Point < T > (X+(T)src.X,Y+(T)src.Y,Z+(T)src.Z);
// 	}

	template <class src_T>
	Point<T> 	operator +  (const Point<src_T>	& src) const
	{
		Point tmp(GetX()+(T)src.GetX(),GetY()+(T)src.GetY(),GetZ()+(T)src.GetZ());
		return tmp;
	}

	template <class src_T>
	Point<T> 	operator -  (const Point<src_T>	& src) const
	{
		Point tmp(GetX()-(T)src.GetX(),GetY()-(T)src.GetY(),GetZ()-(T)src.GetZ());
		return tmp;
	}

// 	template < class PT >
// 	void operator += (const PT & src)
// 	{
// 		X += (T) src.X;
// 		Y += (T) src.Y;
// 		Z += (T) src.Z;
// 	}

	template <class src_T>
	void operator -= (const Point <src_T>  	& src)
	{
		X -= (T) src.X;
		Y -= (T) src.Y;
		Z -= (T) src.Z;
	}

// 	template < class _T >
// 	Point& operator += (const Abstract_Point < _T >& rhs)
// 	{
// 		X += rhs.GetX();
// 		Y += rhs.GetY();
// 		Z += rhs.GetZ();
// 	}

	bool operator > (const Point& rhs) const
	{
		T sqr_lhs((*this)*(*this)), sqr_rhs(rhs*rhs);
		return sqr_lhs > sqr_rhs;
	}

	bool operator < (const Point& rhs) const
	{
		T sqr_lhs((*this)*(*this)), sqr_rhs(rhs*rhs);
		return sqr_lhs < sqr_rhs;
	}

	template < class PT>
	void operator += (const PT& rhs)
	{
		X += T(rhs.GetX());
		Y += T(rhs.GetY());
		Z += T(rhs.GetZ());
	}

	template < class PT >
	void operator -= (const PT& rhs)
	{
		X -= T(rhs.GetX());
		Y -= T(rhs.GetY());
		Z -= T(rhs.GetZ());
	}

	template <class factor_T>
	Point&	operator *= (const factor_T  	& factor)
	{
		X = T((factor_T) X * factor);
		Y = T((factor_T) Y * factor);
		Z = T((factor_T) Z * factor);
		return *this;
	}
	
	template <class factor_T>
	Point&	operator /= (const factor_T  & factor)
	{
		X = T((factor_T) X / factor);
		Y = T((factor_T) Y / factor);
		Z = T((factor_T) Z / factor);
		return *this;
	}

	template <class src_T>
	bool 	operator == (const Point <src_T>  	& src) const
	{
		if (GetX() != (T) src.GetX()) return false;
		if (GetY() != (T) src.GetY()) return false;
		if (GetZ() != (T) src.GetZ()) return false;
		return true;
	}

	template <class src_T>
	bool 	operator != (const Point <src_T>  	& src) const
	{
		return !(*this == src);
	}

//   Point operator - (const Point& bp)
//   {
//     return Point(-bp.X,-bp.Y,-bp.Z);
//   }


// 	inline void Add_Move(const Move < Point < T > > & move)
// 	{
// 		v[move.GetAxis()] += move.GetSummand().v[move.GetAxis()];
// 	}

// 	inline void operator += (const Move < Point < T > > & move)
// 	{
// 		v[move.GetAxis()] += move.GetSummand().v[move.GetAxis()];
// 	}

	double Length() const
	{
		T 			sqr_norm((*this)*(*this));
		return 	sqrt((double)sqr_norm);
	}

	double Normalize();
};

template <class T>
double Point < T >::Normalize()
{
	double length(Length()), inverse_length(1.0/length);

	if(!isnan(inverse_length))
	{
		X*=inverse_length;
		Y*=inverse_length;
		Z*=inverse_length;
	}

	return length;
}

// template <class T>
// Point <double> Normalize(T & p)
// {
// 	T norm(p*p);
// 	double inverse_sqrt_norm(1.0/sqrt((double)norm));
// 	Point <double> normalized(	(double)p.GetX()*inverse_sqrt_norm,
// 								(double)p.GetY()*inverse_sqrt_norm,
// 								(double)p.GetZ()*inverse_sqrt_norm);
// 	if(isnan(normalized.GetX()))
// 		std::cout << "vector leading to nan: (" << p << ")\n";
// 	return normalized;
// }

template <class T>
T CrossProduct(const T & a, const T b)
{
	T	 prod(	a.Y*b.Z - a.Z*b.Y,
						a.Z*b.X - a.X*b.Z,
						a.X*b.Y - a.Y*b.X	);
	if (isnan(prod.X))
	{
		std::cout << "vectors leading to nan: (" << a << ") (" << b << ")" << std::endl;
		return T(0.0,0.0,0.0);
	}
	else return prod;
}

// template < class T >
// std::ostream& operator << (std::ostream& stream, const Abstract_Point < T > & p)
// {
// 	stream << p.GetX() << " " <<  p.GetY() << " " << p.GetZ();
// 	return stream;
// }
//
// template < class T >
// std::istream& operator >> (std::istream& stream, Abstract_Point < T > & p)
// {
// 	T temp;
// 	stream >> temp; p.SetX(temp);
// 	stream >> temp; p.SetY(temp);
// 	stream >> temp; p.SetZ(temp);
// 	return stream;
// }

template < class T >
Point < double > operator * (const double& f, const Point < T >& p)
{
	Point <double> result(p);
	result *= f;
	return result;
}

template < class T >
Point < double > operator * (const Point < T >& p, const double& f)
{
	Point <double> result(p);
	result *= f;
	return result;
}

template < class T >
std::ostream& operator << (std::ostream& stream, const Point < T > & p)
{
	stream << p.GetX() << Point < T >::getSeparator() <<  p.GetY() << Point < T >::getSeparator() << p.GetZ();
	return stream;
}

template < class T >
std::istream& operator >> (std::istream& stream, Point < T > & p)
{
	T temp;
	stream >> temp; p.SetX(temp);
	stream >> temp; p.SetY(temp);
	stream >> temp; p.SetZ(temp);
	return stream;
}


typedef Point < int  >     IPoint;
typedef Point < uint >     UPoint;
typedef Point < int64_t >  ILPoint;
typedef Point < uint64_t > ULPoint;
typedef Point < float >    FPoint;
typedef Point < double >   DPoint;

// template < class T, class Alloc > class std::vector;

#include <vector>

typedef std::vector < UPoint > UPointVec;
typedef std::vector < IPoint > IPointVec;
typedef std::vector < FPoint > FPointVec;
typedef std::vector < DPoint > DPointVec;

// class LPoint;
// typedef std::vector < LPoint > LPointVec;

// typedef std::vector < DPoint > DPointVec;

#endif

