#ifndef TOOLS_H
#define TOOLS_H
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iomanip>
// #include <pair>
#include "point.h"
// #include "group.h"
#include "sys/stat.h"
#include "my_types.h"
#include "pair.h"
#include <sys/stat.h>
// #include <iosfwd>
#include <streambuf>
// #include <sstream>
// #include <fstream>
// #include <stdlib.h>

// #include "matrix.h"

typedef size_t 																	index_t;

typedef Pair < size_t > 												index_pair;
typedef Pair < int 		> 												signed_index_pair;


// typedef pair < index_t, index_t > 							index_pair;
// float 	my_fmod ( float a,float b );
// double 	my_fmod ( const double& a,const double& b );

// struct zerocopy_istringbuf : public std::stringbuf
// {
// 	zerocopy_istringbuf(std::string const& s);
// };




uint roundedFraction(float fraction, uint N);





template < class T >
T my_fmod ( const T& a, const T& b )
{
	if ( a > 0.0f ) return a - ( ( static_cast<int> ( a/b ) ) *b );
	else return a - ( ( static_cast<int> ( a/b )-1 ) *b );
}

template < class T1, class T2 > 
Point < T1 > my_fmod(const Point < T1 > & a, const Point < T2 > & b)
{
	Point < T1 > tmp;
	
	tmp.SetX( my_fmod(a.GetX() , typename Point < T1 >::value_type( b.GetX() )));
	tmp.SetY( my_fmod(a.GetY() , typename Point < T1 >::value_type( b.GetY() )));
	tmp.SetZ( my_fmod(a.GetZ() , typename Point < T1 >::value_type( b.GetZ() )));
	
	return tmp;
}

template < class T > 
uint is_negative(const T& val)
{
	return uint(val >> (sizeof(T)*8-1)) & 1;
}

template < class T > 
uint is_positive(const T& val)
{
	return is_negative(-val);
}

template < class T > 
uint is_zero(const T& val)
{
	return 1 - is_positive(val) - is_negative(val);
}

template < class T > 
uint not_zero(const T& val)
{
	return is_positive(val) + is_negative(val);
}

template < class T >
T periodic_space_projection(const T& val, const T& period )
{
  		T tmp = val;
 		int n = int (tmp / period);
		if ( tmp > 0 ) return val - T(n)*period;
		else return val - T(n-1)*period;
		
		
// 		T tmp(val);
// 		if (tmp == T(0)) return tmp;
// 		while (tmp >  period) tmp -= period;
// 		while (tmp <  0)   tmp += period;
// 		return tmp;
}

template < class T , class T2>
T periodic_space_distance(const T& a, const T& b, const T2& period, const T& half_period )
{
	T tmp(a-b);
	while (tmp >  half_period) tmp -= period;
	while (tmp < -half_period) tmp += period;
 	return tmp; 
}

template < class T >
T periodic_space_distance(const T& a, const T& b, const T& half_period )
{
	return periodic_space_distance(a,b,2*half_period,half_period);
}



template < class T, T half_period >
T periodic_space_distance(const T& a, const T& b )
{
	return periodic_space_distance(a,b,half_period);
}

template < class T >  Point < T >  periodic_space_distance 
( const Point < T > & a, const Point < T > & b, const Point < T > & half_period)
{
	return Point < T > (
		 periodic_space_distance(a.GetX(),b.GetX(),half_period.GetX())
		,periodic_space_distance(a.GetY(),b.GetY(),half_period.GetY())
		,periodic_space_distance(a.GetZ(),b.GetZ(),half_period.GetZ()));
}

template < class T >  Point < T >  periodic_space_distance 
( const Point < T > & a, const Point < T > & b, const Point < T > period, const Point < T > & half_period)
{
	return Point < T > (
		 periodic_space_distance(a.GetX(),b.GetX(),period.GetX(),half_period.GetX())
		,periodic_space_distance(a.GetY(),b.GetY(),period.GetY(),half_period.GetY())
		,periodic_space_distance(a.GetZ(),b.GetZ(),period.GetZ(),half_period.GetZ()));
}

class mytie
{
private:

  void* first;
  void* second;

public:

  template < class F, class S >
    mytie(F& _first, S& _second)
  :first  (static_cast < void* > (&_first))
  ,second (static_cast < void* > (&_second))
  { }

  template < class F, class S >
  void inline operator=( const std::pair < F,S >& p)
  {
    *(static_cast < F* > (first )) = p.first;
    *(static_cast < S* > (second)) = p.second;
  }
};

template < class Container >
std::pair < typename Container::iterator , typename Container::iterator >
iterators(Container& c)
{
  typedef typename Container::iterator iter_t;
  return std::pair < iter_t, iter_t > (c.begin(),c.end());
}

template < class Container >
std::pair
<
  typename Container::const_iterator ,
  typename Container::const_iterator
>
const_iterators(const Container& c)
{
  typedef typename Container::const_iterator iter_t;
  return std::pair < iter_t, iter_t > (c.begin(),c.end());
}

template < unsigned int l, class T = uint32_t > class l2mask
{
public:
	inline static T get() {return (((T(1u)) << l) - (T(1u)));}
  inline T operator () () const {return (((T(1u)) << l) - (T(1u)));}
};


template < unsigned int l > class l2mask < l, uint32_t >
{
	typedef uint32_t T;
public:
  enum {v = (uint32_t((uint64_t(1u) << l) - uint64_t(1u)))};
	inline static T get() {return (((T(1u)) << l) - (T(1u)));}
  inline T operator () () const {return (((T(1u)) << l) - (T(1u)));}
};


template < class T>
void printint(const T & value, unsigned length = sizeof(T)*8, unsigned unit = 8, std::ostream& stream = std::cout)
{
// 	typedef unsigned long long op_t;

	for (int i = length-1; i >= 0 ; i--)
	{
		stream << (((((T)1) << i) & (T) value) ? "o" : "-");
		if (i % unit == 0 && i) stream << " ";
	}
	stream << " : " << value << std::endl;
}

template < class T >
void find_heq_bit(const T& value, T& bit, T& shift)
{
	bit = 0; shift = 0;

	while (bit < value)
	{
		if (bit)
			bit <<= 1;
		else
			bit   = 1;

		shift++;
	}
}

unsigned lowest_mask(unsigned a);

template<unsigned int A>
struct Log2{
 enum{value = Log2<(A>>1)>::value + 1};
};


template<>
struct Log2<1>{
 enum{value = 0};
};

uint my_log2(uint v);

template < class T >
void Permutate_Sign(std::vector < Point < T > > & vec,const Point < T > & p
	)
{
	vec.push_back(Point < T > ( p.GetX(), p.GetY(), p.GetZ()));

	vec.push_back(Point < T > (-p.GetX(), p.GetY(), p.GetZ()));
	vec.push_back(Point < T > ( p.GetX(),-p.GetY(), p.GetZ()));
	vec.push_back(Point < T > ( p.GetX(), p.GetY(),-p.GetZ()));

	vec.push_back(Point < T > ( p.GetX(),-p.GetY(),-p.GetZ()));
	vec.push_back(Point < T > (-p.GetX(), p.GetY(),-p.GetZ()));
	vec.push_back(Point < T > (-p.GetX(),-p.GetY(), p.GetZ()));

	vec.push_back(Point < T > (-p.GetX(),-p.GetY(),-p.GetZ()));
}


	template < class T >
void Permutate(std::vector < Point < T > > & vec, const Point < T > & p)
{
	Permutate_Sign(vec, Point < T > ( p.GetX(), p.GetY(), p.GetZ()));
	Permutate_Sign(vec, Point < T > ( p.GetX(), p.GetZ(), p.GetY()));

	Permutate_Sign(vec, Point < T > ( p.GetY(), p.GetZ(), p.GetX()));
	Permutate_Sign(vec, Point < T > ( p.GetY(), p.GetX(), p.GetZ()));

	Permutate_Sign(vec, Point < T > ( p.GetZ(), p.GetX(), p.GetY()));
	Permutate_Sign(vec, Point < T > ( p.GetZ(), p.GetY(), p.GetX()));
}


template < class T >
std::ostream& operator << (std::ostream& stream, const std::vector<T> &v)
{
	typename std::vector<T>::const_iterator it;
	for (it = v.begin(); it < v.end(); it++)
	{
		stream << ' ' << *it;
	}
	return stream;
}


// namespace std
// {

template < class T, class Label >
class LabeledObject : public T
{
public:

// 	LabeledObject( const Label& _label, const T& val = T())
// 	:T(val),label(_label){}

// // 	typedef typename T
	typedef  Label label_type;

	LabeledObject(const T& val = T(), const Label& _label = Label(__default))
	:T(val),label(_label){}
	
	virtual ~LabeledObject() {}



	LabeledObject(const LabeledObject& rhs)
	:T(static_cast < T > (rhs)),label(rhs.label){}

	operator T() const {return static_cast < T > (*this);}
	operator T() 			 {return static_cast < T > (*this);}

	const Label& 	GetLabel() const 							{ return label;}
	void					SetLabel(const Label& _label) {label = _label;}

	static void	SetDefault(const Label& _default) {__default = _default;}
	static const Label& GetDefault() {return  __default;}

	LabeledObject& operator = (const LabeledObject& rhs)
	{
		static_cast < T > (*this) = static_cast < T > (rhs);
		label 			= rhs.label;
		return *this;
	}

	T& operator = (const T& rhs)
	{
		static_cast < T > (*this)	= rhs;
		return static_cast < T > (*this);
	}

private:

	Label 				label;
	static Label 	__default;

};

template < class T , class Label >
typename LabeledObject < T,Label >::label_type LabeledObject < T,Label >::__default;


// template < class T > class NamedObject : public LabeledObject < T, std::string >
// {
// private:
// 	typedef LabeledObject < T, std::string > Base;
//
// public:
// 	NamedObject(const std::string& _name, const T& val = T() )
// 	:Base(val,_name)
// 	{}
//
// 	const std::string& GetName() const 									{return 	Base::GetLabel();}
// 	void							 SetName(const std::string& _name){					Base::SetLabel(_name);}
//
// };

namespace std
{

template < class T > class NamedObject : public T
{

public:

	NamedObject(const string& _name, const T& val = T() ):T(val),name(_name){}

	NamedObject(const NamedObject& rhs):T(static_cast < T > (rhs)),name(rhs.name){}

	operator T() const {return static_cast < T > (*this);}
	operator T() 			 {return static_cast < T > (*this);}

	const string& GetName() const { return name;}
	void					SetName(const string& _name) {name = _name;}

	NamedObject& operator = (const NamedObject& rhs)
	{
		static_cast < T > (*this) = static_cast < T > (rhs);
		name 			= rhs.name;
		return *this;
	}

	T& operator = (const T& rhs)
	{
		static_cast < T > (*this)	= rhs;
		return static_cast < T > (*this);
	}

private:

	string name;

};

	template < class T1, class T2 >
	pair < T1,T2 > operator + (const pair < T1,T2 >& p1, const pair < T1,T2 >& p2)
	{
		return pair < T1,T2 > (p1.first+p2.first, p1.second+p2.second);
	}

	template < class T1, class T2 >
	pair < T1,T2 > operator - (const pair < T1,T2 >& p1, const pair < T1,T2 >& p2)
	{
		return pair < T1,T2 > (p1.first-p2.first, p1.second-p2.second);
	}


}
// template < class T > class NamedObject : public T
// {
// public:
// 	NamedObject(const string& _name, const T& val = T() ):T(val),name(_name){}
// 	NamedObject(const NamedObject& rhs):T(static_cast < T > (rhs)),name(rhs.name){}
//
// 	operator T() const {return static_cast < T > (*this);}
// 	operator T() 			 {return static_cast < T > (*this);}
//
// 	const string& GetName() const { return name;}
// 	void					SetName(const string& _name) {name = _name;}
//
// 	NamedObject& operator = (const NamedObject& rhs)
// 	{
// 		static_cast < T > (*this) = static_cast < T > (rhs);
// 		name 			= rhs.name;
// 		return *this;
// 	}
//
// 	T& operator = (const T& rhs)
// 	{
// 		static_cast < T > (*this)	= rhs;
// 		return static_cast < T > (*this);
// 	}
//
// private:
//
// 	string name;
//
// };

// };

// namespace std
// {
//
// template < class T1, class T2 >
// std::ostream& operator << (std::ostream& str, const std::pair < T1, T2 > & p)
// {
//   str << '(' << p.first << ',' << p.second << ')';
//   return str;
// }
//
// }

// template < class STL >
// print_stl(


// class print_container
// {
// public:
// 	template < class T >
// 	void operator() (const T& c, ostream& str = cout, char separator = ' ')
// 	{
// 		typename T::iterator it;
// 		for (it = c.begin(); it != c.end(); it++)
// 		{
// 			str << ' ' << *it;
// 		}
// 	}
// };


// namespace smart
// {
// template <  class Container >
// std::ostream& operator << (std::ostream& stream, Container &c)
// {
// 	typename Container::iterator it;
//
// 	for (it = c.begin(); it < c.end(); it++)
// 	{
// 		stream << ' ' << *it;
// 	}
//
// 	return stream;
// }

// }

template
<
	class T
// 	template < class T > class vector_T
>
int GetIndex(const std::vector < T > & v, const T& element)
{
	return (&element - &v[0]);
}

// template < class T >
// std::ostream& operator << (std::ostream& stream, const std::pair<T,T> & p)
// {
// 	stream << p.first() << " " <<  p.seccond();
// 	return stream;
// }

// template < class First, class Second, class Result >
// class my_binary
// {
// public:
// 	typedef First 	first_argument_type;
// 	typedef Second 	second_argument_type;
// 	typedef Result 	result_type;
// };


template < class T1, class T2 >
class printpair
{
public:

  template < class t1, class t2 >
  printpair (const std::pair<t1,t2>& src):ref(std::pair<T1,T2>(T1(src.first),T2(src.second))){ }

	void operator () (std::ostream& stream) const
	{
		stream  << '(' <<  ref.first << "|" <<  ref.second << ')';
	}

private:

	const std::pair<T1,T2>& ref;

};

template < class T1, class T2 >
std::ostream& operator << (std::ostream& stream, const printpair < T1, T2 >& print_obj)
{
	print_obj(stream); return stream;
}


// template < class point_T >
//     void permutate_sign(const point_T& p, std::set < point_T > & set)
// {
//   set.insert(point_T ( p.X, p.Y, p.Z));
//
//   set.insert(point_T (-p.X, p.Y, p.Z));
//   set.insert(point_T ( p.X,-p.Y, p.Z));
//   set.insert(point_T ( p.X, p.Y,-p.Z));
//
//   set.insert(point_T ( p.X,-p.Y,-p.Z));
//   set.insert(point_T (-p.X, p.Y,-p.Z));
//   set.insert(point_T (-p.X,-p.Y, p.Z));
//
//   set.insert(point_T (-p.X,-p.Y,-p.Z));
// }
//
// template < class point_T >
//     void permutate_point(const point_T & p, std::set < point_T > & set)
// {
//   permutate_sign(set, point_T ( p.X, p.Y, p.Z));
//   permutate_sign(set, point_T ( p.X, p.Z, p.Y));
//
//   permutate_sign(set, point_T ( p.Y, p.Z, p.X));
//   permutate_sign(set, point_T ( p.Y, p.X, p.Z));
//
//   permutate_sign(set, point_T ( p.Z, p.X, p.Y));
//   permutate_sign(set, point_T ( p.Z, p.Y, p.X));
// }

// template < class Map >
// void insert_bond(Map& table, const Bond& bond, cuint value = 1)
// {
//   std::set < Point < int >  >  bonds(48);
//
//   typedef std::set < Point < int >  > bondset;
//
//   Bond dummy;
//
//   permutate_point( Point < int >(bond.X, bond.Y, bond.Z), bonds);
//
//   for ( bondset::iterator it = bonds.begin(); it != bonds.end(); ++it )
//   {
//     dummy = *it;
//     table[dummy.GetPosData()] = value;
//   }
// }


/*
template < class T >
vector < T > arrange(T start, T step, unsigned N)
{
	vector < T > arrangement;
	T actual(start);

	for (unsigned i = 0; i < N ; i++)
	{
		arrangement.push_back(actual);
		actual+=step;
	}

	return arrangement;
}
*/


/*
template < class T >
vector < vector < T > > Transpose(const vector < vector < T > > & v)
{
	unsigned size(v.front().size());
	vector < vector < T > > new_v(size);

	typename vector < vector < T > >::iterator new_row;
	for (typename vector < vector < T > >::const_iterator old_row = v.begin();
			old_row < v.end();
			old_row++)
	{
		if (old_row->size() != size)
		{
			cout << "error: vector transposing - row not aligned.";
		}
		else
		{
			new_row = new_v.begin();
			for (typename vector < T >::const_iterator old_column = old_row->begin(); old_column < old_row->end(); old_column++)
			{
				new_row->push_back(*old_column);
				new_row++;
			}
		}
	}
	return new_v;
}
*/

template < class T >
std::vector < T >& operator += (std::vector < T >& lhs, const std::vector<T> &rhs)
{
	typename std::vector<T>::iterator li(lhs.begin());
	typename std::vector<T>::const_iterator ri(rhs.begin()),le(lhs.end());
	
	while(li != le)
	{*li += *ri; ++li; ++ri;}
	return lhs;
}

template < class T >
std::vector < T >& operator -= (std::vector < T >& lhs, const std::vector<T> &rhs)
{
	typename std::vector<T>::iterator li(lhs.begin());
	typename std::vector<T>::const_iterator ri(rhs.begin()),le(lhs.end());
	
	while(li != le)
	{*li -= *ri; ++li; ++ri;}
	return lhs;
}

template < class T , class T2>
std::vector < T >& operator *= (std::vector < T >& lhs, const T2& val)
{
	typename std::vector<T>::iterator li(lhs.begin());
	typename std::vector<T>::const_iterator le(lhs.end());
	while(li != le){*li *= val;++li;}
	return lhs;
}


bool FileExists(const std::string& strFilename);

bool file_exists(const std::string& name);


void findFile(const std::string& prefix, std::string& out );

#endif


