#ifndef NUMERICS_H
#define NUMERICS_H

#include "point.h"
#include <cmath>
#include <vector>
#include "static_vector.h"
#include <deque>
#include <algorithm>
#include "matrix.h"

typedef Point < double > PD;
typedef MatrixStatic < double , 3, 4 > M3D;

#define TOL 1e-7

class TensorOrder : public M3D
{
  typedef Point < double > input_type; 

public:
  
  TensorOrder();
  TensorOrder& operator+=(const input_type& val);
  void reset();
  
  const DPoint& getDirector()const;
  const double& getOrderParameter()const;
  void Print() const;
  
private:
  
  bool fresh;
  MatrixStatic < double , 3, 4 > sum;
  DPoint director; // fourth, W, component will contain order parameter...
  uint64_t counter;
  void finish();
  void solve();
  
};

template < class Iter, class T >
void sum_up(Iter first, Iter last, T& result)
{
	while(first != last){result += *first; ++first;}
}

template < class Container >
void sum_up_to_1st(Container& c)
{
	sum_up(c.begin() + 1, c.end(), c.front());
// 	while(first != last){result += *first; ++first;}
}

template < class inIter, class outIter >								
void normalize(const inIter& ifirst, const inIter& ilast, const outIter& ofirst, const outIter& olast)
{
	typename iterator_traits < inIter > :: value_type total = 0;
	typedef typename iterator_traits < outIter > :: value_type Out_T;
	inIter in; outIter out;
	
	for (in = ifirst; in != ilast; ++in) total += *in;
	
	double inv_total = 1.0/double(total);
	
	for (in = ifirst, out = ofirst;	 out!= olast && in != ilast;
			 ++out, ++in) 
		*out = Out_T(double(*in)*inv_total);
};


template < class I > struct numeric_type_traits {
  typedef I scalar;
};
	
template < class I > struct numeric_type_traits < Point <I> >
{
  typedef I scalar;
};

template < class I > struct numeric_type_traits < Pair <I> >
{
  typedef I scalar;
};


template < class T = double, 
		   class sqr_T = typename numeric_type_traits<T>::scalar > 
class M2_new
{
	typedef typename numeric_type_traits<T>::scalar Scalar;

	public:

		M2_new():sum(Scalar(0)),sum_sqr(Scalar(0)),count(Scalar(0)){}
		M2_new(const T& val):sum(val),sum_sqr(val*val),count(1){}

		M2_new& operator+=(const T& val)
		{
			sum		+=val;
			sum_sqr	+=val*val;
			++count;
            return *this;
		}
		
		void add_multiple(const T& val, const uint64_t& number)
		{
			sum 	+= val*T(number);
			sum_sqr += val*val*T(number);
			count 	+= number;
		}

		M2_new& operator+=(const M2_new& rhs)
		{
			sum			+=rhs.sum;
			sum_sqr		+=rhs.sum_sqr;
			count		+=rhs.count;
		}

		void reset() {
			sum		= 0;
			sum_sqr = 0;
			count	= 0;
		}

		uint64_t get0() const {return count;}
		double	 get1() const {return double(sum)	 /double(count);}
		double	 get2() const {return double(sum_sqr)/double(count);}

		double	getD2() const {double M1(get1()); return get2()-M1*M1;}

	private:
		T 						sum;
		sqr_T 					sum_sqr;
		uint64_t 				count;
};

template < class T , class ResT = double, class DoubleT = double >
class Displacement
{
	public:

	typedef T	 	input_type;
	typedef ResT 	output_type;
	typedef DoubleT	double_type;

	ResT operator () (const T& v1,const T& v2) const
	{
		T diff = v1 - v2;
		return ResT(diff*diff);
	}
};

template < class T , class ResT = T, class DoubleT = T >
class DisplacementComponents
{
	public:

	typedef T	 	input_type;
	typedef ResT 	output_type;
	typedef DoubleT	double_type;

	ResT operator () (const T& v1,const T& v2) const
	{
		T diff = v1 - v2;
		return ResT(diff.GetX()*diff.GetX(),diff.GetY()*diff.GetY(),diff.GetZ()*diff.GetZ());
	}
};

template < class T , class ResT = T, class DoubleT = T >
class Product
{
	public:
	typedef T	 		input_type;
	typedef ResT 		output_type;
	typedef DoubleT		double_type;

	ResT operator () (const T& v1,const T& v2) const
	{
// 		cout << "direct = " << v1*v2 << ", casted = " << ResT (v1*v2) << "___";
		return ResT (v1*v2);
	}
};


template < int step=1 >
struct LinearIncrementor{
  template < class T > void operator()(T& i){i=i+T(step);}
};

template < int factor=2, int divisor = 1 >
struct MultiplyIncrementor{
  template < class T > void operator()(T& i)
  {
	int old = i; i*=factor; i/=divisor;
	if(i == old) LinearIncrementor<>()(i);
  }
};

template < class Type, class Init > struct spec_init
{
  Type operator() (const Init& val){return Type(val);}
};

template < class T > struct spec_init < Point < T > , uint >
{
  Point < T > operator()(uint val){return Point < T >(val,val,val);}
  
};

template <class Func, class Incrementor = LinearIncrementor<1> >
class Correlator
{
  
	typedef typename Func::input_type 	T;
	typedef typename Func::output_type 	ResT;
	typedef typename Func::double_type 	DoubleT;

	const size_t window;
	const size_t offset;
	
	public:
	  
	Correlator(uint _win = -1, uint _offset = 0)
	:window(_win),offset(_offset)
	{
	  cout << "window = " << window << ", offset = " << offset << endl;
	}

	DoubleT operator [] (size_t j) const
	{
// 		return double(msd[j-1])/(double(msd.size()-(j)));
		
		DoubleT tmp(sum[j]); tmp *= 1.0/double(counter[j]);	
		cout << j << "  \t" << sum[j] << "  \t" << tmp << endl;
		return tmp;
	}

	size_t size () const
	{
		return sum.size();
	}

	void resetHistory()
	{
		history.clear();
	}

	void reset()
	{
		resetHistory();
		moments.reset();
		counter.clear();
		sum.clear();
	}
	

	
	

	void operator() (const T& val)
	{
		history.push_front(val);
		moments += val;
		
		if (history.size() > 2*window) history.resize(window);

		size_t max_diff	= min(window,history.size()-1);
		size_t sum_idx 	= 0;
		
		for (size_t diff = offset; diff <= max_diff; Incrementor()(diff), ++sum_idx )
		{
			if (sum_idx>=sum.size())
			{
			  sum	 .resize(max(2*sum_idx,size_t(64)), 0u );//spec_init<ResT,uint>()(0));
			  counter.resize(max(2*sum_idx,size_t(64)), 0u );
			}
			sum[sum_idx] += Func()(val,history[diff]); ++counter[sum_idx];
		}
	}

	Correlator& operator+=(const T& val)
	{
		this->operator()(val);
		return *this;
	}


	Correlator& operator += (const Correlator& rhs)
	{
		cout << "#sum of " << size() << " vs. " << rhs.size() << " elements\n";

		size_t max_size = max(size(),rhs.size());
		cout << "#max size = " << max_size<<endl;

		sum.resize(max_size,0u);counter.resize(max_size,0);
		for (uint i = 0; i < size(); ++i)
		{
			sum[i] 			+= rhs.sum[i];
			counter[i] 	    += rhs.counter[i];
		}
		moments += rhs.moments;
		return *this;
	}


	void Print(ostream& str, const double&  factor = 1.0)//, uint skip = 0 )
	{
		if (size())
		{
			size_t max_diff	= min(window,history.size()-1);
			size_t sum_idx 	= 0;
		
			for (size_t diff = offset; diff <= max_diff; Incrementor()(diff), ++sum_idx )
			{
				str << factor*double(diff) << ' '
	//					<< this->operator[](i)
						<< operator[](sum_idx) << ' '
	//					<< sum[i] << ' '
	//					<< counter[i]
						<< endl;
			}
		}
	}

	void PrintNormalized(ostream& str, const double&  factor = 1.0)//, uint skip = 0 )
	{
		double m1_sqr = moments.get1()*moments.get1();

		if (size()) 
		{  
		 	size_t max_diff	= min(window,history.size()-1);
			size_t sum_idx 	= 0;
		
			for (size_t diff = offset; diff <= max_diff; Incrementor()(diff), ++sum_idx )
			{ 
			  str << factor*double(diff) << ' '
					  << (operator[](sum_idx) - m1_sqr)/(moments.get2() - m1_sqr) << ' '
					  << endl;
			}
		}
	}

	const M2_new < T >& GetMoments() const {return moments;}

private:

		deque < T > 		history;
		vector < ResT > 	sum;
		vector < uint64_t > counter;
		M2_new < T >		moments;

};

template < class T > struct Approx
{
private: T compare;
public:
	Approx(const T& c):compare(c){}
	bool operator()(const T& val) const {return compare == val;}
};

template <> struct Approx < double >
{
private: double compare;
public:
	Approx(const double& c):compare(c){}
	bool inline operator()(const double& val)const
	{
		double d (compare - val);	return d < TOL && d >	-TOL;
	}
};

template < class T >
inline bool approx(const T& val1, const T& val2)
{
	double d ((val1 - val2)/val1);	return d < TOL && d >	-TOL;
// 	return Approx<T>(val1)(val2);
}

template < class T >
inline bool approx_0(const T& val)
{
	return (val < TOL && val > -TOL);
// 	return approx(val,T(0));
}

PD roots3(const double& a, const double& b, const double& c);
PD roots3(const double& a, const double& b, const double& c, const double& d);

template < class T >
PD roots3(const Point < T >& p )
{
	return roots3(double(p.X),double(p.Y),double(p.Z));
}

template < class M >
void c_poly_sym3(const M& m, PD& p)
{
	double a(m(0,0)),b(m(1,1)),c(m(2,2)),
				 d(m(0,1)),e(m(0,2)),f(m(1,2)),
				 d2(d*d),e2(e*e),f2(f*f);

	p.X = -(a+b+c);
	p.Y = (b+c)*a + b*c - d2 - e2 - f2;
	p.Z = a*f2+b*e2+c*d2-2.0*d*e*f-a*b*c;

}


template < class M , class VC >
bool eigenvalues3(const M& m, PD& e_vals, VC & e_vecs)
{
	PD cp; c_poly_sym3(m,cp);

	e_vals = roots3(cp);//.X,cp.Y,cp.Z);

	sort (e_vals.v,e_vals.v+3);

// 	cout << "cp = " << cp << ", e_vals = " << e_vals << endl;
//
	M m2;

	int success(1);

	m2 = m;

// 	double od[3] = {m(0,0),m(1,1),m(2,2)};

	m2(0,0) = m(0,0) - e_vals.X;
	m2(1,1) = m(1,1) - e_vals.X;
	m2(2,2) = m(2,2) - e_vals.X;

	success &= ( int )SpecialHom(m2,e_vecs[0]);

// 	cout << m2 << endl << "ev = " << ev << endl;
	m2 = m;
	m2(0,0) = m(0,0) - e_vals.Y;
	m2(1,1) = m(1,1) - e_vals.Y;
	m2(2,2) = m(2,2) - e_vals.Y;

	success &= ( int )SpecialHom(m2,e_vecs[1]);

// 	cout << m2 << endl << "ev = " << ev << endl;

	m2 = m;
	m2(0,0) = m(0,0) - e_vals.Z;
	m2(1,1) = m(1,1) - e_vals.Z;
	m2(2,2) = m(2,2) - e_vals.Z;

	success &= ( int )SpecialHom(m2,e_vecs[2]);

// 	cout << m2 << endl << "ev = " << ev << endl;


	return bool(success);

}

// template < class M, class V > bool solve_hom(M& m, V& result)
// {
// 	vector <
//
//
// }

template < class M > bool gauss(M& m)
{

	for (int p = 0; p != m.dim1()-1; ++p)
	{
		if (m(p,p))
		{
			for (int r = p+1; r != m.dim1(); ++r)
			{
				if (m(r,p))
				{
					double factor = m(r,p)/m(p,p);
					for (int c = p; c != m.dim2(); ++c)
					{
						m(r,c) -= factor*m(p,c);
					}
				}
			}
		}
		else return false;
	}

	for (int p = m.dim1()-1; p != 0; --p)
	{
		if (m(p,p))
		{
			for (int r = p-1; r != -1; --r)
			{
				if (m(r,p))
				{
					double factor = m(r,p)/m(p,p);
					for (int c = p; c != m.dim2(); ++c)
					{
						m(r,c) -= factor*m(p,c);
					}
				}
			}
		}
		else return false;
	}

	for (int p = 0; p != m.dim1(); ++p)
	{
		double factor = 1.0/m(p,p);
		m(p,p) = 1.0;
		for (int c = m.dim1(); c != m.dim2(); ++c)
		{
			m(p,c) *= factor;
		}
	}
	return true;
}

template < class M > void ReducedEchelon(M& m)
{
	int col = 0, r = 0;

// 	typedef typename M::value_type T;

	do
	{
		int i = r;
		while( approx_0(m(i,col)) )
		{
// 			m(i,col) = 0;
			++i;
			if (i == m.dim1())
			{
				++col; i = r;
				if (col == m.dim2()) goto STOP;
			}
		}
		// swap rows r and i and set leading coeff to 1
		double factor = 1.0/m(i,col);

		for (int j = col; j != m.dim2(); ++j)
		{
			swap(m(i,j),m(r,j));
			m(r,j) *= factor;
		}

// 		for (int i = 0; i != m.dim1()*m.dim2(); ++i)
// 		{
// 			if (approx_0(m[i])) m[i] = 0;
// 		}

		// set all other leading coeffs of the col to 0
		for (int i = 0; i != m.dim1(); ++i)
		{
			if (i != r)
			{
				factor = m(i,col);
				for (int j = 0; j != m.dim2(); ++j)
				{
					if (j != col)
					{
						m(i,j) -= factor*m(r,j);
						if (approx_0(m(i,j))) m(i,j) = 0;
					}
					else
						m(i,j) = 0;
				}
			}
		}
		++col; ++r;
// 		cout << m << endl;
	}
	while(col != m.dim2() && r != m.dim1());

	STOP:

// 	for (int i = 0; i != m.dim1()*m.dim2(); ++i)
// 	{
// 		if (approx_0(m[i])) m[i] = 0;
// 	}

	return;

}

template < class M >
bool SpecialHom( M& m, PD& p)
{
// 	cout << m << endl;
	ReducedEchelon(m);
// 	cout << m << endl;

	if (approx(1.0,m(0,0)) && approx(1.0,m(1,1)) && approx_0(m(2,2)))
	{
		p.Z = 1.0;
		p.Y = -m(1,2);
		p.X = -m(0,2);
		p.Normalize();
		return true;
	}
	else
	{
// 		cout << "matrix not in exspected form!\n";
// 		string dummy;
// 		getline(cin,dummy);
// 	if (!success)
// 		cout << m << endl;
		return false;
// 		cout << m << endl;
	}

}

template < class T > struct operation_traits
{
	typedef T 								square;
	typedef double						to_double;
};

template < class T > struct operation_traits < Point < T > >
{
	typedef T 								square;
	typedef Point < double > 	to_double;
};

template < class T > struct operation_traits < Pair < T > >
{
	typedef T 								square;
	typedef Pair < double > 	to_double;
};


/// maps _______________________________________________________________________

template < class T > class distance_map
{
private:
	const T center;
				T	distance;

public:
	typedef T 				value_type;
	typedef T& 				reference;
	typedef const T& 	const_reference;
	typedef T* 				pointer;
	typedef T					key_type;

// 	member_map	(M C::* _offset):offset(_offset){ }
	distance_map(const T& val):center(val){}

	const T& operator[] (const T& val) {
		distance  = val;
		distance -= center;
		return distance;
	}
};

template < class T > class squared_distance_map
{
public:

	typedef typename operation_traits < T >::square	value_type;
	typedef value_type& 														reference;
	typedef const value_type& 											const_reference;
	typedef value_type* 														pointer;
	typedef T																				key_type;

	squared_distance_map(const T& val):center(val){}

	value_type operator[] (const T& val) {
		T 			temp(val);
		temp -= center;
		return temp*temp;
	}

private:

	const 	T 	center;

};

template < class T > struct square_of
{
	typedef typename operation_traits < T >::square	value_type;
	typedef value_type& 														reference;
	typedef const value_type& 											const_reference;
	typedef value_type* 														pointer;
	typedef T																				key_type;
// // 	s(const T& val):center(val){}
	value_type operator[](const key_type& key) const {return key*key;}
	static value_type get(const key_type& key) {return key*key;}

private:
};


/// SUM ________________________________________________________________________

template < class Iter >
size_t
sum_counting(Iter first, Iter last, typename Iter::value_type& result)
{
	size_t count = 0;
	while(first != last) { result += *first; ++first; ++count;}
	return count;
}

template < class Iter >
void
sum(Iter first, Iter last, typename Iter::value_type& result)
{
	while(first != last) result += *first++;
}

template < class Iter >
typename Iter::value_type29182616
sum(Iter first, Iter last)
{
	typename Iter::value_type result;
	sum (first,last,result);
	return result;
}

/// MEAN _______________________________________________________________________

template < class Iter, class Result >
void
mean(Iter first, Iter last, Result& result)
{
	typename Iter::value_type temp(0);
	size_t count(sum_counting(first,last,temp));
	result 	= temp;
	double norm (1.0/count);
	result *=norm;
}

template < class Iter >
typename operation_traits < typename Iter::value_type > :: to_double
mean(Iter first, Iter last)
{
	typename operation_traits < typename Iter::value_type > :: to_double result;
	mean(first,last,result);
	return result;
}

/// MEAN_DYNAMIC _______________________________________________________________

template < class Iter, class Result >
size_t
mean_dynamic(Iter first, Iter last, Result& result)
{
	static double inverse_size(0.0);
	static size_t	last_count(0);

	typename Iter::value_type temp(0);

	size_t count(sum_counting(first,last,temp));

	if (last_count != count)
	{
		last_count 		= count;
		inverse_size 	= 1.0/count;
	}
	result  = temp;
	result *= inverse_size;

	return count;
}

template < class Iter >
typename operation_traits < typename Iter::value_type > :: to_double
mean_dynamic(Iter first, Iter last)
{
	typename operation_traits < typename Iter::value_type > :: to_double result;
	mean_dynamic(first,last,result);
	return result;
}

/// MEAN_KNOWN_SIZE_____________________________________________________________

template < class Iter, class Result >
void
mean_known_size(Iter first,Iter last,const double& inverse_size, Result& result)
{
	typename Iter::value_type temp(0);
	sum(first,last,temp);
	result  = temp;
	result *= inverse_size;
};

template < class Iter >
typename operation_traits < typename Iter::value_type > :: to_double
mean_known_size(Iter first, Iter last, const double& inverse_size )
{
	typename operation_traits < typename Iter::value_type > :: to_double result;
	mean_dynamic(first,last,inverse_size,result);
	return result;
};

template < class iter_T >
Point < double > CenterOfMass(const iter_T& begin, const iter_T& end)
{
	Point < double > result(0.0,0.0,0.0), temp;
	unsigned count(0);
	double inv_count;

	iter_T it(begin);
	do
	{
		temp = *it;
		result += temp;
		++it; ++count;
	}
	while(it != end);

	inv_count = 1.0/(double(count));

	result.X *= inv_count;
	result.Y *= inv_count;
	result.Z *= inv_count;

// 	cout << "center of mass = " << result << endl;

	return result;

}

template < class iter_T >
double RadiusOfGyration(const iter_T& begin, const iter_T& end, const DPoint& com)
{
	uint64_t cnt = 0;
	double result = 0.0;
	iter_T it(begin);
	while(it != end)
	{
		DPoint 	diff = DPoint(*it);
		diff 		-= com;
		result 		+= diff*diff;
		++it; ++cnt;
	}
	return sqrt(result/double(cnt));
}

template < class iter_T >
double SqrRadiusOfGyration(const iter_T& begin, const iter_T& end, const DPoint& com)
{
	uint64_t cnt = 0;
	double result = 0.0;
	iter_T it(begin);
	while(it != end)
	{
		DPoint 	diff = DPoint(*it);
		diff 		-= com;
		result 		+= diff*diff;
		++it; ++cnt;
	}
	return result/double(cnt);
}

template < class iter_T >
DPoint SqrRadiusOfGyrationComponents(const iter_T& begin, const iter_T& end, const DPoint& com)
{
	uint64_t cnt = 0;
	DPoint result;
	iter_T it(begin);
	while(it != end)
	{
		DPoint 	diff = DPoint(*it);
		diff 		-= com;
		result 		+= DPoint(diff.GetX()*diff.GetX(),
							  diff.GetY()*diff.GetY(),
							  diff.GetZ()*diff.GetZ());
		++it; ++cnt;
	}
	result /= double(cnt);
	return result;
}

template < class iter_T >
double RadiusOfGyration(const iter_T& begin, const iter_T& end)
{
	return RadiusOfGyration(begin,end,CenterOfMass(begin,end));
}

template < class Iter  >
double radius_of_gyration(Iter first, Iter last, const size_t& size, const double& inv_size)
{
	M2_new < typename Iter::value_type > moments(first,last,size,inv_size);
	return sqrt(moments.getD2());
}

template < class Iter >
double radius_of_gyration(Iter first, Iter last)
{
	M2_new < typename Iter::value_type > moments(first,last);
	return sqrt(moments.getD2());
}

template < class Iter >
double end_to_end_distance(Iter first, Iter last)
{
	typedef typename Iter::value_type iter_value;
	iter_value diff(*first);
	--last;
	diff -= *last;
	return sqrt(square_of<iter_value>::get(diff));
}


template < class iter_T >
double Rg(const iter_T& begin, const iter_T& end, const Point < double > &CoM)
{
	double 						result(0.0);
	unsigned 					count(0);
	Point < double >	diff;

	iter_T it(begin);

	do
	{
		diff  	= (*it)->GetPosition();
		diff 	 -=	CoM;
		result += diff*diff;
		++it,count++;
	}
	while(it != end);

	result /=(double(count));

// 	cout << "< sqr(Rg) > = " << result << endl;

	return result;
}

template < class iter_T >
double Rg(const iter_T& begin, const iter_T& end)
{
	Point < double > CoM = CenterOfMass(begin,end);
	return Rg(begin,end,CoM);
}


#endif


