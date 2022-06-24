#include "numerics.h"


PD roots3(const double& a, const double& b, const double& c)
{
	static const double

					_1o3		= 1.0 /  3.0,
					_1o27		= 1.0 / 27.0,
					_2o27		= 2.0 / 27.0,

					angle_1 = 2.0*3.141592653589793238462643383*_1o3,
					angle_2 = 2.0*angle_1;

	double 	p 			= b - _1o3*a*a,
				 	q 			= _2o27*a*a*a -_1o3*a*b + c,

				 	po3 		= _1o3*p,
				 	qo2 		=  0.5*q,

				 	Dp 			= po3*po3*po3,
				 	Dq 			= qo2*qo2,

				 	D  			= Dp + Dq;

	if (D < 0.0)
	{

		double 	r 				= 	sqrt( - Dp ),
						phi 			= 	_1o3*acos( - 0.5 * q / r ),
						prefactor = 	2.0*sqrt( - po3 );

		PD result ( prefactor * cos ( phi ) ,
								prefactor * cos ( phi + angle_1 ),
								prefactor * cos ( phi + angle_2 ));

		result.X -= _1o3*a;
		result.Y -= _1o3*a;
		result.Z -= _1o3*a;

		return result;
	}
	else
	{
		cout << " D >= 0 \n";
		return PD(0.0);
	}

};

PD roots3(const double& a, const double& b, const double& c, const double& d)
{
	double f = 1.0/a;
	return roots3(b*f,c*f,d*f);
}

void TensorOrder::solve() { 
	vector < DPoint > e_vecs(3);
	DPoint e_vals;

	eigenvalues3(*this,e_vals,e_vecs);

	uint idx_max = 0; double max = e_vals[0];

	if (e_vals[1] > max){ idx_max = 1; max = e_vals[1];}
	if (e_vals[2] > max){ idx_max = 2; }

	director = e_vecs[idx_max];
	director.SetW(e_vals[idx_max]);

  }
void TensorOrder::Print() const {
cout << sum << endl;
cout << counter << endl;
 }
const DPoint& TensorOrder::getDirector() const {
	const_cast < TensorOrder& > (*this).finish();
	return director;
  }
void TensorOrder::reset() {
M3D::operator=(0);
sum 		 = 0;
fresh 		 = false;
counter 	 = 0;
 }
TensorOrder& TensorOrder::operator+=(const TensorOrder::input_type& val) {
  if ( ! ( val[0] != val[0] || val[1] != val[1] || val[2] != val[2] ))
  {
	fresh = true;
	double sqr_norm_factor = 1.0/(val*val);
	// cout << "adding " << val << " to TensorOrder\n";
	for (uint a = 0; a != 3; ++a)
	{
	  for (uint b = a; b != 3; ++b)
		sum(a,b) += val[a]*val[b]*sqr_norm_factor;
	}
	++counter;
  }
  return *this;
}

TensorOrder::TensorOrder() :sum(0),counter(0),fresh(false){M3D::operator=(0);}
const double& TensorOrder::getOrderParameter() const {return this->getDirector().GetW();}
void TensorOrder::finish() {
	if (fresh)
	{
	  fresh = false;
	  
	  double norm_fact = 1.0/counter;
	  
	  // diagonals:
	  M3D::operator()(0,0)=1.5*sum(0,0)*norm_fact-0.5;
	  M3D::operator()(1,1)=1.5*sum(1,1)*norm_fact-0.5;
	  M3D::operator()(2,2)=1.5*sum(2,2)*norm_fact-0.5;
	  
	  // symmetric non-diagonals:
	  M3D::operator()(0,1)=1.5*sum(0,1)*norm_fact;
	  M3D::operator()(1,0)=1.5*sum(0,1)*norm_fact;
	  
	  M3D::operator()(0,2)=1.5*sum(0,2)*norm_fact;
	  M3D::operator()(2,0)=1.5*sum(0,2)*norm_fact;
	  
	  M3D::operator()(1,2)=1.5*sum(1,2)*norm_fact;
	  M3D::operator()(2,1)=1.5*sum(1,2)*norm_fact;
	  solve();
	}
  }
