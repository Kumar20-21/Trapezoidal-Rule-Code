

#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include <iomanip> 
#include <string>


using std:: vector ;
using namespace std;

/*
 ::  Open the terminal and first compile the programe using the following 
 :: 
 ::  g++ TrapezoidalRule.cc   
 :: Then Run it by typing
 :: ./a.out

 */

/*
 *  Function values $f(x_{j})$ at the equidistant quadrature 
 *  points will be computed under this function. This is just declaration 
 *  of the function and it is defined in the below. You can give any function 
 * whatever you want. Currently it is taken for exponential in the 
 * interval  [a,b].
 */
void get_data(vector<double>& , double , double , int );

/*
 *  Function Declaration for the trapezoidal rule. This function will
 *  will return the integral of function $f$ in the interval [a,b] using
 *  composite trapezoidal rule
 *
 */
double Integral_by_TtapezoidalRule( const vector<double>& , double , double, int  );


int main(){

  /*
   * a:: start point of the integration interval 
   */
  int a=1.0;
  /*
   *  b:: end point of the integration interval 
   */
  int b =2.0;  // end point of the interval

  /*
   *  maximum_grid_size : this is to compute order of convergence
   */


  int maximum_grid_size =512;
  
  /*
   *  maxError :: error incurred in the numerical intefgration
                  this is obtained by comparing the true integral
                  with approximate integral
   */
  double maxError(0.0);

  /*
   *  relError_at_N  :: relative error in the approximation
   */
  double relError_at_N(0.0);

  /*
   *  exact_integ :: exact value of integral, currently written for 
   *                 exponetial, but if u change f, for comparision
   *                 should change this
   */
  double exact_integ = exp(b)-exp(a);

  /*
   *  Following piece of code is just written to produce results in a nice 
   *  table format
   */
  std::cout<<"\n \n \n";
  
  std::cout<<"==========================Error Report for Trapzoidal Rule====================="<<std::endl;

  std::cout<<"\n \n ";

  std::cout<<"------------------------------------------------------------------------------"<<std::endl;
  std::cout<<" max Error   "<< "|     relative Error "<<" |     Order of Convergence"<< "   | N    "<<std::endl;
 std::cout<<"\n ";
  
 std::cout<<"------------------------------------------------------------------------------"<<std::endl;
  /*
   *  To varify the accuracy as well as the order of convergence, we compute the   *  integral of different values of N. For illustartion, say, start with N=4
   *  then keep this error, and compute again for N=8, take its log of ratios      *  this is eactly the order of convergence. Do it for several N and verify 
   *  the order of convergence of the method
   */  
 
  for(int ref_factor=1; ref_factor<=maximum_grid_size ;ref_factor*=2){

    /*
     *  N :: number of quadrature point
     */
    int N=4*ref_factor; ;   

    /*
     * f(N+1, 0.0) : vector of size N+1. this will store function values
     *               at the quadrature points 
     */

  vector<double> f(N+1,0.0); 
  /*
   *   get_data(f,a,b, N); by calling this function, function subject to the int   *   egration will be computed at the quadrature points $x_j$, and stored in f
   */
  get_data(f,a,b, N);  // get data at quadrature point

  /*
   * approx_integ: approximate integral obtained by trapezoidal rule 
   *
   *  Integral_by_TtapezoidalRule(f,a,b,N): this function will compute 
   *  the integral over interval [a,b]. As a input you need to give,
   *  f: function values at the grid point
   *  a: starting point of the integration interval,
   *  b: end point of the integration interval
   *  N: Number of quadrature point
   *  
   */
 
  double approx_integ =  Integral_by_TtapezoidalRule(f,a,b,N);

  /*
   * error : error between true integral and approximate integral
   */
    
  double error =  exact_integ-approx_integ;
  // std::cout<<"error ="<<error<<std::endl;
 
  /*
   *  This piece of code is to compute the order of convergence
   */
   double relError_at_2N =error/fabs(exact_integ);
   double order_of_convergence(0.0) ;
 
   if(ref_factor!=1){
    order_of_convergence =log(relError_at_N/relError_at_2N)/(log(2.0));  
   }

   /*
   *   Again, following piece of code is just written to 
   *   produce results in a nice 
   *   table format
   */
  
   std::cout<<std::scientific;
   std::cout<<setprecision(2)<<std::endl;
   std::cout<< error <<"    |"<<"      "<<relError_at_2N <<"     "<<" |"<<"       " <<order_of_convergence<<"             | "<< N <<std::endl;
     relError_at_N=relError_at_2N;  
  }
  std::cout<<"============================================================================="<<std::endl;
  
}

  /*
   *   get_data(f,a,b, N); by calling this function, function subject 
   *   to the integration will be computed at the quadrature points $x_j$, 
   *    and stored in f
   */
void get_data(vector<double>& f,  double a, double b, int N)
{

  double h =(b-a)/(double)N ;

  for(int jj=0;jj<=N; jj++){
    double x =a+jj*h ;
    f[jj]=exp(x);
  }
}

double Integral_by_TtapezoidalRule(const vector<double>& f, double a, double b, int N)
{
  /*
   * Will give the integral of $f$ in [a,b] using trapezoidal rule
   */
 double h =(b-a)/(double )N;

  double integ(0.0);

  //integ =0.5*(f[0]+f[N]);

  for(int jj=0; jj<=N; jj++){
    if(jj==0 || jj==N){
      integ +=(0.5*f[jj]);
    }
    else{  // here we have mistake during live session
    integ += f[jj];
    }
  }
  integ *= h;
  return (integ);
}

