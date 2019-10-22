/**************************************************************/
/*      Finite difference code for knotted Hopf Skyrmions     */
/*      created April 2019: Jack Binysh, Gareth Alexander     */
/**************************************************************/

#include "SolidAngle.h"    
#include "Geometry.h"
#include "InputOutput.h"
#include "Constants.h"    
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
//#include <complex>
#include <string>
#include <sstream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// ============================================================

const int Nmax = 5000;       // Number of timesteps
const int stepskip = 5000;    // print director file every stepskip timesteps

// ============================================================

/* functions */
void initialise(void);
void startconfig(void);
void update(void);
void writeVTKFiles(void);
// modified numerical recipes routine
#define n 3     
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot);
#undef n

/* global variables */
int n,LL;
double koverk2,q0,Gamma,dt;
double *nx,*ny,*nz;
double *hx,*hy,*hz;

#define BC 1 // Dirichlet boundary conditions on (1) or off (0)

//int main(int argc, char** argv) 
int main (int argc, char*argv[])
{
  float doverp=(float) atof(args[1]);
  koverk2 = atof(args[2]); 
  Nx=(int) atoi(args[3]);
  Ny=(int) atoi(args[4]);
  Nz=(int) atoi(args[5]);
  output_directory = str(args[6]);
  knot_filename= str(args[7]);
  
  q0=(2.f*M_PI/(float)(Nz))*doverp;
  
  Gamma = 0.65;     // relaxation constant
  dt = 1.0;        // integration timestep


  // knot stuff
  if(argc<2)
  {
    std::cerr << "Please specify a folder name - see README" << endl;
    std::cerr << "e.g. SolidAngle Unknot if your executable is called SolidAngle and there's a folder within the knots folder called Unknot" << endl;
    return 1;
  }

  knot_filename = argv[1];
  // end knot stuff

  int step=0;

  InitialiseSystemParameters();
  LL=Nx*Ny*Nz;
  initialise();
  startconfig();

  cout << "starting simulation" << endl;

  for (n=0; n<=Nmax; n++)
  {
    if (step==stepskip)
    {
      writeVTKFiles();  // output VTK files for use with ParaView 
      step=0;
    }
    step++;
    update();
  }
} // end main


/**********************************************************************/
void initialise(void)
{
  nx=new double[LL];
  ny=new double[LL];
  nz=new double[LL];

  hx=new double[LL];
  hy=new double[LL];
  hz=new double[LL];
}

/**********************************************************************/
void startconfig(void)
{
  int j,k,l,m;
  double omega,alpha,R,rho;
  // set initialisation type

  // initialise the knot
  Link Curve;
  if(0==InitialiseFromFile(Curve))
  {
    cout << "Filling in the Geometry of the Input Curve \n";
    ComputeGeometry(Curve);
    OutputScaledKnot(Curve); 
  }

  k=l=m=0;  
  omega=0.0; // overly cautious!!

  int c; // need to adapt to number of components
  c=0; // no real reason for this    
  R = TubeRadius; // radius for the tubular neighbourhood of the curve 
  alpha=0.0; // overly cautious!!
  viewpoint Point;

  double LL=Nx*Ny*Nz;
  for (j=0; j<LL; j++) 
  {    
    double tempnx = 0.0;
    double tempny = 0.0;
    double tempnz = 1.0;

    Point.xcoord = x(k); //1.0*k-Nx/2.0+0.5;
    Point.ycoord = y(l); //1.0*l-Ny/2.0+0.5;
    Point.zcoord = z(m); //1.0*m-Nz/2.0+0.5;

    // calculate the distance to the curve 
    rho = ComputeDistanceOnePoint(Curve,Point);
    // put a Skyrmion texture inside the tubular neighbourhood
    if (rho < R)
    {
      // here is the solid angle
      omega = ComputeSolidAngleOnePoint(Curve,Point);
      // inelegant method for getting longitudinal phase and component
      alpha = ComputeLongitudinalPhase(Curve,Point);
      c = WhichComponent(Curve,Point);  
      // set the director  
      // Chirality +1 is right handed.
      // Chirality -1 is left handed
      tempnx = sin(M_PI*rho/R)*cos(0.5*omega-degrees[c]*alpha);
      tempny = sin(M_PI*rho/R)*sin(0.5*omega-degrees[c]*alpha);
      tempnz = -cos(M_PI*rho/R);
    }

#if BC // Dirichlet boundary conditions 
    if (m==0) 
    {
      tempnx = 0.0;
      tempny = 0.0;
      tempnz = 1.0;
    }
    if (m==Nz-1)
    { 
      tempnx = 0.0;
      tempny = 0.0;
      tempnz = 1.0;
    }
#endif
    nx[j]=tempnx;
    ny[j]=tempny;
    nz[j]=tempnz;
    k++;
    if (k==Nx) {l++; k=0;}
    if (l==Ny) {m++; l=0;}
  }

} // end startconfig

/**********************************************************************/
void update(void)
{
  int j,k,l,m,jj;
  int xup,xdwn,yup,ydwn,zup,zdwn;
  double Dxnx,Dynx,Dznx,Dxxnx,Dyynx,Dzznx;
  double Dxny,Dyny,Dzny,Dxxny,Dyyny,Dzzny;
  double Dxnz,Dynz,Dznz,Dxxnz,Dyynz,Dzznz;
  double hdotn,sqrtndotn;

  k=l=m=0;

  /* Calculate derivatives, molecular field and do update */
  for (j=0; j<LL; j++) {
    // define neighbouring nodes
    xup=j+1; xdwn=j-1; yup=j+Nx; ydwn=j-Nx; zup=j+Nx*Ny; zdwn=j-Nx*Ny;
    // correct for periodic boundaries
    if (k==0) {xdwn+=Nx;}
    if (k==Nx-1) {xup-=Nx;}
    if (l==0) {ydwn+=Nx*Ny;}
    if (l==Ny-1) {yup-=Nx*Ny;}
    if (m==0) {zdwn+=LL;}
    if (m==Nz-1) {zup-=LL;}

#if BC // cheap fix for the Dirichlet boundary conditions
    if (m==0) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
    if (m==Nz-1) {xup=j; xdwn=j; yup=j; ydwn=j; zup=j; zdwn=j;}
#endif

    // calculate first order derivatives
    Dxnx = (nx[xup]-nx[xdwn])/2.0;
    Dynx = (nx[yup]-nx[ydwn])/2.0;
    Dznx = (nx[zup]-nx[zdwn])/2.0;

    Dxny = (ny[xup]-ny[xdwn])/2.0;
    Dyny = (ny[yup]-ny[ydwn])/2.0;
    Dzny = (ny[zup]-ny[zdwn])/2.0;

    Dxnz = (nz[xup]-nz[xdwn])/2.0;
    Dynz = (nz[yup]-nz[ydwn])/2.0;
    Dznz = (nz[zup]-nz[zdwn])/2.0;

    // calculate second order derivatives
    Dxxnx = nx[xup]-2.0*nx[j]+nx[xdwn];
    Dyynx = nx[yup]-2.0*nx[j]+nx[ydwn];
    Dzznx = nx[zup]-2.0*nx[j]+nx[zdwn];

    Dxxny = ny[xup]-2.0*ny[j]+ny[xdwn];
    Dyyny = ny[yup]-2.0*ny[j]+ny[ydwn];
    Dzzny = ny[zup]-2.0*ny[j]+ny[zdwn];

    Dxxnz = nz[xup]-2.0*nz[j]+nz[xdwn];
    Dyynz = nz[yup]-2.0*nz[j]+nz[ydwn];
    Dzznz = nz[zup]-2.0*nz[j]+nz[zdwn];

    // calculate molecular field
    hx[j] = K*(Dxxnx+Dyynx+Dzznx) - 2.0*K*q0*(Dynz-Dzny);
    hy[j] = K*(Dxxny+Dyyny+Dzzny) - 2.0*K*q0*(Dznx-Dxnz);
    hz[j] = K*(Dxxnz+Dyynz+Dzznz) - 2.0*K*q0*(Dxny-Dynx);

    hdotn = nx[j]*hx[j] + ny[j]*hy[j] + nz[j]*hz[j];
    hx[j] -= nx[j]*hdotn;
    hy[j] -= ny[j]*hdotn;
    hz[j] -= nz[j]*hdotn;

#if BC // Dirichlet boundary conditions along z
    if (m==0) {hx[j] = 0.0; hy[j] = 0.0; hz[j] = 0.0;}
    if (m==Nz-1) {hx[j] = 0.0; hy[j] = 0.0; hz[j] = 0.0;}
#endif

    // keep track of boundaries -- fix me!!
    k++;
    if (k==Nx) {l++; k=0;}
    if (l==Ny) {m++; l=0;}
  }

  // do the update -- first order Euler
  for (j=0; j<LL; j++) {
    // director
    nx[j] += Gamma*hx[j]*dt;
    ny[j] += Gamma*hy[j]*dt;
    nz[j] += Gamma*hz[j]*dt;
    //normalise
    sqrtndotn = sqrt(nx[j]*nx[j] + ny[j]*ny[j] + nz[j]*nz[j]);
    nx[j] /= sqrtndotn;
    ny[j] /= sqrtndotn;
    nz[j] /= sqrtndotn;
  }

} // end update

/**********************************************************************/
void writeVTKFiles(void)
{
  int j;

  string fn = output_dir + "/vtk_director.vtk";
  ofstream Aout (fn.c_str());
  Aout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
  Aout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
  Aout << "ORIGIN " << x(0) << ' ' << y(0) << ' ' << z(0) << '\n';
  Aout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
  Aout << "POINT_DATA " << Nx*Ny*Nz << '\n';
  Aout << "VECTORS n double\n";

  for (j=0; j<LL; j++) {
    Aout << nx[j] << " " << ny[j] << " " << nz[j] << endl;
  }

  Aout.close();
} // end writeVTKFiles


/**************************************************************************/
/*  The following routines are based on those given in Numerical Recipes  */
/*   and are here copied from Alex and Davide's lattice Boltzmann code    */
/**************************************************************************/

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
                            a[k][l]=h+s*(g-h*tau);
#define n 3
void jacobi(double (*a)[n], double d[], double (*v)[n], int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c;
  double b[n],z[n];

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip< n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
            && (fabs(d[iq])+g) == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((fabs(h)+g) == fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0;j<=ip-1;j++) {
            ROTATE(a,j,ip,j,iq)
          }
          for (j=ip+1;j<=iq-1;j++) {
            ROTATE(a,ip,j,j,iq)
          }
          for (j=iq+1;j<n;j++) {
            ROTATE(a,ip,j,iq,j)
          }
          for (j=0;j<n;j++) {
            ROTATE(v,j,ip,j,iq)
          }
          ++(*nrot);
        }
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  cout << "Too many iterations in routine jacobi" << endl;
  exit(0);
}
#undef n
#undef ROTATE


