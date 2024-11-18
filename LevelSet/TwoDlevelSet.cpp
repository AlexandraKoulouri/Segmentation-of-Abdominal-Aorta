///////////////////////////////////////////////////////////////////////////////////////////////
/*Level set implementation according to  Chunming Li, Chenyang Xu, Changfeng Gui, and  Martin D. Fox, 
//“Level Set Evolution Without Re-initialization: A New Variational Formulation”, 
//IEEE International Conference on Computer Vision and Pattern Recognition (CVPR)
//vol. 1, pp. 430-436, San Diego, 2005*/ 
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <mex.h>
#include <mat.h>
#include <matrix.h>

#define cimg_plugin "cimgmatlab.h"

#include "CImg.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace cimg_library;
using namespace std;

//globa values
const double  epsilon=1.5; // the papamater smooth Dirac function (default value 1.5);
const float precision=0.1f;//precision of the error estimation

//functions
CImg<double> DiracU( CImg<double>& u0) ;
CImg<double> Heaviside(CImg<double>& u0);
CImg<float> GradientVectorFlow(const CImg<>& Img);

//-----------------
// Main-MexFunction
//-----------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   if (nrhs < 8) mexErrMsgTxt("No enough input arguments.");
   if (nrhs >8) mexErrMsgTxt("Too many input arguments.");
   if (nrhs == 8){
         
       //Input Parameters (7 inputs)
       CImg<double>  u(prhs[0],true);    //Initial level set function u(distance function)
       CImg<double>  g(prhs[1],true);    //image
       CImg<double>  Img(prhs[2],true);  //g->edge indicator function
      
       
       const double dt= mxGetScalar( prhs[3]);    //time-step of iteration
       const double lambda = mxGetScalar(prhs[4]);//coefficient of the weighted length term L(u)
       const double mu = mxGetScalar(prhs[5]);    // coefficient of the internal (penalizing) energy term P(u)
       const double v = mxGetScalar(prhs[6]);     //coefficient of the weighted area term A(u)

       const int nb_iter=mxGetScalar(prhs[7]);
       //////////////////////////////////////////////////////////////////////////////////////////////////
       
       //Initialization
       CImg<double> dg(g.dimx(),g.dimy(),2);    //derivatives of g function
      // CImg<double> Laplac_u(u.dimx(),u.dimy());//Laplacian operator for level set function u
       //CImg<double> K(u.dimx(),u.dimy());    //Curvature of u
       CImg<double> du(g.dimx(),g.dimy(),2); //derivatives of u
       CImg<double> veloc(u.dimx(),u.dimy());
        
       //Normalize gradient of level set function u.
       CImg<double> N_u(u.dimx(),u.dimy(),2); 
       
       float E=1e20f;//Initial error
       //const unsigned int nb_iter = 400;//Number of iterations  
       
        //////////////////////////////////////////////////////////////////////////////////////////////

        //GVF=>gradient Flow estimation
       // CImg<double> flow = GradientVectorFlow(Img);
       
       //Estimate the gradient of g
       // Compute 2D gradient (spatial derivatives).
       cimg_for3XY(g,x,y){
               dg(x,y,0)=0.5*(g(_n1x,y)-g(_p1x,y)), 
               dg(x,y,1)=0.5*(g(x,_n1y)-g(x,_p1y));
        }
                
        for (unsigned int iter=0; iter<= nb_iter; iter++) {
                
              const float Eold = E;
              E = 0;               
              cimg_for3XY(u,x,y) {
                 const double
                      ux=0.5*(u(_n1x,y)-u(_p1x,y)),
                      uy=0.5*(u(x,_n1y)-u(x,_p1y));  //derivatives of u

                  const double Mag_du= sqrt(pow(ux,2)+ pow(uy,2)+1e-10); //magnitude of grad(u)
                      N_u(x,y,0)=ux/Mag_du;
                      N_u(x,y,1)=uy/Mag_du;
              }
                                   
             
             CImg<double> diracF=DiracU(u);
             CImg<double>HeavisideF=Heaviside(u);
           
              
             cimg_for3XY( N_u,x,y) {
                     
                  const double
                     dN_ux=0.5*(N_u(_n1x,y,0)-N_u(_p1x,y,0)),
                     dN_uy=0.5*(N_u(x,_n1y,1)-N_u(x,_p1y,1)),
                     Laplac_u=(u(_n1x,y) + u(_p1x,y) + u(x,_n1y) + u(x,_p1y))-4*u(x,y),
                     K=dN_ux+dN_uy;
                  const double
                     uxx=(u(_n1x,y)+u(_p1x,y)-2*u(x,y)),
                     uyy=(u(x,_n1y)+u(x,_p1y)-2*u(x,y)),
                     ux =0.5*(u(_n1x,y)-u(_p1x,y)),
                     uy =0.5*(u(x,_n1y)-u(x,_p1y)),//derivatives of u
                     Mag_du = sqrt(pow(ux,2)+ pow(uy,2)+1e-10); //magnitude of grad(u)
                     
                  veloc(x,y)=lambda* diracF(x,y)*( dg(x,y,0)* N_u(x,y,0) +dg(x,y,1)* N_u(x,y,1) + g(x,y)*K)+mu*(Laplac_u-K)+v*g(x,y)*diracF(x,y);//-0.08*diracF(x,y)*(flow(x,y,0)*uxx+flow(x,y,1)*uyy);//+0.025* vel(x,y) *diracF(x,y)); 
                    
                  E+=lambda*g(x,y)*diracF(x,y)* Mag_du+1/2*pow( Mag_du-1,2)+HeavisideF(x,y)*g(x,y);
             }
             double m, M=veloc.maxmin(m);
             const double xdt=dt/max(abs(m),abs(M));
             u+=dt*veloc;
             if (abs(Eold-E)<0.1f) break;                            
    }
 
        
   
  
     
 /* for (unsigned int iter=0; iter<160; iter++) {
        const float Eold = E;
              E = 0;               
           
        CImg<double> diracF=DiracU(u);
        CImg<double>HeavisideF=Heaviside(u);
        cimg_for3XY(u,x,y)  {
            const double
                ux = 0.5*(u(_n1x,y)-u(_p1x,y)),
                uy = 0.5*(u(x,_n1y)-u(x,_p1y)),
                uxx =( u(_n1x,y)+u(_p1x,y)-2*u(x,y)),
                uyy = (u(x,_n1y)+u(x,_p1y)-2*u(x,y)),
                uxy = 0.25*(u(_p1x,_p1y)+u(_n1x,_n1y)-u(_n1x,_p1y)-u(_p1x,_n1y)),
                ngrad =ux*ux+uy*uy+1e-10, 
                kk=(uxx*uy*uy-2*uy*ux*uxy+uyy*ux*ux)/pow(sqrt(ngrad),3),
                Laplace_u1=uxx+uyy,
                Nux=ux/sqrt(ngrad),
                Nuy=uy/sqrt(ngrad);
           // term1(x,y)=lambda*diracF(x,y)*( dg(x,y,0)* Nux+ dg(x,y,1)* Nuy+g(x,y)*k);
             
            veloc(x,y)=4*g(x,y)*diracF(x,y)+lambda*diracF(x,y)*(g(x,y)* kk+dg(x,y,0)*  Nux+dg(x,y,1)* Nuy)+mu*(Laplace_u1-kk);    
            E+=lambda*g(x,y)*diracF(x,y)* sqrt(ngrad)+1/2*pow( sqrt(ngrad)-1,2)+HeavisideF(x,y)*g(x,y);
        }
        double m, M = veloc.maxmin(m);
        const double xdt = dt/max(abs(m),abs(M));
        u+=xdt*veloc;
        if (abs(Eold-E)<0.01f) break;      
        
   }*/
         
   plhs[0]= u.toMatlab();
                   
  } 
      
   return;    
}    
  


   
CImg<double> DiracU(CImg<double>& u0) {

  CImg<double> u(u0.dimx(),u0.dimy());
  u.fill(0);
   cimg_forXY(u0,x,y) { 
                  
       if (u0(x,y)<=epsilon && u0(x,y)>=-epsilon){
       u(x,y)=(double)1/(2*epsilon)*(1+cos(3.14*u0(x,y)/epsilon));
       }
  
   }
    return u;
}
     
CImg<double> Heaviside(CImg<double>& u0) {

  CImg<double> u(u0.dimx(),u0.dimy());
  u.fill(0);
   cimg_forXY(u0,x,y) { 
                  
       if (u0(x,y)<=epsilon && u0(x,y)>=-epsilon){
       u(x,y)=(double) 1/2-u0(x,y)/(2*epsilon)+1/(2*3.14)*sin(-3.14*u0(x,y)/epsilon);
       }
       if (u0(x,y)<-epsilon) u(x,y)=1;
  
   }
    return u;
}
       

CImg<float> GradientVectorFlow(const CImg<>& Img){

    //const CImg<> Img  = src.get_pointwise_norm(1).resize(src.dimx(),src.dimy()).normalize(0,1);
    CImg<float> u(Img.dimx(),Img.dimy(),1,2),temp(Img.dimx(),Img.dimy(),1,2);
    u.fill(0);temp.fill(0);
    float  K=0.08f;
     float dt=0.20f;
     float E=1e20f;
  
    CImg<float> blurredImg = Img.get_blur(1);//.normalize(0,1);
    //CImg<float> blurredImg = ApplyGaussian(Img,sigma);
    CImg<float> MagGrad(Img),f(Img), fx(Img), fy(Img),p(Img);
    fx.fill(0);fy.fill(0);
    //CImg<float>pp(Img);
    //define the EDGE MAP f(x,y)
    cimg_for3XY(f,x,y){
           const float 
                Ix = 0.5f*(blurredImg(_n1x,y)-blurredImg(_p1x,y)),
                Iy = 0.5f*(blurredImg(x,_n1y)-blurredImg(x,_p1y));
           f(x,y) = sqrt(Ix*Ix+Iy*Iy);
     }
     f.normalize(0,1);
     
     //Estimate the fx,fy derivatives
       cimg_for3XY(f,x,y){
               
          fx(x,y)= 0.5f*(f(_n1x,y)-f(_p1x,y));
          fy(x,y)= 0.5f*(f(x,_n1y)-f(x,_p1y));
          MagGrad(x,y) = sqrt(pow(fx(x,y),2)+pow(fy(x,y),2)); 
          p(x,y)=exp((float) (-MagGrad(x,y)/K));
       }   
       
       float m, M = p.maxmin(m);
      const double xdt = dt/max(abs(m),abs(M));
   
      for (unsigned int iter=0; iter<100; iter++) {
      
          const float Eold = E;
          E = 0;    
          cimg_for3XY(u,x,y){
             float tmf=0;
              const float
                //p=exp((float) (-MagGrad(x,y)/K)),
                q=1-p(x,y);
               temp(x,y,0) = p(x,y)*(u(_n1x,y,0)+u(x,_n1y,0)+u(_p1x,y,0)+u(x,_p1y,0)-4*u(x,y,0))-q*(u(x,y,0)-fx(x,y));
               temp(x,y,1) = p(x,y)*(u(_n1x,y,1)+u(x,_n1y,1)+u(_p1x,y,1)+u(x,_p1y,1)-4*u(x,y,1))-q*(u(x,y,1)-fy(x,y));
              const float
                 ux  = 0.5f*(u(_n1x,y,0)-u(_p1x,y,0)),
                 uy  = 0.5f*(u(x,_n1y,0)-u(x,_p1y,0)),
                 vx  = 0.5f*(u(_n1x,y,1)-u(_p1x,y,1)),
                 vy  = 0.5f*(u(x,_n1y,1)-u(x,_p1y,1));
                 
               E+=p(x,y)*(ux*ux+uy*uy+vx*vx+vy*vy)+q*(sqrt(pow(u(x,y,0)-fx(x,y),2)+pow(u(x,y,1)-fy(x,y),2)));
          }
        
        u+=xdt*temp;
        if (abs(Eold-E)<0.9f) break;
        //std::fprintf(stderr,"\dt - Iteration %d ",xdt); std::fflush(stderr);
     }


    return u;
}




