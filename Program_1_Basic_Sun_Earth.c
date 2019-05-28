//  Program 1 Basic Sun Earth RK Method

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Runge Kutta function
double rk(double*,double*,double*,double*,double,double,double,double,double);

int main(void) {
    setvbuf(stdout,NULL,_IONBF,0); setvbuf(stderr,NULL,_IONBF,0);
    
    double xsun,ysun,xplanet, vxplanet, yplanet, vyplanet,r;    // planet & sun coordinates
    double t,dt=1000, orbit, day=86400, year=365.25*day;        // time & orbit settings
    double msun=1.98892e30, au=1.495978707e11, g=6.67408e-11;   // constants
    
    FILE * finout; //opening data file
    finout = fopen("data.txt","w");
    
    //setting initial conditions, planet is placed along the x axis with velocity pointing in y direction
    xsun=0, ysun=0;
    r=1*au, xplanet=r,yplanet=0, vxplanet=0, vyplanet=sqrt((g*msun/r));
    
    // loop to cycle through 1 year
    for (t=dt; t<year; t+=dt){
        
        // sending the variables to RK function
        orbit=rk(&xplanet,&vxplanet,&yplanet,&vyplanet,xsun,ysun,dt,g,msun);
       
        // printf ("%g\t%g\n",xplanet, yplanet);
        fprintf (finout,"%g\t%g\n",xplanet/au,yplanet/au);
    }
    
    // closing file
    fclose(finout);
    
    return 0;
    
}

double rk(double* xplanet, double* vxplanet, double* yplanet, double* vyplanet, double xsun, double ysun, double dt, double g, double msun){
    //  function calculates a values for both x and y, then b values etc. r changes each time to
    //  adjust for the change in step using the RK method.
    
    double r, ax, bx, cx, dx, ay, by, cy, dy;
    double avx, bvx, cvx, dvx, avy, bvy, cvy, dvy;
    
    r=sqrt(pow(xsun-*xplanet,2)+(pow(ysun-*yplanet,2)));
    ax=*vxplanet;               avx=(g*msun*(xsun-*xplanet)/(r*r*r));
    ay=*vyplanet;               avy=(g*msun)*(ysun-*yplanet)/(r*r*r);
    
    r=sqrt(pow(((xsun-*xplanet)+(dt/2)*ax),2)+pow(((ysun-*yplanet)+(dt/2)*ax),2));
    bx=*vxplanet+(dt/2)*avx;    bvx=((g*msun*(xsun-(*xplanet+(dt/2)*ax)))/(r*r*r));
    by=*vyplanet+(dt/2)*avy;    bvy=((g*msun*(ysun-(*yplanet+(dt/2)*ay)))/(r*r*r));
    
    r=sqrt(pow(((xsun-*xplanet)+(dt/2)*bx),2)+pow(((ysun-*yplanet)+(dt/2)*bx),2));
    cx=*vxplanet+(dt/2)*bvx;    cvx=((g*msun*(xsun-(*xplanet+(dt/2)*bx)))/(r*r*r));
    cy=*vyplanet+(dt/2)*bvy;    cvy=((g*msun*(ysun-(*yplanet+(dt/2)*by)))/(r*r*r));
    
    r=sqrt(pow(((xsun-*xplanet)+(dt/2)*cx),2)+pow(((ysun-*yplanet)+(dt/2)*cx),2));
    dx=*vxplanet+(dt*cvx);      dvx=((g*msun*(xsun-(*xplanet+(dt*cx))))/(r*r*r));
    dy=*vyplanet+(dt*cvy);      dvy=((g*msun*(ysun-(*yplanet+(dt*cy))))/(r*r*r));


    *vxplanet = *vxplanet+(dt/6)*(avx+(2*bvx)+(2*cvx)+dvx);
    *xplanet = *xplanet+(dt/6)*(ax+(2*bx)+(2*cx)+dx);
    
    *vyplanet = *vyplanet+(dt/6)*(avy+(2*bvy)+(2*cvy)+dvy);
    *yplanet = *yplanet+(dt/6)*(ay+(2*by)+(2*cy)+dy);
    
    
    return 0;
}
