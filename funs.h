void PN1(int n,double c,double m,double coes[2],double r[6])
{
    double x=r[0],vx=r[1],y=r[2],vy=r[3],z=r[4],vz=r[5];
    double v=sqrt(vx*vx+vy*vy+vz*vz),R=sqrt(x*x+y*y+z*z);
    double dRdt=(x*vx+y*vy+z*vz)/R;
    double A1=4*m/(c*c*R)-v*v/(c*c);
    double A2=-9*pow(m/(c*c*R),2)+2*m*dRdt*dRdt/(c*c*c*c*R);
    double A3=16*pow(m/(c*c*R),3)-pow(m/(c*c*R),2)*dRdt*dRdt/(c*c);
    double B1=4*dRdt/c;
    double B2=-2*dRdt*m/(c*c*c*R);
    double B3=4*pow(m/(c*c*R),2)*dRdt/c;
    switch(n)
    {
        case 1 :
            coes[0]=A1,coes[1]=B1;
            break;
        case 2 :
            coes[0]=A1+A2,coes[1]=B1+B2;
            break;
        case 3 :
            coes[0]=A1+A2+A3,coes[1]=B1+B2+B3;
            break;
    }
    
}
void PN2(int n,double c,double m,double coes[4], double r3[6],double r4[6])
{
    //coes[0]=A3,coes[1]=B3,coes[2]=A4,coes[3]=B4;
    double x3=r3[0],y3=r3[2],z3=r3[4],vx3=r3[1],vy3=r3[3],vz3=r3[5];
    double x4=r4[0],y4=r4[2],z4=r4[4],vx4=r4[1],vy4=r4[3],vz4=r4[5];
    double v3=sqrt(vx3*vx3+vy3*vy3+vz3*vz3);
    double v4=sqrt(vx4*vx4+vy4*vy4+vz4*vz4);
    double R3=sqrt(x3*x3+y3*y3+z3*z3),R4=sqrt(x4*x4+y4*y4+z4*z4);
    double dR3dt=(x3*vx3+y3*vy3+z3*vz3)/R3;
    double dR4dt=(x4*vx4+y4*vy4+z4*vz4)/R4;
    double A31=4*m/(c*c*R3)-v3*v3/(c*c);
    double A32=-9*pow(m/(c*c*R3),2)+2*m*dR3dt*dR3dt/(c*c*c*c*R3);
    double A33=16*pow(m/(c*c*R3),3)-pow(m/(c*c*R3),2)*dR3dt*dR3dt/(c*c);
    double B31=4*dR3dt/c;
    double B32=-2*dR3dt*m/(c*c*c*R3);
    double B33=4*pow(m/(c*c*R3),2)*dR3dt/c;
    double A41=4*m/(c*c*R4)-v4*v4/(c*c);
    double A42=-9*pow(m/(c*c*R4),2)+2*m*dR4dt*dR4dt/(c*c*c*c*R4);
    double A43=16*pow(m/(c*c*R4),3)-pow(m/(c*c*R4),2)*dR4dt*dR4dt/(c*c);
    double B41=4*dR4dt/c;
    double B42=-2*dR4dt*m/(c*c*c*R4);
    double B43=4*pow(m/(c*c*R4),2)*dR4dt/c;
    switch(n)
    {
        case 1 :
            coes[0]=A31,coes[1]=B31,coes[2]=A41,coes[3]=B41;
            break;
        case 2 :
            coes[0]=A31+A32,coes[1]=B31+B32,coes[2]=A41+A42,coes[3]=B41+B42;
            break;
        case 3 :
            coes[0]=A31+A32+A33,coes[1]=B31+B32+B33,coes[2]=A41+A42+A43,coes[3]=B41+B42+B43;
            break;
    }

}
struct rhs1
{
    const VecDoub m;
    const double ai;
    const double af;
    const double tau;
    const double c;
    double m1,m2;
    bool PN;
    const int orderofPN;
    rhs1(const VecDoub &mm,double aii,double aff,double tauu,double cc,bool PNN,const int n) :
        m(mm),ai(aii),af(aff),tau(tauu),c(cc),PN(PNN),orderofPN(n)
        {m1=m[0],m2=m[1];}
    void operator() (const double xx, VecDoub &yy, VecDoub &dydx)
    {
        double t=yy[5],phi=yy[4];
		double a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
		if(yy[4]>2.0*M_PI)
            yy[4]=yy[4]-2.0*M_PI;
        double x12=a*cos(phi),y12=a*sin(phi),dadt=(af-ai)/tau*exp(-t/tau);
		double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
        double x=yy[0],y=yy[2],z=yy[6],vx=yy[1],vy=yy[3],vz=yy[7];
        double ex=yy[8],ey=yy[9],ez=yy[10],E2=yy[11];
        double R1=sqrt(x*x+y*y+z*z),R2=sqrt((x-x12)*(x-x12)+(y-y12)*(y-y12)+z*z);
        double gx=m2*(x12-x)/pow(R2,3)-m2*x12/pow(a,3);
        double gy=m2*(y12-y)/pow(R2,3)-m2*y12/pow(a,3);
        double gz=-m2*z/pow(R2,3);
        
        if(PN==true)
        {
            double r[6]={x,vx/R1,y,vy/R1,z,vz/R1},coes[2];
            PN1(orderofPN,c,m1,coes,r);
            double nx=x/R1,ny=y/R1,nz=z/R1;
            double fx_pn=m1/pow(R1,2)*(coes[0]*nx+coes[1]*vx/(c*R1));
            double fy_pn=m1/pow(R1,2)*(coes[0]*ny+coes[1]*vy/(c*R1));
            double fz_pn=m1/pow(R1,2)*(coes[0]*nz+coes[1]*vz/(c*R1));
            gx=gx+fx_pn,gy=gy+fy_pn,gz=gz+fz_pn;
			double r2[6]={x-x12,vx/R1-vx12,y-y12,vy/R1-vy12,z,vz/R1};
            PN1(orderofPN,c,m2,coes,r2);
			nx=(x-x12)/R2,ny=(y-y12)/R2,nz=z/R2;
            fx_pn=m2/pow(R2,2)*(coes[0]*nx+coes[1]*(vx/R1-vx12)/c);
            fy_pn=m2/pow(R2,2)*(coes[0]*ny+coes[1]*(vy/R1-vy12)/c);
            fz_pn=m2/pow(R2,2)*(coes[0]*nz+coes[1]*vz/R1/c);
			gx=gx+fx_pn,gy=gy+fy_pn,gz=gz+fz_pn;
        }
        double temp1=vx*gx+vy*gy+vz*gz,temp2=x*gx+y*gy+z*gz,temp3=vx*x+vy*y+vz*z;
        double efx=2*x*temp1-vx*temp2-gx*temp3;
        double efy=2*y*temp1-vy*temp2-gy*temp3;
        double efz=2*z*temp1-vz*temp2-gz*temp3;
        dydx[0]=yy[1];
        dydx[1]=2*E2*x-ex+R1*R1*gx;
        dydx[2]=yy[3];
        dydx[3]=2*E2*y-ey+R1*R1*gy;
		dydx[4]=R1*w;
		dydx[5]=R1;
        dydx[6]=yy[7];
        dydx[7]=2*E2*z-ez+R1*R1*gz;
        dydx[8]=efx;
        dydx[9]=efy;
        dydx[10]=efz;
        dydx[11]=temp1;
    }
};
struct rhs2
{
    const VecDoub m;
    const double ai;
    const double af;
    const double tau;
    const double c;
    double m1,m2,m3,m4,m34;
    bool PN;
    const int orderofPN;
    rhs2(VecDoub mm,double aii,double aff,double tauu,double cc,bool PNN,const int n) :
    m(mm),ai(aii),af(aff),tau(tauu),c(cc),PN(PNN),orderofPN(n)
    {m1=m[0],m2=m[1],m3=m[2],m4=m[3],m34=m3+m4;}
    void operator() (const double xx, VecDoub &yy, VecDoub &dydx)
    {

        double t=yy[5];
        double a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
        if(yy[4]>2*M_PI)
            yy[4]=yy[4]-2.0*M_PI;
        double phi=yy[4],dadt=(af-ai)/tau*exp(-t/tau);
        double x12=a*cos(phi),y12=a*sin(phi);
        double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
        double X=yy[0],Y=yy[2],Z=yy[14],x=yy[6],y=yy[8],z=yy[16];
        double VX=yy[1],VY=yy[3],VZ=yy[15],vx=yy[7],vy=yy[9],vz=yy[17];
        double x13=X-m4/m34*x,y13=Y-m4/m34*y,z13=Z-m4/m34*z;
        double x14=X+m3/m34*x,y14=Y+m3/m34*y,z14=Z+m3/m34*z;
        double vx13=VX-m4/m34*vx,vy13=VY-m4/m34*vy,vz13=VZ-m4/m34*vz;
        double vx14=VX+m3/m34*vx,vy14=VY+m3/m34*vy,vz14=VZ+m3/m34*vz;
        double r13=sqrt(x13*x13+y13*y13+z13*z13);
        double r23=sqrt((x13-x12)*(x13-x12)+(y13-y12)*(y13-y12)+z13*z13);
        double r14=sqrt(x14*x14+y14*y14+z14*z14);
        double r24=sqrt((x14-x12)*(x14-x12)+(y14-y12)*(y14-y12)+(z14*z14));
        double r=sqrt(x*x+y*y+z*z),drdtau=(x*vx+y*vy+z*vz)/r;
        double Fx31=0.0,Fy31=0.0,Fz31=0.0,Fx41=0.0,Fy41=0.0,Fz41=0.0;
        double Fx32=0.0,Fy32=0.0,Fz32=0.0,Fx42=0.0,Fy42=0.0,Fz42=0.0;
        if(PN==true)
        {
            double r31[6]={x13,vx13/r,y13,vy13/r,z13,vz13/r};
            double r41[6]={x14,vx14/r,y14,vy14/r,z14,vz14/r},coes[4];
            PN2(orderofPN,c,m1,coes,r31,r41);
            double A3=coes[0],B3=coes[1],A4=coes[2],B4=coes[3];
            double nx13=x13/r13,ny13=y13/r13,nz13=z13/r13,nx14=x14/r14,ny14=y14/r14,nz14=z14/r14;
            Fx31=m1*(A3*nx13+B3*vx13/(c*r))/pow(r13,2);
            Fy31=m1*(A3*ny13+B3*vy13/(c*r))/pow(r13,2);
            Fz31=m1*(A3*nz13+B3*vz13/(c*r))/pow(r13,2);
            Fx41=m1*(A4*nx14+B4*vx14/(c*r))/pow(r14,2);
            Fy41=m1*(A4*ny14+B4*vy14/(c*r))/pow(r14,2);
            Fz41=m1*(A4*nz14+B4*vz14/(c*r))/pow(r14,2);

            double r32[6]={x13-x12,vx13/r-vx12,y13-y12,vy13/r-vy12,z13,vz13/r};
            double r42[6]={x14-x12,vx14/r-vx12,y14-y12,vy14/r-vy12,z14,vz14/r};
            PN2(orderofPN,c,m2,coes,r32,r42);
            A3=coes[0],B3=coes[1],A4=coes[2],B4=coes[3];
            double nx23=(x13-x12)/r23,ny23=(y13-y12)/r23,nz23=z13/r23;
            double nx24=(x14-x12)/r24,ny24=(y14-y12)/r24,nz24=z14/r24;
            Fx32=m2*(A3*nx23+B3*(vx13/r-vx12)/c)/pow(r23,2);
            Fy32=m2*(A3*ny23+B3*(vy13/r-vy12)/c)/pow(r23,2);
            Fz32=m2*(A3*nz23+B3*vz13/(c*r))/pow(r23,2);
            Fx42=m2*(A4*nx24+B4*(vx14/r-vx12)/c)/pow(r24,2);
            Fy42=m2*(A4*ny24+B4*(vy14/r-vy12)/c)/pow(r24,2);
            Fz42=m2*(A4*nz24+B4*vz14/(c*r))/pow(r24,2);
        }

        double ex=yy[11],ey=yy[12],ez=yy[13],E2=yy[10];
        double fx=-m3*m1/m34*x13/pow(r13,3)-m4*m1/m34*x14/pow(r14,3)+m2*m3/m34*(x12-x13)/pow(r23,3)
                +m2*m4/m34*(x12-x14)/pow(r24,3)-m2*x12/pow(a,3)+m3/m34*(Fx31+Fx32)+m4/m34*(Fx41+Fx42);
        double fy=-m3*m1/m34*y13/pow(r13,3)-m4*m1/m34*y14/pow(r14,3)+m2*m3/m34*(y12-y13)/pow(r23,3)
                +m2*m4/m34*(y12-y14)/pow(r24,3)-m2*y12/pow(a,3)+m3/m34*(Fy31+Fy32)+m4/m34*(Fy41+Fy42);
        double fz=-m3*m1/m34*z13/pow(r13,3)-m4*m1/m34*z14/pow(r14,3)-m2*m3/m34*z13/pow(r23,3)
                -m2*m4/m34*z14/pow(r24,3)+m3/m34*(Fz31+Fz32)+m4/m34*(Fz41+Fz42);
        double gx=m1*(x13/pow(r13,3)-x14/pow(r14,3))+m2*((x12-x14)/pow(r24,3)-(x12-x13)/pow(r23,3))+Fx41-Fx31+Fx42-Fx32;
        double gy=m1*(y13/pow(r13,3)-y14/pow(r14,3))+m2*((y12-y14)/pow(r24,3)-(y12-y13)/pow(r23,3))+Fy41-Fy31+Fy42-Fy42;
        double gz=m1*(z13/pow(r13,3)-z14/pow(r14,3))+m2*(-z14/pow(r24,3)+z13/pow(r23,3))+Fz41-Fz31+Fz42-Fz32;
        double temp1=vx*gx+vy*gy+vz*gz,temp2=x*gx+y*gy+z*gz,temp3=vx*x+vy*y+vz*z;
        double efx=2.0*x*temp1-vx*temp2-gx*temp3;
        double efy=2.0*y*temp1-vy*temp2-gy*temp3;
        double efz=2.0*z*temp1-vz*temp2-gz*temp3;
        dydx[0]=yy[1];
        dydx[1]=r*r*fx+drdtau/r*yy[1];
        dydx[2]=yy[3];
        dydx[3]=r*r*fy+drdtau/r*yy[3];
        dydx[4]=r*w;
        dydx[5]=r;
        dydx[6]=yy[7];
        dydx[7]=2.0*E2*x-ex+r*r*gx;
        dydx[8]=yy[9];
        dydx[9]=2.0*E2*y-ey+r*r*gy;
        dydx[10]=temp1;
        dydx[11]=efx;
        dydx[12]=efy;
        dydx[13]=efz;
        dydx[14]=yy[15];
        dydx[15]=r*r*fz+drdtau/r*yy[15];
        dydx[16]=yy[17];
        dydx[17]=2.0*E2*z-ez+r*r*gz;

    }

};
struct rhs3
{
    const VecDoub m;
    const double ai;
    const double af;
    const double tau;
    const double c;
    double m1,m2;
    bool PN;
    const int orderofPN;
    rhs3(const VecDoub &mm,double aii,double aff,double tauu,double cc,bool PNN,const int n) :
        m(mm),ai(aii),af(aff),tau(tauu),c(cc),PN(PNN),orderofPN(n)
        {m1=m[0],m2=m[1];}
    void operator() (const double xx, VecDoub &yy, VecDoub &dydx)
    {
        double t=yy[5];
        double a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
        if(yy[4]>2*M_PI)
            yy[4]=yy[4]-2.0*M_PI;
        double phi=yy[4],m12=m1+m2,dadt=(af-ai)/tau*exp(-t/tau);
        double x12=a*cos(phi),y12=a*sin(phi),vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
        double x1=-m2/m12*x12,y1=-m2/m12*y12,vx1=-m2/m12*vx12,vy1=-m2/m12*vy12;
        double x2=m1/m12*x12,y2=m1/m12*y12,vx2=m1/m12*vx12,vy2=m1/m12*vy12;
        //m1=0,x1=0,y1=0,vx1=0,vy1=0,m2=4*M_PI*M_PI,x2=100,y2=0,vx2=0,vy2=0;
        double x3=yy[0],y3=yy[2],z3=yy[6];
        double r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
        double vx3=yy[1]/r13,vy3=yy[3]/r13,vz3=yy[7]/r13;
        double r23=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+z3*z3);
        double dr13dt=((x3-x1)*(vx3-vx1)+(y3-y1)*(vy3-vy1)+z3*vz3)/r13;
        double fx=m1/pow(r13,3)*(x1-x3)+m2/pow(r23,3)*(x2-x3);
        double fy=m1/pow(r13,3)*(y1-y3)+m2/pow(r23,3)*(y2-y3);
        double fz=-m1/pow(r13,3)*z3-m2/pow(r23,3)*z3;
        if(PN==true)
        {
            double r[6]={x3-x1,vx3-vx1,y3-y1,vy3-vy1,z3,vz3},coes[2];
            PN1(orderofPN,c,m1,coes,r);
            double nx=(x3-x1)/r13,ny=(y3-y1)/r13,nz=z3/r13;
            double fx_pn=m1/pow(r13,2)*(coes[0]*nx+coes[1]*(vx3-vx1)/c);
            double fy_pn=m1/pow(r13,2)*(coes[0]*ny+coes[1]*(vy3-vy1)/c);
            double fz_pn=m1/pow(r13,2)*(coes[0]*nz+coes[1]*vz3/c);
            fx=fx+fx_pn,fy=fy+fy_pn,fz=fz+fz_pn;
            double r2[6]={x3-x2,vx3-vx2,y3-y2,vy3-vy2,z3,vz3};
            PN1(orderofPN,c,m2,coes,r2);
            nx=(x3-x2)/r23,ny=(y3-y2)/r23,nz=z3/r23;
            fx_pn=m2/pow(r23,2)*(coes[0]*nx+coes[1]*(vx3-vx2)/c);
            fy_pn=m2/pow(r23,2)*(coes[0]*ny+coes[1]*(vy3-vy2)/c);
            fz_pn=m2/pow(r23,2)*(coes[0]*nz+coes[1]*vz3/c);
            fx=fx+fx_pn,fy=fy+fy_pn,fz=fz+fz_pn;
        }
        dydx[0]=yy[1];
        dydx[1]=r13*r13*fx+yy[1]*dr13dt;
        dydx[2]=yy[3];
        dydx[3]=r13*r13*fy+yy[3]*dr13dt;
        dydx[4]=w*r13;
        dydx[5]=r13;
        dydx[6]=yy[7];
        dydx[7]=r13*r13*fz+yy[7]*dr13dt;
    }
};
struct rhs4
{
    const VecDoub m;
    const double ai;
    const double af;
    const double tau;
    const double c;
    double m1,m2,m3,m4,m34,m12;
    bool PN;
    const int orderofPN;
    rhs4(VecDoub mm,double aii,double aff,double tauu,double cc,bool PNN,const int n) :
    m(mm),ai(aii),af(aff),tau(tauu),c(cc),PN(PNN),orderofPN(n)
    {m1=m[0],m2=m[1],m3=m[2],m4=m[3],m34=m3+m4,m12=m1+m2;}
    void operator() (const double xx, VecDoub &yy, VecDoub &dydx)
    {

        double t=yy[5];
        double a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
        if(yy[4]>2*M_PI)
            yy[4]=yy[4]-2*M_PI;
        double phi=yy[4],dadt=(af-ai)/tau*exp(-t/tau);
        double x12=a*cos(phi),y12=a*sin(phi),vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
        double x1=-m2/m12*x12,y1=-m2/m12*y12,vx1=-m2/m12*vx12,vy1=-m2/m12*vy12;
        double x2=m1/m12*x12,y2=m1/m12*y12,vx2=m1/m12*vx12,vy2=m1/m12*vy12;
        double X=yy[0],Y=yy[2],Z=yy[14],x=yy[6],y=yy[8],z=yy[16];
        double VX=yy[1],VY=yy[3],VZ=yy[15],vx=yy[7],vy=yy[9],vz=yy[17];
        double x3=X-m4/m34*x,y3=Y-m4/m34*y,z3=Z-m4/m34*z;
        double x4=X+m3/m34*x,y4=Y+m3/m34*y,z4=Z+m3/m34*z;
        double vx3=VX-m4/m34*vx,vy3=VY-m4/m34*vy,vz3=VZ-m4/m34*vz;
        double vx4=VX+m3/m34*vx,vy4=VY+m3/m34*vy,vz4=VZ+m3/m34*vz;
        double r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
        double r23=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+z3*z3);
        double r14=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+z4*z4);
        double r24=sqrt((x4-x2)*(x4-x2)+(y4-y2)*(y4-y2)+z4*z4);
        double r=sqrt(x*x+y*y+z*z),drdtau=(x*vx+y*vy+z*vz)/r;
        double Fx31=0.0,Fy31=0.0,Fz31=0.0,Fx41=0.0,Fy41=0.0,Fz41=0.0;
        double Fx32=0.0,Fy32=0.0,Fz32=0.0,Fx42=0.0,Fy42=0.0,Fz42=0.0;
        if(PN==true)
        {
            double r31[6]={x3-x1,vx3/r-vx1,y3-y1,vy3/r-vy1,z3,vz3/r};
            double r41[6]={x4-x1,vx4/r-vx1,y4-y1,vy4/r-vy1,z4,vz4/r},coes[4];
            PN2(orderofPN,c,m1,coes,r31,r41);
            double nx31=(x3-x1)/r13,ny31=(y3-y1)/r13,nz31=z3/r13,nx41=(x4-x1)/r14,ny41=(y4-y1)/r14,nz41=z4/r14;
            Fx31=m1*(coes[0]*nx31+coes[1]*(vx3/r-vx1)/c)/pow(r13,2);
            Fy31=m1*(coes[0]*ny31+coes[1]*(vy3/r-vy1)/c)/pow(r13,2);
            Fz31=m1*(coes[0]*nz31+coes[1]*vz3/(c*r))/pow(r13,2);
            Fx41=m1*(coes[2]*nx41+coes[3]*(vx4/r-vx1)/c)/pow(r14,2);
            Fy41=m1*(coes[2]*ny41+coes[3]*(vy4/r-vy1)/c)/pow(r14,2);
            Fz41=m1*(coes[2]*nz41+coes[3]*vz4/(c*r))/pow(r14,2);
            double r32[6]={x3-x2,vx3/r-vx2,y3-y2,vy3/r-vy2,z3,vz3/r};
            double r42[6]={x4-x2,vx4/r-vx2,y4-y2,vy4/r-vy2,z4,vz4/r};
            PN2(orderofPN,c,m2,coes,r32,r42);
            double nx32=(x3-x2)/r23,ny32=(y3-y2)/r23,nz32=z3/r23,nx42=(x4-x2)/r24,ny42=(y4-y2)/r24,nz42=z4/r24;
            Fx32=m2*(coes[0]*nx32+coes[1]*(vx3/r-vx2)/c)/pow(r23,2);
            Fy32=m2*(coes[0]*ny32+coes[1]*(vy3/r-vy2)/c)/pow(r23,2);
            Fz32=m2*(coes[0]*nz32+coes[1]*vz3/(c*r))/pow(r23,2);
            Fx42=m2*(coes[2]*nx42+coes[3]*(vx4/r-vx2)/c)/pow(r24,2);
            Fy42=m2*(coes[2]*ny42+coes[3]*(vy4/r-vy2)/c)/pow(r24,2);
            Fz42=m2*(coes[2]*nz42+coes[3]*vz4/(c*r))/pow(r24,2);
        }

        double ex=yy[11],ey=yy[12],ez=yy[13],E2=yy[10];
        double fx=-m3*m1/m34*(x3-x1)/pow(r13,3)-m4*m1/m34*(x4-x1)/pow(r14,3)+m2*m3/m34*(x2-x3)/pow(r23,3)
                +m2*m4/m34*(x2-x4)/pow(r24,3)+m3/m34*(Fx31+Fx32)+m4/m34*(Fx41+Fx42);
        double fy=-m3*m1/m34*(y3-y1)/pow(r13,3)-m4*m1/m34*(y4-y1)/pow(r14,3)+m2*m3/m34*(y2-y3)/pow(r23,3)
                +m2*m4/m34*(y2-y4)/pow(r24,3)+m3/m34*(Fy31+Fy32)+m4/m34*(Fy41+Fy42);
        double fz=-m3*m1/m34*z3/pow(r13,3)-m4*m1/m34*z4/pow(r14,3)-m2*m3/m34*z3/pow(r23,3)
                -m2*m4/m34*z4/pow(r24,3)+m3/m34*(Fz31+Fz32)+m4/m34*(Fz41+Fz42);
        double gx=m1*((x3-x1)/pow(r13,3)-(x4-x1)/pow(r14,3))+m2*((x2-x4)/pow(r24,3)-(x2-x3)/pow(r23,3))+Fx41-Fx31+Fx42-Fx32;
        double gy=m1*((y3-y1)/pow(r13,3)-(y4-y1)/pow(r14,3))+m2*((y2-y4)/pow(r24,3)-(y2-y3)/pow(r23,3))+Fy41-Fy31+Fy42-Fy42;
        double gz=m1*(z3/pow(r13,3)-z4/pow(r14,3))+m2*(-z4/pow(r24,3)+z3/pow(r23,3))+Fz41-Fz31+Fz42-Fz32;
        double temp1=vx*gx+vy*gy+vz*gz,temp2=x*gx+y*gy+z*gz,temp3=vx*x+vy*y+vz*z;
        double efx=2.0*x*temp1-vx*temp2-gx*temp3;
        double efy=2.0*y*temp1-vy*temp2-gy*temp3;
        double efz=2.0*z*temp1-vz*temp2-gz*temp3;
        dydx[0]=yy[1];
        dydx[1]=r*r*fx+drdtau/r*yy[1];
        dydx[2]=yy[3];
        dydx[3]=r*r*fy+drdtau/r*yy[3];
        dydx[4]=r*w;
        dydx[5]=r;
        dydx[6]=yy[7];
        dydx[7]=2.0*E2*x-ex+r*r*gx;
        dydx[8]=yy[9];
        dydx[9]=2.0*E2*y-ey+r*r*gy;
        dydx[10]=temp1;
        dydx[11]=efx;
        dydx[12]=efy;
        dydx[13]=efz;
        dydx[14]=yy[15];
        dydx[15]=r*r*fz+drdtau/r*yy[15];
        dydx[16]=yy[17];
        dydx[17]=2.0*E2*z-ez+r*r*gz;

    }

};
template <class T>
Doub rtnewt(T &funcd, const Doub x1, const Doub x2, const Doub xacc) {
	const Int JMAX=20;
	Doub rtn=0.5*(x1+x2);
	for (Int j=0;j<JMAX;j++) {
		Doub f=funcd(rtn);
		Doub df=funcd.df(rtn);
		Doub dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			throw("Jumped out of brackets in rtnewt");
		if (abs(dx) < xacc) return rtn;
	}
	throw("Maximum number of iterations exceeded in rtnewt");
}
struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64();}
	inline Uint int32() { return (Uint)int64(); }
};
struct Normaldev : Ran {
	Doub mu,sig;
	Normaldev(Doub mmu, Doub ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig){}
	Doub dev() {
		Doub u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = abs(v) + 0.386595;
			q = SQR(x) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597
			&& (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
		return mu + sig*v/u;
	}
};
struct powerlaw
{
    const double alpha;
    powerlaw(double a) : alpha(a) {}
    double operator()(double x)
    {
        return pow(x,alpha);
    }
};
void pow_dis(double left,double right,double alpha,int seed, vector<double> &pd,int flag, vector<double> M)//alpha !=-1
{
	powerlaw myf(alpha);
	double a=qromb(myf,left,right,1e-10);
	double x_min=pow(left,alpha+1)/(alpha+1)/a;
	Ran myrandom(seed);
    if(flag==0)
    {
        double temp=0.0;
        for(int i=0;i<pd.size();i++)
	    {
		    temp=myrandom.doub()+x_min;
		    pd[i]=round(pow((alpha+1)*a*temp,1.0/(alpha+1))*1e6)/1e6;
	    }
    }
    else
    {
        double m_min=0.1;
        double temp=0.0;
        for(int i=0;i<pd.size();i++)
	    {
            for(int j=0;j<20;j++)
            {
                temp=myrandom.doub()+x_min;
                if(pow((alpha+1)*a*temp,1.0/(alpha+1))*M[i]<m_min)
                   continue;
                else
                   {pd[i]=round(pow((alpha+1)*a*temp,1.0/(alpha+1))*1e6)/1e6;break;}
                if(j==20)
                    throw("check mass ratio!");
            }
	    }
    }
}
void semimajor_axes(int seed, vector<double> &pd, vector<double> M)
{
    powerlaw myf(-1);
    Ran myrandom(seed);
    for(int i=0;i<pd.size();i++)
    {
        double mp=M[i];
        double left=2.5*pow(mp,0.6);
        double right=2000;
        double a=qromb(myf,left,right,1e-10);
        double x_min=log(left)/a;
        double temp=myrandom.doub()+x_min;
        pd[i]=round(exp(a*temp)*1e6)/1e6;
    }
}

void transform_matrix1(double phi,double w,double I,double TM[3][3])//tranform from (x,y,z) to (X,Y,Z)
{
    double P1[3][3]={{cos(w),-sin(w),0.0},{sin(w),cos(w),0.0},{0.0,0.0,1.0}};
    double P2[3][3]={{1.0,0.0,0.0},{0.0,cos(I),-sin(I)},{0.0,sin(I),cos(I)}};
    double P3[3][3]={{cos(phi),-sin(phi),0.0},{sin(phi),cos(phi),0.0},{0.0,0.0,1.0}};
    double temp[3][3];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
    {
        temp[i][j]=P2[i][0]*P1[0][j]+P2[i][1]*P1[1][j]+P2[i][2]*P1[2][j];
    }
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
    {
        TM[i][j]=P3[i][0]*temp[0][j]+P3[i][1]*temp[1][j]+P3[i][2]*temp[2][j];
    }
}
void transform_matrix2(double phi,double w,double I,double TM[3][3])//tranform from (X,Y,Z) to (x,y,z)
{
    double rP1[3][3]={{cos(w),sin(w),0.0},{-sin(w),cos(w),0.0},{0.0,0.0,1.0}};
    double rP2[3][3]={{1.0,0.0,0.0},{0.0,cos(I),sin(I)},{0.0,-sin(I),cos(I)}};
    double rP3[3][3]={{cos(phi),sin(phi),0.0},{-sin(phi),cos(phi),0.0},{0.0,0.0,1.0}};
    double temp[3][3];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
    {
        temp[i][j]=rP2[i][0]*rP3[0][j]+rP2[i][1]*rP3[1][j]+rP2[i][2]*rP3[2][j];
    }
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
    {
        TM[i][j]=rP1[i][0]*temp[0][j]+rP1[i][1]*temp[1][j]+rP1[i][2]*temp[2][j];
    }
}
void Kepler_Sol(double m,double kp[6],double r[6])
{
    double a=kp[0],e=kp[1],I=kp[2],phi=kp[3],w=kp[4],f=kp[5];
    double n=sqrt(m/(a*a*a));
    double x=a*(1-e*e)/(1+e*cos(f))*cos(f),y=a*(1-e*e)/(1+e*cos(f))*sin(f);
    double vx=-n*a/sqrt(1-e*e)*sin(f),vy=n*a/sqrt(1-e*e)*(e+cos(f));
    double TM[3][3];
    transform_matrix1(phi,w,I,TM);
    r[0]=TM[0][0]*x+TM[0][1]*y,r[2]=TM[1][0]*x+TM[1][1]*y,r[4]=TM[2][0]*x+TM[2][1]*y;
    r[1]=TM[0][0]*vx+TM[0][1]*vy,r[3]=TM[1][0]*vx+TM[1][1]*vy,r[5]=TM[2][0]*vx+TM[2][1]*vy;
}
struct Fundcd
{
    double n;
    double t;
    double e;
    double tau;
    void change_parameters(double nn,double tt,double ee,double tauu)
    {
        n=nn;
        t=tt;
        e=ee;
        tau=tauu;
    }
    double operator()(const double E)
    {
        return E-e*sin(E)-n*(t-tau);
    }
    double df(const double E)
    {
        return 1-e*cos(E);
    }
};
void Kepler_Solver(double t,double P,const VecDoub &kp,Fundcd &fx,double r[6])
{
    double tempt=t,a=kp[0],e=kp[1],n=2*M_PI/P,tau=kp[6],phi=kp[7],w=kp[4],I=kp[2];
    if((t-tau)>P)
        tempt=t-(ceil((t-tau)/P)-1)*P;
    fx.change_parameters(n,tempt,e,tau);
    double E=rtnewt(fx,-0.5*M_PI,2.5*M_PI,1e-12),dEdt=n/(1-e*cos(E));
    double x=a*(cos(E)-e),y=a*sqrt(1-e*e)*sin(E);
    double vx=-a*sin(E)*dEdt,vy=a*sqrt(1-e*e)*cos(E)*dEdt;
    double TM[3][3];
    transform_matrix1(phi,w,I,TM);
    r[0]=TM[0][0]*x+TM[0][1]*y,r[2]=TM[1][0]*x+TM[1][1]*y,r[4]=TM[2][0]*x+TM[2][1]*y;
    r[1]=TM[0][0]*vx+TM[0][1]*vy,r[3]=TM[1][0]*vx+TM[1][1]*vy,r[5]=TM[2][0]*vx+TM[2][1]*vy;

}
void Kepler_Parameters_2D(double yy[6],VecDoub &kp,double m,double t)
{
    double x=yy[0],vx=yy[1],y=yy[2],vy=yy[3];
    double r=sqrt(x*x+y*y);
    double v=sqrt(vx*vx+vy*vy);
    double drdt=(x*vx+y*vy)/r;
    double hz=x*vy-y*vx;
    double h_M=abs(hz);
    double energy=-m/r+v*v/2.0;
    double a=-m/(2.0*energy);
    double e=sqrt(abs(1-h_M*h_M/(a*m)));
    double I=acos(hz/h_M);
    double f=0.0,w=0.0,half_E=0.0;
    if(e<1e-6)
        f=0.0;
    else
    {
        f=acos(1.0/e*(a*(1-e*e)/r-1));

        if(1.0/e*(a*(1-e*e)/r-1)-1.0>0.0)
            f=0.0;
        if(1.0/e*(a*(1-e*e)/r-1)+1.0<0.0)
            f=M_PI;
        if(drdt<0)
            f=2*M_PI-f;
    }
    double w_f=acos(x/r);
    if(y<0)
        w_f=2*M_PI-w_f;
    if(hz>0)
    {
        w=w_f-f;
        if(w<0)
            w=w+2*M_PI;
        half_E=atan(sqrt((1-e)/(1+e))*tan(f/2.0));
        if(half_E<0)
            half_E=half_E+M_PI;
    }
    else
    {
        w=w_f+f;
        if(w>2*M_PI)
            w=w-2*M_PI;
        half_E=atan(sqrt((1-e)/(1+e))*tan(f/2.0));
        if(half_E<0)
        half_E=half_E+M_PI;
    }
    double E=2.0*half_E;
    kp[0]=a,kp[1]=e,kp[2]=I,kp[3]=f,kp[4]=0.0,kp[5]=E;
    double n=sqrt(m/(a*a*a));
    double tau=t-(E-e*sin(E))/n;
    kp[6]=tau,kp[7]=w;
}
bool Kepler_Parameters_3D(double yy[6],VecDoub &kp,double m,double t)
{

    double x=yy[0],vx=yy[1],y=yy[2],vy=yy[3],z=yy[4],vz=yy[5];
    double r=sqrt(x*x+y*y+z*z);
    double v=sqrt(vx*vx+vy*vy+vz*vz);
    double drdt=(x*vx+y*vy+z*vz)/r;
    double h[3]={y*vz-z*vy,z*vx-x*vz,x*vy-y*vx};
    double h_M=sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
    double energy=-m/r+v*v/2.0;
    if(energy>0)
        return false;
    double a=-m/(2.0*energy);
    double e=sqrt(abs(1-h_M*h_M/(a*m)));
    double I=acos(h[2]/h_M);
    if(I<1e-7*M_PI || I>(1.0-1e-7)*M_PI)
    {
        Kepler_Parameters_2D(yy,kp,m,t);
        return true;
    }
    double f=0.0,phi=0.0,w=0.0;
    if(e<1e-6)
        f=0.0;
    else
    {
        f=acos(1.0/e*(a*(1-e*e)/r-1));

        if(1.0/e*(a*(1-e*e)/r-1)-1.0>0.0)
            f=0.0;
        if(1.0/e*(a*(1-e*e)/r-1)+1.0<0.0)
            f=M_PI;
        if(drdt<0)
            f=2*M_PI-f;

    }
    if(h[0]>0)
    {
        phi=acos(-h[1]/(h_M*sin(I)));
        if(-h[1]/(h_M*sin(I))-1.0>0)
            phi=0.0;
        if(-h[1]/(h_M*sin(I))+1.0<0)
            phi=M_PI;
        double w_f=0.0;
        if(abs(cos(phi))<1e-6)
            {
                w_f=acos(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I)));

                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }
        else
            {
                w_f=acos(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I)));

                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }

        if(z/sin(I)<0)
            w_f=2*M_PI-w_f;
        w=w_f-f;
    }
    else
    {
        phi=acos(-h[1]/(h_M*sin(I)));
        if(-h[1]/(h_M*sin(I))-1.0>0)
            phi=0.0;
        if(-h[1]/(h_M*sin(I))+1.0<0)
            phi=M_PI;
        phi=2*M_PI-phi;
        double w_f=0.0;
        if(cos(phi)<1e-6)
            {
                w_f=acos(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I)));
                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }
        else
            {
                w_f=acos(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I)));
                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }

        if(z/sin(I)<0)
            w_f=2*M_PI-w_f;
        w=w_f-f;
    }

    double half_E=atan(sqrt((1-e)/(1+e))*tan(f/2.0));
    if(half_E<0)
        half_E=half_E+M_PI;
    double E=2.0*half_E;
    double n=sqrt(m/(a*a*a));
    double tau=t-(E-e*sin(E))/n;

    kp[0]=a,kp[1]=e,kp[2]=I,kp[3]=f,kp[4]=w,kp[5]=E,kp[6]=tau,kp[7]=phi;
    return true;
}
void Kepler_2D(double m,double yy[6],double kp[6])
{

    double x=yy[0],vx=yy[1],y=yy[2],vy=yy[3];
    double r=sqrt(x*x+y*y);
    double v=sqrt(vx*vx+vy*vy);
    double drdt=(x*vx+y*vy)/r;
    double hz=x*vy-y*vx;
    double h_M=abs(hz);
    double energy=-m/r+v*v/2.0;
    double a=-m/(2.0*energy);
    double e=sqrt(abs(1-h_M*h_M/(a*m)));
    double I=acos(hz/h_M);
    double f=0.0,w=0.0;
    if(e<1e-6)
        f=0.0;
    else
    {
        f=acos(1.0/e*(a*(1-e*e)/r-1));

        if(1.0/e*(a*(1-e*e)/r-1)-1.0>0.0)
            f=0.0;
        if(1.0/e*(a*(1-e*e)/r-1)+1.0<0.0)
            f=M_PI;
        if(drdt<0)
            f=2*M_PI-f;
    }
    double w_f=acos(x/r);
    if(y<0)
        w_f=2*M_PI-w_f;
    if(hz>0)
    {
        w=w_f-f;
        if(w<0)
            w=w+2*M_PI;
    }
    else
    {
        w=w_f+f;
        if(w>2*M_PI)
            w=w-2*M_PI;

    }
    kp[0]=a,kp[1]=e,kp[2]=I,kp[3]=w,kp[4]=0.0,kp[5]=f;
    
}
void Kepler_3D(double m,double yy[6],double kp[6])
{

    double x=yy[0],vx=yy[1],y=yy[2],vy=yy[3],z=yy[4],vz=yy[5];
    double r=sqrt(x*x+y*y+z*z);
    double v=sqrt(vx*vx+vy*vy+vz*vz);
    double drdt=(x*vx+y*vy+z*vz)/r;
    double h[3]={y*vz-z*vy,z*vx-x*vz,x*vy-y*vx};
    double h_M=sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
    double energy=-m/r+v*v/2.0;
    double a=-m/(2.0*energy);
    double e=sqrt(abs(1-h_M*h_M/(a*m)));
    double I=acos(h[2]/h_M);
    double f=0.0,phi=0.0,w=0.0;
    if(e<1e-6)
        f=0.0;
    else
    {
        f=acos(1.0/e*(a*(1-e*e)/r-1));

        if(1.0/e*(a*(1-e*e)/r-1)-1.0>0.0)
            f=0.0;
        if(1.0/e*(a*(1-e*e)/r-1)+1.0<0.0)
            f=M_PI;
        if(drdt<0)
            f=2*M_PI-f;

    }
    if(h[0]>0)
    {
        phi=acos(-h[1]/(h_M*sin(I)));
        if(-h[1]/(h_M*sin(I))-1.0>0)
            phi=0.0;
        if(-h[1]/(h_M*sin(I))+1.0<0)
            phi=M_PI;
        double w_f=0.0;
        if(abs(cos(phi))<1e-6)
            {
                w_f=acos(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I)));

                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }
        else
            {
                w_f=acos(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I)));

                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }

        if(z/sin(I)<0)
            w_f=2*M_PI-w_f;
        w=w_f-f;
    }
    else
    {
        phi=acos(-h[1]/(h_M*sin(I)));
        if(-h[1]/(h_M*sin(I))-1.0>0)
            phi=0.0;
        if(-h[1]/(h_M*sin(I))+1.0<0)
            phi=M_PI;
        phi=2*M_PI-phi;
        double w_f=0.0;
        if(cos(phi)<1e-6)
            {
                w_f=acos(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I)));
                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/sin(phi)*(y/r-cos(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }
        else
            {
                w_f=acos(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I)));
                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))-1.0>0.0)
                    w_f=0.0;
                if(1.0/cos(phi)*(x/r+sin(phi)*z/(r*sin(I))*cos(I))+1.0<0.0)
                    w_f=M_PI;
            }

        if(z/sin(I)<0)
            w_f=2*M_PI-w_f;
        w=w_f-f;
    }
    if(w<0)
       w=w+2.0*M_PI;
    kp[0]=a,kp[1]=e,kp[2]=I,kp[3]=phi,kp[4]=w,kp[5]=f;
}
void var_transform_smbh(int n,const VecDoub ystart,VecDoub &ystart1,const VecDoub &m)//transform location and velocities to integral varibles
{
	if(n==1)
	{
        double m1=m[0],x=ystart[0],vx=ystart[1],y=ystart[2],vy=ystart[3],z=ystart[10],vz=ystart[11];
        double R=sqrt(x*x+y*y+z*z);
		double E2=(vx*vx+vy*vy+vz*vz)/2.0-m1/R;
		vx=vx*R,vy=vy*R,vz=vz*R;
		double v2=vx*vx+vy*vy+vz*vz,dRdtau=(x*vx+y*vy+z*vz)/R;
        double ex=v2*x/pow(R,2)-dRdtau/R*vx-m1*x/R;
        double ey=v2*y/pow(R,2)-dRdtau/R*vy-m1*y/R;
        double ez=v2*z/pow(R,2)-dRdtau/R*vz-m1*z/R;
        
		ystart1[0]=x,ystart1[1]=vx,ystart1[2]=y,ystart1[3]=vy,ystart1[4]=ystart[4],ystart1[5]=ystart[5];
		ystart1[6]=z,ystart1[7]=vz,ystart1[8]=ex,ystart1[9]=ey,ystart1[10]=ez,ystart1[11]=E2;
	}
	else
	{
		double m34=m[2]+m[3];
        double X=ystart[0],Y=ystart[2],Z=ystart[10],VX=ystart[1],VY=ystart[3],VZ=ystart[11],t=ystart[5],phi=ystart[4];
        double x=ystart[6],y=ystart[8],z=ystart[12],vx=ystart[7],vy=ystart[9],vz=ystart[13];
        double r=sqrt(x*x+y*y+z*z);
        double E2=(vx*vx+vy*vy+vz*vz)/2.0-m34/r;
        vx=vx*r,vy=vy*r,vz=vz*r;
        double v2=vx*vx+vy*vy+vz*vz,drdtau=(x*vx+y*vy+z*vz)/r;
        double ex=v2*x/pow(r,2)-drdtau/r*vx-m34*x/r;
        double ey=v2*y/pow(r,2)-drdtau/r*vy-m34*y/r;
        double ez=v2*z/pow(r,2)-drdtau/r*vz-m34*z/r;
        ystart1[0]=X,ystart1[1]=VX*r,ystart1[2]=Y,ystart1[3]=VY*r;
        ystart1[4]=phi,ystart1[5]=t,ystart1[6]=x,ystart1[7]=vx,ystart1[8]=y,ystart1[9]=vy;
        ystart1[10]=E2,ystart1[11]=ex,ystart1[12]=ey,ystart1[13]=ez,ystart1[14]=Z,ystart1[15]=VZ*r,ystart1[16]=z,ystart1[17]=vz;
	}
}
void frame_transform(VecDoub &ystart,const VecDoub &m,const double ai,const double af,const double tau)
{
    double m1=m[0],m2=m[1];
    double t=ystart[5],phi=ystart[4];
    double a=af+(ai-af)*exp(-t/tau),dadt=(af-ai)/tau*exp(-t/tau);
    double x12=a*cos(phi),y12=a*sin(phi),w=sqrt((m1+m2)/pow(a,3));
    double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
    double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
    double vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
    ystart[0]=ystart[0]+x1,ystart[1]=ystart[1]+vx1,ystart[2]=ystart[2]+y1,ystart[3]=ystart[3]+vy1;
}
void var_transform_cm(int n,const VecDoub y1,VecDoub &y2,const VecDoub &m,const double ai, const double af,const double tau)//transform location and velocities to integral varibles
{
    double m1=m[0],m2=m[1];
    if(n==1)
    {
        double t=y1[5],phi=y1[4];
        double a=af+(ai-af)*exp(-t/tau);
        double x12=a*cos(phi),y12=a*sin(phi);
        double x11=-m2/(m1+m2)*x12,y11=-m2/(m1+m2)*y12;
        double R1=sqrt((y1[0]-x11)*(y1[0]-x11)+(y1[2]-y11)*(y1[2]-y11)+y1[10]*y1[10]);

        y2[0]=y1[0],y2[1]=y1[1]*R1,y2[2]=y1[2],y2[3]=y1[3]*R1,y2[4]=y1[4],y2[5]=y1[5],y2[6]=y1[10],y2[7]=y1[11]*R1;
    }
    if(n==2)
    {
         double m34=m[2]+m[3];
         double X=y1[0],Y=y1[2],Z=y1[10],VX=y1[1],VY=y1[3],VZ=y1[11],t=y1[5],phi=y1[4];
         double x=y1[6],y=y1[8],z=y1[12],vx=y1[7],vy=y1[9],vz=y1[13];
         double r=sqrt(x*x+y*y+z*z);
         double E2=(vx*vx+vy*vy+vz*vz)/2.0-m34/r;
         vx=vx*r,vy=vy*r,vz=vz*r;
         double v2=vx*vx+vy*vy+vz*vz,drdtau=(x*vx+y*vy+z*vz)/r;
         double ex=v2*x/pow(r,2)-drdtau/r*vx-m34*x/r;
         double ey=v2*y/pow(r,2)-drdtau/r*vy-m34*y/r;
         double ez=v2*z/pow(r,2)-drdtau/r*vz-m34*z/r;
         y2[0]=X,y2[1]=VX*r,y2[2]=Y,y2[3]=VY*r;
         y2[4]=phi,y2[5]=t,y2[6]=x,y2[7]=vx,y2[8]=y,y2[9]=vy;
         y2[10]=E2,y2[11]=ex,y2[12]=ey,y2[13]=ez,y2[14]=Z,y2[15]=VZ*r,y2[16]=z,y2[17]=vz;
    }
}
void breakup(double m,double p[6],VecDoub &ystart)
{
    double x=p[0],vx=p[1],y=p[2],vy=p[3],z=p[4],vz=p[5];
    double R=sqrt(x*x+y*y+z*z);
	double E2=(vx*vx+vy*vy+vz*vz)/2.0-m/R;
	vx=vx*R,vy=vy*R,vz=vz*R;
	double v2=vx*vx+vy*vy+vz*vz,dRdtau=(x*vx+y*vy+z*vz)/R;
    double ex=v2*x/pow(R,2)-dRdtau/R*vx-m*x/R;
    double ey=v2*y/pow(R,2)-dRdtau/R*vy-m*y/R;
    double ez=v2*z/pow(R,2)-dRdtau/R*vz-m*z/R;
    ystart[0]=x,ystart[1]=vx,ystart[2]=y,ystart[3]=vy;
	ystart[6]=z,ystart[7]=vz,ystart[8]=ex,ystart[9]=ey,ystart[10]=ez,ystart[11]=E2;
}
void Infor_Tide(double p3[6],double p4[6],double infor[8],double M)
{
    double x3=p3[0],vx3=p3[1],y3=p3[2],vy3=p3[3],z3=p3[4],vz3=p3[5];
    double x4=p4[0],vx4=p4[1],y4=p4[2],vy4=p4[3],z4=p4[4],vz4=p4[5];
    double E3=(vx3*vx3+vy3*vy3+vz3*vz3)/2.0-M/sqrt(x3*x3+y3*y3+z3*z3);
    double E4=(vx4*vx4+vy4*vy4+vz4*vz4)/2.0-M/sqrt(x4*x4+y4*y4+z4*z4);
    double kp[6];
    if(E3<0)
    {
        Kepler_3D(M,p3,kp);
        infor[0]=0,infor[1]=kp[0],infor[2]=kp[1],infor[3]=kp[2];
    }
    else
    {
        double I=acos(z3/sqrt(x3*x3+y3*y3+z3*z3));
        double theta=acos(x3/sqrt(x3*x3+y3*y3));
        infor[0]=1,infor[1]=sqrt(vx3*vx3+vy3*vy3+vz3*vz3),infor[2]=I,infor[3]=theta;
    }
    if(E4<0)
    {
        Kepler_3D(M,p4,kp);
        infor[4]=0,infor[5]=kp[0],infor[6]=kp[1],infor[7]=kp[2];
    }
    else
    {
        double I=acos(z4/sqrt(x4*x4+y4*y4+z4*z4));
        double theta=acos(x4/sqrt(x4*x4+y4*y4));
        infor[4]=1,infor[5]=sqrt(vx4*vx4+vy4*vy4+vz4*vz4),infor[6]=I,infor[7]=theta;
    }
}
void TDE_output(int ending,double ai,double af, double tau,VecDoub m,int removed_particle,double tde_par[3],VecDoub ystart,VecDoub &ystart1)
{
    double m1=m[0],m2=m[1],m3=m[2],m4=m[3];
    double t=ystart[5],X=ystart[0],Y=ystart[2],Z=ystart[10],x=ystart[6],y=ystart[8],z=ystart[12];
    double VX=ystart[1],VY=ystart[3],VZ=ystart[11],vx=ystart[7],vy=ystart[9],vz=ystart[13],phi=ystart[4];
    double a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t/tau);
    double x12=a*cos(phi),y12=a*sin(phi),vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
    double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
    double x2=m1/(m1+m2)*x12,y2=m1/(m1+m2)*y12;
    double x3=X-m4/(m3+m4)*x,y3=Y-m4/(m3+m4)*y,z3=Z-m4/(m3+m4)*z;
    double x4=X+m3/(m3+m4)*x,y4=Y+m3/(m3+m4)*y,z4=Z+m3/(m3+m4)*z;
    double r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
    double r23=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+z3*z3);
    double r14=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+z4*z4);
    double r24=sqrt((x4-x2)*(x4-x2)+(y4-y2)*(y4-y2)+z4*z4);
    double vx3=VX-m4/(m3+m4)*vx,vy3=VY-m4/(m3+m4)*vy,vz3=VZ-m4/(m3+m4)*vz;
    double vx4=VX+m3/(m3+m4)*vx,vy4=VY+m3/(m3+m4)*vy,vz4=VZ+m3/(m3+m4)*vz;
    if(ending==30)
    {
        removed_particle=3,tde_par[0]=t,tde_par[1]=r13,tde_par[2]=r23;
        x=x4-x1,vx=vx4-vx1,y=y4-y1,vy=vy4-vy1,z=z4,vz=vz4;
        double R=sqrt(x*x+y*y+z*z);
		double E2=(vx*vx+vy*vy+vz*vz)/2.0-m1/R;
		vx=vx*R,vy=vy*R,vz=vz*R;
		double v2=vx*vx+vy*vy+vz*vz,dRdtau=(x*vx+y*vy+z*vz)/R;
        double ex=v2*x/pow(R,2)-dRdtau/R*vx-m1*x/R;
        double ey=v2*y/pow(R,2)-dRdtau/R*vy-m1*y/R;
        double ez=v2*z/pow(R,2)-dRdtau/R*vz-m1*z/R;
        ystart1[0]=x,ystart1[1]=vx,ystart1[2]=y,ystart1[3]=vy,ystart1[4]=ystart[4],ystart1[5]=ystart[5];
		ystart1[6]=z,ystart1[7]=vz,ystart1[8]=ex,ystart1[9]=ey,ystart1[10]=ez,ystart1[11]=E2;
    }
    else
    {
        removed_particle=4,tde_par[0]=t,tde_par[1]=r14,tde_par[2]=r24;
        x=x3-x1,vx=vx3-vx1,y=y3-y1,vy=vy3-vy1,z=z3,vz=vz3;
        double R=sqrt(x*x+y*y+z*z);
		double E2=(vx*vx+vy*vy+vz*vz)/2.0-m1/R;
		vx=vx*R,vy=vy*R,vz=vz*R;
		double v2=vx*vx+vy*vy+vz*vz,dRdtau=(x*vx+y*vy+z*vz)/R;
        double ex=v2*x/pow(R,2)-dRdtau/R*vx-m1*x/R;
        double ey=v2*y/pow(R,2)-dRdtau/R*vy-m1*y/R;
        double ez=v2*z/pow(R,2)-dRdtau/R*vz-m1*z/R;
        ystart1[0]=x,ystart1[1]=vx,ystart1[2]=y,ystart1[3]=vy,ystart1[4]=ystart[4],ystart1[5]=ystart[5];
		ystart1[6]=z,ystart1[7]=vz,ystart1[8]=ex,ystart1[9]=ey,ystart1[10]=ez,ystart1[11]=E2;
    }
}
void Output_CM_resize(int number,int values_count,vector<vector<double> > &ysave)
{
    vector<vector<double> > tempmat(ysave);
    for(int i=0;i<8;i++)
        ysave[i].resize(number);
    for (int i=0; i<8; i++)
        for (int k=0; k<values_count; k++)
            ysave[i][k]=tempmat[i][k];
}
void save1_CM(double ai,double af,double tau,VecDoub m,int &values_count,Output &out,vector<vector<double> > &ysave)
{
    double m1=m[0],m2=m[1];
    if(ysave[0].size()<(values_count+out.count))
        Output_CM_resize(out.count+values_count,values_count,ysave);
    for(int i=0;i<out.count;i++)
    { 
        double X=out.ysave[0][i],VX=out.ysave[1][i],Y=out.ysave[2][i],VY=out.ysave[3][i],phi=out.ysave[4][i],t=out.ysave[5][i],Z=out.ysave[6][i],VZ=out.ysave[7][i];
        double R1=sqrt(X*X+Y*Y+Z*Z);
        double a=af+(ai-af)*exp(-t/tau),dadt=(af-ai)/tau*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
        double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
        double x1=-m2/(m1+m2)*a*cos(phi),y1=-m2/(m1+m2)*a*sin(phi),vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
        ysave[0][values_count+i]=X+x1,ysave[1][values_count+i]=VX/R1+vx1,ysave[2][values_count+i]=Y+y1,ysave[3][values_count+i]=VY/R1+vy1,ysave[4][values_count+i]=phi;
        ysave[5][values_count+i]=t,ysave[6][values_count+i]=Z,ysave[7][values_count+i]=VZ/R1;
    }
    values_count=values_count+out.count;
}
void save2_CM(double ai,double af,double tau,VecDoub m,int &values_count,Output &out,vector<vector<double> > &ysave)
{
    if(ysave[0].size()<(values_count+out.count))
        Output_CM_resize(out.count+values_count,values_count,ysave);
    double m1=m[0],m2=m[1];
    for(int i=0;i<out.count;i++)
    {
        double t=out.ysave[5][i],phi=out.ysave[4][i];
        double a=af+(ai-af)*exp(-t/tau),dadt=(af-ai)/tau*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
        double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
            double x1=-m2/(m1+m2)*a*cos(phi),y1=-m2/(m1+m2)*a*sin(phi),vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
        double x=out.ysave[6][i],y=out.ysave[8][i],z=out.ysave[16][i],r=sqrt(x*x+y*y+z*z);
        ysave[0][values_count+i]=out.ysave[0][i]+x1,ysave[1][values_count+i]=out.ysave[1][i]/r+vx1,ysave[2][values_count+i]=out.ysave[2][i]+y1,ysave[3][values_count+i]=out.ysave[3][i]/r+vy1;
        ysave[4][values_count+i]=out.ysave[4][i],ysave[5][values_count+i]=out.ysave[5][i],ysave[6][values_count+i]=out.ysave[14][i],ysave[7][values_count+i]=out.ysave[15][i]/r;
    }
    values_count=values_count+out.count;
}
void Output_resize(int number,int values_count,vector<vector<double> > &ysave)
{
    vector<vector<double> > tempmat(ysave);
    for(int i=0;i<14;i++)
        ysave[i].resize(number);
    for (int i=0; i<14; i++)
        for (int k=0; k<values_count; k++)
            ysave[i][k]=tempmat[i][k];
}
void save2(double ai,double af,double tau,VecDoub m,int &values_count,Output &out,vector<vector<double> > &ysave)
{
    if(ysave[0].size()<(values_count+out.count))
        Output_resize(out.count+values_count,values_count,ysave);
    double m1=m[0],m2=m[1];
    for(int i=0;i<out.count;i++)
    {
        double t=out.ysave[5][i],phi=out.ysave[4][i];
        double a=af+(ai-af)*exp(-t/tau),dadt=(af-ai)/tau*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
        double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
        double x1=-m2/(m1+m2)*a*cos(phi),y1=-m2/(m1+m2)*a*sin(phi),vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
        double x=out.ysave[6][i],y=out.ysave[8][i],z=out.ysave[16][i],r=sqrt(x*x+y*y+z*z);
        ysave[0][values_count+i]=out.ysave[0][i]+x1,ysave[1][values_count+i]=out.ysave[1][i]/r+vx1,ysave[2][values_count+i]=out.ysave[2][i]+y1,ysave[3][values_count+i]=out.ysave[3][i]/r+vy1;
        ysave[4][values_count+i]=out.ysave[4][i],ysave[5][values_count+i]=out.ysave[5][i],ysave[6][values_count+i]=out.ysave[6][i],ysave[7][values_count+i]=out.ysave[7][i]/r;
        ysave[8][values_count+i]=out.ysave[8][i],ysave[9][values_count+i]=out.ysave[9][i]/r,ysave[10][values_count+i]=out.ysave[14][i],ysave[11][values_count+i]=out.ysave[15][i]/r;
        ysave[12][values_count+i]=out.ysave[16][i],ysave[13][values_count+i]=out.ysave[17][i]/r;
    }
    values_count=values_count+out.count;
}
void save3_CM(double ai,double af,double tau,VecDoub m,int &values_count,Output &out,vector<vector<double> > &ysave)
{
    double m2=m[1],m12=m[0]+m[1];
    if(ysave[0].size()<(values_count+out.count))
        Output_CM_resize(out.count+values_count,values_count,ysave);
    for(int i=0;i<out.count;i++)
    { 
        double X=out.ysave[0][i],VX=out.ysave[1][i],Y=out.ysave[2][i],VY=out.ysave[3][i],phi=out.ysave[4][i],t=out.ysave[5][i],Z=out.ysave[6][i],VZ=out.ysave[7][i];
        double a=af+(ai-af)*exp(-t/tau);
        double x1=-m2/m12*a*cos(phi),y1=-m2/m12*a*sin(phi);
        double R1=sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1)+Z*Z);
        ysave[0][values_count+i]=X,ysave[1][values_count+i]=VX/R1,ysave[2][values_count+i]=Y,ysave[3][values_count+i]=VY/R1,ysave[4][values_count+i]=phi;
        ysave[5][values_count+i]=t,ysave[6][values_count+i]=Z,ysave[7][values_count+i]=VZ/R1;
    }
    values_count=values_count+out.count;
}
void save4_CM(int &values_count,Output &out,vector<vector<double> > &ysave)
{
    if(ysave[0].size()<(values_count+out.count))
        Output_CM_resize(out.count+values_count,values_count,ysave);
   for(int i=0;i<out.count;i++)
    {
        double x=out.ysave[6][i],y=out.ysave[8][i],z=out.ysave[16][i],r=sqrt(x*x+y*y+z*z);
        ysave[0][values_count+i]=out.ysave[0][i],ysave[1][values_count+i]=out.ysave[1][i]/r,ysave[2][values_count+i]=out.ysave[2][i],ysave[3][values_count+i]=out.ysave[3][i]/r;
        ysave[4][values_count+i]=out.ysave[4][i],ysave[5][values_count+i]=out.ysave[5][i],ysave[6][values_count+i]=out.ysave[14][i],ysave[7][values_count+i]=out.ysave[15][i]/r;
    }
    values_count=values_count+out.count;
}
void save4(int &values_count,Output out,vector<vector<double> > &ysave)
{
    if(ysave[0].size()<(values_count+out.count))
        Output_resize(out.count+values_count,values_count,ysave);
    
    for(int i=0;i<out.count;i++)
    {
        double x=out.ysave[6][i],y=out.ysave[8][i],z=out.ysave[16][i],r=sqrt(x*x+y*y+z*z);
        ysave[0][values_count+i]=out.ysave[0][i],ysave[1][values_count+i]=out.ysave[1][i]/r,ysave[2][values_count+i]=out.ysave[2][i],ysave[3][values_count+i]=out.ysave[3][i]/r;
        ysave[4][values_count+i]=out.ysave[4][i],ysave[5][values_count+i]=out.ysave[5][i],ysave[6][values_count+i]=out.ysave[6][i],ysave[7][values_count+i]=out.ysave[7][i]/r;
        ysave[8][values_count+i]=out.ysave[8][i],ysave[9][values_count+i]=out.ysave[9][i]/r,ysave[10][values_count+i]=out.ysave[14][i],ysave[11][values_count+i]=out.ysave[15][i]/r; 
        ysave[12][values_count+i]=out.ysave[16][i],ysave[13][values_count+i]=out.ysave[17][i]/r;
    }
    //X,VX,Y,VY,phi,t,x,vx,y,vy,Z,VZ,z,vz
    values_count=values_count+out.count;
}
double get_mean(vector<double> f,double left,double right)
{
    double sum=0;
    for(int i=left;i<right;i++)
        sum=sum+f[i];
    return sum/(right-left);
}
void find_resonance(vector<double> T,vector<double> fr,vector<double> a,vector<double> e,vector<double> I,double inf[7],double catelogy[16],double results[2],bool &res)
{   //inf[7]: release time,a0,e0,I0,af,ef,If
    int temp=int(round(T[T.size()-1]/10)+1);
    double x[temp];
    for(int i=0;i<temp+1;i++)
        x[i]=T[0]+i*10;
    int flag=0,j=0;
    double res_fr=0,res_time=0;
    vector<double> f(temp);
    double t[temp],dis_a[temp],dis_e[temp],dis_I[temp];
    for(int i=1;i<T.size();i++)
    {
        //if(abs(T[i]-x[j])<0.1 && fr[i]>1)
        if(abs(T[i]-x[j])<0.1 || (T[i]>x[j] && T[i-1]<x[j]))
        {
            f[j]=fr[i],t[j]=T[i],dis_a[j]=a[i],dis_e[j]=e[i],dis_I[j]=I[i];
            
            j++;
            
        }
    }
    int start_resonance=0,duration=7;
    double tolerence=0.05;
    for(int i=0;i<j-duration+1;i++)
    {
        double mean=get_mean(f,i,i+duration);
        flag=1;
        for(int k=i;k<i+duration;k++)
        {
            if(abs(f[k]-mean)/mean>tolerence)
               {flag=0;break;}
        }
        if(flag==1)
            {
                res_time=t[i],res_fr=mean,res=true;
                inf[1]=dis_a[i],inf[2]=dis_e[i],inf[3]=dis_I[i];
                start_resonance=i+duration-1;
                break;
             }
    }
    
    if(res)
    {
        for(int i=0;i<16;i++)
        {
            if(abs(catelogy[i]-res_fr)/catelogy[i]<tolerence)
                {results[0]=res_time,results[1]=catelogy[i];break;}
            else if(i==15)
                {results[0]=res_time,results[1]=res_fr;break;}
            else
                continue;
        }
       
        for(int k=start_resonance;k<j;k++)
            {
                if(abs(f[k]-results[1])/results[1]>tolerence)
                    {inf[0]=t[k-1],inf[4]=dis_a[k-1],inf[5]=dis_e[k-1],inf[6]=dis_I[k-1],inf[7]=0;break;}
                else
                    inf[0]=t[k],inf[4]=dis_a[k],inf[5]=dis_e[k],inf[6]=dis_I[k],inf[7]=1;
            }
    }
    
}