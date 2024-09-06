struct Output {
	Int kmax;
	Int nvar;
	Int nsave;
	bool dense;
	Int count;
	int steps;
	Doub x1,x2,xout,dxout;
	VecDoub xsave;
	MatDoub ysave;
	Doub minseparation,mindistance1,mindistance2,temp_mindistance1;
	Output() : kmax(-1),dense(false),count(0) {}
	Output(const Int nsavee) : kmax(500),nsave(nsavee),count(0),xsave(kmax) {
		dense = nsave > 0 ? true : false;
	}
	void init(const Int neqn, const Doub xlo, const Doub xhi) {
		nvar=neqn;
		if (kmax == -1) return;
		ysave.resize(nvar,kmax);
		if (dense) {
			x1=xlo;
			x2=xhi;
			xout=x1;
			dxout=(x2-x1)/nsave;
		}
	}
	void resize() {
		Int kold=kmax;
		kmax *= 2;
		VecDoub tempvec(xsave);
		xsave.resize(kmax);
		for (Int k=0; k<kold; k++)
			xsave[k]=tempvec[k];
		MatDoub tempmat(ysave);
		ysave.resize(nvar,kmax);
		for (Int i=0; i<nvar; i++)
			for (Int k=0; k<kold; k++)
				ysave[i][k]=tempmat[i][k];
	}
	void save(double h,const Doub x, VecDoub_I &y) {
		if (kmax <= 0) return;
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=y[i];
		xsave[count++]=x;
	}
};
struct phi_integ
{
    const double m12;
    const double ai;
    const double af;
    const double tau;
    phi_integ(double mm,double aii,double aff,double tauu) : 
    m12(mm),ai(aii),af(aff),tau(tauu){}
    double operator()(double t)
    {
        double a=af+(ai-af)*exp(-t/tau);
        return sqrt(m12/pow(a,3));
    }
};
template<class Stepper>
struct Odeint_one {
	static const long int MAXSTP=1e8;
	double EPS;
	int nbad;
	int nvar;
	bool dense;
	VecDoub frag;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;
    int nok;
    double x1,x2,hmin,P,R1_criterion,R2_criterion,c;
    double separation,distance1,distance2;
	Odeint_one(VecDoub &ystartt,const double xx1,const double xx2,const int real_nvarr,
		const double atol,const double rtol,const double h1,const VecDoub &fragg,
		const double hminn,Output &outt,typename Stepper::Dtype &derivss);
    int judgement1();
	int integrate1();

};
template<class Stepper>
Odeint_one<Stepper>::Odeint_one(VecDoub &ystartt, const double xx1, const double xx2,const int real_nvarr,
	const double atol, const double rtol, const double h1,const VecDoub &fragg,const double hminn,
	Output &outt,typename Stepper::Dtype &derivss) : nvar(ystartt.size()),dense(outt.dense),
	y(nvar),dydx(nvar),ystart(ystartt),out(outt),derivs(derivss),x(xx1),nok(0),nbad(0),frag(fragg),
	x1(xx1),x2(xx2),hmin(hminn),
	s(y,dydx,x,atol,real_nvarr,rtol,dense) {
    R1_criterion=frag[10];
    R2_criterion=frag[11];
    distance1=frag[12];
    distance2=frag[13];
    separation=frag[14];
    c=frag[7];
	EPS=numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}
template<class Stepper>
int Odeint_one<Stepper>::integrate1() {
	derivs(x,y,dydx);
	out.save(h,x,y);
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			{throw("x2 is too small");}
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
        if(y[5]/frag[8]>frag[19])
            {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;return -1;}
        int ending=judgement1();
		if(ending !=3)
            {out.save(s.hdid,x,y);out.steps=nstp,out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;return ending;}
        if(y[5]>8.0*frag[6])
           {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;return 2;}
        if(nstp%5==0 && y[5]<4.0*frag[6])
            out.save(s.hdid,x,y);
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
    out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;
	throw("Too many steps in routine Odeint");
}
template<class Stepper>
int Odeint_one<Stepper>::judgement1()
{
    double m1=frag[0],m2=frag[1],m12=frag[0]+frag[1],ai=frag[4],af=frag[5],tau=frag[6];
    double t=y[5],phi=y[4];
    double a=af+(ai-af)*exp(-t/tau);
    double X=y[0],Y=y[2],Z=y[6],R1=sqrt(X*X+Y*Y+Z*Z);
    double VX=y[1]/R1,VY=y[3]/R1,VZ=y[7]/R1;
    if(sqrt(VX*VX+VY*VY+VZ*VZ)>c)
        {out.save(s.hdid,x,y);cout<<s.hdid<<endl;throw("case1, excess speed of light!");}
    double x12=a*cos(phi),y12=a*sin(phi);
    double R2=sqrt((X-x12)*(X-x12)+(Y-y12)*(Y-y12)+Z*Z);
    if(distance1>R1)
        distance1=R1;
    if(distance2>R2)
        distance2=R2;
    
    if(R1>a)//criterion of changing frame
        return 1;
    else if(R1<R1_criterion || R2<R2_criterion)
        return 0;
    else
        return 3;
}
template<class Stepper>
struct Odeint_two {
	static const long int MAXSTP=2e8;
	double EPS;
	int nok;
	int nbad;
	int nvar;
	double x1,x2,hmin,P,R1_criterion,R2_criterion,c,collision_radius,TDE_radius13,TDE_radius23,TDE_radius14,TDE_radius24;
	double separation,distance1,distance2,temp_mindistance1=100;
	bool dense;
	VecDoub frag;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;
	Odeint_two(VecDoub &ystartt,const double xx1,const double xx2,const int real_nvarr,
		const double atol,const double rtol,const double h1,const VecDoub &fragg,
		const double hminn,Output &outt,typename Stepper::Dtype &derivss);
	int judgement2();
	int integrate2();

};

template<class Stepper>
Odeint_two<Stepper>::Odeint_two(VecDoub &ystartt, const double xx1, const double xx2,const int real_nvarr,
	const double atol, const double rtol, const double h1,const VecDoub &fragg,const double hminn,
	Output &outt,typename Stepper::Dtype &derivss) : nvar(ystartt.size()),
	y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),frag(fragg),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,real_nvarr,rtol,dense) {
    R1_criterion=frag[10];
    R2_criterion=frag[11];
    distance1=frag[12];
    distance2=frag[13];
    separation=frag[14];
    c=frag[7];
    collision_radius=frag[9];
    TDE_radius13=frag[15];
    TDE_radius23=frag[16];
    TDE_radius14=frag[17];
    TDE_radius24=frag[18];
    P=frag[8];
	EPS=numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}
template<class Stepper>
int Odeint_two<Stepper>::integrate2() {
	derivs(x,y,dydx);
	out.save(h,x,y);
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			{throw("x2 is too small");}
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
        if(y[5]/frag[8]>frag[19])
            {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2,out.temp_mindistance1=temp_mindistance1;return -1;}
        int ending=judgement2();
        if(ending !=5)
        {
            out.save(s.hdid,x,y);
            out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2,out.temp_mindistance1=temp_mindistance1;
            return ending;
        }
        //if(nstp%50==0 && y[5]<4.0*frag[6]) out.save(s.hdid,x,y);
        if(nstp%100==0) out.save(s.hdid,x,y);
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
    out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2,out.temp_mindistance1=temp_mindistance1;
	throw("Too many steps in routine Odeint");
}
template<class Stepper>
int Odeint_two<Stepper>::judgement2()
{
    double m1=frag[0],m3=frag[2],m4=frag[3];
    double ai=frag[4],af=frag[5],tau=frag[6];
    double t=y[5],phi=y[4],a=af+(ai-af)*exp(-t/tau);
    double x12=a*cos(phi),y12=a*sin(phi);
    double X=y[0],Y=y[2],Z=y[14];
    double R1=sqrt(X*X+Y*Y+Z*Z);
    double R2=sqrt((X-x12)*(X-x12)+(Y-y12)*(Y-y12)+Z*Z);
    if(distance1>R1)
        distance1=R1;
    if(distance2>R2)
        distance2=R2;
    if(temp_mindistance1>R1)
        temp_mindistance1=R1;
    double rx=y[6],ry=y[8],rz=y[16];
    double x13=X-m4/(m3+m4)*rx,y13=Y-m4/(m3+m4)*ry,z13=Z-m4/(m3+m4)*rz;
    double x14=X+m3/(m3+m4)*rx,y14=Y+m3/(m3+m4)*ry,z14=Z+m3/(m3+m4)*rz;
    double r13=sqrt(x13*x13+y13*y13+z13*z13);
    double r14=sqrt(x14*x14+y14*y14+z14*z14);
    double r34=sqrt(rx*rx+ry*ry+rz*rz);
    if(separation>r34)
        separation=r34;
    double r23=sqrt((x13-x12)*(x13-x12)+(y13-y12)*(y13-y12)+z13*z13);
    double r24=sqrt((x14-x12)*(x14-x12)+(y14-y12)*(y14-y12)+z14*z14);
    if(r13<TDE_radius13 || r23<TDE_radius23)
        return 30;
    if(r14<TDE_radius14 || r24<TDE_radius24)
        return 40;
    if(r34>10.0*r13 || r34>10.0*r14 || r34>10.0)
        return 4;
    if(r34<collision_radius)
        {cout<<r34<<" "<<collision_radius<<endl;return 0;}
    else if(R1>R1_criterion && R2>R2_criterion)
        return 1;
    else
        return 5;
}
template<class Stepper>
struct Odeint_three {
	static const long int MAXSTP=1e8;
	double EPS;
	int nok;
	int nbad;
	int nvar;
	double x1,x2,hmin,P,R1_criterion,R2_criterion,c;
	double separation,distance1,distance2;
	bool dense;
	VecDoub frag;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	phi_integ &f_phi;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;
	Odeint_three(VecDoub &ystartt,const double xx1,const double xx2,const int real_nvarr,
		const double atol,const double rtol,const double h1,const VecDoub &fragg,phi_integ &f,
		const double hminn,Output &outt,typename Stepper::Dtype &derivss);
    int judgement3();
	int integrate3();
	template <class T>
    Doub qromb(T &func, Doub a, Doub b, const Doub eps);
    void transform_matrix(double phi,double w,double I,double TM[3][3]);
    void Kepler_2D(double yy[6],double kp[8],double m,double t);
    void Kepler_3D(double yy[6],double kp[8],double m,double t);

};
template<class Stepper>
Odeint_three<Stepper>::Odeint_three(VecDoub &ystartt, const double xx1, const double xx2,const int real_nvarr,
	const double atol, const double rtol, const double h1,const VecDoub &fragg,phi_integ &f,const double hminn,
	Output &outt,typename Stepper::Dtype &derivss) : nvar(ystartt.size()),
	y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),frag(fragg),f_phi(f),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,real_nvarr,rtol,dense) {
    R1_criterion=frag[10];
    R2_criterion=frag[11];
    distance1=frag[12];
    distance2=frag[13];
    separation=frag[14];
    c=frag[7];
    P=frag[8];
	EPS=numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}
template<class Stepper>
int Odeint_three<Stepper>::integrate3() {
	derivs(x,y,dydx);
	out.save(h,x,y);
	int ending=3;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			{throw("x2 is too small");}
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
        if(y[5]/frag[8]>frag[19])
            {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;return -1;}
		ending=judgement3();
		if(ending !=3)
            {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;return ending;}
        if(y[5]>8.0*frag[6])
           {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;return 2;}
        if(nstp%5==0 && y[5]<4.0*frag[6])
            out.save(s.hdid,x,y);
		
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
    out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;
	throw("Too many steps in routine Odeint");
}
template<class Stepper>
int Odeint_three<Stepper>::judgement3()
{
    double m1=frag[0],m2=frag[1],m12=frag[0]+frag[1],ai=frag[4],af=frag[5],tau=frag[6];
    double t=y[5],phi=y[4];
    double a=af+(ai-af)*exp(-t/tau);
    double x12=a*cos(phi),y12=a*sin(phi);
    double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,x2=m1/(m1+m2)*x12,y2=m1/(m1+m2)*y12;
    double X=y[0],Y=y[2],Z=y[6];
    double R=sqrt(X*X+Y*Y+Z*Z);
    double R1=sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1)+Z*Z);
    double R2=sqrt((X-x2)*(X-x2)+(Y-y2)*(Y-y2)+Z*Z);
    if(distance1>R1)
        distance1=R1;
    if(distance2>R2)
        distance2=R2;
    double VX=y[1]/R1,VY=y[3]/R1,VZ=y[7]/R1;
	double dRdt=VX*X+VY*Y+VZ*Z;
    if(sqrt(VX*VX+VY*VY+VZ*VZ)>c)
        {out.save(s.hdid,x,y);cout<<s.hdid<<endl;throw("case1, excess speed of light!");}
    double E=(VX*VX+VY*VY+VZ*VZ)/2.0-(m1+m2)/R;
    if(R>60.0*a && E>0)//criterion of binary escape as a whole
        return 1;
    else if(R1<R1_criterion || R2<R2_criterion)
        return 0;
    else if(R>60.0*a && dRdt>0 )
	{
        double r[6]={X,VX,Y,VY,Z,VZ},kp[8];
        Kepler_3D(r,kp,m12,t);
		if(kp[0]>0 && kp[0]*(1-kp[1])>60.0*a)
		    return 2;
        double n=sqrt(m12/pow(kp[0],3)),e=kp[1],f1=kp[3],E1=kp[5];
        double f2=2.0*M_PI-f1;
        double half_E2=atan(sqrt((1-e)/(1+e))*tan(f2/2.0));
        if(half_E2<0)
            half_E2=half_E2+M_PI;
        double E2=2.0*half_E2;
        double deta_t=(E2-E1-e*(sin(E2)-sin(E1)))/n;
        double t2=t+deta_t,dEdt=n/(1-e*cos(E2));
        double deta_phi=qromb(f_phi,t,t2,1e-10);
        double phi2=phi+deta_phi;
        if(phi2>2.0*M_PI)
            phi2=phi2-(ceil(phi2/(2.0*M_PI))-1)*2.0*M_PI;
        a=af+(ai-af)*exp(-t2/tau);
        x12=a*cos(phi2),y12=a*sin(phi2);
        x1=-m2/m12*x12,y1=-m2/m12*y12;
        double nx=kp[0]*(cos(E2)-e),ny=kp[0]*sqrt(1-e*e)*sin(E2);
        double nvx=-kp[0]*sin(E2)*dEdt,nvy=kp[0]*sqrt(1-e*e)*cos(E2)*dEdt;
        double TM[3][3];
        transform_matrix(kp[7],kp[4],kp[2],TM);
        X=TM[0][0]*nx+TM[0][1]*ny,Y=TM[1][0]*nx+TM[1][1]*ny,Z=TM[2][0]*nx+TM[2][1]*ny;
        VX=TM[0][0]*nvx+TM[0][1]*nvy,VY=TM[1][0]*nvx+TM[1][1]*nvy,VZ=TM[2][0]*nvx+TM[2][1]*nvy;
        R1=sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1)+Z*Z);
        y[0]=X,y[1]=VX*R1,y[2]=Y,y[3]=VY*R1,y[4]=phi2,y[5]=t2,y[6]=Z,y[7]=VZ*R1;
		derivs(x,y,dydx);
        return 3;
	}
    else
	   return 3;
}
template<class Stepper>
template <class T>
Doub Odeint_three<Stepper>::qromb(T &func, Doub a, Doub b, const Doub eps)
{
	const Int JMAX=20, JMAXP=JMAX+1, K=5;
	VecDoub s(JMAX),h(JMAXP);
	Poly_interp polint(h,s,K);
	h[0]=1.0;
	Trapzd<T> t(func,a,b);
	for (Int j=1;j<=JMAX;j++) {
		s[j-1]=t.next();
		if (j >= K) {
			Doub ss=polint.rawinterp(j-K,0.0);
			if (abs(polint.dy) <= eps*abs(ss)) return ss;
		}
		h[j]=0.25*h[j-1];
	}
    throw("Too many steps in routine qromb");
}
template<class Stepper>
void Odeint_three<Stepper>::transform_matrix(double phi,double w,double I,double TM[3][3])
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
template<class Stepper>
void Odeint_three<Stepper>::Kepler_2D(double yy[6],double kp[8],double m,double t)
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
template<class Stepper>
void Odeint_three<Stepper>::Kepler_3D(double yy[6],double kp[8],double m,double t)
{
    double x=yy[0],vx=yy[1],y=yy[2],vy=yy[3],z=yy[4],vz=yy[5];
    double r=sqrt(x*x+y*y+z*z);
    double v=sqrt(vx*vx+vy*vy+vz*vz);
    double drdt=(x*vx+y*vy+z*vz)/r;
    double h[3]={y*vz-z*vy,z*vx-x*vz,x*vy-y*vx};
    double h_M=sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
    double energy=-m/r+v*v/2.0;
    if(energy>0)
        throw("test particle escapes");
    double a=-m/(2.0*energy);
    double e=sqrt(abs(1-h_M*h_M/(a*m)));
    double I=acos(h[2]/h_M);
    if(I<1e-7*M_PI || I>(1.0-1e-7)*M_PI)
    {
        Kepler_2D(yy,kp,m,t);
        return;
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
}
template<class Stepper>
struct Odeint_four {
	static const long int MAXSTP=2e8;
	double EPS;
	int nok;
	int nbad;
	int nvar;
	double x1,x2,hmin,P,R1_criterion,R2_criterion,c,collision_radius,TDE_radius13,TDE_radius23,TDE_radius14,TDE_radius24;
	double separation,distance1,distance2,temp_mindistance1=100;
	bool dense;
	VecDoub frag;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;
	Odeint_four(VecDoub &ystartt,const double xx1,const double xx2,const int real_nvarr,
		const double atol,const double rtol,const double h1,const VecDoub &fragg,
		const double hminn,Output &outt,typename Stepper::Dtype &derivss);
	int judgement4();
	int integrate4();

};

template<class Stepper>
Odeint_four<Stepper>::Odeint_four(VecDoub &ystartt, const double xx1, const double xx2,const int real_nvarr,
	const double atol, const double rtol, const double h1,const VecDoub &fragg,const double hminn,
	Output &outt,typename Stepper::Dtype &derivss) : nvar(ystartt.size()),
	y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),frag(fragg),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,real_nvarr,rtol,dense) {
    R1_criterion=frag[10];
    R2_criterion=frag[11];
    distance1=frag[12];
    distance2=frag[13];
    separation=frag[14];
    c=frag[7];
    collision_radius=frag[9];
    TDE_radius13=frag[15];
    TDE_radius23=frag[16];
    TDE_radius14=frag[17];
    TDE_radius24=frag[18];
    P=frag[8];
	EPS=numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}
template<class Stepper>
int Odeint_four<Stepper>::integrate4() {
	derivs(x,y,dydx);
	out.save(h,x,y);
	int ending=3;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			{throw("x2 is too small");}
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
        if(y[5]/frag[8]>frag[19])
            {out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2,out.temp_mindistance1=temp_mindistance1;return -1;}
		ending=judgement4();
        if(ending!=5)
        {
            out.save(s.hdid,x,y);
            out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2,out.temp_mindistance1=temp_mindistance1;
            return ending;
        }
        //if(nstp%50==0 && y[5]<4.0*frag[6]) out.save(s.hdid,x,y);
        if(nstp%100==0)   out.save(s.hdid,x,y);
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
    out.save(s.hdid,x,y);out.steps=nstp;out.minseparation=separation,out.mindistance1=distance1,out.mindistance2=distance2;
	throw("Too many steps in routine Odeint");
}
template<class Stepper>
int Odeint_four<Stepper>::judgement4()
{
    double m1=frag[0],m2=frag[1],m3=frag[2],m4=frag[3];
    double ai=frag[4],af=frag[5],tau=frag[6];
    double t=y[5],phi=y[4],a=af+(ai-af)*exp(-t/tau);
    double x12=a*cos(phi),y12=a*sin(phi);
    double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,x2=m1/(m1+m2)*x12,y2=m1/(m1+m2)*y12;
    double X=y[0],Y=y[2],Z=y[14];
    double R1=sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1)+Z*Z);
    double R2=sqrt((X-x2)*(X-x2)+(Y-y2)*(Y-y2)+Z*Z);
    double rx=y[6],ry=y[8],rz=y[16];
    double x3=X-m4/(m3+m4)*rx,y3=Y-m4/(m3+m4)*ry,z3=Z-m4/(m3+m4)*rz;
    double x4=X+m3/(m3+m4)*rx,y4=Y+m3/(m3+m4)*ry,z4=Z+m3/(m3+m4)*rz;
    double r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
    double r14=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+z4*z4);
    double r34=sqrt(rx*rx+ry*ry+rz*rz);
    if(separation>r34)
        separation=r34;
    if(distance1>R1)
        distance1=R1;
    if(distance2>R2)
        distance2=R2;
    if(temp_mindistance1>R1)
        temp_mindistance1=R1;
    double r23=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+z3*z3);
    double r24=sqrt((x4-x2)*(x4-x2)+(y4-y2)*(y4-y2)+z4*z4);
    if(r13<TDE_radius13 || r23<TDE_radius23)
        return 30;
    if(r14<TDE_radius14 || r24<TDE_radius24)
        return 40;
    if(r34>10.0*r13 || r34>10.0*r14 || r34>10.0)
        return 4;
    if(r34<collision_radius)
        {cout<<r34<<" "<<collision_radius<<endl;return 0;}
    else if(R1>R1_criterion && R2>R2_criterion)
        return 1;
    else
        return 5;
}
template<class Stepper>
struct Odeint_five {
	static const long int MAXSTP=2e8;
	double EPS;
	int nok;
	int nbad;
	int nvar;
	double x1,x2,hmin,P,c,TDE_radius1,TDE_radius2;
	bool dense;
	VecDoub frag;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp,removed_particle;
	double x,h;
	Odeint_five(VecDoub &ystartt,const double xx1,const double xx2,const int real_nvarr,
		const double atol,const double rtol,const double h1,const VecDoub &fragg,
		const double hminn,Output &outt,typename Stepper::Dtype &derivss,int par);
    int judgement5();
	int integrate5();
};

template<class Stepper>
Odeint_five<Stepper>::Odeint_five(VecDoub &ystartt, const double xx1, const double xx2,const int real_nvarr,
	const double atol, const double rtol, const double h1,const VecDoub &fragg,const double hminn,
	Output &outt,typename Stepper::Dtype &derivss,int par) : nvar(ystartt.size()),
	y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),frag(fragg),removed_particle(par),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,real_nvarr,rtol,dense) {
    c=frag[7];
    if(removed_particle==4)
        {TDE_radius1=frag[15],TDE_radius2=frag[16];}
    else
        {TDE_radius1=frag[17],TDE_radius2=frag[18];}
    P=frag[8];
	EPS=numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}
template<class Stepper>
int Odeint_five<Stepper>::integrate5() {
	derivs(x,y,dydx);
	out.save(h,x,y);
    int ending=3;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
        ending=judgement5();
        if(ending!=3)
             {out.steps=nstp;out.save(s.hdid,x,y);return ending;}
	    //if(nstp%10==0)
        //   out.save(s.hdid,x,y);
		if (abs(s.hnext) <= hmin) 
            throw("Step size too small in Odeint");
		h=s.hnext;
	}
    out.save(s.hdid,x,y);
	throw("Too many steps in routine Odeint");
}
template<class Stepper>
int Odeint_five<Stepper>::judgement5()
{
    double t=y[5],phi=y[4],m1=frag[0],m2=frag[1],m12=frag[1]+frag[0];
    double ai=frag[4],af=frag[5],tau=frag[6];
    double a=af+(ai-af)*exp(-t/tau),x12=a*cos(phi),y12=a*sin(phi);
    double x13=y[0],y13=y[2],z13=y[6];
    double r13=sqrt(x13*x13+y13*y13+z13*z13);
	double r23=sqrt((x13-x12)*(x13-x12)+(y13-y12)*(y13-y12)+z13*z13);
    double vx13=y[1]/r13,vy13=y[3]/r13,vz13=y[7]/r13;
    
    if(sqrt(vx13*vx13+vy13*vy13+vz13*vz13)>c)
        throw("excess the speed of light!");
    if(r13>a )
        return 4;
    else if(r13<TDE_radius1 || r23<TDE_radius2)//TDE radius
        {cout<<r13<<" "<<r23<<" "<<TDE_radius1<<" "<<TDE_radius2<<endl;return 0;} 
	else
        return 3;
}
template<class Stepper>
struct Odeint_six {
	static const long int MAXSTP=2e8;
	double EPS;
	int nok;
	int nbad;
	int nvar;
	double x1,x2,hmin,P,c,TDE_radius1,TDE_radius2;
	bool dense;
	VecDoub frag;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	phi_integ &f_phi;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp,removed_particle;
	double x,h;
	Odeint_six(VecDoub &ystartt,const double xx1,const double xx2,const int real_nvarr,
		const double atol,const double rtol,const double h1,const VecDoub &fragg,phi_integ &f,
		const double hminn,Output &outt,typename Stepper::Dtype &derivss,int par);
    int judgement6();
	int integrate6();
	template <class T>
    Doub qromb(T &func, Doub a, Doub b, const Doub eps);
    void transform_matrix(double phi,double w,double I,double TM[3][3]);
    void Kepler_2D(double yy[6],double kp[8],double m,double t);
    void Kepler_3D(double yy[6],double kp[8],double m,double t);
    bool escape=false;
};

template<class Stepper>
Odeint_six<Stepper>::Odeint_six(VecDoub &ystartt, const double xx1, const double xx2,const int real_nvarr,
	const double atol, const double rtol, const double h1,const VecDoub &fragg,phi_integ &f,const double hminn,
	Output &outt,typename Stepper::Dtype &derivss,int par) : nvar(ystartt.size()),
	y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),frag(fragg),removed_particle(par),f_phi(f),
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,real_nvarr,rtol,dense) {
    c=frag[7];
    if(removed_particle==4)
        {TDE_radius1=frag[15],TDE_radius2=frag[16];}
    else
        {TDE_radius1=frag[17],TDE_radius2=frag[18];}
    P=frag[8];
	EPS=numeric_limits<double>::epsilon();
	h=SIGN(h1,x2-x1);
	for (int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}
template<class Stepper>
int Odeint_six<Stepper>::integrate6() {
	derivs(x,y,dydx);
	out.save(h,x,y);
   
    int ending=3;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
        //if(!escape)
        ending=judgement6();
        if(ending!=3)
             {out.steps=nstp;out.save(s.hdid,x,y);return ending;}
        if(y[5]>8.0*frag[6])
            {out.steps=nstp;out.save(s.hdid,x,y);return 2;}
	    //if(nstp%10==0)
        //   out.save(s.hdid,x,y);
		if (abs(s.hnext) <= hmin) 
            throw("Step size too small in Odeint");
		h=s.hnext;
	}
    out.save(s.hdid,x,y);
	throw("Too many steps in routine Odeint");
}
template<class Stepper>
int Odeint_six<Stepper>::judgement6()
{
    double t=y[5],phi=y[4],m1=frag[0],m2=frag[1],m12=frag[1]+frag[0];
    double ai=frag[4],af=frag[5],tau=frag[6];
    double a=af+(ai-af)*exp(-t/tau);
    double x3=y[0],y3=y[2],z3=y[6];
    double r3=sqrt(x3*x3+y3*y3+z3*z3);
    double x12=a*cos(phi),y12=a*sin(phi);
    double x1=-m2/m12*x12,y1=-m2/m12*y12,x2=m1/m12*x12,y2=m1/m12*y12;
    double r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
    double r23=sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+z3*z3);
    double vx3=y[1]/r13,vy3=y[3]/r13,vz3=y[7]/r13;
	double dr3dt=vx3*x3+vy3*y3+vz3*z3;
    if(sqrt(vx3*vx3+vy3*vy3+vz3*vz3)>c)
        throw("excess the speed of light!");
    double E=(vx3*vx3+vy3*vy3+vz3*vz3)/2.0-m12/r3;
    if(r3>60.0*a && E>0)
        //{escape=true;return 3;}
        return 1;
    else if(r13<TDE_radius1 || r23<TDE_radius2)//TDE radius
        return 0;
    else if(r3>60*a && dr3dt>0)
    {
        double r[6]={x3,vx3,y3,vy3,z3,vz3},kp[8];
        Kepler_3D(r,kp,m12,t);
        if(kp[0]>0 && kp[0]*(1-kp[1])>60.0*a)
		    return 2;
        double n=sqrt(m12/pow(kp[0],3)),e=kp[1],f1=kp[3],E1=kp[5];
        double f2=2.0*M_PI-f1;
        double half_E2=atan(sqrt((1-e)/(1+e))*tan(f2/2.0));
        if(half_E2<0)
            half_E2=half_E2+M_PI;
        double E2=2.0*half_E2;
        double deta_t=(E2-E1-e*(sin(E2)-sin(E1)))/n;
        double t2=t+deta_t,dEdt=n/(1-e*cos(E2));
        double deta_phi=qromb(f_phi,t,t2,1e-10);
        double phi2=phi+deta_phi;
        if(phi2>2.0*M_PI)
            phi2=phi2-(ceil(phi2/(2.0*M_PI))-1)*2.0*M_PI;
        a=af+(ai-af)*exp(-t2/tau);
        x12=a*cos(phi2),y12=a*sin(phi2);
        x1=-m2/m12*x12,y1=-m2/m12*y12;
        double nx=kp[0]*(cos(E2)-e),ny=kp[0]*sqrt(1-e*e)*sin(E2);
        double nvx=-kp[0]*sin(E2)*dEdt,nvy=kp[0]*sqrt(1-e*e)*cos(E2)*dEdt;
        double TM[3][3];
        transform_matrix(kp[7],kp[4],kp[2],TM);
        x3=TM[0][0]*nx+TM[0][1]*ny,y3=TM[1][0]*nx+TM[1][1]*ny,z3=TM[2][0]*nx+TM[2][1]*ny;
        vx3=TM[0][0]*nvx+TM[0][1]*nvy,vy3=TM[1][0]*nvx+TM[1][1]*nvy,vz3=TM[2][0]*nvx+TM[2][1]*nvy;
        r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
        y[0]=x3,y[1]=vx3*r13,y[2]=y3,y[3]=vy3*r13,y[4]=phi2,y[5]=t2,y[6]=z3,y[7]=vz3*r13;
		derivs(x,y,dydx);
        return 3;
    }
	else
        return 3;
}
template<class Stepper>
template <class T>
Doub Odeint_six<Stepper>::qromb(T &func, Doub a, Doub b, const Doub eps)
{
	const Int JMAX=20, JMAXP=JMAX+1, K=5;
	VecDoub s(JMAX),h(JMAXP);
	Poly_interp polint(h,s,K);
	h[0]=1.0;
	Trapzd<T> t(func,a,b);
	for (Int j=1;j<=JMAX;j++) {
		s[j-1]=t.next();
		if (j >= K) {
			Doub ss=polint.rawinterp(j-K,0.0);
			if (abs(polint.dy) <= eps*abs(ss)) return ss;
		}
		h[j]=0.25*h[j-1];
	}
    throw("Too many steps in routine qromb");
}
template<class Stepper>
void Odeint_six<Stepper>::transform_matrix(double phi,double w,double I,double TM[3][3])
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
template<class Stepper>
void Odeint_six<Stepper>::Kepler_2D(double yy[6],double kp[8],double m,double t)
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
template<class Stepper>
void Odeint_six<Stepper>::Kepler_3D(double yy[6],double kp[8],double m,double t)
{
    double x=yy[0],vx=yy[1],y=yy[2],vy=yy[3],z=yy[4],vz=yy[5];
    double r=sqrt(x*x+y*y+z*z);
    double v=sqrt(vx*vx+vy*vy+vz*vz);
    double drdt=(x*vx+y*vy+z*vz)/r;
    double h[3]={y*vz-z*vy,z*vx-x*vz,x*vy-y*vx};
    double h_M=sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
    double energy=-m/r+v*v/2.0;
    if(energy>0)
        throw("test particle escapes");
    double a=-m/(2.0*energy);
    double e=sqrt(abs(1-h_M*h_M/(a*m)));
    double I=acos(h[2]/h_M);
    if(I<1e-7*M_PI || I>(1.0-1e-7)*M_PI)
    {
        Kepler_2D(yy,kp,m,t);
        return;
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
}
