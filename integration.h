struct Base_interp
{
	Int n, mm, jsav, cor, dj;
	const Doub *xx, *yy;
	Base_interp(VecDoub_I &x, const Doub *y, Int m)
		: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
		dj = MAX(1,(int)pow((Doub)n,0.25));
	}

	Doub interp(Doub x) {
		Int jlo = cor ? hunt(x) : locate(x);
		return rawinterp(jlo,x);
	}

	Int locate(const Doub x);
	Int hunt(const Doub x);
	
	Doub virtual rawinterp(Int jlo, Doub x) = 0;

};
Int Base_interp::locate(const Doub x)
{
	Int ju,jm,jl;
	if (n < 2 || mm < 2 || mm > n) throw("locate size error");
	Bool ascnd=(xx[n-1] >= xx[0]);
	jl=0;
	ju=n-1;
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if ((x >= xx[jm]) == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}
Int Base_interp::hunt(const Doub x)
{
	Int jl=jsav, jm, ju, inc=1;
	if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
	Bool ascnd=(xx[n-1] >= xx[0]);
	if (jl < 0 || jl > n-1) {
		jl=0;
		ju=n-1;
	} else {
		if ((x >= xx[jl]) == ascnd) {
			for (;;) {
				ju = jl + inc;
				if (ju >= n-1) { ju = n-1; break;}
				else if ((x < xx[ju]) == ascnd) break;
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl;
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break;}
				else if ((x >= xx[jl]) == ascnd) break;
				else {
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if ((x >= xx[jm]) == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}
struct Poly_interp : Base_interp
{
	Doub dy;
	Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
		: Base_interp(xv,&yv[0],m), dy(0.) {}
	Doub rawinterp(Int jl, Doub x);
};
struct Quadrature{
	Int n;
	virtual Doub next() = 0;
};
template<class T>
struct Trapzd : Quadrature {
	Doub a,b,s;
	T &func;
	Trapzd() {};
	Trapzd(T &funcc, const Doub aa, const Doub bb) :
		func(funcc), a(aa), b(bb) {n=0;}
	Doub next() {
		Doub x,tnm,sum,del;
		Int it,j;
		n++;
		if (n == 1) {
			return (s=0.5*(b-a)*(func(a)+func(b)));
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}
};
Doub Poly_interp::rawinterp(Int jl, Doub x)
{
	Int i,m,ns=0;
	Doub y,den,dif,dift,ho,hp,w;
	const Doub *xa = &xx[jl], *ya = &yy[jl];
	VecDoub c(mm),d(mm);
	dif=abs(x-xa[0]);
	for (i=0;i<mm;i++) {
		if ((dift=abs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<mm;m++) {
		for (i=0;i<mm-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ((den=ho-hp) == 0.0) throw("Poly_interp error");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
	}
	return y;
}
template <class T>
Doub qromb(T &func, Doub a, Doub b, const Doub eps=1.0e-10) {
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
