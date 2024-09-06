#include "nr3.h"
#include "integration.h"
#include "odeint.h"
#include "funs.h"
#include "stepperdopr853.h"
#include "ArgumentParser.h"

const double G=6.67e-8;
const double M_sun=1.99e33;
const double pc=3.09e18;
const double yr=3.16e7;

void set_initial_condition(const double ai,const double af,const double tau,double binary_ini[6],double CM_ini[6],int seed[13],VecDoub frag,bool &merger,bool &SN,VecDoub &ystart,VecDoub m,VecDoub &kp,double &f0,double &T_KL)
{
    double collision=frag[9],lifetime3=frag[19];
    Ran myrandom(seed[12]);
    double m34=m[2]+m[3],m1=m[0],m2=m[1];
    double binary_period=2.0*M_PI*sqrt(binary_ini[0]*binary_ini[0]*binary_ini[0]/m34);
    double CM_period=2.0*M_PI*sqrt(CM_ini[0]*CM_ini[0]*CM_ini[0]/m1);
    double start_time=-tau*log(2.0*CM_ini[0]/ai);//orbital phase of m2
    T_KL=CM_period*CM_period/binary_period;
    if(CM_ini[0]>0.5*ai)
        start_time=0.0;
    else
    {
        CM_ini[5]=CM_ini[5]+start_time/CM_period*2.0*M_PI;
        if(CM_ini[5]>2.0*M_PI)
            CM_ini[5]=CM_ini[5]-(ceil(CM_ini[5]/(2.0*M_PI))-1)*2.0*M_PI;
        phi_integ f_phi(m1+m2,ai,af,tau);
        double deta_phi=qromb(f_phi,0,start_time,1e-10);
        f0=f0+deta_phi;
        if(f0>2.0*M_PI)
            f0=f0-(ceil(f0/(2.0*M_PI))-1)*2.0*M_PI;
        double I1=CM_ini[2],phi1=CM_ini[3],w1=CM_ini[4];
        double I2=binary_ini[2],phi2=binary_ini[3],w2=binary_ini[4],f2=binary_ini[5];
        double m=1,a2=1,e2=0.1;
        double kp[6]={a2,e2,I2,phi2,w2,f2},r0[6];
        Kepler_Sol(m,kp,r0);
        double TM[3][3],r1[6];
        transform_matrix2(phi1,w1,I1,TM);
        double x0=r0[0],vx0=r0[1],y0=r0[2],vy0=r0[3],z0=r0[4],vz0=r0[5];
        r1[0]=TM[0][0]*x0+TM[0][1]*y0+TM[0][2]*z0,r1[2]=TM[1][0]*x0+TM[1][1]*y0+TM[1][2]*z0,r1[4]=TM[2][0]*x0+TM[2][1]*y0+TM[2][2]*z0;
        r1[1]=TM[0][0]*vx0+TM[0][1]*vy0+TM[0][2]*vz0,r1[3]=TM[1][0]*vx0+TM[1][1]*vy0+TM[1][2]*vz0,r1[5]=TM[2][0]*vx0+TM[2][1]*vy0+TM[2][2]*vz0;
        Kepler_3D(m,r1,kp);
        if(abs(kp[0]-a2)>1e-10 || abs(kp[1]-e2)>1e-10)
            throw("check ZKL part!");
        double mutual_angle=kp[2]/M_PI*180.0;
        cout<<"mutual angle="<<mutual_angle<<endl;
        double I_b,w_b;
        if(T_KL>start_time || mutual_angle<39.2 || mutual_angle>140.77)
        {
            cout<<"No ZKL Oscillation"<<endl;
            binary_ini[5]=binary_ini[5]+start_time/binary_period*2.0*M_PI;
            if(binary_ini[5]>2.0*M_PI)
                binary_ini[5]=binary_ini[5]-(ceil(binary_ini[5]/(2.0*M_PI))-1)*2.0*M_PI;
        }
        else
        {
            double e_max=sqrt(1-5.0/3.0*cos(kp[2])*cos(kp[2]));
            if(binary_ini[0]*(1-e_max)<collision)
                merger=true;
            else
            {
                double e=1.0;
                while(e>e_max)
                    e=myrandom.doub(); 
                I_b=acos(cos(kp[2])/sqrt(1-e*e)); 
                if(sqrt(2.0/(5.0*sin(I_b)*sin(I_b)))>1)
                    w_b=M_PI/2.0;
                else
                    w_b=asin(sqrt(2.0/(5.0*sin(I_b)*sin(I_b))));
                if(w_b<0)
                    w_b=w_b+2.0*M_PI;
                kp[1]=e,kp[2]=I_b,kp[3]=myrandom.doub()*2.0*M_PI,kp[4]=w_b,kp[5]=myrandom.doub()*2.0*M_PI;
                Kepler_Sol(m,kp,r1);
                transform_matrix1(phi1,w1,I1,TM);
                double x1=r1[0],vx1=r1[1],y1=r1[2],vy1=r1[3],z1=r1[4],vz1=r1[5];
                r0[0]=TM[0][0]*x1+TM[0][1]*y1+TM[0][2]*z1,r0[2]=TM[1][0]*x1+TM[1][1]*y1+TM[1][2]*z1,r0[4]=TM[2][0]*x1+TM[2][1]*y1+TM[2][2]*z1;
                r0[1]=TM[0][0]*vx1+TM[0][1]*vy1+TM[0][2]*vz1,r0[3]=TM[1][0]*vx1+TM[1][1]*vy1+TM[1][2]*vz1,r0[5]=TM[2][0]*vx1+TM[2][1]*vy1+TM[2][2]*vz1;
                Kepler_3D(m,r1,kp);
                if(abs(kp[0]-a2)>1e-10)
                    throw("check ZKL part!");
                for(int k=1;k<6;k++)
                    binary_ini[k]=kp[k];
            }
        }
    }
    for(int k=1;k<6;k++)
        {binary_ini[k]=round(binary_ini[k]*1e6)/1e6,CM_ini[k]=round(CM_ini[k]*1e6)/1e6;}
    if((merger==false && start_time/frag[8]>lifetime3) || (merger==true && lifetime3<T_KL/frag[8]))
        SN=true;
    double r[6],CM[6];
    Kepler_Sol(m34,binary_ini,r);
    double x=r[0],vx=r[1],y=r[2],vy=r[3],z=r[4],vz=r[5];
    if(! Kepler_Parameters_3D(r,kp,m34,0.0) )
        cout<<"Not a bound orbit!"<<endl;
    Kepler_Sol(m1,CM_ini,CM);
    
    double X=CM[0],VX=CM[1],Y=CM[2],VY=CM[3],Z=CM[4],VZ=CM[5];//coordinates relative to m_1
    cout<<"initial conditions of the star binary: ";
    for(int i=0;i<6;i++)
        cout<<binary_ini[i]<<" ";
    cout<<endl;
    cout<<"initial conditions of the center of mass of the star binary: ";
    for(int i=0;i<6;i++)
        cout<<CM_ini[i]<<" ";
    cout<<endl;
    cout<<"T_KL="<<T_KL/frag[8]<<", start time="<<start_time/frag[8]<<endl;
    f0=round(f0*1e6)/1e6;
    cout<<setprecision(8);
    cout<<"f0="<<f0<<endl;
    ystart[0]=X,ystart[1]=VX,ystart[2]=Y,ystart[3]=VY,ystart[4]=f0,ystart[5]=start_time,ystart[6]=x,ystart[7]=vx,ystart[8]=y,ystart[9]=vy;
    ystart[10]=Z,ystart[11]=VZ,ystart[12]=z,ystart[13]=vz;
}

int main(int argc, char *argv[])
{
      ProgramOptions options = parseArguments(argc, argv);
    const int N=1;// amount of cases
    double M_SMBH=4.3e6,L=1.0; // solar mass and parsec, the mass of central body and the distance of ai
    double tidal_force=1e-5; // the criterion of neglecting tidal force of SMBH and IMBH on sterllar binary
    const bool PN=true;     
    const int orderofPN=1;  // oder of post-Newtonian approximation in simulation
  
    double m_min=1.0,m_max=100.0; // the mass range of the primary star
    double m1=4.0*M_PI*M_PI,m2=options.q_MBHB*m1; // the mass of SMBH and IMBH in simulation
    if(options.q_MBHB==-1.0 )
        m2=1e-3*m1; 
    else if(options.q_MBHB>0 && options.q_MBHB<1)
        m2=options.q_MBHB*m1;
    else throw("Wrong input of mass ratio between IMBH and SMBH!");
    const double ai=100.0,af=0.0; // a_i represents L parsec in simulation, the initial and final semi-major axis of IMBH revolving around SMBH
    double decay_time=options.tau; // define the variable 'decay_time' so than we can define a constant 'tau' later
    if(options.tau==-1.0)
        decay_time=300.0;
    else if(options.tau>0)
        decay_time=options.tau;
    else throw("Wrong input of decay time!");

    const double Period=2.0*M_PI*sqrt(ai*ai*ai/(m1+m2)),tau=decay_time*Period; // assume a=a_f+(a_i-a_f)*exp[-t/tau], tau quantifies the speed of infalling rate of IMBH
    double R_sun=6.96/3.09e8*ai/L; // The solar radius in simulation
    const double c=2.99e10*sqrt(m1/ai)/sqrt(G*M_SMBH*M_sun/(L*pc)); // G=1 in simulation, and derive the speed of light in simulation
    const double P12=2.0*M_PI*sqrt(pc*pc*pc/G/(M_SMBH+0.001*M_SMBH)/M_sun)/yr; // Initial period of IMBH revolving around SMBH
    int ini_seed=options.seed;
    if(options.seed==-1)
    {
        auto now = std::chrono::system_clock::now();
        auto now_c = std::chrono::system_clock::to_time_t(now); // use current time as initial seed
        ini_seed=now_c;
    }
    Ran myrandom(ini_seed);
    
    int seed[13];
    for(int i=0;i<14;i++)
        seed[i]=myrandom.int32();

    ofstream outfile;
    outfile.open("ending.txt");
    
    const int nvar1=12,nvar2=18,nvar3=8;
    const int real_nvar1=12,real_nvar2=18,real_nvar3=8;
    const double atol_1=1e-12,rtol_1=1e-10,atol_2=1e-16,rtol_2=1e-12,atol_3=1e-12,rtol_3=1e-10;
    const double h1=0.01,hmin=0.0,x11=0.0,x22=1e9;
    double R1_criterion=0.0,R2_criterion=0.0,binary_initial[6],CM_initial[6];
    VecDoub ystart(14),ystart1(nvar1),ystart2(nvar2),ystart3(nvar3),m(4),frag(21),kp_binary(8);
    
    
    vector<double> Mp(N),mass_ratio(N);// the mass and mass ratio of the stella binary
    pow_dis(m_min,m_max,-1.7,seed[0],Mp,0,Mp);
    pow_dis(0,1,0,seed[1],mass_ratio,1,Mp);
    m[0]=m1,m[1]=m2;
    
    //generate random direction of inital angular momentum of the small binary
    //random semi-majox axis, random initial inclination, argument of pericenter for center of mass
    vector<double> a_dis(N),I_dis(N),phi_dis(N),w_dis(N),f_dis(N);
    vector<double> a_CM_dis(N),I_CM_dis(N),phi_CM_dis(N),w_CM_dis(N),f_CM_dis(N);
    semimajor_axes(seed[2],a_dis,Mp);
    pow_dis(log(0.1),log(0.7),0,seed[3],a_CM_dis,0,Mp);
    Ran ran1(seed[4]),ran2(seed[5]),ran3(seed[6]),ran4(seed[7]),ran5(seed[8]),ran6(seed[9]),ran7(seed[10]);
    double angle_mean=0.0,angle_dispersion=3;
    Normaldev norran_I(angle_mean,angle_dispersion,seed[11]);
    Normaldev norran_phi(0,3.5,seed[13]); //The height of the stellar disk
    Ran ran8(seed[13]);
    for(int i=0;i<N;i++)
    {
        I_dis[i]=round(acos(2.0*ran1.doub()-1)*1e6)/1e6;
        phi_dis[i]=round(ran2.doub()*1e6)/1e6;
        w_dis[i]=round(ran3.doub()*1e6)/1e6;
        f_dis[i]=round(ran4.doub()*1e6)/1e6;
        I_CM_dis[i]=round(norran_I.dev()*1e6)/1e6;// The unit is degree rather than radian
        while(I_CM_dis[i]<0)
            I_CM_dis[i]=round(norran_I.dev()*1e6)/1e6;
        phi_CM_dis[i]=round(ran5.doub()*360*1e6)/1e6;
        w_CM_dis[i]=round(ran6.doub()*1e6)/1e6;
        f_CM_dis[i]=round(ran7.doub()*1e6)/1e6;
    }
    double r_binary[6],Period_binary;
    phi_integ f_phi(m1+m2,ai,af,tau);
    rhs1 d1(m,ai,af,tau,c,PN,orderofPN);
    rhs3 d3(m,ai,af,tau,c,PN,orderofPN);
    
    if(options.a_b ==-1.0)
        ;
    else if (options.a_b>100 && options.a_b<2000)
        a_dis[0]=options.a_b;
    else
        throw("The initial semi-major axis of the stellar binary is out of range!");
    if(options.a_CM ==-1.0)
        ;
    else if(options.a_CM>0.1 && options.a_CM<0.7)
        a_CM_dis[0]=log(options.a_CM);
    else
        throw("The initial semi-major axis of the center of mass is out of range!");
    if(options.I_CM ==-1.0)
        ;
    else if(options.I_CM>0 && options.I_CM<90)
        I_CM_dis[0]=options.I_CM;  
    else
        throw("The initial inclination of the center of mass is out of range!");
    
    for(int st=0;st<N;st++)
    {
        bool break_up=false,merger=false,binary_escape=false,binary_trap=false,TDE=false,SN=false;
        double m3=Mp[st]/M_SMBH*m1,m4=m3*mass_ratio[st]; // the mass of the primary star and the mass of the secondary star
        if(options.Mp ==-1.0)
            ;
        else if( options.Mp>1 && options.Mp<100)
            m3=options.Mp/M_SMBH*m1; 
        else
            {throw("The mass of the primary star is out of range!");}
        if(options.Ms ==-1.0)
            ;
        else if( options.Ms>0.1 && options.Ms<100)
            m4=options.Ms/M_SMBH*m1; 
        else
            {throw("The mass of the secondary star is out of range!");}
        double m34=m3+m4;
        double lifetime3=(2500+670*pow(Mp[st],2.5)+pow(Mp[st],4.5))/(0.033*pow(Mp[st],1.5)+0.35*pow(Mp[st],4.5))*1e6/P12; 
        double lifetime4=(2500+670*pow(Mp[st]*mass_ratio[st],2.5)+pow(Mp[st]*mass_ratio[st],4.5))/(0.033*pow(Mp[st]*mass_ratio[st],1.5)+0.35*pow(Mp[st]*mass_ratio[st],4.5))*1e6/P12;
         
        
        outfile<<"Initial seed: "<<ini_seed<<endl;
        outfile <<scientific << std::setprecision(2);
        outfile<<scientific<<"M_IMBH/M_SMBH: "<<m2/m1<<endl;
        outfile << fixed << setprecision(2);
        outfile<<"Mp: "<<Mp[0]<<" solar mass"<<endl;
        outfile<<"Ms: "<<Mp[0]*mass_ratio[0]<<" solar mass"<<endl;
        outfile<<"a_CM: "<<exp(a_CM_dis[0])<<" a_{MBHM,i}"<<endl;
        outfile<<"a_b: "<<a_dis[0]<<" R_sun"<<endl;
        outfile<<"I_CM: "<<I_CM_dis[0]<<" degree"<<endl;
        outfile<<"tau: "<<decay_time<<" P_{MBHB,i}"<<endl;
        outfile<<endl;
        
        m[2]=m3,m[3]=m4;
        double TDE_radius13=pow(Mp[st],0.6)*pow(m1/m3,1.0/3.0)*R_sun;
        double TDE_radius23=pow(Mp[st],0.6)*pow(m2/m3,1.0/3.0)*R_sun;
        double TDE_radius14=0.0,TDE_radius24=0.0,collision_radius=pow(Mp[st],0.6)*R_sun;
        if(Mp[st]*mass_ratio[st]<1.0)
        {
            TDE_radius14=pow(Mp[st]*mass_ratio[st],0.8)*pow(m1/m4,1.0/3.0)*R_sun;
            TDE_radius24=pow(Mp[st]*mass_ratio[st],0.8)*pow(m2/m4,1.0/3.0)*R_sun;
            collision_radius=collision_radius+pow(Mp[st]*mass_ratio[st],0.8)*R_sun;
        }
        else
        {
            TDE_radius14=pow(Mp[st]*mass_ratio[st],0.6)*pow(m1/m4,1.0/3.0)*R_sun;
            TDE_radius24=pow(Mp[st]*mass_ratio[st],0.6)*pow(m2/m4,1.0/3.0)*R_sun;
            collision_radius=collision_radius+pow(Mp[st]*mass_ratio[st],0.6)*R_sun;
        }
        frag[0]=m1,frag[1]=m2,frag[2]=m3,frag[3]=m4,frag[4]=ai,frag[5]=af,frag[6]=tau;
        frag[7]=c,frag[8]=Period,frag[9]=collision_radius,frag[15]=TDE_radius13,frag[16]=TDE_radius23,frag[17]=TDE_radius14,frag[18]=TDE_radius24;
        frag[19]=lifetime3,frag[20]=lifetime4;
        //set initial conditions for the stellar binary and its mass center
        binary_initial[0]=a_dis[st]*R_sun,binary_initial[1]=0.0,binary_initial[2]=I_dis[st],binary_initial[3]=phi_dis[st]*2.0*M_PI;
        binary_initial[4]=w_dis[st]*2.0*M_PI,binary_initial[5]=f_dis[st]*2.0*M_PI;
        CM_initial[0]=exp(a_CM_dis[st])*ai;
        CM_initial[3]=phi_CM_dis[st]/180.0*M_PI;CM_initial[2]=I_CM_dis[st]/180.0*M_PI;
        CM_initial[1]=0.0,CM_initial[4]=w_CM_dis[st]*2.0*M_PI,CM_initial[5]=f_CM_dis[st]*2.0*M_PI;
        seed[12]=myrandom.int32();
        double f0=myrandom.doub()*2.0*M_PI;//initial true anomaly of MBHB
        double T_KL=0.0;
        set_initial_condition(ai,af,tau,binary_initial,CM_initial,seed,frag,merger,SN,ystart,m,kp_binary,f0,T_KL);
        
        if(SN)
        {
            outfile<<"The primary star (m_3) evolves off the main sequence at time t_sim/tau="<<lifetime3*Period/tau<<endl;
            continue;
        }
        double minseparation=binary_initial[0],mindistance1=ai,mindistance2=ai;
        R1_criterion=pow(8*kp_binary[0]*kp_binary[0]*kp_binary[0]*m1/(m3+m4)/tidal_force,1.0/3.0);
        R2_criterion=pow(8*kp_binary[0]*kp_binary[0]*kp_binary[0]*m2/(m3+m4)/tidal_force,1.0/3.0);
        cout<<" R1_criterion="<<R1_criterion<<" R2_criterion="<<R2_criterion<<endl;
        frag[10]=R1_criterion,frag[11]=R2_criterion,frag[12]=mindistance1,frag[13]=mindistance2,frag[14]=minseparation;
        if(merger)
        {
            outfile<<"The two stars collide during  Zeipel-Lidov-Kozai oscillation at time t_sim/tau="<<T_KL/tau<<endl;
            continue;
        }
        
        Fundcd funcd;
        int total_steps=0;
        double t=0.0,R1=0.0,R2=0.0;
        int binary_state=3;
        bool stage2=false;
        rhs2 d2(m,ai,af,tau,c,PN,orderofPN);
        rhs4 d4(m,ai,af,tau,c,PN,orderofPN);
        //First, calculate in the frame centered on m_1
        int values_count=0;
        vector<vector<double> > My_Output;
        My_Output.resize(8);
        for(int i=0;i<8;i++)
            My_Output[i].resize(5000);
        double a_CM_last=0,e_CM_last=0,I_CM_last=0,a_last=0,e_last=0,I_last=0,R1_last=0;
        for(int i=0;i<500000;i++)
        {
            double X=ystart[0],Y=ystart[2],Z=ystart[10];
            R1=sqrt(X*X+Y*Y+Z*Z),t=ystart[5];
            double phi=ystart[4],a=af+(ai-af)*exp(-t/tau);
            double x12=a*cos(phi),y12=a*sin(phi);
            R2=sqrt((X-x12)*(X-x12)+(Y-y12)*(Y-y12)+Z*Z);
            if(R1>R1_criterion && R2>R2_criterion)
            {
                t=ystart[5];
                r_binary[0]=ystart[6],r_binary[1]=ystart[7],r_binary[2]=ystart[8],r_binary[3]=ystart[9],r_binary[4]=ystart[12],r_binary[5]=ystart[13];
                Kepler_Parameters_3D(r_binary,kp_binary,m34,t);
                frag[10]=R1_criterion,frag[11]=R2_criterion;
                frag[12]=mindistance1,frag[13]=mindistance2,frag[14]=minseparation;
                Period_binary=2.0*M_PI*sqrt(kp_binary[0]*kp_binary[0]*kp_binary[0]/m34);
                var_transform_smbh(1,ystart,ystart1,m);
                Output out1(-1);
                Odeint_one<StepperDopr853<rhs1> >  ode1(ystart1,x11,x22,real_nvar1,atol_1,rtol_1,h1,frag,hmin,out1,d1);
                try
                {
                    binary_state=ode1.integrate1();
                }catch(const char* msg)
                {
                    cout<<msg<<endl;
                    binary_state=3;
                    save1_CM(ai,af,tau,m,values_count,out1,My_Output);
                    cout<<"case 1"<<" i="<<i<<" "<<out1.ysave[5][out1.count-1]<<endl;
                    break;
                }
                save1_CM(ai,af,tau,m,values_count,out1,My_Output);
                mindistance1=out1.mindistance1,mindistance2=out1.mindistance2;
                total_steps=total_steps+out1.steps; 
                int n=out1.count-1;
                double X=out1.ysave[0][n],Y=out1.ysave[2][n],Z=out1.ysave[6][n],R=sqrt(X*X+Y*Y+Z*Z);
                double VX=out1.ysave[1][n]/R,VY=out1.ysave[3][n]/R,VZ=out1.ysave[7][n]/R,phi=out1.ysave[4][n];
                t=out1.ysave[5][n]; 
                ystart[0]=X,ystart[1]=VX,ystart[2]=Y,ystart[3]=VY,ystart[4]=phi,ystart[5]=t,ystart[10]=Z,ystart[11]=VZ;
                Kepler_Solver(t,Period_binary,kp_binary,funcd,r_binary);
                ystart[6]=r_binary[0],ystart[7]=r_binary[1],ystart[8]=r_binary[2],ystart[9]=r_binary[3];
                ystart[12]=r_binary[4],ystart[13]=r_binary[5];
                if(binary_state==1)
                {
                    cout<<"need to change frame "<<" i="<<i<<" "<<R<<endl;stage2=true;
                    break;  
                }
                else if(binary_state==2)
                {
                    outfile<<"The stellar binary is trapped at time t_sim/tau="<<t/tau<<endl;
                    binary_trap=true;
                    break;
                }
                else if(binary_state==-1)
                {
                    outfile<<"The primary star evolves off MS at time t_sim/tau="<<t/tau<<endl;
                    SN=true;
                    break;
                }
                else ;
            }
            else
            {
                double r_CM[6],kp[6];
                r_CM[0]=ystart[0],r_CM[1]=ystart[1],r_CM[2]=ystart[2],r_CM[3]=ystart[3],r_CM[4]=ystart[10],r_CM[5]=ystart[11];
                Kepler_3D(m1,r_CM,kp);
                a_CM_last=kp[0],e_CM_last=kp[1],I_CM_last=kp[2];
                r_binary[0]=ystart[6],r_binary[1]=ystart[7],r_binary[2]=ystart[8],r_binary[3]=ystart[9],r_binary[4]=ystart[12],r_binary[5]=ystart[13];
                Kepler_3D(m34,r_binary,kp);
                a_last=kp[0],e_last=kp[1],I_last=kp[2];

                frag[10]=R1_criterion,frag[11]=R2_criterion;
                frag[12]=mindistance1,frag[13]=mindistance2,frag[14]=minseparation;
                var_transform_smbh(2,ystart,ystart2,m);
                Output out2(-1);
                Odeint_two<StepperDopr853<rhs2> >  ode2(ystart2,x11,x22,real_nvar2,atol_2,rtol_2,h1,frag,hmin,out2,d2);
                try
                {
                    binary_state=ode2.integrate2();
                }catch(const char* msg)
                {
                    cout<<msg<<endl;
                    binary_state=3;
                    save2_CM(ai,af,tau,m,values_count,out2,My_Output);
                    cout<<"case 2"<<" i="<<i<<" "<<out2.ysave[5][out2.count-1]/Period<<endl;
                    break;
                }
                R1_last=out2.temp_mindistance1;
                save2_CM(ai,af,tau,m,values_count,out2,My_Output);
                total_steps=total_steps+out2.steps;
                minseparation=out2.minseparation,mindistance1=out2.mindistance1,mindistance2=out2.mindistance2;
                int n=out2.count-1;
                double t=out2.ysave[5][n],phi=out2.ysave[4][n];
                double X=out2.ysave[0][n],Y=out2.ysave[2][n],Z=out2.ysave[14][n];
                double x=out2.ysave[6][n],y=out2.ysave[8][n],z=out2.ysave[16][n];
                double r=sqrt(x*x+y*y+z*z);
                double VX=out2.ysave[1][n]/r,VY=out2.ysave[3][n]/r,VZ=out2.ysave[15][n]/r;
                double vx=out2.ysave[7][n]/r,vy=out2.ysave[9][n]/r,vz=out2.ysave[17][n]/r;
                ystart[0]=X,ystart[1]=VX,ystart[2]=Y,ystart[3]=VY,ystart[4]=phi,ystart[5]=t;
                ystart[6]=x,ystart[7]=vx,ystart[8]=y,ystart[9]=vy,ystart[10]=Z,ystart[11]=VZ,ystart[12]=z,ystart[13]=vz;
                if(binary_state==-1)
                {
                    outfile<<"The primary star evolves off MS at time t_sim/tau="<<t/tau<<endl;
                    SN=true;
                    break;
                }
                else if(binary_state==30)
                {
                    outfile<<"The primary star is tidal disrupted at time t_sim/tau="<<t/tau<<endl;
                    TDE=true;
                    break;
                }
                else if(binary_state==40)
                {
                    outfile<<"The secondary star is tidal disrupted at time t_sim/tau="<<t/tau<<endl;
                    TDE=true;
                    break;
                }
                else if(binary_state==4)
                {
                    outfile<<"The stellar binary is broken up at time: t_sim/tau="<<t/tau<<endl;
                    break_up=true;
                    break;
                }
                else if(binary_state==0)
                {
                    outfile<<"The two stars collide at time t_sim/tau="<<t/tau<<endl;
                    merger=true;
                    break;
                }
                else ;            
                r_binary[0]=x,r_binary[1]=vx,r_binary[2]=y,r_binary[3]=vy,r_binary[4]=z,r_binary[5]=vz;
                double E_int=(vx*vx+vy*vy+vz*vz)/2.0-m34/sqrt(x*x+y*y+z*z);
                if(E_int>0.0)
                    {outfile<<"The stellar binary is broken up at time t_sim/tau="<<t/tau<<endl;
                    break_up=true;break;}
                double a_bin=-m34/(2.0*E_int);
                double temp1=pow(8*a_bin*a_bin*a_bin*m1/(m3+m4)/tidal_force,1.0/3.0);
                double temp2=pow(8*a_bin*a_bin*a_bin*m2/(m3+m4)/tidal_force,1.0/3.0);
                if(abs(temp1-R1_criterion)>0.01*R1_criterion)
                    R1_criterion=temp1;
                if(abs(temp2-R2_criterion)>0.01*R2_criterion)
                    R2_criterion=temp2;
            }
        }

        

        frame_transform(ystart,m,ai,af,tau);
        if(stage2)
        {
            for(int i=0;i<500000;i++)
            {
                double X=ystart[0],Y=ystart[2],Z=ystart[10];
                t=ystart[5];double phi=ystart[4],a=af+(ai-af)*exp(-t/tau);
                double x12=a*cos(phi),y12=a*sin(phi);
                double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,x2=m1/(m1+m2)*x12,y2=m1/(m1+m2)*y12;
                R1=sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1)+Z*Z);
                R2=sqrt((X-x2)*(X-x2)+(Y-y2)*(Y-y2)+Z*Z);
                if(R1>R1_criterion && R2>R2_criterion)
                {
                    t=ystart[5];
                    r_binary[0]=ystart[6],r_binary[1]=ystart[7],r_binary[2]=ystart[8],r_binary[3]=ystart[9],r_binary[4]=ystart[12],r_binary[5]=ystart[13];
                    Kepler_Parameters_3D(r_binary,kp_binary,m34,t);
                    frag[10]=R1_criterion,frag[11]=R2_criterion;
                    frag[12]=mindistance1,frag[13]=mindistance2,frag[14]=minseparation;
                    Period_binary=2*M_PI*sqrt(pow(kp_binary[0],3)/m34);
                    var_transform_cm(1,ystart,ystart3,m,ai,af,tau);
                    Output out3(-1);
                    
                    Odeint_three<StepperDopr853<rhs3> >  ode3(ystart3,x11,x22,real_nvar3,atol_3,rtol_3,h1,frag,f_phi,hmin,out3,d3);
                    try
                    {
                        binary_state=ode3.integrate3();
                    }catch(const char* msg)
                    {
                        cout<<msg<<endl;
                        binary_state=3;
                        save3_CM(ai,af,tau,m,values_count,out3,My_Output);
                        cout<<"case 3"<<" i="<<i<<" "<<out3.ysave[5][out3.count-1]/Period<<endl;
                        break;
                    }
                    save3_CM(ai,af,tau,m,values_count,out3,My_Output);
                    mindistance1=out3.mindistance1,mindistance2=out3.mindistance2;
                    total_steps=total_steps+out3.steps; 
                    int n=out3.count-1;
                    t=out3.ysave[5][n],phi=out3.ysave[4][n]; 
                    a=af+(ai-af)*exp(-t/tau);
                    x12=a*cos(phi),y12=a*sin(phi);
                    x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
                    X=out3.ysave[0][n],Y=out3.ysave[2][n],Z=out3.ysave[6][n];
                    R1=sqrt((X-x1)*(X-x1)+(Y-y1)*(Y-y1)+Z*Z);
                    double VX=out3.ysave[1][n]/R1,VY=out3.ysave[3][n]/R1,VZ=out3.ysave[7][n]/R1;
                    ystart[0]=X,ystart[1]=VX,ystart[2]=Y,ystart[3]=VY,ystart[4]=phi,ystart[5]=t,ystart[10]=Z,ystart[11]=VZ;
                    Kepler_Solver(t,Period_binary,kp_binary,funcd,r_binary);
                    ystart[6]=r_binary[0],ystart[7]=r_binary[1],ystart[8]=r_binary[2],ystart[9]=r_binary[3];
                    ystart[12]=r_binary[4],ystart[13]=r_binary[5];
                    if(binary_state==1)
                    {
                        outfile<<"The stellar binary runs away at time t_sim/tau= "<<t/tau<<endl;
                        binary_escape=true;
                        break;
                    }
                    else if(binary_state==2)
                    {
                        outfile<<"The stellar binary is trapped at time t_sim/tau= "<<t/tau<<endl;
                        binary_trap=true;
                        break;
                    }
                    else if(binary_state==-1)
                    {
                        outfile<<"The primary star evolves off MS at time t_sim/tau= "<<t/tau<<endl;
                        SN=true;
                        break;
                    }
                    else ;
                    
                    
                }
                else
                {
                    double dadt=(af-ai)/tau*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3));
                    double vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
                    double vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
                    double r_CM[6],kp[6];
                    r_CM[0]=ystart[0]-x1,r_CM[1]=ystart[1]-vx1,r_CM[2]=ystart[2]-y1,r_CM[3]=ystart[3]-vy1,r_CM[4]=ystart[10],r_CM[5]=ystart[11];
                    Kepler_3D(m1,r_CM,kp);
                    a_CM_last=kp[0],e_CM_last=kp[1],I_CM_last=kp[2];
                    r_binary[0]=ystart[6],r_binary[1]=ystart[7],r_binary[2]=ystart[8],r_binary[3]=ystart[9],r_binary[4]=ystart[12],r_binary[5]=ystart[13];
                    Kepler_3D(m34,r_binary,kp);
                    a_last=kp[0],e_last=kp[1],I_last=kp[2];

                    frag[10]=R1_criterion,frag[11]=R2_criterion;
                    frag[12]=mindistance1,frag[13]=mindistance2,frag[14]=minseparation;
                    var_transform_cm(2,ystart,ystart2,m,ai,af,tau);
                    Output out4(-1);
                    rhs4 d4(m,ai,af,tau,c,PN,orderofPN);
                    Odeint_four<StepperDopr853<rhs4> >  ode4(ystart2,x11,x22,real_nvar2,atol_2,rtol_2,h1,frag,hmin,out4,d4);
                    try
                    {
                        binary_state=ode4.integrate4();
                    }catch(const char* msg)
                    {
                        cout<<msg<<endl;
                        binary_state=3;
                        save4_CM(values_count,out4,My_Output);
                        cout<<"case 4"<<" i="<<i<<" "<<out4.ysave[5][out4.count-1]/Period<<endl;
                        break;
                    }
                    R1_last=out4.temp_mindistance1;
                    save4_CM(values_count,out4,My_Output);
                    total_steps=total_steps+out4.steps;
                    minseparation=out4.minseparation,mindistance1=out4.mindistance1,mindistance2=out4.mindistance2;
                    int n=out4.count-1;
                    t=out4.ysave[5][n],phi=out4.ysave[4][n];
                    X=out4.ysave[0][n],Y=out4.ysave[2][n],Z=out4.ysave[14][n];
                    double x=out4.ysave[6][n],y=out4.ysave[8][n],z=out4.ysave[16][n];
                    double r=sqrt(x*x+y*y+z*z);
                    double VX=out4.ysave[1][n]/r,VY=out4.ysave[3][n]/r,VZ=out4.ysave[15][n]/r;
                    double vx=out4.ysave[7][n]/r,vy=out4.ysave[9][n]/r,vz=out4.ysave[17][n]/r;
                    ystart[0]=X,ystart[1]=VX,ystart[2]=Y,ystart[3]=VY,ystart[4]=phi,ystart[5]=t;
                    ystart[6]=x,ystart[7]=vx,ystart[8]=y,ystart[9]=vy,ystart[10]=Z,ystart[11]=VZ,ystart[12]=z,ystart[13]=vz;
                    if(binary_state==-1)
                    {
                        outfile<<"The primary star evolves off MS at time t_sim/tau= "<<t/tau<<endl;
                        SN=true;
                        break;
                    }
                    else if(binary_state==30)
                    {
                        outfile<<"The primary star is tidal disrupted at time t_sim/tau= "<<t/tau<<endl;
                        TDE=true;
                        break;
                    }
                    else if(binary_state==40)
                    {
                        outfile<<"The secondary star is tidal disrupted at time t_sim/tau= "<<t/tau<<endl;
                        TDE=true;
                        break;
                    }
                    else if(binary_state==4)
                    {
                        outfile<<"The stellar binary is broken up at time t_sim/tau= "<<t/tau<<endl;
                        break_up=true;
                        break;
                    }
                    else if(binary_state==0)
                    {
                        outfile<<"The two stars collide at time t_sim/tau= "<<t/tau<<endl;
                        merger=true;
                        break;
                    }
                    else ;
                    r_binary[0]=x,r_binary[1]=vx,r_binary[2]=y,r_binary[3]=vy,r_binary[4]=z,r_binary[5]=vz;
                    double E_int=(vx*vx+vy*vy+vz*vz)/2.0-m34/r;
                    if(E_int>0.0)
                        {outfile<<"The stellar binary is broken up at time t_sim/tau= "<<t/tau<<endl;break_up=true;break;}
                    double a_bin=-m34/(2.0*E_int);
                    double temp1=pow(8*a_bin*a_bin*a_bin*m1/(m3+m4)/tidal_force,1.0/3.0);
                    double temp2=pow(8*a_bin*a_bin*a_bin*m2/(m3+m4)/tidal_force,1.0/3.0);
                    if(abs(temp1-R1_criterion)>0.01*R1_criterion)
                        R1_criterion=temp1;
                    if(abs(temp2-R2_criterion)>0.01*R2_criterion)
                        R2_criterion=temp2;
                }
            }
        }  
        
        ofstream NumPyFile;
        NumPyFile.open("data.npy", ios::binary | ios::out);
        NumPyFile.write("\x93NUMPY\x01\x00\x76\x00", 10);
        NumPyFile << "{'descr': '<f8', 'fortran_order': False, 'shape': (                                                                  \n";
        int rows=0,columns=6;
        vector<double> myoutput(columns);
        for(int i=0;i<values_count;i++)
            {
                double X=My_Output[0][i],Y=My_Output[2][i],Z=My_Output[6][i];
                double phi=My_Output[4][i],t=My_Output[5][i];
                double a=af+(ai-af)*exp(-t/tau);
                double x2=m1/(m1+m2)*a*cos(phi),y2=m1/(m1+m2)*a*sin(phi);
                myoutput[0]=t/Period,myoutput[1]=X,myoutput[2]=Y,myoutput[3]=Z;
                myoutput[4]=x2,myoutput[5]=y2;
                NumPyFile.write((char *)(myoutput.data()), sizeof(double) * columns);
                rows=rows+1;
            }
         NumPyFile.seekp(61);
         NumPyFile << rows<< ", " << columns << "), }";
         NumPyFile.close();  

        int fate_of_binary=-2;
        if(SN)
            fate_of_binary=-1;
        else if (merger)
            fate_of_binary=0;
        else if(binary_escape)
            fate_of_binary=1;
        else if(binary_trap)
            fate_of_binary=2;
        else if(TDE)
            fate_of_binary=3;
        else if(break_up)
            fate_of_binary=4;
        else
            fate_of_binary=5;
        cout<<"fate_of_binary="<<fate_of_binary<<endl;
        if(break_up)
        {
            t=ystart[5];
            int particle3_ending=3,particle4_ending=3;
            double X=ystart[0],Y=ystart[2],Z=ystart[10],x=ystart[6],y=ystart[8],z=ystart[12];
            double VX=ystart[1],VY=ystart[3],VZ=ystart[11],vx=ystart[7],vy=ystart[9],vz=ystart[13],phi=ystart[4];
            double a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t/tau);
            double x12=a*cos(phi),y12=a*sin(phi),vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
            double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
            double x3=X-m4/(m3+m4)*x,y3=Y-m4/(m3+m4)*y,z3=Z-m4/(m3+m4)*z;
            double x4=X+m3/(m3+m4)*x,y4=Y+m3/(m3+m4)*y,z4=Z+m3/(m3+m4)*z;
            double r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
            double r14=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+z4*z4);
            double vx3=VX-m4/(m3+m4)*vx,vy3=VY-m4/(m3+m4)*vy,vz3=VZ-m4/(m3+m4)*vz;
            double vx4=VX+m3/(m3+m4)*vx,vy4=VY+m3/(m3+m4)*vy,vz4=VZ+m3/(m3+m4)*vz;
            double position3[6]={x3-x1,vx3-vx1,y3-y1,vy3-vy1,z3,vz3},position4[6]={x4-x1,vx4-vx1,y4-y1,vy4-vy1,z4,vz4};
            double t3=t,t4=t;
            double infor_tide[8];
            Infor_Tide(position3,position4,infor_tide,m1);
            breakup(m1,position3,ystart1);
            ystart1[4]=phi,ystart1[5]=t3;
            int removed_particle=4;
            Output out5(-1);
            //In the frame centered on m1
            Odeint_five<StepperDopr853<rhs1> >  ode5(ystart1,x11,x22,real_nvar1,atol_1,rtol_1,h1,frag,hmin,out5,d1,removed_particle);
            try
            {
                particle3_ending=ode5.integrate5();
            }catch(const char* msg)
            {
                outfile<<"There is something wrong in the simulation of the first particle!"<<endl;
                cerr<<msg<<" case5"<<endl;
            }
            cout<<out5.count<<" "<<particle3_ending<<endl;
            
            int n_data3=out5.count;
            vector<vector<double> > Inf_particle3,temp_particle3;
            Inf_particle3.resize(8),temp_particle3.resize(8);
            for(int i=0;i<8;i++)
                {Inf_particle3[i].resize(out5.count),temp_particle3[i].resize(out5.count);}
            for(int i=0;i<out5.count;i++)
            {
                t=out5.ysave[5][i],phi=out5.ysave[4][i];
                a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t/tau);
                x12=a*cos(phi),y12=a*sin(phi),vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
                x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
                double x13=out5.ysave[0][i],y13=out5.ysave[2][i],z13=out5.ysave[6][i];
                r13=sqrt(x13*x13+y13*y13+z13*z13);
                double vx13=out5.ysave[1][i]/r13,vy13=out5.ysave[3][i]/r13,vz13=out5.ysave[7][i]/r13;
                x3=x13+x1,vx3=vx13+vx1,y3=y13+y1,vy3=vy13+vy1,z3=z13,vz3=vz13;
                Inf_particle3[0][i]=t,Inf_particle3[1][i]=phi,Inf_particle3[2][i]=x3,Inf_particle3[3][i]=y3;
                Inf_particle3[4][i]=z3,Inf_particle3[5][i]=vx3,Inf_particle3[6][i]=vy3,Inf_particle3[7][i]=vz3;
                temp_particle3[0][i]=t,temp_particle3[1][i]=phi,temp_particle3[2][i]=x3,temp_particle3[3][i]=y3;
                temp_particle3[4][i]=z3,temp_particle3[5][i]=vx3,temp_particle3[6][i]=vy3,temp_particle3[7][i]=vz3;
            }

            int n=out5.count-1;
            t3=out5.ysave[5][n];
            phi=out5.ysave[4][n];
            double x13=out5.ysave[0][n],y13=out5.ysave[2][n],z13=out5.ysave[6][n];
            r13=sqrt(x13*x13+y13*y13+z13*z13);
            double vx13=out5.ysave[1][n]/r13,vy13=out5.ysave[3][n]/r13,vz13=out5.ysave[7][n]/r13;
            a=af+(ai-af)*exp(-t3/tau);
            x12=a*cos(phi),y12=a*sin(phi);
            x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
            double r23=sqrt((x13-x12)*(x13-x12)+(y13-y12)*(y13-y12)+z13*z13);
            if(particle3_ending==4)
            {
                w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t3/tau);
                vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
                vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
                ystart3[0]=x13+x1,ystart3[1]=(vx13+vx1)*r13,ystart3[2]=y13+y1,ystart3[3]=(vy13+vy1)*r13,ystart3[4]=phi,ystart3[5]=t3,ystart3[6]=z13,ystart3[7]=vz13*r13;
                Output out6(-1);
                Odeint_six<StepperDopr853<rhs3> >  ode6(ystart3,x11,x22,real_nvar3,atol_3,rtol_3,h1,frag,f_phi,hmin,out6,d3,removed_particle);//the center of mass of the BHB is origin
                try
                {
                    particle3_ending=ode6.integrate6();
                }catch(const char* msg)
                {
                    cout<<msg<<"  case6"<<endl;
                }
                cout<<out6.count<<" "<<particle3_ending<<endl;
                
                n_data3=n_data3+out6.count;
                for(int i=0;i<8;i++)
                    Inf_particle3[i].resize(out5.count+out6.count);
                for(int i=0;i<out5.count;i++)
                {
                    Inf_particle3[0][i]=temp_particle3[0][i],Inf_particle3[1][i]=temp_particle3[1][i],Inf_particle3[2][i]=temp_particle3[2][i],Inf_particle3[3][i]=temp_particle3[3][i];
                    Inf_particle3[4][i]=temp_particle3[4][i],Inf_particle3[5][i]=temp_particle3[5][i],Inf_particle3[6][i]=temp_particle3[6][i],Inf_particle3[7][i]=temp_particle3[7][i];
                }  
                for(int i=0;i<out6.count;i++)
                {
                    
                    t=out6.ysave[5][i],phi=out6.ysave[4][i];
                    x3=out6.ysave[0][i],y3=out6.ysave[2][i],z3=out6.ysave[6][i];
                    a=af+(ai-af)*exp(-t/tau);
                    x12=a*cos(phi),y12=a*sin(phi);
                    x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
                    r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
                    vx3=out6.ysave[1][i]/r13,vy3=out6.ysave[3][i]/r13,vz3=out6.ysave[7][i]/r13;
                   
                    Inf_particle3[0][out5.count+i]=t,Inf_particle3[1][out5.count+i]=phi,Inf_particle3[2][out5.count+i]=x3,Inf_particle3[3][out5.count+i]=y3;
                    Inf_particle3[4][out5.count+i]=z3,Inf_particle3[5][out5.count+i]=vx3,Inf_particle3[6][out5.count+i]=vy3,Inf_particle3[7][out5.count+i]=vz3;
                }

                n=out6.count-1;
                t3=out6.ysave[5][n],phi=out6.ysave[4][n];
                a=af+(ai-af)*exp(-t3/tau);
                x12=a*cos(phi),y12=a*sin(phi);
                x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
                x3=out6.ysave[0][n],y3=out6.ysave[2][n],z3=out6.ysave[6][n];
                r13=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+z3*z3);
                r23=sqrt((x13-x12)*(x13-x12)+(y13-y12)*(y13-y12)+z13*z13);
                double vx3=out6.ysave[1][n]/r13,vy3=out6.ysave[3][n]/r13,vz3=out6.ysave[7][n]/r13;
                position3[0]=x3,position3[1]=vx3,position3[2]=y3,position3[3]=vy3,position3[4]=z3,position3[5]=vz3;
            }
    
            if(particle3_ending !=0 && particle3_ending !=1 && particle3_ending !=2)
            {
                outfile<<"There is something wrong in the simulation of the first particle!"<<endl;
                cout<<33<<","<<particle3_ending<<","<<t3/Period<<","<<33<<","<<33<<","<<33<<","<<33<<","<<33<<","<<33<<","<<33;
                
                continue;
            }

            Output out7(-1);
            removed_particle=3;
            breakup(m1,position4,ystart1);
            Odeint_five<StepperDopr853<rhs1> >  ode7(ystart1,x11,x22,real_nvar1,atol_1,rtol_1,h1,frag,hmin,out7,d1,removed_particle);
             try
            {
                particle4_ending=ode7.integrate5();
            }catch(const char* msg)
            {
                cout<<msg<<" case 7"<<endl;
            }
            cout<<out7.count<<" "<<particle4_ending<<endl;
            
            int n_data4=out7.count;
            vector<vector<double> > Inf_particle4,temp_particle4;
            Inf_particle4.resize(8),temp_particle4.resize(8);
            for(int i=0;i<8;i++)
                {Inf_particle4[i].resize(out7.count),temp_particle4[i].resize(out7.count);}
         
            n=out7.count-1;
            t4=out7.ysave[5][n];
            phi=out7.ysave[4][n];
            double x14=out7.ysave[0][n],y14=out7.ysave[2][n],z14=out7.ysave[6][n];
            r14=sqrt(x14*x14+y14*y14+z14*z14);
            double vx14=out7.ysave[1][n]/r14,vy14=out7.ysave[3][n]/r14,vz14=out7.ysave[7][n]/r14;
            a=af+(ai-af)*exp(-t4/tau);
            x12=a*cos(phi),y12=a*sin(phi);
            x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
            
            double r24=sqrt((x14-x12)*(x14-x12)+(y14-y12)*(y14-y12)+z14*z14);
            if(particle4_ending==4)
            {
                w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t4/tau);
                vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
                vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
                ystart3[0]=x14+x1,ystart3[1]=(vx14+vx1)*r14,ystart3[2]=y14+y1,ystart3[3]=(vy14+vy1)*r14,ystart3[4]=phi,ystart3[5]=t4,ystart3[6]=z14,ystart3[7]=vz14*r14;
                Output out8(-1);
                Odeint_six<StepperDopr853<rhs3> >  ode8(ystart3,x11,x22,real_nvar3,atol_3,rtol_3,h1,frag,f_phi,hmin,out8,d3,removed_particle);//the center of mass of the BHB is origin
                try
                {
                    particle4_ending=ode8.integrate6();
                }catch(const char* msg)
                {
                    cout<<msg<<"  case8"<<endl;
                }
                cout<<out8.count<<" "<<particle4_ending<<endl;
                
                n_data4=n_data4+out8.count;
                for(int i=0;i<8;i++)
                    Inf_particle4[i].resize(out7.count+out8.count);
                for(int i=0;i<out7.count;i++)
                {
                    Inf_particle4[0][i]=temp_particle4[0][i],Inf_particle4[1][i]=temp_particle4[1][i],Inf_particle4[2][i]=temp_particle4[2][i],Inf_particle4[3][i]=temp_particle4[3][i];
                    Inf_particle4[4][i]=temp_particle4[4][i],Inf_particle4[5][i]=temp_particle4[5][i],Inf_particle4[6][i]=temp_particle4[6][i],Inf_particle4[7][i]=temp_particle4[7][i];
                }  
                for(int i=0;i<out8.count;i++)
                {
                    
                    t=out8.ysave[5][i],phi=out8.ysave[4][i];
                    x4=out8.ysave[0][i],y4=out8.ysave[2][i],z4=out8.ysave[6][i];
                    a=af+(ai-af)*exp(-t/tau),w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t/tau);
                    x12=a*cos(phi),y12=a*sin(phi);
                    x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
                    r14=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+z4*z4);
                    vx4=out8.ysave[1][i]/r14,vy4=out8.ysave[3][i]/r14,vz4=out8.ysave[7][i]/r14;
                    Inf_particle4[0][out7.count+i]=t,Inf_particle4[1][out7.count+i]=phi,Inf_particle4[2][out7.count+i]=x4,Inf_particle4[3][out7.count+i]=y4;
                    Inf_particle4[4][out7.count+i]=z4,Inf_particle4[5][out7.count+i]=vx4,Inf_particle4[6][out7.count+i]=vy4,Inf_particle4[7][out7.count+i]=vz4;
                }

                n=out8.count-1;
                t4=out8.ysave[5][n],phi=out8.ysave[4][n];
                a=af+(ai-af)*exp(-t4/tau);
                x12=a*cos(phi),y12=a*sin(phi);
                x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
                x4=out8.ysave[0][n],y4=out8.ysave[2][n],z4=out8.ysave[6][n];
                r14=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+z4*z4);
                r24=sqrt((x14-x12)*(x14-x12)+(y14-y12)*(y14-y12)+z14*z14);
                double vx4=out8.ysave[1][n]/r14,vy4=out8.ysave[3][n]/r14,vz4=out8.ysave[7][n]/r14;
                position4[0]=x4,position4[1]=vx4,position4[2]=y4,position4[3]=vy4,position4[4]=z4,position4[5]=vz4;
            }
           
            double kp3[6],kp4[6];
            if(particle3_ending==0)
            {
                if(particle4_ending==0)
                   outfile<<"The primary star is tidal disrupted at time t_sim/tau="<<t3/tau
                            <<", and the secondary star is tidal disrupted at time t_sim/tau="<<t4/tau<<endl;
                else if(particle4_ending==1)
                {
                    outfile<<"The primary star is tidal disrupted at time t_sim/tau="<<t3/tau
                            <<", and the secondary star runs away at time t_sim/tau="<<t4/tau<<endl;
                }
                else if(particle4_ending==2)
                {
                    outfile<<"The primary star is tidal disrupted at time t_sim/tau="<<t3/tau
                            <<", and the secondary star is trapped at time t_sim/tau="<<t4/tau<<endl;
                }
                else
                {
                    outfile<<"The primary star is tidal disrupted at time t_sim/tau="<<t3/tau
                            <<", and there is something wrong in the simulation of the secondary star"<<endl;
                    Kepler_3D(m1+m2,position4,kp4);
                    cout<<30<<","<<3<<","<<t3/Period<<","<<r13<<","<<r23<<","<<t4/Period<<","<<kp4[0]<<","<<kp4[1]<<","<<kp4[2]<<","<<kp4[3];
                }
            }
            else if(particle3_ending==1)
            {
                double r3=sqrt(position3[0]*position3[0]+position3[2]*position3[2]+position3[4]*position3[4]);
                Kepler_3D(m1+m2,position3,kp3);
                double E3=-(m1+m2)/(2.0*kp3[0]);
                double v3_10=sqrt(2.0*(E3+(m1+m2)/10/ai));
                
                if(particle4_ending==0)
                    outfile<<"The primary star runs away at time t_sim/tau="<<t3/tau
                            <<", and the secondary star is tidal disrupted at time t_sim/tau="<<t4/tau<<endl;
                else if(particle4_ending==1)
                {
                    outfile<<"The primary star runs away at time t_sim/tau="<<t3/tau
                            <<", and the secondary star also runs away at time t_sim/tau="<<t4/tau<<endl;
                }
                else if(particle4_ending==2)
                {
                    outfile<<"The primary star runs away at time t_sim/tau="<<t3/tau
                            <<", and the secondary star is trapped at time t_sim/tau="<<t4/tau<<endl;
                }
                else
                {
                    outfile<<"The primary star runs away at time t_sim/tau="<<t3/tau
                            <<", and there is something wrong in the simulation of the secondary star"<<endl;
                    Kepler_3D(m1+m2,position4,kp4);
                    cout<<31<<","<<t3/Period<<","<<v3_10<<","<<r3<<","<<kp3[0]<<","<<kp3[1]<<","<<kp3[2]<<","<<kp3[3]<<","<<kp3[4]<<","<<t4/Period<<","<<kp4[0]<<","<<kp4[1]<<","<<kp4[2]<<","<<kp4[3];
                }
            }
            else//particle3_ending==2
            {
                Kepler_3D(m1+m2,position3,kp3);
                if(particle4_ending==0)
                    outfile<<"The primary star is trapped at time t_sim/tau="<<t3/tau
                            <<", and the secondary star is tidal disrupted at time t_sim/tau="<<t4/tau<<endl;
                else if(particle4_ending==1)
                {
                    outfile<<"The primary star is trapped at time t_sim/tau="<<t3/tau
                            <<", and the secondary star also runs away at time t_sim/tau="<<t4/tau<<endl;
                }
                else if(particle4_ending==2)
                {
                   outfile<<"The primary star is trapped at time t_sim/tau="<<t3/tau
                            <<", and the secondary star is also trapped at time t_sim/tau="<<t4/tau<<endl;
                }
                else
                {
                    outfile<<"The primary star is trapped at time t_sim/tau="<<t3/tau
                            <<", and there is something wrong in the simulation of the secondary star"<<endl;
                    Kepler_3D(m1+m2,position4,kp4);
                    cout<<32<<","<<kp3[0]<<","<<kp3[1]<<","<<kp3[2]<<","<<kp3[3]<<","<<t4/Period<<","<<kp4[0]<<","<<kp4[1]<<","<<kp4[2]<<","<<kp4[3];
                }
            }
            
        }
        else if(TDE)
        {
            t=ystart[5];
            double tde_particle[3];
            int removed_particle=0,particle_ending=3;
            TDE_output(binary_state,ai,af,tau,m,removed_particle,tde_particle,ystart,ystart1);
            double position[6]={0,0,0,0,0,0};
            Output out(-1);
            Odeint_five<StepperDopr853<rhs1> >  ode(ystart1,x11,x22,real_nvar1,atol_1,rtol_1,h1,frag,hmin,out,d1,removed_particle);
            try
            {
                particle_ending=ode.integrate5();
            }catch(const char* msg)
            {
                outfile<<"There is something wrong in the simulation after the TDE event!"<<endl;
               cout<<msg<<" TDE1"<<endl;
            }
            int n=out.count-1;
            double t=out.ysave[5][n],phi=out.ysave[4][n];
            double x=out.ysave[0][n],y=out.ysave[2][n],z=out.ysave[6][n];
            double a=af+(ai-af)*exp(-t/tau);
            double x12=a*cos(phi),y12=a*sin(phi);
            double r1=sqrt(x*x+y*y+z*z);
            double r2=sqrt((x-x12)*(x-x12)+(y-y12)*(y-y12)+z*z);
            if(particle_ending==4)
            {
                double w=sqrt((m1+m2)/pow(a,3)),dadt=(af-ai)/tau*exp(-t/tau);
                double x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12,vx12=dadt*cos(phi)-w*a*sin(phi),vy12=dadt*sin(phi)+w*a*cos(phi);
                double vx1=-m2/(m1+m2)*vx12,vy1=-m2/(m1+m2)*vy12;
                double vx=out.ysave[1][n]/r1,vy=out.ysave[3][n]/r1,vz=out.ysave[7][n]/r1;
                ystart3[0]=x+x1,ystart3[1]=(vx+vx1)*r1,ystart3[2]=y+y1,ystart3[3]=(vy+vy1)*r1,ystart3[4]=phi,ystart3[5]=t,ystart3[6]=z,ystart3[7]=vz*r1;
                Output out(-1);
                Odeint_six<StepperDopr853<rhs3> >  ode(ystart3,x11,x22,real_nvar3,atol_3,rtol_3,h1,frag,f_phi,hmin,out,d3,removed_particle);
                try
                {
                    particle_ending=ode.integrate6();
                }catch(const char* msg)
                {
                    outfile<<"There is something wrong in the simulation after the TDE event!"<<endl;
                   cout<<msg<<" TDE2"<<endl;
                }
                n=out.count-1;
                t=out.ysave[5][n],phi=out.ysave[4][n];
                a=af+(ai-af)*exp(-t/tau);
                x12=a*cos(phi),y12=a*sin(phi);
                x1=-m2/(m1+m2)*x12,y1=-m2/(m1+m2)*y12;
                double x2=m1/(m1+m2)*x12,y2=m1/(m1+m2)*y12;;
                x=out.ysave[0][n],y=out.ysave[2][n],z=out.ysave[6][n];
                r1=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+z*z);
                r2=sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+z*z);
                vx=out.ysave[1][n]/r1,vy=out.ysave[3][n]/r1,vz=out.ysave[7][n]/r1;
                position[0]=x,position[1]=vx,position[2]=y,position[3]=vy,position[4]=z,position[5]=vz;
            }
            if(particle_ending==0)
            {
                if(binary_state==30)
                    outfile<<"The secondary star is tidal disrupted at time t_sim/tau="<<t/tau<<endl;
                else
                    outfile<<"The primary star is tidal disrupted at time t_sim/tau="<<t/tau<<endl;
            }
            else if(particle_ending==1)
            {
                if(binary_state==30)
                    outfile<<"The secondary star runs away at time t_sim/tau="<<t/tau<<endl;
                else
                    outfile<<"The primary star runs away at time t_sim/tau="<<t/tau<<endl;
            }
            else if(particle_ending==2)
            {
                if(binary_state==30)
                    outfile<<"The secondary star is trapped at time t_sim/tau="<<t/tau<<endl;
                else
                    outfile<<"The primary star is trapped at time t_sim/tau="<<t/tau<<endl;
            }
            else
            {
                outfile<<"There is something wrong in the simulation after the TDE event!"<<endl;
                cout<<300<<","<<removed_particle<<","<<tde_particle[0]/Period<<","<<tde_particle[1]<<","<<tde_particle[2]
                <<","<<t/Period<<","<<300<<","<<300<<","<<300<<","<<300;
            }
        } 
        else
        {
            t=ystart[5];
            cout<<a_dis[st]<<","<<binary_initial[1]<<","<<binary_initial[2]<<","<<binary_initial[3]
            <<","<<binary_initial[4]<<","<<binary_initial[5]<<","<<a_CM_dis[st]<<","<<CM_initial[1]
            <<","<<CM_initial[2]<<","<<CM_initial[3]<<","<<CM_initial[4]<<","<<CM_initial[5]
            <<","<<Mp[st]<<","<<mass_ratio[st]<<","<<minseparation<<","<<mindistance1<<","<<mindistance2
            <<","<<3<<","<<3<<","<<3<<","<<3<<","<<3<<","<<3<<","<<3<<","<<3<<","<<t/Period<<",";
            double X=ystart[0],Y=ystart[2],Z=ystart[10],x=ystart[6],y=ystart[8],z=ystart[12];
            double VX=ystart[1],VY=ystart[3],VZ=ystart[11],vx=ystart[7],vy=ystart[9],vz=ystart[13];
            double r[6],kp[6];
            r[0]=X,r[1]=VX,r[2]=Y,r[3]=VY,r[4]=Z,r[5]=VZ;
            Kepler_3D(m1+m2,r,kp);
            cout<<43<<","<<2<<","<<kp[0]<<","<<kp[1]<<","<<kp[2]<<","<<kp[3]<<",";
            r[0]=x,r[2]=y,r[1]=vx,r[3]=vy,r[4]=z,r[5]=vz;
            Kepler_3D(m34,r,kp);
            cout<<kp[0]<<","<<kp[1]<<","<<kp[2]<<","<<kp[3];
        }
     
    }
    outfile.close();
    return 0;
}