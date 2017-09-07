%This code was transferred from Mikhail Shneider's S7 code for RBS. The
%full code can be found at the end of this .m file. That Fortran code was
%transferred to Matlab syntax by Emanuel Stockman. Any time this code (or
%transcript of it) is used in work related to a publication, it MUST be
%mentioned that the S7 code was developed at Princeton and should reference
%the original paper:
%1. X.P. Pan, M.N. Shneider, R.B. Miles, Coherent Rayleigh-Brillouin
% Scattering in Molecular Gases, Phys. Rev. A, Vol 69, 033814, 1-16 (2004).

%All of the control parameters are at the beginning of the program
%INPUTS:
%   Pinputatm = gas pressure in atm
%   Tinput = gas temperature in Kelvin
%   gastypeinput = 1-N2, 2-CO2, 3-Ar, 4-C3H8 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAYLEIGH PROFILES FOR CO2, Ar, and C3H8 ARE NOT EXPERIMENTALLY VERIFIED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT:
%   nusHz = frequency axis in Hz
%   signalNorm = RBS signal, normalized to 1 with respect to frequency
%   y = y parameter
%SAMPLE CALL:
%   [nusHz,signalNorm,y]=Mikhail_S7(1.0,273,1)
%   This statement will provide the Rayleigh profile for Nitrogen at 1 atm
%   and 273 K.
%   Typically, this statement is executed in a loop to provide a series of
%   profiles at different temperatures and pressures and then stored in a
%   look up table database. It must be clear that modeling mixtures of
%   gases is not as straight forward as multiplying the relative profile by
%   the molar concentration weighted Rayleigh cross section. 

%%
function [nusHz,signalNorm,y]=Mikhail_S7(PinputAtm,Tinput,gastypeinput)

WN_0 = 1 / 532.26E-7; % 18788.4; % Laser wave number (set your laser wavelength here), i.e., 1/lambda [1/cm]
d_nu=5e9;          %half-range of frequency axis [Hz]
n_omegas=1000;      %number of points in the grid

i_case=0; %0-SpontaneousRBS and 1-CoherentRBS (note: original setting was 0)
% Spontaneous, e.g., gas in test cell
% Coherent , e.g., moving grating formed by crossing laser beams

n=7;                    %number of momemtns
nm=14;                  % number of moments *2???

%Gas properties
gastype=gastypeinput;      %1-N2, 2-CO2, 3-Ar, 4-C3H8
T_curr=Tinput;
P_currAtm =PinputAtm;  %pressure in atm
P_currPa = P_currAtm * 1.01325e5;
P_currTorr = P_currPa * .00750061683;
ng=3.295e22 * P_currTorr * 293 / T_curr;  % gas molecule density (per 1 m^3)

%Constants
c=2.9979e8; % Speed of light [m/s]
kb=1.38d-23; % Boltzmann constant [J/K]
R = 8.14472;    % Ideal gas constant J/mol/K
Navog = 6.0221353e23 ;   %Avogadro's const [1/mol]

%Laser properties
%nu_0=WN_0*c/2/pi;  %Laser frequency [?]
lam_0=1/WN_0/100;  % Laser wavelength [m]
theta = 90*(pi/180);
% I think theta is the crossing angle of the pump beams that generate the optical grating
k=(4*pi/lam_0)*sin(theta/2); % Magnitude of scattering wavevector (why is theta = pi?) [1/m]

range_omegas=2*pi*d_nu;%k*v0*3;
omegas=linspace(-range_omegas,range_omegas,2*n_omegas);
d_omega=mean(diff(omegas));

%%

if gastype==1
    %%%%%Air
    m_m=(1.66053886e-27)*28.9645; % N2
    c_int=1;
    c_tr=3/2;
    viscosity=1.e-6*(0.20882 + 0.07255*T_curr - 5.1641e-5*(T_curr^2) + 1.93243E-8*(T_curr^3));
    heat_cond=1.e-3*(0.30778 + 0.09855*T_curr - 5.25911E-5*(T_curr^2)+ 1.91616E-8*(T_curr^3));
    if  T_curr>1000
        viscosity=40.6e-6*sqrt(T_curr/1000);
        heat_cond=65.3e-3*sqrt(T_curr/1000);
    end
    bulk_vis=(128/175)*viscosity;

elseif gastype==2
    %%%	 CO2
    m_m=1.66e-27*44.; % CO2
    c_int=1.;
    c_tr=3./2.;
    viscosity=1.e-6*(-0.50373+0.05687*T_curr-1.72921E-5*(T_curr^2));
    bulk_vis=0.25*viscosity;
    heat_cond=0.25*(1.4*9-5.)*viscosity*(c_int+c_tr)*188.93;  % R=188.93 J/kg*K
  %      	heat_cond=1.e-3*(-1.02474+0.03545*T_curr+9.68252E-5*(T_curr^2)-*6.67982E-8*(T_curr^3))  ! low frequency limit!^^
  elseif gastype==3
      %% Argon
      viscosity=1.e-6*(0.30075+0.0867*T_curr-4.02559E-5*(T_curr^2)); % Ar, 200<T<500 K	 (Pa*sec)
      m_m=1.66e-27*39.948; % Ar
      c_int=1.e-33;
      c_tr=3./2.;
      bulk_vis=1.e-33;
      heat_cond=0.0175;
elseif gastype==4
    %%Propane information taken from NIST Paper:
    %
    %"A Reference Multiparameter Viscosity Equation for Propane with an
    %Optimized Functional Form" by Scalabrin, Marchi, and Span
    %Conductivity%
    %An improved empirical correlation for the thermal conductivity of
    %propane." By Ramires, de Castro, Perkins in Intl. Journ. of
    %Thermophysics
%%
    m_m = 1.66e-27 * 44.098;
    MW = 0.044098;  %kg/mol
    c_int = 1;
    c_tr = 3/2;
    Tc = 369.825; Pc=4.248*1e6; %critical values
    Tr = T_curr/Tc;
    Hc = MW^.5 * Pc^(2/3) / R^(1/6) / Navog^(1/3) / Tc^(1/6);
    viscosity= Hc * (0.5947468 * Tr + -0.01062820*Tr^4.5 + 0.002142285*Tr^6);
    heat_cond= -1.25e-2 + 8.42e-5*T_curr + 6.02e-8*T_curr^2;      %W/m/K   Dilute gas limit
    bulk_vis=0.25*viscosity;
%%
end
%Match Mikhail's variable nomenclature
tem=T_curr;
pressure_atm=P_currAtm;
pressure_pa=pressure_atm*1.01325e5;
pressure_torr=pressure_pa*0.00750061683;
unit=i;   %(0.d0,1.d0)
zero=0;   %(0.d0,0.d0)
one=1;    %(1.d0,0.d0)

accel=1;

v0=sqrt(2*kb*T_curr/m_m); % Most probably speed in 3D

y7=(3./2.)*ng*kb*tem/(k*v0*viscosity);
y=2.*y7/3;

if i_case==0
    display(['SRBS: 7 moments  y = ' num2str(y) ', p = ' num2str(P_currAtm) ' atm, T = ' num2str(T_curr), ' K, lambda_0 = ' num2str(lam_0), ' m']);
else
    display(['CRBS: 7 moments  y = ' num2str(y) ', p = ' num2str(P_currAtm) ' atm, T = ' num2str(T_curr), ' K, lambda_0 = ' num2str(lam_0), ' m']);
end


%%
%J's coefficients
gamma_int=c_int/(c_tr+c_int);
j020=-ng*kb*tem/viscosity;
j100=-(2/3)*ng*(kb*tem/bulk_vis)*((c_int/(c_tr+c_int))^2);
j030=(3./2.)*j020;
%	j110=-ng*kb*tem*(2./(3.*viscosity)-5.*(gamma_int^2)/(9.*bulk_vis))
j110=5.*j100/6.+2.*j020/3;
%	j011110= -ng*sqrt(5./(18.*c_int))*(gamma_int^2)*kb*tem/bulk_vis
j011110=sqrt(5./(8.*c_int))*j100;

aj1=-(2./3.)*gamma_int*kb*tem/((c_tr+c_int)*bulk_vis);
aj2=(2.*bulk_vis/(5.*viscosity))*((c_tr+c_int)^2)+c_int*(1.+c_int/3.)+(gamma_int^2)*m_m*heat_cond/(6.*kb*bulk_vis);
aj3=-1.+4.*m_m*heat_cond/(15.*kb*viscosity)+2.*(gamma_int^2)*m_m*heat_cond/(9.*kb*bulk_vis);

j011=ng*aj1*aj2/aj3;

area=0;
for hh=1:(2*n_omegas)
    sk=omegas(hh)/(k*v0);
    z= sk + i*y7;

    x_l=sk-15;
	x_r=sk+15;
    
%   dt=0.001;
%     w0=zero;
%     for t=x_l:dt:x_r
% 		w0=w0+exp(-t*t)*dt/(z-t);
%     end

     tt=linspace(x_l,x_r,10000);
     w0=trapz(tt,exp(-tt.^2)./(z-tt));
     clear t;
    
    w1=-sqrt(pi)+z*w0;
    w2=z*w1;
    w3=-0.5*sqrt(pi)+z*w2;
    w4=z*w3;
    w5=-3.*sqrt(pi)/4.+z*w4;
    w6=z*w5;


    %	i - coefficients: I_b1_b2_u1_u2

    factor=unit/(k*v0);
    i0000=factor*w0/(sqrt(pi));
    i0100=factor*(z*w0-sqrt(pi))*sqrt(2./pi);
    i0001=i0100;
    i0010=factor*(2.*w2-w0)/(sqrt(6.*pi));
    i1000=i0010;
    i0011=factor*(2.*w3-3.*w1)/(sqrt(5.*pi));
    i1100=i0011;
    i0101=factor*2.*w2/sqrt(pi);
    i0110=factor*(-w1+2*w3)/sqrt(3.*pi);
    i1001=i0110;
    i0111=factor*(-3.*w2+2.*w4)*sqrt(2./(5.*pi));
    i1101=i0111;
    i1111=factor*(13.*w2-12.*w4+4.*w6)/(5.*sqrt(pi));
    i0002=factor*(-w0+2.*w2)/sqrt(3.*pi);
    i0200=i0002;
    i0211=factor*(-w1+8.*w3-4.*w5)/sqrt(15.*pi);
    i1102=i0211;
    i0202=factor*2.*(w0-2.*w2+2.*w4)/(3.*sqrt(pi));
    i0210=factor*(w0+4.*w2-4.*w4)/(3.*sqrt(2.*pi));
    i1002=i0210;

    i0102=factor*(-w1+2.*w3)*sqrt(2./(3.*pi));
    i0201=i0102;

    i1010=factor*(5.*w0-4.*w2+4.*w4)/(6.*sqrt(pi));
    i1110=factor*(7.*w1-8.*w3+4.*w5)/sqrt(30.*pi);
    i1011=i1110;

    %	matrix elements
    %	matrix A
    a(1,1)=-j030*i0000-one;
    a(2,1)=-j030*i0001;
    a(3,1)=-j030*i0011;
    a(4,1)=-j030*i0002;
    a(5,1)=-j030*i0010;
    a(6,1)=zero;
    a(7,1)=zero;

    a(1,2)=-j030*i0100;
    a(2,2)=-j030*i0101-one;
    a(3,2)=-j030*i0111;
    a(4,2)=-j030*i0102;
    a(5,2)=-j030*i0110;
    a(6,2)=zero;
    a(7,2)=zero;

    a(1,3)=(j030-j110)*i1100;
    a(2,3)=(j030-j110)*i1101;
    a(3,3)=(j030-j110)*i1111+one;
    a(4,3)=(j030-j110)*i1102;
    a(5,3)=(j030-j110)*i1110;
    a(6,3)=-j011110*i0100;
    a(7,3)=-j011110*i0101;

    a(1,4)=(j020-j030)*i0200;
    a(2,4)=(j020-j030)*i0201;
    a(3,4)=(j020-j030)*i0211;
    a(4,4)=(j020-j030)*i0202-(3./2.)*one;
    a(5,4)=(j020-j030)*i0210;
    a(6,4)=zero;
    a(7,4)=zero;

    a(1,5)=(j030-j100)*i1000;
    a(2,5)=(j030-j100)*i1001;
    a(3,5)=(j030-j100)*i1011;
    a(4,5)=(j030-j100)*i1002;
    a(5,5)=(j030-j100)*i1010+one;
    a(6,5)=-j100*i0000*sqrt(c_tr/c_int);
    a(7,5)=-j100*i0001*sqrt(c_tr/c_int);

    a(1,6)=j100*i1000;
    a(2,6)=j100*i1001;
    a(3,6)=j100*i1011;
    a(4,6)=j100*i1002;
    a(5,6)=j100*i1010;
    a(6,6)=((j100*c_tr/c_int-j030)*i0000-one)*sqrt(c_int/c_tr);
    a(7,6)=(j100*c_tr/c_int-j030)*i0001*sqrt(c_int/c_tr);

    a(1,7)=j011110*i1100;
    a(2,7)=j011110*i1101;
    a(3,7)=j011110*i1111;
    a(4,7)=j011110*i1102;
    a(5,7)=j011110*i1110;
    a(6,7)=(j011-j030)*i0100;
    a(7,7)=(j011-j030)*i0101-one;

    b_factor=m_m*accel/sqrt(kb*tem);
    if i_case==1
        %	CRBS
        b(1)=-b_factor*i0100;
        b(2)=-b_factor*i0101;
        b(3)=-b_factor*i0111;
        b(4)=-b_factor*i0102;
        b(5)=-b_factor*i0110;
        b(6)=zero;
        b(7)=zero;
    else
        %	Tenti's case (spontaneous)
        b(1)=-b_factor*i0000;
        b(2)=-b_factor*i0001;
        b(3)=-b_factor*i0011;
        b(4)=-b_factor*i0002;
        b(5)=-b_factor*i0010;
        b(6)=zero;
        b(7)=zero;
    end
    for ii=1:nm
        for jj=1:nm
            if (ii<=n) && (jj<=n)
                am(ii,jj)=real(a(ii,jj));
            end
            if (ii>n) && (jj<=n)
                am(ii,jj)=imag(a(ii-n,jj));
            end
            if (ii<=n) && (jj>n)
                am(ii,jj)=-imag(a(ii,jj-n));
            end
            if (ii>n) && (jj>n)
                am(ii,jj)=real(a(ii-n,jj-n));
            end
        end
    end
    bm=zeros(nm,1);             %bm should be column matrix
    for kk=1:nm
        if(kk<=n)
            bm(kk)=real(b(kk));
        else
            bm(kk)=imag(b(kk-n));
        end
    end

    %	ludcmp replaces a given matrix a() by the LU decomposition
    %	of a rowwise permutation of itself
    [LfromA,UfromA]=lu(am);    %call ludcmp(am,nm,indx)
    %	lubksb	solves the set of linear equations A*X=B
    xm=(am^-1)*bm;                         %call lubksb(am,nm,indx,bm,xm)

    x(1)=xm(1) + i*xm(nm-n+1);
    
    if (i_case==0)
        signal(hh)=real(x(1)); % SRBS
    else
        signal(hh)=x(1)*conj(x(1));    % CRBS
    end
    area = area + signal(hh)*d_omega/(k*v0);   % area under the signal plot
end

nusHz=omegas/2/pi;
normalization=trapz(nusHz,signal);
signalNorm=signal/normalization;
% signalNorm=signal;

% nusHz=linspace(-10e9,0,2000);
% signalNorm=interp1(signal(nus<0),nus(nus<0),nusHz);
% signalNorm(isnan(signalNorm))=0;


end












































%%ORIGINAL FORTRAN CODE

%       parameter (n=7,n_omega=300,nm=14)	! 7 moments
%       complex*16 a(n,n), b(n), x(n)
%       complex*16 w0,w1,w2,w3,w4,w5,w6,unit,factor,zero,z,one,
% 	*i0000,i0100,i0001,i0010,i1000,i0011,i1100,i0101,i0110,i1001,
%      *i1101,i0111,i1111,i0002,i1102,i0211,i0202,i0210,i1002,i1010,
% 	*i0201,i0200,i1110,i0102,i1011 
% 	real*8 k,kb,lambda,m_m,ng,x_l,x_r,integr1,integr2
% 	real*8 j020,j100,j030,j110,j011110,j011,area,
%      *	signal(-n_omega:n_omega),omega(-n_omega:n_omega),
% 	*	am(nm,nm),bm(nm),xm(nm)
% 	integer i_case ! 0 - spontaneous; 1 - coherent
%       EXTERNAL func,func1
% 	real*8 func,func1
% 	common sk,y7										    
% 	open(1,file='292K.dat')
% 
% 	i_case=1	! 0 - SRBS ; any other, CRBS
% 
% 		  unit=(0.d0,1.d0)						    
% 		  zero=(0.d0,0.d0)
% 		  one=(1.d0,0.d0)	
% 
% 	omega(-n_omega)=-3.d10
% 	d_omega=dfloat(abs(omega(-n_omega))/n_omega)
% 	do i=-n_omega+1,n_omega
% 	omega(i)=omega(i-1)+d_omega
% 	enddo
% 
% 	     pi=3.141592654
% 	     dt=0.001
% 
% 		x_l=sk-15.
% 		x_r=sk+15.
% 
% 
% 	tem=292.
% 	pressure_atm=5.	   
% 
% 	pressure_pa=pressure_atm*1.01325e5
% 	pressure_torr=pressure_pa*0.00750061683
% 
% 	ng=3.295e22*pressure_torr*293./tem  ! gas molecule density (per 1 m^3) 
% 
% c	Nitrogen
%       m_m=1.66d-27*28. ! N2
% 		c_int=1.
% 		c_tr=3./2.
% 	viscosity=1.e-6*(0.20882+0.07255*tem-5.1641E-5*(tem**2)+
% 	*1.93243E-8*(tem**3))  
%      	heat_cond=1.e-3*(0.30778+0.09855*tem-5.25911E-5*(tem**2)+
%      *1.91616E-8*(tem**3))
%     	if(tem.gt.1000.) then
% 	viscosity=40.6e-6*sqrt(tem/1000.)
% 	heat_cond=65.3e-3*sqrt(tem/1000.)
% 	endif
% 	bulk_vis=(128./175.)*viscosity
% c	Nitrogen
% 
% c	 CO2
% c      m_m=1.66e-27*44. ! CO2
% c		c_int=1.
% c		c_tr=3./2.
% c	viscosity=1.e-6*(-0.50373+0.05687*tem-1.72921E-5*(tem**2))
% c	bulk_vis=0.25*viscosity
% c	heat_cond=0.25*(1.4*9-5.)*viscosity*(c_int+c_tr)*188.93	! R=188.93 J/kg*K
% 
% c	heat_cond=1.e-3*(-1.02474+0.03545*tem+9.68252E-5*(tem**2)-
% c	*6.67982E-8*(tem**3))  ! low frequency limit!****
% c	 CO2
% 
% 
% c     Argon
% c      viscosity=1.e-6*(0.30075+0.0867*tem-4.02559E-5*(tem**2)) ! Ar, 200<T<500 K	 (Pa*sec)
% c      m_m=1.66e-27*39.948 ! Ar
% c		c_int=1.e-33
% c		c_tr=3./2.
% c      bulk_vis=1.e-33
% c	heat_cond=0.0175
% 
% 	 
% 	accel=1.
% 	
% 	kb=1.38d-23
% 	lambda=532.26d-9  !	 1060.d-9
% 	k=(4*pi/lambda)*sin(pi/4.) 
% 
% 	
% 									v0=sqrt(2.*kb*tem/m_m)
% 
% 	omega(-n_omega)=-k*v0*3. !3.d10
% 	d_omega=dfloat(abs(omega(-n_omega))/n_omega)
% 						do i=-n_omega+1,n_omega
% 						omega(i)=omega(i-1)+d_omega
% 						enddo
% 
% 
% 
% 		y7=(3./2.)*ng*kb*tem/(k*v0*viscosity)
% 		y=2.*y7/3.
% 
% 			if(i_case.eq.0) then
% 		write(*,*) 'SRBS: 7 moments  ','y=',y,
%      *'   p=',pressure_atm,'   T=',tem  
% 		else
% 		write(*,*) 'CRBS: 7 moments  ','y=',y,'   p= ',pressure_atm
% 			endif
% c	 J's coefficients 
% 
%        	gamma_int=c_int/(c_tr+c_int)
% 		j020=-ng*kb*tem/viscosity
% 		j100=-(2./3.)*ng*(kb*tem/bulk_vis)*((c_int/(c_tr+c_int))**2)
% 		j030=(3./2.)*j020  
% c		j110=-ng*kb*tem*(2./(3.*viscosity)-
% c     *5.*(gamma_int**2)/(9.*bulk_vis))
% 		j110=5.*j100/6.+2.*j020/3.
% c	j011110=
% c	*-ng*sqrt(5./(18.*c_int))*(gamma_int**2)*kb*tem/bulk_vis
% 		j011110=sqrt(5./(8.*c_int))*j100
% 	aj1=-(2./3.)*gamma_int*kb*tem/((c_tr+c_int)*bulk_vis)
% 	aj2=(2.*bulk_vis/(5.*viscosity))*((c_tr+c_int)**2)+
%      *c_int*(1.+c_int/3.)+(gamma_int**2)*m_m*heat_cond/(6.*kb*bulk_vis)
% 	aj3=-1.+4.*m_m*heat_cond/(15.*kb*viscosity)+
% 	*2.*(gamma_int**2)*m_m*heat_cond/(9.*kb*bulk_vis)
% 		j011=ng*aj1*aj2/aj3
% 
% 			area=0.
% 			do 100 i=-n_omega,n_omega     
% 			sk=omega(i)/(k*v0)
% 
% 		z=dcmplx(sk,y7)
% 
%       call qsimp(func,x_l,x_r,integr1)
%       call qsimp(func1,x_l,x_r,integr2)
% 
%       w0=dcmplx(integr1,-integr2)
% 
% c		w0=zero
% c
% c		do t=x_l,x_r,dt
% c		w0=w0+exp(-t*t)*dt/(z-t)
% c		enddo
% 
% 			
% 			w1=-sqrt(pi)+z*w0
% 			w2=z*w1
% 			w3=-0.5*sqrt(pi)+z*w2
% 			w4=z*w3
% 			w5=-3.*sqrt(pi)/4.+z*w4
% 			w6=z*w5 
% 
% c	i - coefficients: I_b1_b2_u1_u2
% 
% 		factor=unit/(k*v0)
%       i0000=factor*w0/(sqrt(pi))
% 	i0100=factor*(z*w0-sqrt(pi))*sqrt(2./pi)
% 	i0001=i0100
% 	i0010=factor*(2.*w2-w0)/(sqrt(6.*pi))
% 	i1000=i0010
% 	i0011=factor*(2.*w3-3.*w1)/(sqrt(5.*pi))
% 	i1100=i0011
% 	i0101=factor*2.*w2/sqrt(pi)
% 	i0110=factor*(-w1+2*w3)/sqrt(3.*pi)
% 	i1001=i0110
% 	i0111=factor*(-3.*w2+2.*w4)*sqrt(2./(5.*pi))
% 	i1101=i0111
% 	i1111=factor*(13.*w2-12.*w4+4.*w6)/(5.*sqrt(pi))
%       i0002=factor*(-w0+2.*w2)/sqrt(3.*pi)
% 	i0200=i0002
% 	i0211=factor*(-w1+8.*w3-4.*w5)/sqrt(15.*pi)
% 	i1102=i0211
% 	i0202=factor*2.*(w0-2.*w2+2.*w4)/(3.*sqrt(pi))
% 	i0210=factor*(w0+4.*w2-4.*w4)/(3.*sqrt(2.*pi))
% 	i1002=i0210
% 
% 	i0102=factor*(-w1+2.*w3)*sqrt(2./(3.*pi))
% 	i0201=i0102
% 
%       i1010=factor*(5.*w0-4.*w2+4.*w4)/(6.*sqrt(pi))
% 	i1110=factor*(7.*w1-8.*w3+4.*w5)/sqrt(30.*pi)
% 	i1011=i1110
% 
% c	matrix elements
% c	matrix A
% 	a(1,1)=-j030*i0000-one
% 	a(2,1)=-j030*i0001
% 	a(3,1)=-j030*i0011
% 	a(4,1)=-j030*i0002
% 	a(5,1)=-j030*i0010
% 	a(6,1)=zero
% 	a(7,1)=zero
% 
% 	a(1,2)=-j030*i0100
% 	a(2,2)=-j030*i0101-one
% 	a(3,2)=-j030*i0111
% 	a(4,2)=-j030*i0102
% 	a(5,2)=-j030*i0110
% 	a(6,2)=zero
% 	a(7,2)=zero
% 
% 	a(1,3)=(j030-j110)*i1100
% 	a(2,3)=(j030-j110)*i1101
% 	a(3,3)=(j030-j110)*i1111+one
% 	a(4,3)=(j030-j110)*i1102
% 	a(5,3)=(j030-j110)*i1110
% 	a(6,3)=-j011110*i0100
% 	a(7,3)=-j011110*i0101
% 
% 	a(1,4)=(j020-j030)*i0200
% 	a(2,4)=(j020-j030)*i0201
% 	a(3,4)=(j020-j030)*i0211
% 	a(4,4)=(j020-j030)*i0202-(3./2.)*one
% 	a(5,4)=(j020-j030)*i0210
% 	a(6,4)=zero
% 	a(7,4)=zero
% 
% 	a(1,5)=(j030-j100)*i1000
% 	a(2,5)=(j030-j100)*i1001
% 	a(3,5)=(j030-j100)*i1011
% 	a(4,5)=(j030-j100)*i1002
% 	a(5,5)=(j030-j100)*i1010+one
% 	a(6,5)=-j100*i0000*sqrt(c_tr/c_int)
% 	a(7,5)=-j100*i0001*sqrt(c_tr/c_int)
% 
% 	a(1,6)=j100*i1000
% 	a(2,6)=j100*i1001
% 	a(3,6)=j100*i1011
% 	a(4,6)=j100*i1002
% 	a(5,6)=j100*i1010
% 	a(6,6)=((j100*c_tr/c_int-j030)*i0000-one)*sqrt(c_int/c_tr)
% 	a(7,6)=(j100*c_tr/c_int-j030)*i0001*sqrt(c_int/c_tr)
% 
% 	a(1,7)=j011110*i1100
% 	a(2,7)=j011110*i1101
% 	a(3,7)=j011110*i1111
% 	a(4,7)=j011110*i1102
% 	a(5,7)=j011110*i1110
% 	a(6,7)=(j011-j030)*i0100
% 	a(7,7)=(j011-j030)*i0101-one
% 
% 	b_factor=m_m*accel/sqrt(kb*tem)
% 	if(i_case.eq.0) goto 500
% c	CRBS
% 		b(1)=-b_factor*i0100
% 		b(2)=-b_factor*i0101
% 		b(3)=-b_factor*i0111
% 		b(4)=-b_factor*i0102
% 		b(5)=-b_factor*i0110
% 		b(6)=zero
% 		b(7)=zero
% 		goto 510
% 
% c	Tenti's case (spontaneous)
% 500      b(1)=-b_factor*i0000
% 	b(2)=-b_factor*i0001
% 	b(3)=-b_factor*i0011
% 	b(4)=-b_factor*i0002
% 	b(5)=-b_factor*i0010
% 	b(6)=zero
% 	b(7)=zero
% 
%  510  continue
% 	do ii=1,nm
% 	do j=1,nm
% 	if(ii.le.n.and.j.le.n)	am(ii,j)=dreal(a(ii,j))
% 	if(ii.gt.n.and.j.le.n)	am(ii,j)=dimag(a(ii-n,j))
% 	if(ii.le.n.and.j.gt.n)	am(ii,j)=-dimag(a(ii,j-n))
% 	if(ii.gt.n.and.j.gt.n)	am(ii,j)=dreal(a(ii-n,j-n))
% 		enddo
% 		enddo
% 	     do j=1,nm
% 		 if(j.le.n) then
% 		 bm(j)=dreal(b(j))
% 		 else
% 		 bm(j)=dimag(b(j-n))
% 		endif
% 		enddo
% 
% 
% c	ludcmp replaces a given matrix a() by the LU decomposition 
% c	of a rowwise permutation of itself
%       call ludcmp(am,nm,indx)
% c	lubksb	solves the set of linear equations A*X=B
% 	call lubksb(am,nm,indx,bm,xm)
% 	x(1)=dcmplx(xm(1),xm(nm-n+1))
% 		if(i_case.eq.0)	then
% 	 	signal(i)=real(x(1)) ! SRBS
%       	else
% 		signal(i)=x(1)*conjg(x(1))    ! CRBS
% 		endif     
% 		area=area+signal(i)*d_omega/(k*v0)   ! area under the signal plot
% c	signal_t=x(5)*conjg(x(5))/(b_factor**2)
% 100	continue
% 		do i=-n_omega,n_omega
% 		phase_velocity=omega(i)/k
% 	gauss=exp(-0.5*m_m*(phase_velocity**2)/kb/tem)
% 	write(1,10) omega(i)/(2*pi)/1.e9,signal(i)/area,
%      *	phase_velocity,gauss !    signal(0), !
% omega(i)/(k*v0), signal(i)/area,  !    signal(0), !
% 
% 		end do
% 
%    10 format(4(2x,e12.5))
%   	stop
% 	end
% 
% 	real*8 function func(t)
% 	real*8 t 
% 	common x,y
% 		func=(x-t)*dexp(-t**2)/((x-t)**2+y**2)
% 	return
% 	end
% 		real*8 function func1(t)
% 		real*8 t 
% 	common x,y
% 		func1=y*dexp(-t**2)/((x-t)**2+y**2)
% 	return
% 	end
% 
% 
