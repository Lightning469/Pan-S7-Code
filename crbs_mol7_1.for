      parameter (n=7,n_omega=300,nm=14)	! 7 moments
      complex*16 a(n,n), b(n), x(n)
      complex*16 w0,w1,w2,w3,w4,w5,w6,unit,factor,zero,z,one,
	*i0000,i0100,i0001,i0010,i1000,i0011,i1100,i0101,i0110,i1001,
     *i1101,i0111,i1111,i0002,i1102,i0211,i0202,i0210,i1002,i1010,
	*i0201,i0200,i1110,i0102,i1011 
	real*8 k,kb,lambda,m_m,ng,x_l,x_r,integr1,integr2
	real*8 j020,j100,j030,j110,j011110,j011,area,
     *	signal(-n_omega:n_omega),omega(-n_omega:n_omega),
	*	am(nm,nm),bm(nm),xm(nm)
	integer i_case ! 0 - spontaneous; 1 - coherent
      EXTERNAL func,func1
	real*8 func,func1
	common sk,y7										    
	open(1,file='signal.dat')

	i_case=0	! 0 - SRBS ; any other, CRBS

		  unit=(0.d0,1.d0)						    
		  zero=(0.d0,0.d0)
		  one=(1.d0,0.d0)	

	omega(-n_omega)=-3.d10
	d_omega=dfloat(abs(omega(-n_omega))/n_omega)
	do i=-n_omega+1,n_omega
	omega(i)=omega(i-1)+d_omega
	enddo

	     pi=3.141592654
c	     dt=0.001

		x_l=sk-15.
		x_r=sk+15.


	tem=292.
	pressure_atm=1.5877		   

	pressure_pa=pressure_atm*1.01325e5
	pressure_torr=pressure_pa*0.00750061683

	ng=3.295e22*pressure_torr*293./tem  ! gas molecule density (per 1 m^3) 

c	Nitrogen
      m_m=1.66d-27*28. ! N2
		c_int=1.
		c_tr=3./2.
	viscosity=1.e-6*(0.20882+0.07255*tem-5.1641E-5*(tem**2)+
	*1.93243E-8*(tem**3))
	bulk_vis=(128./175.)*viscosity
	heat_cond=1.e-3*(0.30778+0.09855*tem-5.25911E-5*(tem**2)+
	*1.91616E-8*(tem**3))
c	Nitrogen

c	 CO2
c      m_m=1.66e-27*44. ! CO2
c		c_int=1.
c		c_tr=3./2.
c	viscosity=1.e-6*(-0.50373+0.05687*tem-1.72921E-5*(tem**2))
c	bulk_vis=0.25*viscosity
c	heat_cond=0.25*(1.4*9-5.)*viscosity*(c_int+c_tr)*188.93	! R=188.93 J/kg*K

c	heat_cond=1.e-3*(-1.02474+0.03545*tem+9.68252E-5*(tem**2)-
c	*6.67982E-8*(tem**3))  ! low frequency limit!****
c	 CO2


c     Argon
c      viscosity=1.e-6*(0.30075+0.0867*tem-4.02559E-5*(tem**2)) ! Ar, 200<T<500 K	 (Pa*sec)
c      m_m=1.66e-27*39.948 ! Ar
c		c_int=1.e-33
c		c_tr=3./2.
c      bulk_vis=1.e-33
c	heat_cond=0.0175

	 
	accel=1.
	
	kb=1.38d-23
	lambda=532.d-9  !	 1060.d-9
	k=(4*pi/lambda) !*sin(89*pi/180)
	
									v0=sqrt(2.*kb*tem/m_m)

	omega(-n_omega)=-k*v0*3. !3.d10
	d_omega=dfloat(abs(omega(-n_omega))/n_omega)
						do i=-n_omega+1,n_omega
						omega(i)=omega(i-1)+d_omega
						enddo



		y7=(3./2.)*ng*kb*tem/(k*v0*viscosity)
		y=2.*y7/3.

			if(i_case.eq.0) then
		write(*,*) 'SRBS: 7 moments  ','y=',y,'   p=',pressure_atm
		else
		write(*,*) 'CRBS: 7 moments  ','y=',y,'   p= ',pressure_atm
			endif

c	 J's coefficients 

       	gamma_int=c_int/(c_tr+c_int)
		j020=-ng*kb*tem/viscosity
		j100=-(2./3.)*ng*(kb*tem/bulk_vis)*((c_int/(c_tr+c_int))**2)
		j030=(3./2.)*j020  
c		j110=-ng*kb*tem*(2./(3.*viscosity)-
c     *5.*(gamma_int**2)/(9.*bulk_vis))
		j110=5.*j100/6.+2.*j020/3.
c	j011110=
c	*-ng*sqrt(5./(18.*c_int))*(gamma_int**2)*kb*tem/bulk_vis
		j011110=sqrt(5./(8.*c_int))*j100
	aj1=-(2./3.)*gamma_int*kb*tem/((c_tr+c_int)*bulk_vis)
	aj2=(2.*bulk_vis/(5.*viscosity))*((c_tr+c_int)**2)+
     *c_int*(1.+c_int/3.)+(gamma_int**2)*m_m*heat_cond/(6.*kb*bulk_vis)
	aj3=-1.+4.*m_m*heat_cond/(15.*kb*viscosity)+
	*2.*(gamma_int**2)*m_m*heat_cond/(9.*kb*bulk_vis)
		j011=ng*aj1*aj2/aj3

			area=0.
			do 100 i=-n_omega,n_omega     
			sk=omega(i)/(k*v0)

		z=dcmplx(sk,y7)

      call qsimp(func,x_l,x_r,integr1)
      call qsimp(func1,x_l,x_r,integr2)

      w0=dcmplx(integr1,-integr2)

c		w0=zero
c
c		do t=x_l,x_r,dt
c		w0=w0+exp(-t*t)*dt/(z-t)
c		enddo

			
			w1=-sqrt(pi)+z*w0
			w2=z*w1
			w3=-0.5*sqrt(pi)+z*w2
			w4=z*w3
			w5=-3.*sqrt(pi)/4.+z*w4
			w6=z*w5 

c	i - coefficients: I_b1_b2_u1_u2

		factor=unit/(k*v0)
      i0000=factor*w0/(sqrt(pi))
	i0100=factor*(z*w0-sqrt(pi))*sqrt(2./pi)
	i0001=i0100
	i0010=factor*(2.*w2-w0)/(sqrt(6.*pi))
	i1000=i0010
	i0011=factor*(2.*w3-3.*w1)/(sqrt(5.*pi))
	i1100=i0011
	i0101=factor*2.*w2/sqrt(pi)
	i0110=factor*(-w1+2*w3)/sqrt(3.*pi)
	i1001=i0110
	i0111=factor*(-3.*w2+2.*w4)*sqrt(2./(5.*pi))
	i1101=i0111
	i1111=factor*(13.*w2-12.*w4+4.*w6)/(5.*sqrt(pi))
      i0002=factor*(-w0+2.*w2)/sqrt(3.*pi)
	i0200=i0002
	i0211=factor*(-w1+8.*w3-4.*w5)/sqrt(15.*pi)
	i1102=i0211
	i0202=factor*2.*(w0-2.*w2+2.*w4)/(3.*sqrt(pi))
	i0210=factor*(w0+4.*w2-4.*w4)/(3.*sqrt(2.*pi))
	i1002=i0210

	i0102=factor*(-w1+2.*w3)*sqrt(2./(3.*pi))
	i0201=i0102

      i1010=factor*(5.*w0-4.*w2+4.*w4)/(6.*sqrt(pi))
	i1110=factor*(7.*w1-8.*w3+4.*w5)/sqrt(30.*pi)
	i1011=i1110

c	matrix elements
c	matrix A
	a(1,1)=-j030*i0000-one
	a(2,1)=-j030*i0001
	a(3,1)=-j030*i0011
	a(4,1)=-j030*i0002
	a(5,1)=-j030*i0010
	a(6,1)=zero
	a(7,1)=zero

	a(1,2)=-j030*i0100
	a(2,2)=-j030*i0101-one
	a(3,2)=-j030*i0111
	a(4,2)=-j030*i0102
	a(5,2)=-j030*i0110
	a(6,2)=zero
	a(7,2)=zero

	a(1,3)=(j030-j110)*i1100
	a(2,3)=(j030-j110)*i1101
	a(3,3)=(j030-j110)*i1111+one
	a(4,3)=(j030-j110)*i1102
	a(5,3)=(j030-j110)*i1110
	a(6,3)=-j011110*i0100
	a(7,3)=-j011110*i0101

	a(1,4)=(j020-j030)*i0200
	a(2,4)=(j020-j030)*i0201
	a(3,4)=(j020-j030)*i0211
	a(4,4)=(j020-j030)*i0202-(3./2.)*one
	a(5,4)=(j020-j030)*i0210
	a(6,4)=zero
	a(7,4)=zero

	a(1,5)=(j030-j100)*i1000
	a(2,5)=(j030-j100)*i1001
	a(3,5)=(j030-j100)*i1011
	a(4,5)=(j030-j100)*i1002
	a(5,5)=(j030-j100)*i1010+one
	a(6,5)=-j100*i0000*sqrt(c_tr/c_int)
	a(7,5)=-j100*i0001*sqrt(c_tr/c_int)

	a(1,6)=j100*i1000
	a(2,6)=j100*i1001
	a(3,6)=j100*i1011
	a(4,6)=j100*i1002
	a(5,6)=j100*i1010
	a(6,6)=((j100*c_tr/c_int-j030)*i0000-one)*sqrt(c_int/c_tr)
	a(7,6)=(j100*c_tr/c_int-j030)*i0001*sqrt(c_int/c_tr)

	a(1,7)=j011110*i1100
	a(2,7)=j011110*i1101
	a(3,7)=j011110*i1111
	a(4,7)=j011110*i1102
	a(5,7)=j011110*i1110
	a(6,7)=(j011-j030)*i0100
	a(7,7)=(j011-j030)*i0101-one

	b_factor=m_m*accel/sqrt(kb*tem)
	if(i_case.eq.0) goto 500
c	CRBS
		b(1)=-b_factor*i0100
		b(2)=-b_factor*i0101
		b(3)=-b_factor*i0111
		b(4)=-b_factor*i0102
		b(5)=-b_factor*i0110
		b(6)=zero
		b(7)=zero
		goto 510

c	Tenti's case (spontaneous)
500      b(1)=-b_factor*i0000
	b(2)=-b_factor*i0001
	b(3)=-b_factor*i0011
	b(4)=-b_factor*i0002
	b(5)=-b_factor*i0010
	b(6)=zero
	b(7)=zero

 510  continue
	do ii=1,nm
	do j=1,nm
	if(ii.le.n.and.j.le.n)	am(ii,j)=dreal(a(ii,j))
	if(ii.gt.n.and.j.le.n)	am(ii,j)=dimag(a(ii-n,j))
	if(ii.le.n.and.j.gt.n)	am(ii,j)=-dimag(a(ii,j-n))
	if(ii.gt.n.and.j.gt.n)	am(ii,j)=dreal(a(ii-n,j-n))
		enddo
		enddo
	     do j=1,nm
		 if(j.le.n) then
		 bm(j)=dreal(b(j))
		 else
		 bm(j)=dimag(b(j-n))
		endif
		enddo


c	ludcmp replaces a given matrix a() by the LU decomposition 
c	of a rowwise permutation of itself
      call ludcmp(am,nm,indx)
c	lubksb	solves the set of linear equations A*X=B
	call lubksb(am,nm,indx,bm,xm)
	x(1)=dcmplx(xm(1),xm(nm-n+1))
		if(i_case.eq.0)	then
	 	signal(i)=real(x(1)) ! SRBS
      	else
		signal(i)=x(1)*conjg(x(1))    ! CRBS
		endif     
		area=area+signal(i)*d_omega/(k*v0)   ! area under the signal plot
c	signal_t=x(5)*conjg(x(5))/(b_factor**2)
100	continue
		do i=-n_omega,n_omega
		phase_velocity=omega(i)/k
	gauss=exp(-0.5*m_m*(phase_velocity**2)/kb/tem)
	write(1,10) omega(i)/(2*pi)/1.e9,signal(i)/area,
     *	phase_velocity,gauss !    signal(0), !
omega(i)/(k*v0), signal(i)/area,  !    signal(0), !

		end do

   10 format(4(2x,e12.5))
  	stop
	end

	real*8 function func(t)
	real*8 t 
	common x,y
		func=(x-t)*dexp(-t**2)/((x-t)**2+y**2)
	return
	end
		real*8 function func1(t)
		real*8 t 
	common x,y
		func1=y*dexp(-t**2)/((x-t)**2+y**2)
	return
	end


