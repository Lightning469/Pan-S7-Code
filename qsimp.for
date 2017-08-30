      SUBROUTINE qsimp(func,a,b,s)
      INTEGER JMAX
      REAL*8 a,b,func,s,EPS
c      EXTERNAL func
      PARAMETER (EPS=1.d-6, JMAX=20)
CU    USES trapzd
      INTEGER j
      REAL*8 os,ost,st
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j)
        s=(4.*st-ost)/3.
        if (dabs(s-os).lt.EPS*dabs(os).or.
     *s.eq.0..and.os.eq.0.) return
        os=s
        ost=st
11    continue
      pause 'too many steps in qsimp'
      END
