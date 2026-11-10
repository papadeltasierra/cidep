      program mode1
c
c.....Purpose : to simulate a mode1 esr signal + a misti mode1 esr
c.....signal c
c
c.....The equation to be solved is
c
c             .
c             M(t) = LM(t) + Peq.F/T1
c
c             L = | -1/T2    del        0      |
c                 | -del    -1/T2       omega1 |
c                 |  0      -omega1    -1/T1   |
c
c             F = | 0 |
c                 | 0 |
c                 | 1 |
c
c.....Initially, M(0) = Pi.F
c
c.....r(1) = T1
c.....r(2) = T2
c.....r(3) = Pi/Peq
c.....r(4) = omega1
c.....r(5) = offset from resonance
c.....r(6) = timespan
c.....r(7) = delay after flash
c
c.....nzero = npt/10
c.....  t = 0.0   corresponds to data(nzero-1)
c.....  t = r(6)  corresponds to data(npt)
c
c.....If you wish to see just the contribution of the initial term, input
c.....Pi/Peq = 1.0 and ask to exclude Peq.
c
c.....If you wish to see just the contribution of the eqm term, input
c.....Pi/Peq = 1.0 and ask to exclude Pi.
c
c     ------------------------------------------------------------------
c
      common/shape/data
c
      dimension data(2048,2)
      dimension r(7)
      character ans*1
c
      data r,npt/6.7,0.66,-50.0,1.068,0.0,20.48,0.0,1024/
c
10    continue
      peq=1.0
      if( r(3) .eq. 0.0 ) r(3)=1.0
      call input(r,npt,peq)
      call calc(r,npt,peq)
      nspec = 1
      if(r(7).ne.0.0) then
        nspec = 2
        call mcalc(r,npt,peq)
      endif
      call splot( data, npt, nspec )
      call clear
      print 27
      print 28
27    format(//)
28    format(8X,'Another go?(Y) ',\)
26    format(I2)
      read 25,ans
25    format(A1)
      if(ans.ne.'N'.and.ans.ne.'n') goto 10
      call clear
      stop
      end
c     ------------------------------------------------------------------
      subroutine calc(r,npt,peq)
c
      common/shape/data
c
      dimension data(2048,2)
      dimension r(7),ar(7),rs(6),rss(6)
      complex a(3),as(3),s(3),ss(3),c0,c1,c2,ct2,cw1,ctime,tinc
c
      equivalence ( rs, s ), ( rss, ss )
c
c.....cw1   = cmplx(r(3),0.0)
c.....ctime = cmplx(r(6),0.0)
c.....tinc  = cmplx(r(7),0.0)        time between samples
c.....t1    = 1.0/t1                 inverse of T1
c.....ct2   = cmplx(1.0/t2,0.0)      complex form of T2
c
      nzero=npt/10
      ncalc=npt-nzero+1
      do 10 i=1,nzero
10      data(i,1)=0.0
c
      t2=1.0/r(2)
      w1=r(4)
      w1sq=w1*w1
      offset=17.60844567*r(5)
c
c.....r(5)*17.61 = field sweep in 1e6 radHz
c.....r(5)*[ 2.pi.g.bohr/planck ], g = 2.002322
c
      offset=offset*offset
      t1=1.0/r(1)
      ct2=cmplx(t2,0.0)
      ar(3)=t1+t2+t2
      ar(2)=t2*t2+2.0*t1*t2+w1sq+offset
      ar(1)=t1*t2*t2+w1sq*t2+offset*t1
      cw1=cmplx(w1,0.0)
      tinc=cmplx(r(6)/(float(npt)-1.0),0.0)
c
c     call cubrt(ar,s,ind)
c     call denom(s,ss,ind)
      call cubrt(ar,rs,ind)
      call denom(rs,rss,ind)
c
c.....a(k)  = w1*(s(k)+t2)*ss(k)
c.....ss(k) = 1.0/[(s(1)-s(2))*(s(1)-s(3))]
c.....as(k) = [1.0/t1*s(k)]*a(k)
c.....feq   = -(sum over the 3 as(i) terms)
c
      c1 = ( 0.0, 0.0 )
      do 22 k=1,3
        a(k)=cw1*(s(k)+ct2)*ss(k)
        feq=t1/(cabs(s(k))**2.0)
        as(k)=cmplx(feq,0.0)*a(k)*conjg(s(k))
        c1=c1+as(k)
22    continue
      feq=-real(c1)
c
      ctime=tinc
      do 25 k=nzero,npt
        c0=(0.0,0.0)
        c1=(0.0,0.0)
c
c.....c0  = sum over the 3 terms of the type a(i)exp(s(i)t)
c.....c1  = sum over the 3 terms of the type [a(i)/t1*s(i)]exp(s(i)t)
c.....feq  = sum over the 3 terms of the type a(i)/[t1*s(i)]
c
        do 26 l=1,3
          c2=cexp(ctime*s(l))
          c0=c0+a(l)*c2
26        c1=c1+as(l)*c2
c
c.....My(t) = [p(i)/peq]*a(i)+[a(i)/t1*s(i)](exp(s(i)t)-1)
c.....         summation signs having been assumed
c
        data(k,1)=r(3)*real(c0)+peq*(real(c1)+feq)
c
c.....data(k)=r(3)*real(c0)+peq*(real(c1)+feq)
c.....total signal
c
c.....data(k)=real(c1)+feq
c.....signal in absence of initial term
c
c.....data(k)=r(3)*real(c0)
c.....signal in absence of eqm term
c
c.....print *,data(k),real(ctime),k
c
        ctime=ctime+tinc
25    continue
      return
      end
c     ------------------------------------------------------------------
      subroutine denom(rcub,rs,indcub)
c
c.....calculate rs(), input the roots rcub() and indcub(= 0 if any of
c.....the roots coincide).
c
c.....(rs(1),rs(2)) = 1/[(s(1)-s(2))*(s(1)-s(3))]
c
c.....the other elements of rs() are obtained by cyclic permutation of
c.....the eigenvalue subscripts.
c
c.....unlike the QUARTIC routine this is all in single precision.
c
      dimension rs(6),rcub(6)
      if(indcub.eq.0) goto 200
      do 6 k=1,3
        c0=1.00
        c1=0.00
        k1=(k-1)*2+1
        k2=k1+1
        do 7 l=1,2
          kl=k+l
          if(kl.gt.3) kl=kl-3
          l1=(kl-1)*2+1
          l2=l1+1
          a0=rcub(k1)-rcub(l1)
          a1=rcub(k2)-rcub(l2)
          c2=c0*a0-c1*a1
          c3=c0*a1+c1*a0
          c0=c2
          c1=c3
7       continue
        a2=1.00/(c0*c0+c1*c1)
        rs(k1)=a2*c0
        rs(k2)=-a2*c1
6     continue
      return
200   print 300
      return
300   format(/,1X,'STOP-at least one pair of identical roots')
      end
c     -----------------------------------------------------------------
      subroutine cubrt(a,rr,ind)
c
c.....Routine to solve for roots of a cubic - in single precision.
c
c.....Input data are as follows:
c.....         coefficient of s cubed   = 1.0
c.....             ''      of s squared = a(3)
c.....             ''      of s         = a(2)
c.....            constant term         = a(1)
c
c.....Output array s(6)  <=rr(6) is complex
c
      dimension a(3),rr(6)
c.....print *,rr
      c1=1.0/3.0
      c2=sqrt(0.75)
      do 1,j=1,6
1       rr(j)=0.0
      a0=(a(1))
      a1=(a(2))
      a2=(a(3))
      z1=a2*c1
      q=(a1*c1)-(z1*z1)
      r=((a1*a2-3.0*a0)*c1*0.5)-z1*z1*z1
      z3=r*r
      z1=q*q*q+z3
      if (z1) 100,200,300
c
100   ind=-1
      a0=sqrt(z3+abs(z1))
      if(a0.eq.0.0) goto 10
      a0=alog(abs(a0))
      a0=c1*a0
      a0=exp(a0)
10    a1=sqrt(-z1)
      a1=atan2(a1,r)
      a1=a1*c1
      z1=a0*cos(a1)
      z2=a0*sin(a1)
      a2=-c1*a2
      z3=2.00*c2*z2
      rr(3)=(-z1-z3+a2)
      rr(5)=(-z1+z3+a2)
      rr(1)=-a(3)-rr(3)-rr(5)
      return
c
200   ind=0
      if(r.eq.0.0) goto 2
      s1=alog(abs(r))
      s1=c1*s1
      s1=exp(s1)
      if(r.lt.0.0) s1=-s1
2     s1=r
      a2=-a2*c1
      rr(3)=(-s1+a2)
      rr(5)=rr(3)
      rr(1)=-a(3)-2.0*rr(3)
      return
c
300   ind=+1
      z1=sqrt(z1)
      s1=r+z1
      a0=0.0
      if(s1.eq.0.0) goto 3
      a0=alog(abs(s1))
      a0=c1*a0
      a0=exp(a0)
      if(s1.lt.0.0) a0=-a0
3     s1=a0
      s2=r-z1
      a0=0.0
      if(s2.eq.0.0) goto 4
      a0=alog(abs(s2))
      a0=a0*c1
      a0=exp(a0)
      if(s2.lt.0.0) a0=-a0
4     s2=a0
      a2=-a2*c1
      rr(3)=(-0.5*(s1+s2)+a2)
      rr(4)=(c2*(s1-s2))
      rr(5)=rr(3)
      rr(6)=-rr(4)
      rr(1)=-a(3)-2*rr(3)
c.....print *,rr
      return
      end
c     ------------------------------------------------------------------
      subroutine input(r,npt,peq)
c
c.....Purpose : To input the lineshape parameters
c
      dimension r(7)
      character esc*1,sb*1,csi*2,pfa*1
c
      esc=char(27)
      sb=char(91)
      csi=esc//sb
c
100   print 7
      call vtc(1,1)
      call clear
      print 1,csi,csi,(r(i),i=1,7),npt
c
      do 20 j=1,7
        call vtc(j+3,38)
        print 6
        call vtc(j+3,38)
        read(*,2,err=3)fa
        if(fa.eq.0)go to 21
        r(j)=fa
21      print 7
        call vtc(j+3,38)
        print 6
        call vtc(j+3,38)
        print 8,r(j)
20    continue
      go to 101
c
3     r(j)=0.0
101   continue
      call vtc(11,38)
      print 6
      call vtc(11,38)
      read(*,12)ifa
      if(ifa.ne.0) npt=ifa
      call vtc(11,38)
      print 11,npt
      if(r(3).ne.1.0) goto 150
      print 501
      print 500
      read 550,pfa
      if(pfa.eq.'E'.or.pfa.eq.'e') peq=0.0
      if(pfa.eq.'I'.or.pfa.eq.'i') peq=1.0
      if(pfa.eq.'I'.or.pfa.eq.'i') r(3)=0.0
150   print 99
      return
c
500   format(2X,'Want to exclude the equilibrium(E) or',1X,
     &          'initial(I) term?(N) ',\)
501   format(/)
550   format(1A1)
1     format(A2,'32m',31X,'Mode 1 simulation',A2,'0m',/,
     &  //,3X,'Spin-lattice relaxation time, T1',T38,': ',
     &  G12.6,T52,' microseconds',/,3X,'Spin-spin relaxation time,
     &  T2',T38,': ',G12.6,T57,'"',/
     &  ,3X,'P(init)/P(eq)',T38,': ',
     &  G12.6,/,3X,'Omega 1, w1',T38,': ',
     &  G12.6,T52,' radian megahertz',/,3X,'Offset from resonance',T38,
     &  ': ',G12.6,T52,' Gauss',/,3X,'Timespan ',T38,': ',
     &  G12.6,T52,' microseconds',/,3X,'Delay after flash',T38,': ',
     &  G12.6,T52,' microseconds',/
     &  ,3X,'No. of points ',T38,':     ',I4,'      ')
2     format(G12.0)
6     format ('            ',\)
7     format(1X,\)
8     format(1X,G12.6,\)
11    format(5X,I4,\)
12    format(I4)
99    format(/)
c
      end
c     ------------------------------------------------------------------
      subroutine mcalc(r,npt,peq)
c
      common/shape/data
c
      dimension data(2048,2)
      dimension r(7),ar(7),rs(6),rss(6)
      complex a(3),as(3),s(3),ss(3),c0,c1,c2,ct2,cw1,ctime,tinc
c
      equivalence ( rs, s ), ( rss, ss )
c
c.....cw1   = cmplx(r(3),0.0)
c.....ctime = cmplx(r(6),0.0)
c.....tinc  = cmplx(r(7),0.0)        time between samples
c.....t1    = 1.0/t1                 inverse of T1
c.....ct2   = cmplx(1.0/t2,0.0)      complex form of T2
c
      nzero=npt/10
      mcal=npt-nzero+1
      idelay=int(((float(npt))/r(6))*r(7))
      nzero=nzero+idelay
      ncalc=npt-nzero+1
      sfc=exp(-r(7)/r(1))
c
      do 10 i=1,nzero
10      data(i,2)=0.0
c
      t2=1.0/r(2)
      w1=r(4)
      w1sq=w1*w1
      offset=17.60844567*r(5)
c
c.....r(5)*17.61 = field sweep in 1e6 radHz
c.....r(5)*[ 2.pi.g.bohr/planck ], g = 2.002322
c
      offset=offset*offset
      t1=1.0/r(1)
      ct2=cmplx(t2,0.0)
      ar(3)=t1+t2+t2
      ar(2)=t2*t2+2.0*t1*t2+w1sq+offset
      ar(1)=t1*t2*t2+w1sq*t2+offset*t1
      cw1=cmplx(w1,0.0)
      tinc=cmplx(r(6)/(float(npt)-1.0),0.0)
c
c     call cubrt(ar,s,ind)
c     call denom(s,ss,ind)
      call cubrt(ar,rs,ind)
      call denom(rs,rss,ind)
c
c.....a(k)  = w1*(s(k)+t2)*ss(k)
c.....ss(k) = 1.0/[(s(1)-s(2))*(s(1)-s(3))]
c.....as(k) = [1.0/t1*s(k)]*a(k)
c.....feq   = -(sum over the 3 as(i) terms)
c
      c1=(0.0,0.0)
      do 22 k=1,3
        a(k)=cw1*(s(k)+ct2)*ss(k)
        feq=t1/(cabs(s(k))**2.0)
        as(k)=cmplx(feq,0.0)*a(k)*conjg(s(k))
        c1=c1+as(k)
22    continue
      feq=-real(c1)
c
      ctime=tinc
      do 25 k=nzero,npt
        c0=(0.0,0.0)
        c1=(0.0,0.0)
c
c.....c0   = sum over the 3 terms of the type a(i)exp(s(i)t)
c.....c1   = sum over the 3 terms of the type [a(i)/t1*s(i)]exp(s(i)t)
c.....feq  = sum over the 3 terms of the type a(i)/[t1*s(i)]
c
        do 26 l=1,3
          c2=cexp(ctime*s(l))
          c0=c0+a(l)*c2
26        c1=c1+as(l)*c2
c
c.....My(t) = ([(p(i)/peq)-1]*exp(-tau/t1)+1)*a(i))
c.....      + [a(i)/t1*s(i)](exp(s(i)t)-1)
c
c.....              summation signs having been assumed
c
        data(k,2)=(sfc*(r(3)-1.0)+1.0)*real(c0) + peq*(real(c1)+feq)
c
c.....data(k)=r(3)*real(c0)+peq*(real(c1)+feq)
c.....total signal
c
c.....data(k)=real(c1)+feq
c.....signal in absence of initial term
c
c.....data(k)=r(3)*real(c0)
c.....signal in absence of eqm term
c
c.....print *,zdata(k),real(ctime),k
c
        ctime=ctime+tinc
25    continue
      return
      end
c.....-----------------------------------------------------------------
      subroutine splot( data, npts, nspec )
c
      dimension data(2048,2)
      character esc*1
c
      esc = char( 27 )
c
      call gmode
      call clrscr
c
c.....Set display colour and no line wrapping.
c
      print 10,esc,esc,esc
10    format(' ',a1,'[25;10H',a1,'[7l',a1,'[37;1m',
     &            'Mode 1 simulation',\)
c
      call scrnbox
      call scrnplot ( data, npts, nspec, .TRUE. )
      call scrncurs ( data, npts, keyin )
      call tmode
c
c.....Reset colour and turn line wrapping back on.
c
      print 20, esc, esc
20    format (a1,'[7h',a1,'[37;0m')
c
      return
      end

