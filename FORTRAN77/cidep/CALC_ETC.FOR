$LARGE
      subroutine feqm
c
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /factor/g(2),gg(2),fill2(2)
      common /t/rfd,rep
      common /fe/sorbu(2048)
c
      dimension vu(1000)
      character rfd*1,rep*1
c
      nli=nlines(1)
      do 100 lx=1,nli
100     vu(lx)=float(itotal(lx,1))
      do 200 j=1,npt
200     sorbu(j)=0.0
c
      finc=swp/float(npt-1)
c.....field increment
      var1=-finc*0.5
      bf=714.4844244*freq/g(1)
c.....resonant field/Gauss, freq/GHz
      if(rfd.eq.'N') goto 915
      stfld=714.4844244*freq/g(1)
c.....move centre field to resonant
915   f=stfld-bf-(swp/2.0)
c.....field of radical 1
c
      do 920 j=1,npt
        do 24 lx=1,nli
          if(abs(f+sumami(lx,1))+var1.gt.0.0) goto 24
          sorbu(j)=sorbu(j)+vu(lx)
24      continue
920   f=f+finc
c
      return
      end
c
      subroutine denom(rcub,rs,indcub)
c
c     calculate rs(), input the roots rcub() and indcub(= 0 if any of
c     the roots coincide).
c
c     (rs(1),rs(2)) = 1/[(s(1)-s(2))*(s(1)-s(3))]
c
c     the other elements of rs() are obtained by cyclic permutation of
c     the eigenvalue subscripts.
c     Unlike the QUARTIC routine this is all in single precision.
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
c
      return
c
200   print 300
      return
c
300   format(' ',/,1X,'STOP-at least one pair of identical roots')
c
      end
c     ------------------------------------------------------------------
      subroutine cubrt(a,rr,ind)
c
c......Routine to solve for roots of a cubic - in single precision.
c......Input data are as follows:
c......         coefficient of s cubed   = 1.0
c......             ''      of s squared = a(3)
c......             ''      of s         = a(2)
c......            constant term         = a(1)
c......Output array s(6)  <=rr(6) is complex
c
      dimension a(3),rr(6)
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
c
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
c
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
c
      return
      end
c     ------------------------------------------------------------------
      subroutine nozero
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /sense/specp(4096),spece(4096),deg(1000),pol(1000),
     &              ipos(1000),shapep(1024),shapee(1024),swp,stfld,freq1
c
      complex chan
c
      nlines=0
      do 925 i=1,nstick
        if(spece(i).eq.0.0) goto 925
        nlines=nlines+1
        deg(nlines)=spece(i)
        pol(nlines)=specp(i)
        ipos(nlines)=i
925   continue
c
      return
      end
c     ------------------------------------------------------------------
      subroutine const(r)
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /sense/specp(4096),spece(4096),deg(1000),pol(1000),
     &              ipos(1000),shapep(1024),shapee(1024),swp,stfld,freq1
      common /area/finc,t1,ct2,arr(3),cw1,ctime,tinc
c
      dimension r(5)
      complex ct2,chan,ctime,tinc,cw1
c
      t2=1.0/r(2)
      w1=r(3)
      w1sq=w1*w1
c
      nhalf=nshape/2
      nspec=nshape+nstick
      finc=17.60844567*swp/float(nstick-1)
c
      t1=1.0/r(1)
      ct2=cmplx(t2,0.0)
      arr(3)=t1+t2+t2
      arr(2)=t2*t2+2.0*t1*t2+w1sq
      arr(1)=t1*t2*t2+w1sq*t2
      cw1=cmplx(w1,0.0)
      ctime=cmplx(r(4),0.0)
      tinc=cmplx(r(5),0.0)
      fchan=float(nchan)
      chan=cmplx(fchan,0.0)
c
      return
      end
c     ------------------------------------------------------------------
      subroutine loop
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /sense/specp(4096),spece(4096),deg(1000),pol(1000),
     &              ipos(1000),shapep(1024),shapee(1024),swp,stfld,freq1
      common /area/finc,t1,ct2,arr(3),cw1,ctime,tinc
c
      complex a(3),as(3),s(3),c0,c1,ct2,cw1,ctime,tinc,chan
c
      delta=0.0
      k=nhalf
      do 30 i=1,nhalf
        call calc(a,as,s,feq,delta)
        c0=(0.0,0.0)
        c1=(0.0,0.0)
c
        do 26 l=1,3
          a(l)=a(l)*s(l)
          as(l)=as(l)*s(l)
          c0=c0+a(l)
26        c1=c1+as(l)
c
        shapep(k)=real(c0)
        shapee(k)=real(c1)+feq
c
        k=k-1
30      delta=delta+finc
c
      return
      end
c     ------------------------------------------------------------------
      subroutine calc(a,as,s,feq,delta)
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /area/finc,t1,ct2,arr(3),cw1,ctime,tinc
c
      dimension r(5),ar(3)
      complex a(3),as(3),s(3),ss(3),c1,c2,ct2,cw1,ctime,tinc,chan
c
      del=delta*delta
      ar(3)=arr(3)
      ar(2)=arr(2)+del
      ar(1)=arr(1)+del*t1
      call cubrt(ar,s,ind)
      call denom(s,ss,ind)
c
      c1=(0.0,0.0)
      do 22 k=1,3
        a(k)=cw1*(s(k)+ct2)*ss(k)
        feq=t1/(cabs(s(k))**2.0)
        as(k)=cmplx(feq,0.0)*a(k)*conjg(s(k))
        c1=c1+as(k)
22    continue
      feq=-real(c1)
      if(nchan-1) 23,23,30
c
23    do 24 k=1,3
24      s(k)=cexp(ctime*s(k))
      goto 100
c
30    feq=feq*fchan
c
45    do 50 k=1,3
        c1=(1.0,0.0)-cexp(tinc*s(k))
        s(k)=cexp(ctime*s(k))*conjg(c1)*((1.0,0.0)
     &      -cexp(chan*tinc*s(k)))/(cabs(c1)**2)
50    continue
c
100   return
      end
c     ------------------------------------------------------------------
      subroutine convol
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /sense/specp(4096),spece(4096),deg(1000),pol(1000),
     &              ipos(1000),shapep(1024),shapee(1024),swp,stfld,freq1
c
      complex chan
c
c      do 10 i=1,nspec
c        specp(i)=0.0
c10      spece(i)=0.0
c      do 30 j=1,nlines
c        k1=ipos(j)
c        k2=nshape+k1
c        do 20 i=1,nhalf
c          k2=k2-1
c          sp=shapep(i)*pol(j)
c          specp(k2)=specp(k2)+sp
c          specp(k1)=specp(k1)+sp
c          sp=shapee(i)*deg(j)
c          spece(k2)=spece(k2)+sp
c          spece(k1)=spece(k1)+sp
c20        k1=k1+1
c30    continue
c
      do 10 i=1,nspec
        specp(i)=0.0
10      spece(i)=0.0
      do 30 j=1,nlines
        k1=ipos(j)
        k2=nshape+k1
        do 20 i=1,nhalf
          k2=k2-1
          sp=shapep(i)*pol(j)
          specp(k2)=specp(k2)+sp
          specp(k1)=specp(k1)+sp
          sp=shapee(i)*deg(j)
          spece(k2)=spece(k2)+sp
          spece(k1)=spece(k1)+sp
20        k1=k1+1
30    continue
c
      return
      end
c     ------------------------------------------------------------------
      subroutine group
c
c Purpose : to separate the lines of the stick spectrum into sub-spectra of equal degenaracy
c
c IPOS  :array of line positions in groups of equal degeneracy, in order of equal degeneracy
c
c IWORK :a work array
c
c NLINES:array of the number of lines in each group
c
c NGROUP:number of different degeneracies
c
c SIZE:array of degeneracies in order of increasing degeneracy
c
c SORB:array containing the NSTICK-point stick spectrum
c
c
c
      common/field/swp,stfld,freq,finc
      common/nliste/nstick,nspec,nhalf,nshape,nchan,fchan,
     &              chan,nc,ngroup
      common/exch/sorb(2048),fres(2048),awt(2048),fill(130)
      common/exch1/spec(2560),sub(2560),size(20),ipos(100),
     &             iwork(100),nlines(20),shape(1024)
c
      complex chan
c
c.....remove the 0'sfrom the stick spectrum and count the number of
c.....lines (NL) nl=0
c
      do 925,j=1,nstick
        if(sorb(j).eq.0.0) go to 925
        nl=nl+1
        sorb(nl)=sorb(j)
        iwork(nl)=j
925   continue
c
      smax=0.0
      do 955,j=1,nl
955     if(smax.lt.sorb(j)) smax=sorb(j)
c
      tmin=0.0
      ngroup=0
950   smin=smax
c
      do 935,j=1,nl
935     if(smin.ge.sorb(j).and.sorb(j).gt.tmin) smin=sorb(j)
c
      ngroup=ngroup+1
      size(ngroup)=smin
      tmin=smin
      if(smin.lt.smax) goto 950
c
      l=1
      do 970,i=1,ngroup
        nlines(i)=0
        si=size(i)
        do 960,j=1,nl
          if(sorb(j).ne.si) go to 960
          nlines(i)=nlines(i)+1
          ipos(l)=iwork(j)
          l=l+1
960     continue
970   continue
c
      return
      end
c     ------------------------------------------------------------------
      subroutine conste(r,fdegi)
c
      common/exch2/q(4),x,y,t1,calpha,cw1,gamma2,ctime,tinc
      common/field/swp,stfld,freq,finc
      common/nliste/nstick,nspec,nhalf,nshape,nchan,fchan,
     &              chan,nc,ngroup
c
      dimension r(7)
      complex calpha,cw1,ct1,gamma2,ctime,tinc,chan
c
      t1=1.0/r(1)
      t2=1.0/r(2)
      w3=1.0/r(3)
      if((100.0*r(1)).lt.r(3)) w3=0.0
      w2=(1.0-fdegi)*w3
      w1=r(5)
      alpha=t2+w2
      beta=t1+w2
      gamma=t1+w3-w2
      delta=-w2*(w3-w2)
      x=beta+gamma
      y=beta*gamma+delta
      q(4)=x+2.0*alpha
      q(3)=y+alpha*alpha+2.*alpha*x+w1*w1
      q(2)=2.0*alpha*y+x*alpha*alpha+w1*w1*(alpha+gamma)
      q(1)=alpha*alpha*y+w1*w1*alpha*gamma
      cw1=cmplx(w1,0.0)
      gamma2=cmplx(t1+w3,0.0)
      ctime=cmplx(r(6),0.0)
      tinc=cmplx(r(7),0.0)
      calpha=cmplx(alpha,0.0)
      fchan=float(nchan)
      chan=cmplx(fchan,0.0)
      nc=nchan-1
c
      return
      end
c     ------------------------------------------------------------------
      subroutine loope(r)
c
c.....Purpose : To calculate the low field half of the lineshape
c
c.....Input : (a) contents of AREA
c.....        (b) r(4) - polarisation
c.....        (c) swp,nstick,nshape
c
c.....Output : SHAPE() - array containing half the lineshape
c
c.....Subroutine used : CALC to solve the quartic at each field
c.....position
c
      common/exch2/q(4),x,y,t1,calpha,cw1,gamma2,ctime,tinc
      common/field/swp,stfld,freq,finc
      common/nliste/nstick,nspec,nhalf,nshape,nchan,fchan,
     &              chan,nc,ngroup
      common/exch/sorb(2048),fres(2048),awt(2048),fill(130)
      common/exch1/spec(2560),sub(2560),size(20),ipos(100),
     &             iwork(100),nlines(20),shape(1024)
      dimension r(7)
      complex a(4),as(4),s(4),c0,c1,calpha,cw1,gamma2,
     &        ctime,tinc,chan
c
      delta=0.0
      mpt=nshape/2
      k=mpt
      do 30 i=1,mpt
        call calce(a,as,s,feq,delta)
        c0=(0.0,0.0)
        c1=(0.0,0.0)
c
c.....c0 = sum over the 4 terms of the type a(i)exp(s(i)t)
c.....c1 = sum over the 4 terms of the type a(i)/[t1*s(i)exp(s(i)t]
c
        do 26 l=1,4
          a(l)=a(l)*s(l)
          as(l)=as(l)*s(l)
          c0=c0+a(l)
26        c1=c1+as(l)
c
c.....My(t) = [p(i)/p(eq)*a(i) + [a(i)/t1*s(i)](exp(s(i)t) - 1
c.....the summations signs have been assumed
c
        shape(k)=r(4)*real(c0)+real(c1)+feq
c.....total signal
c.....  shape(k)=real(c0)
c.....                ;%signal in absence of eqm term
c.....  shape(mpt+i)=shape(k) ;%high field = low field side
        k=k-1
30      delta=delta+finc
c.....INCREASE the offset
c
      return
      end
c     ------------------------------------------------------------------
      subroutine convole(j,istart)
c
c.....Purpose : To calculate the convolution
c
c.....SHAPE -array of dimension NSHAPE , symmetrical ie.
c.....      SHAPE(I) = SHAPE(NA-I+1). Thus only the 1st half of the
c.....      array is required ; NHALF = NSHAPE/2
c
c.....B     -array of dimension NSTICK, B() consists of  0's. The
c.....      other elements = SJ
c.....      Eg.  B = (0 SJ 0 0 0 SJ 0 0 0 0), NSTICK = 10.
c.....      B() itself is NOT REQUIRED
c
c.....SUB    -array of dimension NSPEC = NSHAPE + NSTICK -1
c
c.....IPOS   -array of dimension NL(no. of lines in this sub-spectrum).
c.....        IPOS is the array containing the positions of the non-zero
c.....        values of B(), eg. in the above case IPOS = (2 6) and
c.....        NL = 2. The lines of the ith sub-spectrum start at
c.....        ipos(j), j=Sum{nlines(k), k=1,j-1}.
c
c.....SPEC    -array of dimension NSPEC, is the sum of arrays SUB()
c.....        due to many such convolutions
c
      common/exch2/q(4),x,y,t1,calpha,cw1,gamma2,ctime,tinc
      common/field/swp,stfld,freq,finc
      common/nliste/nstick,nspec,nhalf,nshape,nchan,fchan,
     &              chan,nc,ngroup
      common/exch/sorb(2048),fres(2048),awt(2048),fill(130)
      common/exch1/spec(2560),sub(2560),size(20),ipos(100),
     &             iwork(100),nlines(20),shape(1024)
c
        complex calpha,cw1,gamma2,ctime,tinc,chan
c
      do 10 i=1,nspec
10      sub(i)=0.0
      nl=nlines(j)
      kk=istart
      do 30 k=1,nl
        k1=ipos(kk)
        k2=nshape+k1
        do 20 i=1,nhalf
          k2=k2-1
          sub(k2)=sub(k2)+shape(i)
          sub(k1)=sub(k1)+shape(i)
20        k1=k1+1
30    kk=kk+1
c.....Add this (the jth) sub-spectrum to the total spectrum
      sj=size(j)
c.....degeneracy of the lines in this sub-spectrum
      do 40 i=1,nspec
40      spec(i)=spec(i)+sj*sub(i)
c
      return
      end
c     ------------------------------------------------------------------
      subroutine calce(a,as,s,feq,delta)
c
c.....Input : Contents of AREA, delta, nchan, fchan, chan
c
c.....Output : a(),as(),s(),feq
c
      common/exch2/q(4),x,y,t1,calpha,cw1,gamma2,ctime,tinc
      common/field/swp,stfld,freq,finc
      common/nliste/nstick,nspec,nhalf,nshape,nchan,fchan,
     &              chan,nc,ngroup
      common/exch/sorb(2048),fres(2048),awt(2048),fill(130)
      common/exch1/spec(2560),sub(2560),size(20),ipos(100),
     &             iwork(100),nlines(20),shape(1024)
c
      complex a(4),as(4),s(4),ss(4),c0,c1,calpha,cw1,gamma2,
     &        ctime,geo,sumint,tinc,chan
        dimension quad(4)
c
      del=delta*delta
      quad(4)=q(4)
      quad(3)=q(3)+del
      quad(2)=q(2)+x*del
      quad(1)=q(1)+y*del
      call quart(quad,s,ss)
c
c.....a(k) =w1*(s(k)+t2+w2)*(s(k)+t2+w3)/ss(k)
c.....ss(k)=1.0/(s(1)-s(2))*(s(2)-s(3))*(s(3)-s(4))*(s(4)-s(1))
c.....as(k)=[t(1)/s(k)]*a(k)
c.....feq  =-(sum over the 4as(i) terms)
c
      c1=(0.0,0.0)
      do 22 k=1,4
        a(k)=cw1*(s(k)+calpha)*(s(k)+gamma2)*ss(k)
        feq=t1/(cabs(s(k))**2.0)
        as(k)=cmplx(feq,0.0)*a(k)*conjg(s(k))
        c1=c1+as(k)
22    continue
      feq=-real(c1)
c.....integrate ?
      if(nchan-1) 23,23,30
c.....Don't "integrate"
23    do 24 k=1,4
24      s(k)=cexp(ctime*s(k))
      goto 100
c
c....."integrate" by expressing the sum in terms of a geometric series
c.....2 methods for calculating the sum of this series :
c.....                (a) numerical -lines 35 - 40
c.....                (b) analytical-lines 45 - 50
c.....The accuracy for the 2 methods is the same (within the limits of
c.....the error due to the use of single precision arithmetic) while the
c.....second is clearly the faster.
c
30      feq=feq*fchan
c35   do 40 k=1,4
c.....sumint=cexp(ctime*s(k))
c.....geo=ctime
c.....do 37 l=1,nc   ;%nc=nchan-1
c.....geo=geo+tinc
c37   sumint=cexp(geo*s(k))+sumint
c40   s(k)=sumint
c.....go to 100
c use the analytical expression for the sum of a series
45    do 50 k=1,4
      sumint=(1.0,0.0)-cexp(tinc*s(k))
      s(k)=cexp(ctime*s(k))*conjg(sumint)*
     &     ((1.0,0.0)-cexp(chan*tinc*s(k)))/(cabs(sumint)**2)
50    continue
c
100   return
      end
c     ------------------------------------------------------------------
      subroutine quart(a,rr,rs)
c
c.....routine used for solving slow exchange quartic. July 1984
c
      dimension a(4),rr(8),rs(8)
      double precision c0,c1,c2,c3,c4,c5,a0,a1,a2,a3,rcub(8)
      double precision acub(3),c(2),try1,try2
c
      do 100 j=1,8
100     rr(j)=0.0
      a0=dble(a(1))
      a1=dble(a(2))
      a2=dble(a(3))
      a3=dble(a(4))
      acub(3)=-a2
      acub(2)=(a1*a3-4.0D00*a0)
      acub(1)=-(a1*a1+a0*a3*a3-4.0D00*a0*a2)
      call dcubrt(acub,rcub,indcub)
      i=-1
1     i=i+2
      if(i.gt.6) stop 'no solns'
      c5=rcub(i)
      c0=a3*a3*0.25D00+c5-a2
      c1=c5*c5*0.25D00-a0
      if(indcub.eq.1)goto 2
      if(c0.lt.0.00D0.or.c1.lt.0.D00) goto 1
c
2     if(c0.lt.0.00D0.or.c1.lt.0.D00) print *,' negative!',c0,c1
      c3=dsqrt(c0)
      c4=dsqrt(c1)
      i=1
      c1=a3*0.5d00
      c2=c5*0.5d00
      try1=(c1-c3)*(c2+c4)+(c1+c3)*(c2-c4)
      try2=(c1-c3)*(c2-c4)+(c1+c3)*(c2+c4)
      if(dabs(try1-a1).lt.dabs(try2-a1))goto 23
      c3=-c3
23    c(1)=c1+c3
      c(2)=c2+c4
3     call quad(c,rcub(i),indq)
      if(i.ne.1)goto 5
      i=5
      c(1)=c1-c3
      c(2)=c2-c4
      go to 3
5     do 10,i=1,8
10    rr(i)=real(rcub(i))
      call vector(rcub,rs,indq)
c
      return
      end
c     ------------------------------------------------------------------
      subroutine vector(rcub,rs,indq)
c
      dimension rs(8)
      double precision c0,c1,c2,c3,a0,a1,a2,rcub(8)
c
      if(indq.eq.0)go to 200
      do 6 k=1,4
        c0=1.00d00
        c1=0.00d00
        k1=(k-1)*2+1
        k2=k1+1
        do 7 l=1,3
          kl=k+l
          if(kl.gt.4) kl=kl-4
          l1=(kl-1)*2+1
          l2=l1+1
          a0=rcub(k1)-rcub(l1)
          a1=rcub(k2)-rcub(l2)
          c2=c0*a0-c1*a1
          c3=c0*a1+c1*a0
          c0=c2
          c1=c3
7       continue
        a2=1.00d00/(c0*c0+c1*c1)
        rs(k1)=sngl(a2*c0)
        rs(k2)=sngl(-a2*c1)
6     continue
c
      return
c
200   print 300
c
      return
c
300   format(' ',/,1X,'STOP-at least one pair of identical roots   ')
c
      end
c     ------------------------------------------------------------------
      subroutine dcubrt(a,rr,ind)
c
      double precision r,q,s1,s2,z1,z2,z3,a0,a1,a2,c1,c2,rr(6),a(3)
c
      c1=1.0d00/3.0d00
      c2=dsqrt(0.75d00)
      a0=dble(a(1))
      a1=dble(a(2))
      a2=dble(a(3))
      z1=a2*c1
      do 1 j=1,5,2
        rr(j)=-z1
1       rr(j+1)=0.0d00
      q=(a1*c1)-(z1*z1)
      r=(a1*a2*c1-a0)*0.5d0-z1*z1*z1
      z3=r*r
      z1=q*q*q+z3
      if(z1)100,200,300
100   ind=-1
      a0=dsqrt(z3+dabs(z1))
      if(a0.eq.0.0d00)goto 10
      a0=a0**(c1)
10    a1=dsqrt(-z1)
      a1=datan2(a1,r)
      a1=a1*c1
      z1=a0*dcos(a1)
      z2=a0*dsin(a1)
      rr(1)=(2.0D00*z1)+rr(1)
      z3=2.00d00*c2*z2
      rr(3)=(-z1-z3)+rr(3)
      rr(5)=(-z1+z3)+rr(5)
      return
c
200   ind=0
      if(r.eq.0.0d00)goto 2
      s1=dabs(r)**(c1)
      if(r.lt.0.0d00)s1=-s1
2     s1=r
      rr(1)=(2.0d00*s1)+rr(1)
      rr(3)=-s1+rr(3)
      rr(5)=rr(3)
      return
c
300   ind=+1
      z1=dsqrt(z1)
      s1=r+z1
      a0=0.0d00
      if(s1.eq.0.0d00)goto 3
      a0=dabs(s1)**(c1)
      if(s1.lt.0.0d00)a0=-a0
3     s1=a0
      s2=r-z1
      a0=0.0d00
      if(s2.eq.0.0d00)goto 4
      a0=dabs(s2)**(c1)
      if(s2.lt.0.0d00)a0=-a0
4     s2=a0
      rr(1)=s1+s2+rr(1)
      return
c
      end
c     ------------------------------------------------------------------
      subroutine quad(a,s,ind)
c
      double precision a(2),s(4),b,det
c
      b=-a(1)*0.5d0
      do 5 i=1,3,2
        s(i)=b
5       s(i+1)=0.0d0
      det=b*b-a(2)
      ind=int(sngl((dsign(1.0d0,det))))
      if(ind) 10,50,30
10    det=sqrt(-det)
      s(2)=det
      s(4)=-det
      return
c
30    det=sqrt(det)
      s(1)=s(1)+det
      s(3)=s(3)-det
50    return
c
      end

