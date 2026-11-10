$LARGE
      subroutine exchin(r)
c
c.....Purpose : To input the lineshape parameters
c
      common/exch2/q(4),x,y,t1,calpha,cw1,gamma2,ctime,tinc
      common/nliste/nstick,nspec,nhalf,nshape,nchan,fchan,
     &              chan,nc,ngroup
c
      complex calpha,cw1,gamma2,ctime,tinc,chan
      dimension r(7)
c
100   print 7
      call clear
      call vtc(1,1)
      print 1,(r(i),i=1,7),nshape,nchan
      print 7
c
      do 20 j=1,7
        call vtc(j+2,38)
        print 6
        call vtc(j+2,38)
        read(*,2,err=3)fa
        if(fa.eq.0)go to 21
        r(j)=fa
21      print 7
        call vtc(j+2,38)
        print 6
        call vtc(j+2,38)
        print 8,r(j)
20    continue
      go to 101
3     r(j)=0.0
c
101   call vtc(10,38)
      print 6
      call vtc(10,38)
      read(*,12,err=9001) ifa
      if(ifa.ne.0) nshape=ifa
9001  call vtc(10,38)
      print 11,nshape
c
      call vtc(11,38)
      print 6
      call vtc(11,38)
      read(*,12,err=9002)ifa
      if(ifa.ne.0) nchan=ifa
9002  call vtc(11,38)
      print 11,nchan
c
      return
c
1     format(/,3X,'Spin-lattice relaxation time, T1',T38,'  ',
     &  G12.6,T52,' microseconds',/,3X,'Spin-spin relaxation time,
     &  T2',T38,'  ',
     &  G12.6,T57,'"',/,3X,'Electron exchange time, (1./k[R])',T38,'  ',
     &  G12.6,T57,'"',/,3X,'Polarisation ratio, P(i)/Peq',T38,'  ',
     &  G12.6,T57,/,3X,'Microwave power, w1',T38,'  ',
     &  G12.6,T52,' radian megahertz',/,3X,'Time after flash',T38,'  ',
     &  G12.6,T52,' microseconds',/,3X,'Sample interval',T38,'  ',
     &  G12.6,T52,' microseconds',/,3X,'Npts in lineshape',T38,'  ',
     &  I4,'      ',/,3X,'No. of samples ',T38,'  ',I4,'      ')
2     format(G12.0)
6     format(' ',13X,\)
7     format(' ',1X,\)
8     format(' ',1X,G12.6,\)
11    format(' ',1X,I4,\)
12    format(I4)
      end
c     ------------------------------------------------------------------
c
      subroutine spece
c
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
      common / exch2 / q(4),x,y,t1,calpha,cw1,gamma2,ctime,tinc
      common / field / swp,stfld,freq,finc
      common / nliste / nstick,nspec,nhalf,nshape,nchan,fchan,
     &                 chan,nc,ngroup
      common / exch / sorb(2048),fres(2048),awt(2048),fill(130)
      common / exch1 / spec(2560),sub(2560),size(20),ipos(100),
     &                 iwork(100),nlines(20),shape(1024)
c
      complex calpha,cw1,gamma2,ctime,tinc,chan
      dimension fnam (3),r(7),rs(5),ispec(2560)
      character fnam1*1,esc*1,sb*1,csi*2,ans*1
c
c.....set default parameters
c
      r(1)   = 6.7
      r(2)   = 0.66
      r(3)   = 5.5
      r(4)   = -10.0
      r(5)   = 0.1
      r(6)   = 5.0
      r(7)   = 0.01
      nshape = 512
      nchan  = 100
c
      if(mflag(3).ne.1) then
        print 899
        call getscan( keyin, 0 )
        return
899     format(1X,/,'      CREATE T.M. STICK SPECTRUM FIRST     ',/
     &         /,5X,'Hit Any Key to Continue',/)
      endif
c
      esc = char(27)
      sb  = char(91)
      csi = esc//sb
c
c.....set default parameters
c
      r(1)   = 6.7
      r(2)   = 0.66
      r(3)   = 5.5
      r(4)   = -10.0
      r(5)   = 0.1
      r(6)   = 5.0
      r(7)   = 0.01
      nshape = 512
      nchan  = 100
c
      nstick=cnpt
      swp=cswp
      stfld=cstfld
      freq=cfreq1
      do 37 kk=1,nstick
37      sorb(kk)=cstke(kk)
c
      call clear
      print 1,csi,csi
1     format(' ',A2,'32m',
     &  20X,'E.S.R. simulation by direct convolution',/,
     &  1X,A2,'0m',/,8X,'( for systems undergoing slow ',
     &  'degenerate electron exchange )',///)
c
      print 2
2     format(' ',5X,'Want to see stick spectra ?(N) ',\)
      read(*,3,err=9001) ans
3     format(A1)
c.....call regbox
      if(ans.eq.'Y'.or.ans.eq.'y') call splots(sorb,nstick,8)
9001  call group
c.....separate the lines into groups of equal degeneracy
40    call exchin(r)
c.....input lineshape parameters
      nspec = nshape + nstick
c.....no. of points in the final spectrum
      nhalf = nshape / 2
c.....clear array that will contain the final spectrum.
      do 90 i=1,nspec
90      spec(i)=0.0
      istart=1
c.....pointer to the start of the ith sub-spectrum
c.....calculate the total degeneracy
      totdeg=0.0
      do 95 i=1,ngroup
95      totdeg=totdeg+float(nlines(i))*size(i)
      if(totdeg.ne.ctotdeg) then
        print 898
        totdeg=ctotdeg
      end if
898   format(1X,//,5X,'N.B. only part of spectrum')
c.....change field units to rad/MHz
      finc = 17.60844567 * swp / ( float( nstick - 1 ) )
c.....offset * 17.61 = offset in rad/(1.e-6s)
c.....ie. offset*( 2.pi.g.bohr/planck ), g = 2.002322
c
c.....main calculation cycle
      print 150,int(totdeg)
      do 100 j=1,ngroup
        print 200,j,nlines(j),int(size(j))
c.....calculate the lineshape
        call conste(r,(size(j)/totdeg))
        call loope(r)
c.....do the convolution
        call convole(j,istart)
c.....Display the ith sub-spectrum
c.....call splot(sub,nspec,11,nshape,0)
        istart=istart+nlines(j)
c.....move pointer to next sub-spectrum
100   continue
c
c.....Plot the total spectrum
c
      call splot(spec,nspec,11,nshape,0)
c
      call clear
      call vtc(1,1)
      print 560
      read(*,3,err=9998) fnam1
      rs(1)=r(1)
      rs(2)=r(2)
      rs(3)=r(5)
      rs(4)=r(6)
      rs(5)=r(7)
      if(fnam1.eq.'y'.or.fnam1.eq.'Y') call vfile(spec,nspec,nshape,
     &  swp,stfld,freq,rs,nchan)
9998  print 550
      read(*,3,err=9999) fnam1
      if(fnam1.ne.'m'.and.fnam1.ne.'M') goto 40
9999  return
c
150   format(//,1X,'  Total degeneracy = ',I4,//)
200   format(' ',1X,'Group ',I2,' with ',I3,' lines and degeneracy ',
     &  I3)
450   format(3A4)
550   format(' ',' Fresh lineshape (L), or return to Master Program ',
     &  '(M) ?(L) ',\)
560   format(' ',' Do you want this spectrum on file ?(N) ',\)
c
      end

