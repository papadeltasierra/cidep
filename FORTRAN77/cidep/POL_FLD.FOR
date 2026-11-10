$LARGE
        subroutine pol(same,itm,irnc)
c
c.....Purpose : to calculate the radical pair polarisation on each line,
c.....          barring some constant multiplying terms which 'option'
c.....          will put in.
c
c.....Input   : itotal(1024,2), sumami(1024,2)
c
c.....Output  : v(i,k) and itot(k), the total no. of nuclear spin states of
c.....          kth radical.
c
      common /a/v(1024,2),total(1024,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
c
      character same*1,ans*1
c
      gg(1)=stfld*(g(1)-g(2))/2.00232
      gg(2)=-gg(1)
      do 20 i=1,2
        nli=nlines(i)
        if(same.ne.'N') gg(i)=0.0
c.....net effect.
        do 10 j=1,nli
          total(j,i)=float(itotal(j,i))
10        v(j,i)=0.0
20    continue
c
      ans='N'
      if(itm.eq.0) then
        print 22
        read(*,23,err=9005) ans
        goto 9006
9005    ans='N'
9006  end if
      irnc=0
      if(ans.eq.'Y'.or.ans.eq.'y') then
        irnc=1
        print 36,rnc
        read(*,33,err=26) rncD
        if(rncD.ne.0.0) rnc=rncD
        if(abs(rnc).lt.1.0e-6) rnc=1.0E-6
26    end if
      if(ans.ne.'Y'.and.ans.ne.'y'.and.itm.eq.0) then
        print 32,ca
        read(*,33,err=34) caD
        if(caD.ne.0.0) ca=caD
34    end if
22    format(' ',1X,'change Q1/2 dependence to Qn dependence ?(N)',
     &           17X,\)
32    format(' ',1X,'coeff. a in linear term ?(',F6.4,' 1/sqrt(G))',
     &           18X,\)
33    format(F8.0)
36    format(' ',1X,'coeff. n ?(',F6.4,')                        ',
     &           19X,\)
23    format(A1)
      if(ans.eq.'y') ans='Y'
c
      i=1
      ii=2
25    itot(i)=0.0
      nli=nlines(i)
      nlii=nlines(ii)
      do 40 j=1,nli
        sumji=sumami(j,i)+gg(i)
        do 30 k=1,nlii
          q=sumji-sumami(k,ii)
          s=sign(1.0,q)
          if(ans.ne.'Y') then
            v(j,i)=v(j,i)+s*(sqrt(s*q)-ca*(s*q))*total(k,ii)
          else
            v(j,i)=v(j,i)+s*((s*q)**rnc)*total(k,ii)
          end if
30      continue
40      itot(i)=itot(i)+itotal(j,i)
      if(i.eq.1) itot(2)=itot(1)
      if(i.eq.2) goto 50
      i=2
      ii=1
      if(same.eq.'N')  goto 25
50    return
      end
c
      subroutine option(pair,same)
c
c.....Purpose : to scale the polarisations by a constant which is dependant,
c.....          largely, on the theory employed.
c
c.....Input   : nlines(2), v(1024,2), same
c
c.....Output  : vrel(1024,2), v(1024,2)
c
c
c.....OPTION : allows the use of different RPM equations, all of which incorporate
c.....         the square root dependance on Q.--only the one from ref 1 is incorporated, however.
c
c.....       references:
c.....                   1. F.J.Adrian, J.C.P. 54,3918,(1971).
c
c
c.....1Gauss ~ g.beta/h  Hertz ; for g = 2.0023
c
c.....pol = 0.021*sign(Q)*sqrt(tau*Q)/pol(eq)
c
c.....Q(rad) = Q(Gauss)*2.80247e6*2*pi
c.....tau = 10ps
c.....pol(eq) = 0.7e-3
c.....Q = (Bohr*B0*deltag+sumami(radical1)-sumami(radical2))/2
c
c.....The "2" in the above expression is added in at this stage so the calculated
c.....polarisations are scaled by:
c
c.....         0.021*sqrt[(2.8e6*2*pi*10e-12)/2] / 0.7e-3  = 0.281 for F-pairs.
c.....    or   0.085*sqrt[(2.8e6*2*pi*10e-12)/2] / 0.7e-3  = 1.139 for S/T pairs.
c
c
      common /a/v(1024,2),total(1024,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
c
      character same*1,pair*1
c
      scale=-0.281368326
c.....             ;% ref 1, eqn 20, Kr=1.0
      if(pair.eq.'S') scale=1.138871796
c.....ref 1, eqn 17
      if(pair.eq.'T') scale=-1.138871796
c.....ref 1, eqn 17
      sgn=-1.0
      if(pair.eq.'S') sgn=1.0
c
      i=1
12    vmin=abs(v(1,i))
      do 11 j=1,nlines(i)
11      if(vmin.gt.abs(v(j,i)).and.v(j,i).ne.0.0) vmin=abs(v(j,i))
      if(i.eq.2) goto 9
      if(same.eq.'N') then
        i=2
        goto 12
      end if
c
9     i=1
      ii=2
c
10    vmin=1.0/vmin
      do 25 j=1,nlines(i)
        vrel(j,i)=v(j,i)*vmin*sgn
25      v(j,i)=v(j,i)*scale
c
      if(i.eq.2) goto 30
      i=2
      ii=1
      if(same.eq.'N') goto 10
c
c.....check that the sum of the polarisations over the 2 radicals is zero.
c
c30      ss=0.0
c.....   sum=0.0
c.....   do 40 i=1,nlines(1)
c.....   ss=ss+abs(v(l,1))*float(itotal(i,1))
c40      sum=sum+v(i,1)*float(itotal(i,1))
c.....   do 50 i=1,nlines(2)
c.....   ss=ss+abs(v(i,2))*float(itotal(i,2))
c50      sum=sum+v(i,2)*float(itotal(i,2))
c.....   print 100,sum,ss,sum/ss
30      return
c100     format(/,1X,'Sum of pol = ',G12.6,/,1X,'Total = ',G12.6,
c.....&/,1X,'Fractional sum = ',G12.6,/)
c
      end
c
      subroutine sort(same)
c
c.....Purpose : to rearrange the contents of sumami and itotal so that they are in
c.....          order of increasing A*Mi.
c
c
      common /a/v(1024,2),total(1024,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
c
      dimension dv(1024,2),dvrel(1024,2),dsumami(1024,2)
      dimension ndtotal(1024,2)
      character same*1
      integer dj
c
      i=0
100   i=i+1
      nli=nlines(i)
      nlii=nli-1
      do 500 j=1,nlii
        k=j+1
        sumji=sumami(j,i)
        do 400 jj=k,nli
          if(sumji-sumami(jj,i).le.0.0) goto 400
c
          buffs=sumji
          buffv=v(j,i)
          buffr=vrel(j,i)
          ibuffj=itotal(j,i)
c
          sumji=sumami(jj,i)
          v(j,i)=v(jj,i)
          vrel(j,i)=vrel(jj,i)
          itotal(j,i)=itotal(jj,i)
c
          sumami(jj,i) = buffs
          v(jj,i)      = buffv
          vrel(jj,i)   = buffr
          itotal(jj,i) = ibuffj
c
400     continue
        sumami(j,i)=sumji
500   continue
c
      if(i.eq.2) goto 751
      if(same.eq.'N') goto 100
      do 600 j=1,nli
        sumami(j,2)=sumami(j,1)
600     itotal(j,2)=itotal(j,1)
c
751   do 750 i=1,2
        do 750 j=1,nlines(i)
          dv(j,i)=v(j,i)
          dvrel(j,i)=vrel(j,i)
          dsumami(j,i)=sumami(j,i)
          ndtotal(j,i)=itotal(j,i)
750   continue
      do 760 i=1,2
        do 760 j=1,nlines(i)
          dj=nlines(i)-j+1
          v(j,i)=dv(dj,i)
          vrel(j,i)=dvrel(dj,i)
          sumami(j,i)=dsumami(dj,i)
          itotal(j,i)=ndtotal(dj,i)
760   continue
      return
      end

      subroutine field2(i,itm,tm)
c
c.....Purpose : to find the field position of each line.
c
c.....Input   : nlines(2), itot(2), itotal(1024,2), swp, npt, freq, g(2), sumami(1024,2)
c
c.....Output  : sorb(2048), containing the stick spectrum
c
c
c.....nat     : no. of lines for the present radical
c
c.....(Planck's constant / Bohr magneton) = 0.7144844244 e-6 s/Gauss
c
c.....               The polarisations are divided by the total degeneracy of the radical pair
c.....(i.e. the product of the total degeneracies of each radical). This is to take in to
c.....account the actual occurence of each electron spin state in a given ensemble of radicals.
c
c.....N.B. The sign convention +ve sumami = lowfield is used.
c
      common /a/v(1024,2),total(1024,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /s/sorb(2048),fill(352)
      common /t/rfd,rep
      character rfd*1,rfdD*1,rep*1,repD*1
      character same*1
c
      if(i.ne.1) goto 781
      do 780 k=1,npt
780     sorb(k)=0.0
c
781   nli=nlines(i)
      tot12=1.0/(float(itot(1))*float(itot(2)))
      if(istm(i).eq.1) goto 905
c
      if(itm.eq.1) then
        print 9992,tm
        read(*,9993,err=9994) tmD
        if(tmD.ne.0.0) tm=tmD
        goto 9994
      end if
      tm=1.0
      print 90,rep
      read(*,95,err=879) repD
      if(repD.eq.'u') repD='U'
      if(repD.eq.'p') repD='P'
      if(repD.ne.'U'.and.repD.ne.'P') goto 879
      rep=repD
879   if(rep.ne.'U') goto 890
c.....calc an unpolarised spectrum
9994  do 880 lx=1,nli
880     v(lx,i)=tm*float(itotal(lx,i))
      goto 905
c.....calc polarised spectrum
890     do 900 lx=1,nli
900       v(lx,i)=v(lx,i)*tot12*float(itotal(lx,i))
c
905   relc1=0.0
      if(i.eq.1) goto 908
      print 99
      read(*,100,err=9002) relc1D
      if(relc1D.ne.0.0) relc1=relc1D
9002  if(relc1.lt.1.0e-20) relc1=1.0
908   do 910 j=1,npt
910     sorb(j)=sorb(j)*relc1
c
      finc=swp/float(npt-1)
c.....field increment
      var1=-finc*0.5
      bf=714.4844244*freq/g(i)
c.....resonant field/Gauss, freq/GHz
      if(i.eq.2) goto 915
      print 110,rfd
      read(*,95,err=914) rfdD
      if(rfdD.eq.'n') rfdD='N'
      if(rfdD.eq.'y') rfdD='Y'
      if(rfdD.ne.'N'.and.rfdD.ne.'Y') goto 914
      rfd=rfdD
914   if(rfd.eq.'N') goto 915
      stfld=714.4844244*freq/g(1)
c.....move centre field to resonant
915   f=stfld-bf-(swp/2.0)
c.....field of radical 1
c
      do 920 j=1,npt
        do 24 lx=1,nli
          if(abs(f+sumami(lx,i))+var1.gt.0.0) goto 24
          sorb(j)=sorb(j)+v(lx,i)
24      continue
920   f=f+finc
      return
c
9992  format(' ',1X,'T.M. polarisation ratio ? (',F6.2,')',20X,
     &          '    ',\)
9993  format(F6.0)
90    format(' ',1X,'Do you want the polarised(P) or unpolarised(U)',
     &           1X,'spectrum ?(',A1,') ',\)
95    format(A1)
99    format(' ',1X,'Relative conc. of radical 1 ?(1.0) ',\)
100   format(G16.4)
110   format(' ',1X,'Want spectrum centred on resonant field of',
     &           1X,'radical 1 ?(',A1,')    ',\)
      end

