$LARGE
      subroutine mclear
c
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
c.....mflag 1 - sto in cstkp(1,)
c.....      2 - st_ in cstkp(2,)
c.....      3 - tm  in cstkp(3,)
c.....      4 - total - 0 - nothing in cstkt()
c.....                  1 - sto contribution to cstkt()
c.....                  2 - st_ contribution to cstkt()
c.....                  4 - tm  contribution to cstkt()
c
c.....         e.g.     7 - sto + st_ + tm contributions to cstkt()
c
      do 10 i=1,4
10      mflag(i)=0
c
      cmix(1)=0.0
      cmix(2)=0.0
      cmix(3)=0.0
c
      return
      end
c     ------------------------------------------------------------------
      subroutine rpm(itm)
c
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
C
      common /a/v(1024,2),total(1024,2)
      common /b/quant(20,2),a(20,2),nnuc(20,2),lines(20,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /f/ideg(20,20,2),nbuff(30),nk(10),jzcr(20)
     &         ,sumstm(1024,2),stmnet(2)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /t/rfd,rep
c
      dimension sorb(2048)
      dimension gD(2),ngroupD(2),nnucD(20,2),aD(20,2),quantD(20,2)
      equivalence (sorb(1),vrel(1,1))
c
      character same*1,sameD*1,sec*1,secD*1,pair*1,pairD*1
      character ans*1,rfd*1,rep*1,reply*1,fle*1,pha*1
      character esc*1,sb*1,csi*2
c
c.....Set default spectrum data.
c
      g(1)       = 2.0025
      ngroup(1)  = 2
      nnuc(1,1)  = 1
      a(1,1)     = 21.94
      quant(1,1) = 0.5
      nnuc(2,1)  = 6
      a(2,1)     = 24.69
      quant(2,1) = 0.5
c
      g(2)       = 2.0007
      ngroup(2)  = 1
      nnuc(1,2)  = 1
      a(1,2)     = 1.0
      quant(1,2) = 0.5
c
      swp        = 180.0
      stfld      = 3483.0
      freq       = 9.777
      npt        = 1024
c
      ca         = 0.0
      rnc        = 1.0
c
      same       = 'Y'
      pair       = 'T'
      sec        = 'Y'
      rep        = 'P'
      rfd        = 'Y'
      tm         = -50.0
c
      esc=char(27)
      sb=char(91)
      csi=esc//sb
      istf=1
10    nflag=0
c
      call clear
      if(itm.eq.0) print 1000,csi,csi
      if(itm.eq.1) print 9991,csi,csi
      if(itm.eq.1) then
        same='Y'
        istf=3
        goto 51
      end if
      print 1025,same
      read(*,1120,err=9001) sameD
c.....same = 'N' if radicals not identical
      if(sameD.eq.'n') sameD='N'
      if(sameD.eq.'y') sameD='Y'
      if(sameD.ne.'Y'.and.sameD.ne.'N') goto 51
9001  same=sameD
51    i=1
      istm(1)=0
      istm(2)=0
50    print 1026
      print 1030,i
      read(*,1121,err=9002) gD(i)
c.....g-value
      if(gD(i).ne.0.0) g(i)=gD(i)
9002  print *,esc,'[2A'
      print 1031,i,g(i)
      read(*,1122,err=9003) ngroupD(i)
c.....no. of groups of nuclei
      if(ngroupD(i).ne.0) ngroup(i)=ngroupD(i)
9003  print *,esc,'[2A'
      print 1032,i,g(i),ngroup(i)
      print *,'    '
      nsec(i)=0
      if(ngroup(i).eq.0) then
        ngroup(i)=1
        nnuc(1,i)=1
        a(1,i)=1.0
        quant(1,i)=0.0
        goto 61
      end if
      do 60 j=1,ngroup(i)
        if(j.gt.1) goto 99
        print 1040,j
        goto 991
99      print 1039,j
991     read(*,1122,err=9004) nnucD(j,i)
c.....no. of nuclei of group(j), radical(i)
        if(nnucD(j,i).ne.0) nnuc(j,i)=nnucD(j,i)
9004    print *,esc,'[2A'
        print 1041,j,nnuc(j,i)
        read(*,1123,err=9005) aD(j,i)
c.....hyp coupling constant
        if(aD(j,i).ne.0.0) a(j,i)=aD(j,i)
9005    print *,esc,'[2A'
        print 1042,j,nnuc(j,i),a(j,i)
        read(*,1124,err=9006) quantD(j,i)
c.....nuclear spin qu.no.
        if(quantD(j,i).ne.0.0) quant(j,i)=quantD(j,i)
9006    print *,esc,'[2A'
        print 1043,j,nnuc(j,i),a(j,i),quant(j,i)
60      if(abs(a(j,i)).gt.10.0) nsec(i)=1.0
c.....consider 2nd order
61    if(i.eq.2) goto 900
c.....if a > 0.3%
      if(same.ne.'N') goto 200
      i=2
      goto 50
200   ngroup(2)=ngroup(1)
      g(2)=g(1)
      nsec(2)=nsec(1)
      do 220 j=1,ngroup(1)
        nnuc(j,2)=nnuc(j,1)
        a(j,2)=a(j,1)
220     quant(j,2)=quant(j,1)
500   continue
900   continue
c
      call first(same,nflag)
c.....first order
      if(nflag.eq.1) goto 9999
      print 1026
      if(nsec(1).eq.0.and.nsec(2).eq.0) goto 905
      print 1060,sec
c.....want second order effects
      read(*,1120,err=9007) secD
      print 1026
      if(secD.eq.'y') secD='Y'
      if(secD.eq.'n') secD='N'
      if(secD.ne.'Y'.and.secD.ne.'N') goto 904
9007  sec=secD
904   if(sec.ne.'Y') goto 905
      nsec(1)=nsec(1)+1
      nsec(2)=nsec(2)+1
      call second(same)
c.....second order
905   ef=0
      call ami(same,ef)
c.....calc stick spectrum
      if(ef.eq.1) return
935   call pol(same,itm,irnc)
c.....calc polarisations
c
      if(itm.eq.1) goto 990
      print 1050,pair
      read(*,1120,err=9008) pairD
      if(pairD.eq.'t') pairD='T'
      if(pairD.eq.'s') pairD='S'
      if(pairD.eq.'f') pairD='F'
      if(pairD.ne.'T'.and.pairD.ne.'S'.and.pairD.ne.'F') goto 52
9008  pair=pairD
52    call option(pair,same)
c
      print 1230
c.....want to see the numerical results
      read(*,1120,err=990) ans
      if(ans.eq.'y') ans='Y'
      if(ans.ne.'Y') goto 990
c
c.....rearrange the table in order of increasing field.
c.....if(nlines(1).lt.100.and.nlines(2).lt.100) goto 945
c.....with >100 lines on either radical query the need for reordering.
c
c.....print 1070
c.....read(*,1120,945) ans
c.....if(ans.eq.'n') ans='N'
c.....if(ans.eq.'N') goto 948
c.....No
945   call sort(same)
c.....Yes
948   print 1231
      read(*,1120,err=949) fle
      if(fle.eq.'p') fle='P'
      if(fle.ne.'P') goto 949
      call facit(same,irnc)
      goto 990
949   print 1021
      i=1
950   print 1200,i
      print 1210
      nli=nlines(i)
      do 970 l=1,nli
970     print 1220,l,v(l,i),vrel(l,i),float(itotal(l,i)*v(l,i)),
     &             itotal(l,i),sumami(l,i)
      print *,'    '
975   if(same.ne.'N') goto 990
      if(i.eq.2) goto 990
      i=2
      goto 950
990   do 995 i=1,2
        if(i.eq.2.and.same.ne.'N') goto 995
        call field2(i,itm,tm)
c.....calc the spectrum for the ith radical
        call rpmplt(i,itm,tm,1)
c.....plot the spectrum of the ith radical
c.....call hold(1,2)
c.....wait for 1 second
        if(i.eq.2) goto 995
c.....for 2 radicals, 2nd plot will be sum
        call clear
        print *,'    '
        if(rep.eq.'U') goto 997
        if(itm.eq.0) print 1272
        if(itm.eq.1) print 1273
        read(*,1120,err=9009) reply
        if(reply.eq.'y') reply='Y'
9009    if(reply.ne.'Y') call sfile(1,istf)
        if(reply.eq.'Y') call sfile(0,istf)
c
        if(itm.eq.1) then
          ctotdeg=0.0
          do 899 itd=1,nlines(1)
899         ctotdeg=ctotdeg+float(itotal(itd,1))
        end if
c
995   continue
      if(reply.eq.'Y'.and.same.ne.'N') goto 997
      call clear
997   print 1280
      if(itm.eq.0) print 1250
      if(itm.eq.1) print 1251
      read(*,1120,err=9999) ans
      if(ans.ne.'M'.and.ans.ne.'m') goto 10
9999  return
c
9991  format(' ',//,' ',A2,'32m',
     &  31X,'Triplet Mechanism',/,
     &  1X,A2,'0m')
1000  format(' ',//,' ',A2,'32m',
     &  25X,'STo - Radical Pair Mechanism',/,
     &  1X,A2,'0m')
1021  format(/,1X,'Sign convention : A(+ve), Lowfield(+ve sumami)')
1025  format(' ',1X,'Are the two radicals identical ?(',A1,') ',\)
1026  format('     ')
1030  format(' ',1X,'Radical(',I1,') g-value ',\)
1031  format(' ',1X,'Radical(',I1,') g-value ',F8.6,T30,
     &  '; no. of groups of nuclei ',\)
1032  format(' ',1X,'Radical(',I1,') g-value ',F8.6,T30,
     &  '; no. of groups of nuclei ',I2,\)
1039  format(' ',3X,'group(',I2,') : no. of nuclei(<16) ',\)
1040  format(' ',3X,'group(',I2,') : no. of nuclei(<16) ',\)
1041  format(' ',3X,'group(',I2,') : no. of nuclei(<16) ',
     &  I2,'; hyp coupling(Gauss) ',\)
1042  format(' ',3X,'group(',I2,') : no. of nuclei(<16) ',
     &  I2,'; hyp coupling(Gauss) ',F7.3,'; nuc spin ',\)
1043  format(' ',3X,'group(',I2,') : no. of nuclei(<16) ',
     &  I2,'; hyp coupling(Gauss) ',F7.3,'; nuc spin ',F3.1,\)
1060  format(' ',1X,'Want to see the effect of second order',
     &  ' splittings ?(',A1,')       ',\)
1050  format(' ',1X,'S,T or F pair ?(',A1,')',43X,\)
1070  format(' ',1X,'Do you really want the lines reordered ?(Y)',
     &  17X,\)
1200  format(/,1X,'for radical(',I1,')')
1210  format(1X,'Line no.',T15,'pol',T31,'rel pol',T47,'size',
     &  T62,'deg',T69,'Sum of A*Mi')
1220  format(1X,T2,I4,T10,G14.6,T27,G12.6,T45,G13.6,T59,I5,T68,F10.3)
1230  format(' ',1X,'Table of polarisations?(N)',35X,\)
1231  format(' ',1X,'Output file ( printer(P) or terminal(T) ) ?(T)'
     &  ,15X,\)
1250  format(' ',1X,'Return to Master Program(M) or RPM(R) ?(R)',8X,\)
1251  format(' ',1X,'Return to Master Program(M) or TM(T) ?(T)',8X,\)
1272  format(' ',1X,'Do you want this STo stick spectrum on file ?'
     &  ,'(N) ',6X,\)
1273  format(' ',1X,'Do you want this TM stick spectrum on file ?'
     &  ,'(N) ',7X,\)
1280  format(/)
c
1120  format(A1)
1121  format(F8.6)
1122  format(I2)
1123  format(F7.0)
1124  format(F3.1)
      end

