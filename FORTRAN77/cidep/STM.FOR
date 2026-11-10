$LARGE
      subroutine stm(itm)
c
      common /a/v(1024,2),total(1024,2)
      common /b/quant(20,2),a(20,2),nnuc(20,2),lines(20,2)
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /d/ngroup(2),nlines(2),nsec(2),itot(2),istm(2)
      common /e/qj(20,20,2),qm(20)
      common /f/ideg(20,20,2),nbuff(30),nk(10),jzcr(20),
     &          sumstm(1024,2),stmnet(2)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /t/rfd,rep
c
      dimension sorb(2048),rtot(2)
      dimension gD(2),ngroupD(2),nnucD(20,2),aD(20,2),quantD(20,2)
      equivalence (sorb(1),vrel(1,1))
c
      character same*1,sameD*1,sec*1,secD*1,pair*1,pairD*1
      character ans*1,rfd*1,rep*1,reply*1,fle*1,pha*1,esc*1,sb*1,csi*2
      character anet*1,anetD*1
c
      data same,pair,sec,anet /'Y','T','Y','Y'/
c
c.....STM calculates the form of the polarisation resulting from mixing
c.....between the S and T-minus states of a radical pair which are
c.....coupled by off-diagonal elements of the electron-nuclear hyperfine
c.....interaction. The ST-minus contribution to a given line in the
c.....spectrum of one radical consists of two terms, one due to nuclear
c.....spin flips from the particular nuclear spin state of the radical
c.....(multiplet contribution), the other due to nuclear spin flips on
c.....the counter radical (net contribution).
c
c.....stmpol(k,j,i) is the contribution to the multiplet effect generated
c.....by the kth component of the jth multiplet in radical i.
c
c.....sumstm(l,i) is the multiplet polarisation line l in the spectrum
c.....of radical i.
c
c.....stmnet(i) is the net effect on radical i which is the sum of the
c.....polarisations on the counter radical divided by the total number
c.....of nuclear spin functions of the counter radical.
c
c
      rep = 'P'
      rft = 'Y'
      esc=char(27)
      sb=char(91)
      csi=esc//sb
      istf=2
c
10    nflag=0
      call clear
      print 1000,csi,csi
      print 1025,same
      read(*,1120,err=9001) sameD
c.....same = 'N' if radicals not identical
      if(sameD.eq.'n') sameD='N'
      if(sameD.eq.'y') sameD='Y'
      if(sameD.ne.'Y'.and.sameD.ne.'N') goto 51
9001  same=sameD
51    i=1
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
60    continue
61    if(i.eq.2) goto 900
c.....if a > 0.3%
      if(same.ne.'N') goto 200
      i=2
      goto 50
200   ngroup(2)=ngroup(1)
      g(2)=g(1)
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
904   continue
c
      istm(1)=1
      istm(2)=1
c.....instructs subroutines to carry out
c.....st-minus operations
c
      call second(same)
c.....second order and st-minus components
905   ef=0
      call ami(same,ef)
c.....calc stick spectrum and st-minus polarisns
      if(ef.eq.1) return
c
      stmnet(1)=0.0
      stmnet(2)=0.0
      i=1
906   itot(i)=0
c.....itot(i) is found
      do 907 jg=1,nlines(i)
c.....here since pol is
c.....not called
        itot(i)=itot(i)+itotal(jg,i)
        stmnet(i)=stmnet(i)+sumstm(jg,i)
907   continue
      if(same.eq.'Y') itot(2)=itot(1)
      stmnet(i)=stmnet(i)/float(itot(i))
      if(same.eq.'Y') stmnet(2)=stmnet(1)
      if(same.eq.'N') i=i+1
      if(i.eq.2) goto 906
      rtot(1)=float(itot(1))
      rtot(2)=float(itot(2))
c
      if(same.eq.'Y') goto 908
c.....net effect on radical 1 derives from radical 2 and vice versa.
      sm=stmnet(2)
      stmnet(2)=stmnet(1)
      stmnet(1)=sm
c
908   i=1
      print 800,stmnet(1)
      print 801,stmnet(2)
800   format(1X,' Net effect on radical 1 (from radical 2) = -',G12.6)
801   format(1X,' Net effect on radical 2 (from radical 1) = -',G12.6
     &          ,/)
      do 909 iln=1,nlines(i)
        v(iln,1)=0.0
909     v(iln,2)=0.0
      ans='Y'
      print 1300,anet
c.....st-minus net effect ?
      read(*,1120,err=920) anetD
      if(anetD.eq.'n') anetD='N'
      if(anetD.eq.'y') anetD='Y'
      if(anetD.eq.'Y'.or.anetD.eq.'N') anet=anetD
920   if(anet.eq.'N') goto 930
c.....multiplet effect only
      do 925 iln=1,nlines(i)
925     v(iln,i)=v(iln,i)-stmnet(i)*itotal(iln,i)
930   do 935 iln=1,nlines(i)
        v(iln,i)=v(iln,i)-sumstm(iln,i)
935     v(iln,i)=v(iln,i)/rtot(i)
c.....same scale for both radicals
      do 940 iln=1,nlines(i)
940     vrel(iln,i)=v(iln,i)
      if(same.eq.'N') i=i+1
      if(i.eq.2) goto 920
c
      print 1230
c.....want to see the numerical results
      read(*,1120,err=990) ans
      if(ans.eq.'y') ans='Y'
      if(ans.ne.'Y') goto 990
c
c.....rearrange the table in order of increasing field.
c
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
        call rpmplt(i,itm,tm,2)
c.....plot the spectrum of the ith radical
c       call hold(1,2)
c.....wait for 1 second
        if(i.eq.2) goto 99
c
c.....for 2 radicals, 2nd plot will be sum
        call clear
        print *,'    '
        if(rep.eq.'U') goto 997
        print 1272
        read(*,1120,err=9007) reply
        if(reply.eq.'y') reply='Y'
9007    if(reply.ne.'Y') call sfile(1,istf)
        if(reply.eq.'Y') call sfile(0,istf)
995   continue
      if(reply.eq.'Y'.and.same.ne.'N') goto 997
      call clear
997   print 1280
      print 1250
      read(*,1120,err=9999) ans
      if(ans.ne.'m'.and.ans.ne.'M') goto 10
9999  return
c
1000  format(' ',A2,'32m',
     &  25X,'ST_ - Radical Pair Mechanism',/,
     &  1X,A2,'0m',/,30X,'( for triplet pairs )',//)
1021  format(/,1X,'Sign convention : A (+ve), Lowfield (+ve sumami)')
1025  format(' ',1X,'Are the two radicals identical?(',A1,') ',\)
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
     &  I2,'; hyp coupling(Gauss) ',F7.3,'; nuc spin ',F3.1,/,\)
1200  format(/,1X,'for radical(',I1,')')
1210  format(1X,'Line no.',T15,'pol',T31,'rel pol',T47,'size',
     &  T62,'deg',T69,'Sum of A*Mi')
1220  format(1X,T2,I4,T10,G14.6,T27,F12.6,T45,G13.6,T59,I5,T68,F10.4)
1230  format(' ',1X,'Table of polarisations?(N)',35X,\)
1231  format(' ',1X,'Output file ( printer(P) or terminal(T) )?(T)'
     &,15X,\)
1250  format(' ',1X,'Return to Master Program(M) or ST_(S) ?(S)',8X,\)
1272  format(' ',1X,'Do you want this ST_ stick spectrum on file ?'
     &  ,'(N)',7X,\)
1280  format(/)
1300  format(' ',1X,'Want to include ST_ net effect? (',A1,')',28X,\)
c
1120  format(A1)
1121  format(F8.6)
1122  format(I2)
1123  format(F7.0)
1124  format(F3.1)
      end

