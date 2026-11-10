$LARGE
      subroutine sfile(istor,istf)
c
      common /s/sorb(2048),fill(352)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /fe/sorbu(2048)
C
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
      character fnam*40
      character esc*1
c
      esc=char(27)
      call feqm
      if(istor.eq.1) goto 98
      print *,'     '
151   print 1
      goto 15
10    print *,esc,'[2A'
      print 14
      print *,esc,'[2A'
      goto 151
15    read(*,2,err=10) fnam
c
      open( 3, form='binary', file=fnam, status='new', err=10 )
      write(3)npt
      write(3)swp,stfld,freq
      write(3)(sorb(i),i=1,npt)
      write(3)(sorbu(i),i=1,npt)
      endfile(3)
      close(3)
      print 3,fnam
C
98    if(cnpt.ne.npt.or.cswp.ne.swp.or.
     &   cstfld.ne.stfld.or.cfreq1.ne.freq1) then
        mflag(1)=0
        mflag(2)=0
        mflag(3)=0
        mflag(4)=0
      end if
c
      mflag(istf)=1
      cnpt=npt
      cswp=swp
      cstfld=stfld
      cfreq1=freq1
c.....scale data so that range is 100.0
c.....this makes the scaling in OPTION redundant
      smin=sorb(1)
      smax=sorb(1)
      do 88 jj=1,npt
        smin = amin0( smin, sorb(jj) )
88      smax = amax0( smax, sorb(jj) )
      scal = 100.0 / ( smax - smin )
c
      do 90 jj=1,npt
C
c.....  DO NOT SCALE DATA TO A RANGE OF 100.0 IF ISFLAG=1
c.....  THEREBY LEAVING DATA AS SCALED BY OPTION
C
        if(isflag.eq.1) then
          cstkp(istf,jj)=sorb(jj)
        else
          cstkp(istf,jj)=sorb(jj)*scal
        endif
90    cstke(jj)=sorbu(jj)
c
      if(mflag(4).eq.0) then
        do 91 jj=1,npt
91        cstkt(jj)=sorb(jj)
        mflag(4)=int(2**(istf-1))
        cmix(istf)=1.0
      end if
      if(istf.eq.1) cmix(1)=1.0
c
      return
1     format(' ',' Input filename..(n.b. radical 1 only)...... ',\)
14    format(1X,70X)
2     format(A40)
3     format(1X,' File : ',A27,' complete.',/,5X,
     &  '(n.b. POL and EQM stick spectra stored.) ',//)
c
        end
c.....------------------------------------------------------------------
      subroutine tfile
c
      common /mlist/nstick,nspec,stfld,freq1,swp,rat2
      common /tense/specp(2048),spece(2048),specpt(2048)
      character fnam*40
      character esc*1
c
      esc = char(27)
      print *,'        '
      goto 7
6     print *,esc,'[3A'
      print 31
      print *,esc,'[A'
7     print 1
1     format(' ',' Input filename : ',\)
      read(*,2,err=6) fnam
2     format(A40)
31    format('  Input filename....',
     &       '...                  ')
      open( 3, form='binary', file=fnam, status='new', err=6)
      print 20,fnam
      print *,'    '
      print *,'    '
      write(3)nstick
      write(3)swp,stfld,freq1
      write(3)(specpt(i),i=1,nstick)
      write(3)(spece(i),i=1,nstick)
      endfile(3)
      close(3)
20    format(5X,'writing to file : ',A27)
      return
      end
c.....------------------------------------------------------------------
      subroutine vfile(specp,nspec,nshape,swp,stfld,freq1,sr,nchan)
c
      dimension specp(4096),sr(5)
      character fnam*40
      character esc*1
c
      esc = char(27)
      print *,'     '
      goto 15
10    print *,esc,'[2A'
      print 14
      print *,esc,'[2A'
15    print 1
      read(*,2,err=10) fnam
c
      ist=(nshape/2)+1
      ifn=nspec-(nshape/2)
      npt=ifn-ist+1
      open( 3, form='binary', file=fnam, status='new',err=10)
      write(3)npt
      write(3)swp,stfld,freq1,sr(3),sr(4),sr(5),nchan
      write(3)(specp(i),i=ist,ifn)
      endfile(3)
      close(3)
78    format(I6)
      print 3,fnam
      return
1     format(' ',' Input filename :  ',\)
14    format('                                                   ')
2     format(A40)
3     format(1X,' File :  ',A27,' is complete.',//)
c
      end
c.....------------------------------------------------------------------
      subroutine facit(same,irnc)
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
      common /ch/mon
c
      character ref*25,same*1,mon*3
      integer cldat(7)
c
      call gettim ( ihr, imin, isec, i100th )
      call getdat ( iyr, imon, iday )
c
      print 10
      read(*,20,err=9001) ref
9001  open(5,file='PRN')
c
      write(5,1) ihr, imin, isec, iday ,mon, iday ,ref
      if(istm(1).eq.0) then
        write(5,2)
        if(irnc.eq.1) then
          write(5,7) rnc
        else
          if(ca.eq.0.0) write(5,3)
          if(ca.ne.0.0) write(5,8) ca
        end if
      end if
      if(istm(1).eq.1) write(5,6)
      do 22 i=1,2
        if(i.eq.2.and.same.ne.'N') goto 44
        write(5,4)i,g(i),ngroup(i)
        do 33 j=1,ngroup(i)
          write(5,5)j,nnuc(j,i),a(j,i),quant(j,i)
33      continue
44      continue
22    continue
c
c
      write(5,1019)
      write(5,1021)
      i=1
950   write(5,1200)i
      write(5,1210)
      nli=nlines(i)
      do 970 l=1,nli
970     write(5,1220) l,v(l,i),vrel(l,i),float(itotal(l,i)*v(l,i)),
     &                itotal(l,i),sumami(l,i)
      write(5,1000)
975   if(same.ne.'N') goto 990
      if(i.eq.2) goto 990
      i=2
      goto 950
990   write(5,1230)
      close(5)
      return
c
10    format(' ',1X,'Reference..(<=25characters)     ',\)
20    format(A25)
1     format(5X,I2.2,':',I2.2,':',I2.2,3X,I2.2,'-',A3,'-',I4,28X
     &  ,'ref: ',A25)
2     format(/,22X,'STo - Radical Pair Mechanism simulation',/
     &  ,22X,'=======================================')
3     format(/,5X,'n.b. "Q**1/2"-dependance used',/)
4     format(/,5X,'radical(',I1,') : g-value ',F8.6,' ; no. of groups'
     &  ,1X,I2)
5     format(10X,'group(',I2,') : no. of nuclei ',I2,'; hyp coupling'
     &  ,'(Gauss) ',F6.3,'; nuc spin ',F3.1)
6     format(/,22X,'ST_ - Radical Pair Mechanism simulation',/
     &  ,22X,'=======================================',/)
7     format(/,5X,'n.b. "Q**',F4.2,'"-dependance used',/)
8     format(/,5X,'n.b. "Q**1/2 - ',F6.4,'*Q"-dependance used',/)
1000  format(/)
1019  format(//,10X,'Table of polarisations ',/,10X,
     &  '----------------------')
1021  format(/,5X,'Sign convention : A(+ve), Lowfield(+ve sumami)')
1200  format(/,5X,'for radical(',I1,')')
1210  format(1X,'Line no.',T15,'pol',T32,'rel pol',T48,'size',
     &  T62,'deg',T70,'Sum of A*Mi')
1220  format(5X,T2,I4,T10,G14.6,T27,F12.6,T45,G13.6,T59,I5,T69,F10.4)
1230  format(//,25X,'---------------------------------')
c
        end
c     -----------------------------------------------------------------
        subroutine gendat(IST,IER)
c
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
      common /mlist/nstick,nspec,stfld,freq1,swp,rat2
      common /tense/specp(2048),spece(2048),specpt(2048)
      character fnam*40
      character esc*1
c
      esc=char(27)
      ier=0
      print 4
4     format(/)
      goto 7
6     print *,esc,'[3A'
      ier = ier + 1
      if(ier.gt.4) return
      print 31
      print *,esc,'[A'
7     print 1
1     format(' ',' Input filename..(stick spectrum) : ',\)
      read(*,2,err=6) fnam
2     format(A40)
31    format('  Input filename..(stick spectrum)..',
     &  '...                  ')
      open( 3, form='binary', file=fnam, status='old',err=6)
      print 20,fnam
      print *,'    '
      print *,'    '
      read(3)nstick
      read(3)swp,stfld,freq1
      read(3)(specp(i),i=1,nstick)
      read(3)(spece(i),i=1,nstick)
      close(3)
20    format(5X,'Reading file : ',A27)
      ier = 0
C
      if(ist.eq.1) then
        do 77 i=1,nstick
          cstkp(1,i)=specp(i)
77        cstke(i)=spece(i)
        cnpt=nstick
        cswp=swp
        cstfld=stfld
        cfreq1=freq1
        cmix(1)=1.0
      else
        do 78 i=1,nstick
          cstkp(2,i)=specp(i)
78        cstke(i)=spece(i)
        if(cnpt.ne.npt.or.cswp.ne.swp.or.
     &     cstfld.ne.stfld.or.cfreq1.ne.freq1) then
          print 72
          call getscan(keyin)
        end if
        cnpt=nstick
        cswp=swp
        cstfld=stfld
        cfreq1=freq1
        cmix(2)=1.0
      end if
c
      return
c
72    format(/,5X,'WARNING - PARAMETERS DO NOT MATCH',/
     &       /,5X,'Press Any Key To Contine',/)
c
      end
c
c.....------------------------------------------------------------------
      subroutine gendat1(ier)
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /sense/specp(4096),spece(4096),deg(1000),pol(1000),
     &              ipos(1000),shapep(1024),shapee(1024),swp,stfld,freq1
      character fnam*40
      character esc*1
c
      esc=char(27)
      ier = 0
      print *,'   '
      goto 7
6     print *,esc,'[3A'
      ier = ier + 1
      if ( ier .gt. 4 ) return
      print 31
      print *,esc,'[A'
7     print 1
1     format(' ',' Input filename..(stick spectrum) : ',\)
      read(*,2,err=6) fnam
2     format(A40)
31    format('  Input filename..(stick spectrum)..',
     &  '...                  ')
      open( 3, form='binary', file=fnam, status='old', err=6 )
      print 20,fnam
      print *,'    '
      print *,'    '
      read(3)nstick
      read(3)swp,stfld,freq1
      read(3)(specp(i),i=1,nstick)
      read(3)(spece(i),i=1,nstick)
      close(3)
20    format(5X,'Reading file : ',A27)
      IER=0
      return
      end

