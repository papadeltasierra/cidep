$LARGE
      subroutine stmfix
c
c.....Purpose : To add an ST_ contribution to a polarised stick spectrum.
c
c.....Input   : File (:KAM3)*:DATR from RPM
c.....        : Phase of TM
c.....        : Relative magnitude of TM wrt RPM
c
c.....TO CHANGE CMIX(2)
c
      common /mlist/nstick,nspec,stfld,freq1,swp,rat2
      common /tense/specp(2048),spece(2048),specpt(2048)
C
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
      character ans*1,rep2*1,esc*1,sb*1,csi*2
c
      esc=char(27)
      sb=char(91)
      csi=esc//sb
      if(mflag(1).eq.0.and.mflag(3).eq.0) then
        print 76
        call getscan( keyin, 0)
        return
      end if
10    call clear
      print 1,csi,csi
      print 2
      if(mflag(4).eq.1.or.mflag(4).eq.3) then
        ip=1
        ip1=5
      end if
c
      if(mflag(4).eq.5.or.mflag(4).eq.7) then
        ip=4
        ip1=6
      end if
      if(mflag(4).eq.4.or.mflag(4).eq.6) then
        ip=3
        ip1=13
      end if
76    format(1x,/,4x,'  no rpm or tm spectrum created - either create'
     &           ,' or read with tmfix ',/,
     &           5X,'Hit Any Key to Continue',/)
c
      if(mflag(2).eq.0) then
        print 77
        read(1,4,err=75) ans
        goto 73
75      ans='n'
73      if(ans.eq.'y'.or.ans.eq.'y') then
          call gendat(2,ier)
          if(ier.ne.0) return
        else
          print 74
          call getscan( keyin, 0 )
          return
        end if
        ip=7
        ip1=7
74      format(1x,/,5x,'create st_ - rpm polarised stick spectrum',/
     &              5x,'Hit Any Key to Continue',/)
77      format(' ',5x,'do you have an st_ spectrum on file ?(n) ',\)
      else
        nstick=cnpt
        swp=cswp
        stfld=cstfld
        freq1=cfreq1
        do 88 ii=1,cnpt
88        specp(ii)=cmix(1)*cstkp(1,ii)+cmix(3)*cstkp(3,ii)
        do 89 ii=1,cnpt
89        spece(ii)=cstkp(2,ii)
      end if
c
      print 3
      read(*,4,err=6) ans
      if(ans.eq.'y') ans='y'
      if(ans.ne.'y') goto 6
c
      call gmode
c     call unigo
      call scrnbox
c     call unibox
      call scrnstck( specp, nstick, 1 )
c     call unistk ( specp, nstick, fmin, ff )
      call getscan( keyin, 0 )
c     x = float ( nstick / 2 )
c     y = specp ( nstick / 2 )
c     call unicrs ( x, y, key )
c
      call clrscr
c     call gclear
      call scrnbox
c     call unibox
      call scrnstck ( spece, nstick, 1 )
c     call unistk ( spece, nstick, fmin, ff )
      call getscan( keyin, 0 )
c
c     x = float ( nstick / 2 )
c     y = spece ( nstick / 2 )
c     call unicrs ( x, y, key )
      call tmode
c     call uniend
c
c     if ( key.eq.112.or.key.eq.80 )
c    &   call tekgraph ( spece, nstick, 1, 1)
c
6     continue
      iflag=1
      call sfix(iflag,cmix,IP1)
      if(iflag.eq.1) goto 10
      call clear
c
      if(cmix(2).eq.0.0) then
        mflag(4)=iand( mflag(4), 5 )
      else
        mflag(4)=ior( mflag(4), 2 )
        do 78 ii=1,nstick
78        cstkt(ii)=specpt(ii)
      end if
c
      return
c
1     format(' ',A2,'32m',
     &           33X,'ST_ - RPM Fix',/,
     &           1X,A2,'0m',/,5X,
     &  '( adds ST_ contributions to polarised terms from STo',
     &  ' and/or TM )')
2     format(/)
3     format(' ',' Do you want to see stick spectra ?(N) ',8X,\)
4     format(A1)
c
      end
c
      subroutine sfix(iflag,cmix,ip)
c
      common /mlist/nstick,nspec,stfld,freq1,swp,rat2
      common /tense/specp(2048),spece(2048),specpt(2048)
      dimension cmix(3)
      character esc*1,ans*1
c
      esc=char(27)
      call clear
100   print 1
      print 2,cmix(2)
1     format(//)
2     format(' ',5X,'Input ST_ : STo (+/or TM) ratio (',F8.5,')   ',\)
      read(*,3,err=9001) ratio
3     format(F12.0)
      if(ratio.ne.0.0) cmix(2)=ratio
9001  print *,esc,'[3A'
      print 4,cmix(2)
4     format(7X,'ST_ : STo (+/or TM) ratio    ',F8.5,15X)
c
      do 5 i=1,nstick
        if(spece(i).eq.0.0.and.specp(i).eq.0.0) goto 5
        specpt(i)=(cmix(2)*spece(i))+specp(i)
5     continue
c
      call splots(specpt,nstick,ip)
      print 1
      print 6
6     format(' ',5X,'Do you want this (STo +/or ST_ + T.M. ) ',
     &              'spectrum on file ?(N) ',\)
      read(*,7,err=9002) ans
7     format(A1)
      if(ans.eq.'y') ans='Y'
      if(ans.eq.'Y') call tfile
9002  print 18
18    format(' ',5X,'Want to change ratio(C) or return to Master '
     &             ,'program ? (M) ? ',\)
      rat2=ratio
      read(*,7,err=9999) ans
      if(ans.eq.'c') ans='C'
      if(ans.ne.'C') iflag=-1
      if(ans.eq.'C') goto 100
9999  return
      end

