$LARGE
      subroutine specl
c
c.....Purpose : To add a full lineshape ( both pol. and equilibrium terms )
c.....to a polarised  stick spectrum.
c
c.....Input   : File containing stick spectrum.
c.....        : Lineshape parameters.
c
      common /nlist/nstick,nspec,nhalf,nshape,nchan,fchan,chan,
     &              nlines
      common /sense/specp(4096),spece(4096),deg(1000),pol(1000),
     &              ipos(1000),shapep(1024),shapee(1024),swp,stfld,freq1
c
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     &                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
      dimension r(5),sr(5)
      complex chan
      character ans*1,rep2*1,rp3*1,esc*1,sb*1,csi*2
c
      data r / 5.0, 0.6, 0.1, 1.0, 0.01 /
c
      nchan = 150
      esc=char(27)
      sb=char(91)
      csi=esc//sb
10    call clear
      print 1,csi,csi
      print 2
      if(mflag(4).eq.0) then
        call gendat1(IER)
        if(ier.ne.0) return
      else
        nstick=cnpt
        swp=cswp
        stfld=cstfld
        freq1=cfreq1
        do 36 kk=1,nstick
36        specp(kk)=cstkt(kk)
39      do 35 kk=1,nstick
35        spece(kk)=cstke(kk)
      end if
      print 3
      read(*,4,err=6) ans
      if(ans.eq.'y') ans='Y'
      if(ans.ne.'Y') goto 6
c
      call gmode
      call scrnbox
      call scrnstck ( specp, nstick, 1 )
      call getscan( keyin )
c
      call clrscr
      call scrnbox
      call scrnstck ( spece, nstick, 1 )
      call getscan ( keyin )
      call tmode
c
c     if ( keyin .eq. #?? ) call tekgraph ( spece, nstick, 1, 1)
c
6     continue
      call nozero
40    call clear
      call input(r)
      do 93 l=1,5
93      sr(l)=r(l)
c
      call const(r)
      call loop
      call convol
c
      inpts = nspec - nshape - 1
      ifred = nshape/2 + 1
c
      call gmode
      call scrnbox
      call scrnplot( specp(ifred), inpts, 1, .FALSE. )
      call scrncurs( specp(ifred), inpts, keyin )
      if ( keyin .eq. #01 ) goto 7099
c
      call clrscr
      call scrnbox
      call scrnplot( spece(ifred), inpts, 1, .FALSE. )
      call scrncurs( spece(ifred), inpts, keyin )
      if ( keyin .eq. #01 ) goto 7099
c
c48    do 50 i=1, nspec
      iend = (inpts + ifred - 1)
48    do 50 i=ifred, iend
50      specp(i)=specp(i)+spece(i)
c
      call clrscr
      call scrnbox
      call scrnplot( specp(ifred), inpts, 1, .FALSE. )
      call scrncurs( specp(ifred), inpts, keyin )
      call tmode
c
c.....iabf=0
c.....call splot(specp,nspec,9,nshape,iabf)
c.....if(iabf.eq.1) goto 70
c.....call splot(spece,nspec,10,nshape,iabf)
c.....if(iabf.eq.1) goto 70
c
c.....call splot(specp,nspec,11,nshape,iabf)
c
c
      print 199
      print 200,pol(1),deg(1)
      print 198
      print 201
      read(*,101,err=70) polmag
      if(polmag.eq.0.0) goto 70
      print *,esc,'[3A'
      print 202,polmag
      print 198
      polmag = polmag / pol(1)
      do 60 i=1,nlines
60      pol(i)=pol(i)*polmag
      call convol
      call gmode
      goto 48
7099  call tmode
c     call uniend
70    finc=swp/(float(nstick)-1.0)
      fswp=finc*(float(nspec)-1.0)
      if(iabf.eq.1) then
        call clear
        print 198
        goto 66
      end if
      print 198
c.....print 300,finc,fswp
      print 301,nstick
      print 198
      print 350
      read(*,450,err=66) rp3
      if(rp3.eq.'y') rp3='Y'
      if(rp3.eq.'Y') call vfile(specp,nspec,nshape,swp,stfld,freq1,
     &                          sr,nchan)
66    print 198
      print 400
      read(*,450,err=11) rep2
      if(rep2.ne.'M'.and.rep2.ne.'m') goto 40
      goto 11
c
1     format(' ',//,' ',A2,'32m',
     &  20X,'E.S.R. simulation by direct convolution',/,
     &  1X,A2,'0m',/,8X,
     &  '(n.b. incorporates both polarised and equilibrium terms.)')
2     format(/)
3     format(' ',' Do you want to see stick spectra ?(N) ',8X,\)
4     format(A1)
c
199   format(///)
198   format(' ',/)
200   format(6X,'Polarisation of 1st line : ',G12.6,/,
     &  6X,'degeneracy of 1st line   : ',G12.6,/)
201   format(' ',5X,'New Polarization (0 to exit) ? <0.0> ',\)
202   format(6X,'New Polarization           ',G12.6,/)
300   format(6X,'Field increment is        : ',G12.6,' (same as '
     &  ,'stick spectrum.)',/,6X,'Field sweep, however, now : ',G12.6
     &  ,' Gauss',/)
301   format(6X,'n.b. no. of points in the final spectrum : ',I4,/
     &  ,6X,'( i.e. same as in stick spectrum. )')
350   format(' ',5X,'Do you want spectrum on file  ?(N) ',\)
400   format(' ',5X,'Fresh lineshape(L) or return to Master program'
     &  ,'(M) ? (L) ',\)
101   format(G12.6)
450   format(A1)
c
11    call clear
      return
      end

