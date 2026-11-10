$LARGE
      subroutine rpmplt(i,itm,tm,istf)
c
      common /c/sumami(1024,2),itotal(1024,2),vrel(1024,2)
      common /e/qj(20,20,2),qm(20)
      common /factor/g(2),gg(2),fill2(2)
      common /fac/finc,freq,stfld,fscan,npt,swp,ca,rnc
      common /g/buff(18)
      common /s/sorb(2048),fill(352)
c
      fmin=sorb(1)
      fmax=sorb(1)
      do 28 j=1,npt
        if(fmin.gt.sorb(j)) fmin=sorb(j)
28      if(fmax.lt.sorb(j)) fmax=sorb(j)
      call splotr(fmax,fmin,npt,i,itm,tm,istf)
      return
      end
c     ------------------------------------------------------------------
      subroutine splotr(fmax,fmin,npt,irad,itm,tm,istf)
c
      common /s/sorb(2048),fill(352)
      dimension psorb(2048)
      character ans*1
c
      ff=fmax-fmin
      if(ff.eq.0.0) ff=1.0
c
      call gmode
c     call unigo
      call scrnbox
c     call unibox
      call scrnstck ( sorb, npt, 1 )
c     call unistk ( sorb, npt, fmin, ff )
c     call gcharj ( 1 )
c     call gcharf ( 'hard' )
c     call gchar  ( 'Radical 1 only$', ( xsizmm/2.0), 2.0, 4.0 )
c
      call getscan( keyin, 0 )
c     x = npt / 2
c     y = sorb ( npt / 2 )
c     call unicrs ( x, y, key )
c
      call tmode
c     call uniend
c
c     if(itm.eq.0.and.irad.eq.1) print 81
c     if(itm.eq.0.and.irad.eq.2) print 82
c     if(itm.eq.1) then
c              call message(3)
c              call vtc(4,53)
c              print 84,tm
c     else
c              if(istf.eq.1) call message(1)
c              if(istf.eq.2) call message(2)
c     end if
c81   format('$','Radical 1 only ')
c82   format('$','Radicals 1 & 2 ')
c84   format('$','P(I)/Peq = ',F6.2)

      end
c     ------------------------------------------------------------------
      subroutine splots (specl,nstick,ia)
c
      character ans*1
      dimension specl(2048),splt(2048)
c
      call gmode
c     call unigo
      call scrnbox
c     call unibox
      call scrnstck ( specl, nstick, 1 )
c     call unistk ( specl, nstick, fmin, ff )
      call getscan( keyin, 0 )
c
c     x = float ( nstick / 2 )
c     y = specl ( nstick / 2 )
c     call unicrs ( x, y, key )
c
      call tmode
c     call uniend
c
9999  return
      end
c     ------------------------------------------------------------------
      subroutine splot(specp,nspec,ib,nshape,iabf)
c
      dimension specp(4096),spt(2048)
      character ans*1
c
c     call message(ib)
      istart=1+(nshape/2)
      iend=(nspec)-(nshape/2)
      inn=(nspec-nshape)
c
      do 89 ik=1,inn
89      spt(ik)=specp(istart-1+ik)
c
      call gmode
c     call unigo
      call scrnbox
c     call unibox
      call scrnplot( spt, inn, 1, .FALSE. )
c     call uniplt ( spt, inn, 1, 1)
      call scrncurs( spt, inn, keyin )
c     call unicurs( spt, inn )
      call tmode
c
c      CALL VT220
c38    CALL VTC (3,2)
c      CALL CURSOR ( NX , NY , ICH )
c      IF (ICH.EQ.3) GOTO 37
c      IF (ICH.EQ.13.OR.ICH.EQ.10) GOTO 399
c      IDX = ( NX * NNPTS ) / 4095
c      DY = SPECP(IDX+(NSHAPE/2))
c      IF (DY.LT.999.0.AND.DY.GT.0.0001) THEN
c               PRINT 39,IDX,DY
c      ELSE
c                 PRINT 35,IDX,DY
c        ENDIF
c        GOTO 38
cC
c399     CONTINUE
cC
c        CALL VTC(3,2)
c        print 4
c        read(*,5,err=36) ans
c        if(ans.eq.'y') ans='Y'
c        IF(ANS.EQ.'Y') THEN
c                 call hplot(specp,nspec,ib,nshape)
cC                 CALL REGBOX
c        END IF
c4       format('$','    Hard plot ?(N) ')
c        GOTO 36
c37      IABF=1
cC        CALL REGIS(0)
c        RETURN
c5       format(A1)
c36      if(ib.ge.11.and.ans.ne.'Y') THEN
cC                 call REGIS(0)
c        ENDIF
c
      return
      end
