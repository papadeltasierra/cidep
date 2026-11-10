$LARGE
      program cidep
c
c     esr master simulation program              august 1986
c
      common / master / mflag(4),cnpt,cswp,cstfld,cfreq1,
     +                  cstkp(3,2048),cstke(2048),cstkt(2048),
     &                  cmix(3),ctotdeg,isflag
c
      character esc*1,sb*1,csi*2
      esc=char(27)
      sb=char(91)
      csi=esc//sb
c
100   call clear
      print 110,csi,csi,csi,csi
C
110   format (' ',//,24x,A2,'32m',
     &          'CIDEP master simulation program',
     &          a2,'0m',//,
     &12x,'1  -  create STo - RPM polarised stick spectrum',/,
     &12x,'2  -  add a TM distortion to a polarised stick spectrum',/,
     &12x,'3  -  create pure TM polarised stick spectrum',/,
     &12x,'4  -  convolute lineshape to stick spectrum',/,
     &12x,'5  -  create ST_ - RPM polarised stick spectrum',/,
     &12x,'6  -  add ST_ contribution to polarised stick spectrum' ,/,
     &12x,'7  -  convolute slow exchange lineshape to TM stick spectrum'
     &,/,
     &12x,'8  -  clear current data',/,
     &12x,'9  -  set flag to prevent stick spectrum normalization',/,
     &12x,'0  -  exit',//,36x,a2,'32m','Option ? ',a2,'0m'\)
c
120   call getscan ( keyin, 0 )
      keyin = keyin - #01
c
      if (keyin .lt. 1 .or. keyin .gt. 10 ) goto 120
      goto (20,30,40,50,60,70,80,90,130,10) keyin
c
c.....0 pressed
10    call clear
      stop
c
c.....1 pressed
20    call mclear
      call rpm(0)
      goto 100
c
c.....2 pressed
30    call tmfix
      goto 100
c
c.....3 pressed
40    call mclear
      call rpm(1)
      goto 100
c
c.....4 pressed
50    call specl
      goto 100
c
c.....5 pressed
60    call stm(0)
      goto 100
c
c.....6 pressed
70    call stmfix
      goto 100
c
c.....7 pressed
80    call spece
      goto 100
c
c.....8 pressed
90    call mclear
      isflag=0
      goto 100
c
c.....9 pressed
130   isflag=1
      goto 100
c
      end
