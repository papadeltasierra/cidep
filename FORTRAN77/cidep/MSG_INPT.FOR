$LARGE
      subroutine message(imess)
c
      character esc*1,sb*1,csi*2
c
      esc=char(27)
      sb=char(91)
      csi=esc//sb
c
      call vtc(32,1)
      print 1001,csi
      if(imess.lt.1.or.imess.gt.13) goto 1000
c
      goto (1,2,3,4,5,6,7,8,9,10,11,12) imess
c
1000  print 9999
      goto 99
1     print 101
      goto 99
2     print 102
      goto 99
3     print 103
      goto 99
4     print 104
      goto 99
5     print 105
      goto 99
6     print 106
      goto 99
7     print 107
      goto 99
8     print 108
      goto 99
9     print 109
      goto 99
10    print 110
      goto 99
11    print 111
      goto 99
12    print 112
      goto 99
13    print 113
      goto 99
c
99    print 1002,csi
      return
c
1001  format(A2,'7l',17X,\)
1002  format(16X,'A^',2X,A2,'0m',\)
c
101   format('       STo - RPM polarised  spectrum       ',\)
102   format('       ST_ - RPM polarised  spectrum       ',\)
103   format('         TM - polarised spectrum           ',\)
104   format('       STo + TM polarised spectrum         ',\)
105   format('       STo + ST_ polarised spectrum        ',\)
106   format('     STo + ST_ + TM polarised spectrum     ',\)
107   format('           Polarised spectrum              ',\)
108   format('          Unpolarised spectrum             ',\)
109   format('          Polarised terms only             ',\)
110   format('         Unpolarised terms only            ',\)
111   format('        Polarised and Eqm. terms           ',\)
112   format('Polarised and Eqm. terms with slow exchange',\)
113   format('        ST_ + TM polarised spectrum        ',\)
9999  format(80X)
      end
c     ------------------------------------------------------------------
      subroutine mont(imon)
c
      common /ch/mon
      character mon*3
      if(imon.eq.1)  mon='Jan'
      if(imon.eq.2)  mon='Feb'
      if(imon.eq.3)  mon='Mar'
      if(imon.eq.4)  mon='Apr'
      if(imon.eq.5)  mon='May'
      if(imon.eq.6)  mon='Jun'
      if(imon.eq.7)  mon='Jul'
      if(imon.eq.8)  mon='Aug'
      if(imon.eq.9)  mon='Sep'
      if(imon.eq.10) mon='Oct'
      if(imon.eq.11) mon='Nov'
      if(imon.eq.12) mon='Dec'
      return
      end
c     ------------------------------------------------------------------
      subroutine input(r)
c
      common /nlist/nstick,nspec,nhalf,nr(2),fchan,chan,nlines
c
      dimension r(5)
      character esc*1
c
      esc = char( 27 )
c
100   call clear
      print 1,esc,esc,(r(i),i=1,5),nr(2)
c
      do 20 j=1,5
c       call vtc(j+4,45)
c       print 6
        call vtc(j+4,47)
        read(*,2,err=3)fa
        if(fa.eq.0)go to 21
        r(j)=fa
c21      print 7
c       call vtc(j+4,45)
c       print 6
21      call vtc(j+4,47)
        print 8,r(j)
        goto 20
3       fa=0.0
        goto 21
20    continue
c
101   do 150,i=2,2
c       call vtc(10,45)
c       print 6
        call vtc(10,47)
        read(*,12,err=17)ifa
        goto 18
17      ifa=0
18      if(ifa.ne.0) nr(i)=ifa
        call vtc(10,47)
        print 11,nr(i)
150   continue
      print 89
      nr(1)=nstick
c
      return
c
1     format(' ',//,28x,a1,'[32mEnter Lineshape Parameters',a1,'[0m',//,
     &  10X,'T1                < ',G12.6,' > :  ',12x,' microseconds',/,
     &  10X,'T2                < ',G12.6,' > :  ',12x,' microseconds',/,
     &  10X,'w1                < ',G12.6,' > :  ',12x,' radian MHz  ',/,
     &  10X,'delay after flash < ',G12.6,' > :  ',12x,' microseconds',/,
     &  10X,'sample interval   < ',G12.6,' > :  ',12x,' microseconds',/,
     &  10X,'no. of samples            < ',I4,' > :  ',12x)
2     format(G12.0)
c6     format(' ',12X,\)
c7     format(' ',1X,\)
8     format(G12.6,\)
11    format(I4,\)
12    format(I4)
89    format(///,10X,'type ARROWS with SHIFT and CTRL to move cursor',/,
     &  10X,'type RETURN to exit cursor mode',/,
     &  10X,'type ESC in cursor mode to skip remaining plots',//)
      end
c
