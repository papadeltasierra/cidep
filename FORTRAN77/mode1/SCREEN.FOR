      subroutine clear
c
      character esc*1
c
      esc = char(27)
      print 10, esc, esc
10    format (' ',a1,'[2J',a1,'[0;0H',\)
c
      return
      end
c     -----------------------------------------------------------------
      subroutine vtc (irow, icol)
c
      character esc*1
c
      esc = char(27)
      print 10, esc, irow, icol
10    format (' ',a1,'[',i2.2,';',i2.2,'H',\)
c
      return
      end
