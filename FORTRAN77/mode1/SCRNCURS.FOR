c     -----------------------------------------------------------------
c     Subroutine to display a cursor, move the cursor, and display the
c     position of the cursor within the data set and the value of the
c     data point.
c     -----------------------------------------------------------------
c
      subroutine scrncurs( data, npts, keyin )
c
      parameter ( idepth  = 334 )
      parameter ( iwidth  = 639 )
c
c     The print line contains values which are hardcoded for the
c     particular screen format. Change these as well as the parameters
c     when porting to a different screen.
c
      dimension data(2048)
      character esc*1
      real minimum, maximum
      integer*4 ix, ixnew
c
      esc = char( 27 )
      npts_m = npts - 1
c
      if ( npts_m .eq. 0 ) then
        return
      endif
c
      call colxor ( 1 )
      call colour ( 15 )
c
      idatidx = npts / 2
      ix = ( int4(idatidx - 1) * iwidth ) / npts_m
      call putpt ( ix, 1 )
      call dline ( ix, idepth )
c
110   call getscan( keyin, 0 )
c
      if ( keyin .eq. #4b ) then
        idatidx = idatidx - 1
      elseif ( keyin .eq. #4d ) then
        idatidx = idatidx + 1
      elseif ( keyin .eq. #73 ) then
        idatidx = idatidx - 16
      elseif ( keyin .eq. #74 ) then
        idatidx = idatidx + 16
      elseif ( keyin .eq. #9b ) then
        idatidx = idatidx - 256
      elseif ( keyin .eq. #9d ) then
        idatidx = idatidx + 256
      elseif ( keyin .eq. #e0 .or. keyin .eq. #1c ) then
        goto 120
      else
        goto 110
      endif
c
130   if ( idatidx .gt. npts ) then
        idatidx = idatidx - npts
        goto 130
      else if ( idatidx .le. 0 ) then
        idatidx = idatidx + npts
        goto 130
      endif
c
      ixnew = ( int4(idatidx - 1) * iwidth ) / npts_m
c
      if ( ixnew .ne. ix ) then
        call dline ( ix, 1 )
        ix = ixnew
        call putpt ( ix, 1 )
        call dline ( ix, idepth )
      endif
c
      print 140, esc, idatidx, esc, data(idatidx)
140   format (a1,'[25;60H',i4,a1,'[25;70',f8.2,\)
      goto 110
c
120   continue
      call dline ( ix, ioffset )
      call colxor ( 0 )
c
      return
      end
