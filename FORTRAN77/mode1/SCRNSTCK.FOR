c     -----------------------------------------------------------------
c     Subroutine to plot stick data to screen. The data will be
c     presented to the subroutine in real form which is not ideal but
c     can't be helped.
c     -----------------------------------------------------------------
c
      subroutine scrnstck( data, npts, nspec )
c
      parameter ( depth   = 332.0 )
      parameter ( iwidth  = 639 )
      parameter ( ioffset = 333 )
c
      dimension data(2048,*)
      logical isame
      real minimum, maximum
      integer*4 ix
      integer*2 icolour(6)
c
      npts_m = npts - 1
c
      if ( npts_m .eq. 0 ) then
        return
      endif
c
      call colour ( 2 )
c
      maximum = data(1,1)
      minimum = data(1,1)
      do 100 i = 1, nspec
        do 110 j = 1, npts
          maximum = amax0 ( maximum, data(j,i) )
          minimum = amin0 ( minimum, data(j,i) )
110     continue
100   continue
      minimum = amin0 ( minimum, 0 )
      if ( maximum .eq. minimum ) then
        scale = 1.0
      else
        scale = depth / ( maximum - minimum )
      endif
c
      ibase = ioffset + ( minimum * scale )
c
      call putpt ( 0, ibase )
      call dline ( iwidth, ibase )
c
      do 120 i = 1, nspec
c
        do 130 j = 1, npts
          if ( data(j, i) .ne. 0.0 ) then
            ix = ( int4(j - 1) * iwidth ) / npts_m
            iy = ioffset - nint( ( data(j,i) - minimum ) * scale )
            call putpt ( ix, ibase )
            call dline ( ix, iy )
          endif
130     continue
c
120   continue
c
      return
      end
