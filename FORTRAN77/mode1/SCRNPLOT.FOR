c     -----------------------------------------------------------------
c     Subroutine to plot data to screen. The data will be presented to
c     the subroutine in real form which is not ideal but can't be
c     helped.
c     -----------------------------------------------------------------
c
      subroutine scrnplot( data, npts, nspec, isame )
c
      parameter ( depth   = 331.0 )
      parameter ( iwidth  = 639 )
      parameter ( ioffset = 333 )
c
      dimension data(2048,*)
      logical isame
      real minimum, maximum
      integer*4 ix
      integer*2 icolour(6)
c
      data icolour / 13, 2, 12, 9, 14, 11 /
c
      npts_m = npts - 1
c
      if ( npts_m .eq. 0 ) then
        return
      endif
c
      if ( isame ) then
        maximum = data(1,1)
        minimum = data(1,1)
        do 100 i = 1, nspec
          do 110 j = 1, npts
            maximum = amax0 ( maximum, data(j,i) )
            minimum = amin0 ( minimum, data(j,i) )
110       continue
100     continue
        if ( maximum .eq. minimum ) then
          scale = 1.0
        else
          scale = depth / ( maximum - minimum )
        endif
      endif
c
      do 120 i = 1, nspec
c
        if ( .not.isame ) then
          maximum = data(1,i)
          minimum = data(1,i)
          do 130 j = 1, npts
            maximum = amax0 ( maximum, data(j,i) )
            minimum = amin0 ( minimum, data(j,i) )
130       continue
          if ( maximum .eq. minimum ) then
            scale = 1.0
          else
            scale = depth / ( maximum - minimum )
          endif
        endif
c
        jcolour = 1 + mod(i, 6)
        kcolour = icolour( jcolour )
        call colour ( kcolour )
c
        ix = 0
        iy = ioffset - nint( ( data(1,i) - minimum ) * scale )
c
        call putpt ( ix, iy )
c
        do 140 j = 2, npts
          ix = ( int4(j - 1) * iwidth ) / npts_m
          iy = ioffset - nint( ( data(j,i) - minimum ) * scale )
          call dline ( ix, iy )
140     continue
c
120   continue
c
      return
      end
