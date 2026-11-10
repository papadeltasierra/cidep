c     -----------------------------------------------------------------
c     Subroutine to draw a box into which a spectrum may be plotted.
c     -----------------------------------------------------------------
      subroutine scrnbox
c
      parameter ( idepth  = 335 )
      parameter ( iwidth  = 639 )
c
      call colour ( 15 )
      call putpt  ( 0, 0 )
      call dline  ( 0, idepth )
      call dline  ( iwidth, idepth )
      call dline  ( iwidth, 0 )
      call dline  ( 0, 0 )
c
      return
      end
