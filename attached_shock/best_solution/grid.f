
      program ramp_grid
c *****************************************************************
c ****  Generates grid to fit a specified shock geometry       ****
c *****************************************************************

      Parameter(imax=500,jmax=500)

      Double precision xcor,dxc,xplt,xrmp,ang,dyw,ymax,ydis
      Double precision f,con,xp,xm,xc,yj,fc,gc
      Double precision x(jmax,imax),y(jmax,imax)
      Integer nx1,nx2,ny,i,j,ii,jj

      ang  = 20.0     ! ramp angle
      angt = 0.0     ! angle of top bdy
      ymax = 0.8    ! max y value at i=1

      nx1 = 100
      nx2 = 100
      nx  = nx1+nx2
      ny  = 50

      xcor = 0        ! x loc of corner
      xplt = 0.5      ! length of plate before the corner
      xrmp = 0.5     ! length of ramp after the corner

      ang  = ang *4.*atan2(1.,1.)/180.
      angt = angt*4.*atan2(1.,1.)/180.
c
c *** put points on j=2 line
c
      j = 2
      do i = 1,nx1+1
            f = float(i-1)/float(nx1)
         ii = nx1+3-i
         x(j,ii) = xcor - f*xplt
         y(j,ii) = 0.
      enddo
      x(j,1) = 2.*x(j,2) - x(j,3)
      y(j,1) = 0.

      do i = 1,nx2+1
            f = float(i-1)/float(nx2)
         ii = nx1+2 + i-1
         x(j,ii) = xcor + f*xrmp
         y(j,ii) = min(0.3,(x(j,ii) - xcor) * tan(ang))
      enddo
      x(j,nx1+nx2+3) = 2.*x(j,nx1+nx2+2) - x(j,nx1+nx2+1)
      y(j,nx1+nx2+3) = 0.3!(x(j,nx1+nx2+3) - xcor) * tan(ang)
c
c *** stretching in y-dir (i=1)
c
      do i = 1,nx+3
         ytop = ymax + max(0.,(x(2,i)-xcor)*tan(angt))
         ybot = y(2,i)
         ydis = ytop - ybot
         do j = 1,ny+1
            f = float(j-1)/float(ny)

            jj = j+1
            y(jj,i) = y(2,i) + f*ydis
         enddo
         y(1,i) = 2.*y(2,i) - y(3,i)
         y(ny+3,i) = 2.*y(ny+2,i) - y(ny+1,i)
      enddo
c
c *** form grid by shifting j=2 line up and right
c
      do i = 1,nx+3
         do j = 1,ny+3
            x(j,i) = x(2,i)
         enddo
      enddo
c
c *** Write tecplot file
c      
      open(12,file='grid.dat')
      write(12,*) 'variables = x,y'
      write(12,*) 'zone t=wedge i=',nx+3,', j=',ny+3,' f=block'
      write(12,102) ((x(j,i),i=1,nx+3),j=1,ny+3)
      write(12,102) ((y(j,i),i=1,nx+3),j=1,ny+3)

      xstr = 1.0
      do i = 2,nx+2
         dxrat = (x(1,i+1)-x(1,i))/(x(1,i)-x(1,i-1))
         if (xstr .lt. dxrat) xstr=dxrat
      enddo

      ystr = 1.0
      do j=2,ny+2
         dyrat = (y(j+1,1)-y(j,1))/(y(j,1)-y(j-1,1))
         if (ystr .lt. dyrat) ystr=dyrat
      enddo

      print *,'x-stretching = ',xstr
      print *,'y-stretching = ',ystr

101   format(i4,7E16.8)
102   format(7E16.8)

      stop
      end
