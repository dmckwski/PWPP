      module fft
      implicit none

      contains

      subroutine fft1d(nblk,ndipole,datin,datout,isign)
      implicit none
      integer :: nblk,ndipole,isign,oldndipole,ix
      save :: oldndipole
      integer, parameter :: mxtrig=10000
      real(8), save :: trig(mxtrig)
      real(8) :: rdat(nblk,0:ndipole-1),idat(nblk,0:ndipole-1),pi
      complex(8) :: datin(nblk,0:ndipole-1)
      complex(8) :: datout(nblk,0:ndipole-1)
      data oldndipole/0/
      data pi/3.141592653589793/

      if(ndipole.ne.oldndipole) then
        oldndipole=ndipole
!write(*,*) 'fft',ndipole,nblk,isign
        call setgpfa(trig(:),ndipole)
      endif
      do ix=0,ndipole-1
         rdat(:,ix)=dble(datin(:,ix))
         idat(:,ix)=dimag(datin(:,ix))
      enddo
      call cgpfa(rdat,idat,trig,nblk,ndipole,isign)
      do ix=0,ndipole-1
         datout(:,ix)=dcmplx(rdat(:,ix),idat(:,ix))
      enddo
      datout=datout/sqrt(dble(ndipole))
      end subroutine fft1d

      subroutine fft2d(nblk,ndx,ndy,wx,wy,s0x,s0y,datin,datout,isign)
      implicit none
      integer :: nblk,ndx,ndy,isign,oldndx,oldndy,ix,iy, &
                 i1x,i2x,i1y,i2y
      save :: oldndx,oldndy
      integer, parameter :: mxtrig=10000
      real(8), save :: trigx(mxtrig),trigy(mxtrig)
      real(8) :: rdat(nblk,-max(ndx,ndy)/2:(max(ndx,ndy)-1)/2), &
                 idat(nblk,-max(ndx,ndy)/2:(max(ndx,ndy)-1)/2), &
                 s0x,s0y,dx,dy,pi,wx,wy
      complex(8) :: datin(nblk,-ndx/2:(ndx-1)/2,-ndy/2:(ndy-1)/2), &
                    datout(nblk,-ndx/2:(ndx-1)/2,-ndy/2:(ndy-1)/2),ci
      complex(8) :: pfx(-ndx/2:(ndx-1)/2), &
                    pfy(-ndy/2:(ndy-1)/2), &
                    pf0x(-ndx/2:(ndx-1)/2), &
                    pf0y(-ndy/2:(ndy-1)/2), &
                    pfin(-max(ndx,ndy)/2:(max(ndx,ndy)-1)/2), &
                    pfout(-max(ndx,ndy)/2:(max(ndx,ndy)-1)/2)
      data oldndx,oldndy,ci/0,0,(0.d0,1.d0)/
      data pi/3.141592653589793/
      i1x=-ndx/2
      i2x=(ndx-1)/2
      i1y=-ndy/2
      i2y=(ndy-1)/2
      if(ndx.ne.oldndx) then
        call setgpfa(trigx(:),ndx)
      endif
      if(ndy.ne.oldndy) then
        call setgpfa(trigy(:),ndy)
      endif
      oldndx=ndx
      oldndy=ndy
      dx=wx/dble(ndx)
      dy=wy/dble(ndy)
      do ix=i1x,i2x
         pfx(ix)=cdexp(isign*ci*(s0x*dx*(dble(ix))+2.*pi*dble((ix-i1x)*i1x)/dble(ndx)))
         pf0x(ix)=cdexp(isign*ci*(2.*pi*dble((ix-i1x)*i1x)/dble(ndx)))
      enddo
      do iy=i1y,i2y
         pfy(iy)=cdexp(isign*ci*(s0y*dy*(dble(iy))+2.*pi*dble((iy-i1y)*i1y)/dble(ndy)))
         pf0y(iy)=cdexp(isign*ci*(2.*pi*dble((iy-i1y)*i1y)/dble(ndy)))
      enddo

      if(isign.eq.-1) then
         pfin(i1y:i2y)=pfy(i1y:i2y)
         pfout(i1y:i2y)=pf0y(i1y:i2y)
      else
         pfin(i1y:i2y)=pf0y(i1y:i2y)
         pfout(i1y:i2y)=pfy(i1y:i2y)
      endif
      do ix=i1x,i2x
         rdat=0.d0
         idat=0.d0
         do iy=i1y,i2y
            rdat(:,iy)=dble(datin(:,ix,iy)*pfin(iy))
            idat(:,iy)=dimag(datin(:,ix,iy)*pfin(iy))
         enddo
         call cgpfa(rdat(1:nblk,i1y:i2y),idat(1:nblk,i1y:i2y),trigy,nblk,ndy,isign)
         do iy=i1y,i2y
            datout(:,ix,iy)=dcmplx(rdat(:,iy),idat(:,iy))*pfout(iy)
         enddo
      enddo
      if(isign.eq.-1) then
         pfin(i1x:i2x)=pfx(i1x:i2x)
         pfout(i1x:i2x)=pf0x(i1x:i2x)
      else
         pfin(i1x:i2x)=pf0x(i1x:i2x)
         pfout(i1x:i2x)=pfx(i1x:i2x)
      endif
      do iy=i1y,i2y
         rdat=0.d0
         idat=0.d0
         do ix=i1x,i2x
            rdat(:,ix)=dble(datout(:,ix,iy)*pfin(ix))
            idat(:,ix)=dimag(datout(:,ix,iy)*pfin(ix))
         enddo
         call cgpfa(rdat(1:nblk,i1x:i2x),idat(1:nblk,i1x:i2x),trigx,nblk,ndx,isign)
         do ix=i1x,i2x
            datout(:,ix,iy)=dcmplx(rdat(:,ix),idat(:,ix))*pfout(ix)
         enddo
      enddo
      datout=datout*cdexp(isign*ci*2.d0*pi*(dble(i1x*i1x)/dble(ndx)+dble(i1y*i1y)/dble(ndy)))/sqrt(dble(ndx*ndy))
      end subroutine fft2d

      subroutine cgpfa(cr,ci,trig,nblk,m,isign)
!      use iso_c_binding
      implicit none
      integer :: m,isign,nblk,n,i,inc
      real(8) :: trig(*),cr(nblk*m),ci(nblk*m)
!      real(8), pointer :: cr(:)
!      real(8) :: cr(2*nblk*m)
!      complex(8), target :: c(nblk*m)
!      call C_F_POINTER(C_LOC(c), cr, [2*nblk*m])
!      cr(1:2*nblk*m-1:2)=dble(c(1:nblk*m))
!      cr(2:2*nblk*m:2)=dimag(c(1:nblk*m))
!      inc=2*nblk
      inc=nblk
      do n=1,nblk
!         i=2*n-1
         i=n
         call gpfa(cr(i:),ci(i:),trig,inc,1,m,1,isign)
      enddo
!      c(1:nblk*m)=dcmplx(cr(1:2*nblk*m-1:2),cr(2:2*nblk*m:2))
      end subroutine cgpfa


!*********************************************************************
!                                                                    *
!     gpfapack - fortran implementation of the self-sorting          *
!     in-place generalized prime factor (complex) fft [gpfa]         *
!                                                                    *
!     written by clive temperton                                     *
!     recherche en prevision numerique / ecmwf                       *
!                                                                    *
!     the package consists of the setup routine setgpfa, together    *
!     with the routines gpfa, gpfa2f, gpfa3f, gpfa5f                 *
!                                                                    *
!*********************************************************************
!
!        subroutine 'setgpfa'
!        setup routine for self-sorting in-place
!            generalized prime factor (complex) fft [gpfa]
!
!        call setgpfa(trigs,n)
!
!        input :
!        -----
!        n is the length of the transforms. n must be of the form:
!          -----------------------------------
!            n = (2**ip) * (3**iq) * (5**ir)
!          -----------------------------------
!
!        output:
!        ------
!        trigs is a table of twiddle factors,
!          of length 2*ipqr (real) words, where:
!          --------------------------------------
!            ipqr = (2**ip) + (3**iq) + (5**ir)
!          --------------------------------------
!
!        written by clive temperton 1990
!
!----------------------------------------------------------------------
!
      subroutine setgpfa(trigs,n)
      implicit none
      integer :: n,nn,ifac,ll,kk,nj(3),ip,iq,ir,ni,irot,kink,k,i
      real(8) :: trigs(*),twopi,del,angle
!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do ll = 1 , 3
         kk = 0
         do while (mod(nn,ifac).eq.0)
            kk = kk + 1
            nn = nn / ifac
         enddo
         nj(ll) = kk
         ifac = ifac + ll
      enddo
!
      if (nn.ne.1) then
         write(6,40) n
   40    format(' *** warning!!!',i10,' is not a legal value of n ***')
         return
      endif
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute list of rotated twiddle factors
!     ---------------------------------------
      nj(1) = 2**ip
      nj(2) = 3**iq
      nj(3) = 5**ir
!
      twopi = 8.0d0 *datan(1.0d0)
      i = 1
!
      do ll = 1 , 3
         ni = nj(ll)
         if (ni.eq.1) cycle
!
         del = twopi / dble(ni)
         irot = n / ni
         kink = mod(irot,ni)
         kk = 0
!
         do k = 1 , ni
            angle = dble(kk) * del
            trigs(i) = cos(angle)
            trigs(i+1) = sin(angle)
            i = i + 2
            kk = kk + kink
            if (kk.gt.ni) kk = kk - ni
         enddo
      enddo
      end subroutine setgpfa
!        subroutine 'gpfa'
!        self-sorting in-place generalized prime factor (complex) fft
!
!        *** this is the all-fortran version ***
!            -------------------------------
!
!        call gpfa(a,b,trigs,inc,jump,n,lot,isign)
!
!        a is first real input/output vector
!        b is first imaginary input/output vector
!        trigs is a table of twiddle factors, precalculated
!              by calling subroutine 'setgpfa'
!        inc is the increment within each data vector
!        jump is the increment between data vectors
!        n is the length of the transforms:
!          -----------------------------------
!            n = (2**ip) * (3**iq) * (5**ir)
!          -----------------------------------
!        lot is the number of transforms
!        isign = +1 for forward transform
!              = -1 for inverse transform
!
!        written by clive temperton
!        recherche en prevision numerique
!        atmospheric environment service, canada
!
!----------------------------------------------------------------------
!
!        definition of transform
!        -----------------------
!
!        x(j) = sum(k=0,...,n-1)(c(k)*exp(isign*2*i*j*k*pi/n))
!
!---------------------------------------------------------------------
!
!        for a mathematical development of the algorithm used,
!        see:
!
!        c temperton : "a generalized prime factor fft algorithm
!          for any n = (2**p)(3**q)(5**r)",
!          siam j. sci. stat. comp., may 1992.
!
!----------------------------------------------------------------------
!
      subroutine gpfa(a,b,trigs,inc,jump,n,lot,isign)
      implicit none
      integer :: inc,jump,n,lot,isign,nn,ifac,ll,kk,nj(3),ip,iq,ir,i
      real(8) :: a(*),b(*),trigs(*)
!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do ll = 1 , 3
         kk = 0
         do while (mod(nn,ifac).eq.0)
            kk = kk + 1
            nn = nn / ifac
         enddo
         nj(ll) = kk
         ifac = ifac + ll
      enddo
!
      if (nn.ne.1) then
         write(6,40) n
   40    format(' *** warning!!!',i10,' is not a legal value of n ***')
         return
      endif
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute the transform
!     ---------------------
      i = 1
      if (ip.gt.0) then
         call gpfa2f(a,b,trigs,inc,jump,n,ip,lot,isign)
         i = i + 2 * ( 2**ip)
      endif
      if (iq.gt.0) then
         call gpfa3f(a,b,trigs(i),inc,jump,n,iq,lot,isign)
         i = i + 2 * (3**iq)
      endif
      if (ir.gt.0) then
         call gpfa5f(a,b,trigs(i),inc,jump,n,ir,lot,isign)
      endif
!
      end subroutine gpfa
!     fortran version of *gpfa2* -
!     radix-2 section of self-sorting, in-place, generalized pfa
!     central radix-2 and radix-8 passes included
!      so that transform length can be any power of 2
!
!-------------------------------------------------------------------
!
      subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)
      implicit none
      integer :: inc,jump,n,mm,lot,isign,lvr,n2,inq,jstepx,ninc,ink, &
                 m2,m8,m,mh,nblox,left,nb,nvex,la,mu,ipass,jstep,jstepl, &
                 jjj,ja,nu,jb,jc,jd,j,l,kk,k,je,jf,jg,jh,laincl,ji,jj,jk,&
                 jl,jm,jn,jo,jp,istart,ll
      real(8) :: a(*),b(*),trigs(*),s,ss,t0,t2,t1,t3,u0,u2,u1,u3,co1,si1, &
                 co2,si2,co3,si3,c1,c2,c3,co4,si4,co5,si5,co6,si6,co7,si7
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n2 = 2**mm
      inq = n/n2
      jstepx = (n2-n) * inc
      ninc = n * inc
      ink = inc * inq
!
      m2 = 0
      m8 = 0
      if (mod(mm,2).eq.0) then
         m = mm/2
      else if (mod(mm,4).eq.1) then
         m = (mm-1)/2
         m2 = 1
      else if (mod(mm,4).eq.3) then
         m = (mm-3)/2
         m8 = 1
      endif
      mh = (m+1)/2
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = dble(isign)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
         if (left.le.lvr) then
            nvex = left
         else if (left.lt.(2*lvr)) then
            nvex = left/2
            nvex = nvex + mod(nvex,2)
         else
            nvex = lvr
         endif
         left = left - nvex
!
         la = 1
!
!  loop on type i radix-4 passes
!  -----------------------------
         mu = mod(inq,4)
         if (isign.eq.-1) mu = 4 - mu
         ss = 1.0
         if (mu.eq.3) ss = -1.0
!
         if (mh.eq.0) go to 200
!
         do 160 ipass = 1 , mh
            jstep = (n*inc) / (4*la)
            jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
            do 120 jjj = 0 , (n-1)*inc , 4*jstep
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do 115 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc.lt.istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd.lt.istart) jd = jd + ninc
                  j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                  do 110 l = 1 , nvex
                     t0 = a(ja+j) + a(jc+j)
                     t2 = a(ja+j) - a(jc+j)
                     t1 = a(jb+j) + a(jd+j)
                     t3 = ss * ( a(jb+j) - a(jd+j) )
                     u0 = b(ja+j) + b(jc+j)
                     u2 = b(ja+j) - b(jc+j)
                     u1 = b(jb+j) + b(jd+j)
                     u3 = ss * ( b(jb+j) - b(jd+j) )
                     a(ja+j) = t0 + t1
                     a(jc+j) = t0 - t1
                     b(ja+j) = u0 + u1
                     b(jc+j) = u0 - u1
                     a(jb+j) = t2 - u3
                     a(jd+j) = t2 + u3
                     b(jb+j) = u2 + t3
                     b(jd+j) = u2 - t3
                     j = j + jump
  110             continue
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
  115          continue
  120       continue
!
!  finished if n2 = 4
!  ------------------

            if (n2.eq.4) go to 490
            kk = 2 * la
!
!  loop on nonzero k
!  -----------------
            do 150 k = ink , jstep-ink , ink
               co1 = trigs(kk+1)
               si1 = s*trigs(kk+2)
               co2 = trigs(2*kk+1)
               si2 = s*trigs(2*kk+2)
               co3 = trigs(3*kk+1)
               si3 = s*trigs(3*kk+2)
!
!  loop along transform
!  --------------------
               do 140 jjj = k , (n-1)*inc , 4*jstep
                  ja = istart + jjj
!
!     "transverse" loop
!     -----------------
                  do 135 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd.lt.istart) jd = jd + ninc
                     j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep,shortloop
                     do 130 l = 1 , nvex
                        t0 = a(ja+j) + a(jc+j)
                        t2 = a(ja+j) - a(jc+j)
                        t1 = a(jb+j) + a(jd+j)
                        t3 = ss * ( a(jb+j) - a(jd+j ) )
                        u0 = b(ja+j) + b(jc+j)
                        u2 = b(ja+j) - b(jc+j)
                        u1 = b(jb+j) + b(jd+j)
                        u3 = ss * ( b(jb+j) - b(jd+j) )
                        a(ja+j) = t0 + t1
                        b(ja+j) = u0 + u1
                        a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                        b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                        a(jc+j) = co2*(t0-t1) - si2*(u0-u1)
                        b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
                        a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
                        b(jd+j) = si3*(t2+u3) + co3*(u2-t3)
                        j = j + jump
  130                continue
!-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
  135             continue
  140          continue
!-----( end of loop along transforms )
               kk = kk + 2*la
  150       continue
!-----( end of loop on nonzero k )
            la = 4*la
  160    continue
!-----( end of loop on type i radix-4 passes)
!
!  central radix-2 pass
!  --------------------
  200 continue
      if (m2.eq.0) go to 300
!
      jstep = (n*inc) / (2*la)
      jstepl = jstep - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 220 jjj = 0 , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 215 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 210 l = 1 , nvex
      t0 = a(ja+j) - a(jb+j)
      a(ja+j) = a(ja+j) + a(jb+j)
      a(jb+j) = t0
      u0 = b(ja+j) - b(jb+j)
      b(ja+j) = b(ja+j) + b(jb+j)
      b(jb+j) = u0
      j = j + jump
  210 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  215 continue
  220 continue
!
!  finished if n2=2
!  ----------------
      if (n2.eq.2) go to 490
!
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 260 k = ink , jstep - ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
!
!  loop along transforms
!  ---------------------
      do 250 jjj = k , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 245 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (kk.eq.n2/2) then
!cdir$ ivdep, shortloop
      do 230 l = 1 , nvex
      t0 = ss * ( a(ja+j) - a(jb+j) )
      a(ja+j) = a(ja+j) + a(jb+j)
      a(jb+j) = ss * ( b(jb+j) - b(ja+j) )
      b(ja+j) = b(ja+j) + b(jb+j)
      b(jb+j) = t0
      j = j + jump
  230 continue
!
      else
!
!cdir$ ivdep, shortloop
      do 240 l = 1 , nvex
      t0 = a(ja+j) - a(jb+j)
      a(ja+j) = a(ja+j) + a(jb+j)
      u0 = b(ja+j) - b(jb+j)
      b(ja+j) = b(ja+j) + b(jb+j)
      a(jb+j) = co1*t0 - si1*u0
      b(jb+j) = si1*t0 + co1*u0
      j = j + jump
  240 continue
!
      endif
!
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  245 continue
  250 continue
!-----(end of loop along transforms)
      kk = kk + 2 * la
  260 continue
!-----(end of loop on nonzero k)
!-----(end of radix-2 pass)
!
      la = 2 * la
      go to 400
!
!  central radix-8 pass
!  --------------------

  300 continue
      if (m8.eq.0) go to 400
      jstep = (n*inc) / (8*la)
      jstepl = jstep - ninc
      mu = mod(inq,8)
      if (isign.eq.-1) mu = 8 - mu
      c1 = 1.0
      if (mu.eq.3.or.mu.eq.7) c1 = -1.0
      c2 = sqrt(0.5)
      if (mu.eq.3.or.mu.eq.5) c2 = -c2
      c3 = c1 * c2
!
!  stage 1
!  -------
      do 320 k = 0 , jstep - ink , ink
      do 315 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 312 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
!cdir$ ivdep, shortloop
      do 310 l = 1 , nvex
      t0 = a(ja+j) - a(je+j)
      a(ja+j) = a(ja+j) + a(je+j)
      t1 = c1 * ( a(jc+j) - a(jg+j) )
      a(je+j) = a(jc+j) + a(jg+j)
      t2 = a(jb+j) - a(jf+j)
      a(jc+j) = a(jb+j) + a(jf+j)
      t3 = a(jd+j) - a(jh+j)
      a(jg+j) = a(jd+j) + a(jh+j)
      a(jb+j) = t0
      a(jf+j) = t1
      a(jd+j) = c2 * ( t2 - t3 )
      a(jh+j) = c3 * ( t2 + t3 )
      u0 = b(ja+j) - b(je+j)
      b(ja+j) = b(ja+j) + b(je+j)
      u1 = c1 * ( b(jc+j) - b(jg+j) )
      b(je+j) = b(jc+j) + b(jg+j)
      u2 = b(jb+j) - b(jf+j)
      b(jc+j) = b(jb+j) + b(jf+j)
      u3 = b(jd+j) - b(jh+j)
      b(jg+j) = b(jd+j) + b(jh+j)
      b(jb+j) = u0
      b(jf+j) = u1
      b(jd+j) = c2 * ( u2 - u3 )
      b(jh+j) = c3 * ( u2 + u3 )
      j = j + jump
  310 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  312 continue
  315 continue
  320 continue
!
!  stage 2
!  -------
!
!  k=0 (no twiddle factors)
!  ------------------------
      do 330 jjj = 0 , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 328 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
!cdir$ ivdep, shortloop
      do 325 l = 1 , nvex
      t0 = a(ja+j) + a(je+j)
      t2 = a(ja+j) - a(je+j)
      t1 = a(jc+j) + a(jg+j)
      t3 = c1 * ( a(jc+j) - a(jg+j) )
      u0 = b(ja+j) + b(je+j)
      u2 = b(ja+j) - b(je+j)
      u1 = b(jc+j) + b(jg+j)
      u3 = c1 * ( b(jc+j) - b(jg+j ) )
      a(ja+j) = t0 + t1
      a(je+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(je+j) = u0 - u1
      a(jc+j) = t2 - u3
      a(jg+j) = t2 + u3
      b(jc+j) = u2 + t3
      b(jg+j) = u2 - t3
      t0 = a(jb+j) + a(jd+j)
      t2 = a(jb+j) - a(jd+j)
      t1 = a(jf+j) - a(jh+j)
      t3 = a(jf+j) + a(jh+j)
      u0 = b(jb+j) + b(jd+j)
      u2 = b(jb+j) - b(jd+j)
      u1 = b(jf+j) - b(jh+j)
      u3 = b(jf+j) + b(jh+j)
      a(jb+j) = t0 - u3
      a(jh+j) = t0 + u3
      b(jb+j) = u0 + t3
      b(jh+j) = u0 - t3
      a(jd+j) = t2 + u1
      a(jf+j) = t2 - u1
      b(jd+j) = u2 - t1
      b(jf+j) = u2 + t1
      j = j + jump
  325 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  328 continue
  330 continue
!
      if (n2.eq.8) go to 490
!
!  loop on nonzero k
!  -----------------
      kk = 2 * la
!
      do 350 k = ink , jstep - ink , ink
!
      co1 = trigs(kk+1)
      si1 = s * trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s * trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s * trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s * trigs(4*kk+2)
      co5 = trigs(5*kk+1)
      si5 = s * trigs(5*kk+2)
      co6 = trigs(6*kk+1)
      si6 = s * trigs(6*kk+2)
      co7 = trigs(7*kk+1)
      si7 = s * trigs(7*kk+2)
!
      do 345 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 342 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
!cdir$ ivdep, shortloop
      do 340 l = 1 , nvex
      t0 = a(ja+j) + a(je+j)
      t2 = a(ja+j) - a(je+j)
      t1 = a(jc+j) + a(jg+j)
      t3 = c1 * ( a(jc+j) - a(jg+j) )
      u0 = b(ja+j) + b(je+j)
      u2 = b(ja+j) - b(je+j)
      u1 = b(jc+j) + b(jg+j)
      u3 = c1 * ( b(jc+j) - b(jg+j ) )
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(je+j) = co4*(t0-t1) - si4*(u0-u1)
      b(je+j) = si4*(t0-t1) + co4*(u0-u1)
      a(jc+j) = co2*(t2-u3) - si2*(u2+t3)
      b(jc+j) = si2*(t2-u3) + co2*(u2+t3)
      a(jg+j) = co6*(t2+u3) - si6*(u2-t3)
      b(jg+j) = si6*(t2+u3) + co6*(u2-t3)
      t0 = a(jb+j) + a(jd+j)
      t2 = a(jb+j) - a(jd+j)
      t1 = a(jf+j) - a(jh+j)
      t3 = a(jf+j) + a(jh+j)
      u0 = b(jb+j) + b(jd+j)
      u2 = b(jb+j) - b(jd+j)
      u1 = b(jf+j) - b(jh+j)
      u3 = b(jf+j) + b(jh+j)
      a(jb+j) = co1*(t0-u3) - si1*(u0+t3)
      b(jb+j) = si1*(t0-u3) + co1*(u0+t3)
      a(jh+j) = co7*(t0+u3) - si7*(u0-t3)
      b(jh+j) = si7*(t0+u3) + co7*(u0-t3)
      a(jd+j) = co3*(t2+u1) - si3*(u2-t1)
      b(jd+j) = si3*(t2+u1) + co3*(u2-t1)
      a(jf+j) = co5*(t2-u1) - si5*(u2+t1)
      b(jf+j) = si5*(t2-u1) + co5*(u2+t1)
      j = j + jump
  340 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  342 continue
  345 continue
      kk = kk + 2 * la
  350 continue
!
      la = 8 * la
!
!  loop on type ii radix-4 passes
!  ------------------------------
  400 continue
      mu = mod(inq,4)
      if (isign.eq.-1) mu = 4 - mu
      ss = 1.0
      if (mu.eq.3) ss = -1.0
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 4*jstep
!
      do 420 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      a(jc+j) = a(ji+j)
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(jb+j) = a(je+j)
      a(ja+j) = t0 + t1
      a(ji+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(jc+j) = u0 - u1
      b(jd+j) = b(jm+j)
      a(je+j) = t2 - u3
      a(jd+j) = t2 + u3
      b(jb+j) = u2 + t3
      b(jm+j) = u2 - t3
!----------------------
      t0 = a(jb+j) + a(jg+j)
      t2 = a(jb+j) - a(jg+j)
      t1 = a(jf+j) + a(jh+j)
      t3 = ss * ( a(jf+j) - a(jh+j) )
      a(jg+j) = a(jj+j)
      u0 = b(je+j) + b(jg+j)
      u2 = b(je+j) - b(jg+j)
      u1 = b(jf+j) + b(jh+j)
      u3 = ss * ( b(jf+j) - b(jh+j) )
      b(je+j) = b(jb+j)
      a(jb+j) = t0 + t1
      a(jj+j) = t0 - t1
      b(jg+j) = b(jj+j)
      b(jb+j) = u0 + u1
      b(jj+j) = u0 - u1
      a(jf+j) = t2 - u3
      a(jh+j) = t2 + u3
      b(jf+j) = u2 + t3
      b(jh+j) = u2 - t3
!----------------------
      t0 = a(jc+j) + a(jk+j)
      t2 = a(jc+j) - a(jk+j)
      t1 = a(jg+j) + a(jl+j)
      t3 = ss * ( a(jg+j) - a(jl+j) )
      u0 = b(ji+j) + b(jk+j)
      u2 = b(ji+j) - b(jk+j)
      a(jl+j) = a(jo+j)
      u1 = b(jg+j) + b(jl+j)
      u3 = ss * ( b(jg+j) - b(jl+j) )
      b(ji+j) = b(jc+j)
      a(jc+j) = t0 + t1
      a(jk+j) = t0 - t1
      b(jl+j) = b(jo+j)
      b(jc+j) = u0 + u1
      b(jk+j) = u0 - u1
      a(jg+j) = t2 - u3
      a(jo+j) = t2 + u3
      b(jg+j) = u2 + t3
      b(jo+j) = u2 - t3
!----------------------
      t0 = a(jm+j) + a(jl+j)
      t2 = a(jm+j) - a(jl+j)
      t1 = a(jn+j) + a(jp+j)
      t3 = ss * ( a(jn+j) - a(jp+j) )
      a(jm+j) = a(jd+j)
      u0 = b(jd+j) + b(jl+j)
      u2 = b(jd+j) - b(jl+j)
      u1 = b(jn+j) + b(jp+j)
      u3 = ss * ( b(jn+j) - b(jp+j) )
      a(jn+j) = a(jh+j)
      a(jd+j) = t0 + t1
      a(jl+j) = t0 - t1
      b(jd+j) = u0 + u1
      b(jl+j) = u0 - u1
      b(jn+j) = b(jh+j)
      a(jh+j) = t2 - u3
      a(jp+j) = t2 + u3
      b(jh+j) = u2 + t3
      b(jp+j) = u2 - t3
      j = j + jump
  410 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass.eq.m) go to 490
!
      kk = 2*la
!
!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 4*jstep
!
      do 450 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      a(jc+j) = a(ji+j)
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(jb+j) = a(je+j)
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(je+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      b(jd+j) = b(jm+j)
      a(ji+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jm+j) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      t0 = a(jb+j) + a(jg+j)
      t2 = a(jb+j) - a(jg+j)
      t1 = a(jf+j) + a(jh+j)
      t3 = ss * ( a(jf+j) - a(jh+j) )
      a(jg+j) = a(jj+j)
      u0 = b(je+j) + b(jg+j)
      u2 = b(je+j) - b(jg+j)
      u1 = b(jf+j) + b(jh+j)
      u3 = ss * ( b(jf+j) - b(jh+j) )
      b(je+j) = b(jb+j)
      a(jb+j) = t0 + t1
      b(jb+j) = u0 + u1
      b(jg+j) = b(jj+j)
      a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jj+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jj+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jh+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jh+j) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      t0 = a(jc+j) + a(jk+j)
      t2 = a(jc+j) - a(jk+j)
      t1 = a(jg+j) + a(jl+j)
      t3 = ss * ( a(jg+j) - a(jl+j) )
      u0 = b(ji+j) + b(jk+j)
      u2 = b(ji+j) - b(jk+j)
      a(jl+j) = a(jo+j)
      u1 = b(jg+j) + b(jl+j)
      u3 = ss * ( b(jg+j) - b(jl+j) )
      b(ji+j) = b(jc+j)
      a(jc+j) = t0 + t1
      b(jc+j) = u0 + u1
      b(jl+j) = b(jo+j)
      a(jg+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jg+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jk+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jk+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jo+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jo+j) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      t0 = a(jm+j) + a(jl+j)
      t2 = a(jm+j) - a(jl+j)
      t1 = a(jn+j) + a(jp+j)
      t3 = ss * ( a(jn+j) - a(jp+j) )
      a(jm+j) = a(jd+j)
      u0 = b(jd+j) + b(jl+j)
      u2 = b(jd+j) - b(jl+j)
      a(jn+j) = a(jh+j)
      u1 = b(jn+j) + b(jp+j)
      u3 = ss * ( b(jn+j) - b(jp+j) )
      b(jn+j) = b(jh+j)
      a(jd+j) = t0 + t1
      b(jd+j) = u0 + u1
      a(jh+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jh+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jl+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jl+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jp+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jp+j) = si3*(t2+u3) + co3*(u2-t3)
      j = j + jump
  440 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 4*la
  480 continue
!-----( end of loop on type ii radix-4 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa2f
!     fortran version of *gpfa3* -
!     radix-3 section of self-sorting, in-place
!        generalized pfa
!
!-------------------------------------------------------------------
!
      subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)
      implicit none
      integer :: inc,jump,n,mm,lot,isign,lvr,inq,jstepx,ninc,ink, &
                 m,mh,nblox,left,nb,nvex,la,mu,ipass,jstep,jstepl, &
                 jjj,ja,nu,jb,jc,jd,j,l,kk,k,je,jf,jg,jh,laincl,ji, &
                 n3,istart,ll
      real(8) :: a(*),b(*),trigs(*),s,t2,t1,t3,u2,u1,u3,co1,si1, &
                 co2,si2,c1, &
                 sin60
      data sin60/0.866025403784437/
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n3 = 3**mm
      inq = n/n3
      jstepx = (n3-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,3)
      if (isign.eq.-1) mu = 3-mu
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = sin60
      if (mu.eq.2) c1 = -c1
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = float(isign)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type i radix-3 passes
!  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 3*jstep
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 115 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jb+j) = t2 - u3
      b(jb+j) = u2 + t3
      a(jc+j) = t2 + u3
      b(jc+j) = u2 - t3
      j = j + jump
  110 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  115 continue
  120 continue
!
!  finished if n3 = 3
!  ------------------
      if (n3.eq.3) go to 490
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 150 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 3*jstep
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  130 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 3*la
  160 continue
!-----( end of loop on type i radix-3 passes)
!
!  loop on type ii radix-3 passes
!  ------------------------------
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
      laincl = la*ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 3*jstep
!
      do 420 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      a(jb+j) = a(jd+j)
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      b(jb+j) = b(jd+j)
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jd+j) = t2 - u3
      b(jd+j) = u2 + t3
      a(jc+j) = t2 + u3
      b(jc+j) = u2 - t3
!----------------------
      t1 = a(je+j) + a(jf+j)
      t2 = a(jb+j) - 0.5 * t1
      t3 = c1 * ( a(je+j) - a(jf+j) )
      a(jf+j) = a(jh+j)
      u1 = b(je+j) + b(jf+j)
      u2 = b(jb+j) - 0.5 * u1
      u3 = c1 * ( b(je+j) - b(jf+j) )
      b(jf+j) = b(jh+j)
      a(jb+j) = a(jb+j) + t1
      b(jb+j) = b(jb+j) + u1
      a(je+j) = t2 - u3
      b(je+j) = u2 + t3
      a(jh+j) = t2 + u3
      b(jh+j) = u2 - t3
!----------------------
      t1 = a(jf+j) + a(ji+j)
      t2 = a(jg+j) - 0.5 * t1
      t3 = c1 * ( a(jf+j) - a(ji+j) )
      t1 = a(jg+j) + t1
      a(jg+j) = a(jc+j)
      u1 = b(jf+j) + b(ji+j)
      u2 = b(jg+j) - 0.5 * u1
      u3 = c1 * ( b(jf+j) - b(ji+j) )
      u1 = b(jg+j) + u1
      b(jg+j) = b(jc+j)
      a(jc+j) = t1
      b(jc+j) = u1
      a(jf+j) = t2 - u3
      b(jf+j) = u2 + t3
      a(ji+j) = t2 + u3
      b(ji+j) = u2 - t3
      j = j + jump
  410 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass.eq.m) go to 490
!
      kk = 2*la
!
!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 3*jstep
!
      do 450 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      a(jb+j) = a(jd+j)
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      b(jb+j) = b(jd+j)
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jd+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jd+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
!----------------------
      t1 = a(je+j) + a(jf+j)
      t2 = a(jb+j) - 0.5 * t1
      t3 = c1 * ( a(je+j) - a(jf+j) )
      a(jf+j) = a(jh+j)
      u1 = b(je+j) + b(jf+j)
      u2 = b(jb+j) - 0.5 * u1
      u3 = c1 * ( b(je+j) - b(jf+j) )
      b(jf+j) = b(jh+j)
      a(jb+j) = a(jb+j) + t1
      b(jb+j) = b(jb+j) + u1
      a(je+j) = co1*(t2-u3) - si1*(u2+t3)
      b(je+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jh+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jh+j) = si2*(t2+u3) + co2*(u2-t3)
!----------------------
      t1 = a(jf+j) + a(ji+j)
      t2 = a(jg+j) - 0.5 * t1
      t3 = c1 * ( a(jf+j) - a(ji+j) )
      t1 = a(jg+j) + t1
      a(jg+j) = a(jc+j)
      u1 = b(jf+j) + b(ji+j)
      u2 = b(jg+j) - 0.5 * u1
      u3 = c1 * ( b(jf+j) - b(ji+j) )
      u1 = b(jg+j) + u1
      b(jg+j) = b(jc+j)
      a(jc+j) = t1
      b(jc+j) = u1
      a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
      a(ji+j) = co2*(t2+u3) - si2*(u2-t3)
      b(ji+j) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  440 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 3*la
  480 continue
!-----( end of loop on type ii radix-3 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa3f
!     fortran version of *gpfa5* -
!     radix-5 section of self-sorting, in-place,
!        generalized pfa
!
!-------------------------------------------------------------------
!
      subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)
      implicit none
      integer :: inc,jump,n,mm,lot,isign,lvr,inq,jstepx,ninc,ink, &
                 m,mh,nblox,left,nb,nvex,la,mu,ipass,jstep,jstepl, &
                 jjj,ja,nu,jb,jc,jd,j,l,kk,k,je,jf,jg,jh,laincl,ji,jj,jk,&
                 jl,jm,jn,jo,jp,n5,jq,jr,js,jt,ju,jv,jw,jx,jy,istart,ll
      real(8) :: a(*),b(*),trigs(*),s,t2,t1,t3,u2,u1,u3,co1,si1, &
                 co2,si2,co3,si3,c1,c2,c3,co4,si4, &
                 sin36,sin72,qrt5,t4,t5,t6,t7,t8,t9,t10,t11, &
                 ax,bx,u4,u5,u6,u7,u8,u9,u10,u11
      data sin36/0.587785252292473/, sin72/0.951056516295154/, &
           qrt5/0.559016994374947/
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n5 = 5 ** mm
      inq = n / n5
      jstepx = (n5-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,5)
      if (isign.eq.-1) mu = 5 - mu
!
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = qrt5
      c2 = sin72
      c3 = sin36
      if (mu.eq.2.or.mu.eq.3) then
         c1 = -c1
         c2 = sin36
         c3 = sin72
      endif
      if (mu.eq.3.or.mu.eq.4) c2 = -c2
      if (mu.eq.2.or.mu.eq.4) c3 = -c3
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = float(isign)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type i radix-5 passes
!  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      kk = 0
!
!  loop on k
!  ---------
      do 150 k = 0 , jstep-ink , ink
!
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 5*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (k.eq.0) then
!
!cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jb+j) = t8 - u11
      b(jb+j) = u8 + t11
      a(je+j) = t8 + u11
      b(je+j) = u8 - t11
      a(jc+j) = t9 - u10
      b(jc+j) = u9 + t10
      a(jd+j) = t9 + u10
      b(jd+j) = u9 - t10
      j = j + jump
  110 continue
!
      else
!
!cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jb+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jb+j) = si1*(t8-u11) + co1*(u8+t11)
      a(je+j) = co4*(t8+u11) - si4*(u8-t11)
      b(je+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jc+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jc+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  130 continue
!
      endif
!
!-----( end of loop across transforms )
!
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 5*la
  160 continue
!-----( end of loop on type i radix-5 passes)
!
      if (n.eq.5) go to 490
!
!  loop on type ii radix-5 passes
!  ------------------------------
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
      kk = 0
!
!     loop on k
!     ---------
      do 470 k = 0 , jstep-ink , ink
!
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 5*jstep
!
      do 450 jjj = ll , (n-1)*inc , 5*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = ja + laincl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jf + laincl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = jl + jstepl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jk + laincl
      if (jp.lt.istart) jp = jp + ninc
      jq = jp + jstepl
      if (jq.lt.istart) jq = jq + ninc
      jr = jq + jstepl
      if (jr.lt.istart) jr = jr + ninc
      js = jr + jstepl
      if (js.lt.istart) js = js + ninc
      jt = js + jstepl
      if (jt.lt.istart) jt = jt + ninc
      ju = jp + laincl
      if (ju.lt.istart) ju = ju + ninc
      jv = ju + jstepl
      if (jv.lt.istart) jv = jv + ninc
      jw = jv + jstepl
      if (jw.lt.istart) jw = jw + ninc
      jx = jw + jstepl
      if (jx.lt.istart) jx = jx + ninc
      jy = jx + jstepl
      if (jy.lt.istart) jy = jy + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (k.eq.0) then
!
!cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      a(jb+j) = a(jf+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jc+j) = a(jk+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      b(jb+j) = b(jf+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jc+j) = b(jk+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jf+j) = t8 - u11
      b(jf+j) = u8 + t11
      a(je+j) = t8 + u11
      b(je+j) = u8 - t11
      a(jk+j) = t9 - u10
      b(jk+j) = u9 + t10
      a(jd+j) = t9 + u10
      b(jd+j) = u9 - t10
!----------------------
      t1 = a(jg+j) + a(jj+j)
      t2 = a(jh+j) + a(ji+j)
      t3 = a(jg+j) - a(jj+j)
      t4 = a(jh+j) - a(ji+j)
      a(jh+j) = a(jl+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jb+j) - 0.25 * t5
      a(jb+j) = a(jb+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ji+j) = a(jq+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jg+j) + b(jj+j)
      u2 = b(jh+j) + b(ji+j)
      u3 = b(jg+j) - b(jj+j)
      u4 = b(jh+j) - b(ji+j)
      b(jh+j) = b(jl+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jb+j) - 0.25 * u5
      b(jb+j) = b(jb+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ji+j) = b(jq+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jg+j) = t8 - u11
      b(jg+j) = u8 + t11
      a(jj+j) = t8 + u11
      b(jj+j) = u8 - t11
      a(jl+j) = t9 - u10
      b(jl+j) = u9 + t10
      a(jq+j) = t9 + u10
      b(jq+j) = u9 - t10
!----------------------
      t1 = a(jh+j) + a(jo+j)
      t2 = a(jm+j) + a(jn+j)
      t3 = a(jh+j) - a(jo+j)
      t4 = a(jm+j) - a(jn+j)
      a(jn+j) = a(jr+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jc+j) - 0.25 * t5
      a(jc+j) = a(jc+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jo+j) = a(jw+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jh+j) + b(jo+j)
      u2 = b(jm+j) + b(jn+j)
      u3 = b(jh+j) - b(jo+j)
      u4 = b(jm+j) - b(jn+j)
      b(jn+j) = b(jr+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jc+j) - 0.25 * u5
      b(jc+j) = b(jc+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jo+j) = b(jw+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jh+j) = t8 - u11
      b(jh+j) = u8 + t11
      a(jw+j) = t8 + u11
      b(jw+j) = u8 - t11
      a(jm+j) = t9 - u10
      b(jm+j) = u9 + t10
      a(jr+j) = t9 + u10
      b(jr+j) = u9 - t10
!----------------------
      t1 = a(ji+j) + a(jt+j)
      t2 = a(jn+j) + a(js+j)
      t3 = a(ji+j) - a(jt+j)
      t4 = a(jn+j) - a(js+j)
      a(jt+j) = a(jx+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jp+j) - 0.25 * t5
      ax = a(jp+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jp+j) = a(jd+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(jd+j) = ax
      u1 = b(ji+j) + b(jt+j)
      u2 = b(jn+j) + b(js+j)
      u3 = b(ji+j) - b(jt+j)
      u4 = b(jn+j) - b(js+j)
      b(jt+j) = b(jx+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jp+j) - 0.25 * u5
      bx = b(jp+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jp+j) = b(jd+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(jd+j) = bx
      a(ji+j) = t8 - u11
      b(ji+j) = u8 + t11
      a(jx+j) = t8 + u11
      b(jx+j) = u8 - t11
      a(jn+j) = t9 - u10
      b(jn+j) = u9 + t10
      a(js+j) = t9 + u10
      b(js+j) = u9 - t10
!----------------------
      t1 = a(jv+j) + a(jy+j)
      t2 = a(jo+j) + a(jt+j)
      t3 = a(jv+j) - a(jy+j)
      t4 = a(jo+j) - a(jt+j)
      a(jv+j) = a(jj+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ju+j) - 0.25 * t5
      ax = a(ju+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ju+j) = a(je+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(je+j) = ax
      u1 = b(jv+j) + b(jy+j)
      u2 = b(jo+j) + b(jt+j)
      u3 = b(jv+j) - b(jy+j)
      u4 = b(jo+j) - b(jt+j)
      b(jv+j) = b(jj+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ju+j) - 0.25 * u5
      bx = b(ju+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ju+j) = b(je+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(je+j) = bx
      a(jj+j) = t8 - u11
      b(jj+j) = u8 + t11
      a(jy+j) = t8 + u11
      b(jy+j) = u8 - t11
      a(jo+j) = t9 - u10
      b(jo+j) = u9 + t10
      a(jt+j) = t9 + u10
      b(jt+j) = u9 - t10
      j = j + jump
  410 continue
!
      else
!
!cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      a(jb+j) = a(jf+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jc+j) = a(jk+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      b(jb+j) = b(jf+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jc+j) = b(jk+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jf+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jf+j) = si1*(t8-u11) + co1*(u8+t11)
      a(je+j) = co4*(t8+u11) - si4*(u8-t11)
      b(je+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jk+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jk+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(jg+j) + a(jj+j)
      t2 = a(jh+j) + a(ji+j)
      t3 = a(jg+j) - a(jj+j)
      t4 = a(jh+j) - a(ji+j)
      a(jh+j) = a(jl+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jb+j) - 0.25 * t5
      a(jb+j) = a(jb+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ji+j) = a(jq+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jg+j) + b(jj+j)
      u2 = b(jh+j) + b(ji+j)
      u3 = b(jg+j) - b(jj+j)
      u4 = b(jh+j) - b(ji+j)
      b(jh+j) = b(jl+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jb+j) - 0.25 * u5
      b(jb+j) = b(jb+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ji+j) = b(jq+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jg+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jg+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jj+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jj+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jl+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jl+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jq+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jq+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(jh+j) + a(jo+j)
      t2 = a(jm+j) + a(jn+j)
      t3 = a(jh+j) - a(jo+j)
      t4 = a(jm+j) - a(jn+j)
      a(jn+j) = a(jr+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jc+j) - 0.25 * t5
      a(jc+j) = a(jc+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jo+j) = a(jw+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jh+j) + b(jo+j)
      u2 = b(jm+j) + b(jn+j)
      u3 = b(jh+j) - b(jo+j)
      u4 = b(jm+j) - b(jn+j)
      b(jn+j) = b(jr+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jc+j) - 0.25 * u5
      b(jc+j) = b(jc+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jo+j) = b(jw+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jh+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jh+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jw+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jw+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jm+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jm+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jr+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jr+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(ji+j) + a(jt+j)
      t2 = a(jn+j) + a(js+j)
      t3 = a(ji+j) - a(jt+j)
      t4 = a(jn+j) - a(js+j)
      a(jt+j) = a(jx+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jp+j) - 0.25 * t5
      ax = a(jp+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jp+j) = a(jd+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(jd+j) = ax
      u1 = b(ji+j) + b(jt+j)
      u2 = b(jn+j) + b(js+j)
      u3 = b(ji+j) - b(jt+j)
      u4 = b(jn+j) - b(js+j)
      b(jt+j) = b(jx+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jp+j) - 0.25 * u5
      bx = b(jp+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jp+j) = b(jd+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(jd+j) = bx
      a(ji+j) = co1*(t8-u11) - si1*(u8+t11)
      b(ji+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jx+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jx+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jn+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jn+j) = si2*(t9-u10) + co2*(u9+t10)
      a(js+j) = co3*(t9+u10) - si3*(u9-t10)
      b(js+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(jv+j) + a(jy+j)
      t2 = a(jo+j) + a(jt+j)
      t3 = a(jv+j) - a(jy+j)
      t4 = a(jo+j) - a(jt+j)
      a(jv+j) = a(jj+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ju+j) - 0.25 * t5
      ax = a(ju+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ju+j) = a(je+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(je+j) = ax
      u1 = b(jv+j) + b(jy+j)
      u2 = b(jo+j) + b(jt+j)
      u3 = b(jv+j) - b(jy+j)
      u4 = b(jo+j) - b(jt+j)
      b(jv+j) = b(jj+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ju+j) - 0.25 * u5
      bx = b(ju+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ju+j) = b(je+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(je+j) = bx
      a(jj+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jj+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jy+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jy+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jo+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jo+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jt+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jt+j) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  440 continue
!
      endif
!
!-----(end of loop across transforms)
!
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 5*la
  480 continue
!-----( end of loop on type ii radix-5 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa5f
      end module fft
