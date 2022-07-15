!
!  special function for the multiple sphere problem
!
      module specialfuncs
      implicit none
      real(8), private :: pi=3.141592653589793d0
      contains

         subroutine timewrite(iunit,char1,time,line_break)
         use intrinsics
         implicit none
         integer :: iunit
         real(8) :: time,time2
         logical :: linebreak
         logical, optional :: line_break
         character(*) :: char1
         if(present(line_break)) then
            linebreak=line_break
         else
            linebreak=.true.
         endif
         if(time.gt.3600.d0) then
            time2=time/3600.d0
            if(linebreak) then
               write(iunit,'(a,f9.3,'' hours'')') char1,time2
            else
               write(iunit,'(a,f9.3,'' hours'',$)') char1,time2
            endif
         elseif(time.gt.60.d0) then
            time2=time/60.d0
            if(linebreak) then
               write(iunit,'(a,f9.2,'' min'')') char1,time2
            else
               write(iunit,'(a,f9.2,'' min'',$)') char1,time2
            endif
         else
            if(linebreak) then
               write(iunit,'(a,f9.2,'' sec'')') char1,time
            else
               write(iunit,'(a,f9.2,'' sec'',$)') char1,time
            endif
         endif
         if(linebreak) call flush(iunit)
         end subroutine timewrite

         subroutine mieext(x,refrel,qext,qsca,gsca)
         implicit none
         real(8) ::  gsca,qext,qsca,x
         complex(8) :: refrel
         integer :: n,nstop,nmx,nn
         real(8) ::  chi,chi0,chi1,dx,en,fn,psi,psi0,psi1, &
                     xstop,ymod
         complex(8) :: an,an1,bn,bn1,drefrl,xi,xi1,y
         complex(8), allocatable :: d(:)

         dx=x
         drefrl=refrel
         y=x*drefrl
         ymod=abs(y)
         xstop=x+4.*x**0.3333+2.
         nmx=nint(max(xstop,ymod))+15
         allocate(d(nmx))
         nstop=nint(xstop)
         d(nmx)=(0.,0.)
         nn=nmx-1
         do n=1,nn
            en=nmx-n+1
            d(nmx-n)=(en/y)-(1./(d(nmx-n+1)+en/y))
         enddo
         psi0=cos(dx)
         psi1=sin(dx)
         chi0=-sin(dx)
         chi1=cos(dx)
         xi1=dcmplx(psi1,-chi1)
         qsca=0.
         gsca=0.
         qext=0.d0
         do n=1,nstop
            en=n
            fn=(2.e0*en+1.)/(en*(en+1.))
            psi=(2.e0*en-1.)*psi1/dx-psi0
            chi=(2.e0*en-1.)*chi1/dx-chi0
            xi=dcmplx(psi,-chi)
            if(n.gt.1)then
               an1=an
               bn1=bn
            endif
            an=(d(n)/drefrl+en/dx)*psi-psi1
            an=an/((d(n)/drefrl+en/dx)*xi-xi1)
            bn=(drefrl*d(n)+en/dx)*psi-psi1
            bn=bn/((drefrl*d(n)+en/dx)*xi-xi1)
            qsca=qsca+dble((2.*en+1.)*(abs(an)**2+abs(bn)**2))
            qext=qext+dble((2.*en+1.)*(an+bn))
            gsca=gsca+dble(((2.*en+1.)/(en*(en+1.))) &
                 *(dble(an)*dble(bn)+dimag(an)*dimag(bn)))
            if(n.gt.1) then
               gsca=gsca+dble(((en-1.)*(en+1.)/en) &
                 *(dble(an1)*dble(an)+dimag(an1)*dimag(an) &
                 +dble(bn1)*dble(bn)+dimag(bn1)*dimag(bn)))
            endif
            psi0=psi1
            psi1=psi
            chi0=chi1
            chi1=chi
            xi1=dcmplx(psi1,-chi1)
         enddo
         gsca=dble(2.d0*gsca/qsca)
         qsca=dble((2.d0/(dx*dx))*qsca)
         qext=dble((2.d0/(dx*dx))*qext)
         deallocate(d)
         end subroutine mieext
!
! cartosphere takes the cartesian point (x,y,z) = xp(1), xp(2), xp(3)
! and converts to polar form: r: radius, ct: cos(theta), ep = exp(i phi)
!
         subroutine cartosphere(xp,r,ct,ep)
         implicit none
         real(8) :: xp(3),r,ct
         complex(8) :: ep
         r=xp(1)*xp(1)+xp(2)*xp(2)+xp(3)*xp(3)
         if(r.eq.0.d0) then
            ct=1.d0
            ep=(1.d0,0.d0)
            return
         endif
         r=sqrt(r)
         ct=xp(3)/r
         if(xp(1).eq.0.d0.and.xp(2).eq.0.d0) then
            ep=(1.d0,0.d0)
         else
            ep=dcmplx(xp(1),xp(2))/sqrt(xp(1)*xp(1)+xp(2)*xp(2))
         endif
         return
         end subroutine cartosphere
!
! euler rotation of a point (x,y,z) = xp(1), xp(2), xp(3)
! November 2012
!
         subroutine eulerrotation(xp,eulerangf,dir,xprot)
         implicit none
         integer :: dir
         real(8) :: xp(3),eulerangf(3),eulerang(3),cang(3),sang(3), &
                    mat1(3,3),mat2(3,3),mat3(3,3),xprot(3),xpt(3)
         xpt=xp
         if(dir.eq.1) then
            eulerang=eulerangf
         else
            eulerang(1:3)=-eulerangf(3:1:-1)
         endif
         cang=cos(eulerang)
         sang=sin(eulerang)
         mat1(1,:) = (/cang(1),sang(1),0.d0/)
         mat1(2,:) = (/-sang(1),cang(1),0.d0/)
         mat1(3,:) = (/0.d0,0.d0,1.d0/)
         mat2(1,:) = (/cang(2),0.d0,-sang(2)/)
         mat2(2,:) = (/0.d0,1.d0,0.d0/)
         mat2(3,:) = (/sang(2),0.d0,cang(2)/)
         mat3(1,:) = (/cang(3),sang(3),0.d0/)
         mat3(2,:) = (/-sang(3),cang(3),0.d0/)
         mat3(3,:) = (/0.d0,0.d0,1.d0/)
         xpt=matmul(mat1,xpt)
         xpt=matmul(mat2,xpt)
         xpt=matmul(mat3,xpt)
         xprot=xpt
         end subroutine eulerrotation
!
!
! inverse of a 2 X 2 complex matrix.
! March 2013
!
         subroutine twobytwoinverse(mat,imat)
         implicit none
         integer :: s,t,ss,st
         complex(8) :: mat(2,2),imat(2,2),tmat(2,2),det
         tmat=mat
         det=mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)
         do s=1,2
            ss=(-1)**s
            do t=1,2
               st=(-1)**t
               imat(s,t)=ss*st*tmat(3-t,3-s)/det
            enddo
         enddo
         end subroutine twobytwoinverse

         subroutine effectiverefractiveindex(ndat,edat,d,rieff, &
                    e0)
         use mpidefs
         implicit none
         integer :: ndat,i,rank
         real(8) :: d,xdat(ndat),phase(ndat),amplitude(ndat), &
                    phaseslope,phaseintercept, &
                    ampslope,ampintercept,oldphase, &
                    newphase
         complex(8) :: edat(ndat),rieff,e0
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)

         oldphase=-1.d10
         do i=1,ndat
            xdat(i)=d*dble(i-1)
            phase(i)=datan2(dimag(edat(i)),dble(edat(i)))
            amplitude(i)=dlog(cdabs(edat(i)))
         enddo
         oldphase=phase(1)
         do i=2,ndat
            newphase=phase(i)
            do while(abs(newphase-oldphase).gt.pi)
               if(newphase.gt.oldphase) then
                  newphase=newphase-2*pi
               else
                  newphase=newphase+2*pi
               endif
            enddo
            phase(i)=newphase
            oldphase=newphase
         enddo
         call linearregression(ndat,phase,xdat, &
              phaseslope,phaseintercept)
         call linearregression(ndat,amplitude,xdat, &
              ampslope,ampintercept)
         rieff=dcmplx(phaseslope,-ampslope)
         e0=dexp(ampintercept)*cdexp((0.d0,1.d0)*phaseintercept)
         end subroutine effectiverefractiveindex

         subroutine linearregression(ndat,fdat,xdat,a,b)
         implicit none
         integer :: ndat
         real(8) :: fdat(ndat),xdat(ndat),a,b,xbar,x2bar, &
                    fbar,xfbar
         fbar=sum(fdat(1:ndat))/dble(ndat)
         xbar=sum(xdat(1:ndat))/dble(ndat)
         xfbar=sum(xdat(1:ndat)*fdat(1:ndat))/dble(ndat)
         x2bar=sum(xdat(1:ndat)*xdat(1:ndat))/dble(ndat)
         a=(xfbar-xbar*fbar)/(x2bar-xbar*xbar)
         b=(fbar*x2bar-xbar*xfbar)/(x2bar-xbar*xbar)
         end subroutine linearregression

         subroutine groupfilename(firststring,number,laststring,newstring)
         implicit none
         integer :: number
         character*128 :: firststring,laststring,newstring,sform,intfile
         if(number.lt.10) then
            sform='(a,i1,a,a)'
         elseif(number.lt.100) then
            sform='(a,i2,a,a)'
         else
            sform='(a,i3,a,a)'
         endif
         write(intfile,fmt=sform) trim(firststring),number,'_',trim(laststring)
         read(intfile,'(a)') newstring
         end subroutine groupfilename

         real(8) function circumscribingsphere(shape,aspect_ratio_y, &
             aspect_ratio_z)
         implicit none
         integer :: shape
         real(8) :: ary,arz,chex,c
         real(8), optional :: aspect_ratio_y,aspect_ratio_z
         data chex/0.8660254/
         if(present(aspect_ratio_y)) then
            ary=aspect_ratio_y
         else
            ary=1.d0
         endif
         if(present(aspect_ratio_z)) then
            arz=aspect_ratio_z
         else
            arz=1.d0
         endif
         go to (10,20,30,40) shape
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     prolate spheroid
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
    10   c=1.d0/(ary*arz)**(1.d0/3.d0)
         circumscribingsphere=c*max(1.d0,ary,arz)
         return
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     rectangular solid
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c

    20   c=(pi/(6.d0*ary*arz))**(1.d0/3.d0)
         circumscribingsphere=c*sqrt(1.d0+ary*ary+arz*arz)
         return
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     circular cylinder
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
    30   circumscribingsphere=sqrt(1.d0+arz*arz)
         return
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     hexagon cylinder
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c

    40   circumscribingsphere=sqrt(1.d0+chex*chex+arz*arz)
         return
         end function circumscribingsphere

         logical function insurf(pos,scale,shape,aspect_ratio_y, &
             aspect_ratio_z)
         implicit none
         integer :: shape
         real(8) :: pos(3),scale,ary,arz,x0,y1,z0,z1,dz,xc,zc,chex, &
            x,y,z,r,x1,yedge,c
         real(8), optional :: aspect_ratio_y,aspect_ratio_z
         data x0,x1,y1,z0,z1,dz/-0.57735,1.1547,1.,-0.408248,1.22474,1.63299/
         data xc,zc/-0.422650,-0.591752/
         data chex/0.8660254/

         if(present(aspect_ratio_y)) then
            ary=aspect_ratio_y
         else
            ary=1.d0
         endif
         if(present(aspect_ratio_z)) then
            arz=aspect_ratio_z
         else
            arz=1.d0
         endif
         x=pos(1)
         y=pos(2)
         z=pos(3)

         go to (10,20,30,40) shape
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     prolate spheroid
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
    10   c=1.d0/(ary*arz)**(1.d0/3.d0)
         if (abs(x).gt.c*scale .or. abs(y).gt.c*scale*ary .or. &
              abs(z).gt.c*scale*arz) then
            insurf = .false.
            return
         end if
         r = sqrt(x*x+y*y/ary/ary+z*z/arz/arz)
         if (r.gt.c*scale) then
            insurf = .false.
         else
            insurf = .true.
         end if
         return
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     rectangular solid
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c

    20   c=(pi/(6.d0*ary*arz))**(1.d0/3.d0)
         if (abs(x).gt.c*scale .or. abs(y).gt.c*scale*ary .or. &
              abs(z).gt.c*scale*arz) then
            insurf = .false.
         else
            insurf = .true.
         end if
         return
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     circular cylinder
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
    30   if (abs(x).gt.scale .or. abs(y).gt.scale .or. &
              abs(z).gt.scale*arz) then
            insurf = .false.
            return
         end if
         r = sqrt(x*x+y*y)
         if (r.gt.scale) then
            insurf = .false.
         else
            insurf = .true.
         end if
         return
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c
!   c     hexagon cylinder
!   c
!   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c

    40   if (abs(x).gt.scale .or. abs(y).gt.scale*chex .or. &
              abs(z).gt.scale*arz) then
            insurf = .false.
            return
         end if
         if (abs(x).le..5*scale) then
            insurf = .true.
            return
         end if
         yedge = chex* ((scale-abs(x))/.5)
         if (abs(y).le.yedge) then
            insurf = .true.
         else
            insurf = .false.
         end if
         return
         end function insurf


         subroutine spherical_particle_dipole_map(npart,spos,radius, &
            ndx,ndy,ml,dz,dlat,dipolemap)
         implicit none
         integer :: npart,ndx,ndy,ml,j1,j2,i,nz,iz,nx,ix,ny,iy,ipos(3),ixp,iyp,izp,i1x,i2x,i1y,i2y
         integer(4) :: dipolemap(-ndx/2:(ndx-1)/2,-ndy/2:(ndy-1)/2,-ml/2:(ml-1)/2)
         real(8) :: spos(3,npart),radius,dz,dlat, &
                    wx,wy,thick,rmean(3),rpos(3),rsamp(3),z,xmax,x,rho,ymax,y,r

         wx=dble(ndx)*dlat
         wy=dble(ndy)*dlat
         thick=dble(ml)*dz
         i1x=-ndx/2
         i2x=(ndx-1)/2
         i1y=-ndy/2
         i2y=(ndy-1)/2
         j1=-ml/2
         j2=ml+j1-1
         do i=1,3
            rmean(i)=sum(spos(i,:))/dble(npart)
            spos(i,:)=spos(i,:)-rmean(i)
         enddo
         dipolemap=0

         do i=1,npart
            nz=ceiling(radius/dz)
            rsamp=spos(:,i)
            do iz=-nz-1,nz
               z=dble(iz)*dz+dz/2.d0
               izp=1+floor((z+rsamp(3))/dz)
               xmax=sqrt(radius*radius-z*z)
               nx=ceiling(xmax/dlat)
               do ix=-nx-1,nx
                  x=dble(ix)*dlat+dlat/2.d0
                  rho=sqrt(x*x+z*z)
                  if(rho.gt.radius) cycle
                  ymax=sqrt(radius*radius-rho*rho)
                  ny=ceiling(ymax/dlat)
                  do iy=-ny-1,ny
                     y=dble(iy)*dlat+dlat/2.d0
                     r=sqrt(rho*rho+y*y)
                     if(r.gt.radius) cycle
                     rpos=(/x,y,z/)
                     ipos=1+floor((rpos(:)+rsamp(:))/(/dlat,dlat,dz/))
                     ixp=ipos(1)
                     iyp=ipos(2)
                     izp=ipos(3)
                     if(ixp.lt.i1x) then
                        ixp=ndx+ixp
                     elseif(ixp.gt.i2x) then
                        ixp=ixp-ndx
                     endif
                     if(iyp.lt.i1y) then
                        iyp=ndy+iyp
                     elseif(iyp.gt.i2y) then
                        iyp=iyp-ndy
                     endif
                     if(izp.lt.j1.or.izp.gt.j2) cycle
                     dipolemap(ixp,iyp,izp)=1
                  enddo
               enddo
            enddo
         enddo
         end subroutine spherical_particle_dipole_map

         subroutine ipairswap(n1,n2)
         implicit none
         integer :: n1,n2,n
         n=n1
         n1=n2
         n2=n
         end subroutine ipairswap

         subroutine fpairswap(n1,n2)
         implicit none
         real(8) :: n1,n2,n
         n=n1
         n1=n2
         n2=n
         end subroutine fpairswap

         subroutine lpairswap(n1,n2)
         implicit none
         logical :: n1,n2,n
         n=n1
         n1=n2
         n2=n
         end subroutine lpairswap

         logical function checkn235(n)
         implicit none
         integer :: n,nn,ifac,ll,kk
         nn = n
         ifac = 2
         do ll = 1 , 3
            kk = 0
            do while (mod(nn,ifac).eq.0)
               kk = kk + 1
               nn = nn / ifac
            enddo
            ifac = ifac + ll
         enddo
         checkn235=nn.eq.1
         end function checkn235

         integer function correctn235(n)
         implicit none
         integer :: n,n1
         n1=n
         do while(.not.checkn235(n1))
            n1=n1+1
         enddo
         correctn235=n1
         end function correctn235

         real(8) function mymod(a,b)
         real(8) a,b
         mymod=mod(a,b)
         if (mymod.lt.0.d0) mymod=mymod+b
         end function mymod

         subroutine psdsamp(sigma,maxradius,x)
         implicit none
         integer :: i
         real(8) :: sigma,maxradius,r2pi,f1,fd,x,fmax,s2,xmax, &
                    t1
   !      real(4), external :: ran3
         data r2pi/2.5066282746310002d0/
         if(sigma.eq.0.d0) then
            x=1.d0
            return
         endif
         s2=sigma*sigma
         f1=1.d0
         fd=0.d0
         xmax=exp(-2.5d0*s2)
         t1=(log(xmax)+1.5d0*s2)
         fmax=exp(-t1*t1/(2.d0*s2))/r2pi/xmax/sigma
         i=0
         do while(f1.gt.fd)
            i=i+1
            x=maxradius*ran3(1)
            f1=fmax*ran3(1)
            t1=(log(x)+1.5d0*s2)
            fd=exp(-t1*t1/(2.d0*s2))/r2pi/x/sigma
         enddo
         end subroutine psdsamp

         real(4) function ran3(idum)
         integer :: idum
         integer :: mbig,mseed,mz
         real(4) ::  fac
         parameter(mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
         integer :: i,ii,k
         integer :: mj,mk
         integer, save:: iff,inext,inextp,ma(55)
         data iff/0/
         if(idum.lt.0.or.iff.eq.0) then
            iff=1
            mj=mseed-iabs(idum)
            mj=mod(mj,mbig)
            ma(55)=mj
            mk=1
            do i=1,54
               ii=mod(21*i,55)
               ma(ii)=mk
               mk=mj-mk
               if(mk.lt.mz) mk=mk+mbig
               mj=ma(ii)
            enddo
            do k=1,4
               do i=1,55
                  ma(i)=ma(i)-ma(1+mod(i+30,55))
                  if(ma(i).lt.mz) ma(i)=ma(i)+mbig
               enddo
            enddo
            inext=0
            inextp=31
            idum=1
         endif
20       inext=inext+1
         if(inext.eq.56) inext=1
         inextp=inextp+1
         if(inextp.eq.56) inextp=1
         mj=ma(inext)-ma(inextp)
         if(mj.lt.mz) mj=mj+mbig
         ma(inext)=mj
         ran3=mj*fac
         if(ran3.gt.1.e-9.and.ran3.lt.0.99999999) return
         goto 20
         end function ran3

         complex(8) function bruggeman_mixing_rule(fv,ri1,ri0)
         implicit none
         real(8) :: fv
         complex(8) :: ri1,ri0,e1,e0,ea
         e0=ri0*ri0
         e1=ri1*ri1
         ea=.25d0*(2.d0*e0-e1-3.d0*fv*(e0-e1) &
           +cdsqrt(8.d0*e0*e1+(e1*(1.d0-3.d0*fv)+e0*(3.d0*fv-2.d0))**2.))
         bruggeman_mixing_rule=cdsqrt(ea)
         end function bruggeman_mixing_rule

         complex(8) function bruggeman_mixing_rule_multicomp(ncomp,fv,ri)
         implicit none
         integer :: ncomp,n,imax,itemp(1)
         real(8) :: fv(ncomp)
         complex(8) :: ri(ncomp),e(ncomp),ea,et,eaold
         itemp=maxloc(fv)
         imax=itemp(1)
         e=ri*ri
         ea=sum(fv(1:ncomp)*e(1:ncomp))
         eaold=ea
         do
            et=0.d0
            do n=1,ncomp
               if(n.ne.imax) then
                  et=et-(e(n)-ea)/(e(n)+2.d0*ea)*fv(n)
               endif
            enddo
            ea=e(imax)*(fv(imax)-et)/(fv(imax)+2.d0*et)
            if(cdabs(ea-eaold).lt.1.d-5) exit
            eaold=ea
         enddo
         bruggeman_mixing_rule_multicomp=cdsqrt(ea)
         end function bruggeman_mixing_rule_multicomp

      end module specialfuncs

!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!
!  Multidimensional (1D-6D) B-spline interpolation of data on a regular grid.
!  Basic pure subroutine interface.
!
!### Notes
!
!  This module is based on the B-spline and spline routines from [1].
!  The original Fortran 77 routines were converted to free-form source.
!  Some of them are relatively unchanged from the originals, but some have
!  been extensively refactored. In addition, new routines for
!  1d, 4d, 5d, and 6d interpolation were also created (these are simply
!  extensions of the same algorithm into higher dimensions).
!
!### See also
!  * An object-oriented interface can be found in [[bspline_oo_module]].
!
!### References
!
!  1. DBSPLIN and DTENSBS from the
!     [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).
!     Original code is public domain.
!  2. Carl de Boor, "A Practical Guide to Splines",
!     Springer-Verlag, New York, 1978.
!  3. Carl de Boor, [Efficient Computer Manipulation of Tensor
!     Products](http://dl.acm.org/citation.cfm?id=355831),
!     ACM Transactions on Mathematical Software,
!     Vol. 5 (1979), p. 173-182.
!  4. D.E. Amos, "Computation with Splines and B-Splines",
!     SAND78-1968, Sandia Laboratories, March, 1979.
!  5. Carl de Boor,
!     [Package for calculating with B-splines](http://epubs.siam.org/doi/abs/10.1137/0714026),
!     SIAM Journal on Numerical Analysis 14, 3 (June 1977), p. 441-472.

    module bspline_sub_module

    implicit none

    private

    !Spline function order (order = polynomial degree + 1)
    integer,parameter,public :: bspline_order_quadratic = 3
    integer,parameter,public :: bspline_order_cubic     = 4
    integer,parameter,public :: bspline_order_quartic   = 5
    integer,parameter,public :: bspline_order_quintic   = 6

      integer :: order
      integer :: inbvx_dat(32),inbvy_dat(32),iloy_dat(32),nx_dat,ny_dat,n_dat,nz_dat
      real(8) :: delkx,delky
      real(8), allocatable :: kx_dat(:),ky_dat(:),tx_dat(:,:),ty_dat(:,:), &
         bcoef_dat(:,:,:)


    !main routines:
    public :: db2ink, db2val

    public :: smatsplinecalc,smatsplineeval
    public :: get_status_message

    contains


      subroutine smatsplinecalc(nx1,nx2,ny1,ny2,smat, &
         k0x,k0y,wx,wy,spline_order,auto_zero_s12)
      implicit none
      logical :: autozero
      logical, optional :: auto_zero_s12
      integer :: nx1,nx2,ny1,ny2,sx,sy,s,i,iflag,px,py
      integer, optional :: spline_order
      real(8) :: wx,wy,k0x,k0y,pi,smat(32,nx1:nx2,ny1:ny2), &
         tdat(nx1:nx2,ny1:ny2,32), &
         cp,sp,kx,ky,kr,ct
      if(present(spline_order)) then
         order=spline_order
      else
         order=4
      endif
      if(present(auto_zero_s12)) then
         autozero=auto_zero_s12
      else
         autozero=.false.
      endif
      pi=4.d0*datan(1.d0)
      delkx=2.d0*pi/wx
      delky=2.d0*pi/wy
      nx_dat=nx2-nx1+1
      ny_dat=ny2-ny1+1
      if(allocated(kx_dat)) then
         deallocate(kx_dat,ky_dat,bcoef_dat,tx_dat,ty_dat)
      endif
      allocate(kx_dat(nx_dat),ky_dat(ny_dat), &
         bcoef_dat(nx_dat,ny_dat,32), &
         tx_dat(nx_dat+order,32), &
         ty_dat(ny_dat+order,32))
      do i=1,32
         inbvx_dat(i)=1
         inbvy_dat(i)=1
         iloy_dat(i)=1
      enddo
      s=0
      do sx=nx1,nx2
         s=s+1
         kx_dat(s)=k0x+delkx*dble(sx)
      enddo
      s=0
      do sy=ny1,ny2
         s=s+1
         ky_dat(s)=k0y+delky*dble(sy)
      enddo
      do sy=ny1,ny2
         ky=ky_dat(sy-ny1+1)
         do sx=nx1,nx2
            kx=kx_dat(sx-nx1+1)
            kr=sqrt(kx*kx+ky*ky)
            ct=sqrt(1.d0-min(1.d0,kr*kr))
            if(kr.eq.0.d0) then
               cp=1.d0
               sp=0.d0
            else
               cp=kx/kr
               sp=ky/kr
            endif
            if(kr.ge.1.d0) then
               px=int((cp-k0x)/delkx)
               py=int((sp-k0y)/delky)
            else
               px=sx
               py=sy
            endif
!if(smat(1,px,py).eq.0.d0) write(*,'(4i5,e13.5)') px,py,sx,sy,kr
            tdat(sx,sy,1:32)=smat(1:32,px,py)*ct
         enddo
      enddo
      tdat(0,0,17:32)=.25*(tdat(-1,0,17:32)+tdat(1,0,17:32)+tdat(0,-1,17:32)+tdat(0,1,17:32))
      if(autozero) then
         tdat(0,0,2)=0.d0
         tdat(0,0,5)=0.d0
         tdat(0,0,12)=0.d0
         tdat(0,0,15)=0.d0
         tdat(0,0,18)=0.d0
         tdat(0,0,21)=0.d0
         tdat(0,0,28)=0.d0
         tdat(0,0,31)=0.d0
      endif

      do i=1,32
         call db2ink(kx_dat,nx_dat,ky_dat,ny_dat, &
            tdat(nx1:nx2,ny1:ny2,i), &
            order,order,0,tx_dat(:,i), &
            ty_dat(:,i),bcoef_dat(:,:,i),iflag)
      enddo
      end subroutine smatsplinecalc

      subroutine smatsplineeval(kx,ky,smat,i_start,i_end)
      implicit none
      integer :: i,iflag,istart,iend,i1
      integer, optional :: i_start,i_end
      real(8) :: kx,ky,smat(*),kxt,kyt,kr,ct
      if(present(i_start)) then
         istart=i_start
      else
         istart=1
      endif
      if(present(i_end)) then
         iend=i_end
      else
         iend=32
      endif
      kxt=max(kx,kx_dat(1))
      kxt=min(kx,kx_dat(nx_dat))
      kyt=max(ky,ky_dat(1))
      kyt=min(ky,ky_dat(ny_dat))
      kr=sqrt(kxt*kxt+kyt*kyt)
      ct=sqrt(1.d0-min(1.d0,kr*kr))
      if(kr.gt.1.d0) then
         kxt=(1.d0)/kr*kxt
         kyt=(1.d0)/kr*kyt
      endif

      i1=0
      do i=istart,iend
         i1=i1+1
         call db2val(kxt,kyt,0,0,tx_dat(:,i),ty_dat(:,i),nx_dat,ny_dat, &
            order,order,bcoef_dat(:,:,i),smat(i1),iflag, &
            inbvx_dat(i),inbvy_dat(i),iloy_dat(i),.true.)
            smat(i1)=smat(i1)/max(sqrt(delkx*delky)/2.d0,ct)
!         if(i.eq.1) then
!            smat(i1)=(10.d0**smat(i1))
!         endif
!         if(i.eq.17) then
!            smat(i1)=(10.d0**smat(i1))
!         endif
      enddo

      end subroutine smatsplineeval


!*****************************************************************************************
!>
!  Determines the parameters of a function that interpolates
!  the two-dimensional gridded data
!  $$ [x(i),y(j),\mathrm{fcn}(i,j)] ~\mathrm{for}~ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y $$
!  The interpolating function and its derivatives may
!  subsequently be evaluated by the function [[db2val]].
!
!  The interpolating function is a piecewise polynomial function
!  represented as a tensor product of one-dimensional b-splines. the
!  form of this function is
!
!  $$ s(x,y) = \sum_{i=1}^{n_x} \sum_{j=1}^{n_y} a_{ij} u_i(x) v_j(y) $$
!
!  where the functions \(u_i\) and \(v_j\) are one-dimensional b-spline
!  basis functions. the coefficients \( a_{ij} \) are chosen so that
!
!  $$ s(x(i),y(j)) = \mathrm{fcn}(i,j) ~\mathrm{for}~ i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y $$
!
!  Note that for each fixed value of \(y\), \( s(x,y) \) is a piecewise
!  polynomial function of \(x\) alone, and for each fixed value of \(x\), \( s(x,y) \)
!  is a piecewise polynomial function of \(y\) alone. in one dimension
!  a piecewise polynomial may be created by partitioning a given
!  interval into subintervals and defining a distinct polynomial piece
!  on each one. the points where adjacent subintervals meet are called
!  knots. each of the functions \(u_i\) and \(v_j\) above is a piecewise
!  polynomial.
!
!  Users of [[db2ink]] choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the \(x\) and
!  \(y\) directions (`kx` and `ky`). users also may define their own knot
!  sequence in \(x\) and \(y\) separately (`tx` and `ty`). if `iflag=0`, however,
!  [[db2ink]] will choose sequences of knots that result in a piecewise
!  polynomial interpolant with `kx-2` continuous partial derivatives in
!  \(x\) and `ky-2` continuous partial derivatives in \(y\). (`kx` knots are taken
!  near each endpoint in the \(x\) direction, not-a-knot end conditions
!  are used, and the remaining knots are placed at data points if `kx`
!  is even or at midpoints between data points if `kx` is odd. the \(y\)
!  direction is treated similarly.)
!
!  After a call to [[db2ink]], all information necessary to define the
!  interpolating function are contained in the parameters `nx`, `ny`, `kx`,
!  `ky`, `tx`, `ty`, and `bcoef`. These quantities should not be altered until
!  after the last call of the evaluation routine [[db2val]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db2ink(x,nx,y,ny,fcn,kx,ky,iknot,tx,ty,bcoef,iflag)

    implicit none

    integer,intent(in)                      :: nx     !! Number of \(x\) abcissae
    integer,intent(in)                      :: ny     !! Number of \(y\) abcissae
    integer,intent(in)                      :: kx     !! The order of spline pieces in \(x\)
                                                      !! ( \( 2 \le k_x < n_x \) )
                                                      !! (order = polynomial degree + 1)
    integer,intent(in)                      :: ky     !! The order of spline pieces in \(y\)
                                                      !! ( \( 2 \le k_y < n_y \) )
                                                      !! (order = polynomial degree + 1)
    real(8),dimension(:),intent(in)        :: x      !! `(nx)` array of \(x\) abcissae. Must be strictly increasing.
    real(8),dimension(:),intent(in)        :: y      !! `(ny)` array of \(y\) abcissae. Must be strictly increasing.
    real(8),dimension(:,:),intent(in)      :: fcn    !! `(nx,ny)` matrix of function values to interpolate.
                                                      !! `fcn(i,j)` should contain the function value at the
                                                      !! point (`x(i)`,`y(j)`)
    integer,intent(in)                      :: iknot  !! knot sequence flag:
                                                      !!
                                                      !! * 0 = knot sequence chosen by [[db1ink]].
                                                      !! * 1 = knot sequence chosen by user.
    real(8),dimension(:),intent(inout)     :: tx     !! The `(nx+kx)` knots in the \(x\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(8),dimension(:),intent(inout)     :: ty     !! The `(ny+ky)` knots in the \(y\) direction for the spline interpolant.
                                                      !!
                                                      !! * If `iknot=0` these are chosen by [[db2ink]].
                                                      !! * If `iknot=1` these are specified by the user.
                                                      !!
                                                      !! Must be non-decreasing.
    real(8),dimension(:,:),intent(out)     :: bcoef  !! `(nx,ny)` matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                     :: iflag  !! *  0 = successful execution.
                                                      !! *  2 = `iknot` out of range.
                                                      !! *  3 = `nx` out of range.
                                                      !! *  4 = `kx` out of range.
                                                      !! *  5 = `x` not strictly increasing.
                                                      !! *  6 = `tx` not non-decreasing.
                                                      !! *  7 = `ny` out of range.
                                                      !! *  8 = `ky` out of range.
                                                      !! *  9 = `y` not strictly increasing.
                                                      !! * 10 = `ty` not non-decreasing.
                                                      !! * 700 = `size(x)`  \( \ne \) `size(fcn,1)`
                                                      !! * 701 = `size(y)`  \( \ne \) `size(fcn,2)`
                                                      !! * 706 = `size(x)`  \( \ne \) `nx`
                                                      !! * 707 = `size(y)`  \( \ne \) `ny`
                                                      !! * 712 = `size(tx)` \( \ne \) `nx+kx`
                                                      !! * 713 = `size(ty)` \( \ne \) `ny+ky`
                                                      !! * 800 = `size(x)`  \( \ne \) `size(bcoef,1)`
                                                      !! * 801 = `size(y)`  \( \ne \) `size(bcoef,2)`

    real(8),dimension(nx*ny) :: temp
    real(8),dimension(max(2*kx*(nx+1),2*ky*(ny+1))) :: work
    logical :: status_ok

    !check validity of inputs

    call check_inputs('db2ink',&
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,&
                        kx=kx,ky=ky,&
                        x=x,y=y,&
                        tx=tx,ty=ty,&
                        f2=fcn,&
                        bcoef2=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        !choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
        end if

        !construct b-spline coefficients

                      call dbtpcf(x,nx,fcn, nx,ny,tx,kx,temp, work,iflag)
        if (iflag==0) call dbtpcf(y,ny,temp,ny,nx,ty,ky,bcoef,work,iflag)

    end if

    end subroutine db2ink
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db2ink]] or one of its
!  derivatives at the point (`xval`,`yval`).
!
!  To evaluate the interpolant
!  itself, set `idx=idy=0`, to evaluate the first partial with respect
!  to `x`, set `idx=1,idy=0`, and so on.
!
!  [[db2val]] returns 0.0 if `(xval,yval)` is out of range. that is, if
!```fortran
!   xval < tx(1) .or. xval > tx(nx+kx) .or.
!   yval < ty(1) .or. yval > ty(ny+ky)
!```
!  if the knots tx and ty were chosen by [[db2ink]], then this is equivalent to:
!```fortran
!   xval < x(1) .or. xval > x(nx)+epsx .or.
!   yval < y(1) .or. yval > y(ny)+epsy
!```
!  where
!```fortran
!   epsx = 0.1*(x(nx)-x(nx-1))
!   epsy = 0.1*(y(ny)-y(ny-1))
!```
!
!  The input quantities `tx`, `ty`, `nx`, `ny`, `kx`, `ky`, and `bcoef` should be
!  unchanged since the last call of [[db2ink]].
!
!### History
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag,inbvx,inbvy,iloy,extrap)

    implicit none

    integer,intent(in)                   :: idx      !! \(x\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: idy      !! \(y\) derivative of piecewise polynomial to evaluate.
    integer,intent(in)                   :: nx       !! the number of interpolation points in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer,intent(in)                   :: ny       !! the number of interpolation points in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    integer,intent(in)                   :: kx       !! order of polynomial pieces in \(x\).
                                                     !! (same as in last call to [[db2ink]])
    integer,intent(in)                   :: ky       !! order of polynomial pieces in \(y\).
                                                     !! (same as in last call to [[db2ink]])
    real(8),intent(in)                  :: xval     !! \(x\) coordinate of evaluation point.
    real(8),intent(in)                  :: yval     !! \(y\) coordinate of evaluation point.
    real(8),dimension(nx+kx),intent(in) :: tx       !! sequence of knots defining the piecewise polynomial
                                                     !! in the \(x\) direction.
                                                     !! (same as in last call to [[db2ink]])
    real(8),dimension(ny+ky),intent(in) :: ty       !! sequence of knots defining the piecewise
                                                     !! polynomial in the \(y\) direction.
                                                     !! (same as in last call to [[db2ink]])
    real(8),dimension(nx,ny),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db2ink]].
    real(8),intent(out)                 :: f        !! interpolated value
    integer,intent(out)                  :: iflag    !! status flag:
                                                     !!
                                                     !! * \( = 0 \)   : no errors
                                                     !! * \( \ne 0 \) : error
    integer,intent(inout)                :: inbvx    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer,intent(inout)                :: inbvy    !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    integer,intent(inout)                :: iloy     !! initialization parameter which must be set to 1
                                                     !! the first time this routine is called,
                                                     !! and must not be changed by the user.
    logical,intent(in),optional          :: extrap   !! if extrapolation is allowed
                                                     !! (if not present, default is False)

    integer :: k, lefty, mflag, kcol
    real(8),dimension(ky) :: temp
    real(8),dimension(3*max(kx,ky)) :: work

    f = 0.0d0

    iflag = check_value(xval,tx,1,extrap); if (iflag/=0) return
    iflag = check_value(yval,ty,2,extrap); if (iflag/=0) return

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag,extrap); if (mflag /= 0) return

    kcol = lefty - ky
    do k=1,ky
        kcol = kcol + 1
        call dbvalu(tx,bcoef(:,kcol),nx,kx,idx,xval,inbvx,work,iflag,temp(k),extrap)
        if (iflag/=0) return !error
    end do

    kcol = lefty - ky + 1
    call dbvalu(ty(kcol:),temp,ky,ky,idy,yval,inbvy,work,iflag,f,extrap)

    end subroutine db2val
!*****************************************************************************************

!*****************************************************************************************
!>
!  Checks if the value is withing the range of the knot vectors.
!  This is called by the various `db*val` routines.

    pure function check_value(x,t,i,extrap) result(iflag)

    implicit none

    integer :: iflag  !! returns 0 if value is OK, otherwise returns `600+i`
    real(8),intent(in) :: x !! the value to check
    integer,intent(in) :: i !! 1=x, 2=y, 3=z, 4=q, 5=r, 6=s
    real(8),dimension(:),intent(in) :: t  !! the knot vector
    logical,intent(in),optional :: extrap  !! if extrapolation is allowed
                                           !! (if not present, default is False)

    logical :: allow_extrapolation  !! if extrapolation is allowed

    if (present(extrap)) then
        allow_extrapolation = extrap
    else
        allow_extrapolation = .false.
    end if

    if (allow_extrapolation) then
        ! in this case all values are OK
        iflag = 0
    else
        if (x<t(1) .or. x>t(size(t))) then
            iflag = 600 + i  ! value out of bounds (601, 602, etc.)
        else
            iflag = 0
        end if
    end if

    end function check_value
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the validity of the inputs to the `db*ink` routines.
!  Prints warning message if there is an error,
!  and also sets iflag and status_ok.
!
!  Supports up to 6D: `x`,`y`,`z`,`q`,`r`,`s`
!
!### Notes
!
!  The code is new, but the logic is based on the original
!  logic in the CMLIB routines `db2ink` and `db3ink`.
!
!### History
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine check_inputs(routine,&
                            iknot,&
                            iflag,&
                            nx,ny,nz,nq,nr,ns,&
                            kx,ky,kz,kq,kr,ks,&
                            x,y,z,q,r,s,&
                            tx,ty,tz,tq,tr,ts,&
                            f1,f2,f3,f4,f5,f6,&
                            bcoef1,bcoef2,bcoef3,bcoef4,bcoef5,bcoef6,&
                            status_ok)

    implicit none

    character(len=*),intent(in)                         :: routine
    integer,intent(in)                                  :: iknot !! = 0 if the `INK` routine is computing the knots.
    integer,intent(out)                                 :: iflag
    integer,intent(in),optional                         :: nx,ny,nz,nq,nr,ns
    integer,intent(in),optional                         :: kx,ky,kz,kq,kr,ks
    real(8),dimension(:),intent(in),optional           :: x,y,z,q,r,s
    real(8),dimension(:),intent(in),optional           :: tx,ty,tz,tq,tr,ts
    real(8),dimension(:),intent(in),optional           :: f1,bcoef1
    real(8),dimension(:,:),intent(in),optional         :: f2,bcoef2
    real(8),dimension(:,:,:),intent(in),optional       :: f3,bcoef3
    real(8),dimension(:,:,:,:),intent(in),optional     :: f4,bcoef4
    real(8),dimension(:,:,:,:,:),intent(in),optional   :: f5,bcoef5
    real(8),dimension(:,:,:,:,:,:),intent(in),optional :: f6,bcoef6
    logical,intent(out)                                 :: status_ok

    logical :: error

    status_ok = .false.

    if ((iknot < 0) .or. (iknot > 1)) then

        !write(error_unit,'(A,1X,I5)') &
        !    trim(routine)//' - iknot is out of range: ',iflag
        iflag = 2

    else

        call check('x',nx,kx,x,tx,[3,  4, 5, 6,706,712],iflag,error); if (error) return
        call check('y',ny,ky,y,ty,[7,  8, 9,10,707,713],iflag,error); if (error) return
        call check('z',nz,kz,z,tz,[11,12,13,14,708,714],iflag,error); if (error) return
        call check('q',nq,kq,q,tq,[15,16,17,18,709,715],iflag,error); if (error) return
        call check('r',nr,kr,r,tr,[19,20,21,22,710,716],iflag,error); if (error) return
        call check('s',ns,ks,s,ts,[23,24,25,26,711,717],iflag,error); if (error) return

        if (present(x) .and. present(f1) .and. present(bcoef1)) then
            if (size(x)/=size(f1,1))     then; iflag = 700; return; end if
            if (size(x)/=size(bcoef1,1)) then; iflag = 800; return; end if
        end if
        if (present(x) .and. present(y) .and. present(f2) .and. present(bcoef2)) then
            if (size(x)/=size(f2,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f2,2))     then; iflag = 701; return; end if
            if (size(x)/=size(bcoef2,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef2,2)) then; iflag = 801; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(f3) .and. &
            present(bcoef3)) then
            if (size(x)/=size(f3,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f3,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f3,3))     then; iflag = 702; return; end if
            if (size(x)/=size(bcoef3,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef3,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef3,3)) then; iflag = 802; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(f4) .and. present(bcoef4)) then
            if (size(x)/=size(f4,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f4,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f4,3))     then; iflag = 702; return; end if
            if (size(q)/=size(f4,4))     then; iflag = 703; return; end if
            if (size(x)/=size(bcoef4,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef4,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef4,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef4,4)) then; iflag = 803; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(f5) .and. present(bcoef5)) then
            if (size(x)/=size(f5,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f5,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f5,3))     then; iflag = 702; return; end if
            if (size(q)/=size(f5,4))     then; iflag = 703; return; end if
            if (size(r)/=size(f5,5))     then; iflag = 704; return; end if
            if (size(x)/=size(bcoef5,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef5,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef5,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef5,4)) then; iflag = 803; return; end if
            if (size(r)/=size(bcoef5,5)) then; iflag = 804; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(s) .and. present(f6) .and. present(bcoef6)) then
            if (size(x)/=size(f6,1))     then; iflag = 700; return; end if
            if (size(y)/=size(f6,2))     then; iflag = 701; return; end if
            if (size(z)/=size(f6,3))     then; iflag = 702; return; end if
            if (size(q)/=size(f6,4))     then; iflag = 703; return; end if
            if (size(r)/=size(f6,5))     then; iflag = 704; return; end if
            if (size(s)/=size(f6,6))     then; iflag = 705; return; end if
            if (size(x)/=size(bcoef6,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef6,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef6,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef6,4)) then; iflag = 803; return; end if
            if (size(r)/=size(bcoef6,5)) then; iflag = 804; return; end if
            if (size(s)/=size(bcoef6,6)) then; iflag = 805; return; end if
        end if

        status_ok = .true.
        iflag = 0

    end if

    contains

        pure subroutine check(s,n,k,x,t,ierrs,iflag,error)  !! check `t`,`x`,`n`,`k` for validity

        implicit none

        character(len=1),intent(in)               :: s     !! coordinate string: 'x','y','z','q','r','s'
        integer,intent(in)              ,optional :: n     !! size of `x`
        integer,intent(in)              ,optional :: k     !! order
        real(8),dimension(:),intent(in),optional :: x     !! abcissae vector
        real(8),dimension(:),intent(in),optional :: t     !! knot vector `size(n+k)`
        integer,dimension(:),intent(in)           :: ierrs !! int error codes for `n`,`k`,`x`,`t`,
                                                           !! `size(x)`,`size(t)` checks
        integer,intent(out)                       :: iflag !! status return code
        logical,intent(out)                       :: error !! true if there was an error

        if (present(n) .and. present(k) .and. present(x) .and. present(t)) then
            call check_n('n'//s,n,x,[ierrs(1),ierrs(5)],iflag,error); if (error) return
            call check_k('k'//s,k,n,ierrs(2),iflag,error); if (error) return
            call check_x(s,n,x,ierrs(3),iflag,error); if (error) return
            if (iknot /= 0) then
                call check_t('t'//s,n,k,t,[ierrs(4),ierrs(6)],iflag,error); if (error) return
            end if
        end if

        end subroutine check

        pure subroutine check_n(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)     :: s
        integer,intent(in)              :: n
        real(8),dimension(:),intent(in):: x     !! abcissae vector
        integer,dimension(2),intent(in) :: ierr  !! [n<3 check, size(x)==n check]
        integer,intent(out)             :: iflag !! status return code
        logical,intent(out)             :: error

        if (n < 3) then
            !write(error_unit,'(A,1X,I5)') &
            !    trim(routine)//' - '//trim(s)//' is out of range: ',n
            iflag = ierr(1)
            error = .true.
        else
            if (size(x)/=n) then
                !write(error_unit,'(A,1X,I5)') &
                !    trim(routine)//' - '//trim(s)//' is not abscissa vector size'
                iflag = ierr(2)
                error = .true.
            else
                error = .false.
            end if
        end if

        end subroutine check_n

        pure subroutine check_k(s,k,n,ierr,iflag,error)

        implicit none

        character(len=*),intent(in) :: s
        integer,intent(in)          :: k
        integer,intent(in)          :: n
        integer,intent(in)          :: ierr
        integer,intent(out)         :: iflag !! status return code
        logical,intent(out)         :: error

        if ((k < 2) .or. (k >= n)) then
            !write(error_unit,'(A,1X,I5)') &
            !    trim(routine)//' - '//trim(s)//' is out of range: ',k
            iflag = ierr
            error = .true.
        else
            error = .false.
        end if

        end subroutine check_k

        pure subroutine check_x(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        real(8),dimension(:),intent(in)  :: x
        integer,intent(in)                :: ierr
        integer,intent(out)               :: iflag !! status return code
        logical,intent(out)               :: error

        integer :: i

        error = .true.
        do i=2,n
            if (x(i) <= x(i-1)) then
                iflag = ierr
                !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
                !            ' array must be strictly increasing'
                return
            end if
        end do
        error = .false.

        end subroutine check_x

        pure subroutine check_t(s,n,k,t,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        integer,intent(in)                :: k
        real(8),dimension(:),intent(in)  :: t
        integer,dimension(2),intent(in)   :: ierr  !! [non-decreasing check, size check]
        integer,intent(out)               :: iflag !! status return code
        logical,intent(out)               :: error

        integer :: i

        error = .true.

        if (size(t)/=(n+k)) then
            !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
            !            ' array is not the correct size'
            iflag = ierr(2)
            return
        end if

        do i=2,n + k
            if (t(i) < t(i-1))  then
                iflag = ierr(1)
                !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
                !            ' array must be non-decreasing'
                return
            end if
        end do
        error = .false.

        end subroutine check_t

    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbknot chooses a knot sequence for interpolation of order k at the
!  data points x(i), i=1,..,n.  the n+k knots are placed in the array
!  t.  k knots are placed at each endpoint and not-a-knot end
!  conditions are used.  the remaining knots are placed at data points
!  if n is even and between data points if n is odd.  the rightmost
!  knot is shifted slightly to the right to insure proper interpolation
!  at x(n) (see page 350 of the reference).
!
!### History
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbknot(x,n,k,t)

    implicit none

    integer,intent(in)                 :: n
    integer,intent(in)                 :: k
    real(8),dimension(n),intent(in)   :: x
    real(8),dimension(:),intent(out)  :: t

    integer  :: i, j, ipj, npj, ip1, jstrt
    real(8) :: rnot

    !put k knots at each endpoint
    !(shift right endpoints slightly -- see pg 350 of reference)
    rnot = x(n) + 0.1d0*( x(n)-x(n-1) )
    do j=1,k
        t(j)   = x(1)
        npj    = n + j
        t(npj) = rnot
    end do

    !distribute remaining knots

    if (mod(k,2) == 1)  then

        !case of odd k --  knots between data points

        i = (k-1)/2 - k
        ip1 = i + 1
        jstrt = k + 1
        do j=jstrt,n
            ipj = i + j
            t(j) = 0.5d0*( x(ipj) + x(ipj+1) )
        end do

    else

        !case of even k --  knots at data points

        i = (k/2) - k
        jstrt = k+1
        do j=jstrt,n
            ipj = i + j
            t(j) = x(ipj)
        end do

    end if

    end subroutine dbknot
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbtpcf computes b-spline interpolation coefficients for nf sets
!  of data stored in the columns of the array fcn. the b-spline
!  coefficients are stored in the rows of bcoef however.
!  each interpolation is based on the n abcissa stored in the
!  array x, and the n+k knots stored in the array t. the order
!  of each interpolation is k.
!
!### History
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag)

    integer,intent(in)                    :: n
    integer,intent(in)                    :: nf
    integer,intent(in)                    :: ldf
    integer,intent(in)                    :: k
    real(8),dimension(n),intent(in)      :: x
    real(8),dimension(ldf,nf),intent(in) :: fcn
    real(8),dimension(*),intent(in)      :: t
    real(8),dimension(nf,n),intent(out)  :: bcoef
    real(8),dimension(*),intent(out)     :: work   !! work array of size >= `2*k*(n+1)`
    integer,intent(out)                   :: iflag  !!   0: no errors
                                                    !! 301: n should be >0

    integer :: i, j, m1, m2, iq, iw

    ! check for null input

    if (nf > 0)  then

        ! partition work array
        m1 = k - 1
        m2 = m1 + k
        iq = 1 + n
        iw = iq + m2*n+1

        ! compute b-spline coefficients

        ! first data set

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw),iflag)
        if (iflag == 0) then
            do i=1,n
                bcoef(1,i) = work(i)
            end do

            !  all remaining data sets by back-substitution

            if (nf == 1)  return
            do j=2,nf
                do i=1,n
                    work(i) = fcn(i,j)
                end do
                call dbnslv(work(iq),m2,n,m1,m1,work)
                do i=1,n
                    bcoef(j,i) = work(i)
                end do
            end do
        end if

    else
        !write(error_unit,'(A)') 'dbtpcf - n should be >0'
        iflag = 301
    end if

    end subroutine dbtpcf
!*****************************************************************************************

!*****************************************************************************************
!>
!  dbintk produces the b-spline coefficients, bcoef, of the
!  b-spline of order k with knots t(i), i=1,...,n+k, which
!  takes on the value y(i) at x(i), i=1,...,n.  the spline or
!  any of its derivatives can be evaluated by calls to [[dbvalu]].
!
!  the i-th equation of the linear system a*bcoef = b for the
!  coefficients of the interpolant enforces interpolation at
!  x(i), i=1,...,n.  hence, b(i) = y(i), for all i, and a is
!  a band matrix with 2k-1 bands if a is invertible.  the matrix
!  a is generated row by row and stored, diagonal by diagonal,
!  in the rows of q, with the main diagonal going into row k.
!  the banded system is then solved by a call to dbnfac (which
!  constructs the triangular factorization for a and stores it
!  again in q), followed by a call to dbnslv (which then
!  obtains the solution bcoef by substitution).  dbnfac does no
!  pivoting, since the total positivity of the matrix a makes
!  this unnecessary.  the linear system to be solved is
!  (theoretically) invertible if and only if
!          t(i) < x(i) < t(i+k),        for all i.
!  equality is permitted on the left for i=1 and on the right
!  for i=n when k knots are used at x(1) or x(n).  otherwise,
!  violation of this condition is certain to lead to an error.
!
!# Error conditions
!
!  * improper input
!  * singular system of equations
!
!### History
!  * splint written by carl de boor [5]
!  * dbintk author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations. (jec)
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbintk(x,y,t,n,k,bcoef,q,work,iflag)

    implicit none

    integer,intent(in)                :: n      !! number of data points, n >= k
    real(8),dimension(n),intent(in)  :: x      !! vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real(8),dimension(n),intent(in)  :: y      !! corresponding vector of length n containing data
                                                !! point ordinates.
    real(8),dimension(*),intent(in)  :: t      !! knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer,intent(in)                :: k      !! order of the spline, k >= 1
    real(8),dimension(n),intent(out) :: bcoef  !! a vector of length n containing the b-spline coefficients
    real(8),dimension(*),intent(out) :: q      !! a work vector of length (2*k-1)*n, containing
                                                !! the triangular factorization of the coefficient
                                                !! matrix of the linear system being solved.  the
                                                !! coefficients for the interpolant of an
                                                !! additional data set (x(i),yy(i)), i=1,...,n
                                                !! with the same abscissa can be obtained by loading
                                                !! yy into bcoef and then executing
                                                !! call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
    real(8),dimension(*),intent(out) :: work   !! work vector of length 2*k
    integer,intent(out)               :: iflag  !! *   0: no errors.
                                                !! * 100: k does not satisfy k>=1.
                                                !! * 101: n does not satisfy n>=k.
                                                !! * 102: x(i) does not satisfy x(i)<x(i+1) for some i.
                                                !! * 103: some abscissa was not in the support of the.
                                                !! corresponding basis function and the system is singular.
                                                !! * 104: the system of solver detects a singular system.
                                                !! although the theoretical conditions for a solution were satisfied.

    integer :: iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real(8) :: xi
    logical :: found

    if (k<1) then
        !write(error_unit,'(A)') 'dbintk - k does not satisfy k>=1'
        iflag = 100
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbintk - n does not satisfy n>=k'
        iflag = 101
        return
    end if

    jj = n - 1
    if (jj/=0) then
        do i=1,jj
            if (x(i)>=x(i+1)) then
                !write(error_unit,'(A)') 'dbintk - x(i) does not satisfy x(i)<x(i+1) for some i'
                iflag = 102
                return
            end if
        end do
    end if

    np1 = n + 1
    km1 = k - 1
    kpkm2 = 2*km1
    left = k
    ! zero out all entries of q
    lenq = n*(k+km1)
    do i=1,lenq
        q(i) = 0.0d0
    end do

    ! loop over i to construct the n interpolation equations
    do i=1,n

        xi = x(i)
        ilp1mx = min(i+k,np1)
        ! find left in the closed interval (i,i+k-1) such that
        !         t(left) <= x(i) < t(left+1)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
            !             ' corresponding basis function and the system is singular'
            iflag = 103
            return
        end if
        found = .false.
        do
            found = (xi<t(left+1))
            if (found) exit
            left = left + 1
            if (left>=ilp1mx) exit
        end do
        if (.not. found) then
            left = left - 1
            if (xi>t(left+1)) then
                !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
                !             ' corresponding basis function and the system is singular'
                iflag = 103
                return
            end if
        end if
        ! the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
        call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork, iflag)
        if (iflag/=0) return

        ! we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
        ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
        ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
        ! as a two-dim. array , with  2*k-1  rows (see comments in
        ! dbnfac). in the present program, we treat  q  as an equivalent
        ! one-dimensional array (because of fortran restrictions on
        ! dimension statements) . we therefore want  bcoef(j) to go into
        ! entry
        !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
        !            = i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
        ! of q.
        jj = i - left + 1 + (left-k)*(k+km1)
        do j=1,k
            jj = jj + kpkm2
            q(jj) = bcoef(j)
        end do

    end do

    ! obtain factorization of a, stored again in q.
    call dbnfac(q, k+km1, n, km1, km1, iflag)

    if (iflag==1) then !success
        ! solve  a*bcoef = y  by backsubstitution
        do i=1,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
        iflag = 0
    else  !failure
        !write(error_unit,'(A)') 'dbintk - the system of solver detects a singular system'//&
        !             ' although the theoretical conditions for a solution were satisfied'
        iflag = 104
    end if

    end subroutine dbintk
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns in w the LU-factorization (without pivoting) of the banded
!  matrix a of order nrow with (nbandl + 1 + nbandu) bands or diagonals
!  in the work array w .
!
!  gauss elimination without pivoting is used. the routine is
!  intended for use with matrices a which do not require row inter-
!  changes during factorization, especially for the totally
!  positive matrices which occur in spline calculations.
!  the routine should not be used for an arbitrary banded matrix.
!
!### Work array
!
! **Input**
!
!        w array of size (nroww,nrow) contains the interesting
!        part of a banded matrix  a , with the diagonals or bands of  a
!        stored in the rows of  w , while columns of  a  correspond to
!        columns of  w . this is the storage mode used in  linpack  and
!        results in efficient innermost loops.
!           explicitly,  a  has  nbandl  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   nbandu  bands above the diagonal
!        and thus, with    middle = nbandu + 1,
!          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
!                                              j=1,...,nrow .
!        for example, the interesting entries of a (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  w
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        all other entries of  w  not identified in this way with an en-
!        try of  a  are never referenced .
!
! **Output**
!
!  * if  iflag = 1, then
!        w contains the lu-factorization of  a  into a unit lower triangu-
!        lar matrix  l  and an upper triangular matrix  u (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  a . this makes it possible to solve any particular linear
!        system  a*x = b  for  x  by a
!              call dbnslv ( w, nroww, nrow, nbandl, nbandu, b )
!        with the solution x  contained in  b  on return .
!  * if  iflag = 2, then
!        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  a  does not have an lu-factorization. this implies that
!        a  is singular in case it is totally positive .
!
!### History
!  * banfac written by carl de boor [5]
!  * dbnfac from CMLIB [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)

    integer,intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
    integer,intent(in) :: nrow    !! matrix order
    integer,intent(in) :: nbandl  !! number of bands of a below the main diagonal
    integer,intent(in) :: nbandu  !! number of bands of a above the main diagonal
    integer,intent(out) :: iflag  !! indicating success(=1) or failure (=2)
    real(8),dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

    integer :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
    real(8) :: factor, pivot

    iflag = 1
    middle = nbandu + 1   ! w(middle,.) contains the main diagonal of a.
    nrowm1 = nrow - 1

    if (nrowm1 < 0) then
        iflag = 2
        return
    elseif (nrowm1 == 0) then
        if (w(middle,nrow)==0.0d0) iflag = 2
        return
    end if

    if (nbandl<=0) then
        ! a is upper triangular. check that diagonal is nonzero .
        do i=1,nrowm1
            if (w(middle,i)==0.0d0) then
                iflag = 2
                return
            end if
        end do
        if (w(middle,nrow)==0.0d0) iflag = 2
        return
    end if

    if (nbandu<=0) then
        ! a is lower triangular. check that diagonal is nonzero and
        ! divide each column by its diagonal.
        do i=1,nrowm1
            pivot = w(middle,i)
            if (pivot==0.0d0) then
                iflag = 2
                return
            end if
            jmax = min(nbandl,nrow-i)
            do j=1,jmax
                w(middle+j,i) = w(middle+j,i)/pivot
            end do
        end do
        return
    end if

    ! a is not just a triangular matrix. construct lu factorization
    do i=1,nrowm1
        ! w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0d0) then
            iflag = 2
            return
        end if
        ! jmax is the number of (nonzero) entries in column i
        ! below the diagonal.
        jmax = min(nbandl,nrow-i)
        ! divide each entry in column i below diagonal by pivot.
        do j=1,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
        end do
        ! kmax is the number of (nonzero) entries in row i to
        ! the right of the diagonal.
        kmax = min(nbandu,nrow-i)
        ! subtract a(i,i+k)*(i-th column) from (i+k)-th column
        ! (below row i).
        do k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1,jmax
                w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
        end do
    end do

    ! check the last diagonal entry.
    if (w(middle,nrow)==0.0d0) iflag = 2

    end subroutine dbnfac
!*****************************************************************************************

!*****************************************************************************************
!>
!  Companion routine to [[dbnfac]]. it returns the solution x of the
!  linear system a*x = b in place of b, given the lu-factorization
!  for a in the work array w from dbnfac.
!
!  (with \( a = l*u \), as stored in w), the unit lower triangular system
!  \( l(u*x) = b \) is solved for \( y = u*x \), and y stored in b. then the
!  upper triangular system \(u*x = y \) is solved for x. the calculations
!  are so arranged that the innermost loops stay within columns.
!
!### History
!  * banslv written by carl de boor [5]
!  * dbnslv from SLATEC library [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)

    integer,intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    integer,intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    integer,intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    integer,intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    real(8),dimension(nroww,nrow),intent(in) :: w    !! describes the lu-factorization of a banded matrix a of order `nrow` as constructed in [[dbnfac]].
    real(8),dimension(nrow),intent(inout) :: b  !! * **in**: right side of the system to be solved
                                                 !! * **out**: the solution x, of order nrow

    integer :: i, j, jmax, middle, nrowm1

    middle = nbandu + 1
    if (nrow/=1) then

        nrowm1 = nrow - 1
        if (nbandl/=0) then

            ! forward pass
            ! for i=1,2,...,nrow-1, subtract right side(i)*(i-th column of l)
            !                       from right side (below i-th row).
            do i=1,nrowm1
                jmax = min(nbandl,nrow-i)
                do j=1,jmax
                    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
                end do
            end do

        end if

        ! backward pass
        ! for i=nrow,nrow-1,...,1, divide right side(i) by i-th diagonal
        !                          entry of u, then subtract right side(i)*(i-th column
        !                          of u) from right side (above i-th row).
        if (nbandu<=0) then
            ! a is lower triangular.
            do i=1,nrow
                b(i) = b(i)/w(1,i)
            end do
            return
        end if

        i = nrow
        do
            b(i) = b(i)/w(middle,i)
            jmax = min(nbandu,i-1)
            do j=1,jmax
                b(i-j) = b(i-j) - b(i)*w(middle-j,i)
            end do
            i = i - 1
            if (i<=1) exit
        end do

    end if

    b(1) = b(1)/w(middle,1)

    end subroutine dbnslv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculates the value of all (possibly) nonzero basis
!  functions at x of order max(jhigh,(j+1)*(index-1)), where t(k)
!  <= x <= t(n+1) and j=iwork is set inside the routine on
!  the first call when index=1.  ileft is such that t(ileft) <=
!  x < t(ileft+1).  a call to dintrv(t,n+1,x,ilo,ileft,mflag)
!  produces the proper ileft.  dbspvn calculates using the basic
!  algorithm needed in dbspvd.  if only basis functions are
!  desired, setting jhigh=k and index=1 can be faster than
!  calling dbspvd, but extra coding is required for derivatives
!  (index=2) and dbspvd is set up for this purpose.
!
!  left limiting values are set up as described in dbspvd.
!
!### Error Conditions
!
!  * improper input
!
!### History
!  * bsplvn written by carl de boor [5]
!  * dbspvn author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork,iflag)

    implicit none

    real(8),dimension(*),intent(in)  :: t        !! knot vector of length n+k, where
                                                  !! n = number of b-spline basis functions
                                                  !! n = sum of knot multiplicities-k
                                                  !! dimension t(ileft+jhigh)
    integer,intent(in)                :: jhigh    !! order of b-spline, 1 <= jhigh <= k
    integer,intent(in)                :: k        !! highest possible order
    integer,intent(in)                :: index    !! index = 1 gives basis functions of order jhigh
                                                  !!       = 2 denotes previous entry with work, iwork
                                                  !!         values saved for subsequent calls to
                                                  !!         dbspvn.
    real(8),intent(in)               :: x        !! argument of basis functions, t(k) <= x <= t(n+1)
    integer,intent(in)                :: ileft    !! largest integer such that t(ileft) <= x < t(ileft+1)
    real(8),dimension(k),intent(out) :: vnikx    !! vector of length k for spline values.
    real(8),dimension(*),intent(out) :: work     !! a work vector of length 2*k
    integer,intent(out)               :: iwork    !! a work parameter.  both work and iwork contain
                                                  !! information necessary to continue for index = 2.
                                                  !! when index = 1 exclusively, these are scratch
                                                  !! variables and can be used for other purposes.
    integer,intent(out)               :: iflag    !!   0: no errors
                                                  !! 201: k does not satisfy k>=1
                                                  !! 202: jhigh does not satisfy 1<=jhigh<=k
                                                  !! 203: index is not 1 or 2
                                                  !! 204: x does not satisfy t(ileft)<=x<=t(ileft+1)

    integer :: imjp1, ipj, jp1, jp1ml, l
    real(8) :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i),
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1) then
        !write(error_unit,'(A)') 'dbspvn - k does not satisfy k>=1'
        iflag = 201
        return
    end if
    if (jhigh>k .or. jhigh<1) then
        !write(error_unit,'(A)') 'dbspvn - jhigh does not satisfy 1<=jhigh<=k'
        iflag = 202
        return
    end if
    if (index<1 .or. index>2) then
        !write(error_unit,'(A)') 'dbspvn - index is not 1 or 2'
        iflag = 203
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1)) then
        !write(error_unit,'(A)') 'dbspvn - x does not satisfy t(ileft)<=x<=t(ileft+1)'
        iflag = 204
        return
    end if

    iflag = 0

    if (index==1) then
        iwork = 1
        vnikx(1) = 1.0d0
        if (iwork>=jhigh) return
    end if

    do
        ipj = ileft + iwork
        work(iwork) = t(ipj) - x
        imjp1 = ileft - iwork + 1
        work(k+iwork) = x - t(imjp1)
        vmprev = 0.0d0
        jp1 = iwork + 1
        do l=1,iwork
            jp1ml = jp1 - l
            vm = vnikx(l)/(work(l)+work(k+jp1ml))
            vnikx(l) = vm*work(l) + vmprev
            vmprev = vm*work(k+jp1ml)
        end do
        vnikx(jp1) = vmprev
        iwork = jp1
        if (iwork>=jhigh) exit
    end do

    end subroutine dbspvn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluates the b-representation (`t`,`a`,`n`,`k`) of a b-spline
!  at `x` for the function value on `ideriv=0` or any of its
!  derivatives on `ideriv=1,2,...,k-1`.  right limiting values
!  (right derivatives) are returned except at the right end
!  point `x=t(n+1)` where left limiting values are computed.  the
!  spline is defined on `t(k)` \( \le \) `x` \( \le \) `t(n+1)`.
!  dbvalu returns a fatal error message when `x` is outside of this
!  interval.
!
!  To compute left derivatives or left limiting values at a
!  knot `t(i)`, replace `n` by `i-1` and set `x=t(i), i=k+1,n+1`.
!
!### Error Conditions
!
!  * improper input
!
!### History
!  * bvalue written by carl de boor [5]
!  * dbvalu author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbvalu(t,a,n,k,ideriv,x,inbv,work,iflag,val,extrap)

    implicit none

    real(8),intent(out)             :: val     !! the interpolated value
    integer,intent(in)               :: n       !! number of b-spline coefficients.
                                                !! (sum of knot multiplicities-`k`)
    real(8),dimension(:),intent(in) :: t       !! knot vector of length `n+k`
    real(8),dimension(n),intent(in) :: a       !! b-spline coefficient vector of length `n`
    integer,intent(in)               :: k       !! order of the b-spline, `k >= 1`
    integer,intent(in)               :: ideriv  !! order of the derivative, `0 <= ideriv <= k-1`.
                                                !! `ideriv = 0` returns the b-spline value
    real(8),intent(in)              :: x       !! argument, `t(k) <= x <= t(n+1)`
    integer,intent(inout)            :: inbv    !! an initialization parameter which must be set
                                                !! to 1 the first time [[dbvalu]] is called.
                                                !! `inbv` contains information for efficient processing
                                                !! after the initial call and `inbv` must not
                                                !! be changed by the user.  distinct splines require
                                                !! distinct `inbv` parameters.
    real(8),dimension(:),intent(inout) :: work !! work vector of length at least `3*k`
    integer,intent(out)              :: iflag   !! status flag:
                                                !!
                                                !! * 0: no errors
                                                !! * 401: `k` does not satisfy `k` \( \ge \) 1
                                                !! * 402: `n` does not satisfy `n` \( \ge \) `k`
                                                !! * 403: `ideriv` does not satisfy 0 \( \le \) `ideriv` \(<\) `k`
                                                !! * 404: `x` is not greater than or equal to `t(k)`
                                                !! * 405: `x` is not less than or equal to `t(n+1)`
                                                !! * 406: a left limiting value cannot be obtained at `t(k)`
    logical,intent(in),optional :: extrap   !! if extrapolation is allowed
                                            !! (if not present, default is False)

    integer :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real(8) :: fkmj
    real(8) :: xt
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    ! make a temp copy of x (for computing the
    ! interval) in case extrapolation is allowed
    if (extrapolation_allowed) then
        if (x<t(1)) then
            xt = t(1)
        else if(x>t(n+k)) then
            xt = t(n+k)
        else
            xt = x
        end if
    else
        xt = x
    end if

    val = 0.0d0

    if (k<1) then
        iflag = 401  ! dbvalu - k does not satisfy k>=1
        return
    end if

    if (n<k) then
        iflag = 402  ! dbvalu - n does not satisfy n>=k
        return
    end if

    if (ideriv<0 .or. ideriv>=k) then
        iflag = 403  ! dbvalu - ideriv does not satisfy 0<=ideriv<k
        return
    end if

    kmider = k - ideriv

    ! find *i* in (k,n) such that t(i) <= x < t(i+1)
    ! (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1
    call dintrv(t, n+1, xt, inbv, i, mflag)
    if (xt<t(k)) then
        iflag = 404  ! dbvalu - x is not greater than or equal to t(k)
        return
    end if

    if (mflag/=0) then

        if (xt>t(i)) then
            iflag = 405  ! dbvalu - x is not less than or equal to t(n+1)
            return
        end if

        do
            if (i==k) then
                iflag = 406  ! dbvalu - a left limiting value cannot be obtained at t(k)
                return
            end if
            i = i - 1
            if (xt/=t(i)) exit
        end do

    end if

    ! difference the coefficients *ideriv* times
    ! work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

    imk = i - k
    do j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
    end do

    if (ideriv/=0) then
        do j=1,ideriv
            kmj = k - j
            fkmj = dble(kmj)
            do jj=1,kmj
                ihi = i + jj
                ihmkmj = ihi - kmj
                work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
            end do
        end do
    end if

    ! compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    ! given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

    if (ideriv/=km1) then
        ip1 = i + 1
        kpk = k + k
        j1 = k + 1
        j2 = kpk + 1
        do j=1,kmider
            ipj = i + j
            work(j1) = t(ipj) - x
            ip1mj = ip1 - j
            work(j2) = x - t(ip1mj)
            j1 = j1 + 1
            j2 = j2 + 1
        end do
        iderp1 = ideriv + 1
        do j=iderp1,km1
            kmj = k - j
            ilo = kmj
            do jj=1,kmj
                work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)*&
                            work(k+jj))/(work(kpk+ilo)+work(k+jj))
                ilo = ilo - 1
            end do
        end do
    end if

    iflag = 0
    val = work(1)

    end subroutine dbvalu
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the largest integer `ileft` in 1 \( \le \) `ileft` \( \le \) `lxt`
!  such that `xt(ileft)` \( \le \) `x` where `xt(*)` is a subdivision of
!  the `x` interval.
!  precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   mflag=0
!         if xt(lxt) <= x           then ileft=lxt, mflag=1
!```
!
!  that is, when multiplicities are present in the break point
!  to the left of `x`, the largest index is taken for `ileft`.
!
!### History
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).
!  * Jacob Williams, 3/4/2017 : added extrapolation option.

    pure subroutine dintrv(xt,lxt,xx,ilo,ileft,mflag,extrap)

    implicit none

    integer,intent(in)                 :: lxt    !! length of the `xt` vector
    real(8),dimension(lxt),intent(in) :: xt     !! a knot or break point vector of length `lxt`
    real(8),intent(in)                :: xx     !! argument
    integer,intent(inout)              :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  distinct splines
                                                 !! require distinct `ilo` parameters.
    integer,intent(out)                :: ileft  !! largest integer satisfying `xt(ileft)` \( \le \) `x`
    integer,intent(out)                :: mflag  !! signals when `x` lies out of bounds
    logical,intent(in),optional        :: extrap !! if extrapolation is allowed
                                                 !! (if not present, default is False)

    integer :: ihi, istep, middle
    real(8) :: x

    x = get_temp_x_for_extrap(xx,xt,extrap)

    ihi = ilo + 1
    if ( ihi>=lxt ) then
        if ( x>=xt(lxt) ) then
            mflag = 1
            ileft = lxt
            return
        end if
        if ( lxt<=1 ) then
            mflag = -1
            ileft = 1
            return
        end if
        ilo = lxt - 1
        ihi = lxt
    endif

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=lxt ) then
                if ( x>=xt(lxt) ) then
                    mflag = 1
                    ileft = lxt
                    return
                end if
                ihi = lxt
            elseif ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    return
                end if
            elseif ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    endif

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        middle = (ilo+ihi)/2
        if ( middle==ilo ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! note. it is assumed that middle = ilo in case ihi = ilo+1
        if ( x<xt(middle) ) then
            ihi = middle
        else
            ilo = middle
        endif
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the value of `x` to use for computing the interval
!  in `t`, depending on if extrapolation is allowed or not.
!
!  If extrapolation is allowed and x is > t(1) or < t(n), then either
!  t(1) or t(n) is returned. Otherwise, `x` is returned.

    pure function get_temp_x_for_extrap(x,t,extrap) result(xt)

    implicit none

    real(8),intent(in)              :: x       !! variable value
    real(8),dimension(:),intent(in) :: t       !! knot vector for b-splines
    real(8)                         :: xt      !! The value returned (it will either
                                                !! be `t(1)`, `x`, or `t(n)`)
    logical,intent(in),optional      :: extrap  !! if extrapolation is allowed
                                                !! (if not present, default is False)

    integer :: n  !! size of `t`
    logical :: extrapolation_allowed  !! if extrapolation is allowed

    if (present(extrap)) then
        extrapolation_allowed = extrap
    else
        extrapolation_allowed = .false.
    end if

    n = size(t)

    if (extrapolation_allowed) then
        if (x<t(1)) then
            xt = t(1)
        elseif (x>t(n)) then
            xt = t(n)
        else
            xt = x
        end if
    else
        xt = x
    end if

    end function get_temp_x_for_extrap
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a message string associated with the status code.

    pure function get_status_message(iflag) result(msg)

    implicit none

    integer,intent(in)           :: iflag  !! return code from one of the routines
    character,allocatable :: msg(:)    !! status message associated with the flag

    character(len=10) :: istr   !! for integer to string conversion
    integer           :: istat  !! for write statement

    select case (iflag)

    case(  0); msg='Successful execution'

    case(  1); msg='Error in evaluate_*d: class is not initialized'

    case(  2); msg='Error in db*ink: iknot out of range'
    case(  3); msg='Error in db*ink: nx out of range'
    case(  4); msg='Error in db*ink: kx out of range'
    case(  5); msg='Error in db*ink: x not strictly increasing'
    case(  6); msg='Error in db*ink: tx not non-decreasing'
    case(  7); msg='Error in db*ink: ny out of range'
    case(  8); msg='Error in db*ink: ky out of range'
    case(  9); msg='Error in db*ink: y not strictly increasing'
    case( 10); msg='Error in db*ink: ty not non-decreasing'
    case( 11); msg='Error in db*ink: nz out of range'
    case( 12); msg='Error in db*ink: kz out of range'
    case( 13); msg='Error in db*ink: z not strictly increasing'
    case( 14); msg='Error in db*ink: tz not non-decreasing'
    case( 15); msg='Error in db*ink: nq out of range'
    case( 16); msg='Error in db*ink: kq out of range'
    case( 17); msg='Error in db*ink: q not strictly increasing'
    case( 18); msg='Error in db*ink: tq not non-decreasing'
    case( 19); msg='Error in db*ink: nr out of range'
    case( 20); msg='Error in db*ink: kr out of range'
    case( 21); msg='Error in db*ink: r not strictly increasing'
    case( 22); msg='Error in db*ink: tr not non-decreasing'
    case( 23); msg='Error in db*ink: ns out of range'
    case( 24); msg='Error in db*ink: ks out of range'
    case( 25); msg='Error in db*ink: s not strictly increasing'
    case( 26); msg='Error in db*ink: ts not non-decreasing'
    case(700); msg='Error in db*ink: size(x) /= size(fcn,1)'
    case(701); msg='Error in db*ink: size(y) /= size(fcn,2)'
    case(702); msg='Error in db*ink: size(z) /= size(fcn,3)'
    case(703); msg='Error in db*ink: size(q) /= size(fcn,4)'
    case(704); msg='Error in db*ink: size(r) /= size(fcn,5)'
    case(705); msg='Error in db*ink: size(s) /= size(fcn,6)'
    case(706); msg='Error in db*ink: size(x) /= nx'
    case(707); msg='Error in db*ink: size(y) /= ny'
    case(708); msg='Error in db*ink: size(z) /= nz'
    case(709); msg='Error in db*ink: size(q) /= nq'
    case(710); msg='Error in db*ink: size(r) /= nr'
    case(711); msg='Error in db*ink: size(s) /= ns'
    case(712); msg='Error in db*ink: size(tx) /= nx+kx'
    case(713); msg='Error in db*ink: size(ty) /= ny+ky'
    case(714); msg='Error in db*ink: size(tz) /= nz+kz'
    case(715); msg='Error in db*ink: size(tq) /= nq+kq'
    case(716); msg='Error in db*ink: size(tr) /= nr+kr'
    case(717); msg='Error in db*ink: size(ts) /= ns+ks'
    case(800); msg='Error in db*ink: size(x) /= size(bcoef,1)'
    case(801); msg='Error in db*ink: size(y) /= size(bcoef,2)'
    case(802); msg='Error in db*ink: size(z) /= size(bcoef,3)'
    case(803); msg='Error in db*ink: size(q) /= size(bcoef,4)'
    case(804); msg='Error in db*ink: size(r) /= size(bcoef,5)'
    case(805); msg='Error in db*ink: size(s) /= size(bcoef,6)'

    case(100); msg='Error in dbintk: k does not satisfy k>=1'
    case(101); msg='Error in dbintk: n does not satisfy n>=k'
    case(102); msg='Error in dbintk: x(i) does not satisfy x(i)<x(i+1) for some i'
    case(103); msg='Error in dbintk: some abscissa was not in the support of the '//&
                    'corresponding basis function and the system is singular'
    case(104); msg='Error in dbintk: the system of solver detects a singular system '//&
                   'although the theoretical conditions for a solution were satisfied'

    case(201); msg='Error in dbspvn: k does not satisfy k>=1'
    case(202); msg='Error in dbspvn: jhigh does not satisfy 1<=jhigh<=k'
    case(203); msg='Error in dbspvn: index is not 1 or 2'
    case(204); msg='Error in dbspvn: x does not satisfy t(ileft)<=x<=t(ileft+1)'

    case(301); msg='Error in dbtpcf: n should be > 0'

    case(401); msg='Error in dbvalu: k does not satisfy k>=1'
    case(402); msg='Error in dbvalu: n does not satisfy n>=k'
    case(403); msg='Error in dbvalu: ideriv does not satisfy 0<=ideriv<k'
    case(404); msg='Error in dbvalu: x is not greater than or equal to t(k)'
    case(405); msg='Error in dbvalu: x is not less than or equal to t(n+1)'
    case(406); msg='Error in dbvalu: a left limiting value cannot be obtained at t(k)'

    case(501); msg='Error in initialize_*d_specify_knots: tx is not the correct size (kx+nx)'
    case(502); msg='Error in initialize_*d_specify_knots: ty is not the correct size (ky+ny)'
    case(503); msg='Error in initialize_*d_specify_knots: tz is not the correct size (kz+nz)'
    case(504); msg='Error in initialize_*d_specify_knots: tq is not the correct size (kq+nq)'
    case(505); msg='Error in initialize_*d_specify_knots: tr is not the correct size (kr+nr)'
    case(506); msg='Error in initialize_*d_specify_knots: ts is not the correct size (ks+ns)'

    case(601); msg='Error in db*val: x value out of bounds'
    case(602); msg='Error in db*val: y value out of bounds'
    case(603); msg='Error in db*val: z value out of bounds'
    case(604); msg='Error in db*val: q value out of bounds'
    case(605); msg='Error in db*val: r value out of bounds'
    case(606); msg='Error in db*val: s value out of bounds'

    case default
        write(istr,fmt='(I10)',iostat=istat) iflag
        msg = 'Unknown status flag: '//trim(adjustl(istr))
    end select

    end function get_status_message
!*****************************************************************************************

!*****************************************************************************************
    end module bspline_sub_module
!*****************************************************************************************
