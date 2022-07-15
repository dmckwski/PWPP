      module solver
      logical :: initialize_constants,every_iteration_cv
      integer :: ndx_cv,ndy_cv,pol_cv,run_print_unit_cv
      integer, public :: world_numprocs,world_rank
      real(8) :: k0x_cv,k0y_cv,wx_cv,wy_cv,h_cv,dw_cv,d_cv
      real(8), private :: pi
      complex(8) :: rimedium_cv,ri0_cv,ris_cv,rpar0_cv, &
          rper0_cv,rpars_cv,rpers_cv
!      complex(8), allocatable :: epsilon_cv(:,:)
      complex(4), allocatable :: epsilon_cv(:,:)
      data pi/3.141592653589793/
      data run_print_unit_cv/6/
      data world_rank/0/
      data every_iteration_cv/.false./
      data initialize_constants/.true./

      contains

         subroutine epsilon_mult(ndsx,ndsy,ml,gcoef,acoef, &
             eps,con_tran,kp_x,kp_y)
         use fft
         implicit none
         logical, optional :: con_tran
         logical :: contran
         integer :: ndsx,ndsy,ml,is1x,is2x,is1y,is2y,ndpx,ndpy,ip1x,ip2x,ip1y,ip2y,nds2,ndp2,j,m,isign
         real(8) :: kpx,kpy
         real(8), optional :: kp_x,kp_y
         complex(8) :: gcoef(3,ndsx*ndsy,ml),acoef(3,ndsx*ndsy,ml)
         complex(4) :: eps(min(ndx_cv*ndy_cv,4*ndsx*ndsy),ml)
         complex(8), allocatable :: gt(:,:,:),gt2(:,:)
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif
         if(present(kp_x)) then
            kpx=kp_x
         else
            kpx=0.d0
         endif
         if(present(kp_y)) then
            kpy=kp_y
         else
            kpy=0.d0
         endif
         if((ndx_cv.eq.1.and.ndy_cv.eq.1).or.(ndsx.eq.1.and.ndsy.eq.1)) then
            if(contran) then
               do j=1,ml
                  acoef(:,1,j)=conjg(eps(1,j))*gcoef(:,1,j)
               enddo
            else
               do j=1,ml
                  acoef(:,1,j)=eps(1,j)*gcoef(:,1,j)
               enddo
            endif
            return
         endif
         isign=1
         is1x=-ndsx/2
         is2x=(ndsx-1)/2
         is1y=-ndsy/2
         is2y=(ndsy-1)/2
         ndpx=min(2*ndsx,ndx_cv)
         ndpy=min(2*ndsy,ndy_cv)
         ip1x=-ndpx/2
         ip2x=(ndpx-1)/2
         ip1y=-ndpy/2
         ip2y=(ndpy-1)/2
         nds2=ndsx*ndsy
         ndp2=ndpx*ndpy
         allocate(gt(3,ip1x:ip2x,ip1y:ip2y),gt2(3,ndp2))
         do j=1,ml
            gt=0.d0
            gt(:,is1x:is2x,is1y:is2y)=reshape(gcoef(1:3,1:nds2,j),(/3,ndsx,ndsy/))
            call fft2d(3,ndpx,ndpy,wx_cv,wy_cv,0.d0,0.d0,gt,gt2,isign)
            if(contran) then
               do m=1,3
                  gt2(m,1:ndp2)=conjg(eps(1:ndp2,j))*gt2(m,1:ndp2)
               enddo
            else
               do m=1,3
                  gt2(m,1:ndp2)=eps(1:ndp2,j)*gt2(m,1:ndp2)
               enddo
            endif
            call fft2d(3,ndpx,ndpy,wx_cv,wy_cv,kpx,kpy,gt2,gt,-isign)
            acoef(:,1:nds2,j)=reshape(gt(1:3,is1x:is2x,is1y:is2y),(/3,nds2/))
         enddo
         deallocate(gt,gt2)
         end subroutine epsilon_mult

         subroutine gfunc_mult(ndsx,ndsy,ml,acoef,g_coef,con_tran, &
            p_coef)
         implicit none
         logical, optional :: con_tran
         logical :: contran
         integer :: ndsx,ndsy,ml,is1x,is2x,is1y,is2y,s,sy,sx,j,nds2
         real(8) :: ky,kx,kr2,kr,kxh,kyh,rmat(2,2),h,maxfunc
         complex(8) :: acoef(3,ndsx*ndsy,ml), &
                       kz,skz,ci,m,m0,ms,kz0,kzs, &
                       kz2,ckzmh,skzmh,f1den,f2den,expd2, &
                       expd,expj,upar(ml+1),uper(ml+1), &
                       dpar(0:ml),dper(0:ml),gc(3,ml),ac(3,ml), &
                       f0par,df0par,uparj,uperj,dparj,dperj, &
                       fhpar,dfhpar
         complex(8), optional :: g_coef(3,ndsx*ndsy,ml), &
            p_coef(2,ndsx*ndsy,2)
         complex(8), allocatable, save :: vf0par(:,:,:), &
            vfhpar(:,:,:),f0per(:,:),fhper(:,:)
         data ci,maxfunc/(0.d0,1.d0),0.d0/
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif
         is1x=-ndsx/2
         is2x=(ndsx-1)/2
         is1y=-ndsy/2
         is2y=(ndsy-1)/2
         nds2=ndsx*ndsy
         m=rimedium_cv
         m0=ri0_cv
         ms=ris_cv
         h=h_cv

         if(initialize_constants) then
            if(allocated(vf0par)) then
               deallocate(vf0par,vfhpar,f0per,fhper)
            endif
            allocate(vf0par(2,0:ml,nds2), &
               vfhpar(2,ml+1,nds2), &
               f0per(0:ml,nds2), &
               fhper(ml+1,nds2))
            s=0
            do sy=is1y,is2y
               ky=k0y_cv+2.d0*pi*dble(sy)/wy_cv
               do sx=is1x,is2x
                  s=s+1
                  kx=k0x_cv+2.d0*pi*dble(sx)/wx_cv
                  kr2=kx*kx+ky*ky
                  kr=sqrt(kr2)
                  kz=cdsqrt((1.d0,0.d0)-kr2/m/m)
                  kz0=cdsqrt((1.d0,0.d0)-kr2/m0/m0)
                  kzs=cdsqrt((1.d0,0.d0)-kr2/ms/ms)
                  skz=cdsqrt(1.d0-kz*kz)
                  kz2=kz*kz

                  ckzmh=.5d0*(1.d0+cdexp(2.d0*(0.d0,1.d0)*kz*m*h))
                  skzmh=.5d0*(0.d0,1.d0)*(1.d0-cdexp(2.d0*(0.d0,1.d0)*kz*m*h))
                  f1den=cdsqrt(d_cv*kz/(m*((0.d0,1.d0)*kz*m*(kzs*m0 + kz0*ms)*ckzmh &
                     + (kz0*m*kzs*m + kz*m0*kz*ms)*skzmh)))
                  f2den=cdsqrt(d_cv/(kz*m*((0.d0,1.d0)*kz*m*(kz0*m0 + kzs*ms)*ckzmh &
                    + (kz*kz*m*m + kz0*kzs*m0*ms)*skzmh)))

                  expd2=cdexp((0.d0,1.d0)*d_cv*kz*m/2.d0)
                  expd=expd2*expd2
                  expj=expd2
                  do j=1,ml
                     f0par=.5d0*(0.d0,1.d0) &
                        *(m*kz0+m0*kz+expj*expj*(m*kz0-m0*kz))
                     df0par=.5d0*m*kz &
                        *(m*kz0+m0*kz+expj*expj*(-m*kz0+m0*kz))
                     vf0par(1:2,j,s)=(/f0par, &
                        -(0.d0,1.d0)*skz*df0par/(m*kz2)/)*f1den
                     f0per(j,s)=.5d0*(m*kz+m0*kz0+expj*expj*(m*kz-m0*kz0))*f2den
                     expj=expj*expd
                  enddo
                  expj=expd2
                  do j=ml,1,-1
                     fhpar=.5d0*(0.d0,1.d0) &
                        *(m*kzs+ms*kz+expj*expj*(m*kzs-ms*kz))
                     dfhpar=.5d0*m*kz &
                        *(m*kzs+ms*kz+expj*expj*(-m*kzs+ms*kz))
                     vfhpar(1:2,j,s)=(/fhpar, &
                        -(0.d0,1.d0)*skz*dfhpar/(m*kz2)/)*f1den
                     fhper(j,s)=-.5d0*(m*kz+ms*kzs+expj*expj*(m*kz-ms*kzs))*f2den
                     expj=expj*expd
                  enddo
                  vfhpar(1,ml+1,s)=(0.d0,1.d0)*m*kzs*f1den
                  fhper(ml+1,s)=-m*kz*f2den
                  vf0par(1,0,s)=(0.d0,1.d0)*m*kz0*f1den
                  f0per(0,s)=m*kz*f2den
               enddo
            enddo
            initialize_constants=.false.
!write(*,'(e13.5)') maxfunc
         endif

         s=0
         do sy=is1y,is2y
            ky=k0y_cv+2.d0*pi*dble(sy)/wy_cv
            do sx=is1x,is2x
               s=s+1
               kx=k0x_cv+2.d0*pi*dble(sx)/wx_cv
               kr2=kx*kx+ky*ky
               kr=sqrt(kr2)
               if(kr2.eq.0.d0) then
                  kxh=1.d0
                  kyh=0.d0
               else
                  kxh=kx/kr
                  kyh=ky/kr
               endif
               kz=cdsqrt((1.d0,0.d0)-kr2/m/m)

               rmat(1,:)=(/kxh,kyh/)
               rmat(2,:)=(/kyh,-kxh/)
               do j=1,ml
                  gc((/1,2/),j)=matmul(rmat,acoef((/1,2/),s,j))
                  gc(3,j)=acoef(3,s,j)
               enddo
               if(contran) then
                  gc=dconjg(gc)
                  gc(3,:)=-gc(3,:)
               endif

               expd2=cdexp((0.d0,.5d0)*d_cv*kz*m)
               expd=expd2*expd2
               upar(1)=0.d0
               uper(1)=0.d0
               do j=2,ml
                  upar(j)=expd*(vf0par(1,j-1,s)*gc(1,j-1) &
                     +vf0par(2,j-1,s)*gc(3,j-1) &
                     +upar(j-1))
                  uper(j)=expd*(f0per(j-1,s)*gc(2,j-1) &
                     +uper(j-1))
               enddo
               dpar(ml)=0.d0
               dper(ml)=0.d0
               do j=ml-1,1,-1
                  dpar(j)=expd*(vfhpar(1,j+1,s)*gc(1,j+1) &
                     -vfhpar(2,j+1,s)*gc(3,j+1) &
                     +dpar(j+1))
                  dper(j)=expd*(fhper(j+1,s)*gc(2,j+1) &
                     +dper(j+1))
               enddo

               if(present(g_coef)) then
                  do j=1,ml
                     uparj=upar(j)+.5d0*(vf0par(1,j,s)*gc(1,j) &
                        +vf0par(2,j,s)*gc(3,j))
                     uperj=uper(j)+.5d0*(f0per(j,s)*gc(2,j))
                     dparj=dpar(j)+.5d0*(vfhpar(1,j,s)*gc(1,j) &
                        -vfhpar(2,j,s)*gc(3,j))
                     dperj=dper(j)+.5d0*(fhper(j,s)*gc(2,j))
                     ac(1,j)=vfhpar(1,j,s)*uparj &
                         +vf0par(1,j,s)*dparj
                     ac(2,j)=fhper(j,s)*uperj &
                         +f0per(j,s)*dperj
                     ac(3,j)=vfhpar(2,j,s)*uparj &
                         -vf0par(2,j,s)*dparj &
                         -gc(3,j)/m/m
                  enddo
                  if(contran) then
                     ac=dconjg(ac)
                     ac(3,:)=-ac(3,:)
                  endif
                  do j=1,ml
                     g_coef((/1,2/),s,j)=matmul(rmat,ac((/1,2/),j))
                     g_coef(3,s,j)=ac(3,j)
                  enddo
               endif
               if(present(p_coef)) then
                  kz0=cdsqrt((1.d0,0.d0)-kr2/m0/m0)
                  kzs=cdsqrt((1.d0,0.d0)-kr2/ms/ms)
                  upar(ml)=expd2*(upar(ml)+vf0par(1,ml,s)*gc(1,ml) &
                     +vf0par(2,ml,s)*gc(3,ml))
                  uper(ml)=expd2*(uper(ml)+f0per(ml,s)*gc(2,ml))
                  dpar(1)=expd2*(dpar(1)+vfhpar(1,1,s)*gc(1,1) &
                     -vfhpar(2,1,s)*gc(3,1))
                  dper(1)=expd2*(dper(1)+fhper(1,s)*gc(2,1))
                  ac(1,ml)=vfhpar(1,ml+1,s)*upar(ml)
                  ac(2,ml)=fhper(ml+1,s)*uper(ml)
                  p_coef(1,s,2)=ac(1,ml)/kzs
                  p_coef(2,s,2)=ac(2,ml)
                  ac(1,1)=vf0par(1,0,s)*dpar(1)
                  ac(2,1)=f0per(0,s)*dper(1)
                  p_coef(1,s,1)=ac(1,1)/kz0
                  p_coef(2,s,1)=ac(2,1)
               endif
            enddo
         enddo
         end subroutine gfunc_mult

         subroutine interaction_matrix_mult(ndsx,ndsy,ml,gcoef,acoef,con_tran)
         use mpidefs
         implicit none
         logical, optional ::con_tran
         logical :: contran
         integer :: ndsx,ndsy,ml
!         real(8) :: starttime
         complex(8) :: acoef(3,ndsx*ndsy,ml),gcoef(3,ndsx*ndsy,ml), &
             acoeft(3,ndsx*ndsy,ml)
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif

         if(contran) then
            call epsilon_mult(ndsx,ndsy,ml,gcoef,acoeft,epsilon_cv,con_tran=contran)
            call gfunc_mult(ndsx,ndsy,ml,acoeft,g_coef=acoef,con_tran=contran)
         else
!            starttime=mstm_mpi_wtime()
            call gfunc_mult(ndsx,ndsy,ml,gcoef,g_coef=acoeft,con_tran=contran)
!            time1=mstm_mpi_wtime()
            call epsilon_mult(ndsx,ndsy,ml,acoeft,acoef,epsilon_cv,con_tran=contran)
!            time2=mstm_mpi_wtime()
!write(*,'(2e13.5)') (time1-starttime)/(time2-starttime),(time2-time1)/(time2-starttime)
         endif
         end subroutine interaction_matrix_mult

         subroutine incident_internal_field(z,pol,einc,m_0,m_slab,m_h)
         implicit none
         integer :: pol
         real(8) :: z,h,kr2,kr,kxh,kyh,rmat(2,2),ppar,pper
         complex(8) :: einc(3),ci,m,m0,ms,kz0,kzs,kz,ckzmh, &
            skzmh,f1den,f2den,ckzmhz,skzmhz,skz
         complex(8), optional :: m_0,m_slab,m_h
         data ci/(0.d0,1.d0)/
         ppar=dble(2-pol)
         pper=dble(pol-1)
         if(present(m_slab)) then
            m=m_slab
         else
            m=rimedium_cv
         endif
         if(present(m_0)) then
            m0=m_0
         else
            m0=ri0_cv
         endif
         if(present(m_h)) then
            ms=m_h
         else
            ms=ris_cv
         endif
         h=h_cv
         kr2=k0x_cv*k0x_cv+k0y_cv*k0y_cv
         kr=sqrt(kr2)
         if(kr2.eq.0.d0) then
            kxh=1.d0
            kyh=0.d0
         else
            kxh=k0x_cv/kr
            kyh=k0y_cv/kr
         endif
         rmat(1,:)=(/kxh,kyh/)
         rmat(2,:)=(/kyh,-kxh/)
         kz0=cdsqrt(1.d0-kr2/m0/m0)
         kzs=cdsqrt(1.d0-kr2/ms/ms)
         kz=cdsqrt(1.d0-kr2/m/m)
         skz=cdsqrt(1.d0-kz*kz)
         ckzmh=cdcos(kz*m*h)
         skzmh=cdsin(kz*m*h)
         ckzmhz=cdcos(kz*m*(h-z))
         skzmhz=cdsin(kz*m*(h-z))
         f1den=1.d0/(m*(ci*kz*m*(kzs*m0 + kz0*ms)*ckzmh &
            + (kz0*m*kzs*m + kz*m0*kz*ms)*skzmh))
         f2den=1.d0/(kz*m*(ci*kz*m*(kz0*m0 + kzs*ms)*ckzmh &
           + (kz*kz*m*m + kz0*kzs*m0*ms)*skzmh))
         einc(1)=2*f1den*kz*kz0*m*m0*(ci*kzs*m*ckzmhz + kz*ms*skzmhz)*ppar
         einc(2)=2*f2den*kz*kz0*m*m0*((-ci)*kz*m*ckzmhz - kzs*ms*skzmhz)*pper
         einc(3)=2*f1den*skz*kz0*m*m0*((-ci)*kz*ms*ckzmhz - kzs*m*skzmhz)*ppar
         einc(1:2)=matmul(rmat,einc(1:2))
         end subroutine incident_internal_field

         subroutine incident_external_field(eref,etra)
         implicit none
         real(8) :: h,kr2,kr,kxh,kyh,rmat(2,2),ppar,pper
         complex(8) :: eref(3),etra(3),ci,m,m0,ms,kz0,kzs,kz,ckzmh, &
            skzmh,f1den,f2den
         data ci/(0.d0,1.d0)/
         ppar=1
         pper=1
         m=rimedium_cv
         m0=ri0_cv
         ms=ris_cv
         h=h_cv
         kr2=k0x_cv*k0x_cv+k0y_cv*k0y_cv
         kr=sqrt(kr2)
         if(kr2.eq.0.d0) then
            kxh=1.d0
            kyh=0.d0
         else
            kxh=k0x_cv/kr
            kyh=k0y_cv/kr
         endif
         rmat(1,:)=(/kxh,kyh/)
         rmat(2,:)=(/kyh,-kxh/)
         kz0=cdsqrt(1.d0-kr2/m0/m0)
         kzs=cdsqrt(1.d0-kr2/ms/ms)
         kz=cdsqrt(1.d0-kr2/m/m)
         ckzmh=cdcos(kz*m*h)
         skzmh=cdsin(kz*m*h)
         f1den=1.d0/(m*(ci*kz*m*(kzs*m0 + kz0*ms)*ckzmh &
            + (kz0*m*kzs*m + kz*m0*kz*ms)*skzmh))
         f2den=1.d0/(kz*m*(ci*kz*m*(kz0*m0 + kzs*ms)*ckzmh &
           + (kz*kz*m*m + kz0*kzs*m0*ms)*skzmh))
         eref(1)=f1den*kz0*m*(ci*kz*m*(kzs*m0 - kz0*ms)*ckzmh &
            + (-(kz0*kzs*m*m) + kz*kz*m0*ms)*skzmh)
         eref(2)=f2den*kz*m*((-ci)*kz*m*(kz0*m0 - kzs*ms)*ckzmh &
            + (kz*kz*m*m - kz0*kzs*m0*ms)*skzmh)
         eref(3)=f1den*cdsqrt(1.d0 - kz0*kz0)*m*(ci*kz*m*(kzs*m0 - kz0*ms)*ckzmh &
            + (-(kz0*kzs*m*m) + kz*kz*m0*ms)*skzmh)
         etra(1)=2.d0*ci*f1den*kz*kz0*kzs*m*m*m0
         etra(2)=(-2.d0*ci)*f2den*kz*kz*kz0*m*m*m0
         etra(3)=(-2.d0*ci)*f1den*kz*kz0*cdsqrt(1.d0 - kzs*kzs)*m*m*m0
         end subroutine incident_external_field

         subroutine hemispherical_flux(ndsx,ndsy,ml,acoef, &
            hflux,cflux)
         implicit none
         integer :: ndsx,ndsy,ml,is1x,is2x,is1y,is2y,sx,sy,s
         real(8) :: hflux(2,2,2),kr2,kzinc,smat(4,4),kz0,kzs,cflux(2),hfluxt(2,2),kx,ky
         complex(8) :: acoef(3,ndsx*ndsy,ml,2),ampmat(2,2), &
            pcoef(2,ndsx*ndsy,2,2),eref(3),etra(3)

         kr2=k0x_cv*k0x_cv+k0y_cv*k0y_cv
         is1x=-ndsx/2
         is2x=(ndsx-1)/2
         is1y=-ndsy/2
         is2y=(ndsy-1)/2
         kzinc=sqrt(1.d0-kr2/(ri0_cv*ri0_cv))
         hflux=0.
         call gfunc_mult(ndsx,ndsy,ml,acoef(:,:,:,1), &
             p_coef=pcoef(:,:,:,1))
         call gfunc_mult(ndsx,ndsy,ml,acoef(:,:,:,2), &
             p_coef=pcoef(:,:,:,2))
         call incident_external_field(eref,etra)

         do sy=is1y,is2y
            ky=k0y_cv+2.d0*pi*dble(sy)/wy_cv
            do sx=is1x,is2x
               kx=k0x_cv+2.d0*pi*dble(sx)/wx_cv
               kr2=kx*kx+ky*ky
               s=sx-is1x+1+(sy-is1y)*ndsx
               if(kr2/dble(ri0_cv*ri0_cv).lt.1.d0) then
                  ampmat(:,:)=pcoef(:,s,1,:)
                  kz0=sqrt(1.d0-kr2/dble(ri0_cv*ri0_cv))
                  if(sx.eq.0.and.sy.eq.0) then
                     ampmat(1,1)=ampmat(1,1)+eref(1)/kz0
                     ampmat(2,2)=ampmat(2,2)+eref(2)
                  endif
                  call amplitude_to_scatteringmatrix(ampmat(:,:), &
                     smat(:,:))
!                  hfluxt(1,1)=(smat(1,1)+smat(1,2) &
!                     +smat(2,1)+smat(2,2))*ri0_cv*kz0/kzinc/2.d0
!                  hfluxt(2,1)=(smat(1,1)+smat(1,2) &
!                     -smat(2,1)-smat(2,2))*ri0_cv*kz0/kzinc/2.d0
!                  hfluxt(2,2)=(smat(1,1)-smat(1,2) &
!                     -smat(2,1)+smat(2,2))*ri0_cv*kz0/kzinc/2.d0
!                  hfluxt(1,2)=(smat(1,1)-smat(1,2) &
!                     +smat(2,1)-smat(2,2))*ri0_cv*kz0/kzinc/2.d0
                  hfluxt(1,1)=(smat(1,1)+smat(1,2) &
                     +smat(2,1)+smat(2,2))*kz0/kzinc/2.d0
                  hfluxt(2,1)=(smat(1,1)+smat(1,2) &
                     -smat(2,1)-smat(2,2))*kz0/kzinc/2.d0
                  hfluxt(2,2)=(smat(1,1)-smat(1,2) &
                     -smat(2,1)+smat(2,2))*kz0/kzinc/2.d0
                  hfluxt(1,2)=(smat(1,1)-smat(1,2) &
                     +smat(2,1)-smat(2,2))*kz0/kzinc/2.d0
                  hflux(:,:,1)=hflux(:,:,1)+hfluxt
                  if(sx.eq.0.and.sy.eq.0) then
                     cflux(1)=(hfluxt(1,1)+hfluxt(1,2)+hfluxt(2,1)+hfluxt(2,2))/2.d0
                  endif
               endif
               if(kr2/dble(ris_cv*ris_cv).lt.1.d0.and.dimag(ris_cv).eq.0.d0) then
                  ampmat(:,:)=pcoef(:,s,2,:)
                  kzs=sqrt(1.d0-kr2/dble(ris_cv*ris_cv))
                  if(sx.eq.0.and.sy.eq.0) then
                     ampmat(1,1)=ampmat(1,1)+etra(1)/kzs
                     ampmat(2,2)=ampmat(2,2)+etra(2)
                  endif
                  call amplitude_to_scatteringmatrix(ampmat(:,:), &
                      smat(:,:))
                  hfluxt(1,1)=(smat(1,1)+smat(1,2) &
                     +smat(2,1)+smat(2,2))*ris_cv/ri0_cv*kzs/kzinc/2.d0
                  hfluxt(2,1)=(smat(1,1)+smat(1,2) &
                     -smat(2,1)-smat(2,2))*ris_cv/ri0_cv*kzs/kzinc/2.d0
                  hfluxt(1,2)=(smat(1,1)-smat(1,2) &
                     +smat(2,1)-smat(2,2))*ris_cv/ri0_cv*kzs/kzinc/2.d0
                  hfluxt(2,2)=(smat(1,1)-smat(1,2) &
                     -smat(2,1)+smat(2,2))*ris_cv/ri0_cv*kzs/kzinc/2.d0
!                  hfluxt(1,1)=(smat(1,1)+smat(1,2) &
!                     +smat(2,1)+smat(2,2))*kzs/kzinc/2.d0
!                  hfluxt(2,1)=(smat(1,1)+smat(1,2) &
!                     -smat(2,1)-smat(2,2))*kzs/kzinc/2.d0
!                  hfluxt(1,2)=(smat(1,1)-smat(1,2) &
!                     +smat(2,1)-smat(2,2))*kzs/kzinc/2.d0
!                  hfluxt(2,2)=(smat(1,1)-smat(1,2) &
!                     -smat(2,1)+smat(2,2))*kzs/kzinc/2.d0
                  hflux(:,:,2)=hflux(:,:,2)+hfluxt
                  if(sx.eq.0.and.sy.eq.0) then
                     cflux(2)=(hfluxt(1,1)+hfluxt(1,2)+hfluxt(2,1)+hfluxt(2,2))/2.d0
                  endif
               endif
            enddo
         enddo
         end subroutine hemispherical_flux

         subroutine absorption_distribution(ndsx,ndsy,ml,acoef,qabsj)
         implicit none
         integer :: ndsx,ndsy,ml,n,j
         real(8) :: z,qabsj(ml)
         complex(8) :: acoef(3,-ndsx/2:(ndsx-1)/2,-ndsy/2:(ndsy-1)/2,ml,2), &
            pcoef(3,-ndsx/2:(ndsx-1)/2,-ndsy/2:(ndsy-1)/2,ml,2), &
            ac(3,-ndsx/2:(ndsx-1)/2,-ndsy/2:(ndsy-1)/2,ml,2), &
            einc(3,2), &
            a(3,-ndsx/2:(ndsx-1)/2,-ndsy/2:(ndsy-1)/2,2), &
            p(3,-ndsx/2:(ndsx-1)/2,-ndsy/2:(ndsy-1)/2,2)

         call gfunc_mult(ndsx,ndsy,ml,acoef(:,:,:,:,1),g_coef=pcoef(:,:,:,:,1))
         call gfunc_mult(ndsx,ndsy,ml,acoef(:,:,:,:,2),g_coef=pcoef(:,:,:,:,2))

         do j=1,ml
            z=(dble(j)-0.5d0)*d_cv
            call incident_internal_field(z,1,einc(:,1))
            call incident_internal_field(z,2,einc(:,2))
            pcoef(:,0,0,j,1)=pcoef(:,0,0,j,1)+einc(:,1)
            pcoef(:,0,0,j,2)=pcoef(:,0,0,j,2)+einc(:,2)
         enddo

         epsilon_cv=epsilon_cv+rimedium_cv*rimedium_cv

         call epsilon_mult(ndsx,ndsy,ml,pcoef(:,:,:,:,1), &
             ac(:,:,:,:,1),epsilon_cv)
         call epsilon_mult(ndsx,ndsy,ml,pcoef(:,:,:,:,2), &
             ac(:,:,:,:,2),epsilon_cv)

         epsilon_cv=epsilon_cv-rimedium_cv*rimedium_cv

         n=3*ndsx*ndsy*2
         do j=1,ml
            a=ac(:,:,:,j,:)
            p=pcoef(:,:,:,j,:)
            qabsj(j)=-dimag(dot_product( &
               reshape(a(:,:,:,:),(/n/)), &
               reshape(p(:,:,:,:),(/n/)) &
               ))
         enddo
         end subroutine absorption_distribution

         subroutine scattered_field(ndsx,ndsy,ml,nx1,nx2,ny1,ny2, &
            acoef,smat,incident_rotate,alpha_rotate, &
            include_incident)
         implicit none
         logical :: incrot,incinc
         logical, optional :: incident_rotate,include_incident
         integer :: ndsx,ndsy,ml,is1x,is2x,is1y,is2y,s,sy,sx,ss0,nx1,nx2,ny1,ny2
         real(8) :: kr2,kzinc,kx,ky,kr,kxh,kyh,h, &
            smat(4,4,2,nx1:nx2,ny1:ny2),alpha,cp,sp
         real(8), optional :: alpha_rotate
         complex(8) :: acoef(3,ndsx*ndsy,ml,2),m,m0,ms,kz0,kzs, &
             eref(3),etra(3),ci, &
            ealpha,ephi,ampmat(2,2), &
            pcoef(2,ndsx*ndsy,2,2)
         data ci/(0.d0,1.d0)/
         if(present(incident_rotate).and.k0x_cv.eq.0.d0.and.k0y_cv.eq.0.d0) then
            incrot=incident_rotate
         else
            incrot=.false.
         endif
         if(present(alpha_rotate)) then
            alpha=alpha_rotate
         else
            alpha=0.d0
         endif
         if(present(include_incident)) then
            incinc=include_incident
         else
            incinc=.false.
         endif
         ealpha=cdexp(-ci*alpha)
         is1x=-ndsx/2
         is2x=(ndsx-1)/2
         is1y=-ndsy/2
         is2y=(ndsy-1)/2
         ss0=-is1x+1+ndsx*(-is1y)
         m=rimedium_cv
         m0=ri0_cv
         ms=ris_cv
         h=h_cv
         kr2=k0x_cv*k0x_cv+k0y_cv*k0y_cv
         kzinc=sqrt(1.d0-kr2/m0/m0)
         call gfunc_mult(ndsx,ndsy,ml,acoef(:,:,:,1), &
             p_coef=pcoef(:,:,:,1))
         call gfunc_mult(ndsx,ndsy,ml,acoef(:,:,:,2), &
             p_coef=pcoef(:,:,:,2))
         call incident_external_field(eref,etra)
         smat=0.d0
         do sy=ny1,ny2
            ky=k0y_cv+2.d0*pi*dble(sy)/wy_cv
            do sx=nx1,nx2
               kx=k0x_cv+2.d0*pi*dble(sx)/wx_cv
               s=sx-is1x+1+(sy-is1y)*ndsx
               kr2=kx*kx+ky*ky
               kr=sqrt(kr2)
               if(kr2.eq.0.d0) then
                  kxh=1.d0
                  kyh=0.d0
               else
                  kxh=kx/kr
                  kyh=ky/kr
               endif
               ephi=dcmplx(kxh,kyh)
               if(incrot) then
                  cp=dble(ealpha*ephi)
                  sp=dimag(ealpha*ephi)
               else
                  cp=1.d0
                  sp=0.d0
               endif
               if(kr2/dble(ri0_cv*ri0_cv).lt.1.d0) then
                  ampmat(:,1)=pcoef(:,s,1,1)*cp+pcoef(:,s,1,2)*sp
                  ampmat(:,2)=-pcoef(:,s,1,1)*sp+pcoef(:,s,1,2)*cp
                  kz0=sqrt(1.d0-kr2/dble(ri0_cv*ri0_cv))
                  if(sx.eq.0.and.sy.eq.0.and.incinc) then
                     ampmat(1,1)=ampmat(1,1)+eref(1)/kz0
                     ampmat(2,2)=ampmat(2,2)+eref(2)
                  endif
                  call amplitude_to_scatteringmatrix(ampmat(:,:), &
                     smat(:,:,1,sx,sy))
                  smat(:,:,1,sx,sy)=smat(:,:,1,sx,sy)*kz0/kzinc
               endif
               if(kr2/dble(ris_cv*ris_cv).lt.1.d0 &
                     .and.dimag(ris_cv).eq.0.d0) then
                  ampmat(:,1)=pcoef(:,s,2,1)*cp+pcoef(:,s,2,2)*sp
                  ampmat(:,2)=-pcoef(:,s,2,1)*sp+pcoef(:,s,2,2)*cp
                  kzs=sqrt(1.d0-kr2/dble(ris_cv*ris_cv))
                  if(sx.eq.0.and.sy.eq.0.and.incinc) then
                     ampmat(1,1)=ampmat(1,1)+etra(1)/kzs
                     ampmat(2,2)=ampmat(2,2)+etra(2)
                  endif
                  call amplitude_to_scatteringmatrix(ampmat(:,:), &
                     smat(:,:,2,sx,sy))
                  smat(:,:,2,sx,sy)= &
                  smat(:,:,2,sx,sy)*kzs/kzinc*dble(ris_cv)
               endif
            enddo
         enddo
         end subroutine scattered_field

         subroutine amplitude_to_scatteringmatrix(samat,sm,i_max)
         implicit none
         integer :: i,j,imax
         integer, optional :: i_max
         real(8) :: sm(4,4)
         complex(8) :: sa(4),sp(4,4),samat(2,2)
         sa(1)=samat(2,2)
         sa(2)=samat(1,1)
         sa(3)=samat(1,2)
         sa(4)=samat(2,1)
         if(present(i_max)) then
            imax=i_max
         else
            imax=4
         endif
         if(imax.eq.1) then
            sm(1,1)=0.d0
            do i=1,4
               sm(1,1)=sm(1,1)+sa(i)*dconjg(sa(i))
            enddo
            sm(1,1)=sm(1,1)*16.d0
            return
         endif
         do i=1,4
            do j=1,4
               sp(i,j)=sa(i)*dconjg(sa(j))/2.d0
            enddo
         enddo
         sm(1,1)=sp(1,1)+sp(2,2)+sp(3,3)+sp(4,4)
         sm(1,2)=-sp(1,1)+sp(2,2)-sp(3,3)+sp(4,4)
         sm(2,1)=-sp(1,1)+sp(2,2)+sp(3,3)-sp(4,4)
         sm(2,2)=sp(1,1)+sp(2,2)-sp(3,3)-sp(4,4)
         sm(3,3)=2.*(sp(1,2)+sp(3,4))
         sm(3,4)=2.*dimag(sp(2,1)+sp(4,3))
         sm(4,3)=2.*dimag(sp(1,2)-sp(3,4))
         sm(4,4)=2.*(sp(1,2)-sp(3,4))
         sm(1,3)=2.*(sp(2,3)+sp(1,4))
         sm(3,1)=2.*(sp(2,4)+sp(1,3))
         sm(1,4)=-2.*dimag(sp(2,3)-sp(1,4))
         sm(4,1)=-2.*dimag(sp(4,2)+sp(1,3))
         sm(2,3)=-2.*(sp(2,3)-sp(1,4))
         sm(3,2)=-2.*(sp(2,4)-sp(1,3))
         sm(2,4)=-2.*dimag(sp(2,3)+sp(1,4))
         sm(4,2)=-2.*dimag(sp(4,2)-sp(1,3))
         end subroutine amplitude_to_scatteringmatrix

         subroutine ppgf_solver(ndsx,ndsy,ml,niter,eps,acoef,iter, &
            errmax,mpi_comm)
         use mpidefs
         implicit none
         logical :: firstrun
         integer, save :: p1,p2,pcomm,sendrank
         integer :: ndsx,ndsy,ml,niter,iter,p,mpicomm,color, &
                   i1x,i2x,i1y,i2y,ss0,numprocs,rank,j,nsend,nds2,iterwrite
         integer, optional :: mpi_comm
         real(8) :: eps,errmax,z
         complex(8) :: acoef(3,ndsx*ndsy,ml,2),gcoef(3,ndsx*ndsy,ml)
         data firstrun/.true./
         data iterwrite/1/
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size', mpi_size=numprocs, &
              mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank', mpi_rank=rank, &
              mpi_comm=mpicomm)
         if(firstrun) then
            firstrun=.false.
            if(numprocs.eq.4) then
               color=floor(dble(2*rank)/dble(numprocs))+1
               p1=color
               p2=p1
               sendrank=2
               call mstm_mpi(mpi_command='split', &
                  mpi_color=color,mpi_key=rank, &
                  mpi_new_comm=pcomm, &
                  mpi_comm=mpicomm)
            elseif(numprocs.eq.2) then
               color=rank+1
               p1=color
               p2=p1
               sendrank=1
               call mstm_mpi(mpi_command='split', &
                  mpi_color=color,mpi_key=rank, &
                  mpi_new_comm=pcomm, &
                  mpi_comm=mpicomm)
            else
               p1=1
               p2=2
               pcomm=mpicomm
            endif
         endif
         nds2=ndsx*ndsy
         i1x=-ndsx/2
         i2x=(ndsx-1)/2
         i1y=-ndsy/2
         i2y=(ndsy-1)/2
         ss0=-i1x+1+ndsx*(-i1y)
         acoef=0.d0

         do p=p1,p2
            do j=1,ml
               z=d_cv*(dble(j)-0.5d0)
               call incident_internal_field(z,p,acoef(:,ss0,j,p))
            enddo
            call epsilon_mult(ndsx,ndsy,ml,acoef(:,:,:,p),gcoef,epsilon_cv)
            acoef(:,:,:,p)=gcoef
            iter=0
            call cbicg(ndsx,ndsy,ml,niter,eps,gcoef,acoef(:,:,:,p), &
                 iterwrite,iter,errmax,mpi_comm=pcomm)
         enddo
         if(numprocs.gt.1) then
            nsend=3*nds2*ml
            call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=acoef(1:3,1:nds2,1:ml,2), &
                  mpi_number=nsend, &
                  mpi_rank=sendrank, &
                  mpi_comm=mpicomm)
         endif
         end subroutine ppgf_solver

         subroutine cbicg(ndsx,ndsy,ml,niter,eps,pnp,anp,iterwrite,iter, &
                 errmax,mpi_comm)
         use mpidefs
         implicit none
         logical, save :: firstrun
         integer :: ndsx,ndsy,ml,neqns,niter,iter,writetime,iterwrite, &
                    ntot,rank,icon(1),mpicomm,numprocs,rank0
         integer, optional :: mpi_comm
         real(8) :: ptime(5),time1,time2,starttime
         real(8) :: eps,errmax,enorm,minerr
         complex(8)  :: pnp(3*ndsx*ndsy*ml),anp(3*ndsx*ndsy*ml),cak,csk,cbk,csk2
         complex(8) :: cr(3*ndsx*ndsy*ml),cp(3*ndsx*ndsy*ml),cw(3*ndsx*ndsy*ml), &
                       cq(3*ndsx*ndsy*ml),cap(3*ndsx*ndsy*ml),caw(3*ndsx*ndsy*ml),capt(3*ndsx*ndsy*ml), &
                       cawt(3*ndsx*ndsy*ml)
         data writetime,firstrun/0,.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         neqns=3*ndsx*ndsy*ml
         ntot=neqns
         iter=0
         errmax=0.
         minerr=1.d10
         ptime=0.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
!if(rank.eq.0) then
!write(*,'('' rank,numprocs:'',2i5)') rank0,numprocs
!endif

!
! the following is the implementation of the complex biconjugate gradient
! iteration scheme.   Thank you P Flatau.
!
         if(niter.eq.0) then
            anp=pnp
            return
         endif
         enorm=dble(dot_product(pnp(:),pnp(:)))

         if(niter.lt.0) then
            if(rank.eq.0) then
               starttime=mstm_mpi_wtime()
               anp=pnp
               cr=anp
               errmax=1.d0
               do iter=1,abs(niter)
                  time1=mstm_mpi_wtime()
                  call interaction_matrix_mult(ndsx,ndsy,ml,cr,cr)
                  anp=anp+cr
                  errmax=dble(dot_product(cr(:),cr(:)))/enorm
                  if((iterwrite.eq.1.or.firstrun).and.rank0.eq.0) then
                     time2=mstm_mpi_wtime()
                     if(firstrun.or.(time2-starttime.gt.5.d0).or.every_iteration_cv) then
                        write(run_print_unit_cv,'('' iter,err,tpi:'',i5,e12.4,2e12.4)') &
                           iter,errmax,time2-time1
                        call flush(run_print_unit_cv)
                        starttime=time2
                     endif
                     firstrun=.false.
                  endif
               if(errmax.lt.eps) exit
               enddo
            endif
            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
            return
         endif

         cr=0.d0
         icon=0
         if(rank.eq.0) then
            call interaction_matrix_mult(ndsx,ndsy,ml,anp,cr)
         endif

         if(rank.eq.0) then
            cr(:)=pnp(:)-anp(:)+cr(:)
            cq(:)=conjg(cr(:))
            cw(:)=cq(:)
            cp(:)=cr(:)
            csk=dot_product(conjg(cr(:)),cr(:))
            if(cdabs(csk).eq.0.d0) icon=1
         endif
         call mstm_mpi(mpi_command='bcast', &
             mpi_send_buf_i=icon, &
             mpi_number=1, &
             mpi_rank=0, &
             mpi_comm=mpicomm)
         errmax=dot_product(cr(:),cr(:))/enorm
         if(icon(1).eq.1) return
!
!  here starts the main iteration loop
!
         if(rank.eq.0) starttime=mstm_mpi_wtime()
         do iter=1,niter
            if(rank.eq.0) time1=mstm_mpi_wtime()
            cak=(0.d0,0.d0)
            cawt(:)=(0.d0,0.d0)
            capt(:)=(0.d0,0.d0)
            cap(:)=0.d0
            caw(:)=0.d0
            ptime=0.
            if(numprocs.eq.1) then
               call interaction_matrix_mult(ndsx,ndsy,ml,cp,cap)
               call interaction_matrix_mult(ndsx,ndsy,ml,cw,caw, &
                    con_tran=.true.)
            else
               call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=cw, &
                  mpi_number=neqns, &
                  mpi_rank=0, &
                  mpi_comm=mpicomm)
               if(rank.eq.0) then
                  call interaction_matrix_mult(ndsx,ndsy,ml,cp,cap)
               elseif(rank.eq.1) then
                  call interaction_matrix_mult(ndsx,ndsy,ml,cw,caw, &
                       con_tran=.true.)
               endif
               call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=caw, &
                  mpi_number=neqns, &
                  mpi_rank=1, &
                  mpi_comm=mpicomm)
            endif
            if(rank.eq.0) then
               cap(:)=cp(:)-cap(:)
               caw(:)=cw(:)-caw(:)
               cak=dot_product(cw(:),cap(:))
               cak=csk/cak
               anp(:)=anp(:)+cak*cp(:)
               cr(:)=cr(:)-cak*cap(:)
               cq(:)=cq(:)-conjg(cak)*caw(:)
               csk2=dot_product(cq(:),cr(:))
               errmax=dot_product(cr(:),cr(:))
               errmax=errmax/enorm
               minerr=min(minerr,errmax)
               if(errmax.lt.eps.or.cdabs(csk).eq.0.d0) then
                  icon=1
               else
                  cbk=csk2/csk
                  csk=csk2
               endif
               cp(:)=cr(:)+cbk*cp(:)
               cw(:)=cq(:)+conjg(cbk)*cw(:)
               if(((iterwrite.eq.1).or.firstrun.or.every_iteration_cv).and.rank0.eq.0) then
                  time2=mstm_mpi_wtime()
                  if(firstrun.or.(time2-starttime.gt.5.d0).or.every_iteration_cv) then
                     write(run_print_unit_cv,'('' iter,err,min err, tpi:'',i5,2e12.4,e12.4)') &
                        iter,errmax,minerr,time2-time1
                     call flush(run_print_unit_cv)
                     starttime=time2
                  endif
                  firstrun=.false.
               endif
            endif
            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='bcast', &
                mpi_send_buf_i=icon, &
                mpi_number=1, &
                mpi_rank=0, &
                mpi_comm=mpicomm)
            if(icon(1).eq.1) return
         enddo
         end subroutine cbicg

      end module solver
