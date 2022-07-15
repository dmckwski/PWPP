      module target_generation

      type element_list
      integer(2) :: pos(4)
      type(element_list), pointer :: next=>null()
      end type element_list

      type surf_element_list
      integer(2) :: pos(3)
      type(surf_element_list), pointer :: next=>null()
      end type surf_element_list

      type aggregate_type
      integer :: index,size,surf_size,mean_position(3),box(3,2)
      type(element_list), pointer :: elements
      type(surf_element_list), pointer :: surf_elements
      type(aggregate_type), pointer :: next=>null()
      end type aggregate_type

      type(aggregate_type), pointer :: aggregate_data(:)
      integer, allocatable :: aggregate_number(:),root_monomer(:), &
      aggregate_position(:,:),aggregate_size(:)

      contains

         subroutine generate_random_particle_target(ncomp,meanradius,fv,pshape,sigma, &
            dz,dlat,cellmap,cld,cud,npartstart,npart,iflag,aspect_ratio_y,aspect_ratio_z,boundary_model, &
            periodic_boundary,diffusion_steps,max_number_particles, &
            print_unit,sticking_prob,rotation_prob,write_frames,aggregation_frac, &
            shell_thickness)
         use specialfuncs
         use mpidefs
         implicit none
         type(element_list), pointer :: elist,elisttemp
         type(surf_element_list), pointer :: selist,selisttemp
         logical :: noparticleoverlap,particleincell,pb(3)
         logical, optional :: periodic_boundary(3)
         integer :: ncomp,comp,j,i,nparti(ncomp),pshape(ncomp),npart,nd, &
            cld(3),cud(3),izmin,izmax,rsamp(3),maxsamp,boundarymodel,npartleft, &
            fitsamp,izsampmin,izsampmax,ndtot,dsteps,dstep,cellnd, &
            npartstart,iflag,punit,ix,iy,iz,rank0,ndsurf,ndsurftot,writeframes,npart0
         integer, allocatable :: compi(:)
         integer, optional :: boundary_model,diffusion_steps, &
            max_number_particles,print_unit,write_frames
         integer(4) :: cellmap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3))
         integer(2), allocatable :: elems(:,:),selems(:,:)
         real(8) :: meanradius(ncomp),fv(ncomp),sigma(ncomp),ary(ncomp),arz(ncomp), &
                    maxradius,sampradius,dz,dlat,celldim(1:3),shellthick(ncomp), &
                    cellscale(3),ptime(5),stickprob,rotprob,aggfrac,aggfrac0
         real(8), allocatable :: rad(:),eulerang(:,:)
         real(8), optional :: aspect_ratio_y(ncomp),aspect_ratio_z(ncomp),sticking_prob, &
                  rotation_prob,aggregation_frac,shell_thickness(ncomp)
         data maxradius,maxsamp/2.d0,5000/
         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
         iflag=0
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
         if(present(periodic_boundary)) then
            pb=periodic_boundary
         else
            pb=(/.true.,.true.,.false./)
         endif
         if(present(boundary_model)) then
            boundarymodel=boundary_model
         else
            boundarymodel=0
         endif
         if(present(diffusion_steps)) then
            dsteps=diffusion_steps
         else
            dsteps=0
         endif
         if(present(print_unit)) then
            punit=print_unit
         else
            punit=6
         endif
         if(present(sticking_prob)) then
            stickprob=sticking_prob
         else
            stickprob=1.
         endif
         if(present(rotation_prob)) then
            rotprob=rotation_prob
         else
            rotprob=0.d0
         endif
         if(present(write_frames)) then
            writeframes=write_frames
         else
            writeframes=0
         endif
         if(present(aggregation_frac)) then
            aggfrac0=aggregation_frac
         else
            aggfrac0=1.d0
         endif
         if(present(shell_thickness)) then
            shellthick=shell_thickness
         else
            shellthick=0.d0
         endif

         cellscale(1:3)=(/dlat,dlat,dz/)
         celldim(1:3)=dble(cud(1:3)-cld(1:3)+1)*cellscale(1:3)
         cellnd=(cud(1)-cld(1)+1)*(cud(2)-cld(2)+1)*(cud(3)-cld(3)+1)
         nparti(1:ncomp)=nint(fv(1:ncomp)*celldim(1)*celldim(2)*celldim(3) &
             *3./(4.d0*3.1415927*(meanradius(1:ncomp)+shellthick(1:ncomp))**3.))
         if(present(max_number_particles)) then
            npart=sum(nparti(1:ncomp))
            if(npart.gt.max_number_particles) then
               do i=1,ncomp
                  nparti(i)=nint((dble(nparti(i))/dble(npart))*max_number_particles)
               enddo
            endif
         endif
         npart=sum(nparti(1:ncomp))

         allocate(rad(npart),compi(npart),eulerang(3,npart))

         i=0
         do comp=1,ncomp
            do j=1,nparti(comp)
               call psdsamp(sigma(comp),maxradius,sampradius)
               sampradius=sampradius*meanradius(comp)
               i=i+1
               rad(i)=sampradius
               compi(i)=comp
            enddo
         enddo

         do i=1,npart
            do j=1,npart-i
               if(rad(j).lt.rad(j+1)) then
                  sampradius=rad(j)
                  rad(j)=rad(j+1)
                  rad(j+1)=sampradius
                  comp=compi(j)
                  compi(j)=compi(j+1)
                  compi(j+1)=comp
               endif
            enddo
         enddo

         allocate(aggregate_data(npart),aggregate_number(npart), &
            root_monomer(npart),aggregate_position(3,npart), &
            aggregate_size(npart))

         npartleft=npart
         ndtot=0
         ndsurftot=0
         do i=1,npart
!            isamp=i+int(ran3(1)*dble(npartleft))
!            isamp=max(i,isamp)
!            isamp=min(isamp,npart)
!            comp=compi(isamp)
!            sampradius=rad(isamp)
!            trad=rad(i)
!            rad(i)=rad(isamp)
!            rad(isamp)=trad
!            j=compi(i)
!            compi(i)=compi(isamp)
!            compi(isamp)=j
!            npartleft=npartleft-1
            comp=compi(i)
            sampradius=rad(i)
            eulerang(1,i)=2.d0*3.1415927*ran3(1)
            eulerang(2,i)=dacos(-1.d0+2.*ran3(1))
            eulerang(3,i)=2.d0*3.1415927*ran3(1)
            call generate_single_particle_list(sampradius,pshape(comp),ary(comp), &
               arz(comp),eulerang(:,i),dz,dlat,nd,ndsurf,shell_thick=shellthick(comp))
            aggregate_data(i)%index=i
            aggregate_data(i)%size=nd
            aggregate_data(i)%surf_size=ndsurf
            aggregate_data(i)%mean_position=(/0,0,0/)
            allocate(elems(4,nd),selems(3,ndsurf))
            call generate_single_particle_list(sampradius,pshape(comp),ary(comp), &
               arz(comp),eulerang(:,i),dz,dlat,nd,ndsurf,shell_thick=shellthick(comp), &
               particle_dlist=elems(1:4,1:nd), &
               particle_surf_dlist=selems(1:3,1:ndsurf) &
            )
            do j=1,3
               aggregate_data(i)%box(j,:)=(/minval(elems(j,1:nd)),maxval(elems(j,1:nd))/)
            enddo
            allocate(aggregate_data(i)%elements)
            elist=>aggregate_data(i)%elements
            do j=1,nd
               elist%pos(:)=elems(:,j)
               if(j.lt.nd) then
                  allocate(elist%next)
                  elist=>elist%next
               else
                  elist%next=>null()
               endif
            enddo
            allocate(aggregate_data(i)%surf_elements)
            selist=>aggregate_data(i)%surf_elements
            do j=1,ndsurf
               selist%pos(:)=selems(:,j)
               if(j.lt.ndsurf) then
                  allocate(selist%next)
                  selist=>selist%next
               else
                  selist%next=>null()
               endif
            enddo
            aggregate_number(i)=1
            root_monomer(i)=i
            aggregate_position(:,i)=(/0,0,0/)
            aggregate_size(i)=nd
            deallocate(elems,selems)
            ndtot=nd+ndtot
            ndsurftot=ndsurf+ndsurftot
         enddo
         if(rank0.eq.0) then
            write(punit,'('' initial npart, volume fraction:'',i6,e13.5)') npart,dble(ndtot)/dble(cellnd)
            call flush(punit)
         endif

         cellmap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3))=0
         do i=1,npart
            nd=aggregate_data(i)%size
            izmin=aggregate_data(i)%box(3,1)
            izmax=aggregate_data(i)%box(3,2)
            if(boundarymodel.eq.0) then
               izsampmin=cld(3)-izmin
               izsampmax=cud(3)-izmax
            else
               izsampmin=cld(3)-izmax
               izsampmax=cud(3)-izmin
            endif
            do fitsamp=1,maxsamp
               rsamp(1)=floor(dble(cud(1)-cld(1)+1)*ran3(1))+cld(1)
               rsamp(2)=floor(dble(cud(2)-cld(2)+1)*ran3(1))+cld(2)
               rsamp(3)=floor(dble(izsampmax-izsampmin+1)*ran3(1))+izsampmin
               aggregate_position(:,i)=rsamp(:)
               call particle_uc_fit_operation(i, &
                  cellmap,cld,cud,noparticleoverlap,particleincell, &
                  periodic_boundary=pb)
               if((boundarymodel.eq.0).and.(.not.particleincell)) cycle
               if(.not.noparticleoverlap) cycle
               exit
            enddo
            if((.not.noparticleoverlap).or.(.not.particleincell)) then
               if(iflag.eq.0) iflag=1
               write(punit,'('' no fit'',2i6,2l1)') i,fitsamp,noparticleoverlap,particleincell
               call flush(punit)
            endif
            call particle_uc_fit_operation(i, &
               cellmap,cld,cud,noparticleoverlap,particleincell, &
               set_uc_map=.true.,set_uc_map_index=i,periodic_boundary=pb)
         enddo

         npartleft=npart
         npartstart=npart
         if(dsteps.gt.0.and.rank0.eq.0) then
            write(punit,'('' initial number of particles:'',i6)') npart
            call flush(punit)
         endif
         if(writeframes.gt.0.and.rank0.eq.0) then
            open(20,file='atest.dat')
!            call write_cell_map(cellmap,cld,cud,20)
            call write_monomer_positions(npartstart,20)
         endif

         call unit_cell_aggregation(npartstart,npart,cellmap,cld,cud,npartleft, &
            periodic_boundary=pb, &
            diffusion_simulation=.false., &
            process_timings=ptime, &
            sticking_prob=stickprob, &
            rotation_prob=rotprob)
         npart0=npartleft
         npart=npart0

         do dstep=1,dsteps
            npart=npartleft
            call unit_cell_aggregation(npartstart,npart,cellmap,cld,cud,npartleft, &
               periodic_boundary=pb, &
               diffusion_simulation=.true., &
               process_timings=ptime, &
               sticking_prob=stickprob, &
               rotation_prob=rotprob)
            aggfrac=1.d0-dble(npartleft)/dble(npart0)
            if(npart.ne.npartleft.and.rank0.eq.0) then
               ptime(3)=ptime(1)+ptime(2)
               ptime=ptime/ptime(3)
               write(punit,'('' step, number particles, aggregation fraction:'',2i6,2f8.3)') dstep,npartleft,aggfrac
               call flush(punit)
            endif
            if(writeframes.eq.2.and.rank0.eq.0) then
!               call write_cell_map(cellmap,cld,cud,20)
               call write_monomer_positions(npartstart,20)
            endif
            if(npartleft.eq.1) exit
            if(aggfrac.ge.aggfrac0) exit
         enddo

         if(writeframes.eq.1) then
!           call write_cell_map(cellmap,cld,cud,20)
            call write_monomer_positions(npartstart,20)
         endif
         if(writeframes.gt.0) close(20)

         ndtot=0
         do iz=cld(3),cud(3)
            do iy=cld(2),cud(2)
               do ix=cld(1),cud(1)
                  i=cellmap(ix,iy,iz)
                  if(i.ne.0) then
                     ndtot=ndtot+1
                     if(i.gt.0) then
                        cellmap(ix,iy,iz)=compi(i)
                     elseif(i.lt.0) then
                        cellmap(ix,iy,iz)=-compi(-i)
                     endif
                  endif
               enddo
            enddo
         enddo
         if(rank0.eq.0) then
            write(punit,'('' final npart, volume fraction:'',i6,e13.5)') npart,dble(ndtot)/dble(cellnd)
            call flush(punit)
         endif

         do i=1,npartstart
            elist=>aggregate_data(i)%elements
            do while(associated(elist))
               elisttemp=>elist%next
               deallocate(elist)
               nullify(elist)
               elist=>elisttemp
            enddo
            selist=>aggregate_data(i)%surf_elements
            do while(associated(selist))
               selisttemp=>selist%next
               deallocate(selist)
               nullify(selist)
               selist=>selisttemp
            enddo
         enddo

         deallocate(rad,compi,eulerang,root_monomer,aggregate_data, &
            aggregate_number,aggregate_position,aggregate_size)
         end subroutine generate_random_particle_target

         subroutine nearest_neighbor_test(ipart, &
            cellmap,cld,cud,particleincell,nnn,nnarray, &
            periodic_boundary,nn_test,particle_displacement, &
            rotation_matrix)
         implicit none
         type(surf_element_list), pointer :: selist
         type(aggregate_type), pointer :: alist
         logical :: particleincell,pincell,pb(3),nntest,dup
         logical, optional :: periodic_boundary(3),nn_test
         integer :: nd,rpos(3),i,j,r(3),cld(3),cud(3), &
                    nagg,ia,testid,ipart,dpos(3), &
                    nnn,nnarray(*),k
         integer, optional :: particle_displacement(3)
         integer(4) :: cellmap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3))
         integer, optional :: rotation_matrix(3,3)
         if(present(periodic_boundary)) then
            pb=periodic_boundary
         else
            pb=(/.true.,.true.,.false./)
         endif
         if(present(nn_test)) then
            nntest=nn_test
         else
            nntest=.true.
         endif
         if(present(particle_displacement)) then
            dpos(1:3)=particle_displacement(1:3)+aggregate_position(1:3,ipart)
         else
            dpos=aggregate_position(1:3,ipart)
         endif

         particleincell=.true.
         nagg=aggregate_number(ipart)
         nnn=0
         alist=>aggregate_data(ipart)
         do ia=1,nagg
            if(.not.associated(alist)) then
               write(*,'('' error'')')
               stop
            endif
            rpos(1:3)=alist%mean_position(1:3)
            nd=alist%surf_size
            selist=>alist%surf_elements
            if(associated(alist%next)) alist=>alist%next
            do i=1,nd
               r(1:3)=rpos(1:3)+selist%pos(1:3)
               if(present(rotation_matrix)) r(:)=matmul(rotation_matrix,r)
               r(1:3)=r(1:3)+dpos(1:3)
               if(associated(selist%next)) selist=>selist%next
               do j=1,3
                  if(pb(j)) then
                     if(r(j).lt.cld(j)) r(j)=(cud(j)-cld(j)+1)+r(j)
                     if(r(j).gt.cud(j)) r(j)=-(cud(j)-cld(j)+1)+r(j)
                  endif
               enddo
               pincell=(r(1).ge.cld(1)).and.(r(1).le.cud(1)) &
                  .and.(r(2).ge.cld(2)).and.(r(2).le.cud(2)) &
                  .and.(r(3).ge.cld(3)).and.(r(3).le.cud(3))
               if(particleincell) particleincell=pincell
               if(pincell) then
                  testid=abs(cellmap(r(1),r(2),r(3)))
                  if(testid.eq.0) cycle
                  testid=root_monomer(testid)
                  if(testid.le.ipart) cycle
                  dup=.false.
                  do k=1,nnn
                     if(nnarray(k).eq.testid) then
                        dup=.true.
                        exit
                     endif
                  enddo
                  if(dup) cycle
                  nnn=nnn+1
                  nnarray(nnn)=testid
               endif
            enddo
         enddo
         end subroutine nearest_neighbor_test

         subroutine particle_uc_fit_operation(ipart, &
            cellmap,cld,cud,noparticleoverlap,particleincell, &
            periodic_boundary,set_uc_map,set_uc_map_index, &
            self_test,particle_displacement,rotation_matrix)
         implicit none
         type(element_list), pointer :: elist
         type(aggregate_type), pointer :: alist
         logical :: noparticleoverlap,particleincell,pincell,pb(3), &
            setucmap
         logical, optional :: periodic_boundary(3),set_uc_map,self_test
         integer :: nd,rpos(3),i,j,r(3),cld(3),cud(3),nagg,ia, &
                    ipart,dpos(3),itest,sucmi,cellt,comp
         integer, optional :: set_uc_map_index,particle_displacement(3)
         integer(4) :: cellmap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3)),cell
         integer, optional :: rotation_matrix(3,3)
         if(present(periodic_boundary)) then
            pb=periodic_boundary
         else
            pb=(/.true.,.true.,.false./)
         endif
         if(present(set_uc_map)) then
            setucmap=set_uc_map
         else
            setucmap=.false.
         endif
         if(present(set_uc_map_index)) then
            sucmi=set_uc_map_index
         else
            sucmi=1
         endif
         if(present(particle_displacement)) then
            dpos(1:3)=particle_displacement(1:3)+aggregate_position(1:3,ipart)
         else
            dpos(1:3)=aggregate_position(1:3,ipart)
         endif
         if(present(self_test)) then
            if(self_test) then
               itest=root_monomer(ipart)
            else
               itest=0
            endif
         else
            itest=0
         endif

         particleincell=.true.
         noparticleoverlap=.true.
         nagg=aggregate_number(ipart)
         alist=>aggregate_data(ipart)
         do ia=1,nagg
            if(.not.associated(alist)) then
               write(*,'('' error'')')
               stop
            endif
            rpos(1:3)=alist%mean_position(1:3)
            nd=alist%size
            elist=>alist%elements
            if(ia.lt.nagg) alist=>alist%next
            do i=1,nd
               r(1:3)=rpos(1:3)+elist%pos(1:3)
               comp=elist%pos(4)
               if(present(rotation_matrix)) r(:)=matmul(rotation_matrix,r)
               r(1:3)=r(1:3)+dpos(1:3)
               if(i.lt.nd) elist=>elist%next
               do j=1,3
                  if(pb(j)) then
                     if(r(j).lt.cld(j)) r(j)=(cud(j)-cld(j)+1)+r(j)
                     if(r(j).gt.cud(j)) r(j)=-(cud(j)-cld(j)+1)+r(j)
                  endif
               enddo
               pincell=(r(1).ge.cld(1)).and.(r(1).le.cud(1)) &
                  .and.(r(2).ge.cld(2)).and.(r(2).le.cud(2)) &
                  .and.(r(3).ge.cld(3)).and.(r(3).le.cud(3))
               if(particleincell) particleincell=pincell
               if(pincell) then
                  if(setucmap) then
                     if(comp.gt.0) then
                        cellmap(r(1),r(2),r(3))=sucmi
                     else
                        cellmap(r(1),r(2),r(3))=-sucmi
                     endif
                     cycle
                  endif
                  cell=abs(cellmap(r(1),r(2),r(3)))
                  if(cell.eq.0) cycle
                  cellt=root_monomer(cell)
                  if(cellt.eq.itest) cycle
                  noparticleoverlap=.false.
                  return
               endif
            enddo
         enddo
         end subroutine particle_uc_fit_operation

         subroutine generate_single_particle_list(partradius,partshape,partaspectratioy, &
            partaspectratioz,eulerang,dz,dlat,nd,ndsurf,particle_dlist, &
            particle_surf_dlist,shell_thick)
         use specialfuncs
         implicit none
         logical :: inparticle,inshell
         logical, allocatable :: cell(:,:,:)
         integer :: partshape,nd,nz,iz,nx,ix,ny,iy,ndloc,j, &
            cld(3),cud(3),ndsurf,jx,jy,jz,tpos
         integer(2), allocatable :: dlist(:,:)
         integer(2), optional :: particle_dlist(4,*),particle_surf_dlist(3,*)
         real(8) :: partradius,partaspectratioy,partaspectratioz,dz,dlat,csradius, &
            eulerang(3),z,xmax,x,rho,ymax,y,r,rpos(3),rposr(3),dcell(3),shellthick, &
            totpartradius
         real(8), optional :: shell_thick

         if(present(shell_thick)) then
            shellthick=shell_thick
         else
            shellthick=0.d0
         endif
         totpartradius=partradius+shellthick
         csradius=totpartradius*circumscribingsphere(partshape,partaspectratioy,partaspectratioz)
         nz=ceiling(2.*csradius/dz)
         dcell=(/dlat,dlat,dz/)
         ndloc=0
         cld=1000
         cud=-1000
         nd=ceiling(4.*3.14159*csradius**3/(dz*dlat*dlat))
         allocate(dlist(4,nd))
         do iz=-nz/2,nz/2
            z=dble(iz)*dz
            if(abs(z).gt.csradius) cycle
            xmax=sqrt(csradius*csradius-z*z)
            nx=ceiling(2.*xmax/dlat)
            do ix=-nx/2,nx/2
               x=dble(ix)*dlat
               rho=sqrt(x*x+z*z)
               if(rho.gt.csradius) cycle
               ymax=sqrt(csradius*csradius-rho*rho)
               ny=ceiling(2*ymax/dlat)
               do iy=-ny/2,ny/2
                  y=dble(iy)*dlat
                  r=sqrt(rho*rho+y*y)
                  if(r.gt.csradius) cycle
                  rpos=(/x,y,z/)
                  call eulerrotation(rpos,eulerang,1,rposr)
                  inparticle=insurf(rposr,partradius,partshape,partaspectratioy,partaspectratioz)
                  inshell=insurf(rposr,totpartradius,partshape,partaspectratioy,partaspectratioz)
                  if(inshell) then
                     ndloc=ndloc+1
                     dlist(1:3,ndloc)=int2((/ix,iy,iz/))
                     if(inparticle) then
                        dlist(4,ndloc)=1
                     else
                        dlist(4,ndloc)=-1
                     endif
                     cld(1)=min(cld(1),ix-2)
                     cud(1)=max(cud(1),ix+2)
                     cld(2)=min(cld(2),iy-2)
                     cud(2)=max(cud(2),iy+2)
                     cld(3)=min(cld(3),iz-2)
                     cud(3)=max(cud(3),iz+2)
                  endif
               enddo
            enddo
         enddo
         nd=ndloc
         do j=1,3
            tpos=sum(dlist(j,1:nd))/nd
            dlist(j,1:nd)=dlist(j,1:nd)-tpos
            if(present(particle_dlist)) then
               particle_dlist(j,1:nd)=dlist(j,1:nd)
            endif
         enddo
         if(present(particle_dlist)) then
            particle_dlist(4,1:nd)=dlist(4,1:nd)
         endif

         allocate(cell(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3)))
         cell=.true.
         do j=1,nd
            cell(dlist(1,j),dlist(2,j),dlist(3,j))=.false.
         enddo
         ndsurf=0
         do iz=cld(3)+1,cud(3)-1
            do iy=cld(2)+1,cud(2)-1
               do ix=cld(1)+1,cud(1)-1
                  if(.not.cell(ix,iy,iz)) cycle
                  do j=0,26
                     if(j.eq.13) cycle
                     jx=mod(j,3)-1+ix
                     jy=mod(j/3,3)-1+iy
                     jz=mod(j/9,3)-1+iz
                     if(.not.cell(jx,jy,jz)) then
                        ndsurf=ndsurf+1
                        if(present(particle_surf_dlist)) &
                           particle_surf_dlist(1:3,ndsurf)=int2((/ix,iy,iz/))
                        exit
                     endif
                  enddo
               enddo
            enddo
         enddo
         deallocate(cell,dlist)
         end subroutine generate_single_particle_list

         subroutine unit_cell_aggregation(ntotpart,npartin,cellmap,cld,cud,npartout, &
            periodic_boundary,diffusion_simulation,process_timings, &
            sticking_prob,rotation_prob)
         use specialfuncs
         use mpidefs
         implicit none
         logical :: particleincell,noparticleoverlap,pb(3),dsim, &
            sampok,rotate
         logical, optional :: periodic_boundary(3),diffusion_simulation
         integer :: npartin,npartout,cld(3),cud(3),i,j, &
                      i2,nnn,rmat(3,3), &
                      samp,maxsamp,ntotpart,randis(3),nnarray(ntotpart)
         integer(4) :: cellmap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3))
         real(8) :: stickprob,rotprob
         real(8), optional :: process_timings(5),sticking_prob,rotation_prob
         data maxsamp/100/

         if(present(periodic_boundary)) then
            pb=periodic_boundary
         else
            pb=(/.true.,.true.,.false./)
         endif
         if(present(diffusion_simulation)) then
            dsim=diffusion_simulation
         else
            dsim=.false.
         endif
         if(present(sticking_prob)) then
            stickprob=sticking_prob
         else
            stickprob=1.
         endif
         if(present(rotation_prob)) then
            rotprob=rotation_prob
         else
            rotprob=0.d0
         endif

         npartout=npartin
         if(present(process_timings)) process_timings(1)=mstm_mpi_wtime()

         do i=ntotpart,1,-1
            if(i.ne.root_monomer(i)) cycle
            call nearest_neighbor_test(i, &
               cellmap,cld,cud,particleincell,nnn,nnarray(1:), &
               periodic_boundary=pb)
            do j=nnn,1,-1
               if(stickprob.ge.dble(ran3(1))) then
                  i2=nnarray(j)
                  call combine_aggregate_data(i,i2)
                  npartout=npartout-1
               endif
            enddo
         enddo

         if(present(process_timings)) then
            process_timings(1)=mstm_mpi_wtime()-process_timings(1)
            process_timings(2)=mstm_mpi_wtime()
         endif
         if(dsim) then
            do i=1,ntotpart
               if(i.ne.root_monomer(i)) cycle
               sampok=.false.
               samp=0
               do while((.not.sampok).and.(samp.le.maxsamp))
                  samp=samp+1
                  randis(:)=(/-1+int(3*ran3(1)),-1+int(3*ran3(1)),-1+int(3*ran3(1))/)
                  rotate=(dble(ran3(1)).lt.rotprob)
                  if(rotate) then
                     call rotation_matrix(rmat)
                     call particle_uc_fit_operation(i, &
                        cellmap,cld,cud,noparticleoverlap,particleincell, &
                        periodic_boundary=pb, &
                        particle_displacement=randis, &
                        self_test=.true., &
                        rotation_matrix=rmat)
                  else
                     call particle_uc_fit_operation(i, &
                        cellmap,cld,cud,noparticleoverlap,particleincell, &
                        periodic_boundary=pb, &
                        particle_displacement=randis, &
                        self_test=.true.)
                  endif
                  sampok=particleincell.and.noparticleoverlap
               enddo
               if(sampok) then
                  call particle_uc_fit_operation(i, &
                     cellmap,cld,cud,noparticleoverlap,particleincell, &
                     periodic_boundary=pb,set_uc_map=.true.,set_uc_map_index=0)
                  aggregate_position(:,i)=aggregate_position(:,i)+randis(:)
                  if(rotate) call aggregate_rotate(i,rmat)
                  call particle_uc_fit_operation(i, &
                     cellmap,cld,cud,noparticleoverlap,particleincell, &
                     periodic_boundary=pb, &
                     set_uc_map=.true., &
                     set_uc_map_index=i)
               endif
            enddo
         endif
         if(present(process_timings)) then
            process_timings(2)=mstm_mpi_wtime()-process_timings(2)
         endif
         end subroutine unit_cell_aggregation

         subroutine combine_aggregate_data(j1,j2)
         implicit none
         integer :: i1,i2,j1,j2,pos1(3),pos2(3),pos(3),nd1,nd2
         type(aggregate_type), pointer :: alist1,alist2
         i1=min(j1,j2)
         i2=max(j1,j2)
         nd1=aggregate_size(i1)
         nd2=aggregate_size(i2)
         pos1(:)=aggregate_position(:,i1)
         pos2(:)=aggregate_position(:,i2)
         pos(:)=nint(dble(nd1*pos1(:)+nd2*pos2(:))/dble(nd1+nd2))
         pos1(:)=pos1(:)-pos(:)
         pos2(:)=pos2(:)-pos(:)
         alist1=>aggregate_data(i1)
         alist1%mean_position(:)=alist1%mean_position(:)+pos1(:)
         do while(associated(alist1%next))
            alist1=>alist1%next
            alist1%mean_position(:)=alist1%mean_position(:)+pos1(:)
         enddo
         alist2=>aggregate_data(i2)
         alist2%mean_position(:)=alist2%mean_position(:)+pos2(:)
         root_monomer(alist2%index)=i1
         do while(associated(alist2%next))
            alist2=>alist2%next
            alist2%mean_position(:)=alist2%mean_position(:)+pos2(:)
            root_monomer(alist2%index)=i1
         enddo
         alist2%next=>alist1%next
         alist1%next=>aggregate_data(i2)
         aggregate_number(i1)=aggregate_number(i1)+aggregate_number(i2)
         aggregate_position(:,i1)=pos(:)
         aggregate_position(:,i2)=pos(:)
         aggregate_number(i2)=0
         aggregate_size(i1)=nd1+nd2
         aggregate_size(i2)=0
         end subroutine combine_aggregate_data

         subroutine rotation_matrix(rmat)
         use specialfuncs
         implicit none
         integer :: j,jl,jr,cp(3),sp(3),rmat(3,3), &
                rmatt(3,3,3),cs(3),ss(3)
         cs=-1+2*nint((/ran3(1),ran3(1),ran3(1)/))
         ss=-1+2*nint((/ran3(1),ran3(1),ran3(1)/))
         cp=nint((/ran3(1),ran3(1),ran3(1)/))
         sp=1-cp
         cp=cs*cp
         sp=ss*sp
         do j=1,3
            jl=mod(j+1,3)+1
            jr=mod(j,3)+1
            rmatt(j,j,j)=1
            rmatt(j,jl,j)=0
            rmatt(j,jr,j)=0
            rmatt(jl,j,j)=0
            rmatt(jr,j,j)=0
            rmatt(jl,jl,j)=cp(j)
            rmatt(jr,jr,j)=cp(j)
            rmatt(jl,jr,j)=sp(j)
            rmatt(jr,jl,j)=-sp(j)
         enddo
         rmat=matmul(rmatt(:,:,1),matmul(rmatt(:,:,2),rmatt(:,:,3)))
         end subroutine rotation_matrix

         subroutine aggregate_rotate(ipart,rmat)
         implicit none
         type(element_list), pointer :: elist
         type(surf_element_list), pointer :: selist
         type(aggregate_type), pointer :: alist
         integer :: ipart,i,j,rmat(3,3),fpos(3)

         alist=>aggregate_data(ipart)
         do i=1,aggregate_number(ipart)
            fpos(:)=alist%mean_position(:)
            alist%mean_position(:)=matmul(rmat,fpos)
            elist=>alist%elements
            do j=1,alist%size
               fpos(:)=elist%pos(1:3)
               elist%pos(1:3)=matmul(rmat,fpos)
               if(j.lt.alist%size) elist=>elist%next
            enddo
            selist=>alist%surf_elements
            do j=1,alist%surf_size
               fpos(1:3)=selist%pos(1:3)
               selist%pos(:)=matmul(rmat,fpos)
               if(j.lt.alist%surf_size) selist=>selist%next
            enddo
            if(i.lt.aggregate_number(ipart)) alist=>alist%next
         enddo
         end subroutine aggregate_rotate

         subroutine write_cell_map(cellmap,cld,cud,outunit)
         implicit none
         logical, save :: firstrun
         integer :: cld(3),cud(3),nd,ix,iy,iz,outunit,i,j
         integer(4) :: cellmap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3))
         data firstrun/.true./
         nd=0
         if(firstrun) then
            write(outunit,*) cld(1),cud(1),cld(2),cud(2),cld(3),cud(3)
            firstrun=.false.
         endif
         do i=1,2
            if(i.eq.1) then
               nd=0
            else
               write(outunit,*) nd
            endif
            do iz=cld(3),cud(3)
               do iy=cld(2),cud(2)
                  do ix=cld(1),cud(1)
                     if(cellmap(ix,iy,iz).ne.0) then
                        if(i.eq.1) then
                           nd=nd+1
                        else
                           j=cellmap(ix,iy,iz)
                           if(j.ne.0) j=root_monomer(j)
                           write(outunit,'(4i8)') ix,iy,iz,j
                        endif
                     endif
                  enddo
               enddo
            enddo
         enddo
         end subroutine write_cell_map

         subroutine write_monomer_positions(nmonomer,iunit)
         implicit none
         integer :: nmonomer,iunit,i,j
         type(aggregate_type), pointer :: alist
         write(iunit,*) nmonomer
         do i=1,nmonomer
            alist=>aggregate_data(i)
            do j=1,aggregate_number(i)
               write(iunit,'(3i10)') alist%mean_position(:)+aggregate_position(:,i)
               if(j.lt.aggregate_number(i)) then
                  alist=>alist%next
               endif
            enddo
         enddo
         end subroutine write_monomer_positions

      end module target_generation
