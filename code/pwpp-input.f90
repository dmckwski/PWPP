!june 18 original
!27 july: core volume fraction added to output

      module inputinterface
      implicit none
      logical :: loop_job,repeat_run,ml_specified,first_run,ndrs_specified, &
                nd_specified,square_target,configuration_ok
      logical, target :: homogeneous_layer,append_output_file, &
                 print_scattering_matrix,azimuth_average, &
                 copy_input_file,sphere_1_fixed,number_particles_specified, &
                 write_particle_data,auto_ri_medium,d_fixed,equal_d, &
                 ordered_placement,incident_scattering_plane, &
                 print_absorption_distribution,coherent_averaging, &
                 auto_zero_s12,auto_update_ri_medium, &
                 enforce_zero_absorption,circular_window, &
                 rotational_average_epsilon,enforce_symmetry, &
                 print_at_quad_points,read_particle_data, &
                 periodic_sides,ot_specified,static_configuration, &
                 polarizability_form
      integer :: n_nest_loops,i_var_start(5),i_var_stop(5),i_var_step(5), &
                 loop_component_number(5),run_number,n_configuration_groups, &
                 number_particles,run_print_unit,n_angles,n_monomers,n_aggregates
      integer, target :: max_iterations,ndx_width,ndy_width,n_components,part_shape(5), &
                         nd_thickness,boundary_type,nd_width, &
                         s_incident_x,s_incident_y,n_configurations, &
                         max_number_particles,ndrsx_width,ndrsy_width,spline_order, &
                         number_aggregation_cycles,number_diffusion_cycles, &
                         nddsx_width,nddsy_width
      real(8) :: r_var_start(5),r_var_stop(5),r_var_step(5),extinction_coefficient, &
                 set_dipole_spacing,target_volume_fraction,set_lateral_dipole_spacing, &
                 core_volume_fraction
      real(8), private :: pi=3.141592653589793d0
      real (8), target :: incident_beta_deg,incident_alpha_deg,solution_eps,k_cut, &
            target_thickness,part_mean_radius(5),part_sigma(5),part_ary(5), &
            part_arz(5),part_fv(5),dipole_spacing,precon_absorption,target_x_width,target_y_width, &
            target_width,lateral_dipole_spacing,part_shell_thickness(5), &
            interpolation_increment,rotation_probability, &
            sticking_probability,mie_optical_thickness,aggregation_fraction
      complex(8) :: c_var_start(5),c_var_stop(5),c_var_step(5),effective_ri(1)
      complex(8), target :: ri_medium,part_ri(5),ri_back_surface,ri_front_surface, &
                            ri_binder,part_shell_ri(5)
      character*1 :: loop_var_type(5)
      character*20 :: run_date_and_time
      character*128 :: loop_var_label(5),input_file
      character*128, target :: output_file,run_file,dipole_map_input_file, &
            dipole_map_output_file,particle_data_file,particle_data_input_file
      data run_print_unit/6/
      data loop_job,repeat_run,first_run/.false.,.false.,.true./
      data homogeneous_layer/.false./
      data append_output_file/.false./
      data copy_input_file/.false./
      data sphere_1_fixed/.false./
      data number_particles_specified/.false./
      data auto_ri_medium/.true./
      data precon_absorption/0.001d0/
      data print_scattering_matrix/.true./
      data print_absorption_distribution/.false./
      data azimuth_average/.true./
      data n_nest_loops/0/
      data run_number/0/
      data loop_component_number/5*1/
      data i_var_start,i_var_stop,i_var_step/5*0,5*0,5*0/
      data r_var_start,r_var_stop,r_var_step/5*0.d0,5*0.d0,5*0.d0/
      data c_var_start,c_var_stop,c_var_step/5*(0.d0,0.d0),5*(0.d0,0.d0),5*(0.d0,0.d0)/
      data max_iterations/10000/
      data nd_width,ndx_width,ndy_width/500,500,500/
      data nd_thickness/10/
      data boundary_type/0/
      data ml_specified/.false./
      data nd_specified/.true./
      data ndrs_specified/.false./
      data n_components/1/
      data max_number_particles/100000000/
      data n_configurations/1/
      data part_shape/5*1/
      data incident_beta_deg/0.d0/
      data incident_alpha_deg/0.d0/
      data s_incident_x/0/
      data s_incident_y/0/
      data solution_eps/1.d-6/
      data k_cut/2.5d0/
      data target_width,target_x_width,target_y_width/100.d0,100.d0,100.d0/
      data target_thickness/2.d0/
      data dipole_spacing/.2d0/
      data lateral_dipole_spacing/.2d0/
      data part_mean_radius/5*1.d0/
      data part_sigma/5*0.d0/
      data part_ary/5*1.d0/
      data part_arz/5*1.d0/
      data part_fv/.1d0,0.d0,0.d0,0.d0,0.d0/
      data ri_medium/(1.d0,0.d0)/
      data ri_back_surface/(1.d0,0.d0)/
      data ri_front_surface/(1.d0,0.d0)/
      data ri_binder/(1.d0,0.d0)/
      data part_ri/(1.54d0,0.d0),(1.d0,0.d0),(1.d0,0.d0),(1.d0,0.d0),(1.d0,0.d0)/
      data output_file/'ppddtest.dat'/
      data run_file/'run1.dat'/
      data dipole_map_input_file/' '/
      data dipole_map_output_file/' '/
      data write_particle_data/.false./
      data particle_data_file/'particledatatest.dat'/
      data interpolation_increment/.01d0/
      data spline_order/2/
      data d_fixed/.false./
      data equal_d/.true./
      data ordered_placement/.true./
      data incident_scattering_plane/.false./
      data coherent_averaging/.false./
      data auto_zero_s12/.true./
      data auto_update_ri_medium/.false./
      data enforce_zero_absorption/.false./
      data enforce_symmetry/.true./
      data print_at_quad_points/.false./
      data read_particle_data/.false./
      data particle_data_input_file/' '/
      data periodic_sides/.true./
      data circular_window/.false./
      data number_aggregation_cycles/0/
      data number_diffusion_cycles/0/
      data rotation_probability/0.d0/
      data sticking_probability/1.d0/
      data ot_specified/.false./
      data mie_optical_thickness/0.d0/
      data aggregation_fraction/1.d0/
      data part_shell_thickness/5*0.d0/
      data static_configuration/.false./
      data polarizability_form/.false./
      data nddsx_width,nddsy_width/100000,100000/

      contains

         subroutine variable_list_operation(varlabel, &
            var_value,var_type, &
            var_position,var_operation,var_status, &
            i_var_pointer,r_var_pointer,c_var_pointer)
         implicit none
         logical :: operate
         logical, pointer :: lvarvalue
         integer :: varpos,varstatus
         integer, optional :: var_position,var_status
         integer, pointer :: ivarvalue
         integer, optional, pointer :: i_var_pointer
         real(8), pointer :: rvarvalue
         real(8), optional, pointer :: r_var_pointer
         complex(8), pointer :: cvarvalue
         complex(8), optional, pointer :: c_var_pointer
         character*1 :: vartype
         character*1, optional :: var_type
         character*(*), optional :: var_value,var_operation
         character*128 :: varop,sentvarvalue,varlabel
         character*128, pointer :: avarvalue

         if(present(var_operation)) then
            varop=trim(var_operation)
         else
            varop=' '
         endif
         if(present(var_value)) then
            sentvarvalue=trim(var_value)
            operate=.true.
         else
            sentvarvalue=' '
            operate=.false.
         endif
         if(present(var_position)) then
            varpos=var_position
         else
            varpos=1
         endif
         varstatus=0
         vartype='n'

         if(varlabel.eq.'output_file') then
            vartype='a'
            avarvalue=>output_file

         elseif(varlabel.eq.'append_output_file') then
            vartype='l'
            lvarvalue=>append_output_file

         elseif(varlabel.eq.'copy_input_file') then
            vartype='l'
            lvarvalue=>copy_input_file

         elseif(varlabel.eq.'run_file') then
            vartype='a'
            avarvalue=>run_file

         elseif(varlabel.eq.'dipole_map_input_file') then
            vartype='a'
            avarvalue=>dipole_map_input_file

         elseif(varlabel.eq.'dipole_map_output_file') then
            vartype='a'
            avarvalue=>dipole_map_output_file

         elseif(varlabel.eq.'write_particle_data') then
            vartype='l'
            lvarvalue=>write_particle_data

         elseif(varlabel.eq.'particle_data_file') then
            vartype='a'
            avarvalue=>particle_data_file

         elseif(varlabel.eq.'read_particle_data') then
            vartype='l'
            lvarvalue=>read_particle_data

         elseif(varlabel.eq.'particle_data_input_file') then
            vartype='a'
            avarvalue=>particle_data_input_file

         elseif(varlabel.eq.'homogeneous_layer') then
            vartype='l'
            lvarvalue=>homogeneous_layer

         elseif(varlabel.eq.'sphere_1_fixed') then
            vartype='l'
            lvarvalue=>sphere_1_fixed

         elseif(varlabel.eq.'ordered_placement') then
            vartype='l'
            lvarvalue=>ordered_placement

         elseif(varlabel.eq.'circular_window') then
            vartype='l'
            lvarvalue=>circular_window

         elseif(varlabel.eq.'auto_ri_medium') then
            vartype='l'
            lvarvalue=>auto_ri_medium

         elseif(varlabel.eq.'auto_update_ri_medium') then
            vartype='l'
            lvarvalue=>auto_update_ri_medium

         elseif(varlabel.eq.'enforce_zero_absorption') then
            vartype='l'
            lvarvalue=>enforce_zero_absorption

         elseif(varlabel.eq.'enforce_symmetry') then
            vartype='l'
            lvarvalue=>enforce_symmetry

         elseif(varlabel.eq.'polarizability_form') then
            vartype='l'
            lvarvalue=>polarizability_form

         elseif(varlabel.eq.'precon_absorption') then
            vartype='r'
            rvarvalue=>precon_absorption

         elseif(varlabel.eq.'max_iterations') then
            vartype='i'
            ivarvalue=>max_iterations

         elseif(varlabel.eq.'solution_eps') then
            vartype='r'
            rvarvalue=>solution_eps

         elseif(varlabel.eq.'n_configurations') then
            vartype='i'
            ivarvalue=>n_configurations

         elseif(varlabel.eq.'nd_width') then
            vartype='i'
            ivarvalue=>ndx_width
            nd_specified=.true.
            square_target=.true.

         elseif(varlabel.eq.'ndx_width') then
            vartype='i'
            ivarvalue=>ndx_width
            nd_specified=.true.
            square_target=.false.

         elseif(varlabel.eq.'ndy_width') then
            vartype='i'
            ivarvalue=>ndy_width
            nd_specified=.true.
            square_target=.false.

         elseif(varlabel.eq.'ndds_width') then
            vartype='i'
            ivarvalue=>nddsx_width

         elseif(varlabel.eq.'nddsx_width') then
            vartype='i'
            ivarvalue=>nddsx_width

         elseif(varlabel.eq.'nddsy_width') then
            vartype='i'
            ivarvalue=>nddsy_width

         elseif(varlabel.eq.'target_width') then
            vartype='r'
            rvarvalue=>target_width
            nd_specified=.false.
            square_target=.true.

         elseif(varlabel.eq.'target_x_width') then
            vartype='r'
            rvarvalue=>target_x_width
            nd_specified=.false.
            square_target=.false.

         elseif(varlabel.eq.'target_y_width') then
            vartype='r'
            rvarvalue=>target_y_width
            nd_specified=.false.
            square_target=.false.

         elseif(varlabel.eq.'ndrsx_width') then
            vartype='i'
            ivarvalue=>ndrsx_width
            ndrs_specified=.true.

         elseif(varlabel.eq.'ndrsy_width') then
            vartype='i'
            ivarvalue=>ndrsy_width
            ndrs_specified=.true.

         elseif(varlabel.eq.'k_cut') then
            vartype='r'
            rvarvalue=>k_cut
            ndrs_specified=.false.

         elseif(varlabel.eq.'target_thickness') then
            vartype='r'
            rvarvalue=>target_thickness
            ml_specified=.false.

         elseif(varlabel.eq.'nd_thickness') then
            vartype='i'
            ivarvalue=>nd_thickness
            ml_specified=.true.

         elseif(varlabel.eq.'boundary_type') then
            vartype='i'
            ivarvalue=>boundary_type

         elseif(varlabel.eq.'dipole_spacing') then
            vartype='r'
            rvarvalue=>dipole_spacing

         elseif(varlabel.eq.'lateral_dipole_spacing') then
            vartype='r'
            rvarvalue=>lateral_dipole_spacing
            equal_d=.false.

         elseif(varlabel.eq.'d_fixed') then
            vartype='l'
            lvarvalue=>d_fixed

         elseif(varlabel.eq.'equal_d') then
            vartype='l'
            lvarvalue=>equal_d

         elseif(varlabel.eq.'incident_beta_deg') then
            vartype='r'
            rvarvalue=>incident_beta_deg

         elseif(varlabel.eq.'incident_alpha_deg') then
            vartype='r'
            rvarvalue=>incident_alpha_deg

         elseif(varlabel.eq.'s_incident_x') then
            vartype='i'
            ivarvalue=>s_incident_x

         elseif(varlabel.eq.'s_incident_y') then
            vartype='i'
            ivarvalue=>s_incident_y

         elseif(varlabel.eq.'print_scattering_matrix') then
            vartype='l'
            lvarvalue=>print_scattering_matrix

         elseif(varlabel.eq.'print_absorption_distribution') then
            vartype='l'
            lvarvalue=>print_absorption_distribution

         elseif(varlabel.eq.'azimuth_average') then
            vartype='l'
            lvarvalue=>azimuth_average

         elseif(varlabel.eq.'incident_scattering_plane') then
            vartype='l'
            lvarvalue=>incident_scattering_plane

         elseif(varlabel.eq.'max_number_particles') then
            vartype='i'
            ivarvalue=>max_number_particles

         elseif(varlabel.eq.'number_particles_specified') then
            vartype='l'
            lvarvalue=>number_particles_specified

         elseif(varlabel.eq.'part_shape') then
            vartype='i'
            ivarvalue=>part_shape(varpos)

         elseif(varlabel.eq.'part_mean_radius') then
            vartype='r'
            rvarvalue=>part_mean_radius(varpos)

         elseif(varlabel.eq.'part_sigma') then
            vartype='r'
            rvarvalue=>part_sigma(varpos)

         elseif(varlabel.eq.'part_ary') then
            vartype='r'
            rvarvalue=>part_ary(varpos)

         elseif(varlabel.eq.'part_arz') then
            vartype='r'
            rvarvalue=>part_arz(varpos)

         elseif(varlabel.eq.'part_fv') then
            vartype='r'
            rvarvalue=>part_fv(varpos)

         elseif(varlabel.eq.'part_ri') then
            vartype='c'
            cvarvalue=>part_ri(varpos)

         elseif(varlabel.eq.'part_shell_thickness') then
            vartype='r'
            rvarvalue=>part_shell_thickness(varpos)

         elseif(varlabel.eq.'part_shell_ri') then
            vartype='c'
            cvarvalue=>part_shell_ri(varpos)

         elseif(varlabel.eq.'ri_medium') then
            vartype='c'
            cvarvalue=>ri_medium

         elseif(varlabel.eq.'ri_binder') then
            vartype='c'
            cvarvalue=>ri_binder

         elseif(varlabel.eq.'ri_back_surface') then
            vartype='c'
            cvarvalue=>ri_back_surface

         elseif(varlabel.eq.'ri_front_surface') then
            vartype='c'
            cvarvalue=>ri_front_surface

         elseif(varlabel.eq.'interpolation_increment') then
            vartype='r'
            rvarvalue=>interpolation_increment

         elseif(varlabel.eq.'spline_order') then
            vartype='i'
            ivarvalue=>spline_order

         elseif(varlabel.eq.'coherent_averaging') then
            vartype='l'
            lvarvalue=>coherent_averaging

         elseif(varlabel.eq.'auto_zero_s12') then
            vartype='l'
            lvarvalue=>auto_zero_s12

         elseif(varlabel.eq.'print_at_quad_points') then
            vartype='l'
            lvarvalue=>print_at_quad_points

         elseif(varlabel.eq.'periodic_sides') then
            vartype='l'
            lvarvalue=>periodic_sides

         elseif(varlabel.eq.'number_aggregation_cycles') then
            vartype='i'
            ivarvalue=>number_aggregation_cycles

         elseif(varlabel.eq.'number_diffusion_cycles') then
            vartype='i'
            ivarvalue=>number_diffusion_cycles

         elseif(varlabel.eq.'aggregation_fraction') then
            vartype='r'
            rvarvalue=>aggregation_fraction

         elseif(varlabel.eq.'rotation_probability') then
            vartype='r'
            rvarvalue=>rotation_probability

         elseif(varlabel.eq.'sticking_probability') then
            vartype='r'
            rvarvalue=>sticking_probability

         elseif(varlabel.eq.'mie_optical_thickness') then
            vartype='r'
            rvarvalue=>mie_optical_thickness
            ot_specified=.true.

         elseif(varlabel.eq.'ot_specified') then
            vartype='l'
            lvarvalue=>ot_specified

         elseif(varlabel.eq.'static_configuration') then
            vartype='l'
            lvarvalue=>static_configuration

         endif

         if(vartype.eq.'n') then
            varstatus=1
            if(present(var_status)) var_status=varstatus
            return
         endif
         if(present(var_type)) var_type=vartype
         if(present(i_var_pointer)) i_var_pointer=>ivarvalue
         if(present(r_var_pointer)) r_var_pointer=>rvarvalue
         if(present(c_var_pointer)) c_var_pointer=>cvarvalue

         if(operate) then
            if(vartype.eq.'i') then
               call set_string_to_int_variable(sentvarvalue, &
                  ivarvalue,var_operation=varop)
            elseif(vartype.eq.'r') then
               call set_string_to_real_variable(sentvarvalue, &
                  rvarvalue,var_operation=varop)
            elseif(vartype.eq.'c') then
               call set_string_to_cmplx_variable(sentvarvalue, &
                  cvarvalue,var_operation=varop)
            elseif(vartype.eq.'l') then
               call set_string_to_logical_variable(sentvarvalue, &
                  lvarvalue,var_operation=varop)
            elseif(vartype.eq.'a') then
               avarvalue=sentvarvalue
            endif
            if(n_components.lt.5) then
               part_mean_radius(n_components+1:5)=part_mean_radius(n_components)
               part_shape(n_components+1:5)=part_shape(n_components)
               part_fv(n_components+1:5)=part_fv(n_components)
               part_sigma(n_components+1:5)=part_sigma(n_components)
               part_ary(n_components+1:5)=part_ary(n_components)
               part_arz(n_components+1:5)=part_arz(n_components)
               part_ri(n_components+1:5)=part_ri(n_components)
            endif
         endif
         end subroutine variable_list_operation

         subroutine inputdata(inputfiledata,read_status)
         use mpidefs
         implicit none
         integer :: readok,n,ncomp,varstat,rank,stopit
         integer, save :: inputline
         integer, optional :: read_status
         character*128 :: parmid,parmval,varop,inputfiledata(*)
         data inputline/1/

         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         readok=0
         stopit=0
         ncomp=n_components
         do while(readok.eq.0)
            parmid=inputfiledata(inputline)
            inputline=inputline+1
            if(trim(parmid).eq.'run_file') then
               parmval=inputfiledata(inputline)
               inputline=inputline+1
               if(trim(parmval).ne.' ') then
                  run_print_unit=3
                  if(rank.eq.0) then
                     open(3,file=trim(parmval))
                  endif
               endif
               cycle
            endif

            if(parmid(1:1).eq.'!'.or.parmid(1:1).eq.'%') then
               cycle
            endif

            if(trim(parmid).eq.'loop_variable') then
               loop_job=.true.
               n_nest_loops=n_nest_loops+1
               n=n_nest_loops
               parmid=inputfiledata(inputline)
               inputline=inputline+1
               if(trim(parmid).eq.'component_number') then
                  read(inputfiledata(inputline),*) ncomp
                  inputline=inputline+1
                  loop_component_number(n)=ncomp
                  parmid=inputfiledata(inputline)
                  inputline=inputline+1
               else
                  ncomp=1
               endif
               loop_var_label(n)=parmid
               if(trim(parmid).eq.'nd_thickness') then
                  ml_specified=.true.
               elseif(trim(parmid).eq.'target_thickness') then
                  ml_specified=.false.
               endif
               call variable_list_operation(loop_var_label(n), &
                    var_type=loop_var_type(n))
               if(loop_var_type(n).eq.'i') then
                  read(inputfiledata(inputline),*) i_var_start(n),i_var_stop(n),i_var_step(n)
               elseif(loop_var_type(n).eq.'r') then
                  read(inputfiledata(inputline),*) r_var_start(n),r_var_stop(n),r_var_step(n)
               elseif(loop_var_type(n).eq.'c') then
                  read(inputfiledata(inputline),*) c_var_start(n),c_var_stop(n),c_var_step(n)
               endif
               inputline=inputline+1
               cycle

            elseif(trim(parmid).eq.'new_component') then
               n_components=n_components+1
               ncomp=n_components
               cycle

            elseif(trim(parmid).eq.'new_run') then
               repeat_run=.true.
               exit

            elseif(trim(parmid).eq.'end_of_options') then
               repeat_run=.false.
               readok=-1
               exit
            else
               varstat=0
               call variable_list_operation(parmid, &
                   var_status=varstat)
               if(varstat.ne.0) then
                  if(rank.eq.0) then
                     write(run_print_unit,'('' unknown input parameter:'',a)') trim(parmid)
                     call flush(run_print_unit)
                     stopit=1
                  endif
                  cycle
               else
                  parmval=inputfiledata(inputline)
                  inputline=inputline+1
                  if(readok.ne.0) cycle
                  parmval=trim(parmval)
                  varop='assign'
                  call variable_list_operation(parmid,var_value=parmval, &
                      var_position=ncomp,var_operation='assign', &
                      var_status=varstat)
               endif
            endif
         enddo
         if(stopit.eq.1) stop
         if(present(read_status)) read_status=varstat
         end subroutine inputdata

         subroutine main_calling_program()
         use specialfuncs
         use intrinsics
         use mpidefs
         use solver
         use bspline_sub_module
         implicit none
         logical :: writepartdata,coherentfix,incinc,oplace
         integer :: maxiter,ndx,ndy,ndsx,ndsy,ml,mlsamp,nd2,nds2,ndpadx,ndpady, &
                    ndpad2,rank,numprocs,j,iter, &
                    nconfiggroups,configcolor,configgroup, &
                    configcomm,configrank,config0comm,config, &
                    is1x,is2x,is1y,is2y,niter,nsend,ntotpart,ss0,idum(3), &
                    partshape,nang,nax1,nax2,nay1,nay2, &
                    nax,nay,config0,nconfigave, &
                    numprocsperconfig
         integer, save :: ndsold,ndold
         real(8) :: betadeg,alphadeg,thick,d,k0x,k0y,wx,wy,dset,sinb, &
                    extinction(2),hflux(2,2,2),cflux(2),fvc,kx,ky, &
                    starttime,time1,time2,gamma,partdata(9),errmax, &
                    alpharad,dw,hfluxave(2,2,2),cfluxave(2),setthickness, &
                    dtheta,cosb,qext,qsca,g11,pcabs
         real(8), save :: wold,k0xold,k0yold
         real(8), allocatable :: smat(:,:,:,:),smatave(:,:,:,:), &
            qabsj(:),qabsjave(:),scohmat(:,:,:,:)
         complex(8) :: alpha0,e0
         complex(8), allocatable :: acoef(:,:,:,:),amean(:,:,:,:), &
                    acoh(:,:,:),axtemp(:)
         data ndsold,ndold,wold,k0xold,k0yold/0,0,0.d0,0.d0,0.d0/

         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)

         if(square_target) then
            ndx_width=nd_width
            ndy_width=nd_width
            target_x_width=target_width
            target_y_width=target_width
         endif

         if(dipole_map_input_file.ne.' ') then
            if(rank.eq.0) then
               open(4,file=trim(dipole_map_input_file))
               read(4,*) ndx,ndy,ml
               close(4)
               do
                  if(checkn235(ndx)) exit
                  ndx=ndx+1
               enddo
               do
                  if(checkn235(ndy)) exit
                  ndy=ndy+1
               enddo
            endif
            idum=(/ndx,ndy,ml/)
            call mstm_mpi(mpi_command='bcast', &
                 mpi_number=3,mpi_send_buf_i=idum, &
                 mpi_rank=0)
            ndx_width=idum(1)
            ndy_width=idum(2)
            nd_thickness=idum(3)
            ml_specified=.true.
            nd_specified=.true.
            boundary_type=0
         endif

         if(nd_specified) then
            do
               if(checkn235(ndx_width)) exit
               ndx_width=ndx_width+1
            enddo
            do
               if(checkn235(ndy_width)) exit
               ndy_width=ndy_width+1
            enddo
         endif
         maxiter=max_iterations
         ndx=ndx_width
         ndy=ndy_width
         betadeg=incident_beta_deg
         alphadeg=incident_alpha_deg
         alpharad=alphadeg*pi/180.d0
         thick=target_thickness
         ml=nd_thickness
         dset=dipole_spacing
         set_dipole_spacing=dset
         pcabs=precon_absorption
         if(equal_d) then
            lateral_dipole_spacing=dset
         endif
         dw=lateral_dipole_spacing
         set_lateral_dipole_spacing=dw
         if(betadeg.eq.0.d0) then
            sinb=0.d0
            cosb=1.d0
         else
            sinb=sin(betadeg*pi/180.d0)
            cosb=cos(betadeg*pi/180.d0)
         endif
         k0x=sinb*cos(alphadeg*pi/180.d0)
         k0y=sinb*sin(alphadeg*pi/180.d0)
         dtheta=interpolation_increment
         oplace=ordered_placement
         if(sinb.eq.0.d0) then
            gamma=incident_alpha_deg*pi/180.d0
         else
            gamma=0.d0
         endif
         if(maxval(part_fv(1:n_components)).eq.0.d0.and. &
               (.not.number_particles_specified)) then
            homogeneous_layer=.true.
         endif

         if(nd_specified) then
            wx=dw*dble(ndx)
            wy=dw*dble(ndy)
         else
            ndx_width=ceiling(target_x_width/dw)
            ndy_width=ceiling(target_y_width/dw)
            do
               if(checkn235(ndx_width)) exit
               ndx_width=ndx_width+1
            enddo
            do
               if(checkn235(ndy_width)) exit
               ndy_width=ndy_width+1
            enddo
            dw=target_x_width/dble(ndx_width)
            target_y_width=dw*dble(ndy_width)
            ndx=ndx_width
            ndy=ndy_width
            wx=target_x_width
            wy=target_y_width
         endif

         if(equal_d) then
            dw=d
            wx=dble(ndx)*d
            wy=dble(ndy)*d
         endif

         if(number_particles_specified.and.(n_components.eq.1)) then
            boundary_type=0
            if(max_number_particles.eq.1) then
               ml_specified=.true.
!               ml=ceiling(part_mean_radius(1)*2.d0/dset)+2
               thick=2.1d0*part_mean_radius(1)*circumscribingsphere(part_shape(1), &
                   aspect_ratio_y=part_ary(1), &
                   aspect_ratio_z=part_arz(1))
               ml=ceiling(thick/dset)+2
!               sphere_1_fixed=.true.
               target_thickness=thick
            elseif(ml_specified) then
               ml=ceiling(dble(max_number_particles) &
                  *(part_mean_radius(1)**3.)*4.d0*pi/(3.d0*part_fv(1)*dble(ndx*ndy)*dset**3))
               ml=ml+4.*ceiling(part_mean_radius(1)/dset)
            else
!               part_fv(1)=4.d0*pi*(part_mean_radius(1)**3)*dble(max_number_particles) &
!                  /(3.d0*w*w*(target_thickness-2.d0*part_mean_radius(1)))
               part_fv(1)=4.d0*pi*(part_mean_radius(1)**3)*dble(max_number_particles) &
                  /(3.d0*wx*wy*target_thickness)
            endif
            if(part_fv(1).le.0.1) then
               oplace=.false.
            else
               oplace=ordered_placement
            endif
         endif

         if(ml_specified) then
            nd_thickness=ml
            d=dset
            thick=d*dble(ml)
            target_thickness=thick
         else
            ml=ceiling(thick/dset)
            nd_thickness=ml
            if(d_fixed) then
               d=dset
            else
               d=thick/dble(ml)
            endif
         endif

         if(ot_specified) then
            call mieext(part_mean_radius(1)*dble(ri_binder),part_ri(1)/ri_binder,qext,qsca,g11)
            part_fv(1)=4.d0*part_mean_radius(1)*mie_optical_thickness/(3.d0*thick*qext)
         endif

         if(auto_ri_medium.and.pcabs.lt.0.) then
            call mieext(part_mean_radius(1)*dble(ri_binder),part_ri(1)/ri_binder,qext,qsca,g11)
            extinction_coefficient=3.d0*part_fv(1)*qext/4.d0/part_mean_radius(1)
            pcabs=extinction_coefficient/2.d0
            precon_absorption=-pcabs
         endif

         if(boundary_type.eq.0) then
            setthickness=thick-2.d0*part_mean_radius(1)
         else
            setthickness=thick+2.d0*part_mean_radius(1)
         endif

         set_lateral_dipole_spacing=dw
         set_dipole_spacing=d

         kx=k0x+2.d0*pi*s_incident_x/wx
         ky=k0y+2.d0*pi*s_incident_y/wy
         if(ndrs_specified) then
            ndsx=ndrsx_width
            do
               if(checkn235(ndsx).and.checkn235(2*ndsx)) exit
               ndsx=ndsx+1
            enddo
            ndsy=ndrsy_width
            do
               if(checkn235(ndsy).and.checkn235(2*ndsy)) exit
               ndsy=ndsy+1
            enddo
            k_cut=pi/wx*dble(ndsx)
         else
            ndsx=ceiling(k_cut*wx/pi)
            ndsy=ceiling(k_cut*wy/pi)
            do
               if(checkn235(ndsx).and.checkn235(2*ndsx)) exit
               ndsx=ndsx+1
            enddo
            do
               if(checkn235(ndsy).and.checkn235(2*ndsy)) exit
               ndsy=ndsy+1
            enddo
         endif

         ndsx=min(ndsx,ndx)
         ndrsx_width=ndsx
         ndsy=min(ndsy,ndy)
         ndrsy_width=ndsy
         ndpadx=2*ndsx
         ndpadx=min(ndpadx,ndx)
         ndpady=2*ndsy
         ndpady=min(ndpady,ndy)
         nd2=ndx*ndy
         nds2=ndsx*ndsy
         ndpad2=ndpadx*ndpady
         is1x=-(ndsx/2)
         is2x=(ndsx-1)/2
         is1y=-(ndsy/2)
         is2y=(ndsy-1)/2
         ss0=-is1x+1-ndsx*is1y
         target_x_width=wx
         target_y_width=wy
         mlsamp=ml

         ndx_cv=ndx
         ndy_cv=ndy
         k0x_cv=k0x
         k0y_cv=k0y
         wx_cv=wx
         wy_cv=wy
         h_cv=thick
         d_cv=d
         dw_cv=dw
         rimedium_cv=ri_medium
         ri0_cv=ri_front_surface
         ris_cv=ri_back_surface
         run_print_unit_cv=run_print_unit

         if(maxiter.lt.0) then
            numprocsperconfig=2
         else
            numprocsperconfig=4
         endif
         nconfiggroups=numprocs/numprocsperconfig
         nconfiggroups=max(nconfiggroups,1)
         n_configuration_groups=nconfiggroups
         configcolor=floor(dble(nconfiggroups*rank)/dble(numprocs))
         configgroup=configcolor
         call mstm_mpi(mpi_command='split', &
              mpi_color=configcolor,mpi_key=rank, &
              mpi_new_comm=configcomm)
         call mstm_mpi(mpi_command='rank', &
              mpi_rank=configrank, &
              mpi_comm=configcomm)
         configcolor=configrank
         call mstm_mpi(mpi_command='split', &
              mpi_color=configcolor,mpi_key=rank, &
              mpi_new_comm=config0comm)

         coherentfix=(n_configurations*nconfiggroups.gt.1)
         coherentfix=.true.
         incinc=(max_number_particles.eq.1)
         incinc=.false.

         if(allocated(acoef)) deallocate(acoef,axtemp)
         allocate(acoef(3,nds2,ml,2),axtemp(ml))

         if(configrank.eq.0) then
            if(allocated(smat)) deallocate(smat,smatave, &
               amean,qabsj,qabsjave)
            nax1=-int((1.d0+k0x)*wx/(2.d0*pi))
            nay1=-int((1.d0+k0y)*wy/(2.d0*pi))
            nax2=int((1.d0-k0x)*wx/(2.d0*pi))
            nay2=int((1.d0-k0y)*wy/(2.d0*pi))
            nax=(nax2-nax1+1)
            nay=(nay2-nay1+1)
            nang=(nax2-nax1+1)*(nay2-nay1+1)
            allocate(smat(16,2,nax1:nax2,nay1:nay2), &
               smatave(16,2,nax1:nax2,nay1:nay2), &
               amean(3,nds2,ml,2), &
               qabsj(ml),qabsjave(ml))
            smatave=0.
            amean=0.
            hfluxave=0.
            cfluxave=0.
            qabsjave=0.
            nconfigave=0
         endif
         if(rank.eq.0) then
            if(allocated(acoh)) deallocate(acoh)
            allocate(acoh(3,ml,2))
            starttime=mstm_mpi_wtime()
            if(.not.print_at_quad_points) then
               open(2,file=output_file,access='append')
               call print_run_variables(2)
               close(2)
            endif
            call print_run_variables(run_print_unit)
         endif

         initialize_constants=.true.
         if(auto_ri_medium) then
            if(auto_update_ri_medium) then
               config0=-1
            else
               config0=1
            endif
         else
            config0=1
         endif

         do config=config0,max(1,ceiling(dble(n_configurations)/dble(nconfiggroups)))
            if(rank.eq.0) then
               write(run_print_unit,'('' configuration(s) '',i4,''-'',i4,''/'',i4)') &
                  1+nconfiggroups*(config-1),nconfiggroups*config, &
                  n_configurations
               call flush(run_print_unit)
            endif
            time1=mstm_mpi_wtime()

            if(allocated(epsilon_cv)) deallocate(epsilon_cv)
            allocate(epsilon_cv(ndpad2,ml))
            epsilon_cv=1.d0

            if(configrank.eq.0) then
               if(homogeneous_layer) then
                  epsilon_cv=ri_binder*ri_binder
               elseif(max_number_particles.le.1) then
                  call fixed_target_epsilon(ndx,ndy,ndsx,ndsy,ml, &
                     d,dw,thick,wx,wy,part_ri(1), &
                     ri_binder,fvc,ntotpart,epsilon_cv)
               else
                  if(write_particle_data.and.rank.eq.0) then
                     writepartdata=.true.
                     open(10,file=particle_data_file)
                     open(11,file='temppartdat.dat')
                  else
                     writepartdata=.false.
                  endif
                  call random_target_epsilon(ndx,ndy,mlsamp,ndsx,ndsy,ml, &
                                      d,dw,fvc,ntotpart,epsilon_cv)
                  if(rank.eq.0) then
                     write(run_print_unit,'('' npart, fv:'',i8,e13.5)') ntotpart,fvc
                     call flush(run_print_unit)
                  endif
                  if(writepartdata) then
                     rewind(11)
                     write(10,'(i10)') ntotpart
                     do j=1,ntotpart
                        read(11,*) partshape,partdata(:)
                        write(10,'(i5,9e13.5)') partshape,partdata
                     enddo
                     close(11,status='delete')
                     close(10)
                  endif
               endif
            endif

            if(numprocs.gt.1) then
               nsend=ml*ndpad2
               call mstm_mpi(mpi_command='bcast', &
                  mpi_rank=0,mpi_number=nsend, &
                  mpi_send_buf_c=epsilon_cv(:,:), &
                  mpi_comm=configcomm)
            endif

            alpha0=sum(epsilon_cv(:,:))/dble(ml*ndpad2)
            if(auto_ri_medium) then
               if(auto_update_ri_medium) then
                  if(config.eq.config0) then
                     rimedium_cv=cdsqrt(dcmplx(dble(alpha0), &
                        abs(dimag(alpha0))))+pcabs*(0.d0,1.d0)
                  elseif(config.le.1) then
                     rimedium_cv=cdsqrt(dcmplx(dble(alpha0), &
                        abs(dimag(alpha0))))+extinction_coefficient*(0.d0,1.d0)/2.d0
                     initialize_constants=.true.
                  endif
               else
                  rimedium_cv=cdsqrt(dcmplx(dble(alpha0), &
                     abs(dimag(alpha0))))+pcabs*(0.d0,1.d0)
               endif
               if(rank.eq.0) then
                  write(run_print_unit,'('' volume mean ri, ri medium:'', 4e13.5)') &
                     cdsqrt(alpha0),rimedium_cv
               endif
            endif

            epsilon_cv=epsilon_cv-rimedium_cv*rimedium_cv

            niter=maxiter
            call ppgf_solver(ndsx,ndsy,ml,niter,solution_eps,acoef,iter, &
               errmax,mpi_comm=configcomm)
            if(configrank.eq.0) then
               call hemispherical_flux(ndsx,ndsy,ml,acoef, &
                  hflux,cflux)
               if(print_scattering_matrix) then
                  call scattered_field(ndsx,ndsy,ml,nax1,nax2,nay1,nay2, &
                     acoef,smat, &
                     incident_rotate=azimuth_average, &
                     alpha_rotate=alpharad, &
                     include_incident=.false.)
               endif
               if(print_absorption_distribution) then
                  call absorption_distribution(ndsx,ndsy,ml,acoef,qabsj)
               endif
            endif
            if(configrank.eq.0) then
               amean=amean+acoef
               hfluxave=hfluxave+hflux
               cfluxave=cfluxave+cflux
               if(print_scattering_matrix) smatave=smatave+smat
               if(print_absorption_distribution) qabsjave=qabsjave+qabsj
            endif

            call mstm_mpi(mpi_command='barrier')

            if(print_scattering_matrix) then
               nsend=16*2*nang
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=smatave, &
                  mpi_recv_buf_dp=smat, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            nsend=3*nds2*ml*2
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dc=amean, &
               mpi_recv_buf_dc=acoef, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=8
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=hfluxave, &
               mpi_recv_buf_dp=hflux, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=2
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=cfluxave, &
               mpi_recv_buf_dp=cflux, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            if(print_absorption_distribution) then
               nsend=ml
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=qabsjave, &
                  mpi_recv_buf_dp=qabsj, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif

            if(rank.eq.0) then
               nconfigave=nconfigave+nconfiggroups
               acoef=acoef/dble(nconfigave)
               acoh(:,:,:)=acoef(:,ss0,:,:)
               if(ml.ge.4) then
                  axtemp=acoef(1,ss0,:,1)
                  call effectiverefractiveindex(3*ml/4-ml/4+1,axtemp(ml/4:3*ml/4),d, &
                      effective_ri(1),e0)
               endif
               hflux=hflux/dble(nconfigave)
               cflux=cflux/dble(nconfigave)
               extinction_coefficient=-dlog(cflux(2))/thick*cosb
               if(print_scattering_matrix) then
                  smat=smat/dble(nconfigave)
                  if(.not.coherent_averaging) then
                     acoef=0.d0
                     acoef(:,ss0,:,:)=acoh(:,:,:)
                  endif
                  if(.not.((max_number_particles.eq.1.and.number_particles_specified))) then
                     allocate(scohmat(16,2,nax1:nax2,nay1:nay2))
                     call scattered_field(ndsx,ndsy,ml,nax1,nax2,nay1,nay2, &
                        acoef,scohmat, &
                        incident_rotate=azimuth_average, &
                        alpha_rotate=alpharad, &
                        include_incident=.false.)
                     smat=smat-scohmat
                     deallocate(scohmat)
                  endif
                  if(enforce_symmetry.and.alpharad.eq.0) then
                     smat((/1,2,5,6,11,12,15,16/),:,:,nay1:-1) &
                        =0.5d0*(smat((/1,2,5,6,11,12,15,16/),:,:,nay1:-1) &
                        +smat((/1,2,5,6,11,12,15,16/),:,:,nay2:1:-1))
                     smat((/1,2,5,6,11,12,15,16/),:,:,nay2:1:-1) &
                        = smat((/1,2,5,6,11,12,15,16/),:,:,nay1:-1)
                     smat((/3,4,7,8,9,10,13,14/),:,:,nay1:-1) &
                        =0.5d0*(smat((/3,4,7,8,9,10,13,14/),:,:,nay1:-1) &
                        -smat((/3,4,7,8,9,10,13,14/),:,:,nay2:1:-1))
                     smat((/3,4,7,8,9,10,13,14/),:,:,nay2:1:-1) &
                        =-smat((/3,4,7,8,9,10,13,14/),:,:,nay1:-1)
                  endif
                  call smatsplinecalc(nax1,nax2,nay1,nay2, &
                     smat(:,:,nax1:nax2,nay1:nay2),k0x,k0y,wx,wy, &
                     spline_order=spline_order, &
                     auto_zero_s12=auto_zero_s12)
               endif

               if(print_absorption_distribution) then
                  qabsj=qabsj/dble(nconfigave)*d_cv/2.d0
               endif

               time2=mstm_mpi_wtime()
               time1=time2-time1
               time2=time2-starttime
               if(.not.print_at_quad_points) then
                  call print_calculation_results(ml,output_file,config,d,wx,wy, &
                     dtheta,fvc,ntotpart,iter,errmax,hflux,cflux, &
                     qabsj,time1,time2)
               else
                  call print_s11_at_quad_points(output_file,hflux,cflux)
               endif
            endif
            first_run=.false.
            extinction(1)=extinction_coefficient
            call mstm_mpi(mpi_command='bcast', &
               mpi_send_buf_dp=extinction, &
               mpi_number=1, &
               mpi_rank=0)
            extinction_coefficient=extinction(1)
            if(configrank.eq.0.and.config.eq.0) then
               smatave=0.
               amean=0.
               hfluxave=0.
               cfluxave=0.
               qabsjave=0.
               nconfigave=0
            endif

            epsilon_cv=epsilon_cv+rimedium_cv*rimedium_cv

         enddo

         call mstm_mpi(mpi_command='free',mpi_comm=configcomm)
         call mstm_mpi(mpi_command='free',mpi_comm=config0comm)
         end subroutine main_calling_program

         subroutine random_target_epsilon(ndx,ndy,mlsamp,ndsx,ndsy,ml, &
             dz,dlat,fvc,ntotpart,epsilonmap, &
             set_target_thickness,ordered_placement)
         use specialfuncs
         use mpidefs
         use target_generation
         implicit none
         logical, optional :: ordered_placement
         logical :: writedipolemap,readdipolemap,writepartdata,oplace
         integer :: ndx,ndy,mlsamp,ndsx,ndsy,ml,ix,iy,j,ntotpart, &
                    rank0,ndt,mlt,ierr,dmapt,partdataunit, &
                    iflag,cld(3),cud(3)
         integer(4), allocatable, save :: dipolemap(:,:,:)
         real(8) :: dz,dlat,fvc,thick,wx,wy
         real(8), optional :: set_target_thickness
         real(8), allocatable :: spherepos(:,:)
         complex(4) :: epsilonmap(min(ndx*ndy,4*ndsx*ndsy),ml)

         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
         writedipolemap=.false.
         readdipolemap=.false.
         configuration_ok=.true.
         if(rank0.eq.0.and.dipole_map_output_file.ne.' ') then
            writedipolemap=.true.
         endif
         if(rank0.eq.0.and.dipole_map_input_file.ne.' ') then
            readdipolemap=.true.
         endif
         if(rank0.eq.0.and.write_particle_data) then
            writepartdata=.true.
         else
            writepartdata=.false.
         endif
         partdataunit=11
         if(present(set_target_thickness)) then
            thick=set_target_thickness
         else
            thick=mlsamp*dz
         endif
         if(present(ordered_placement)) then
            oplace=ordered_placement
         else
            oplace=.true.
         endif
         wx=dlat*dble(ndx)
         wy=dlat*dble(ndy)
         cld(1)=-ndx/2
         cud(1)=(ndx-1)/2
         cld(2)=-ndy/2
         cud(2)=(ndy-1)/2
         cld(3)=-mlsamp/2
         cud(3)=(mlsamp-1)/2

         if(.not.static_configuration) then
            allocate(dipolemap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3)))
         else
            if(run_number.eq.1) allocate(dipolemap(cld(1):cud(1),cld(2):cud(2),cld(3):cud(3)))
         endif
         fvc=0.d0
         core_volume_fraction=0.d0
         if(readdipolemap) then
            ntotpart=0
            dipolemap=0
            open(4,file=dipole_map_input_file)
            read(4,'(2i6)') ndt,ndt,mlt
            do
               read(4,'(4i6)',iostat=ierr) ix,iy,j,dmapt
               if(ierr.ne.0) exit
               dipolemap(ix,iy,j)=dmapt
               fvc=fvc+1.d0
            enddo
            close(4)
         else
            if(read_particle_data) then
               open(11,file=particle_data_input_file)
               read(11,*) ntotpart
               allocate(spherepos(3,ntotpart))
               do j=1,ntotpart
                  read(11,*) spherepos(:,j)
               enddo
               close(11)
               call spherical_particle_dipole_map(ntotpart,spherepos,part_mean_radius(1), &
                   ndx,ndy,ml,dz,dlat,dipolemap)
               deallocate(spherepos)
            else
               if((.not.static_configuration).or.(run_number.eq.1)) then
                  call generate_random_particle_target(n_components,part_mean_radius, &
                     part_fv,part_shape,part_sigma, &
                     dz,dlat,dipolemap,cld,cud,n_monomers,n_aggregates,iflag, &
                     aspect_ratio_y=part_ary, &
                     aspect_ratio_z=part_arz, &
                     boundary_model=boundary_type, &
                     periodic_boundary=(/periodic_sides,periodic_sides,.false./), &
                     shell_thickness=part_shell_thickness, &
                     diffusion_steps=number_diffusion_cycles, &
                     aggregation_frac=aggregation_fraction, &
                     sticking_prob=sticking_probability, &
                     rotation_prob=rotation_probability, &
                     print_unit=run_print_unit)
                  ntotpart=n_monomers
                  if(iflag.ne.0) then
                     write(run_print_unit,'('' warning: particle sampling algorithm did not converge'')')
                     call flush(run_print_unit)
                     configuration_ok=.false.
                  endif
               endif
            endif
            if(writedipolemap) then
               open(4,file=dipole_map_output_file)
               write(4,'(3i6)') ndx,ndy,ml
            endif
            do j=cld(3),cud(3)
               do iy=cld(2),cud(2)
                  do ix=cld(1),cud(1)
                     if(dipolemap(ix,iy,j).ne.0) then
                        fvc=fvc+1.d0
                        if(writedipolemap) then
                           write(4,'(4i6)') ix,iy,j,dipolemap(ix,iy,j)
                        endif
                        if(dipolemap(ix,iy,j).gt.0) then
                           core_volume_fraction=core_volume_fraction+1.d0
                        endif
                     endif
                  enddo
               enddo
            enddo
            if(writedipolemap) close(4)
         endif
         fvc=fvc/dble(ndx*ndy*ml)
         core_volume_fraction=core_volume_fraction/dble(ndx*ndy*ml)
         call epsilon_filter(ndx,ndy,mlsamp,ndsx,ndsy,ml,wx,wy,dipolemap, &
              n_components,part_ri,ri_binder,epsilonmap,part_shell_ri)
         if(.not.static_configuration) deallocate(dipolemap)
         end subroutine random_target_epsilon

         subroutine epsilon_filter(ndx,ndy,mlt,ndsx,ndsy,ml,wx,wy,dipolemap, &
                       ncomp,ri,ribinder,epsilonmap,ri_shell)
         use specialfuncs
         use fft
         implicit none
         logical :: nonabs
         integer :: ndx,ndy,mlt,ml,ncomp, &
                    ndsx,ndsy,nd2,nds2,ndpadx,ndpady,ndpad2,i1x,i2x,is1x,is2x,ip1x,ip2x, &
                    i1y,i2y,is1y,is2y,ip1y,ip2y, &
                    j1target,j2target,j1,j2,jp,iy,ix,j, &
                    comp,sx,sy,smax
         integer(4) :: dipolemap(-ndx/2:(ndx-1)/2,-ndy/2:(ndy-1)/2,-mlt/2:(mlt-1)/2)
         real(8) :: wx,wy,pi
         complex(8) :: ri(ncomp),ribinder
         complex(8), optional :: ri_shell(ncomp)
         complex(4) :: epsilonmap(min(ndx*ndy,4*ndsx*ndsy),ml)
         complex(8), allocatable :: alphaj(:,:),alphaj2(:,:)
         pi=4.d0*datan(1.d0)
         ndpadx=ndsx+ndsx
         ndpadx=min(ndpadx,ndx)
         ndpady=ndsy+ndsy
         ndpady=min(ndpady,ndy)
         nd2=ndx*ndy
         nds2=ndsx*ndsy
         ndpad2=ndpadx*ndpady
         i1x=-ndx/2
         i2x=(ndx-1)/2
         i1y=-ndy/2
         i2y=(ndy-1)/2
         is1x=-ndsx/2
         is2x=(ndsx-1)/2
         is1y=-ndsy/2
         is2y=(ndsy-1)/2
         ip1x=-ndpadx/2
         ip2x=(ndpadx-1)/2
         ip1y=-ndpady/2
         ip2y=(ndpady-1)/2
         smax=ip1x*ip1y

         j1=-floor(dble(ml)/2.d0)
         j2=j1+ml-1
         epsilonmap=ribinder*ribinder
         nonabs=(abs(dimag(ribinder)).le.1.d-8).and.all(abs(dimag(ri(1:ncomp))).le.1.d-8)
         allocate(alphaj(i1x:i2x,i1y:i2y),alphaj2(ip1x:ip2x,ip1y:ip2y))
         j1target=-floor(dble(mlt)/2.d0)
         j2target=j1target+mlt-1
         do j=j1,j2
            if(j.lt.j1target.or.j.gt.j2target) cycle
            alphaj=ribinder*ribinder
            jp=j-j1+1
            do iy=i1y,i2y
               do ix=i1x,i2x
                  comp=dipolemap(ix,iy,j)
                  if(comp.gt.0) then
                     alphaj(ix,iy)=ri(comp)*ri(comp)
                  elseif(comp.lt.0.and.present(ri_shell)) then
                     alphaj(ix,iy)=ri_shell(-comp)*ri_shell(-comp)
                  endif
               enddo
            enddo
            call fft2d(1,ndx,ndy,wx,wy,0.d0,0.d0,alphaj(:,:),alphaj(:,:),-1)
            alphaj2=0.d0
            if(circular_window) then
               do sy=ip1y,ip2y
                  do sx=ip1x,ip2x
                     if(sx*sx+sy*sy.gt.smax) cycle
                     alphaj2(sx,sy)=alphaj(sx,sy)*sqrt(dble(ndpad2)/dble(nd2))
                  enddo
               enddo
            else
               alphaj2=alphaj(ip1x:ip2x,ip1y:ip2y)*sqrt(dble(ndpad2)/dble(nd2))
            endif
            call fft2d(1,ndpadx,ndpady,wx,wy,0.d0,0.d0,alphaj2,alphaj2,1)
            epsilonmap(:,jp)=reshape(alphaj2,(/ndpad2/))
            if(nonabs.and.enforce_zero_absorption) then
               epsilonmap(:,jp)=dble(epsilonmap(:,jp))
            endif
         enddo
         deallocate(alphaj,alphaj2)
         end subroutine epsilon_filter

         integer function epsilon_func(x,y,z,wx,wy,h)
         use specialfuncs
         implicit none
         real(8) :: x,y,z,wx,wy,h,r1,x0,y0,z0,d,pi,sampradius
         data pi/3.141592653589793d0/
         if(max_number_particles.lt.0) then
            d=part_mean_radius(1)*(4.d0*pi/(3.d0*part_fv(1)))**(1.d0/3.d0)
            x0=d*nint(x/d)
            y0=d*nint(y/d)
            z0=d*nint(z/d)
            call psdsamp(part_sigma(1),d/2.d0,sampradius)
            sampradius=sampradius*part_mean_radius(1)
         else
            x0=0.
            y0=0.
            z0=0.
            sampradius=part_mean_radius(1)
         endif
         r1=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
         if((r1.le.sampradius)) then
            epsilon_func=1
         else
            epsilon_func=0
         endif
         end function epsilon_func

         subroutine fixed_target_epsilon(ndx,ndy,ndsx,ndsy,ml,d,dw,h,wx,wy,ri,riback, &
            fvc,ntotpart,epsilonmap)
         use specialfuncs
         implicit none
         integer :: ndx,ndy,ndsx,ndsy,ml,ndpadx,ndpady,ndpad2,nd2,nds2,i1x,i2x,i1y,i2y, &
              sx,sy,j, &
              jp,ntotpart,nocc,ntot,n
         integer(4), allocatable :: dipolemap(:,:,:)
         real(8) :: x,y,z,d,h,wx,wy,dw,pi, &
            fvc
         complex(8) :: ri,ricomp(1),riback
         complex(4) :: epsilonmap(min(ndx*ndy,4*ndsx*ndsy),ml)

         pi=4.d0*datan(1.d0)
         ricomp=ri
         ndpadx=ndsx+ndsx
         ndpadx=min(ndpadx,ndx)
         ndpady=ndsy+ndsy
         ndpady=min(ndpady,ndy)
         nd2=ndx*ndy
         nds2=ndsx*ndsy
         ndpad2=ndpadx*ndpady
         i1x=-ndx/2
         i2x=(ndx-1)/2
         i1y=-ndy/2
         i2y=(ndy-1)/2
         epsilonmap=0.d0

         allocate(dipolemap(ndx,ndy,-ml/2:(ml-1)/2))
         nocc=0
         ntot=0
         do j=1,ml
            z=(dble(j)-.5d0)*d-h/2.d0
            jp=j-1-ml/2
            do sy=i1y,i2y
               do sx=i1x,i2x
                  x=dble(sx)*dw
                  y=dble(sy)*dw
                  n=epsilon_func(x,y,z,wx,wy,h)
                  dipolemap(sx-i1x+1,sy-i1y+1,jp)=n
                  ntot=ntot+1
                  nocc=nocc+n
               enddo
            enddo
         enddo
         fvc=dble(nocc)/dble(ntot)
         ntotpart=nint(dble(nocc)/(4.d0*pi*part_mean_radius(1)/(dw*dw*d)))
         call epsilon_filter(ndx,ndy,ml,ndsx,ndsy,ml,wx,wy,dipolemap, &
                    1,ricomp,riback,epsilonmap)
         end subroutine fixed_target_epsilon

         recursive subroutine nested_loop(looplevel,rank)
         use solver
         implicit none
         logical :: continueloop
         integer :: looplevel,varposition,loopindex,rank
         integer, pointer :: i_loop_var_pointer
         real(8) :: maxdif,loopdif
         real(8), pointer :: r_loop_var_pointer
         complex(8), pointer :: c_loop_var_pointer
         character*128 :: varlabel
         character*1 :: vartype

         varlabel=loop_var_label(looplevel)
         vartype=loop_var_type(looplevel)
         varposition=loop_component_number(looplevel)
         if(vartype.eq.'i') then
            maxdif=abs(i_var_stop(looplevel)-i_var_start(looplevel))
            call variable_list_operation(varlabel, &
               var_position=varposition, &
               i_var_pointer=i_loop_var_pointer)
               i_loop_var_pointer=i_var_start(looplevel)
         elseif(vartype.eq.'r') then
            maxdif=abs(r_var_stop(looplevel)-r_var_start(looplevel))
            call variable_list_operation(varlabel, &
               var_position=varposition, &
               r_var_pointer=r_loop_var_pointer)
               r_loop_var_pointer=r_var_start(looplevel)
         elseif(vartype.eq.'c') then
            maxdif=cdabs(c_var_stop(looplevel)-c_var_start(looplevel))
            call variable_list_operation(varlabel, &
               var_position=varposition, &
               c_var_pointer=c_loop_var_pointer)
               c_loop_var_pointer=c_var_start(looplevel)
         endif

         loopindex=0
         continueloop=.true.
         do while(continueloop)
            loopindex=loopindex+1
            if(looplevel.eq.n_nest_loops) then
               run_number=run_number+1
               call main_calling_program()
            else
               call nested_loop(looplevel+1,rank)
            endif
            if(vartype.eq.'i') then
               i_loop_var_pointer=i_loop_var_pointer+i_var_step(looplevel)
               loopdif=abs(i_loop_var_pointer-i_var_start(looplevel))
            elseif(vartype.eq.'r') then
               r_loop_var_pointer=r_loop_var_pointer+r_var_step(looplevel)
               loopdif=abs(r_loop_var_pointer-r_var_start(looplevel))
            elseif(vartype.eq.'c') then
               c_loop_var_pointer=c_loop_var_pointer+c_var_step(looplevel)
               loopdif=cdabs(c_loop_var_pointer-c_var_start(looplevel))
            endif
!            if(loopindex.gt.1000) exit
            if(loopdif-maxdif.gt.1.d-6) continueloop=.false.
         enddo
         end subroutine nested_loop

         subroutine output_header(iunit,inputfile)
         implicit none
         integer :: iunit
         character*8 :: rundate
         character*10 :: runtime
         character*128 :: inputfile
         call date_and_time(date=rundate,time=runtime)
         run_date_and_time=trim(rundate)//' '//trim(runtime)
         write(iunit,'(''****************************************************'')')
         write(iunit,'(''****************************************************'')')
         write(iunit,'('' pwpp calculation results'')')
         write(iunit,'('' date, time:'')')
         write(iunit,'(a)') run_date_and_time
         write(iunit,'('' input file:'')')
         write(iunit,'(a)') trim(inputfile)
         end subroutine output_header

         subroutine print_run_variables(iunit)
         use intrinsics
         implicit none
         integer :: iunit,i

         write(iunit,'(''****************************************************'')')
         write(iunit,'('' input variables for run '',i5)') run_number
         write(iunit,'('' max_iterations,solution_eps'')')
         write(iunit,'(i10,e13.5)') max_iterations,solution_eps
         write(iunit,'('' incident alpha,beta(deg)'')')
         write(iunit,'(2e13.5)') incident_alpha_deg,incident_beta_deg
         if(auto_ri_medium) then
            if(precon_absorption.lt.0.) then
               write(iunit,'('' ri medium calculated using Mie optical thickness'')')
               write(iunit,'('' absorption factor'')')
               write(iunit,'(e13.5)') -precon_absorption
            else
               write(iunit,'('' ri medium calculated automatically'')')
               write(iunit,'('' absorption factor'')')
               write(iunit,'(e13.5)') precon_absorption
               if(auto_update_ri_medium) then
                  write(iunit,'('' ri medium updated using extinction coefficient'')')
               endif
            endif
         endif
         write(iunit,'('' ndx_width, ndy_width,,dz,dlat,k_cut'')')
         write(iunit,'(2i5,3e13.5)') ndx_width,ndy_width,dipole_spacing,lateral_dipole_spacing,k_cut
         if(ot_specified) then
            write(iunit,'('' Mie optical thickness specified'')')
            write(iunit,'(e13.5)') mie_optical_thickness
         endif
         write(iunit,'('' number components'')')
         write(iunit,'(i5)') n_components
         write(iunit,*) ' comp   shape rad     sigma   ary     arz    fv           re(m)        im(m) '
         do i=1,n_components
            write(iunit,'(2i6,4f8.3,3e13.5)') i,part_shape(i), &
               part_mean_radius(i),part_sigma(i),part_ary(i),part_arz(i), &
               part_fv(i),part_ri(i)
         enddo
         write(iunit,'('' particle shell properties'')')
         write(iunit,*) ' thickness    re(m)        im(m) '
         do i=1,n_components
            write(iunit,'(f8.2,3e13.5)') part_shell_thickness(i),part_shell_ri(i)
         enddo
         if(read_particle_data) then
            write(iunit,'('' sphere positions read from file:'',a)') trim(particle_data_input_file)
         else
            write(iunit,'('' particle positions generated automatically'')')
            write(iunit,'('' maximum diffusion cycles, sticking, rotation probabilities:'')')
            write(iunit,'(i10,2e13.5)') number_diffusion_cycles,sticking_probability,rotation_probability
         endif
         write(iunit,'('' binder refractive index'')')
         write(iunit,'(2e13.5)') ri_binder
         write(iunit,'('' front medium refractive index'')')
         write(iunit,'(2e13.5)') ri_front_surface
         write(iunit,'('' back medium refractive index'')')
         write(iunit,'(2e13.5)') ri_back_surface
         if(boundary_type.eq.0) then
            write(iunit,'('' natural boundary (particles intact)'')')
         else
            write(iunit,'('' flat boundary (particles sliced)'')')
         endif
         write(iunit,'('' ndsx, ndsy, lateral spacing, sample x, y width '')')
         write(iunit,'(2i6,3e13.5)') ndrsx_width,ndrsy_width,set_lateral_dipole_spacing,target_x_width,target_y_width
         write(iunit,'('' nd thickness, vertical spacing, target thickness'')')
         write(iunit,'(i6,2e13.5)') nd_thickness,set_dipole_spacing,target_thickness

         write(iunit,'(''****************************************************'')')
         write(iunit,'('' calculation results for run '')')
         write(iunit,*) run_number
         call flush(iunit)
         end subroutine print_run_variables

         subroutine azimuth_average_scattering_matrix(ntheta,savemat,s11_normalize)
         use bspline_sub_module
         implicit none
         logical, optional :: s11_normalize
         logical :: s11norm
         integer :: itheta,nphi,iphi,ntheta,i
         real(8) :: smat(32),kx,ky,savemat(32,ntheta), &
                    thetadeg,theta,ct,st,phi,cp,sp
         if(present(s11_normalize)) then
            s11norm=s11_normalize
         else
            s11norm=.false.
         endif
         do itheta=1,ntheta
            thetadeg=dble(itheta-1)*90.d0/dble(max(1,ntheta-1))
            if(itheta.eq.1) thetadeg=0.01d0
            theta=pi*thetadeg/180.d0
            ct=cos(theta)
            st=sin(theta)
            nphi=4*max(1,ceiling(st*dble(ntheta)))
            savemat(:,itheta)=0.d0
            do iphi=0,nphi-1
               phi=2.d0*pi*dble(iphi)/dble(nphi)
               cp=cos(phi)
               sp=sin(phi)
               kx=st*cp
               ky=st*sp
               call smatsplineeval(kx,ky,smat)
               savemat(:,itheta)=savemat(:,itheta)+smat
            enddo
            savemat(:,itheta)=savemat(:,itheta)/dble(nphi)
            if(s11norm) then
               if(savemat(1,itheta).gt.1.d-16) then
                  do i=2,16
                     savemat(i,itheta)=savemat(i,itheta)/savemat(1,itheta)
                     savemat(i+16,itheta)=savemat(i+16,itheta)/savemat(17,itheta)
                  enddo
               else
                  savemat(:,itheta)=0.d0
               endif
            endif
         enddo
         end subroutine azimuth_average_scattering_matrix

         subroutine print_calculation_results(ml,fout,config,d,wx,wy, &
            deltak,fvc,ntotpart,niter,errmax,hflux,cflux, &
            qabsj,time1,time2)
         use solver
         use bspline_sub_module
         implicit none
         integer :: outunit,config,ix,iy,i,j,ntotpart,niter, &
                    ntheta,nthetamax,nax,nang,nay,ml
         real(8) :: hflux(2,2,2),cflux(2), &
            d,wx,wy,kx,ky,time1,time2,fvc,savemat(4,4,2,91),dsca,smats(4,4), &
            smatt(4,4), &
            cohref,cohtra,difref,diftra,qabs,qsca, &
            deltak,kr,qabsj(ml),qabst,errmax
         character*128 :: fout,chartemp
         ntheta=91
         if(fout(1:7).eq.'console') then
            outunit=6
         else
            outunit=2
            open(2,file=trim(fout))
            do
               read(2,'(a)') chartemp
               if(trim(chartemp).eq.run_date_and_time) exit
            enddo
            do
               read(2,'(a)') chartemp
               if(chartemp(1:28).eq.' calculation results for run') then
                  read(2,*) i
                  if(i.eq.run_number) exit
               endif
            enddo
         endif
         cohref=cflux(1)
         cohtra=cflux(2)
         difref=sum(hflux(:,:,1))/2.d0-cohref
         diftra=sum(hflux(:,:,2))/2.d0-cohtra
!         difref=sum(hflux(:,:,1))/2.d0
!         diftra=sum(hflux(:,:,2))/2.d0
         if(homogeneous_layer) then
            dsca=cohref+cohtra+difref+diftra
         else
            dsca=difref+diftra
         endif
         if(read_particle_data) then
            n_monomers=ntotpart
            n_aggregates=ntotpart
         endif
         write(outunit,'('' current, total configurations '')')
         write(outunit,'(2i6)') config*n_configuration_groups,n_configurations
         write(outunit,'('' number iterations, error, time per solution, total time '')')
         write(outunit,'(i6,3e13.5)') niter,errmax,time1, time2
         write(outunit,'('' number monomer particles, number aggregates,total, core volume fractions '')')
         write(outunit,'(2i6,2e13.5)') n_monomers,n_aggregates,fvc,core_volume_fraction
         if(.not.configuration_ok) then
            write(outunit,'('' WARNING: particle configuration algorithm did not sucessfully complete!!!'')')
            write(outunit,'('' results might be in serious error!!!'')')
         endif
         write(outunit,'('' DGF refractive index'')')
         write(outunit,'(2e13.5)') rimedium_cv
         write(outunit,'('' effective refractive index'')')
         write(outunit,'(2e13.5)') effective_ri
         if(number_particles_specified.and.max_number_particles.eq.1) then
            write(outunit,'('' scattering, absorption, extinction efficiency'')')
            qabs=(1.d0-difref-diftra-cohref-cohtra)/pi*target_x_width*target_y_width/part_mean_radius(1)**2
            qsca=dsca/pi*target_x_width*target_y_width/part_mean_radius(1)**2
            write(outunit,'(3e13.5)') qsca,qabs,qsca+qabs
         else
            write(outunit,'('' dimensionless extinction coefficient, optical thickness'')')
            write(outunit,'(3e13.5)') extinction_coefficient,extinction_coefficient*target_thickness
         endif
         write(outunit,'('' par-par r    per-par r    par-par t   per-par t '')')
         write(outunit,'(5e13.5)') hflux(1,1,1),hflux(2,1,1),hflux(1,1,2),hflux(2,1,2)
         write(outunit,'('' par-per r    per-per r    par-per t   per-per t '')')
         write(outunit,'(5e13.5)') hflux(1,2,1),hflux(2,2,1),hflux(1,2,2),hflux(2,2,2)
         write(outunit,'('' unp dif r    unp dif t    unp coh r    unp coh t    unp abs '')')
         write(outunit,'(5e13.5)') difref,diftra,cohref,cohtra, &
            1.d0-difref-diftra-cohref-cohtra
         call flush(outunit)

         if(print_absorption_distribution) then
            write(outunit,'('' layer absorption'')')
            qabst=0.
            do j=1,ml
               write(outunit,'(i5,e13.5)') j,qabsj(j)
               qabst=qabst+qabsj(j)
            enddo
            write(outunit,'('' total absorption'')')
            write(outunit,'(e13.5)') qabst
         endif

         if(print_scattering_matrix) then
            write(outunit,'('' number of directions:'')')
            if(k0x_cv.eq.0.d0.and.k0y_cv.eq.0.d0.and.azimuth_average) then
               call azimuth_average_scattering_matrix(ntheta,savemat, &
                  s11_normalize=.true.)
               nthetamax=91
               savemat(1,1,:,:)=savemat(1,1,:,:)/dsca
               write(outunit,*) nthetamax
               write(outunit,'('' normalized diffuse averaged reflection'')')
               write(outunit,*) ' theta       s11          r12          r13          r14          r21 &
                    &         r22          r23          r24          r31          r32          r33          r34 &
                    &         r41          r42          r43          r44'
               do iy=nthetamax,1,-1
                  write(outunit,'(i5,16e13.5)') 181-iy,((savemat(i,j,1,iy),j=1,4),i=1,4)
               enddo
               write(outunit,'('' normalized diffuse averaged transmission'')')
               write(outunit,*) ' theta       s11          r12          r13          r14          r21 &
                    &         r22          r23          r24          r31          r32          r33          r34 &
                    &         r41          r42          r43          r44'
               do iy=1,nthetamax
                  write(outunit,'(i5,16e13.5)') iy-1,((savemat(i,j,2,iy),j=1,4),i=1,4)
               enddo
            else
               nax=floor(1.d0/deltak)
               if(incident_scattering_plane) then
                  nay=0
                  nang=2*nax+1
               else
                  nay=nax
                  nang=0
                  do iy=-nay,nay
                     ky=dble(iy)*deltak
                     do ix=-nax,nax
                        kx=dble(ix)*deltak
                        if(kx*kx+ky*ky.lt.1.d0) then
                           nang=nang+1
                        endif
                     enddo
                  enddo
               endif
               write(outunit,*) nang
               write(outunit,'('' normalized diffuse reflection'')')
               if(incident_scattering_plane) then
                  write(outunit,*) ' theta       s11          r12          r13          r14          r21 &
                       &         r22          r23          r24          r31          r32          r33          r34 &
                       &         r41          r42          r43          r44'
               else
                  write(outunit,*) ' kx           ky           s11          r12          r13          r14          r21 &
                       &         r22          r23          r24          r31          r32          r33          r34 &
                       &         r41          r42          r43          r44'
               endif
               do iy=-nay,nay
                  do ix=-nax,nax
                     if(incident_scattering_plane) then
                        kr=dble(ix)*deltak
                        kx=kr*cos(incident_alpha_deg*pi/180.d0)
                        ky=kr*sin(incident_alpha_deg*pi/180.d0)
                     else
                        ky=dble(iy)*deltak
                        kx=dble(ix)*deltak
                     endif
                     if(kx*kx+ky*ky.lt.1.d0) then
                        call smatsplineeval(kx,ky,smatt,i_start=1,i_end=16)
                        if(smatt(1,1).gt.1.d-16) then
                           smats(:,:)=smatt(1,1)
!                           smats=1.d0
                           smats(1,1)=1.d0
                           smatt(:,:)=smatt(:,:)/smats(:,:)
                        else
                           smatt=0.d0
                        endif
                        smatt(1,1)=smatt(1,1)/dsca
                        if(incident_scattering_plane) then
                           write(outunit,'(17e13.5)') asin(kr)*180.d0/pi,((smatt(i,j),j=1,4),i=1,4)
                        else
                           write(outunit,'(18e13.5)') kx,ky,((smatt(i,j),j=1,4),i=1,4)
                        endif
                     endif
                  enddo
               enddo
               if(dimag(ris_cv).eq.0.d0) then
                  write(outunit,'('' normalized diffuse transmission'')')
                  if(incident_scattering_plane) then
                     write(outunit,*) ' theta       s11          r12          r13          r14          r21 &
                          &         r22          r23          r24          r31          r32          r33          r34 &
                          &         r41          r42          r43          r44'
                  else
                     write(outunit,*) ' kx           ky           s11          r12          r13          r14          r21 &
                          &         r22          r23          r24          r31          r32          r33          r34 &
                          &         r41          r42          r43          r44'
                  endif
                  do iy=-nay,nay
                     do ix=-nax,nax
                        if(incident_scattering_plane) then
                           kr=dble(ix)*deltak
                           kx=kr*cos(incident_alpha_deg*pi/180.d0)
                           ky=kr*sin(incident_alpha_deg*pi/180.d0)
                        else
                           ky=dble(iy)*deltak
                           kx=dble(ix)*deltak
                        endif
                        if(kx*kx+ky*ky.lt.1.d0) then
                           call smatsplineeval(kx,ky,smatt,i_start=17,i_end=32)
                           if(smatt(1,1).gt.1.d-16) then
                              smats(:,:)=smatt(1,1)
!                              smats=1.d0
                              smats(1,1)=1.d0
                              smatt(:,:)=smatt(:,:)/smats(:,:)
                           else
                              smatt=0.d0
                           endif
                           smatt(1,1)=smatt(1,1)/dsca
                           if(incident_scattering_plane) then
                              write(outunit,'(17e13.5)') asin(kr)*180.d0/pi,((smatt(i,j),j=1,4),i=1,4)
                           else
                              write(outunit,'(18e13.5)') kx,ky,((smatt(i,j),j=1,4),i=1,4)
                           endif
                        endif
                     enddo
                  enddo
               endif
            endif
         endif
         if(outunit.ne.6) close(outunit)
         end subroutine print_calculation_results

         subroutine print_s11_at_quad_points(fout,hflux,cflux)
         use solver
         use bspline_sub_module
         implicit none
         integer :: outunit,i,j,ntheta
         real(8) :: hflux(2,2,2),cflux(2),quadangs(40), &
            dsca,smatt(4,4),s11r,s11t,kx,ky,sb, &
            cohref,cohtra,difref,diftra,phid,st
         character*128 :: fout
         data quadangs/0.0572957,3.8802235,7.1012177,10.290330,13.463632,16.623004, &
            19.767435,22.894965,26.003200,29.089485,32.150969,35.184626,38.187255, &
            41.155475,44.085714,46.974200,49.816942,52.609719,55.348064,58.027253, &
            60.642294,63.187919,65.658580,68.048447,70.351418,72.561132,74.670993, &
            76.674202,78.563804,80.332749,81.973959,83.480423,84.845291,86.061987, &
            87.124333,88.026676,88.764010,89.332102,89.727610,89.948236/
         ntheta=40
         if(fout(1:7).eq.'console') then
            outunit=6
         else
            outunit=2
            if(run_number.eq.1) then
               open(2,file=trim(fout))
            else
               open(2,file=trim(fout),access='append')
            endif
         endif
         cohref=cflux(1)
         cohtra=cflux(2)
         difref=sum(hflux(:,:,1))/2.d0-cohref
         diftra=sum(hflux(:,:,2))/2.d0-cohtra
         dsca=difref+diftra
         sb=sin(pi*incident_beta_deg/180.d0)

         write(outunit,'(e13.5)') incident_beta_deg
         write(outunit,'(4e13.5)') difref,cohref,diftra,cohtra
         do i=1,ntheta
            phid=0.d0
            st=sin(pi*quadangs(i)/180.d0)
            do j=1,80
               kx=st*cos(phid*pi/180.d0)+sb
               ky=st*sin(phid*pi/180.d0)
               call smatsplineeval(kx,ky,smatt,i_start=1,i_end=16)
               s11r=smatt(1,1)/dsca
               call smatsplineeval(kx,ky,smatt,i_start=17,i_end=32)
               s11t=smatt(1,1)/dsca
               write(outunit,'(4e13.5)') quadangs(i),phid,s11r,s11t
               phid=phid+4.5d0
            enddo
         enddo
         if(outunit.ne.6) close(outunit)
         end subroutine print_s11_at_quad_points


         subroutine set_string_to_int_variable(sentvarvalue, &
               ivarvalue,var_operation)
         implicit none
         integer :: itemp
         integer, pointer :: ivarvalue
         character*128 :: sentvarvalue,varop,intfile
         character*(*), optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) itemp
         if(varop(1:6).eq.'assign') then
            ivarvalue=itemp
         elseif(varop(1:3).eq.'add') then
            ivarvalue=ivarvalue+itemp
         endif
         end subroutine set_string_to_int_variable

         subroutine set_string_to_real_variable(sentvarvalue, &
               rvarvalue,var_operation)
         implicit none
         real(8) :: rtemp
         real(8), pointer :: rvarvalue
         character*128 :: sentvarvalue,varop,intfile
         character*128, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) rtemp
         if(varop(1:6).eq.'assign') then
            rvarvalue=rtemp
         elseif(varop(1:3).eq.'add') then
            rvarvalue=rvarvalue+rtemp
         endif
         end subroutine set_string_to_real_variable

         subroutine set_string_to_cmplx_variable(sentvarvalue, &
               cvarvalue,var_operation)
         implicit none
         complex(8) :: ctemp
         complex(8), pointer :: cvarvalue
         character*128 :: sentvarvalue,varop,intfile
         character*128, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) ctemp
         if(varop(1:6).eq.'assign') then
            cvarvalue=ctemp
         elseif(varop(1:3).eq.'add') then
            cvarvalue=cvarvalue+ctemp
         endif
         end subroutine set_string_to_cmplx_variable

         subroutine set_string_to_logical_variable(sentvarvalue, &
               lvarvalue,var_operation)
         implicit none
         logical :: ltemp
         logical, pointer :: lvarvalue
         character*128 :: sentvarvalue,varop,intfile
         character*128, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) ltemp
         if(varop(1:6).eq.'assign') then
            lvarvalue=ltemp
         endif
         end subroutine set_string_to_logical_variable

      end module inputinterface
