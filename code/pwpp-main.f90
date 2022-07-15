      use inputinterface
      use solver
      use mpidefs
      use intrinsics
      use specialfuncs
      implicit none
      integer :: looplevel,rank,numprocs,tdum(3),iseed, &
                 readstat(1),i,numargs,istat,numberinputlines,n
      real(4) :: rannum
      character*128 :: inputfile,inputline,oldoutputfile
      character*8 :: rundate
      character*10 :: runtime
      character*128, allocatable :: inputfiledata(:)
      data oldoutputfile/' '/

      call mstm_mpi(mpi_command='init')
      call mstm_mpi(mpi_command='rank',mpi_rank=rank)
      call mstm_mpi(mpi_command='size',mpi_size=numprocs)

      call itime(tdum)
      iseed =3600*tdum(1)+60*tdum(2)+tdum(3)+rank
      rannum=ran3(-iseed)

      do i=0,numprocs-1
         if(rank.eq.i) then
            numargs=mstm_nargs()
            if(numargs.eq.0) then
               inputfile='ppdda.inp'
            else
               call mstm_getarg(inputfile)
            endif
            input_file=inputfile
            open(2,file=inputfile)
            numberinputlines=0
            istat=0
            do while(istat.eq.0)
               read(2,'(a)',iostat=istat) inputline
               numberinputlines=numberinputlines+1
               if(trim(inputline).eq.'end_of_options') exit
            enddo
            if(istat.ne.0) numberinputlines=numberinputlines+1
            allocate(inputfiledata(numberinputlines))
            rewind(2)
            istat=0
            n=0
            do while(istat.eq.0)
               read(2,'(a)',iostat=istat) inputline
               n=n+1
               inputfiledata(n)=inputline
               if(trim(inputline).eq.'end_of_options') exit
            enddo
            if(istat.ne.0) inputfiledata(numberinputlines)='end_of_options'
            close(2)
         endif
         call mstm_mpi(mpi_command='barrier')
      enddo

      repeat_run=.true.
      first_run=.true.

      do while(repeat_run)
         readstat=0
         do i=0,numprocs-1
            if(i.eq.rank) then
               call inputdata(inputfiledata,read_status=readstat(1))
            endif
            call mstm_mpi(mpi_command='barrier')
         enddo
         if(oldoutputfile.ne.output_file) then
            first_run=.true.
            run_number=0
            oldoutputfile=output_file
         endif
         if(rank.eq.0) then
            if(first_run) then
               if(append_output_file) then
                  open(2,file=output_file,access='append')
               else
                  open(2,file=output_file)
                  close(2,status='delete')
                  open(2,file=output_file)
               endif
               call output_header(2,inputfile)
               close(2)
            endif
         endif

         if(n_nest_loops.eq.0) then
            run_number=run_number+1
            call main_calling_program()
         else
            looplevel=1
            call nested_loop(looplevel,rank)
         endif
         n_nest_loops=0
      enddo
      call mstm_mpi(mpi_command='barrier')
      call mstm_mpi(mpi_command='finalize')
      end
