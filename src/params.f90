! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile,Insect)
    use vars
    use insect_module
    use module_ini_files_parser_mpi
    ! The file we read the PARAMS from
    character(len=strlen),intent(in) :: paramsfile
    ! the insect we initialize here
    type(diptera), intent(inout) :: Insect
    integer :: i
    ! this array contains the entire ascii-params file
    type(inifile):: PARAMS
    real(kind=pr), dimension(1:3) :: defaultvec
    character(len=strlen) :: old_meanflow


    call read_ini_file_mpi(PARAMS, paramsfile, .true.)



    ! Resolution section
    call read_param_mpi(PARAMS,"Resolution","nx",nx, 4)
    call read_param_mpi(PARAMS,"Resolution","ny",ny, 4)
    call read_param_mpi(PARAMS,"Resolution","nz",nz, 4)

    ! Geometry section
    call read_param_mpi(PARAMS,"Geometry","xl",xl, 1.d0)
    call read_param_mpi(PARAMS,"Geometry","yl",yl, 1.d0)
    call read_param_mpi(PARAMS,"Geometry","zl",zl, 1.d0)

    ! lattice spacing is global (since we allow to specify reals in multiples of
    ! grid points, we nedd that value now.)
    dx=xl/dble(nx)
    dy=yl/dble(ny)
    dz=zl/dble(nz)
    pi=4.d0 *datan(1.d0)
    ! scaling for FFTs
    scalex=2.d0*pi/xl
    scaley=2.d0*pi/yl
    scalez=2.d0*pi/zl

    if (nx==1) then
        if (root) write(*,*) "2D run: setting x coordinate accordingly (OVERWRITE!!!)"
        dx = 1.d0
        xl = 1.d0
        if (root) write(*,'("xl=",es12.4," dx=",es12.4)') xl,dx
    endif

    ! Geometry section
    call read_param_mpi(PARAMS,"Geometry","Size",length, 0.d0)

    !-----------------------------------------------------------------------------
    ! Time section
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"Time","nt",nt, 9999999)
    call read_param_mpi(PARAMS,"Time","iTimeMethodFluid",iTimeMethodFluid,"RK4")
    call read_param_mpi(PARAMS,"Time","Tmax",Tmax,1.d9)
    call read_param_mpi(PARAMS,"Time","CFL",cfl,1.d0)
    call read_param_mpi(PARAMS,"Time","dt_max",dt_max,0.d0)
    call read_param_mpi(PARAMS,"Time","dt_fixed",dt_fixed,0.d0)

    !-----------------------------------------------------------------------------
    ! Reynolds number section:
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"ReynoldsNumber","nu",nu,1.d-2)

    !-----------------------------------------------------------------------------
    ! Initial conditions section
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"InitialCondition","inicond",inicond, "none")
    ! if inicond reading from file, which files?
    if (inicond=="infile") then
        call read_param_mpi(PARAMS,"InitialCondition","file_ux",file_ux, "none")
        call read_param_mpi(PARAMS,"InitialCondition","file_uy",file_uy, "none")
        call read_param_mpi(PARAMS,"InitialCondition","file_uz",file_uz, "none")
        call read_param_mpi(PARAMS,"InitialCondition","file_p",file_p, "none")
    endif

    !-----------------------------------------------------------------------------
    ! Saving section
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"Saving","iDoBackup",iDoBackup, 1)
    call read_param_mpi(PARAMS,"Saving","iSaveVelocity",iSaveVelocity, 0)
    call read_param_mpi(PARAMS,"Saving","iSavePress",iSavePress, 0)
    call read_param_mpi(PARAMS,"Saving","iSaveVorticity",iSaveVorticity, 0)
    call read_param_mpi(PARAMS,"Saving","iSaveMask",iSaveMask, 0)
    call read_param_mpi(PARAMS,"Saving","iSaveXMF",iSaveXMF, 0) ! default is no
    call read_param_mpi(PARAMS,"Saving","tsave",tsave, 9.d9)
    call read_param_mpi(PARAMS,"Saving","truntime",truntime, 1.d0)
    call read_param_mpi(PARAMS,"Saving","wtimemax",wtimemax, 8760.d0) ! 1 year
    call read_param_mpi(PARAMS,"Saving","tintegral",tintegral,0.01d0)
    call read_param_mpi(PARAMS,"Saving","tsave_first",tsave_first,0.0d0)
    call read_param_mpi(PARAMS,"Saving","tsave_period",tsave_period,1.0d0)
    call read_param_mpi(PARAMS,"Saving","save_only_one_period", save_only_one_period,"no")
    call read_param_mpi(PARAMS,"Saving","itdrag",itdrag,99999)
    call read_param_mpi(PARAMS,"Saving","itbeam",itbeam,99999)

    !-- dry run, just the mask function
    call read_param_mpi(PARAMS,"DryRun","dry_run_without_fluid",dry_run_without_fluid,"no")
    if (dry_run_without_fluid=="yes") then
        write(*,*) "Attention! This is a dry run without fluid"
        write(*,*) "Deactivating all useless save-switches..."
        idobackup=0
        iSavePress=0
        iSaveVelocity=0
        iSaveVorticity=0
    endif

    !-----------------------------------------------------------------------------
    ! penalization / cavity
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"Penalization","iPenalization",iPenalization, 0)
    call read_param_mpi(PARAMS,"Penalization","iMoving",iMoving, 0)
    call read_param_mpi(PARAMS,"Penalization","iMask",iMask, "none")
    call read_param_mpi(PARAMS,"Penalization","eps",eps, 1.d-2)
    call read_param_mpi(PARAMS,"Penalization","iCavity",iCavity,"no")
    call read_param_mpi(PARAMS,"Penalization","cavity_size",cavity_size,0)
    call read_param_mpi(PARAMS,"Penalization","compute_forces",compute_forces,1)
    call read_param_mpi(PARAMS,"Penalization","unst_corrections",unst_corrections,0)
    call read_param_mpi(PARAMS,"Penalization","iChannel",iChannel,"no")
    if (iChannel=="0") iChannel="no" ! for downward compatibility with older ini files
    if (iChannel=="1") iChannel="xy" ! for downward compatibility with older ini files
    call read_param_mpi(PARAMS,"Penalization","thick_wall",thick_wall,0.2d0)
    call read_param_mpi(PARAMS,"Penalization","pos_wall",pos_wall,0.3d0)

    !-----------------------------------------------------------------------------
    ! Geometry section
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"Geometry","x0",x0, 0.d0)
    call read_param_mpi(PARAMS,"Geometry","y0",y0, 0.d0)
    call read_param_mpi(PARAMS,"Geometry","z0",z0, 0.d0)

    !-----------------------------------------------------------------------------
    ! Saving section
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"Saving","iSaveSolidVelocity",iSaveSolidVelocity,0)

    !-----------------------------------------------------------------------------
    ! MeanFlow section
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow_x",iMeanFlow_x,"free")
    call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow_y",iMeanFlow_y,"free")
    call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow_z",iMeanFlow_z,"free")
    ! for compatibility with old params files:
    call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow",old_meanflow,"unused")
    call read_param_mpi(PARAMS,"MeanFlow","m_fluid",m_fluid, 999.d9)
    call read_param_mpi(PARAMS,"MeanFlow","ux",uxmean, 1.d0)
    call read_param_mpi(PARAMS,"MeanFlow","uy",uymean, 1.d0)
    call read_param_mpi(PARAMS,"MeanFlow","uz",uzmean, 1.d0)
    call read_param_mpi(PARAMS,"MeanFlow","iMeanFlowStartupConditioner", iMeanFlowStartupConditioner,"no")
    call read_param_mpi(PARAMS,"MeanFlow","tau_meanflow",tau_meanflow, 0.d0)
    call read_param_mpi(PARAMS,"MeanFlow","T_release_meanflow",T_release_meanflow,0.d0)

    ! for compatibility with old files:
    if (old_meanflow=="1") then
        iMeanFlow_x="fixed"
        iMeanFlow_y="fixed"
        iMeanFlow_z="fixed"
    endif


    !-----------------------------------------------------------------------------
    ! Incompressibility
    !-----------------------------------------------------------------------------
    call read_param_mpi(PARAMS,"Incompressibility","c_0", c_0, 0.d0)
    call read_param_mpi(PARAMS,"Incompressibility","gamma_p", gamma_p, 1.d0)
    call read_param_mpi(PARAMS,"Incompressibility","method", method,"centered_2nd")


    call clean_ini_file_mpi(PARAMS)

    !-----------------------------------------------------------------------------
    ! DONE..
    !-----------------------------------------------------------------------------
    if (mpirank==0) then
        write (*,*) "*************************************************"
        write (*,'(A,i3)') " *** DONE READING PARAMETERS"
        write (*,*) "*************************************************"
    endif
end subroutine
