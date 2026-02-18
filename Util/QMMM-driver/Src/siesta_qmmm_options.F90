module siesta_qmmm_options
  !! This module has several variables and a single routine in charge of
  !! reading some of the QMMM options.
  implicit none
  private
  public :: read_siesta_qmmm_options
  public :: read_siesta_options

  integer, public :: wricoord    = 1
    !! Write coordinates every wricoord steps.
  integer, public :: mmsteps     = 1
    !! Number of MM steps for every other relaxation step.
  integer, public :: mneigh_freq = 10
    !! Update neighbour lists every mneigh_freq steps.

  logical, public :: launch_siesta_flag    = .false.
    !! Used in fsiesta to indicate whether we launch regular siesta.
  logical, public :: center_qm_system      = .false.
    !! Whether we recenter the center of mass of the QM system as 0.
  logical, public :: siesta_qmmm_usesavexv = .false.
    !! Whether we are restarting from a previous XV file.
  character(len=100), public :: parallel_command = 'mpirun'
    !! Choose the command to run SIESTA in parallel, for cases
    !! when one would like to use mpi_exec or srun.

  logical, public :: opt_idyn         = .false.
    !! Whether we are doing an optimization run.
  logical, public :: fc_idyn          = .false.
    !! Whether we are doing a phonon calculation run.
  logical, public :: therm_idyn       = .false.
    !! Whether we are doing a constant Temperature run.

contains

  subroutine read_siesta_qmmm_options( )
    !! Reads the input QMMM options.
    use fdf, only : fdf_get, leqi
    use sys, only : die

    implicit none
    character(len=80) :: dyntype
    logical           :: default

    ! Write coordinate variable
    wricoord = fdf_get( 'WriteCoordinates', 1 )
    write( 6,'(a,i5,a)' ) &
      'read: Write coordinates each           = ', wricoord, '  steps'

    ! MMxQM steps
    mmsteps = fdf_get( 'MD.MMxQMsteps', 1 )
    write(6,'(a,i5)') &
         'read: MM x QM steps                    = ', mmsteps

    launch_siesta_flag    = fdf_get( 'LaunchSiesta'  , .false. )
    parallel_command      = fdf_get(  'MPICommand'    , 'mpirun')

    ! Options readed here instead of siesta_options
    default               = fdf_get( 'UseSaveData'   , .false. )
    siesta_qmmm_usesavexv = fdf_get( 'MD.UseSaveXV'  , default )
    mneigh_freq           = fdf_get( 'NebListFreq'   , 10      )
    center_qm_system      = fdf_get( 'CenterQmSystem', .false. )
    dyntype               = fdf_get(  'MD.TypeOfRun'  ,'cg'     )

    if ( leqi(dyntype,'cg') .or. leqi(dyntype,'broyden') .or. &
         leqi(dyntype,'fire') ) then
      opt_idyn = .true.
    elseif ( leqi(dyntype,'fc') ) then
      fc_idyn = .true.
    else
      therm_idyn = .true.
    endif

    write(6,'("read: ",73(1h*))')
  end subroutine read_siesta_qmmm_options

  subroutine read_siesta_options( na, ns )
    !! Reads general SIESTA-related options, mainly those partaining to MD.
        ! real*8 charnet           : Net charge (in units of |e|)
    ! real*8 temp              : Temperature for Fermi smearing (Ry)
    ! real*8 etol              : Relative tolerance in CG minimization
    !                            of band structure energy
    ! integer idyn             : Atomic dynamics option:
    !                             0 = Geometry optimization
    !                             1 = Standard MD run (Verlet)
    !                             2 = Nose thermostat MD
    !                             3 = Parrinello-Rahman MD
    !                             4 = Nose thermostat + Parrinello-Rahman MD
    !                             5 = Annealing MD
    !                             6 = Force constants
    !                             8 = Force evaluation
    !                             9 = Explicit set of coordinates
    !                            10 = Lua controlled dynamics
    ! integer istart           : Initial time step for MD
    ! integer ifinal           : Final time step for MD
    ! integer nmove            : Number of steps in *any* MD/optimization
    ! real*8 ftol              : Maximum force for structural optimization
    ! real*8 strtol            : Maximum stress for structural optimization
    ! integer ianneal          : Annealing option for idyn = 5
    !                             1 = Temperature
    !                             2 = Pressure
    !                             3 = Temperature and Pressure
    ! integer iquench          : Quench option: 0 = No; 1 = Yes; 2 = Fire
    ! real*8 dt                : Length of time step (fs)
    ! real*8 dx                : Atomic displacement for Force Constants
    !                             calculation
    ! integer ia1              : First atom to displace for force constants
    ! integer ia2              : Last atom to displace for force constants
    ! real*8 dxmax             : Maximum atomic displacement in one atomic move
    ! real*8 tt                : Target temperature (Kelvin)
    ! real*8 tp                : Target Pressure (Ry/Bohr**3)
    ! real*8 mn                : Mass of Nose variable (Ry/fs**2)
    ! real*8 mpr               : Mass of Parrinello-R. variable (Ry/fs**2)
    ! real*8 bulkm             : Estimate of bulk modulus (Ry/Bohr**3)
    ! real*8 taurelax          : Annealing time to reach targer T and P (fs)
    ! logical usesavecg        : True = try to use continuation CG files
    !                              from disk
    ! real*8 tempinit          : Initial temperature (Kelvin) of the MD simulation
    ! logical varcel           : variable shape for optimization or dynamics
    ! logical broyden_optim    : Use Broyden method for optimization
    use siesta_options
    use precision      , only : dp
    use parallel       , only : IOnode
    use fdf            , only : fdf_get, fdf_defined, fdf_deprecated, leqi
    use files          , only : slabel
    use sys            , only : die
    use units          , only : eV, Ang, Kelvin, GPa
    use siesta_cml     , only : cml_p, mainXML
    use siesta_cml     , only : cmlStartParameterList, cmlEndParameterList, &
                                cmlAddParameter
    use m_target_stress, only : set_target_stress

    implicit none
    integer, intent(in)  :: na
      !! Number of atoms
    integer, intent(in)  :: ns
      !! Number of species

    character(len=22) :: annop, dyntyp
    character(len=20) :: fmt1, fmt3, fmt4, fmt6
    character(len=6)  :: method
    logical           :: qnch, qnch2

    fmt1 = "(a,t53,'= ',2x,l1)"
    fmt3 = "(a,t53,'= ',a)"
    fmt4 = "(a,t53,'= ',i8)"
    fmt6 = "(a,t53,'= ',f10.4,a)"

    if (cml_p) &
      call cmlStartParameterList(mainXML, title='Input Parameters')

    ! for cml output, find the system name & label
    if (cml_p) then
      call cmlAddParameter(xf=mainXML, name='SystemName',             &
                           value=trim(sname), dictref='siesta:sname')
      call cmlAddParameter(xf=mainXML, name='SystemLabel',            &
                           value=trim(slabel), dictref='siesta:slabel')
    endif

    ! Write about Number of species, as before
    if (ionode) write(6,fmt4) 'redata: Number of Atomic Species', ns

    if (ns < 1) &
      call die( 'redata: ERROR: Number of species must be larger than zero.' )

    if (cml_p) &
       call cmlAddParameter( xf=mainXML, title='NumberOfSpecies', &
            value=ns, dictRef='siesta:ns', units="cmlUnits:countable" )

    ! Net charge in the cell-
    charnet = fdf_get('NetCharge',0.0_dp)
    if (ionode) &
       write(6,fmt6) 'redata: Net charge of the system',charnet,' |e|'

    if (cml_p) &
       call cmlAddParameter( xf=mainXML, name='NetCharge', value=charnet, &
            dictRef='siesta:NetCharge', units='siestaUnits:e__')

    ! Use Saved Data
    ! This is duplicated in struct_init because that code requires
    ! UseSaveData to be read before parsing the XV files.
    ! It is duplicated because it can then be used separately.
    usesaveddata = fdf_get('UseSaveData',.false.)
    if (ionode) write(6,fmt1) 'redata: Using Saved Data (generic)', usesaveddata

    ! Dynamics parameters ...
    varcel = fdf_get('MD.VariableCell', .false. )

    ! Type of dynamics
    dyntyp = fdf_get( 'MD.TypeOfRun', 'cg' )
    if (leqi(dyntyp,'cg')) then
      idyn = 0
      usesavecg = fdf_get('MD.UseSaveCG', usesaveddata)
      ! Support the old Broyden switch  for now
      broyden_optim = fdf_get('Optim.Broyden',.false.)
      call fdf_deprecated("Optim.Broyden", "MD.TypeOfRun")

    else if (leqi(dyntyp,'broyden')) then
      idyn = 0
      broyden_optim = .true.
    else if (leqi(dyntyp,'fire')) then
      idyn = 0
      fire_optim = .true.
    else if (leqi(dyntyp,'verlet')) then
      idyn = 1
    else if (leqi(dyntyp,'TDED')) then
      idyn = 1 ! For the time being verlet is used for TDDFT calculations.
      td_elec_dyn = .true.
      rstart_time = 0.00_dp
      totime      = 0.00_dp
    else if (leqi(dyntyp,'nose')) then
      idyn = 2
    else if (leqi(dyntyp,'parrinellorahman')) then
      idyn = 3
    else if (leqi(dyntyp,'noseparrinellorahman')) then
      idyn = 4
    else if (leqi(dyntyp,'anneal')) then
      idyn = 5
    else if (leqi(dyntyp,'fc')) then
      idyn = 6
    else if (leqi(dyntyp,'forces').or.leqi(dyntyp,'master')) then
      idyn = 8
#ifdef NCDF_4
    else if (leqi(dyntyp,'explicit')) then
      idyn = 9
#endif
    else
      call die('Invalid Option selected - value of MD.TypeOfRun not recognised')
    endif

    ! Maximum number of steps in MD/coordinate optimization
    nmove = fdf_get('MD.NumCGsteps',0)
    nmove = fdf_get('MD.Steps',nmove)

    ! Maximum atomic displacement in one step
    dxmax = fdf_get('MD.MaxCGDispl',0.2_dp,'Bohr')
    dxmax = fdf_get('MD.MaxDispl',dxmax,'Bohr')

    ! Tolerance in the maximum atomic force [0.04 eV/Ang]
    ftol = fdf_get('MD.MaxForceTol', 0.00155574_dp, 'Ry/Bohr')

    ! Tolerance in the maximum residual stress (var cell) [1 GPa]
    strtol = fdf_get('MD.MaxStressTol', GPa, 'Ry/Bohr**3')
    strtol = abs(strtol)

    GeometryMustConverge = fdf_get('GeometryMustConverge', .false.)

    if (ionode) then
      select case (idyn)
      case(0)
        if (nmove > 0) then
          if (broyden_optim) then
            write(6,fmt3)'redata: Dynamics option','Broyden coord. optimization'
          elseif (fire_optim) then
            write(6,fmt3) 'redata: Dynamics option', 'FIRE coord. optimization'
          else
            write(6,fmt3) 'redata: Dynamics option','CG coord. optimization'
          endif
          write(6,fmt1) 'redata: Variable cell', varcel
          if (.not. broyden_optim) then
            write(6,fmt1) 'redata: Use continuation files for CG', usesavecg
            write(6,fmt6) 'redata: Max atomic displ per move', dxmax/Ang, ' Ang'
          endif
          write(6,fmt4) 'redata: Maximum number of optimization moves', nmove
          write(6,fmt6) 'redata: Force tolerance', ftol/eV*Ang, ' eV/Ang'
          if (varcel) write(6,fmt6) 'redata: Stress tolerance',strtol/GPa,' GPa'
          if (cml_p) then
            if (broyden_optim) then
              call cmlAddParameter( xf   = mainXML, name = 'MD.TypeOfRun', &
                                    value= 'Broyden' )
            else if (fire_optim) then
              call cmlAddParameter( xf   = mainXML, name = 'MD.TypeOfRun', &
                                    value= 'FIRE' )
            else
              call cmlAddParameter( xf    =mainXML, name  ='MD.TypeOfRun', &
                                    value ='CG' )
              call cmlAddParameter( xf    = mainXML, name  = 'MD.UseSaveCG',&
                                    value = usesavecg )
            endif
            call cmlAddParameter( xf    = mainXML, name  = 'MD.NumCGSteps', &
                                  value = nmove, units = "cmlUnits:countable" )
            call cmlAddParameter( xf    = mainXML, name  = 'MD.Steps', &
                                  value = nmove, units = "cmlUnits:countable" )
            call cmlAddParameter( xf    = mainXML, name  = 'MD.MaxCGDispl',   &
                                  value = dxmax, units = 'siestaUnits:Bohr' )
            call cmlAddParameter( xf    = mainXML, name  = 'MD.MaxDispl',   &
                                  value = dxmax, units = 'siestaUnits:Bohr' )
            call cmlAddParameter( xf=mainXML, name='MD.MaxForceTol',      &
                                  value=ftol, units='siestaUnits:Ry_Bohr' )
            if (varcel) then
              call cmlAddParameter( xf=mainXML, name='MD.MaxStressTol', &
                                  value=strtol, units='siestaUnits:Ry_Bohr__3' )
            endif
          endif
        else
          write(6,fmt3) 'redata: Dynamics option','Single-point calculation'
          if (cml_p) &
            call cmlAddParameter( xf   = mainXML, name = 'MD.TypeOfRun', &
                                  value= 'Single-Point' )
          endif

      case(1)
        write(6,fmt3) 'redata: Dynamics option', 'Verlet MD run'
        if (cml_p) &
          call cmlAddParameter( xf=mainXML, name='MD.TypeOfRun', value='Verlet')

      case(2)
        write(6,fmt3) 'redata: Dynamics option', 'Nose thermostat MD run'
        if (cml_p) &
          call cmlAddParameter( xf=mainXML, name='MD.TypeOfRun', value='Nose')

      case(3)
        write(6,fmt3) 'redata: Dynamics option', 'Parrinello-Rahman MD run'
        if (cml_p) &
          call cmlAddParameter( xf = mainXML, name = 'MD.TypeOfRun', &
                                value = 'Parrinello-Rahman' )

      case(4)
        write(6,fmt3) 'redata: Dynamics option', 'Nose-Parrinello-Rahman MD run'
        if (cml_p) &
          call cmlAddParameter( xf = mainXML, name  = 'MD.TypeOfRun', &
                                value = 'Nose-Parrinello-Rahman' )

      case(5)
        write(6,fmt3) 'redata: Dynamics option', 'Annealing MD run'
        if (cml_p) &
          call cmlAddParameter( xf = mainXML, name  = 'MD.TypeOfRun', &
                                value = 'Annealing' )

      case(6)
        write(6,fmt3) &
          'redata: Dynamics option', 'Force Constants Matrix Calculation'
        if (cml_p) &
          call cmlAddParameter( xf    = mainXML, name  = 'MD.TypeOfRun',   &
                                value = 'Force Constants' )

      case(8)
        write(6,fmt3) 'redata: Dynamics option','Force evaluation'
        if (cml_p) &
          call cmlAddParameter( xf = mainXML, name  = 'MD.TypeOfRun',   &
                                value = 'Force Evaluation' )
#ifdef NCDF_4
      case(9)
        write(6,fmt3) 'redata: Dynamics option','Explicit'
        if (cml_p) &
          call cmlAddParameter( xf    = mainXML, name  = 'MD.TypeOfRun',     &
                                value = 'Explicit' )
#endif
#ifdef SIESTA__FLOOK
      case(10)
        write(6,fmt3) 'redata: Dynamics option','LUA'
        if (cml_p) &
          call cmlAddParameter( xf = mainXML, name  = 'MD.TypeOfRun', &
                                value = 'LUA' )
#endif
      end select
    endif

    ! Initial and final time steps for MD
    istart = fdf_get('MD.InitialTimeStep',1)
    if ( fdf_defined('MD.Steps') ) then
      ifinal = fdf_get('MD.FinalTimeStep',max(1,nmove - istart + 1))
    else
      ifinal = fdf_get('MD.FinalTimeStep',1)
    end if

    ! Length of time step for MD
    dt = fdf_get('MD.LengthTimeStep',1._dp,'fs')

    ! In case of TDDFT, dt is determined from electronic time step.
    if (td_elec_dyn) dt = td_dt * ntded

    ! Quench Option
    qnch  = fdf_get('MD.Quench',.false.)
    qnch2 = fdf_get('MD.FireQuench',.false.)
    if ((qnch .or. qnch2) .and. (idyn==2 .or. idyn==4)) &
      call die( 'redata: ERROR: You cannot quench and '//&
                'use a Nose thermostat simultaneously')

    iquench = 0
    if (qnch)  iquench = 1
    if (qnch2) iquench = 2

    ! Initial Temperature of MD simulation
    ! (draws random velocities from the Maxwell-Boltzmann distribition
    !  at the given temperature)
    tempinit = fdf_get('MD.InitialTemperature',0.0_dp,'K')

    if ( (idyn >= 1) .and. (idyn <= 5) ) then
      if ( ionode ) then
        write(6,fmt4) 'redata: Initial MD time step', istart
        write(6,fmt4) 'redata:   Final MD time step', ifinal
        write(6,fmt6) 'redata: Length of MD time step', dt, ' fs'
        write(6,fmt6) 'redata: Initial Temperature of MD run', tempinit, ' K'
        if ( idyn /= 5 ) then
          if ( qnch2 ) then
            write(6,fmt1) 'redata: Perform a MD Fire quench',qnch2
          else
            write(6,fmt1) 'redata: Perform a MD quench',qnch
          endif
        endif
      endif

      if (cml_p) then
        call cmlAddParameter( xf    = mainXML, name  = 'MD.InitialTimeStep',&
                             value = istart, units = 'cmlUnits:countable' )
        call cmlAddParameter( xf    = mainXML, name  = 'MD.FinalTimeStep',  &
                              value = ifinal, units = 'cmlUnits:countable' )
        call cmlAddParameter( xf=mainXML, name='MD.LengthTimeStep',&
                              value=dt, units='siestaUnits:fs' )
        call cmlAddParameter( xf=mainXML, name='MD.InitialTemperature',&
                              value=tempinit, units='siestaUnits:K' )
        if ( idyn /= 5 ) then
          if ( qnch2 ) then
            call cmlAddParameter(xf=mainXML, name='MD.FireQuench', value=qnch2)
          else
            call cmlAddParameter(xf=mainXML, name='MD.Quench', value=qnch)
          endif
        endif
      endif
    endif

    ! Target Temperature and Pressure
    tt = fdf_get('MD.TargetTemperature',0.0_dp,'K')

    call fdf_deprecated("MD.TargetPressure", "Target.Pressure")
    tp = fdf_get('MD.TargetPressure',0.0_dp,'Ry/Bohr**3')
    tp = fdf_get('Target.Pressure',tp,'Ry/Bohr**3')
    call set_target_stress( tp )

    ! Mass of Nose variable
    mn = fdf_get('MD.NoseMass',100._dp,'Ry*fs**2')

    ! Mass of Parrinello-Rahman variables
    mpr = fdf_get('MD.ParrinelloRahmanMass',100._dp,'Ry*fs**2')

    if ( (idyn == 2) .or. (idyn == 4) ) then
      if ( ionode ) write(6,fmt6) 'redata: Nose mass',mn/eV,' eV*fs**2'
      if ( cml_p ) &
        call cmlAddParameter( xf    = mainXML, name  = 'MD.NoseMass',        &
                              value = mn, units = 'siestaUnits:Ry_fs__2' )
    endif

    if ( (idyn == 3) .or. (idyn == 4) ) then
      if (ionode) &
        write(6,fmt6) 'redata: Parrinello-Rahman mass',mpr/eV,' eV*fs**2'

      if (cml_p) &
        call cmlAddParameter( xf = mainXML, name = 'MD.ParrinelloRahmanMass', &
                              value = mpr, units = 'siestaUnits:Ry_fs__2' )
    endif

    ! Annealing option
    ianneal = 0
    annop   = fdf_get( 'MD.AnnealOption','TemperatureAndPressure' )

    if (idyn == 5) then
      if (leqi(annop,'Temperature')) then
        ianneal = 1
      else if (leqi(annop,'Pressure')) then
        ianneal = 2
      else if (leqi(annop,'TemperatureAndPressure')) then
        ianneal = 3
      else
        call die( 'redata: ERROR: With annealing MD, you must '//&
                  'choose an appropriate value for MD.AnnealOption' )
      endif

      if ( ionode ) then
        select case (ianneal)
          case(1)
            write(6,fmt3) 'redata: Annealing Option', 'Temperature'
            if (cml_p) &
              call cmlAddParameter( xf = mainXML, name  = 'MD.AnnealOption', &
                                    value = 'Temperature' )

          case(2)
            write(6,fmt3) 'redata: Annealing Option', 'Pressure'
            if (cml_p) &
              call cmlAddParameter( xf = mainXML, name  = 'MD.AnnealOption', &
                                    value = 'Pressure' )

          case(3)
            write(6,fmt3) 'redata: Annealing Option', 'Temperature and Pressure'
            if (cml_p) &
              call cmlAddParameter( xf = mainXML, name  = 'MD.AnnealOption',   &
                                    value = 'TemperatureAndPressure')
        end select
      endif
    endif

    if ( (idyn == 2) .or. (idyn == 4) .or. &
         ( (idyn == 5) .and. ((ianneal == 1) .or. (ianneal == 3)) ) ) then
      if (ionode) write(6,fmt6) 'redata: Target Temperature', tt, ' Kelvin'

      if (cml_p) &
        call cmlAddParameter( xf    = mainXML, name  = 'MD.TargetTemperature', &
                              value = tt, units = 'siestaUnits:K' )
    endif

    if ( (idyn == 3) .or. (idyn == 4) .or. &
         ( (idyn == 5) .and. ((ianneal == 2) .or. (ianneal == 3)) ) ) then
      if (ionode) &
        write(6,fmt6) 'redata: Target Pressure', tp/eV*Ang**3, ' eV/Ang**3'

      if (cml_p) &
        call cmlAddParameter( xf=mainXML, name= 'MD.TargetPressure',      &
                              value= tp, units= 'siestaUnits:Ry_Bohr__3' )
    endif

    ! Relaxation Time for Annealing
    taurelax = fdf_get( 'MD.TauRelax', 100._dp, 'fs' )
    if ( idyn == 5 ) then
      if (ionode) &
        write(6,fmt6) 'redata: Annealing Relaxation Time', taurelax,' fs'

      if (cml_p) &
        call cmlAddParameter( xf    = mainXML, name  = 'MD.TauRelax',    &
                              value = taurelax, units = 'siestaUnits:fs' )
    endif

    ! Estimated Bulk modulus (for Pressure annealing) [100 GPa]
    bulkm = fdf_get( 'MD.BulkModulus', 100*Gpa, 'Ry/Bohr**3' )
    if ( ionode ) then
      if ( (idyn == 5) .and. ( (ianneal == 2) .or. (ianneal == 3) ) ) &
        write(6,fmt6) &
          'redata: Approx. Bulk Modulus ', bulkm/eV*Ang**3, ' eV/Ang**3'
    endif

    if (cml_p) &
      call cmlAddParameter( xf = mainXML, name  = 'MD.BulkModulus',  &
                            value = bulkm, units = 'siestaUnits:Ry_Bohr__3' )

    ! Atomic displacement for force constant calculation
    call fdf_deprecated("MD.FCDispl", "FC.Displacement")
    dx = fdf_get('MD.FCDispl',0.04_dp,'Bohr')
    dx = fdf_get('FC.Displacement ',dx,'Bohr')

    ! First and last atoms to displace for calculation of force constants
    call fdf_deprecated("MD.FCFirst", "FC.First")
    ia1 = fdf_get('MD.FCfirst',1) ! Catch old-style keyword (prefer new key)
    ia1 = fdf_get('FC.first',ia1)
    call fdf_deprecated("MD.FCLast", "FC.Last")
    ia2 = fdf_get('MD.FClast',na) ! Catch old-style keyword (prefer new key)
    ia2 = fdf_get('FC.last',ia2)

    ! Check that last atom doesn't exceed total number
    if ( (idyn == 6) .and. (ia2 > na) ) &
      call die( 'redata: ERROR:'//&
                'Last atom index for FC calculation is > number of atoms.' )

    if (idyn == 6) then
      if ( ionode ) then
        write(6,fmt6) 'redata: Atomic displ for force constants', dx/Ang,'  Ang'
        write(6,fmt4) 'redata: First atom to move', ia1
        write(6,fmt4) 'redata: Last atom to move' , ia2
      endif
      if ( cml_p ) then
        call cmlAddParameter( xf    = mainXML, name  = 'MD.FCDispl', &
                              value = dx, units = 'siestaUnits:Bohr' )
        call cmlAddParameter( xf    = mainXML, name  = 'FC.Displacement', &
                              value = dx, units = 'siestaUnits:Bohr' )
        call cmlAddParameter( xf= mainXML, name= 'MD.FCFirst', &
                              value= ia1, units= 'cmlUnits:countable' )
        call cmlAddParameter( xf= mainXML, name= 'FC.First', &
                              value= ia1, units= 'cmlUnits:countable' )
        call cmlAddParameter( xf= mainXML, name= 'MD.FCLast', &
                              value= ia2, units= 'cmlUnits:countable' )
        call cmlAddParameter( xf= mainXML, name= 'FC.Last', &
                              value= ia2, units= 'cmlUnits:countable' )
      endif
    endif

    ! Variable cell shape? Depending on input and type of dynamics
    varcel = varcel .or. (idyn==3) .or. (idyn==4)  &
         .or. (idyn==5 .and. ianneal==1)           &
         .and. (idyn/=1) .and. (idyn/=2)           &
         .and. (idyn/=6) .and. (idyn/=7)           &
         .and. (idyn/=9) &
         .and. (.not. (idyn==5 .and. ianneal/=1) )

    ! First read in, then later correct depending on
    ! the other usages
    use_aux_cell = fdf_get('ForceAuxCell', .false.)

    ! Find some switches
    writek = fdf_get( 'Write.Kpoints', .true. )
    writeF = fdf_get( 'Write.Forces', .true. )
    writec = fdf_get( 'WriteCoorStep', .false. )
    writmd = fdf_get( 'WriteMDhistory', .false. )
    writpx = fdf_get( 'WriteMDXmol', .not. writec )

    ! If no mulliken is requested, set it to false
    rijmin = fdf_get( 'WarningMinimumAtomicDistance', 1.0_dp, 'Bohr' )

    RelaxCellOnly                = fdf_get('MD.RelaxCellOnly', .false.)
    RemoveIntraMolecularPressure = fdf_get( &
         'MD.RemoveIntraMolecularPressure', .false.)

    if ( ionode ) write(6,'(2a)') 'redata: ', repeat('*', 71)
    if ( cml_p ) call cmlEndParameterList( mainXML )
  end subroutine read_siesta_options
end module siesta_qmmm_options

