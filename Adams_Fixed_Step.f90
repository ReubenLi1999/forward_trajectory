!***************************************************************************************************
!> author: Jacob Williams
!  date: 4/7/2018
!> revision author: Hao-si Li
!  date: 27/3/2020
!
!  Programme to calculate the velocities and positions integrated by Adams-Method utilizing EGM2008.
!  txt for spherical harmonies with a given vector consisted of velocity and poistion.
!
!  In this programme, the initial date is 2019.03.01, so its corresponding julian day is 2458543.5,
!  and the residual part can be determined by using seconds divided by 86400. Obviously, once the se
!  cond becomes larger than 86399, the next day, 2020.01.02, is coming.
!
!  According to the reference from IERS, the leapsecond for the whole month of 2019.03 is 37.
!
!
!  Date          Version    Programmer    Description
!  ==========    =======    ==========    ===========
!  2020.03.29    1.0        Hao-si Li     None
!
!
!  Input : 1. a set of known velocity and position
!          2. EGM2008.txt
!          3. sofa subroutines to transfer systems offered by Zi-wang Su
!
!  Output: a dataframe as follows:
!          time(s) x(m) y(m) z(k) vx(m/s) vy(m/s) vz(m/s)
!
!  All the original code are from https://github.com/sausage02/Adams-Bashford
!
!***************************************************************************************************

program ddeabm_fixed_step_test

    use ddeabm_kinds
    use ddeabm_module
    use accl_potl_module

    implicit none

    type, extends(ddeabm_with_event_class) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[ddeabm_class]] to include data used in the deriv routine
        real(wp)           :: mu     = 0.0_wp       !! central body gravitational parameter (m3/s2)
        integer(kind = 4)  :: fevals = 0            !! number of function evaluations
        integer(kind = 4)  :: degree = 10
        logical            :: first  = .true.       !! first point is being exported
        logical            :: w_time = .false.
        double precision   :: date   = 2458543.5_wp !! MJD of 2019.03.01
        real(wp)           :: t0     = 0.0_wp     !! initial time (sec)
        real(wp)           :: tf     = 0
        real(wp)           :: dt     = 10.0_wp      !! output step size (sec)
        integer(kind = 4)  :: days   = 2
        real(wp)           :: final_state(6)
    end type spacecraft

    integer,  parameter    :: n      = 6            !! number of state variables
    real(wp), parameter    :: tol    = 1.0e-14_wp   !! event location tolerance

    type(spacecraft)       :: s
    real(wp), dimension(n) :: x0, xf, x, x02
    real(wp)               :: t, z_target
    real(wp)               :: time_begin, time_end
    integer                :: idid
    integer(kind = 8)      :: file_unit
    character(len = 80)    :: log_file
    integer(kind = 4)      :: index

    log_file = 'log.rpt'
    file_unit = 99

    s%mu = 398600.4415_wp * 1000000000  !earth

    x0 = [&
        -6386965.790572938_wp, -2444576.170399037_wp, -688229.5458564862_wp, &
        -665.4527846303586_wp, -403.3112362290047_wp, 7578.348511165726_wp   &
    ]                    !initial state [r, v] (m, m/s)

    s%fevals = 0

    open(file_unit, file=log_file)

    call cpu_time(time_begin)

    write(file_unit, *) ''
    write(file_unit, *) '------------------------'
    write(file_unit, *) ' ddeabm_fixed_step_test'
    write(file_unit, *) '------------------------'
    write(file_unit, *) ''

    !***************************************************************************

    write(*, *) ''
    write(*, *) '-----------------------'
    write(*, *) 'forward the trajectory:'
    write(*, *) '-----------------------'
    write(*, *) ''

    ! initialize the integrator
    call s%initialize(&
        n, maxnum=500000, df=twobody, rtol=[1.0e-14_wp], atol=[1.0e-14_wp], report=twobody_report &
    )
    
    write(file_unit, '(A/, *(F15.6/))') 'Initial time: ', s%t0
    write(file_unit, '(A/, *(F15.6/))') 'Initial state: ', x0

    s%fevals = 0
    s%first  = .true.
    t = s%t0
    x = x0
    s%tf = 86400.0_wp + s%t0 - s%dt    ! final time   (sec)
    do index = 1, s%days, 1
        call s%first_call()
        call s%integrate(t, x, s%tf, idid=idid, integration_mode=2, tstep=s%dt) 
        x = s%final_state
        t = s%t0
        !forward (report points)
    end do
    xf = x
    write(file_unit, *) ''
    write(file_unit, '(A/, *(I5/))')    'idid: ', idid
    write(file_unit, '(A/, *(F15.6/))') 'Final time:', t
    write(file_unit, '(A/, *(F15.6/))') 'Final state:', xf
    write(file_unit, '(A, I10)')        'Function evaluations:',  s%fevals
    write(file_unit, '(A, I5)')         'Days',  s%days
    write(file_unit, *) ''

    !! ----------------------Preserved-----------------------------------------------!!
    ! t = tf
    ! x = xf
    ! s%fevals = 0
    ! call s%first_call()  !restarting the integration
    ! call s%integrate(t, x, t0, idid=idid, integration_mode=2, tstep=dt)  !backwards
    ! x02 = x
    !
    ! write(file_unit, '(A/,*(I5/))')     'idid: ', idid
    ! write(file_unit, '(A/,*(F15.6/))')  'Initial state:', x02
    ! write(file_unit, *)                 ''
    ! write(file_unit, '(A/,*(E20.12/))') 'Error:', x02 - x0
    ! write(file_unit, '(A,I5)')          'Function evaluations:', s%fevals
    ! write(file_unit, *) ''
    !! ----------------------Preserved-----------------------------------------------!!


    call cpu_time(time_end)
    write(*, *) 'This program costs ', time_end - time_begin, ' s.'

    close(file_unit)

    read(*, *)
    stop

contains
!***************************************************************************************************
    !***********************************************************************************************
    subroutine twobody(me, t_in, x_in, xdot)
    !! derivative routine for two-body orbit propagation
        implicit none
        class(ddeabm_class), intent(inout)      :: me
        real(wp)           , intent(in   )      :: t_in
        real(wp)           , intent(in   )      :: x_in(:)
        real(wp)           , intent(  out)      :: xdot(:)
        double precision                        :: residual_date
        integer(kind = 4)                       :: record_num
        
        residual_date = t_in / 86400.0

        select type (me)
        class is (spacecraft)

            record_num = int((me%tf - me%t0) / me%dt + 1)

            gl%r_inertial = x_in(1: 3)
            gl%v_inertial = x_in(4: 6)

            ! in icrf aka inertial system, calculate the effect of three-body (sun and moon)
            ! and the effect of the relativistic effect

            !allocate(cb09(record_num, 4), cb10(record_num, 4))
            !call RelativisticEffect(gl%v_inertial, gl%r_inertial, gl%acc_rela_inertial)
!
            !call readGM(gm_other)
            !call readCB(record_num, cb09, cb10)
            !call calculate_three_body(&
            !    gl%r_inertial, gl%acc_thre_inertial, (t - me%t0) / 10.0 + 1.0, &
            !    cb09, cb10, gm_other, record_num &
            !)
            !deallocate(cb09, cb10)

            ! description about flag
            ! 2 => earth-fixed-system to inertial system
            ! 1 => inertial system to earth-fixed system
            if (me%w_time) then
                call ITRSandGCRS_P(1, me%date, residual_date, 37, gl%r_inertial, gl%r_earth)
            else
                call rotation(t, gl%r_inertial, gl%r_earth, 1)
            end if

            ! calculate the gravitional accelaration in the earth-fixed system
            call accxyz(gl%r_earth, gl%acc_grav_earth, me%degree)

            if (me%w_time) then
                call ITRSandGCRS_P( &
                    2, me%date, residual_date, 37, gl%acc_grav_earth, gl%acc_grav_inertial &
                )
            else
                call rotation(t, gl%acc_grav_earth, gl%acc_grav_inertial, 2)
            end if

            xdot(1: 3) = gl%v_inertial
            xdot(4: 6) = gl%acc_grav_inertial + gl%acc_rela_inertial ! gl%acc_thre_inertial
            me%fevals = me%fevals + 1
        end select
    end subroutine twobody
    !***********************************************************************************************
    !***********************************************************************************************
    subroutine twobody_report(me, t_in, x_in)
        !! report function - write time, state to console
        implicit none
        class(ddeabm_class), intent(inout)     :: me
        real(wp)           , intent(in   )     :: t_in
        real(wp)           , intent(in   )     :: x_in(:)
        real(wp), allocatable                  :: c_coef(:, :)
        real(wp), allocatable                  :: s_coef(:, :)
        character(len = 80)                    :: output_file
        character(len = 80)                    :: earth_fixed_coor
        character(len = 80)                    :: potential_file
        integer(kind = 8)                      :: file_unit_report
        integer(kind = 8)                      :: file_unit_earth
        integer(kind = 8)                      :: file_unit_potential

        file_unit_report = 33
        output_file = 'output\report.txt'
        file_unit_earth = 67
        earth_fixed_coor = 'output\coordinates_earth_fixed_system.txt'
        file_unit_potential = 96
        potential_file = 'output\potential.txt'

        select type (me)
        class is (spacecraft)

            allocate(c_coef(0: me%degree, 0: me%degree), s_coef(0: me%degree, 0: me%degree))

            if (me%first) then  !print header
                open(file_unit_report, file=output_file, status='replace')
                open(file_unit_earth, file=earth_fixed_coor, status='replace')
                open(file_unit_potential, file=potential_file, status='replace')

                write(file_unit_report, '(*(A15, 1X))') 'time (sec)', 'x (m)', 'y (m)', 'z (m)',&
                                                            'vx (m/s)', 'vy (m/s)', 'vz (m/s)'
                write(file_unit_earth, '(*(A15, 1x))') 'x y z'
                write(file_unit_potential, *) 'potential'

                me%first = .false.

                close(file_unit_report)
                close(file_unit_earth)
                close(file_unit_potential)
            end if

            open(file_unit_report, file=output_file, status='old', access='append')
            open(file_unit_earth, file=earth_fixed_coor, status='old', access='append')
            open(file_unit_potential, file=potential_file, status='old', access='append')

            ! Print t and x vector
            write(file_unit_report, '(*(f30.15, 1x))') t_in, x_in
            me%final_state = x_in

            ! print the positions in
            if (me%w_time) then
                call ITRSandGCRS_P(1, me%date, t_in / 86400.0, 37, x_in(1: 3), gl%r_earth)
            else
                call rotation(t_in, x_in(1: 3), gl%r_earth, 1)
            end if
            write(file_unit_earth, '(*(F30.15, 1X))') gl%r_earth

            ! print potential
            call ReadCS(C_coef, S_coef, me%degree)
            call process(C_coef, S_coef, me%degree + 1, gl%r_earth, 0, me%degree, gl%potential)
            write(file_unit_potential, *) gl%potential

            close(file_unit_report)
            close(file_unit_earth)
            close(file_unit_potential)

            deallocate(c_coef)
            deallocate(s_coef)
        end select
    end subroutine twobody_report
    !***********************************************************************************************
    !***********************************************************************************************
    subroutine rotation(t, position_in, position_out, flag)
        implicit none
        !输入
        !t              时间         单位秒
        !x1,y1,z1       坐标         单位米
        !flag           转换方向 flag=1,空转地；flag=2,地转空

        !计算
        !rz          赫尔默特旋转矩阵
        !w           地球自转角速度
        !pi          圆周率
        !az          地球自转的角度

        !输出
        !t              时间         单位秒
        !x2,y2,z2       坐标         单位米

        real(kind = 8)   , intent(in   )       :: t, position_in(3)
        integer(kind = 4), intent(in   )       :: flag
        real(kind = 8)   , intent(  out)       :: position_out(3)

        real(kind = 8) w, pi
        parameter(pi = 3.1415926535_wp)
        real(kind = 8) p1(3), p2(3)
        real(kind = 8) rz(3, 3), az

        !input
        p1 = position_in

        !calculation
        w = 7.2921158553E-5_wp
        az = w * t
        if(flag == 1)then
            call helmert_rz(rz, az)
            p2 = matmul(rz, p1)
        else if(flag == 2)then
            call helmert_rz(rz, -az)
            p2 = matmul(rz, p1)
        end if

        position_out = p2
    end
    !***********************************************************************************************
    !***********************************************************************************************
    subroutine helmert_rz(rz,az)
        implicit none
    !输入
    !az          欧拉角           单位弧度
    !输出
    !rz          赫尔默特旋转矩阵
        real(kind = 8) rz(3, 3), az
        rz(1, 1) = cos(az)
        rz(1, 2) = sin(az)
        rz(1, 3) = 0.0_wp
        rz(2, 1) = -sin(az)
        rz(2, 2) = cos(az)
        rz(2, 3) = 0.0_wp
        rz(3, 1) = 0.0_wp
        rz(3, 2) = 0.0_wp
        rz(3, 3) = 1.0_wp
    end subroutine

end program ddeabm_fixed_step_test
!*****************************************************************************************
