!***************************************************************************************************
!> author: Jacob Williams
!  date: 4/7/2018
!> revision author: Hao-si Li
!  date: 2/3/2020
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
!  2020.03.02    0.0        Hao-si Li     None
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

    implicit none

    type, extends(ddeabm_with_event_class) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[ddeabm_class]] to include data used in the deriv routine
        real(wp)           :: mu     = 0.0_wp    !! central body gravitational parameter (m3/s2)
        integer(kind = 4)  :: fevals = 0         !! number of function evaluations
        integer(kind = 4)  :: degree = 10
        logical            :: first  = .true.    !! first point is being exported
        logical            :: w_time = .true.
        double precision   :: date   = 2458543.5 !! MJD of 2019.03.01
    end type spacecraft

    integer,  parameter    :: n      = 6          !! number of state variables
    real(wp), parameter    :: tol    = 1.0e-14_wp !! event location tolerance
    real(wp), parameter    :: dt     = 10.0_wp    !! output step size (sec)

    type(spacecraft)       :: s
    real(wp), dimension(n) :: x0, xf, x, x02
    real(wp)               :: t0, tf, t, z_target
    real(wp)               :: time_begin, time_end
    integer                :: idid
    integer(kind = 8)      :: file_unit
    character(len = 80)    :: log_file

    log_file = 'log.rpt'
    file_unit = 99

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
    write(*, *) 'forward/backward test:'
    write(*, *) '-----------------------'
    write(*, *) ''

    !constructor (main body is Earth):

    call s%initialize(&
        n, maxnum=100000, df=twobody, rtol=[1.0e-14_wp], atol=[1.0e-14_wp], report=twobody_report &
    )
    s%mu = 398600.4415_wp * 1000000000  !earth
    s%date = 2458543.5

    !initial conditions:
    !x0 = [&
    !    -6476.414366042598150671_wp, 2351.665002495955675840_wp, 9.768305112335665399_wp, &
    !    -665.4527846303586_wp, -403.3112362290047_wp, 7578.348511165726_wp &
    !]                   !initial state [r, v] (km, km/s)

    x0 = [&
        -6386965.790572938_wp, -2444576.170399037_wp, -688229.5458564862_wp, &
        -665.4527846303586_wp, -403.3112362290047_wp, 7578.348511165726_wp   &
    ]                    !initial state [r, v] (m, m/s)

    t0 = 51.18_wp                   ! initial time (sec)
    tf = 86400.0_wp - dt + t0     ! final time   (sec)
    s%fevals = 0

    write(file_unit, '(A/, *(F15.6/))') 'Initial time: ', t0
    write(file_unit, '(A/, *(F15.6/))') 'Initial state: ', x0
    s%fevals = 0
    s%first  = .true.
    t = t0
    x = x0
    call s%first_call()
    call s%integrate(t, x, tf, idid=idid, integration_mode=2, tstep=dt) !forward (report points)
    xf = x
    write(file_unit, *) ''
    write(file_unit, '(A/, *(I5/))')    'idid: ', idid
    write(file_unit, '(A/, *(F15.6/))') 'Final time:', t
    write(file_unit, '(A/, *(F15.6/))') 'Final state:', xf
    write(file_unit, '(A, I5)')         'Function evaluations:',  s%fevals
    write(file_unit, *) ''

    ! t = tf
    ! x = xf
    ! s%fevals = 0
    ! call s%first_call()  !restarting the integration
    ! call s%integrate(t, x, t0, idid=idid, integration_mode=2, tstep=dt)  !backwards
    ! x02 = x

    ! write(file_unit, '(A/,*(I5/))')     'idid: ', idid
    ! write(file_unit, '(A/,*(F15.6/))')  'Initial state:', x02
    ! write(file_unit, *)                 ''
    ! write(file_unit, '(A/,*(E20.12/))') 'Error:', x02 - x0
    ! write(file_unit, '(A,I5)')          'Function evaluations:', s%fevals
    ! write(file_unit, *) ''


    call cpu_time(time_end)
    write(*, *) 'This program costs ', time_end - time_begin, ' s.'

    close(file_unit)

    read(*, *)
    stop

    contains
!*****************************************************************************************

    !*********************************************************
        subroutine twobody(me, t, x, xdot)

        !! derivative routine for two-body orbit propagation

            implicit none

            class(ddeabm_class), intent(inout)      :: me
            real(wp)           , intent(in   )      :: t
            real(wp)           , intent(in   )      :: x(:)
            real(wp)           , intent(  out)      :: xdot(:)

            real(wp), dimension(3)                  :: r_inertial, v

            double precision                        :: r_earth(3) = 0
            double precision                        :: residual_date
            double precision                        :: acc_earth(3)
            double precision                        :: acc_inertial(3)

            residual_date = t / 86400.0

            ! r_earth = [&
            !     5020587.207190956_wp, 4641854.900745415_wp, -699879.8082799137_wp &
            ! ]
            !
            ! call accxyz(r_earth, acc_earth, 10)
            !
            ! write(*, *) acc_earth

            select type (me)
            class is (spacecraft)

                r_inertial = x(1: 3)
                v = x(4: 6)

                ! description about flag
                ! 2 => earth-fixed-system to inertial system
                ! 1 => inertial system to earth-fixed system
                if (me%w_time) then
                    call ITRSandGCRS_P(1, me%date, residual_date, 37, r_inertial, r_earth)
                else
                    call rotation(t, r_inertial, r_earth, 1)
                end if

                ! calculate the gravitional accelaration in the earth-fixed system
                call accxyz(r_earth, acc_earth, me%degree)

                if (me%w_time) then
                    call ITRSandGCRS_P(2, me%date, residual_date, 37, acc_earth, acc_inertial)
                else
                    call rotation(t, acc_earth, acc_inertial, 2)
                end if

                xdot(1: 3) = v
                xdot(4: 6) = acc_inertial

                me%fevals = me%fevals + 1

            end select

        end subroutine twobody
    !*********************************************************

    !*********************************************************
        subroutine twobody_report(me, t, x)

            !! report function - write time, state to console

            implicit none

            class(ddeabm_class), intent(inout)     :: me
            real(wp)           , intent(in   )     :: t
            real(wp)           , intent(in   )     :: x(:)
            real(wp), allocatable                  :: c(:, :)
            real(wp), allocatable                  :: s(:, :)
            real(wp)                               :: x_earth(3)
            real(wp)                               :: acc_earth(3)
            real(wp)                               :: acc_inertial(3)
            real(kind = 8)                         :: V
            character(len = 80)                    :: output_file
            character(len = 80)                    :: earth_fixed_coor
            character(len = 80)                    :: potential_file
            character(len = 80)                    :: accelaration_file
            integer(kind = 8)                      :: file_unit
            integer(kind = 8)                      :: file_unit_earth
            integer(kind = 8)                      :: file_unit_potential
            integer(kind = 8)                      :: file_unit_acc

            file_unit = 33
            output_file = 'output\report.txt'
            file_unit_earth = 67
            earth_fixed_coor = 'output\coordinates_earth_fixed_system.txt'
            file_unit_potential = 96
            potential_file = 'output\potential.txt'
            file_unit_acc = 54
            accelaration_file = 'output\accelaration.txt'

            select type (me)
            class is (spacecraft)

                allocate(c(0: me%degree, 0: me%degree), s(0: me%degree, 0: me%degree))

                if (me%first) then  !print header
                    open(file_unit, file=output_file, status='replace')
                    open(file_unit_earth, file=earth_fixed_coor, status='replace')
                    open(file_unit_acc, file=accelaration_file, status='replace')
                    open(file_unit_potential, file=potential_file, status='replace')

                    write(file_unit, '(*(A15, 1X))')  'time (sec)', 'x (m)', 'y (m)', 'z (m)', &
                                                'vx (m/s)', 'vy (m/s)', 'vz (m/s)'
                    write(file_unit_earth, '(*(A15, 1x))') 'x y z'
                    write(file_unit_potential, *) 'potential'
                    write(file_unit_acc, *) 'ax', 'ay', 'az'

                    me%first = .false.

                    close(file_unit)
                    close(file_unit_acc)
                    close(file_unit_earth)
                    close(file_unit_potential)
                end if

                open(file_unit, file=output_file, status='old', access='append')
                open(file_unit_acc, file=accelaration_file, status='old', access='append')
                open(file_unit_earth, file=earth_fixed_coor, status='old', access='append')
                open(file_unit_potential, file=potential_file, status='old', access='append')

                ! Print t and x vector
                write(file_unit, '(*(f30.15, 1x))') t, x

                ! print the positions in
                if (me%w_time) then
                    call ITRSandGCRS_P(1, me%date, t / 86400.0, 37, x(1: 3), x_earth)
                else
                    call rotation(t, x(1: 3), x_earth, 1)
                end if
                write(file_unit_earth, '(*(F30.15, 1X))') x_earth

                ! print accelaration
                call accxyz(x_earth, acc_earth, me%degree)
                write(file_unit_acc, *) acc_earth

                ! print potential
                call ReadCS(C, S, me%degree)
                call process(C, S, me%degree + 1, x_earth, 0, me%degree, V)
                write(file_unit_potential, *) v
                write(*, *) v

                close(file_unit)
                close(file_unit_acc)
                close(file_unit_earth)
                close(file_unit_potential)
                deallocate(c)
                deallocate(s)
            end select

        end subroutine twobody_report
    !*********************************************************

    subroutine AccXyz(position, acc, degree)

        Implicit None

        integer          , intent(in   )             :: degree
        double  precision, intent(in   )             :: position(3)
        double  precision, intent(  out)             :: acc(3)

        Real(8) :: C(0: degree, 0: degree),         S(0: degree, 0: degree)
        real(8) :: E(0: degree + 1, 0: degree + 1), F(0: degree + 1, 0: degree + 1)
        !C，S为球谐系数，E，F为递推项，VeloX，VeloY，VeloZ为速度三分量，AccX，AccY，AccZ为加速度三分量
        !Acc为加速度总和，CoordinateX，CoordinateY，CoordinateZ为位置三分量，RSatellite为卫星轨道半径
        !TSatellite为卫星坐标对应的时间
        Real(8) AccXx, AccYy, AccZz, AccXS, AccYS, AccZS
        !AccXx,AccYy,AccZz为理论值给定阶次的加速度值，AccXS,AccYS,AccZS为理论值各阶次加速度值和
        Integer i, j
        Real(8) x, y, z
        !计算理论值时的坐标，以及卫星数据采样间隔
        Real(8) :: mu, REarth

        x = position(1)
        y = position(2)
        z = position(3)


        mu     = 398600.4415_wp * 1000000000.0_wp
        REarth = 6.3781363E+6_wp   !地球半径（单位:m）

        !读取球谐系数文件
        Call ReadCS(C, S, degree)

        Call CalculateEF(x, y, z, degree + 1, E, F)

        AccXS = 0
        AccYS = 0
        AccZS = 0

        do i = 0, degree
            do j = 0, i

                if(j == 0) then
                    AccXx = mu / REarth**2 * (-CalBi(i, j, 1)) * E(i + 1, 1) * C(i, 0)
                else
                    AccXx = mu / 2.0_wp / REarth**2 * &
                        ( CalBi(i, j, 2) * (-E(i + 1, j + 1) * C(i, j) - F(i + 1, j + 1) * S(i, j))&
                        + CalBi(i, j, 3) * ( E(i + 1, j - 1) * C(i, j) + F(i + 1, j - 1) * S(i, j)))
                end if

                if(j == 0) then
                    AccYy = mu / REarth**2 * (-CalBi(i, j, 1)) * F(i + 1, 1) * C(i, 0)
                else
                    AccYy = mu / 2.0_wp / REarth**2 * &
                        ( CalBi(i, j, 2) * (-F(i + 1, j + 1) * C(i, j) + E(i + 1, j + 1) * S(i, j))&
                        + CalBi(i, j, 3) * (-F(i + 1, j - 1) * C(i, j) + E(i + 1, j - 1) * S(i, j)))
                end if

                AccZz = mu / REarth**2 * &
                        ( CalBi(i, j, 4) * (-E(i + 1, j) * C(i, j) - F(i + 1, j) * S(i, j)))

                AccXS = AccXS + AccXx
                AccYS = AccYS + AccYy
                AccZS = AccZS + AccZz

            end do
        end do

        acc(1) = accxs
        acc(2) = accys
        acc(3) = acczs

    end subroutine accxyz

!参数b1到b4的计算
    Real(8) Function CalBi(i, j, m)
        integer i, j, m
        real(8) fac
        if(m == 1) then
            fac   = real((i + 1) * (i + 2) * (2 * i + 1))
            CalBi = sqrt(fac / 2.0_wp / real(2 * i + 3))
        else if(m == 2) then
            fac   = real((i + j + 1) * (i + j + 2) * (2 * i + 1))
            calbi = sqrt(fac / real(2 * i + 3))
        else if(m == 3) then
            fac   = real((i - j + 1) * (i - j + 2) * (2 * i + 1))
            if (j == 1) then
                CalBi = sqrt(fac / real(2 * i + 3) * 2.0_wp)
            else
                calbi = sqrt(fac / real(2 * i + 3))
            end if
        else if(m == 4) then
            fac   = real((i - j + 1) * (i + j + 1) * (2 * i + 1))
            CalBi = sqrt(fac / real(2 * i + 3))
        else
            write(*, *) 'm值输入错误', m
        end if

    End Function CalBi

    !读取C和S
    subroutine ReadCS(C, S, n)
        implicit none
        integer line, col, n
        real(8) C(0: n, 0: n), S(0: n, 0: n)
        real(8) c_read((n + 1) * (n + 2) / 2), s_read((n + 1) * (n + 2) / 2)

        integer i

        C = 0
        S = 0

        open(10, file='input/EGM2008.txt', status='old')
        do i = 1, (n + 1) * (n + 2) / 2
            read(10, *) line, col, c_read(i), s_read(i)
            c(line, col) = c_read(i)
            s(line, col) = s_read(i)
        end do
        close(10)
    end subroutine

    subroutine Process(C, S, degree, position, n1, n2, V)
        implicit none
        real(kind = 8)   , intent(in   ) :: position(3)  ! positions in the earth-fixed system
        real(kind = 8)   , intent(  out) :: V            ! gravitional v
        integer(kind = 4), intent(in   ) :: degree       ! the max degree fo the stocks coefficients
        integer(kind = 4), intent(in   ) :: n1, n2       ! the degree span
        real(kind = 8)   , intent(in   ) :: C(0: n2, 0: n2), S(0: n2, 0: n2)

        real(kind = 8), allocatable           :: E(:, :), F(:, :)
        real(kind = 8)                        :: x, y, z

        allocate(e(0: degree, 0: degree), f(0: degree, 0: degree))

        x = position(1)
        y = position(2)
        z = position(3)


        Call CalculateEF(x, y, z, degree, E, F)
        Call CalculateV(E, F, degree, n1, n2, V, c, s)

        deallocate(e)
        deallocate(f)

    end subroutine Process

    subroutine CalculateEF(x, y, z, n, E, F)
        implicit none
        integer(kind = 4), intent(in   ) :: n
        real(kind = 8)   , intent(in   ) :: x, y, z
        real(kind = 8)   , intent(  out) :: E(0: n, 0: n), F(0: n, 0: n)
        integer(kind = 4)                :: i, j
        real(kind = 8)                   :: r, REarth
        real(kind = 8)                   :: coe1, coe2


        REarth = 6.3781363E+6_wp

        r = sqrt(x**2 + y**2 + z**2)

        E(0, 0) = REarth / r
        E(1, 0) = sqrt(3.0_wp) * z * (REarth**2) / (r**3)
        E(1, 1) = sqrt(3.0_wp) * x * (REarth**2) / (r**3)
        F(0, 0) = 0.0_wp
        F(1, 0) = 0.0_wp
        F(1, 1) = sqrt(3.0_wp) * y * (REarth**2) / (r**3)

        do i = 2, n
            do j = 0, i
                if(i == j) then
                    coe1 = sqrt(real(2 * j + 1) / (2 * j))
                    E(i, j) = coe1 * (x * REarth * E(i - 1, j - 1) / (r**2) - &
                                      y * REarth * F(i - 1, j - 1) / (r**2))
                    F(i, j) = coe1 * (x * REarth * F(i - 1, j - 1) / (r**2) + &
                                      y * REarth * E(i - 1, j - 1) / (r**2))
                else
                    coe1 = sqrt(real(2 * i + 1) * (2 * i - 1) / ((i - j) * (i + j)))
                    coe2 = sqrt((real(2 * i + 1) * (i - j - 1) * (i + j - 1)) / &
                                     ((i - j) * (i + j) * (2 * i - 3)))
                    E(i, j) = coe1 * z * REarth * E(i - 1, j) / (r**2) - &
                              coe2 * REarth**2 * E(i - 2, j) / (r**2)
                    F(i, j) = coe1 * z * REarth * F(i - 1, j) / (r**2) - &
                              coe2 * REarth**2 * F(i - 2, j) / (r**2)
                end if
            end do
        end do

    end subroutine

    subroutine CalculateV(E, F, degree, n1, n2, V, c, s)
        implicit none

        integer(kind = 4), intent(in   )             :: n1, n2
        integer(kind = 4), intent(in   )             :: degree
        real(kind = 8)   , intent(  out)             :: V
        real(kind = 8)   , intent(in   )             :: E(0: Degree, 0: Degree)
        real(kind = 8)   , intent(in   )             :: F(0: Degree, 0: Degree)
        real(kind = 8)   , intent(in   )             :: s(0: n2, 0: n2)
        real(kind = 8)   , intent(in   )             :: c(0: n2, 0: n2)

        integer(kind = 4)                            :: i, j
        real(kind = 8)                               :: REarth, MU

        REarth = 6.3781363E+6_wp
        MU = 3.986004415000000e+14_wp
        V = 0.0_wp

        do i = n1, n2
            do j = 0, i
                V = V + E(i, j) * C(i, j) + F(i, j) * S(i, j)
            end do
        end do

        V = V * MU / REarth
    end subroutine

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
