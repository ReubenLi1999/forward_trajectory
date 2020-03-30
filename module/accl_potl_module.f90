module accl_potl_module
    use ddeabm_kinds
    
    implicit none
    type sate_char
        real(wp)                                :: r_inertial(3)
        double precision                        :: r_earth(3) = 0
        real(wp)                                :: v_inertial(3)
        real(wp)                                :: acc_rela_inertial(3)
        real(wp)                                :: acc_thre_inertial(3)
        real(wp)                                :: potential
        double precision                        :: acc_grav_inertial(3)
        double precision                        :: acc_grav_earth(3)
    end type

    real(wp)                                    :: gm_other(10)
    real(wp), allocatable                       :: cb09(:, :)
    real(wp), allocatable                       :: cb10(:, :)
    
    type(sate_char)                             :: gl

contains
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

    subroutine RelativisticEffect(Velocity, Coordinate, Acceleration)
        implicit none
        real(kind = 8), intent(in   )      :: Velocity(1:3)             
        real(kind = 8), intent(in   )      :: Coordinate(1:3)           
        real(kind = 8), intent(  out)      :: Acceleration(1:3)         

        real(kind = 8)      :: REarth=6.371393E+6_wp, Vlight=2.99792458E+8_wp
        real(kind = 8)      :: MU=0.3986004415E+15_wp, JEarth = 9.8e+8_wp

        Acceleration = 0
        call calculate1(Velocity, Coordinate, Acceleration, Vlight,REarth,MU)
        call calculate2(Velocity, Coordinate, Acceleration, Vlight,REarth,MU)
        call calculate3(Velocity,Coordinate,Acceleration,Vlight,REarth,MU,JEarth)

    end subroutine RelativisticEffect


    subroutine calculate2(Velocity,Coordinate,Acceleration,Vlight,REarth,MU)
        implicit none
        real(8) Velocity(1:3), Coordinate(1:3), Acceleration(1:3)
        real(8) Vlight, REarth, MU


        real(8) V,R,coe,coe1
        integer i,j

        V=0
        R=0

        do i=1,3
            V=sqrt(V**2+Velocity(i)**2)
            R=sqrt(R**2+Coordinate(i)**2)
        end do


        coe = 0
        do i = 1, 3
            coe = coe + Velocity(i) * Coordinate(i)
        end do

        coe1 = MU / ((R**3) * (Vlight**2))

        do i = 1, 3
            Acceleration(i) = Acceleration(i) - 4 * coe1 * coe * Velocity(i)
        end do
    end subroutine


    subroutine calculate1(Velocity,Coordinate,Acceleration,Vlight,REarth,MU)
        implicit none
        real(8) Velocity(1:3),Coordinate(1:3),Acceleration(1:3)
        real(8) Vlight, REarth, MU


        integer i
        real(8) coe, V, R


        V = 0
        R = 0
        do i = 1, 3
            V = sqrt(V**2 + Velocity(i)**2)
            R = sqrt(R**2 + Coordinate(i)**2)
        end do

        coe = MU / ((R**3) * (Vlight**2))
        do i = 1, 3
            Acceleration(i) = Acceleration(i) - coe * Coordinate(i) * (4 * MU / R - V**2)
        end do
    end subroutine

    subroutine calculate3(Velocity,Coordinate,Acceleration,Vlight,REarth,MU,JEarth)
        implicit none
        real(8) Velocity(1:3), Coordinate(1:3), Acceleration(1:3)
        real(8) Vlight, REarth, MU
        real(kind = 8) coe, R, JEarth
        integer(kind = 4) i, j

        R = 0
        do i = 1, 3
            R = sqrt(R**2 + Coordinate(i)**2)
        end do

        Acceleration = Acceleration + 2 * MU / (Vlight**2 * R**3) * (3 * Coordinate(3) * JEarth / R**2 + JEarth * (Velocity(2) - Velocity(1)))
    end subroutine
    
    subroutine A_09MOON(record1,GM,position,CB09,a09,number)

        implicit none
        integer record1,i
        real(kind=8) GM(10),position(3),CB09(record1,4),a09(3)
        real(wp) number
        real(kind=8) Lx,Ly,Lz,Ld

        Lx=position(1)*1.0e-3-CB09(number,1)
        Ly=position(2)*1.0e-3-CB09(number,2)
        Lz=position(3)*1.0e-3-CB09(number,3)
        Ld=sqrt(Lx**2+Ly**2+Lz**2)

        a09(1)=-GM(9)*(Lx/Ld**3+CB09(number,1)/CB09(number,4)**3)*1.0e-6
        a09(2)=-GM(9)*(Ly/Ld**3+CB09(number,2)/CB09(number,4)**3)*1.0e-6
        a09(3)=-GM(9)*(Lz/Ld**3+CB09(number,3)/CB09(number,4)**3)*1.0e-6

    end subroutine

    subroutine A_10SUN(record1,GM,position,CB10,a10,number)

        implicit none
        integer record1,i
        real(kind=8) GM(10),position(3),CB10(record1,4),a10(3)
        real(kind=8) Lx,Ly,Lz,Ld
        real(wp) number

        Lx=position(1)*1.0e-3-CB10(number,1)
        Ly=position(2)*1.0e-3-CB10(number,2)
        Lz=position(3)*1.0e-3-CB10(number,3)
        Ld=sqrt(Lx**2+Ly**2+Lz**2)

        a10(1)=-GM(10)*(Lx/Ld**3+CB10(number,1)/CB10(number,4)**3)*1.0e-6
        a10(2)=-GM(10)*(Ly/Ld**3+CB10(number,2)/CB10(number,4)**3)*1.0e-6
        a10(3)=-GM(10)*(Lz/Ld**3+CB10(number,3)/CB10(number,4)**3)*1.0e-6
    
    end subroutine

    subroutine calculate_three_body(position,acc,number,cb09,cb10,GM, record1)
        implicit none
        real(kind=8) position(3),acc(3),acc09(3),acc10(3)
        integer record1
        real(wp) number
        real(kind=8) cb09(record1,4),cb10(record1,4)
        real(kind=8) GM(10)
        
        call A_09MOON(record1,GM,position,CB09,acc09,number)
        call A_10SUN(record1,GM,position,CB10,acc10,number)
        
        acc=acc09+acc10
    
    end subroutine

    subroutine readCB(record1, CB09, CB10)
        implicit none
        integer record1
        real(kind=8) CB09(record1,4),CB10(record1,4)
        
        character(len = 80) CBfile
        character(len = 80) str
        integer i,j,m

        CBfile = 'input/20190301.txt'
        
        
        open(30,file=CBfile,status='old')
        
        do i=1, record1
            read(30,*)
            read(30,*) CB09(i,1),CB09(i,2),CB09(i,3)                    
            CB09(i,4)=sqrt(CB09(i,1)**2+CB09(i,2)**2+CB09(i,3)**2)
            read(30,*)
            read(30,*) CB10(i,1),CB10(i,2),CB10(i,3)
            CB10(i,4)=sqrt(CB10(i,1)**2+CB10(i,2)**2+CB10(i,3)**2)
        end do
        close(30)
    end subroutine

    subroutine readGM(GM)
        implicit none
        
        character(80) str
        character(len = 80) GMfile
        real(kind=8) GM(*)
        real(kind=8) number1
        integer j

        GMfile = 'input/GM.dat'

        open(50,file=GMfile,status='old')
        read(50,*)
        do j=1,10
            read(50,*) str,number1,GM(j)
        end do
        close(50)
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

end module accl_potl_module