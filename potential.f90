program main 
    implicit none
    character(len=200)          ::FilePar
    character(len=200)          ::FileOut
    character(len=200)          ::FileCoefficient
    character(len=200)          ::FileCoordinateAndVelocity
    character(len=100)          ::str
    
    real(kind=8),allocatable    ::C(:,:), S(:,:)
    integer(kind=4)             ::degree
    integer(kind=4)             ::N
    integer(kind=4)             ::flag
    integer(kind=4)             ::ios=0
    real(kind=8)                ::JulianDay
    real(kind=8)                ::Second
    real(kind=8)                ::x, y, z
    real(kind=8)                ::x_2, y_2, z_2
    real(kind=8)                ::VelocityX, VelocityY, VelocityZ
    real(kind=8)                ::T1
    real(kind=8)                ::T2
    real(kind=8)                ::ConstE
    real(kind=8)                ::P , EK, P2

    FilePar='par.txt'
    ! FileCoordinateAndVelocity='20200101statevec_Earth_gcrs_SCA.txt'
    FileCoordinateAndVelocity='output/report.txt'
    call ReadPar(FilePar, FileCoefficient, degree, N, FileOut)
    allocate(C(0:degree,0:degree), S(0:degree,0:degree))
    call ReadSphericalHarmonicCoefficients(FileCoefficient, C, S, degree)


    open(100, file=FileCoordinateAndVelocity, status='old')
    read(100,*)
    ! read(100,*)
    ! read(100,*)
    ! read(100,*)

    open(110, file='output/potential.txt', status='unknown')
    write(110, *) 'potential'
    
    open(120, file='output/coordinates_earth_fixed_system.txt', status='old')
    read(120,*)


    flag = 0
    write(*,*) ios
    do while(ios==0)

        read(100,*,iostat=ios) Second, x, y, z, VelocityX, VelocityY, VelocityZ
        ! read(100,*,iostat=ios) str, Second, x, y, z
        ! read(100,*,iostat=ios) str, Second, VelocityX, VelocityY, VelocityZ
        ! read(100,*,iostat=ios) 
        ! x = x*1e+3
        ! y = y*1e+3
        ! z = z*1e+3
        ! VelocityX = VelocityX*1e+3
        ! VelocityY = VelocityY*1e+3
        ! VelocityZ = VelocityZ*1e+3
        if(ios/=0) exit
        read(120,*,iostat=ios) x_2, y_2, z_2
        ! x_2 = x_2*1e+3
        ! y_2 = y_2*1e+3
        ! z_2 = z_2*1e+3
        ! x_2 = x
        ! y_2 = y
        ! z_2 = z

        if(ios/=0) exit
        call Calculate(flag, x, y, z, x_2, y_2, z_2, VelocityX, VelocityY, VelocityZ, C, S, degree, N, T1, T2, ConstE, P , EK, P2)
        
        write(110,'(f30.15)') p2
        
        flag = 1
    end do

    close(100)
    close(110)
    close(120)

end program

subroutine ReadSphericalHarmonicCoefficients(FileCoefficient,C,S,n)
    implicit none
    character(len=200)  ::FileCoefficient
    integer(kind=4)     ::n
    real(kind=8)        ::C(0:n,0:n),S(0:n,0:n)

    character(len=200)  ::str
    integer(kind=4)     ::i, j

    C=0
    S=0

    open(10,file='input\EGM2008.txt',status='old')
    read(10,*) str
    do i=0,n
        do j=0,i
            read(10,*) str,str,C(i,j),S(i,j)
        end do
    end do
    close(10)

end subroutine

subroutine ReadPar(FilePar, FileCoefficient, degree, N, FileOut)
    implicit none
    character(len=200)          ::FilePar
    character(len=200)          ::FileOut
    character(len=200)          ::FileCoefficient
    integer(kind=4)             ::degree
    integer(kind=4)             ::N

    character(len=200)          ::str

    open(10, file=FilePar, status='old')
    read(10,*) str, FileCoefficient
    read(10,*) str, degree
    read(10,*) str, N
    read(10,*) str, FileOut
    close(10)

end subroutine

subroutine Process(C, S, degree, position, n1, n2, V)
    implicit none
    real(kind = 8)        , intent(in   ) :: position(3)               ! positions in the earth-fixed system
    real(kind = 8)       , intent(  out) :: V                         ! gravitional v
    integer(kind = 4)     , intent(in   ) :: degree                    ! the max degree fo the stocks coefficients
    real(kind = 8)        , intent(in   ) :: C(0: degree, 0: degree), S(0: degree, 0: degree)
    integer(kind = 4)     , intent(in   ) :: n1, n2                    ! the degree span

    real(kind = 8), allocatable           :: E(:, :), F(:, :)
    real(kind = 8)                        :: x, y, z

    allocate(e(0: degree, 0: degree), f(0: degree, 0: degree))

    x = position(1)
    y = position(2)
    z = position(3)

    Call CalculateEF1(x, y, z, degree, E, F)
    Call CalculateV(E, F, degree, n1, n2, V, c, s)
    deallocate(e)
    deallocate(f)

end subroutine Process

subroutine CalculateEF1(x, y, z, n, E, F)
    implicit none
    integer(kind = 4), intent(in   ) :: n
    real(kind = 8)   , intent(in   ) :: x, y, z
    real(kind = 8)   , intent(  out) :: E(0: n, 0: n), F(0: n, 0: n)
    integer(kind = 4)                :: i, j
    real(kind = 8)                   :: r, REarth
    real(kind = 8)                   :: coe1, coe2


    REarth = 6.3781363E+6

    r = sqrt(x**2 + y**2 + z**2)

    E(0, 0) = REarth / r
    E(1, 0) = sqrt(3.0) * z * (REarth**2) / (r**3)
    E(1, 1) = sqrt(3.0) * x * (REarth**2) / (r**3)
    F(0, 0) = 0.0
    F(1, 0) = 0.0
    F(1, 1) = sqrt(3.0) * y * (REarth**2) / (r**3)

    do i = 2, n
        do j = 0, i
            if(i == j) then
                coe1 = sqrt(real(2 * j + 1) / (2 * j))
                E(i, j) = coe1 * (x * REarth * E(i - 1, j - 1) / (r**2) - y * REarth * F(i - 1, j - 1) / (r**2))
                F(i, j) = coe1 * (x * REarth * F(i - 1, j - 1) / (r**2) + y * REarth * E(i - 1, j - 1) / (r**2))
            else
                coe1 = sqrt(real(2 * i + 1) * (2 * i - 1) / ((i - j) * (i + j)))
                coe2 = sqrt((real(2 * i + 1) * (i - j - 1) * (i + j - 1)) / ((i - j) * (i + j) * (2 * i - 3)))
                E(i, j) = coe1 * z * REarth * E(i - 1, j) / (r**2) - coe2 * (REarth**2) * E(i - 2, j) / (r**2)
                F(i, j) = coe1 * z * REarth * F(i - 1, j) / (r**2) - coe2 * (REarth**2) * F(i - 2, j) / (r**2)
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
    real(kind = 8)   , intent(in   )             :: s(0: Degree, 0: Degree)
    real(kind = 8)   , intent(in   )             :: c(0: Degree, 0: Degree)

    integer(kind = 4)                            :: i, j
    real(kind = 8)                               :: REarth, MU

    REarth = 6.3781363E+6
    MU = 3.986004415000000e+14
    V = 0.0

    do i = n1, n2
        do j = 0, i
            V = V + E(i, j) * C(i, j) + F(i, j) * S(i, j)
        end do
    end do

    V = V * MU / REarth
end subroutine


subroutine CalculateT(x, y, z, V, n1, n2, C, S, degree)
    implicit none
    integer(kind=4)         ::degree
    integer(kind=4)         ::n1, n2
    real(kind=8)            ::V
    real(kind=8)            ::x, y, z
    real(kind=8)            ::C(0:degree,0:degree), S(0:degree,0:degree)


    real(kind=8)            ::REarth=6.3781363E+6, MU=0.3986004415E+15
    real(kind=8),allocatable::E(:,:), F(:,:), coe1(:,:)
    integer(kind=4)         ::i, j, coe2

    coe2 = n2
    if(degree>=2 .and. n2 ==1 ) then
        coe2 = 2
    end if

    allocate(E(0:coe2,0:coe2), F(0:coe2,0:coe2), coe1(1:2,1:3))

    call CalculateEF(x, y, z, REarth, coe2, E, F)

    if(degree>=2) then
        do i=1,3
            coe1(1,i) = C(2,i-1)
        end do
        do i=1,3
            coe1(2,i) = S(2,i-1)
        end do
    
        if(n2==1) then
            ! C(2,0) = -0.48000000000000e-03
            C(2,1) = 0
            C(2,2) = 0
            S(2,0) = 0
            S(2,1) = 0
            S(2,2) = 0
        else
            C(2,0) = 0!-0.004165143790815e-03
        end if
        
    end if


    V=0
    do i=n1,coe2
        do j=0,i
            V = V+E(i,j)*C(i,j)+F(i,j)*S(i,j)
        end do
    end do

    V = V*MU/REarth

    if(degree>=2) then
        do i=1,3
            C(2,i-1)=coe1(1,i)
        end do
        do i=1,3
            S(2,i-1)=coe1(2,i)
        end do
    end if

    deallocate(E,F,coe1)

end subroutine

subroutine CalculateEF(x, y, z, REarth, n, E, F)
    implicit none
    integer(kind=4) ::n
    real(kind=8)    ::x, y, z, REarth
    real(kind=8)    ::E(0:n,0:n), F(0:n,0:n)

    
    integer(kind=4) ::i, j
    real(kind=8)    ::r
    real(kind=8)    ::coe1, coe2

    r=sqrt(x**2+y**2+z**2)
    E(0,0)=REarth/r
    E(1,0)=sqrt(3.0)*z*(REarth**2)/(r**3)
    E(1,1)=sqrt(3.0)*x*(REarth**2)/(r**3)
    F(0,0)=0
    F(1,0)=0
    F(1,1)=sqrt(3.0)*y*(REarth**2)/(r**3)

    do i=2,n
        do j=0,i
            if(i==j) then
                coe1=sqrt(real(2*j+1)/(2*j))
                E(i,j)=coe1*(x*REarth*E(i-1,j-1)/(r**2)-y*REarth*F(i-1,j-1)/(r**2))
                F(i,j)=coe1*(x*REarth*F(i-1,j-1)/(r**2)+y*REarth*E(i-1,j-1)/(r**2))
            else
                coe1=sqrt(real(2*i+1)*(2*i-1)/((i-j)*(i+j)))
                coe2=sqrt((real(2*i+1)*(i-j-1)*(i+j-1))/((i-j)*(i+j)*(2*i-3)))
                E(i,j)=coe1*z*REarth*E(i-1,j)/(r**2)-coe2*(REarth**2)*E(i-2,j)/(r**2)
                F(i,j)=coe1*z*REarth*F(i-1,j)/(r**2)-coe2*(REarth**2)*F(i-2,j)/(r**2)
            end if
        end do
    end do

end subroutine

subroutine Calculate(flag, x, y, z, x_2, y_2, z_2, VelocityX, VelocityY, VelocityZ, C, S, degree, N, T1, T2, ConstE, P , EK, P2)
    implicit none 
    integer(kind=4) ::degree
    integer(kind=4) ::N
    integer(kind=4) ::flag
    real(kind=8)    ::ConstE
    real(kind=8)    ::T1    !观测方程的计算值
    real(kind=8)    ::T2    !EGM2008计算得到的理论值
    real(kind=8)    ::x, y, z
    real(kind=8)    ::x_2, y_2, z_2
    real(kind=8)    ::VelocityX, VelocityY, VelocityZ
    real(kind=8)    ::C(0:degree,0:degree), S(0:degree,0:degree)
    real(kind=8)    ::P, P2
    
    
    real(kind=8)    ::AngularVEarth=7.2921158553E-5     !地球平均自转角速度
    real(kind=8)    ::EK    !动能
    real(kind=8)    ::VW    !离心力位
    real(kind=8)    ::U0    !正常场
    real(kind=8)    ::U1

    real(kind=8)    ::position(3)


    EK=(VelocityX**2+VelocityY**2+VelocityZ**2)/2.0
    VW = AngularVEarth * (x*VelocityY - y*VelocityX)

    if(flag==0) then

        call CalculateT(x_2,y_2,z_2,T2,2,N,C,S,degree)
        call CalculateT(x_2,y_2,z_2,U0,0,1,C,S,degree)
        call CalculateT(x_2,y_2,z_2,U1,N+1,degree,C,S,degree)
        ConstE = EK - U0 - U1 - T2 - VW

    end if
    
    call CalculateT(x_2,y_2,z_2,U0,0,1,C,S,degree)
    call CalculateT(x_2,y_2,z_2,U1,N+1,degree,C,S,degree)
    T1 = EK - U0 - U1 - ConstE - VW
    call CalculateT(x_2,y_2,z_2,T2,2,N,C,S,degree)
    position(1) = x_2
    position(2) = y_2
    position(3) = z_2
    Call Process(C, S, degree, position, 0, N, P2)
    P = U0 + U1 + T2 

end subroutine
