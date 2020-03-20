subroutine Process(C, S, degree, position, n1, n2, V)
    implicit none
    real(kind = 8)             , intent(in   ) :: position(3)               ! positions in the earth-fixed system
    real(kind = 8)             , intent(  out) :: V                         ! gravitional v
    integer(kind = 4)          , intent(in   ) :: degree                    ! the max degree fo the stocks coefficients
    real(kind = 8)             , intent(in   ) :: C(0: degree, 0: degree), S(0: degree, 0: degree)
    integer(kind = 4)          , intent(in   ) :: n1, n2                    ! the degree span

    real(kind = 8)                             :: E(0: degree, 0: degree), F(0: degree, 0: degree)
    real(kind = 8)                             :: x, y, z

    x = position(1)
    y = position(2)
    z = position(3)

    Call CalculateEF(x, y, z, degree, E, F)
    Call CalculateV(E, F, degree, n1, n2, V, c, s)

end subroutine Process

subroutine CalculateEF(x, y, z, n, E, F)
    implicit none
    integer(kind = 4), intent(in   ) :: n
    real(kind = 8)   , intent(in   ) :: x, y, z
    real(kind = 8)   , intent(  out) :: E(0: n, 0: n), F(0: n, 0: n)


    integer(kind = 4) :: i, j
    real(kind = 8)    :: r, REarth
    real(kind = 8)    :: coe1, coe2

    REarth = 6.3781363E+6

    r = sqrt(x**2 + y**2 + z**2)

    E(0, 0) = REarth / r
    E(1, 0) = sqrt(3.0) * z * (REarth**2) / (r**3)
    E(1, 1) = sqrt(3.0) * x * (REarth**2) / (r**3)
    F(0, 0) = 0
    F(1, 0) = 0
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

    integer(kind = 4)             :: i, j
    real(kind = 8)                :: REarth, MU

    REarth = 6.3781363E+6
    MU = 3.986004415000000e+14
    V = 0
    do i = n1, n2
        do j = 0, i
            V = V + E(i, j) * C(i, j) + F(i, j) * S(i, j)
        end do
    end do

    V = V * MU / REarth
end subroutine