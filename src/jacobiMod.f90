!-----------------------------------------------------------------------
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
!
!   This file is part of NEXD 1D.
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but
!   WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with NEXD 1D. If not, see <http://www.gnu.org/licenses/>.
!---------------------------------
module jacobiMod
    ! Module to deal with jacobi polinomials

    use constantsMod
    use rosettaGammaMod
    implicit none

    contains

    subroutine jacobiP(P, x, alpha, beta, N)
        ! Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
        !          (alpha+beta /= -1) at points x for order N and returns
        !          P[1:length(xp)]
        ! Note:    They are normalized to be orthonormal.

        !input variables
        real(kind=custom_real), dimension(:), intent(in) :: x
        real(kind=custom_real), intent(in) :: alpha, beta
        integer, intent(in) :: N
        !output
        real(kind=custom_real), dimension(:), intent(out) :: P
        !local variables
        real(kind=custom_real), dimension(N+1,size(x)) :: PL
        real(kind=custom_real) :: gamma0, gamma1, aold, anew, h1, bnew
        integer :: i, dim

        dim=size(x)
        !initial values for P_0(x) and P_1(x)
        gamma0 = 2**(alpha+beta+1.0)/(alpha+beta+1.0)*lacz_gamma(alpha+1.0)*lacz_gamma(beta+1.0)/lacz_gamma(alpha+beta+1.0)
        !gamma0 = 2**(alpha + beta + 1.0)/(alpha + beta + 1.0) * gamma(alpha + 1.0)*gamma(beta + 1.0)/gamma(alpha + beta + 1.0)
        PL(1,:) = 1.0/sqrt(gamma0)

        if (N == 0) then
            P = PL(1,:)
            return
        endif

        gamma1 = (alpha + 1.0) * (beta + 1.0)/(alpha + beta + 3.0) * gamma0
        PL(2,:) = ((alpha + beta + 2.0) * x(:)/2.0 + (alpha - beta)/2.0)/sqrt(gamma1)

        if (N == 1) then
            P = PL(N+1,:)
            return
        endif

        !repeat value in recurrence
        aold = 2.0/(2.0 + alpha + beta) * sqrt((alpha + 1.0) * (beta + 1.0)/(alpha + beta + 3.0))

        !forward recurrence using the symetry of the recurrence
        do i = 1, N-1
            h1 = 2.0 * i + alpha + beta
            anew = 2.0/(h1+2.0)*sqrt((i+1.0)*(i+1.0+alpha+beta)*(i+1.0+alpha)*(i+1.0+beta)/(h1+1.0)/(h1+3.0))
            bnew = -(alpha**2 - beta**2)/h1/(h1 + 2)
            PL(i + 2,:) = 1.0/anew * (-aold * PL(i,:) + (x(:) - bnew) * PL(i+1,:))
            aold = anew
        enddo

        P = PL(N+1,:)
    end subroutine jacobiP

    subroutine gradJacobiP(dP, x, alpha, beta, N)
        ! Purpose: Evaluate the derivative of the Jacobi polynomial of type
        !          (alpha,beta)>-1, at points r for order N and returns
        !          dP[1:length(r)]

        !input
        integer, intent(in) :: N
        real(kind=custom_real), dimension(:), intent(in) :: x !grid points
        real(kind=custom_real), intent(in) :: alpha, beta
        !output
        real(kind=custom_real), dimension(:), intent(out) :: dP
        !local variables
        real(kind=custom_real), dimension(size(dP)) :: h1

        if (N == 0) then
            dP(:) = 0
        else
            call jacobiP(h1, x, alpha + 1, beta + 1, N - 1)
            dP(:) = sqrt(N * (N + alpha + beta + 1.0)) * h1(:)
        endif
    end subroutine gradJacobiP

end module jacobiMod
