!-----------------------------------------------------------------------
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
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
module riemannflux1Dmod

    use constantsMod
    use materials1DMod

    implicit none

    contains

    subroutine rightRiemannflux(dv, dsigma, zin, zout, vel, fluxR)
        !Input
        real(kind=custom_real), intent(in) :: dv
        real(kind=custom_real), intent(in) :: dsigma
        real(kind=custom_real), intent(in) :: zin  !zk
        real(kind=custom_real), intent(in) :: zout !zk-1
        real(kind=custom_real), intent(in) :: vel
        !output
        real(kind=custom_real), dimension(:), intent(out) :: fluxR
        !local variables
        real(kind=custom_real) :: tmp

        tmp = vel/(zin + zout) * (-zout*dv+dsigma)

        fluxR(1) = tmp              !flux v
        fluxR(2) = -zin * tmp       !flux sigma
    end subroutine rightRiemannflux

    subroutine leftRiemannflux(dv, dsigma, zin, zout, vel, fluxL)
        !Input
        real(kind=custom_real), intent(in) :: dv
        real(kind=custom_real), intent(in) :: dsigma
        real(kind=custom_real), intent(in) :: zin !zk
        real(kind=custom_real), intent(in) :: zout !zk+1
        real(kind=custom_real), intent(in) :: vel
        !output
        real(kind=custom_real), dimension(:), intent(out) :: fluxL
        !local variables
        real(kind=custom_real) :: tmp

        tmp = -vel/(zin + zout) * (zout*dv+dsigma)
        fluxL(1) = tmp              !flux v
        fluxL(2) = zin * tmp        !flux sigma
    end subroutine leftRiemannflux

    subroutine GodunovFlux(UL, Ain, nx, flux)
        !input
        real(kind=custom_real), dimension(:), intent(in) :: UL
        real(kind=custom_real), dimension(:,:), intent(in) :: Ain
        integer, intent(in) :: nx
        !output
        real(kind=custom_real), dimension(:), intent(out) :: flux
        !local variables
        real(kind=custom_real), dimension(:,:), allocatable :: A
        integer, dimension(:,:), allocatable :: T
        integer :: N, i
        integer, dimension(2) :: shapeA

        shapeA=shape(Ain)
        N=shapeA(1)
        allocate(A(N,N))
        allocate(T(N,N))
        A=Ain

        T = 0
        do i=1,(N/2)
            T(2*i-1,2*i-1) = nx**2
            T(2*i  ,2*i  ) = nx
        enddo

        A = matmul(T,matmul(A,T))

        flux = matmul(A,UL)

        if (allocated(A)) deallocate(A)
        if (allocated(T)) deallocate(T)
    end subroutine GodunovFlux

    subroutine rightRiemannFluxLSI(S,vel,zin,zout,fluxR)
        !Input
        real(kind=custom_real), intent(in) :: vel
        real(kind=custom_real), intent(in) :: zin !zk
        real(kind=custom_real), intent(in) :: zout !zk-1
        real(kind=custom_real), intent(in) :: S
        !output
        real(kind=custom_real), dimension(:), intent(out) :: fluxR
        !local variables
        real(kind=custom_real) :: tmp

        tmp = vel/(zin + zout)
        fluxR(1) = tmp * zout * 2. * S
        fluxR(2) = -zin * tmp * zout * 2. * S
    end subroutine rightRiemannFluxLSI

    subroutine leftRiemannFluxLSI(S, vel, zin, zout, fluxL)
        !Input
        real(kind=custom_real), intent(in) :: vel
        real(kind=custom_real), intent(in) :: zin !zk
        real(kind=custom_real), intent(in) :: zout !zk+1
        real(kind=custom_real), intent(in) :: S
        !output
        real(kind=custom_real), dimension(:), intent(out) :: fluxL
        !local variables
        real(kind=custom_real) :: tmp

        tmp = vel/(zin + zout)
        fluxL(1) = tmp * zout * 2. * S
        fluxL(2) = zin * tmp * zout * 2. * S
    end subroutine leftRiemannFluxLSI
end module riemannflux1Dmod
