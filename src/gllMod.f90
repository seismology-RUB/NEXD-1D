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
module gllMod
    ! module to generate gll points
    use constantsMod
    use parameterMod

    implicit none

    contains

    function getGll(mesh)
        ! Function to get gll points for order NGLL-1
        type(meshVar), intent(in) :: mesh
        real(kind=CUSTOM_REAL), dimension(mesh%np) :: x,GetGll
        real(kind=CUSTOM_REAL), dimension(mesh%np,mesh%np):: V
        real(kind=CUSTOM_REAL) :: xold
        integer :: i,k,p

        p = mesh%np-1

        do i=1,mesh%np
            x(i) = cos(pi*(i-1)/p)
        end do

        ! Newton-Raphson iteration
        do i=1,mesh%np
            xold = 2.
            do while (abs(x(i)-xold) > eps)
                xold = x(i)
                V(i,1) = 1.
                V(i,2) = x(i)
                do k=2,p
                    V(i,k+1) = ((2.*k-1.)*x(i)*V(i,k)-(k-1.)*V(i,k-1))/k
                end do
                x(i) = xold-(x(i)*V(i,mesh%np)-V(i,p))/(mesh%np*V(i,mesh%np))
            end do
        end do
        GetGLL=x
    end function getGll
end module gllMod
