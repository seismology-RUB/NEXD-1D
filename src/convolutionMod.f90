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
module convolutionMod
    use parameterMod
    use constantsMod

    implicit none

    type convolution
        real(kind=custom_real)  :: dt
        real(kind=custom_real), dimension(:), pointer :: data => null()
    end type

    contains

    function convolve(h,d) result(res)
        !input
        type(convolution), intent(in) :: d
        real(kind=custom_real), dimension(:), intent(in) :: h
        !output
        !type(convolution) :: res
        real(kind=custom_real), dimension(size(d%data)) :: res
        !local
        integer :: i, j
        integer :: nh, nd

        nh = size(h)
        nd = size(d%data)

        if (nd >= nh) then
            do i=1, max(nh, nd)
                res(i) = 0.
                do j=1,min(i,nh)
                    res(i) = res(i) + h(j)*d%data(i-j+1) !exp(-nu*dt*(i-j))
                enddo
            enddo
        else
            do i=1,max(nh,nd)
                res(i) = 0.
                do j=1,min(i,nd)
                    res(i) = res(i) + d%data(j)*h(i-j+1)
                enddo
            enddo
        endif
        res = res/(nd)
    end function convolve

end module
