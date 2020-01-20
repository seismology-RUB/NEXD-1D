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
module matrixMod

    use constantsMod

    implicit none

    contains
    subroutine print_int_matrix(desc, A)
        character*(*) :: desc
        integer, dimension(:,:) ::  A
        integer :: i,j

        write(*,*)
        write(*,*) DESC
        do i = 1, size(A(:,1))
            write(*,9998) (A(i,j), j = 1, size(A(1,:)))
        enddo
        9998 format( 40(:,1X,i3) )
        return
    end subroutine print_int_matrix

    subroutine print_real_matrix(desc, A)
        character*(*) :: desc
        real(kind=custom_real), dimension(:,:) ::  A
        integer :: i,j

        write(*,*)
        write(*,*) DESC
        do i = 1, size(A(:,1))
            write(*,9998) (A(i,j), j = 1, size(A(1,:)))
        enddo
        9998 format( 11(:,1X,F11.7) )
        return
    end subroutine print_real_matrix

    subroutine print_log_matrix(desc, A)
        character*(*) :: desc
        logical, dimension(:,:) ::  A
        integer :: i,j

        write(*,*)
        write(*,*) DESC
        do i = 1, size(A(:,1))
            write(*,9998) (A(i,j), j = 1, size(A(1,:)))
        enddo
        9998 format( 40(:,1X,l3) )
        return
    end subroutine print_log_matrix

    subroutine print_real_vector(desc, vec)
        character*(*) :: desc
        real(kind=4), dimension(:) ::  vec
        integer :: i

        write(*,*)
        write(*,*) DESC
        do i = 1, size(vec)
            write(*,9998) (vec(i))
        enddo
        9998 format( 11(:,1X,F11.7) )
        return
    end subroutine print_real_vector

    subroutine print_int_vector(desc, vec)
        character*(*) :: desc
        integer, dimension(:) ::  vec
        integer :: i

        write(*,*)
        write(*,*) DESC
        do i = 1, size(vec)
            write(*,9998) (vec(i))
        enddo
        9998 format( 11(:,1X,i3) )
        return
    end subroutine print_int_vector

    subroutine invert(V)
        !invert matrix
        real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: V
        real(kind=CUSTOM_REAL), dimension(2*size(v(:,1))) :: work
        integer, dimension(size(v(:,1))) :: ipvt
        integer :: ierr, iwork
        integer :: N
        N=size(v(:,1))
        iwork=2*N

        !LU trafo
        if (custom_real == size_real) then
            call sgetrf(N,N,V,N,ipvt,ierr)
        elseif (custom_real == size_double) then
            call dgetrf(N,N,V,N,ipvt,ierr)
        else
            if (ierr/=0) write(*,*) "Error LU vdm2d ",ierr
        endif
        ! ivert Pr
        if (custom_real == size_real) then
            call sgetri(N,V,N,ipvt,work,iwork,ierr)
        elseif (custom_real == size_double) then
            call dgetri(N,V,N,ipvt,work,iwork,ierr)
        else
            if (ierr/=0) write(*,*) "Error invert vdm2d ",ierr
        endif


    end subroutine invert

end module
