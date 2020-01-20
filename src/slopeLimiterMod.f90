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
module slopeLimiterMod
    use parameterMod
    use constantsMod
    use errorMessage
    use localOperatorsMod
    implicit none

    type :: limiterVar
        real(kind=custom_real), dimension(:), allocatable :: avg
        real(kind=custom_real), dimension(:), allocatable :: avgr
        real(kind=custom_real), dimension(:), allocatable :: avgl
        real(kind=custom_real), dimension(:,:), allocatable :: qmodal
        real(kind=custom_real), dimension(:,:), allocatable :: qavg
    end type

    contains

    subroutine selectLimiter(par, mesh, vdm, invV, v, s, errmsg)
        !input
        type(meshVar), intent(in) :: mesh
        type(parameterVar), intent(in) :: par
        type (error_message) :: errmsg
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        real(kind=custom_real), dimension(:,:), intent(in) :: invV
        !in/output
        real(kind=custom_real), dimension(:,:), intent(inout) :: v
        real(kind=custom_real), dimension(:,:), intent(inout) :: s
        !local
        character(len=12) :: myname = "slopeLimiter"
        call addTrace(errmsg,myname)

        select case(trim(par%limiter))
            case ("minmod")
                call minmodLimiter(mesh, vdm, invV, v)
                call minmodLimiter(mesh, vdm, invV, s)
            case ("superbee")
                call superbeeLimiter(mesh, vdm, invV, v)
                call superbeeLimiter(mesh, vdm, invV, s)
            case ("MC")
                call mcLimiter(mesh, vdm, invV, v)
                call mcLimiter(mesh, vdm, invV, s)
            case ("none")
                return
            case default
                call add(errmsg,2,'No valid limiter has been chosen',myname)
        end select
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    end subroutine selectLimiter

    subroutine limiterBasis(mesh, lim, vdm, invV, q, start, end)
        !input
        type(meshVar), intent(in) :: mesh
        integer, intent(in) :: start, end
        real(kind=custom_real), dimension(:,:), intent(in) :: q !Solutionmatrix to be limited
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        real(kind=custom_real), dimension(:,:), intent(in) :: invV
        !in-/output
        type(limiterVar) :: lim

        allocate(lim%avg(size(q(1,start:end))))
        allocate(lim%avgr(size(q(1,start:end))))
        allocate(lim%avgl(size(q(1,start:end))))

        !compute modal coefficients
        lim%qmodal = matmul(invV, q)
        !extract cell averages
        lim%qmodal(2:mesh%Np,:) = 0.
        lim%qavg = matmul(vdm,lim%qmodal)
        lim%avg(:) = lim%qavg(1,start:end)
        lim%avgl(:) = (/lim%qavg(1,start),lim%qavg(1,start:end-1)/)
        lim%avgr(:) = (/lim%qavg(1,start+1:end),lim%qavg(1,end)/)
    end subroutine limiterBasis

    subroutine minmodLimiter(mesh, vdm, invV, q)
        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        real(kind=custom_real), dimension(:,:), intent(in) :: invV
        !in-/output
        real(kind=custom_real), dimension(:,:), intent(inout) :: q !Solutionmatrix to be limited
        !local
        type(limiterVar) :: lim
        integer :: i, j
        real(kind=custom_real) :: r,l
        real(kind=custom_real) :: slope
        real(kind=custom_real), dimension(mesh%k) :: x0

        !cell center points
        x0(:) = mesh%x(1,:) + mesh%h(:)/2.

        allocate(lim%qmodal(mesh%np,mesh%k))
        allocate(lim%qavg(mesh%np,mesh%k))
        call limiterBasis(mesh, lim, vdm, invV, q, 1, mesh%k)
        do i = 1, mesh%k
            do j = 1, mesh%np
                r = (lim%avgr(i) - lim%avg(i))/mesh%h(i)
                l = (lim%avg(i) - lim%avgl(i))/mesh%h(i)
                slope = minmod((/r,l/))
                q(j,i) = lim%avg(i) + slope*(mesh%x(j,i) - x0(i))
            enddo
        enddo
    end subroutine minmodLimiter

    subroutine superbeeLimiter(mesh, vdm, invV, q)
        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        real(kind=custom_real), dimension(:,:), intent(in) :: invV
        !in-/output
        real(kind=custom_real), dimension(:,:), intent(inout) :: q !Solutionmatrix to be limited
        !local
        type(limiterVar) :: lim
        integer :: i, j
        real(kind=custom_real) :: r,l
        real(kind=custom_real) :: slope, s1, s2
        real(kind=custom_real), dimension(mesh%k) :: x0

        !cell center points
        x0(:) = mesh%x(1,:) + mesh%h(:)/2.

        allocate(lim%qmodal(mesh%np,mesh%k))
        allocate(lim%qavg(mesh%np,mesh%k))
        call limiterBasis(mesh, lim, vdm, invV, q, 1, mesh%k)

        do i = 1, mesh%k
            do j = 1, mesh%np
                r = (lim%avgr(i) - lim%avg(i))/mesh%h(i)
                l = (lim%avg(i) - lim%avgl(i))/mesh%h(i)
                s1 = minmod((/r, 2*l/))
                s2 = minmod((/2*r, l/))
                slope = maxmod(s1,s2)
                q(j,i) = lim%avg(i) + slope*(mesh%x(j,i) - x0(i))
            enddo
        enddo
    end subroutine superbeeLimiter

    subroutine mcLimiter(mesh, vdm, invV, q)
        !input
        type(meshVar), intent(in) :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        real(kind=custom_real), dimension(:,:), intent(in) :: invV
        !in-/output
        real(kind=custom_real), dimension(:,:), intent(inout) :: q !Solutionmatrix to be limited
        !local
        type(limiterVar) :: lim
        integer :: i, j
        real(kind=custom_real) :: r,l,c
        real(kind=custom_real) :: slope
        real(kind=custom_real), dimension(mesh%k) :: x0

        !cell center points
        x0(:) = mesh%x(1,:) + mesh%h(:)/2.
        allocate(lim%qmodal(mesh%np,mesh%k))
        allocate(lim%qavg(mesh%np,mesh%k))
        call limiterBasis(mesh, lim, vdm, invV, q, 1, mesh%k)

        do i = 1, mesh%k
            do j = 1, mesh%np
                r = (lim%avgr(i) - lim%avg(i))/mesh%h(i)
                l = (lim%avg(i) - lim%avgl(i))/mesh%h(i)
                c = (lim%avgr(i) - lim%avgl(i))/(2*mesh%h(i))
                slope = minmod((/c, 2*l, 2*r/))
                q(j,i) = lim%avg(i) + slope*(mesh%x(j,i) - x0(i))
            enddo
        enddo
    end subroutine mcLimiter

    function minmod(in) result(res)
        !input
        real(kind=custom_real), dimension(:), intent(in) :: in
        !output
        real(kind=custom_real) :: res
        !local
        integer :: m, s
        logical :: ids

        ids = .false.
        res = 0

        m = size(in,1)
        s = sum(signum(in))/m
        !if (abs(s) == 1) ids = .true.
        if (s == -1) then
            res = maxval(in)
        else if (s == 1) then
            res = minval(in)
        else
            res = 0
        endif
    end function minmod

    function maxmod(a,b) result(res)
        !input
        real(kind=custom_real), intent(in) :: a, b
        !output
        real(kind=custom_real) :: res
        !local
        real(kind=custom_real) :: tmp

        res = 0
        tmp = a*b

        if (abs(a) > abs(b) .and. tmp > 0) then
            res = a
        else if(abs(b) > abs(a) .and. tmp > 0) then
            res = b
        else if(tmp <= 0) then
            res = 0
        endif
    endfunction maxmod

    !This function is used to calculate the sign of the entries of a vector
    function signum(v) result(res)
        !input
        real(kind=custom_real), dimension(:), intent(in) :: v
        !output
        integer, dimension(size(v)) :: res
        !local
        integer :: i

        do i = 1, size(v)
            if (v(i) < 0) res(i) = -1
            if (v(i) < epsilon(v(i))) res(i) = 0
            !if (v(i) == 0) res(i) = 0
            if (v(i) > 0) res(i) = 1
        enddo
    end function

    subroutine deallocLimiterArrays(lim)
        type(limiterVar) :: lim
        if (allocated(lim%avg)) deallocate(lim%avg)
        if (allocated(lim%avgr)) deallocate(lim%avgr)
        if (allocated(lim%avgl)) deallocate(lim%avgl)
        if (allocated(lim%qmodal)) deallocate(lim%qmodal)
        if (allocated(lim%qavg)) deallocate(lim%qavg)
    end subroutine
end module
