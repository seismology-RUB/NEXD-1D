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
module receiverMod
    use constantsMod
    use parameterMod
    use errorMessage
    use vandermonde1DMod
    use localOperatorsMod
    use matrixMod
    implicit none

    type :: receiverVar
        integer :: element
        real(kind=custom_real) :: pos
        real(kind=custom_real), dimension(1) :: r
        real(kind=custom_real), dimension(:), allocatable :: data
        real(kind=custom_real), dimension(:), allocatable :: dataSi
        real(kind=custom_real), dimension(:), allocatable :: dataV1
        real(kind=custom_real), dimension(:), allocatable :: dataP1
        real(kind=custom_real), dimension(:), allocatable :: dataV2
        real(kind=custom_real), dimension(:), allocatable :: dataP2
        real(kind=custom_real), dimension(:), allocatable :: interpolation
    end type

    contains

    subroutine initReceiver(par, mesh, rec, vdm, lu, filename, errmsg)
        !input
        type(parameterVar), intent(in) :: par
        type (error_message) :: errmsg
        type(meshVar) :: mesh
        character(len=*), intent(in) :: filename
        integer :: lu
        real(kind=custom_real), dimension(:,:) :: vdm
        !input/output
        type(receiverVar), dimension(:), allocatable :: rec
        !local
        integer :: i
        integer :: ios
        integer :: out
        character(len=22) :: myname = 'readReceiverParameters'
        character(len=255) :: text
        character(len=80) :: outfile


        call addTrace(errmsg,myname)
        open(lu,file = filename, status = 'old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            return
        endif

        text = 'none'
        do while(text /= 'BEGIN')
            read(lu,'(a)') text
        enddo

        allocate(rec(par%recn))

        do i = 1, par%recn
            read(lu,*) rec(i)%pos
            allocate(rec(i)%data(par%tsteps))
            allocate(rec(i)%dataSi(par%tsteps))
            if (par%poroelastic) then
                allocate(rec(i)%dataV1(par%tsteps))
                allocate(rec(i)%dataP1(par%tsteps))
                if (par%fluidn == 2) then
                    allocate(rec(i)%dataV2(par%tsteps))
                    allocate(rec(i)%dataP2(par%tsteps))
                endif
            endif
            if (rec(i)%pos > mesh%xmax .or. rec(i)%pos < mesh%xmin) then
                call add(errmsg,2,'Receiver not on the grid' ,myname)
                return
            endif
            call findReceiver(mesh, rec(i))
            call interpolateReceiver(mesh, vdm, rec(i))
        enddo
        close(lu)

        if (par%log) then
            write (*,"(a80)") "|                              Data on Receiver:                               |"
            write (*,"(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a14, a13, a53)") "|   Number   |","  Position  |", "                                             |"
            write (*,"(a80)") "|------------|------------|----------------------------------------------------|"
            do i = 1, par%recn
               write (*,"(a5, i4, a5, a1, f10.2, a2, a53)") &
               "|    ", i, "    |", " ", rec(i)%pos, " |", "                                             |"
            enddo
            write (*,"(a80)") "|------------------------------------------------------------------------------|"
        endif
        if (par%logfile) then
            outfile = "output/logfile"//trim(par%title)
            outfile = trim(outfile)
            out = 16
            open(out, file = outfile, position="append")
            write (out,"(a80)") "|                              Data on Receiver:                               |"
            write (out,"(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a14, a13, a53)") "|   Number   |","  Position  |", "                                             |"
            write (out,"(a80)") "|------------|------------|----------------------------------------------------|"
            do i = 1, par%recn
               write (out,"(a5, i4, a5, a1, f10.2, a2, a53)") &
               "|    ", i, "    |", " ", rec(i)%pos, " |", "                                             |"
            enddo
            write (out,"(a80)") "|------------------------------------------------------------------------------|"
            close(out)
        endif
    end subroutine

    subroutine findReceiver(mesh, rec)
        ! Find the local coordinate of the receiver (position in the element)
        !input
        type(meshVar), intent(in)   :: mesh
        !in/out
        type(receiverVar) :: rec
        !local
        real(kind=custom_real) :: elementCenter
        real(kind=custom_real) :: diff

        if (rec%pos < epsilon(rec%pos)) then
            rec%pos = rec%pos + eps
        endif
        rec%element = ceiling(((rec%pos+abs(mesh%xmin))*mesh%ncell/(mesh%xmax+abs(mesh%xmin)))+mesh%nghost)
        !Calculate the center of the element the receiver is placed in
        elementCenter = mesh%x(1, rec%element) + (mesh%x(mesh%np, rec%element) - mesh%x(1, rec%element))/2
        !Calculate the position relative to the center of the element
        diff = rec%pos - elementCenter
        !scale to the size of the element
        rec%r(1) = diff/((mesh%x(mesh%np, rec%element) - mesh%x(1, rec%element))/2)
    end subroutine

    subroutine interpolateReceiver(mesh, vdm, rec)
        !input
        type(meshVar), intent(in)   :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        !in/out
        type(receiverVar) :: rec
        !local
        real(kind=custom_real), dimension(mesh%np,mesh%np) :: vdmTinv  !Transposed and inverted vandermonde matrix (V^T)^-1
        real(kind=custom_real), dimension(mesh%np) :: tmp

        allocate(rec%interpolation(mesh%np))
        vdmTinv = transpose(vdm)
        call invert(vdmTinv)
        call vdm1d(mesh, rec%r, tmp)
        rec%interpolation = matmul(vdmTinv, tmp)
    end subroutine

    subroutine recordSeismograms(par, mesh, rec, tstep, solution)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar) :: mesh
        type(receiverVar), dimension(:) :: rec
        integer, intent(in) :: tstep
        type (solutionvector), intent(in) :: solution
        !real(kind=custom_real), dimension(:,:), intent(in) :: V
        !local
        integer :: i, j, r

        do i = 1, mesh%k
            do r = 1, par%recn
                if (i == rec(r)%element) then
                    do j = 1, mesh%np
                        rec(r)%data(tstep) = rec(r)%data(tstep) + solution%v(j,i)*rec(r)%interpolation(j)
                        rec(r)%dataSi(tstep) = rec(r)%dataSi(tstep) + solution%sigma(j,i)*rec(r)%interpolation(j)
                        if (par%poroelastic) then
                            rec(r)%dataV1(tstep) = rec(r)%dataV1(tstep) + solution%v1(j,i)*rec(r)%interpolation(j)
                            rec(r)%dataP1(tstep) = rec(r)%dataP1(tstep) + solution%p1(j,i)*rec(r)%interpolation(j)
                            if (par%fluidn == 2) then
                                rec(r)%dataV2(tstep) = rec(r)%dataV2(tstep) + solution%v2(j,i)*rec(r)%interpolation(j)
                                rec(r)%dataP2(tstep) = rec(r)%dataP2(tstep) + solution%p2(j,i)*rec(r)%interpolation(j)
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo
    end subroutine

end module
