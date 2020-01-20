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
module sourcesMod

    use constantsMod
    use parameterMod
    use genericMaterial
    use errorMessage
    use matrixMod
    use vandermonde1DMod
    use localOperatorsMod

    implicit none

    type :: sourceVar
        character(len=10) :: type
        character(len=1):: direction
        integer :: element
        real(kind=custom_real), dimension(1) :: r                          !Local coordinate of the source in the element
        real(kind=custom_real):: center
        real(kind=custom_real) :: width
        real(kind=custom_real) :: f0
        real(kind=custom_real) :: t0
        real(kind=custom_real) :: amplitude
        real(kind=custom_real), dimension(:), allocatable :: stf
        real(kind=custom_real), dimension(:), allocatable :: diffStf
        real(kind=custom_real), dimension(:), allocatable :: interpolation
    end type

    contains

    subroutine readSourceParameters(par, mesh, src, lu, filename, errmsg)
        !input
        type(parameterVar), intent(in) :: par
        type (error_message) :: errmsg
        type(meshVar) :: mesh
        character(len=*) :: filename
        integer :: lu
        !input/output
        type(sourceVar), dimension(:), allocatable, intent(inout) :: src
        !local
        integer :: i
        integer :: ios
        integer :: out
        character(len=21) :: myname = 'readSourceParameters'
        character(len=255) :: text
        character(len=80) :: outfile
        real(kind=custom_real) :: tmp !temporary storage for the 3rd parameter in thr source file.

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

        allocate(src(par%srcn))

        do i = 1, par%srcn
            read(lu,*) src(i)%type, src(i)%center, tmp, src(i)%amplitude, src(i)%direction
            select case (trim(src(i)%type))
                case ("ricker", "green", "external", "sin3", "dg")
                    if (src(i)%center < mesh%xmin .or. src(i)%center > mesh%xmax) then
                        call add(errmsg,2,'Wavelet not on the grid' ,myname)
                    endif
                    src(i)%f0 = tmp
                    src(i)%t0 = 1.2/src(i)%f0
                case default
                    src(i)%width = tmp
                    if ((src(i)%center-0.5*src(i)%width) < mesh%xmin .or. (src(i)%center+0.5*src(i)%width) > mesh%xmax) then
                        call add(errmsg,2,'Wavelet not on the grid' ,myname)
                    endif
            end select
        enddo
        close(lu)

        if (par%log) then
            write (*,"(a80)") "|                               Data on sources:                               |"
            write (*,"(a80)") "|------------------------------------------------------------------------------|"

            do i = 1, par%srcn
                select case (trim(src(i)%type))
                    case ("ricker", "green", "external", "sin3", "dg")
                        write (*,"(a14, a13, a13, a13, a13, a14)") "|   Number   |","    Type    |", "  Position  |",&
                        " Peak Freq. |", "  Amplitude |", "  Direction  |"
                        write (*,"(a80)") "|------------|------------|------------|------------|------------|-------------|"
                        write (*,"(a5, i4, a5, a1, a10, a2, a1, f10.2, a2, a1, es10.3, a2, a1, es10.3, a2 ,a6, a2, a6)") &
                        "|    ", i, "    |", " ", trim(src(i)%type), " |"," ", src(i)%center," |",&
                        " ", src(i)%f0," |", " ", src(i)%amplitude, " |", "      ",trim(src(i)%direction), "     |"
                    case default
                        write (*,"(a14, a13, a13, a13, a13, a14)") "|   Number   |","    Type    |", "   Center   |",&
                        "    Width   |", "  Amplitude |", "  Direction  |"
                        write (*,"(a80)") "|------------|------------|------------|------------|------------|-------------|"
                        write (*,"(a5, i4, a5, a1, a10, a2, a1, f10.2, a2, a1, f10.2, a2, a1, f10.2, a2 ,a6, a2, a6)") &
                        "|    ", i, "    |", " ", trim(src(i)%type), " |"," ", src(i)%center," |",&
                        " ", src(i)%width," |", " ", src(i)%amplitude, " |", "      ",trim(src(i)%direction), "     |"
                end select
            enddo
            write (*,"(a80)") "|------------------------------------------------------------------------------|"
        endif
        if (par%logfile) then
            outfile = "output/logfile"//trim(par%title)
            outfile = trim(outfile)
            out = 16
            open(out, file = outfile, position="append")
            write (out,"(a80)") "|                               Data on sources:                               |"
            write (out,"(a80)") "|------------------------------------------------------------------------------|"

            do i = 1, par%srcn
                select case (trim(src(i)%type))
                    case ("ricker", "green", "external","sin3","dg")
                        write (out,"(a14, a13, a13, a13, a13, a14)") "|   Number   |","    Type    |", "  Position  |",&
                        " Peak Freq. |", "  Amplitude |", "  Direction  |"
                        write (out,"(a80)") "|------------|------------|------------|------------|------------|-------------|"
                        write (out,"(a5, i4, a5, a1, a10, a2, a1, f10.2, a2, a1, es10.3, a2, a1, es10.3, a2 ,a6, a2, a6)") &
                        "|    ", i, "    |", " ", trim(src(i)%type), " |"," ", src(i)%center," |",&
                        " ", src(i)%f0," |", " ", src(i)%amplitude, " |", "      ",trim(src(i)%direction), "     |"
                    case default
                        write (out,"(a14, a13, a13, a13, a13, a14)") "|   Number   |","    Type    |", "   Center   |",&
                        "    Width   |", "  Amplitude |", "  Direction  |"
                        write (out,"(a80)") "|------------|------------|------------|------------|------------|-------------|"
                        write (out,"(a5, i4, a5, a1, a10, a2, a1, f10.2, a2, a1, f10.2, a2, a1, f10.2, a2 ,a6, a2, a6)") &
                        "|    ", i, "    |", " ", trim(src(i)%type), " |"," ", src(i)%center," |",&
                        " ", src(i)%width," |", " ", src(i)%amplitude, " |", "      ",trim(src(i)%direction), "     |"
                end select
            enddo
            write (out,"(a80)") "|------------------------------------------------------------------------------|"
            close(out)
        endif

    end subroutine

    subroutine init_sources(par, mesh, genmat, vdm, src, solution, dt, errmsg)
        !input
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type (generic_material) :: genmat
        type (error_message) :: errmsg
        real(kind=custom_real) :: dt
        real(kind=custom_real), dimension(:,:) :: vdm
        !input/output
        type(sourceVar), dimension(:), allocatable, intent(inout) :: src
        type (solutionvector), intent(inout) :: solution
        !local
        character(len=30) :: filename, filenr
        character(len=21) :: myname = 'init_sources'
        integer :: i, it
        real(kind=custom_real) :: time

        call addTrace(errmsg,myname)
        call readSourceParameters(par, mesh, src, 1, "input/sources", errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        !Initial conditions
        solution%v = 0.
        solution%sigma = 0.
        if (allocated(solution%v1)) solution%v1 = 0.
        if (allocated(solution%p1)) solution%p1 = 0.
        if (allocated(solution%v2)) solution%v2 = 0.
        if (allocated(solution%p2)) solution%p2 = 0.

        !Initialize Source Wavelet
        do i = 1, par%srcn

            if (src(i)%direction == "R") then
                call addTrace(errmsg,myname)
                select case (trim(src(i)%type))
                    case ("sin2")
                        call initRightpropSineSquared(par, mesh, src(i), genmat, solution)
                    case ("box")
                        call initRightpropBoxcar(par, mesh, src(i), genmat, solution)
                    case ("null")
                        call initNull(par,mesh,solution)
                    case ("test")
                        call initRightpropTestwavelet(mesh, src(i), genmat, solution%v, solution%sigma)
                    case default
                        call add(errmsg,2,'No valid initial solution has been chosen',myname)
                end select
                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
            else if (src(i)%direction == "L") then
                call addTrace(errmsg,myname)
                select case (trim(src(i)%type))
                    case ("sin2")
                        call initLeftpropSineSquared(par, mesh, src(i), genmat, solution)
                    case ("box")
                        call initLeftpropBoxcar(par, mesh, src(i), genmat, solution)
                    case ("null")
                        call initNull(par,mesh,solution)
                    case ("test")
                        call add(errmsg,2,'No left propagting testwavelet defined',myname)
                    case default
                        call add(errmsg,2,'No valid initial solution has been chosen',myname)
                end select
                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
            else if (src(i)%direction == "B") then
                call findSource(mesh, src(i))
                call interpolateSource(mesh, vdm, src(i))
                select case (trim(src(i)%type))
                    case ("ricker", "green", "dg")
                        allocate(src(i)%stf(par%tsteps))
                        allocate(src(i)%diffStf(par%tsteps))
                        write(filenr, "(I3.3)") i
                        filename = "output/stf"//trim(filenr)
                        open(unit = 1, file=trim(filename), status="unknown")
                        do it = 1, par%tsteps
                            time = (float(it) -1.)*dt
                            if (trim(src(i)%type) == "ricker") then
                                src(i)%stf(it) = - stfRicker(time, src(i))
                                src(i)%diffStf(it) = - stfDiffRicker(time, src(i))
                            else if (trim(src(i)%type) == "dg") then
                                src(i)%stf(it) = - stfDiffGauss(time, src(i))
                            else
                                src(i)%stf(it) = - stfGauss(time, src(i))
                            endif
                            write (1,*) time, src(i)%stf(it)
                        enddo
                        close(1)
                    case ("external")
                        allocate(src(i)%stf(par%tsteps))
                        allocate(src(i)%diffStf(par%tsteps))
                        call stfExternal(src(i), dt, par%tsteps, .true.,6,3,trim(par%extwavelet),errmsg)
                        write(filenr, "(I3.3)") i
                        filename = "output/stf"//trim(filenr)
                        open(unit = 1, file=trim(filename), status="unknown")
                        do it = 1, par%tsteps
                            time = (float(it) -1.)*dt
                            write (1,*) time, src(i)%stf(it)
                        enddo
                        close(1)
                    case ("sin3")
                        allocate(src(i)%stf(par%tsteps))
                        allocate(src(i)%diffStf(par%tsteps))
                        write(filenr, "(I3.3)") i
                        filename = "output/stf"//trim(filenr)
                        open(unit = 1, file=trim(filename), status="unknown")
                        do it = 1, par%tsteps
                            time = (float(it) -1.)*dt
                            src(i)%stf(it) = stfSin3(time, src(i))
                            write (1,*) time, src(i)%stf(it)
                        enddo
                        close(1)
                    case default
                        call add(errmsg,2,'No valid initial solution has been chosen',myname)
                end select
            else
                call add(errmsg,2, 'No propagation direction specified. Exiting....', myname)
                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
            endif
        enddo
    end subroutine

    subroutine findSource(mesh, src)
        ! Find the local coordinate of the source (position in the element)
        !input
        type(meshVar), intent(in)   :: mesh
        !in/out
        type(sourceVar) :: src
        !local
        real(kind=custom_real) :: elementCenter
        real(kind=custom_real) :: diff

        if (src%center < epsilon(src%center)) then
            src%center = src%center + eps
        endif
        src%element = ceiling(((src%center+abs(mesh%xmin))*mesh%ncell/(mesh%xmax+abs(mesh%xmin)))+mesh%nghost)
        !Calculate the center of the element the source is placed in
        elementCenter = mesh%x(1, src%element) + (mesh%x(mesh%np, src%element) - mesh%x(1, src%element))/2
        !Calculate the position relative to the center of the element
        diff = src%center - elementCenter
        !scale to the size of the element
        src%r(1) = diff/((mesh%x(mesh%np, src%element) - mesh%x(1, src%element))/2)
    end subroutine

    subroutine interpolateSource(mesh, vdm, src)
        !input
        type(meshVar), intent(in)   :: mesh
        real(kind=custom_real), dimension(:,:), intent(in) :: vdm
        !in/out
        type(sourceVar) :: src
        !local
        real(kind=custom_real), dimension(mesh%np,mesh%np) :: vdmTinv  !Transposed and inverted vandermonde matrix (V^T)^-1
        real(kind=custom_real), dimension(mesh%np) :: tmp

        allocate(src%interpolation(mesh%np))
        vdmTinv = transpose(vdm)
        call invert(vdmTinv)
        call vdm1d(mesh, src%r, tmp)
        src%interpolation = matmul(vdmTinv, tmp)
    end subroutine

    subroutine rhsSource(par, mesh, op, genmat, src, tstep, src_vec)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in)   :: mesh
        type(operator), intent(in) :: op
        type(generic_material), intent(in) :: genmat
        type(sourceVar), dimension(:), intent(in) :: src
        integer, intent(in) :: tstep
        !output
        type(solutionvector) :: src_vec
        !local
        integer :: i, j
        real(kind=custom_real) :: rho
        real(kind=custom_real), dimension(:,:), allocatable :: A

        allocate(A(2+2*par%fluidn,2+2*par%fluidn))

        do i=1, par%srcn
            select case (trim(src(i)%type))
                case ("ricker", "green", "external", "sin3", "dg")
                    if (par%poroelastic) then
                        A = genmat.A.(mesh%etom(src(i)%element))
                        src_vec%v(:,src(i)%element)  = src(i)%interpolation(:) * src(i)%stf(tstep)*(-A(2,1))
                        src_vec%v1(:,src(i)%element) = src(i)%interpolation(:) * src(i)%stf(tstep)*(-A(4,1))
                        if (par%fluidn == 2) then
                            src_vec%v2(:,src(i)%element) = src(i)%interpolation(:) * src(i)%stf(tstep)*(-A(6,1))
                        endif
                    else
                        rho = genmat.rho.(mesh%etom(src(i)%element))
                        src_vec%v(:,src(i)%element) = src(i)%interpolation(:) * src(i)%stf(tstep)/rho
                    endif
                    src_vec%v(:,src(i)%element) = matmul(op%invMass, src_vec%v(:,src(i)%element))
                    do j = 1, mesh%np
                        src_vec%v(j, src(i)%element) = src_vec%v(j,src(i)%element)/op%jacobian(j, src(i)%element)
                    enddo
                    if (par%poroelastic) then
                        src_vec%v1(:,src(i)%element) = matmul(op%invMass, src_vec%v1(:,src(i)%element))
                        do j = 1, mesh%np
                            src_vec%v1(j, src(i)%element) = src_vec%v1(j,src(i)%element)/op%jacobian(j, src(i)%element)
                        enddo
                        if (par%fluidn == 2) then
                            src_vec%v2(:,src(i)%element) = matmul(op%invMass, src_vec%v2(:,src(i)%element))
                            do j = 1, mesh%np
                                src_vec%v2(j, src(i)%element) = src_vec%v2(j,src(i)%element)/op%jacobian(j, src(i)%element)
                            enddo
                        endif
                    endif
                case default
                    src_vec%v = 0.
                    if (par%poroelastic) then
                        src_vec%v1 = 0.
                        if (par%fluidn == 2) then
                            src_vec%v2 = 0.
                        endif
                    endif
            end select
        enddo
    end subroutine

    function stfRicker(t,src)
        !input
        type(sourceVar), intent(in) :: src
        real(kind=custom_real), intent(in) :: t
        !output
        real(kind=custom_real) :: stfRicker
        !local
        real(kind=custom_real) :: aval

        aval = pi*pi*src%f0*src%f0
        stfRicker =  src%amplitude * (1.-2.*aval*(t-src%t0)**2.) * exp(-aval*(t-src%t0)**2.)
    end function stfRicker

    function stfDiffRicker(t,src)
        !input
        type(sourceVar), intent(in) :: src
        real(kind=custom_real), intent(in) :: t
        !output
        real(kind=custom_real) :: stfDiffRicker
        !local
        real(kind=custom_real) :: aval

        aval = pi*pi*src%f0*src%f0
        stfDiffRicker =  2*aval*src%amplitude*(t-src%t0)*exp(-aval*(t-src%t0)**2.)*(2*aval*(t-src%t0)**2-3)
    end function stfDiffRicker

    function stfGauss(t,src)
        !input
        type(sourceVar), intent(in) :: src
        real(kind=custom_real), intent(in) :: t
        !output
        real(kind=custom_real) :: stfGauss
        !local
        real(kind=custom_real) :: aval

        aval = pi*pi*src%f0*src%f0
        stfGauss = -src%amplitude * sqrt(pi) * src%f0 * exp(-aval*(t-src%t0)**2.)
    end function stfGauss

    function stfDiffGauss(t,src)
        !input
        type(sourceVar), intent(in) :: src
        real(kind=custom_real), intent(in) :: t
        !output
        real(kind=custom_real) :: stfDiffGauss
        !local
        real(kind=custom_real) :: aval

        aval = pi*pi*src%f0*src%f0
        stfDiffGauss = -src%amplitude * sqrt(pi) * aval * 2 * (t-src%t0) * src%f0 * exp(-aval*(t-src%t0)**2.)
    end function stfDiffGauss

    function stfSin3(t,src)
        !Input
        type(sourceVar), intent(in) :: src
        real(kind=custom_real), intent(in) :: t
        !output
        real(kind=custom_real) :: stfSin3

        if (t < (1.0/src%f0)) then
            stfSin3 = -src%amplitude * sin(pi*t*src%f0)**3
        else
            stfSin3 = 0.0
        end if
    end function stfSin3

    !subroutine stfExternal(dnew,tnew,dt,nt,filter,w,nl,filename,errmsg)
    subroutine stfExternal(src,dt,nt,filter,w,nl,filename,errmsg)
        !input
        type(error_message), intent(in) :: errmsg
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nt
        integer, intent(in) :: w    ! bedeutung???????
        integer, intent(in) :: nl   ! bedeutung???????
        real(kind=CUSTOM_REAL), intent(in) :: dt
        logical, intent(in) :: filter
        !in/output
        type(sourceVar) :: src
        !local
        character(len=21) :: myname = 'stfExternal'
        integer :: io
        integer :: ns !number of line in external source file
        integer :: i
        integer :: j
        integer :: x0_f !variable to convert x0 (real) to integer
        integer :: fu
        real(kind=CUSTOM_REAL) :: rtemp !temporary variable to read the number of lines in the externals source file
        real(kind=CUSTOM_REAL) :: dt_old
        real(kind=CUSTOM_REAL) :: t0
        real(kind=CUSTOM_REAL) :: x0
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: data,time
        logical ::file_exists

        call addTrace(errmsg,myname)

        !open sourcefile
        inquire(file=trim(filename), exist=file_exists)
        if (.not.file_exists) then
            call add(errmsg,2,'File '//trim(filename)//' does not exist!',myname)
            return
        end if

        ns = 0
        io = 0

        fu = 1

        open(unit=fu,file=trim(filename),action = 'read')
        do while (io >= 0)
            read(fu,*,iostat=io) rtemp,rtemp
            ns = ns+1
        end do
        close(fu)
        ns = ns-1

        allocate(data(ns),time(ns))
        fu = 1

        open(unit=fu,file=trim(filename),action = 'read')
        do i=1,ns
            read(fu,*) time(i),data(i)
        end do
        close(fu)

        dt_old = time(2)-time(1)
        t0 = time(1)

        src%stf(:) = 0.0
        do i=1,nt
            x0 = ((i-1.)*dt) / dt_old +1
            x0_f = int(x0)
            if (i<=nl) then
                do j=1,i+nl
                    src%stf(i) = src%stf(i)+data(j)*lanczos((x0-j),nl)
                end do
            else if (x0 > (ns - nl)) then
                src%stf(i) = data(ns)
            else
                do j=(x0_f-nl),x0_f+nl
                    src%stf(i) = src%stf(i)+data(j)*lanczos((x0-j),nl)
                end do
            end if
        end do
        if (filter) then
            src%stf=movFilter(src%stf,w,nt)
        end if
        src%stf = src%stf - src%stf(1)

        deallocate(data,time)
    end subroutine stfExternal

    subroutine initNull(par, mesh, solution)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in) :: mesh
        !input/output
        type(solutionvector), intent(inout) :: solution
        !local
        integer :: i, j

        do i = 1, mesh%K
            do j = 1, mesh%Np
                solution%sigma(j,i) =      1e-20
                solution%v(j,i)     =      1e-20
                if (par%poroelastic) then
                    solution%p1(j,i)    =      1e-20
                    solution%v1(j,i)    =      1e-20
                    if (par%fluidn == 2) then
                        solution%p2(j,i)    =  1e-20
                        solution%v2(j,i)    =  1e-20
                    endif
                endif
            enddo
        enddo
    end subroutine initNull

    subroutine initRightpropBoxcar(par, mesh, src, genmat, solution)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in) :: mesh
        type(generic_material), intent(in) :: genmat
        type(sourceVar), intent(in) :: src
        !input/output
        type(solutionvector), intent(inout) :: solution
        !local
        real(kind=custom_real), dimension(:,:), allocatable :: multiplier
        integer :: i, j, m, locOnGrid

        locOnGrid = nint((src%center*mesh%ncell/mesh%xmax)+mesh%nghost)

        do i = 1, mesh%K
            do j = 1, mesh%Np
                if (i <= locOnGrid) then
                    m = mesh%ntom(j,i)
                    if (par%poroelastic) then
                        allocate(multiplier(2,2+2*par%fluidn))
                        multiplier = getSourceMultiplier(genmat.A.m)
                        solution%sigma(j,i) =  src%amplitude * multiplier(1,1)
                        solution%v(j,i)     =  src%amplitude * multiplier(1,2)
                        solution%p1(j,i)    =  src%amplitude * multiplier(1,3)
                        solution%v1(j,i)    =  src%amplitude * multiplier(1,4)
                        if (par%fluidn == 2) then
                            solution%p2(j,i)    =  src%amplitude * multiplier(1,5)
                            solution%v2(j,i)    =  src%amplitude * multiplier(1,6)
                        endif
                        deallocate(multiplier)
                    else
                        solution%v(j,i)     =  src%amplitude
                        solution%sigma(j,i) = -src%amplitude * (genmat.imp.m)
                    endif
                endif
            enddo
        enddo
    end subroutine

    subroutine initLeftpropBoxcar(par, mesh, src, genmat, solution)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in) :: mesh
        type(generic_material), intent(in) :: genmat
        type(sourceVar), intent(in) :: src
        !input/output
        type(solutionvector), intent(inout) :: solution
        !local
        real(kind=custom_real), dimension(:,:), allocatable :: multiplier
        integer :: i, j, m

        do i = 1, mesh%K
            do j = 1, mesh%Np
                if (mesh%x(j,i) >= src%center) then
                    m = mesh%ntom(j,i)
                    if (par%poroelastic) then
                        allocate(multiplier(2,2+2*par%fluidn))
                        multiplier = getSourceMultiplier(genmat.A.m)
                        solution%sigma(j,i) =  src%amplitude * multiplier(2,1)
                        solution%v(j,i)     =  src%amplitude * multiplier(2,2)
                        solution%p1(j,i)    =  src%amplitude * multiplier(2,3)
                        solution%v1(j,i)    =  src%amplitude * multiplier(2,4)
                        if (par%fluidn == 2) then
                            solution%p2(j,i)    =  src%amplitude * multiplier(2,5)
                            solution%v2(j,i)    =  src%amplitude * multiplier(2,6)
                        endif
                        deallocate(multiplier)
                    else
                        solution%v(j,i)     = src%amplitude
                        solution%sigma(j,i) = src%amplitude * (genmat.imp.m)
                    endif
                endif
            enddo
        enddo
    end subroutine

    subroutine initRightpropSineSquared(par, mesh, src, genmat, solution)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in) :: mesh
        type(generic_material), intent(in) :: genmat
        type(sourceVar), intent(in) :: src
        !input/output
        type(solutionvector), intent(inout) :: solution
        !local
        integer :: i, j, m
        real(kind=custom_real), dimension(:,:), allocatable :: multiplier
        real(kind=custom_real) :: xc
        real(kind=custom_real) :: w

        xc = src%center
        w = src%width

        do i = 1, mesh%K
            do j = 1, mesh%Np
                if (abs(mesh%x(j,i)-xc) < 0.5*w) then
                    m = mesh%ntom(j,i)
                    if (par%poroelastic) then
                        allocate(multiplier(2,2+2*par%fluidn))
                        multiplier = getSourceMultiplier(genmat.A.m)
                        solution%sigma(j,i) =  src%amplitude * multiplier(1,1) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%v(j,i)     =  src%amplitude * multiplier(1,2) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%p1(j,i)    =  src%amplitude * multiplier(1,3) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%v1(j,i)    =  src%amplitude * multiplier(1,4) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        if (par%fluidn == 2) then
                            solution%p2(j,i)    =  src%amplitude * multiplier(1,5) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                            solution%v2(j,i)    =  src%amplitude * multiplier(1,6) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        endif
                        deallocate(multiplier)
                    else
                        solution%v(j,i)     =  src%amplitude                  * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%sigma(j,i) = -src%amplitude * (genmat.imp.m) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                    endif
                endif
            enddo
        enddo
    end subroutine

    subroutine initLeftpropSineSquared(par, mesh, src, genmat, solution)
        !input
        type(parameterVar), intent(in) :: par
        type(meshVar), intent(in) :: mesh
        type(generic_material), intent(in) :: genmat
        type(sourceVar), intent(in) :: src
        !input/output
        type(solutionvector), intent(inout) :: solution
        !local
        integer :: i, j, m
        real(kind=custom_real), dimension(:,:), allocatable :: multiplier
        real(kind=custom_real) :: xc
        real(kind=custom_real) :: w

        xc = src%center
        w = src%width

        do i = 1, mesh%K
            do j = 1, mesh%Np
                if (abs(mesh%x(j,i)-xc) < 0.5*w) then
                    m = mesh%ntom(j,i)
                    if (par%poroelastic) then
                        allocate(multiplier(2,2+2*par%fluidn))
                        multiplier = getSourceMultiplier(genmat.A.m)
                        solution%sigma(j,i) =  src%amplitude * multiplier(2,1) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%v(j,i)     =  src%amplitude * multiplier(2,2) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%p1(j,i)    =  src%amplitude * multiplier(2,3) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%v1(j,i)    =  src%amplitude * multiplier(2,4) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        if (par%fluidn == 2) then
                            solution%p2(j,i)    =  src%amplitude * multiplier(2,5) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                            solution%v2(j,i)    =  src%amplitude * multiplier(2,6) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        endif
                        deallocate(multiplier)
                    else
                        solution%v(j,i)     =  src%amplitude                  * cos(pi*(mesh%x(j,i)-xc)/w)**2
                        solution%sigma(j,i) = +src%amplitude * (genmat.imp.m) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                    endif
                endif
            enddo
        enddo
    end subroutine

    subroutine initRightpropTestwavelet(mesh, src, genmat, V, Sigma)
        !input
        type(meshVar), intent(in) :: mesh
        type(generic_material), intent(in) :: genmat
        type(sourceVar), intent(in) :: src
        !input/output
        real(kind=custom_real), dimension(:,:), intent(inout) :: V
        real(kind=custom_real), dimension(:,:), intent(inout) :: Sigma
        !local
        integer :: i, j, m
        real(kind=custom_real) :: w
        real(kind=custom_real) :: xc

        xc = src%center
        w = src%width

        do i = 1, mesh%K
            do j = 1, mesh%Np
                if (mesh%x(j,i) < xc) then
                    m = mesh%ntom(j,i)
                    V(j,i)     = +src%amplitude
                    Sigma(j,i) = -src%amplitude * (genmat.imp.m)
                else if (mesh%x(j,i) >= xc .and. mesh%x(j,i) < xc + w/2) then
                    m = mesh%ntom(j,i)
                    V(j,i)     = +src%amplitude                  * cos(pi*(mesh%x(j,i)-xc)/w)**2
                    Sigma(j,i) = -src%amplitude * (genmat.imp.m) * cos(pi*(mesh%x(j,i)-xc)/w)**2
                else
                    V(j,i) = 0
                    Sigma(j,i) = 0
                endif
            enddo
        enddo
    end subroutine

    function getSourceMultiplier(A) result(multiplier)
        !input
        real(kind=custom_real), dimension(:,:), intent(in) :: A
        !output
        real(kind=custom_real), dimension(:,:), allocatable :: multiplier
        !local variables
        type(error_message) :: errmsg
        character(len=19) :: myname = 'getSourceMultiplier'
        character(len=100) :: errstr
        real(kind=custom_real), dimension(:,:), allocatable :: Awork
        real(kind=custom_real), dimension(:), allocatable :: WR, WI    !contain real and imaginary part of eigenvalues, respectively
        real(kind=custom_real), dimension(:,:), allocatable :: VR      !matrix containing the eigenvectors as columns
        real(kind=custom_real), dimension(:,:), allocatable :: VL
        real(kind=custom_real), dimension(:), allocatable :: work
        integer :: N, lwork, info
        integer, dimension(:), allocatable :: ipiv, posmax, posmin
        integer, dimension(2) :: shapeA

        shapeA=shape(A)
        N=shapeA(1)
        allocate(multiplier(2,N))
        allocate(posmax(N))
        allocate(posmin(N))

        allocate(Awork(N,N))
        allocate(ipiv(N))
        allocate(WR(N))
        allocate(WI(N))
        allocate(VR(N,N))
        Awork = A

        !calculate eigenvectors and eigenvalues
        if (custom_real==size_double) then
            !do workspace query
            allocate(work(1)); lwork = -1
            !SYNTAX: dgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )   (LAPACK: http://www.netlib.org/lapack/double/dgeev.f)
            call dgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Workspace query failed: LAPACK routine DGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            !compute eigenvalues and eigenvectors
            call dgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Computation of eigenvalues and -vectors failed: LAPACK routine DGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            deallocate(work)
        elseif (custom_real==size_real) then
            !do workspace query
            allocate(work(1)); lwork = -1
            !SYNTAX: sgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )   (LAPACK: http://www.netlib.org/lapack/double/dgeev.f)
            call sgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Workspace query failed: LAPACK routine SGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            !compute eigenvalues and eigenvectors
            call sgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Computation of eigenvalues and -vectors failed: LAPACK routine SGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            deallocate(work)
        else
            call add(errmsg,2,'CUSTOM_REAL is neither SIZE_REAL nor SIZE_DOUBLE',myname)
            stop
        endif

        multiplier = 0.
        posmax = maxloc(WR)
        posmin = minloc(WR)
        multiplier(1,:) = VR(:,posmax(1))
        multiplier(2,:) = VR(:,posmin(1))
        multiplier(1,:) = multiplier(1,:) / multiplier(1,2) !set entry of velocity v to 1
        multiplier(2,:) = multiplier(2,:) / multiplier(2,2) !set entry of velocity v to 1

        if (allocated(Awork)) deallocate(Awork)
        if (allocated(ipiv)) deallocate(ipiv)
        if (allocated(WR)) deallocate(WR)
        if (allocated(WI)) deallocate(WI)
        if (allocated(VR)) deallocate(VR)
    end function getSourceMultiplier

    function lanczos(x,l)
        !input
        integer, intent(in) :: l
        real(kind=custom_real), intent(in) :: x
        !output
        real(kind=custom_real) :: lanczos

        if (x < epsilon(x)) then
            lanczos = 1.0
        else if (abs(x) < l) then
            lanczos = sin(pi*x)/(pi*x) * sin((pi*x)/l)/((pi*x)/l)
        else
            lanczos= 0.0
        end if
    end function lanczos

    function movFilter(data,w,nt)
        !input
        integer, intent(in) :: nt
        integer, intent(in) :: w
        real(kind=CUSTOM_REAL), dimension(nt), intent(in) :: data
        !output
        real(kind=CUSTOM_REAL) , dimension(nt) :: movFilter
        !local
        integer :: i,j,n,m
        real(kind=CUSTOM_REAL) :: temp

        do i=1,nt
            if (i<=w) then
                n=1
                m=i+w
            else if ( i >= nt-w) then
                n=i-w
                m=nt-1
            else
                n=i-w
                m=i+w
            end if
                temp=0.0
            do j=n,m
                temp=temp+data(j)
            end do
            movFilter(i) = temp / ((2*w)+1)
        end do
    end function movFilter
end module
