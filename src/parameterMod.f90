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
module parameterMod
    use constantsMod
    use errorMessage

    type :: parameterVar
        character(len=80) :: title          !Define the title of the simulation
        character(len=80) :: limiter        !Type of limiter: Minmod, superbee or MC
        integer :: tsteps                   !Total number of timesteps
        integer :: framestep                !intervall in which to create a outputfile
        integer :: tint                     !Time integragtion: Euler, SSP-RK und RK4
        integer :: matn                     !Number of (different) Materials
        integer :: srcn                     !Number of Sources
        character(len=80) :: extwavelet     !File with external source wavelet
        integer :: lsin                     !Number of Interfaces
        integer :: recn                     !Number of Receivers
        integer :: flux_type                !Parameter to chose how the fluxes are calculated: "0" Precalculated for elastic, "1" Godunov, "2" Rusanov
        real(kind=custom_real) :: cfl       !Courrantnumber
        real(kind=custom_real) :: dt        !timestep
        logical :: autodt                   !calculate timestep automatically or not
        logical :: log                      !enable to log on screen
        logical :: logfile                  !write log to file
        logical :: debug                    !enable to set a debugflag (currently not in use)
        logical :: movie                    !select if output is created
        logical :: poroelastic              !Materials are poroelastic if .true. otherwise elastic
        integer :: fluidn                   !Number of immiscible fluids (either 1, i.e. saturated, or 2, i.e. unsaturated/saturated by 2 fluids)
        logical :: calculate_tortuosity     !Tortuosity is calculated according to Berryman (1980): T = 1+r(1-1/phi) (note, that if this is set to .true., in the file porousmaterial r has to be specified instead of T!)
        logical :: physical_coordinates     !Select fortmat of input for the coordinates of the material boundaries and the location of the LSI
        logical :: compare                  !true if numerical solution should be compared to an analytical.
    endtype parameterVar

    type :: meshVar
        character(len=80) :: lbc                                !Boundary condition for the left side
        character(len=80) :: rbc                                !Boundary condition for the right side
        integer :: ncell                                        !Number of regular cells
        integer :: nghost                                       !Number of Ghostelements at each side of the grid
        integer :: nfaces                                       !Number of faces for each element
        integer :: nlsi                                         !Number of Slipinterface nodes
        integer :: K                                            !Total number of elements
        integer :: nv                                           !Number of Vertexpoints
        integer :: n                                            !Order of approximation
        integer :: np                                           !Number of points per element
        integer :: nfp                                          !Number of face nodes
        integer, dimension(:), allocatable :: etom              !Material index at each element -> element to material
        integer, dimension(:), allocatable :: ntof              !entry i gives node index sitting on face i, named Fmask in Hesthaven & Warburton (2008)
        integer, dimension(:,:), allocatable :: nx              !normales at the elements faces
        integer, dimension(:,:), allocatable :: ntom            !material index at nodes -> node to material
        integer, dimension(:,:), allocatable :: etov            !element to vertices map. entry (i,j) gives index of j-th vertex of cell i
        integer, dimension(:,:), allocatable :: etoe            !Mapping of the neighboring elements
        integer, dimension(:,:), allocatable :: etof            !element face to face map, (ne x 2),
                                                                ! entry (i,j) gives local index of face adjacent to element i's j-th face
        real(kind=custom_real) :: xmax                          !Maximum x-value
        real(kind=custom_real) :: xmin                          !Minimum x-value
        real(kind=custom_real), dimension(:), allocatable :: r  !Vector of the gll points
        real(kind=custom_real), dimension(:), allocatable :: vx !array of vertex coordinates (including ghost elements)
        real(kind=custom_real), dimension(:), allocatable :: h  !element width; i-th entry gives width of element i
        real(kind=custom_real), dimension(:,:), allocatable :: x!pysical grid points
    endtype meshVar

    type :: solutionvector
        real(kind=custom_real), dimension(:,:), allocatable :: sigma    !stress of solid
        real(kind=custom_real), dimension(:,:), allocatable :: v        !velocity of solid
        real(kind=custom_real), dimension(:,:), allocatable :: p1       !pressure of fluid 1
        real(kind=custom_real), dimension(:,:), allocatable :: v1       !velocity of fluid 1
        real(kind=custom_real), dimension(:,:), allocatable :: p2       !pressure of fluid 2
        real(kind=custom_real), dimension(:,:), allocatable :: v2       !velocity of fluid 2
    endtype solutionvector

    contains

    subroutine readParfile(this, mesh ,errmsg)
        type (parameterVar) :: this
        type (meshVar) :: mesh
        type (error_message) :: errmsg
        !local variables
        character(len=80) :: filename
        character(len=80) :: outfile
        character(len=20) :: myname = "parameterMod"
        !integer :: i
        integer :: ier
        integer :: out

        call addTrace(errmsg,myname)
        filename=trim('input/parfile')
        !initial test if parfile exists
        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        if (ier /= 0) then
            call add(errmsg,2,'Could not open file: '//trim(filename),myname)
            call print(errmsg)
            return
        endif
        close(19)

        ! Cycle to read the parameters form the parameter file. The order of appearence of the parameters is not important
        call readStringPar(this%title, "title", filename, 0)
        call readStringPar(this%limiter, "limiter", filename, 0)
        call readStringPar(mesh%lbc, "lbc", filename, 0)
        call readStringPar(mesh%rbc, "rbc", filename, 0)
        call readIntPar(this%framestep, "framestep", filename, 0)
        call readIntPar(this%tsteps, "Nsteps", filename, 0)
        call readIntPar(this%tint, "tint", filename, 0)
        call readLogicalPar(this%autodt, "autodt", filename, 0)
        call readFloatPar(this%dt, "dt", filename, 0)
        call readIntPar(mesh%ncell, "ncell", filename, 0)
        call readIntPar(mesh%n, "n", filename, 0)
        call readIntPar(mesh%nghost, "nghost", filename, 0)
        call readIntPar(mesh%nfaces, "nfaces", filename, 0)
        call readIntPar(mesh%nfp, "nfp", filename, 0)
        call readIntPar(this%matn, "matn", filename, 0)
        call readIntPar(this%lsin, "lsin", filename, 0)
        call readIntPar(this%srcn, "srcn", filename, 0)
        call readStringPar(this%extwavelet, "extwavelet", filename, 0)
        call readIntPar(this%recn, "recn", filename, 0)
        call readIntPar(this%flux_type, "flux_type", filename, 0)
        call readFloatPar(this%cfl, "cfl", filename, 0)
        call readFloatPar(mesh%xmin, "xmin", filename, 0)
        call readFloatPar(mesh%xmax, "xmax", filename, 0)
        call readLogicalPar(this%log, "log", filename, 0)
        call readLogicalPar(this%logfile, "logfile", filename, 0)
        call readLogicalPar(this%movie, "movie", filename, 0)
        call readLogicalPar(this%debug, "debug", filename, 0)
        call readLogicalPar(this%poroelastic, "poroelastic", filename, 0)
        if (this%poroelastic) then
            call readIntPar(this%fluidn, "fluidn", filename, 0)
        else
            this%fluidn = 0
        endif
        call readLogicalPar(this%calculate_tortuosity, "calculate_tortuosity", filename, 0)
        call readLogicalPar(this%physical_coordinates, "physical_coordinates", filename, 0)
        call readLogicalPar(this%compare, "compare", filename, 0)

        !test if certain parameters have been read/entered correctly
        if (mesh%nfaces /= 2) call add(errmsg,2,'nfaces must be 2',myname)
        if (mesh%nfp /= 1) call add(errmsg,2,'nfp must be 1',myname)
        if (this%tint /= 1 .and. this%tint /= 5 .and. this%tint /= 3) then
            call add(errmsg,2,'tint must be either 1,3 or 5',myname)
        endif
        if (this%compare .and. this%poroelastic) then
            call add(errmsg,2,'No analytical solution for poroelastic wave propagation available.',myname)
        endif
        if (trim(mesh%lbc) == "periodic" .and. trim(mesh%rbc) /= "periodic") then
            call add(errmsg,2,'Left AND right boundary conditions must be periodic',myname)
        endif
        if (trim(mesh%rbc) == "periodic" .and. trim(mesh%lbc) /= "periodic") then
            call add(errmsg,2,'Left AND right boundary conditions must be periodic',myname)
        endif
        if (trim(mesh%rbc) == "periodic" .and. trim(mesh%lbc) == "periodic") then
            if (mesh%nghost /= 0) then
                call add(errmsg,2,'In case of periodic BC no ghost elements are needed. Should be set to 0!',myname)
            endif
        else if (trim(mesh%rbc) == "absorbing" .or. trim(mesh%lbc) == "absorbing") then
            if (mesh%nghost < 2) then
                call add(errmsg,2,'In case of absorbing BC more than one ghost element is needed. Select 2!',myname)
            endif
        else
            if (mesh%nghost == 0) then
                call add(errmsg,2,'In case of reflecting BC ghost elements are needed. Should be set > 0!',myname)
            endif
        endif
        if (this%poroelastic .and. this%flux_type == 0) then
              call add(errmsg,2,'For poroelastic simulations select flux 1 or 2',myname)
        endif


        !calculate additional mesh parameters from input
        mesh%K = mesh%ncell + mesh%nfaces*mesh%nghost
        mesh%nv = mesh%K + 1
        mesh%np = mesh%n + 1

        if (this%log) then
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, a10, a30)")   "|              Title of the simulation: ", this%title, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Create log: ", this%log, "                             |"
            write (*,"(a40, l10, a30)")   "|                           Debug mode: ", this%debug, "                             |"
            write (*,"(a40, a10, a30)")   "|                        Slope Limiter: ", this%limiter, "                             |"
            write (*,"(a40, l10, a30)")   "|       Compare to analytical solution: ", this%compare, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                 Mesh Parameters                              |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, f10.1, a30)") "|                                 xmin: ", mesh%xmin, "                             |"
            write (*,"(a40, f10.1, a30)") "|                                 xmax: ", mesh%xmax, "                             |"
            write (*,"(a40, i10, a30)")   "|                   Number of elements: ", mesh%ncell, "                             |"
            write (*,"(a40, i10, a30)")   "|               Order of approximation: ", mesh%n, "                             |"
            write (*,"(a40, i10, a30)")   "|                  Ghost-elements/side: ", mesh%nghost, "                             |"
            write (*,"(a40, i10, a30)")   "|              Number of faces/element: ", mesh%nfaces, "                             |"
            write (*,"(a40, i10, a30)")   "|                Number of points/face: ", mesh%nfp, "                             |"
            write (*,"(a40, a13, a27)")   "|       Boundary condition (left side): ", mesh%lbc, "                          |"
            write (*,"(a40, a13, a27)")   "|      Boundary condition (right side): ", mesh%rbc, "                          |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                Movie Parameters                              |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                           Save movie: ", this%movie, "                             |"
            write (*,"(a40, i10, a30)")   "|   Number of timesteps for movieframe: ", this%framestep,&
             "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                   Environment                                |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, l10, a30)")   "|                          Poroelastic: ", this%poroelastic,&
             "                             |"
            if (this%poroelastic) then
                write (*,"(a40, i10, a30)")   "|          Number of immiscible fluids: ", this%fluidn, "                             |"
                write (*,"(a40, l10, a30)")   "|                 Calculate tortuosity: ", this%calculate_tortuosity,&
                 "                             |"
            endif
            if (this%physical_coordinates) then
                write (*,"(a80)")   "|                     Coordinate input:  physical                              |"
            else
                write (*,"(a80)")   "|                     Coordinate input:  Element number                        |"
            endif
            write (*,"(a40, i10, a30)")   "|                            Materials: ", this%matn, "                             |"
            write (*,"(a40, i10, a30)")   "|                              Sources: ", this%srcn, "                             |"
            write (*,"(a40, i10, a30)")   "|                             Receiver: ", this%recn, "                             |"
            write (*,"(a40, i10, a30)")   "|                      Slip Interfaces: ", this%lsin, "                             |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*, "(a80)") "|                                 Timeintegration                              |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
            write (*,"(a40, i10, a30)")   "|                                 Type: ", this%tint, "                             |"
            write (*,"(a40, i10, a30)")   "|                  Number of timesteps: ", this%tsteps, "                             |"
            write (*,"(a40, l10, a30)")   "|                               autodt: ", this%autodt, "                             |"
            if (this%autodt) then
                write (*,"(a40, f10.2, a30)") "|                            cfl value: ", this%cfl, "                             |"
            endif
        end if
        if (this%logfile) then
            outfile = "output/logfile"//trim(this%title)
            outfile = trim(outfile)
            out = 16
            open(out, file = outfile, position="append")
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a40, a10, a30)")   "|              Title of the simulation: ", this%title, "                             |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a40, l10, a30)")   "|                           Create log: ", this%log, "                             |"
            write (out,"(a40, l10, a30)")   "|                           Debug mode: ", this%debug, "                             |"
            write (out,"(a40, a10, a20, a10)")   "|                        Slope Limiter: ", this%limiter, "                   ",&
            "         |"
            write (out,"(a40, l10, a20, a10)")   "|       Compare to analytical solution: ", this%compare, "                   ",&
            "         |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out, "(a80)") "|                                 Mesh Parameters                              |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a40, f10.1, a30)") "|                                 xmin: ", mesh%xmin, "                             |"
            write (out,"(a40, f10.1, a30)") "|                                 xmax: ", mesh%xmax, "                             |"
            write (out,"(a40, i10, a30)")   "|                   Number of elements: ", mesh%ncell, "                             |"
            write (out,"(a40, i10, a30)")   "|               Order of approximation: ", mesh%n, "                             |"
            write (out,"(a40, i10, a20, a10)")   "|                  Ghost-elements/side: ", mesh%nghost, "                    ",&
            "         |"
            write (out,"(a40, i10, a20, a10)")   "|              Number of faces/element: ", mesh%nfaces, "                    ",&
            "         |"
            write (out,"(a40, i10, a30)")   "|                Number of points/face: ", mesh%nfp, "                             |"
            write (out,"(a40, a13, a27)")   "|       Boundary condition (left side): ", mesh%lbc, "                          |"
            write (out,"(a40, a13, a27)")   "|      Boundary condition (right side): ", mesh%rbc, "                          |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out, "(a80)") "|                                Movie Parameters                              |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a40, l10, a30)")   "|                           Save movie: ", this%movie, "                             |"
            write (out,"(a40, i10, a30)")   "|   Number of timesteps for movieframe: ", this%framestep,&
             "                             |"

            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out, "(a80)") "|                                   Environment                                |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a40, l10, a30)")   "|                          Poroelastic: ", this%poroelastic,&
             "                             |"
            if (this%poroelastic) then
                write (out,"(a40, i10, a20, a10)")   "|          Number of immiscible fluids: ", this%fluidn, "                    ",&
                 "         |"
                write (out,"(a40, l10, a30)")   "|                 Calculate tortuosity: ", this%calculate_tortuosity,&
                 "                             |"
            endif
            if (this%physical_coordinates) then
                write (out,"(a80)")   "|                     Coordinate input:  physical                              |"
            else if (.not. this%physical_coordinates ) then
                write (out,"(a80)")   "|                     Coordinate input:  Element number                        |"
            endif
            write (out,"(a40, i10, a30)")   "|                            Materials: ", this%matn, "                             |"
            write (out,"(a40, i10, a30)")   "|                              Sources: ", this%srcn, "                             |"
            write (out,"(a40, i10, a30)")   "|                             Receiver: ", this%recn, "                             |"
            write (out,"(a40, i10, a30)")   "|                      Slip Interfaces: ", this%lsin, "                             |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out, "(a80)") "|                                 Timeintegration                              |"
            write (out, "(a80)") "|------------------------------------------------------------------------------|"
            write (out,"(a40, i10, a30)")   "|                                 Type: ", this%tint, "                             |"
            write (out,"(a40, i10, a20, a10)")   "|                  Number of timesteps: ", this%tsteps, "                    ",&
                       "         |"
            write (out,"(a40, l10, a23, a7)") "|                               autodt: ", this%autodt, "                       ",&
                       "      |"
            if (this%autodt) then
                write (out,"(a40, f10.2, a20, a10)") "|                            cfl value: ", this%cfl, "                    ",&
                           "         |"
            endif
        end if
    end subroutine readParfile

    function getParameterName(line) result(par)
        !Function to call and compare the parameter name to the entry in the file -> find the loocation of a specific parameter!
        character (len=*), intent(in) :: line
        character (len=len(line)) :: tmp
        character (len=50) :: par

        integer, parameter :: tab = 9, space = 32, eq = 61
        integer :: i, ASCII_char, n

        n = 0
        par = " "
        tmp = line(1:index(line, "="))
        do i=1, len(tmp)
            ASCII_char = iachar(tmp(i:i))
            if (ASCII_char /= tab .and. ASCII_char /= space .and. ASCII_char /= eq) then
                n = n + 1
                par(n:n) = tmp(i:i)
            end if
        end do
    end function getParameterName

    subroutine getParameterArray(name, filename, pos, value)
        !Read the parameters value from the input string
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: pos

        character (len=80) :: par, tmp
        character (len=255) :: line
        character (len=80), dimension(:) :: value

        integer :: i, j, l, ier
        integer, parameter :: seek_set = 0
        logical :: found = .false.
        integer :: set_pos

        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        call fseek(19, pos, seek_set, ier)
        set_pos = int(ftell(19))

        do while (.not. is_iostat_end(ier))
            read(19, "(a255)", iostat=ier) line
            line = trim(adjustl(line))
            if (line(1:1) == "#") cycle
            if (line(1:1) == "\n") cycle
            par = getParameterName(line)
            par = trim(par)
            j = scan(line, "=")
            l = scan(line, "#")
            if (par /= name) cycle
            if (par == name) then
                tmp = trim(line(j+1:l-1))
                if (tmp == "") then
                    value(:) = " "
                else
                    read(tmp,*) (value(i),i=1,size(value))
                    do i = 1, size(value)
                        value(i) = trim(value(i))
                    enddo
                    found = .true.
                    exit
                endif
            endif
        enddo
    end subroutine getParameterArray

    subroutine readIntArray(array_to_read, name, filename, pos)
        !function to read interger values from the parameter-file
        !input
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: pos
        !output
        integer, dimension(:), intent(out) :: array_to_read
        !local
        integer :: i
        character (len=80), dimension(size(array_to_read)) :: value
        logical :: emptyval = .false.

        call getParameterArray(name, filename, pos, value)

        do i = 1, size(value)
            if (value(i) /= " ") then
                read(value(i), *) array_to_read(i)
            else
                emptyval = .true.
                cycle
            endif
        enddo
    end subroutine readIntArray

    subroutine readFloatArray(array_to_read, name, filename, pos)
        !function to read interger values from the parameter-file
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: pos
        !output
        real(kind=custom_real), dimension(:), intent(out) :: array_to_read
        !local
        integer :: i
        character (len=80), dimension(size(array_to_read)) :: value

        call getParameterArray(name, filename, pos, value)

        do i = 1, size(value)
            read(value(i), *) array_to_read(i)
        enddo
    end subroutine readFloatArray

    function getParameterValue(name, filename, pos) result(value)
        !Read the parameters value from the input string
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: pos

        character (len=80) :: tmp, par, value
        character (len=255) :: line
        integer :: j,l, ier
        integer, parameter :: seek_set = 0
        integer(kind=8) :: set_pos
        logical :: found = .false.

        open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
        call fseek(19, pos, seek_set, ier)
        set_pos = ftell(19)

        do while (.not. is_iostat_end(ier))
            read(19, "(a255)", iostat=ier) line
            line = trim(adjustl(line))
            if (line(1:1) == "#") cycle
            if (line(1:1) == "\n") cycle
            par = getParameterName(line)
            par = trim(par)
            j = scan(line, "=")
            l = scan(line, "#")
            if (par /= name) cycle
            if (par == name) then
                !print *, name
                if (trim(line(j+1:l-1)) == '' ) then
                    print *, "string is empty"
                    value = " "
                else
                    tmp = trim(line(j+1:l-1))
                    read(tmp, "(a80)") value
                    value = trim(value)
                endif
                found = .true.
            endif
            exit
        enddo
        close(19)
    end function getParameterValue

    subroutine readStringPar(par_to_read, name, filename, pos)
        !Function to read Parametrs as text.
        character (len=*) :: name, par_to_read, filename
        character (len=80) :: value!, tmp
        integer :: pos

        value = getParameterValue(name, filename, pos)

        if (value /= " ") then
            read(value, "(a80)") par_to_read
            par_to_read = trim(adjustl(par_to_read))
        endif
    end subroutine readStringPar

    subroutine readIntPar(par_to_read, name, filename, pos)
        !function to read interger values from the parameter-file
        character (len=*), intent(in) :: name, filename
        integer :: par_to_read, pos
        character (len=80) :: value

        value = getParameterValue(name, filename, pos)
        if (value /= " ") then
            read(value, *) par_to_read
        endif
    end subroutine readIntPar

    subroutine readFloatPar(par_to_read, name, filename, pos)
        !Function to read floating point variables from the parameter-file
        character (len=*), intent(in) :: name, filename
        real(kind=CUSTOM_REAL) :: par_to_read
        character (len=80) :: value
        integer :: pos

        value = getParameterValue(name, filename, pos)
        if (value /= " ") then
            read(value, *) par_to_read
        endif
    end subroutine readFloatPar

    subroutine readLogicalPar(par_to_read, name, filename, pos)
        !Function to read logical variables from the parameter-file
        character (len=*), intent(in) :: name, filename
        logical :: par_to_read
        character (len=80) :: value
        integer :: pos

        value = getParameterValue(name, filename, pos)
        if (value /= " ") then
            read(value, *) par_to_read
        endif
    end subroutine readLogicalPar

    subroutine deallocMeshArrays(mesh)
        type(meshVar) :: mesh
        if (allocated(mesh%etom)) deallocate(mesh%etom)
        if (allocated(mesh%ntom)) deallocate(mesh%ntom)
        if (allocated(mesh%nx)) deallocate(mesh%nx)
        if (allocated(mesh%ntof)) deallocate(mesh%ntof)
        if (allocated(mesh%etoe)) deallocate(mesh%etoe)
        if (allocated(mesh%etof)) deallocate(mesh%etof)
        if (allocated(mesh%etov)) deallocate(mesh%etov)
        if (allocated(mesh%vx)) deallocate(mesh%vx)
        if (allocated(mesh%x)) deallocate(mesh%x)
        if (allocated(mesh%h)) deallocate(mesh%h)
    end subroutine
end module
