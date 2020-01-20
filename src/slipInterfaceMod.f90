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
!--------------------------------------------------------------------
!> \brief Handle slip interface properties characterized by eta and xi
!
module slipInterfaceMod
    use errorMessage
    use constantsMod
    use parameterMod
    use genericMaterial

    implicit none

    type :: interface_spec
        character(len=10) :: type                                   !Interface types: Choices are Elastic, visco and SLS
        real(kind=custom_real) :: nu                                !Relaxation frequency of the fracture
    end type

    type :: lsiVar !für jedes echte interface
        integer :: prop                                             !Propertyindex of Intersfaces
        integer, dimension(2) :: elements                           !Grid elements associated with the interface
        real(kind=custom_real) :: loc                               !Physical location of the interface on the grid
        real(kind=custom_real), dimension(2) :: S                   !Storage variable for the Slipinterfaceactivity
        real(kind=custom_real), dimension(2) :: rhsS                !right hand side of the evolution equation for S
        real(kind=custom_real), dimension(2) :: resS                !Runge-Kutta residual (S)
    end type lsiVar

    contains

    subroutine setup_lsi(par, mesh, genmat, lsi, lsi_spec, V, errmsg)
        !input/output
        type(parameterVar) :: par
        type(meshVar) :: mesh
        type(generic_material) :: genmat
        type(error_message) :: errmsg
        type(lsiVar), dimension(:) :: lsi
        type (interface_spec), dimension(:), allocatable :: lsi_spec
        real(kind=custom_real), dimension(:,:) :: V
        !local
        character(len=15) :: myname = 'setup_lsi'
        integer :: i,j,l,m,nif,nmf
        real(kind=custom_real) :: tmp
        real(kind=custom_real) :: zin
        real(kind=custom_real) :: zout
        real(kind=custom_real) :: nu

        call addTrace(errmsg,myname)
        !Read the specifications of all possible interface types
        call readSlipInterfaceSpecs(lsi_spec, 1, "input/interfaces" ,errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        !read location and property-index of the Slip Interface
        call readSlipInterfaceProperties(par, mesh, lsi, lsi_spec, 1, "input/lsi_eap" , errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        !Set Initial activity of the interfaces
        do l = 1, par%lsin
            do j = 1, mesh%nfaces
                i = lsi(l)%elements(j)    !element
                m = mesh%etoe(i,j)        !neighbouring element
                !node indices on faces
                nmf = mesh%ntof(mesh%etof(i,j))
                nif = mesh%ntof(j)
                select case (trim(lsi_spec(lsi(l)%prop)%type))
                    case ("elastic")
                        lsi(l)%S(j) = (V(nmf, m) - V(nif,i))*mesh%nx(j,i)
                    case default
                        call add(errmsg,2,"No accepted interface-type entered. Abort...",myname)
                        call print(errmsg); return
                end select
                lsi(l)%resS(j) = 0.
                lsi(l)%rhsS(j) = 0.
            enddo
        enddo
    end subroutine

    subroutine readSlipInterfaceProperties(par, mesh, lsi, lsi_spec, lu, filename, errmsg)
        !input
        type (parameterVar) :: par
        type (meshVar) :: mesh
        type (error_message) :: errmsg
        type (interface_spec), dimension(:) :: lsi_spec
        integer :: lu
        character(len=*), intent(in) :: filename
        !in/output
        type(lsiVar), dimension(:) :: lsi
        !local
        integer :: i,j,l,nf, ios, out
        integer :: locOnGrid
        character(len=15) :: myname = 'readLSI'
        character(len=255) :: text
        character(len=10) :: lsii
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
        if (par%logfile) then
            outfile = "output/logfile"//trim(par%title)
            outfile = trim(outfile)
            out = 16
            open(out, file = outfile, position="append")
            write (out,"(a80)") "|                            Data on Interfaces:                               |"
            write (out,"(a80)") "|------------------------------------------------------------------------------|"
        endif
        if (par%log) write (*,"(a80)") "|                            Data on Interfaces:                               |"
        if (par%log) write (*,"(a80)") "|------------------------------------------------------------------------------|"
        do i = 1, par%lsin
            read(lu, *) lsi(i)%loc, lsi(i)%prop

            !check whether material specified in 'mat_bap' is available in 'materials' or not
            if (lsi(i)%prop > size(lsi_spec)) then
                write (lsii,'(i10)') lsi(i)%prop
                call add(errmsg,2,"LSI "//trim(adjustl(lsii))//" specified in 'lsi_eap' is not available in 'interfaces'.",myname)
                return
            endif

            j = lsi(i)%prop
            if (par%physical_coordinates) then
                if (lsi(i)%loc > mesh%xmax) then
                    call add(errmsg,2,'Slip interface not on grid. Abort...',myname)
                    return
                endif
                locOnGrid = nint(((lsi(i)%loc+abs(mesh%xmin))*mesh%ncell/(mesh%xmax+abs(mesh%xmin)))+mesh%nghost)
            else
                locOnGrid = nint(lsi(i)%loc) + mesh%nghost
            endif
            do nf = 1, mesh%nfaces
                ! interface:   Face 2 of element x | face 1 of element x + 1
                lsi(i)%elements(1) = locOnGrid + 1
                lsi(i)%elements(2) = locOnGrid
            enddo

            select case(trim(lsi_spec(j)%type))
                case ("elastic")
                    if (par%log)  then
                        write (*,"(a80)") "|       No.      |         type        |        Loc.       |         nu        |"
                        write (*,"(a80)") "|----------------|---------------------|-------------------|-------------------|"
                        write (*,"(a7,i4.4,a14,a7,a14,f8.2,a12,f8.2,a6 )") &
                        "|      ", i, "      |       ", lsi_spec(j)%type, "       |      ", lsi(i)%loc, "     |      ", lsi_spec(j)%nu,"     |"
                    endif
                    if (par%logfile) then
                        write (out,"(a80)") &
                        "|       No.      |         type        |        Loc.       |         nu        |"
                        write (out,"(a80)") &
                        "|----------------|---------------------|-------------------|-------------------|"
                        write (out,"(a7,i4.4,a14,a7,a14,f8.2,a12,f8.2,a6 )") &
                        "|      ", i, "      |       ", lsi_spec(j)%type, "       |      ", lsi(i)%loc, "     |      ", lsi_spec(j)%nu,"     |"
                    endif
                case default
                    call add(errmsg,2,"No accepted interface-type entered. Abort...",myname)
                    call print(errmsg); return
            end select
            if (par%log) write (*,"(a80)") "|------------------------------------------------------------------------------|"
            if (par%logfile) then
                write (out,"(a80)") "|------------------------------------------------------------------------------|"
                close(out)
            endif
        enddo
        close(lu)
    end subroutine readSlipInterfaceProperties

!---------------------------------------------------------------------
!> \brief Read slip interface specifications from a file
!
    subroutine readSlipInterfaceSpecs(lsi_spec,lu,filename,errmsg)
        !input
        integer :: lu
        character(len=*) :: filename
        type (error_message) :: errmsg
        !output
        type(interface_spec), dimension(:), allocatable :: lsi_spec
        !local
        character(len=27) :: myname = 'readSlipInterfaceSpec'
        character(len=max_length_string) :: text
        integer :: nslip
        integer :: ios
        integer :: i, j
        real(kind=custom_real) :: tmp

        call addTrace(errmsg,myname)
        open(lu,file = filename, status = 'old',iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'Could not open file: '//trim(filename),myname)
            call print(errmsg); return
        endif
        text = 'none'
        do while(text /= 'BEGIN')
            read(lu,'(a)') text
        enddo
        read(lu,*) nslip

        allocate (lsi_spec(nslip))

        do i = 1, nslip
            read(1,*) lsi_spec(i)%type,tmp
            select case(trim(lsi_spec(i)%type))
                case ("elastic")
                    lsi_spec(i)%nu = tmp
                case default
                    call add(errmsg,2,"No accepted interface-type entered. Abort...",myname)
                    call print(errmsg); return
            end select
        enddo
        close(lu)
    end subroutine readSlipInterfaceSpecs
end module slipInterfaceMod
