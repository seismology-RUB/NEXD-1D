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
module materials1DMod

    use parameterMod
    use constantsMod
    use errorMessage
    use matrixMod

    implicit none

    interface setup_materials
        module procedure setup_rvzmaterials
        module procedure setup_poromaterials
    end interface setup_materials

    interface dealloc; module procedure deallocRoVsZimpMaterial; end interface
    interface operator (.rho.); module procedure getDensityRoVsZimpMaterial; end interface
    interface operator (.vs.); module procedure getSVelocityRoVsZimpMaterial; end interface
    interface operator (.imp.); module procedure getImpedanceRoVsZimpMaterial; end interface

    type rovszimp_material
        integer :: nmat
        real(kind=custom_real), dimension(:), pointer :: rho
        real(kind=custom_real), dimension(:), pointer :: vs
        real(kind=custom_real), dimension(:), pointer :: zimp
    end type rovszimp_material

    type porous_material
        integer :: nmat
        !parameters to read in
        integer, dimension(:), pointer :: typ
        real(kind=custom_real), dimension(:), pointer :: rhos   !solid density
        real(kind=custom_real), dimension(:), pointer :: lambda !1st Lamé parameter
        real(kind=custom_real), dimension(:), pointer :: my     !2nd Lamé parameter / shear modulus
        real(kind=custom_real), dimension(:), pointer :: phi    !porosity
        real(kind=custom_real), dimension(:), pointer :: kappa  !permeabilty
        real(kind=custom_real), dimension(:), pointer :: b      !Biot coefficient
        real(kind=custom_real), dimension(:), pointer :: invT   !inverse Tortuosity
        real(kind=custom_real), dimension(:), pointer :: invN   !inverse of Biot modulus
        real(kind=custom_real), dimension(:), pointer :: rho1   !fluid 1 density
        real(kind=custom_real), dimension(:), pointer :: S1     !fluid 1 saturation
        real(kind=custom_real), dimension(:), pointer :: K1     !fluid 1 bulk modulus
        real(kind=custom_real), dimension(:), pointer :: ny1    !fluid 1 viscosity
        real(kind=custom_real), dimension(:), pointer :: rho2   !fluid 2 density
        real(kind=custom_real), dimension(:), pointer :: S2     !fluid 2 saturation
        real(kind=custom_real), dimension(:), pointer :: K2     !fluid 2 bulk modulus
        real(kind=custom_real), dimension(:), pointer :: ny2    !fluid 2 viscosity

        real(kind=custom_real), dimension(:), pointer :: fitting_m   !fitting parameter for van Genuchten model
        real(kind=custom_real), dimension(:), pointer :: fitting_n   !fitting parameter for van Genuchten model
        real(kind=custom_real), dimension(:), pointer :: fitting_chi !fitting parameter for van Genuchten model
        real(kind=custom_real), dimension(:), pointer :: Sr1    !fluid 1 residual saturation
        real(kind=custom_real), dimension(:), pointer :: Sr2    !fluid 2 residual saturation

        !parameters to calculate
        real(kind=custom_real), dimension(:), pointer :: krel1  !fluid 1 relative permeability
        real(kind=custom_real), dimension(:), pointer :: krel2  !fluid 2 relative permeability
        real(kind=custom_real), dimension(:), pointer :: dpcdS1 !capillary pressure term
        real(kind=custom_real), dimension(:), pointer :: M      !modulus

        real(kind=custom_real), dimension(:), pointer :: S1eff  !fluid 1 effective saturation
        real(kind=custom_real), dimension(:), pointer :: S2eff  !fluid 2 effective saturation

        real(kind=custom_real), dimension(:,:,:), pointer :: A  !matrix A
        real(kind=custom_real), dimension(:,:,:), pointer :: AP,AM
        real(kind=custom_real), dimension(:,:,:), pointer :: E  !matrix E

        real(kind=custom_real), dimension(:), pointer :: vmax !maximum velocity calculated from matrix A
        real(kind=custom_real), dimension(:), pointer :: vmin !minimum velocity calculated from matrix A
    end type porous_material

    type :: materialVar
        integer :: types                                            !Number of material types
        real(kind=custom_real), dimension(:), allocatable :: b      !Material boundary
        integer, dimension(:), allocatable :: i                     !index_array for the material property
    endtype materialVar

    contains

    subroutine setup_rvzmaterials(par, mat, mesh, rvzmat, errmsg)
        type(parameterVar) :: par
        type(materialVar) :: mat
        type(meshVar) :: mesh
        type(rovszimp_material) :: rvzmat
        type(error_message) :: errmsg

        !Read Material boundaries and create index Matrix for materialproperties
        call materials1D(par, mat, mesh, 1, "input/mat_bap" , errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        !read Materialtypes and respective properties
        call readMaterialProperties(mat, rvzmat, 1, "input/material", errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    end subroutine

    subroutine setup_poromaterials(par, mat, mesh, poromat, errmsg)
        type(parameterVar) :: par
        type(materialVar) :: mat
        type(meshVar) :: mesh
        type(porous_material) :: poromat
        type(error_message) :: errmsg

        !Read Material boundaries and create index Matrix for materialproperties
        call materials1D(par, mat, mesh, 1, "input/mat_bap" , errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        !read Materialtypes and respective properties
        call readPorousMaterialProperties(par, mat, poromat, 1, "input/porousmaterial", errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    end subroutine

    subroutine readMaterialProperties(mat, rvzmat, lu, filename, errmsg)
        !in/output
        type(materialVar) :: mat
        type(rovszimp_material) :: rvzmat
        !input
        type (error_message) :: errmsg
        integer :: lu
        character(len=*) :: filename
        !local
        integer :: i
        integer :: ios
        integer :: mattypes
        character(len=21) :: myname = 'readMaterialProperies'
        character(len=255) :: text
        character(len=10) :: mati

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
        read(lu,*) mattypes
        !check whether material specified in 'mat_bap' is available in 'materials' or not
        do i = 1, size(mat%i)
            if (mat%i(i) > mattypes) then
                write (mati,'(i10)') mat%i(i)
                call add(errmsg,2,"Material "//trim(adjustl(mati))//&
                 " specified in 'mat_bap' is not available in 'materials'.",myname)
                return
            endif
        enddo
        allocate(rvzmat%rho(mattypes))
        allocate(rvzmat%vs(mattypes))
        allocate(rvzmat%zimp(mattypes))
        do i = 1, mattypes
            read(lu,*) rvzmat%rho(i), rvzmat%vs(i)
            rvzmat%zimp(i) = rvzmat%rho(i)*rvzmat%vs(i)
        enddo

        mat%types = mattypes

        close(lu)
    end subroutine readMaterialProperties

    subroutine readPorousMaterialProperties(par, mat, poromat, lu, filename, errmsg)
        !in/output
        type(materialVar) :: mat
        type(porous_material) :: poromat
        !input
        type (parameterVar), intent(in) :: par
        type (error_message) :: errmsg
        integer :: lu
        character(len=*) :: filename
        !local
        integer :: i!,iii,iij
        integer :: ios
        integer :: mattypes
        character(len=27) :: myname = 'readPorousMaterialProperies'
        character(len=255) :: text
        !character(len=3) :: str
        character(len=10) :: mati
        real(kind=custom_real) :: rho,rhoast,M1,M2,Mti,Mti1,Mti2
        real(kind=custom_real) :: lambdaast1_inv,lambdaast2_inv,lambdati11_inv,lambdati12_inv,lambdati21_inv,lambdati22_inv
        !real(kind=custom_real) :: m1,m2,beta11,beta12,beta21,beta22,rho1,rho22,rho21,rho31,rho32

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
        read(lu,*) mattypes

        !check whether material specified in 'mat_bap' is available in 'porousmaterials' or not
        do i = 1, size(mat%i)
            if (mat%i(i) > mattypes) then
                write (mati,'(i10)') mat%i(i)
                call add(errmsg,2,"Material "//trim(adjustl(mati))//&
                 " specified in 'mat_bap' is not available in 'porousmaterials'.",myname)
                return
            endif
        enddo

        allocate(poromat%typ(mattypes))
        allocate(poromat%rhos(mattypes))
        allocate(poromat%lambda(mattypes))
        allocate(poromat%my(mattypes))
        allocate(poromat%phi(mattypes))
        allocate(poromat%kappa(mattypes))
        allocate(poromat%b(mattypes))
        allocate(poromat%invT(mattypes))
        allocate(poromat%invN(mattypes))
        allocate(poromat%rho1(mattypes))
        allocate(poromat%S1(mattypes))
        allocate(poromat%K1(mattypes))
        allocate(poromat%ny1(mattypes))
        allocate(poromat%rho2(mattypes))
        allocate(poromat%S2(mattypes))
        allocate(poromat%K2(mattypes))
        allocate(poromat%ny2(mattypes))
        allocate(poromat%krel1(mattypes))
        allocate(poromat%krel2(mattypes))
        allocate(poromat%dpcdS1(mattypes))
        allocate(poromat%fitting_m(mattypes))
        allocate(poromat%fitting_n(mattypes))
        allocate(poromat%fitting_chi(mattypes))
        allocate(poromat%Sr1(mattypes))
        allocate(poromat%Sr2(mattypes))
        allocate(poromat%S1eff(mattypes))
        allocate(poromat%S2eff(mattypes))
        allocate(poromat%M(mattypes))
        allocate(poromat%A(2+2*par%fluidn,2+2*par%fluidn,mattypes))
        allocate(poromat%AP(2+2*par%fluidn,2+2*par%fluidn,mattypes))
        allocate(poromat%AM(2+2*par%fluidn,2+2*par%fluidn,mattypes))
        allocate(poromat%E(2+2*par%fluidn,2+2*par%fluidn,mattypes))
        allocate(poromat%vmax(mattypes))
        allocate(poromat%vmin(mattypes))

        do i = 1, mattypes
            read(lu,*) poromat%typ(i), poromat%rhos(i), poromat%lambda(i), poromat%my(i), poromat%phi(i), poromat%kappa(i),&
             poromat%b(i), poromat%invT(i), poromat%invN(i), poromat%rho1(i), poromat%S1(i), poromat%K1(i),&
             poromat%ny1(i), poromat%rho2(i), poromat%S2(i), poromat%K2(i), poromat%ny2(i),&
             poromat%fitting_n(i), poromat%fitting_chi(i), poromat%Sr1(i), poromat%Sr2(i)

            if (par%fluidn == 1) then
                poromat%S1(i) = 1.0
                poromat%S2(i) = 0.0
            endif
            if ((poromat%S1(i) - 1.0) < epsilon(poromat%S1(i))) then
                !It will be required to have rho2=rho1 if one fluid
                !vanishes.
                poromat%rho2(i) = poromat%rho1(i)
                !The second fluid must not have a residual saturation.
                poromat%Sr2(i) = 0.
            endif

            if ((poromat%S1(i) - 1.0) < epsilon(poromat%S1(i)) .and. poromat%S2(i) < epsilon(poromat%S2(i))) then
                !Just to avoid a division by zero. All terms with K2 will become zero anyway in this case.
                poromat%K2(i) = eps
            elseif (poromat%S1(i) < epsilon(poromat%S1(i)) .and. (poromat%S1(i) - 1.0) < epsilon(poromat%S1(i))) then
                !Just to avoid a division by zero. All terms with K1 will become zero anyway in this case.
                poromat%K1(i) = eps
            endif

            !check for correct input
            if ((poromat%phi(i) < 0.) .or. (poromat%phi(i) > 1.)) call add(errmsg,2,'phi must be between 0 and 1',myname)
            if (par%fluidn == 2 .and. (poromat%S1(i) + poromat%S2(i) - 1.0) > epsilon(poromat%S1(i))) call add(errmsg,2,'S1+S2 must be 1',myname)
            if ((poromat%b(i) < 0.) .or. (poromat%b(i) > 1.)) call add(errmsg,2,'b must be between 0 and 1',myname)
            if (poromat%phi(i) > epsilon(poromat%phi(i))) then
                if (par%calculate_tortuosity) then
                    if (poromat%invT(i) > 1. .or. poromat%invT(i) < 0.) call add(errmsg,2,'r must be between 0 and 1',myname)
                else
                    if (poromat%invT(i) > 1. .or. poromat%invT(i) < 0.) call add(errmsg,2,'1/T must be between 0 and 1',myname)
                endif
            endif
            if (poromat%S1(i) > 1-poromat%Sr2(i)) call add(errmsg,2,'S1 <= 1-Sr2 required',myname)
            if (poromat%S2(i) > 1-poromat%Sr1(i)) call add(errmsg,2,'S2 <= 1-Sr1 required',myname)
            if (poromat%S1(i) < poromat%Sr1(i)) call add(errmsg,2,'S1 >= Sr1 required',myname)
            if (poromat%S2(i) < poromat%Sr2(i)) call add(errmsg,2,'S2 >= Sr2 required',myname)

            if (poromat%phi(i) < epsilon(poromat%phi(i))) then
                !When there is no porespace, the skeleton and the
                !material building up the skeleton have the same bulk
                !modulus. This means, that b must be 0.
                poromat%b(i) = 0.
            elseif (poromat%typ(i) == 3 .or. poromat%typ(i) == 4 .or. poromat%typ(i) == 8) then
                !here: K^d = poromat%lambda and K_s = poromat%invN
                poromat%b(i) = 1 - poromat%lambda(i) / poromat%invN(i)
            elseif (poromat%typ(i) == 5 .or. poromat%typ(i) == 6 .or. poromat%typ(i) == 9) then
                !here: K^d = poromat%lambda
                poromat%b(i) = (poromat%phi(i)+1)/2 - &
                 sqrt((poromat%phi(i)+1)**2/4 - poromat%phi(i) - poromat%lambda(i)*poromat%invN(i))
            endif

            if (poromat%typ(i) == 3 .or. poromat%typ(i) == 4 .or. poromat%typ(i) == 8) then
                !here: K_s = poromat%invN before calculation
                poromat%invN(i) = (poromat%b(i)-poromat%phi(i))/poromat%invN(i)
                !now 1/N = poromat%invN!
            endif

            if (par%calculate_tortuosity) then
                !here: r = poromat%invT before calculation
                poromat%invT(i) = 1/ (1 - poromat%invT(i) * ( 1 - 1 / poromat%phi(i) ))
            endif

            poromat%krel1(i) = 1.
            poromat%krel2(i) = 1.
            poromat%dpcdS1(i) = 1.
            if (par%fluidn == 2) then
                poromat%S1eff(i) = (poromat%S1(i)-poromat%Sr1(i))/(1-poromat%Sr2(i)-poromat%Sr1(i))
                poromat%S2eff(i) = (poromat%S2(i)-poromat%Sr2(i))/(1-poromat%Sr1(i)-poromat%Sr2(i))
                select case (poromat%typ(i))
                    case (1,3,5) ! van Genuchten model
                        !calculate more properties using van Genuchten (1980) model
                        poromat%fitting_m(i) = 1 - 1/poromat%fitting_n(i)
                        !calculate relative permeability
                        poromat%krel1(i) = sqrt(poromat%S1eff(i))*(1-(1-poromat%S1eff(i)**(1/poromat%fitting_m(i)))**&
                         poromat%fitting_m(i))**2
                        poromat%krel2(i) = sqrt(1-poromat%S1eff(i))*(1-poromat%S1eff(i)**(1/poromat%fitting_m(i)))**&
                         (2*poromat%fitting_m(i))
                        !calculate capillary pressure influence
                        poromat%dpcdS1(i) = - ((poromat%rho1(i)*g)/((poromat%fitting_n(i)-1)*poromat%fitting_chi(i))) * &
                         ((poromat%S1(i)**(-1/poromat%fitting_m(i))-1)**(1/poromat%fitting_n(i)-1)) * &
                         (poromat%S1(i)**(-(1/poromat%fitting_m(i)+1)))
                    case (2,4,6) ! Brooks & Corey model
                        !here p_b = fitting_n and \lambda_{BC} = fitting_chi
                        !calculate relative permeability
                        poromat%krel1(i) = poromat%S1eff(i)**((2+3*poromat%fitting_chi(i))/poromat%fitting_chi(i))
                        poromat%krel2(i) = (1-poromat%S1eff(i))**2 * &
                         (1-poromat%S1eff(i)**((2+poromat%fitting_chi(i))/poromat%fitting_chi(i)))
                        !calculate capillary pressure influence
                        poromat%dpcdS1(i) = - poromat%fitting_n(i)/(poromat%fitting_chi(i)*(poromat%Sr2(i)-poromat%Sr1(i)))*&
                         poromat%S1(i)**(-(1+poromat%fitting_chi(i))/poromat%fitting_chi(i))
                    case (7,8,9) ! Douglas Jr. et al. model
                        !here A = fitting_n
                        !calculate relative permeability
                        poromat%krel1(i) = (1-poromat%S1(i)/(1-poromat%Sr2(i)))**2
                        poromat%krel2(i) = ( (poromat%S1(i) - poromat%Sr1(i)) / (1-poromat%Sr1(i)) )**2
                        !calculate capillary pressure influence
                        poromat%dpcdS1(i) = - 2 * poromat%fitting_n(i) * ((poromat%S1(i)-poromat%Sr1(i))**(-3) + &
                         poromat%Sr2(i)**2 / ((1-poromat%S1(i))**2 * (1-poromat%Sr1(i)-poromat%Sr2(i))**2 ) )
                end select
            endif

            !calculate M
            poromat%M(i) = 1/(poromat%invN(i) + poromat%phi(i)*(poromat%S1(i)/poromat%K1(i)+&
             poromat%S2(i)/poromat%K2(i)) - ( poromat%S1(i)*poromat%S2(i)/(poromat%K1(i)*poromat%K2(i))*poromat%phi(i)&
             + poromat%S1(i)*poromat%S2(i)*poromat%invN(i)*(poromat%S1(i)/poromat%K2(i)+poromat%S2(i)/poromat%K1(i)))*&
             poromat%dpcdS1(i))

            if (poromat%typ(i) == 3 .or. poromat%typ(i) == 4 .or. poromat%typ(i) == 5 .or. poromat%typ(i) == 6 &
             .or. poromat%typ(i) == 8 .or. poromat%typ(i) == 9) then
                !here: on the rhs K^d = poromat%lambda
                poromat%lambda(i) = poromat%lambda(i) - 2./3. * poromat%my(i) + poromat%M(i)*(poromat%b(i)**2) * &
                 (1-poromat%S1(i)*poromat%S2(i)*(poromat%S1(i)/poromat%K2(i)+poromat%S2(i)/poromat%K1(i))*&
                 poromat%dpcdS1(i))
            endif

            rho    = (1-poromat%phi(i))*poromat%rhos(i) + poromat%phi(i)*poromat%S1(i)*poromat%rho1(i) &
                                                        + poromat%phi(i)*poromat%S2(i)*poromat%rho2(i)
            rhoast = (1-poromat%phi(i))*poromat%rhos(i) + poromat%phi(i)*poromat%S1(i)*poromat%rho1(i) &
                                                        + poromat%phi(i)*poromat%S2(i)*poromat%rho2(i)
            if (poromat%phi(i) > epsilon(poromat%phi(i))) then
                rhoast = rhoast - poromat%invT(i)*poromat%phi(i)*poromat%S1(i)*poromat%rho1(i) &
                                - poromat%invT(i)*poromat%phi(i)*poromat%S2(i)*poromat%rho2(i)
            endif

            lambdaast1_inv = (poromat%phi(i) * poromat%S1(i) * poromat%ny1(i)) / (poromat%kappa(i) * poromat%krel1(i)) !/C1
            lambdaast2_inv = (poromat%phi(i) * poromat%S2(i) * poromat%ny2(i)) / (poromat%kappa(i) * poromat%krel2(i)) !/C2
            lambdati11_inv = (poromat%phi(i) * poromat%S1(i) * (poromat%invT(i)**2-poromat%invT(i)) * lambdaast1_inv)/rhoast
            lambdati12_inv = (poromat%phi(i) * poromat%S2(i) * (poromat%invT(i)**2-poromat%invT(i)) * lambdaast2_inv)/rhoast
            lambdati21_inv = lambdaast1_inv * poromat%invT(i) / poromat%rho1(i) + lambdati11_inv
            lambdati22_inv = lambdaast2_inv * poromat%invT(i) / poromat%rho2(i) + lambdati12_inv

            M1   = poromat%M(i)*(1-poromat%S1(i)*poromat%S2(i)/poromat%K1(i)*poromat%dpcdS1(i))
            M2   = poromat%M(i)*(1-poromat%S1(i)*poromat%S2(i)/poromat%K2(i)*poromat%dpcdS1(i))
            Mti  = poromat%M(i)*(1+poromat%S1(i)*poromat%S2(i)*poromat%invN(i)/poromat%phi(i)*poromat%dpcdS1(i))
            Mti1 = poromat%M(i)* &
             (1-(poromat%S1(i)**2 * poromat%invN(i)/poromat%phi(i) + poromat%S1(i)/poromat%K1(i)) * poromat%dpcdS1(i))
            Mti2 = poromat%M(i)* &
             (1-(poromat%S2(i)**2 * poromat%invN(i)/poromat%phi(i) + poromat%S2(i)/poromat%K2(i)) * poromat%dpcdS1(i))

            poromat%A(:,:,i) =  0.
            poromat%A(1,2,i) = -(poromat%lambda(i) + 2*poromat%my(i))
            poromat%A(2,1,i) = -1/rhoast
            poromat%E(:,:,i) =  0.

            if (poromat%phi(i) > epsilon(poromat%phi(i))) then
                poromat%A(1,2,i) = poromat%A(1,2,i) &
                    + poromat%b(i)*poromat%phi(i)*(poromat%S1(i)*M2 + poromat%S2(i)*M1)
                poromat%A(1,4,i) = -poromat%phi(i)*poromat%S1(i)*M2*poromat%b(i)
                poromat%A(2,3,i) = -poromat%phi(i)*poromat%S1(i)*poromat%invT(i)/rhoast
                poromat%A(3,2,i) =  M2*poromat%b(i) - poromat%phi(i)*(poromat%S1(i)*Mti2 + poromat%S2(i)*Mti)
                poromat%A(3,4,i) =  poromat%phi(i)*poromat%S1(i)*Mti2
                poromat%A(4,1,i) =  (poromat%invT(i)-1)/rhoast
                poromat%A(4,3,i) =  poromat%invT(i)/poromat%rho1(i) &
                    + poromat%phi(i)*poromat%S1(i)*(poromat%invT(i)**2-poromat%invT(i))/rhoast
                if (par%fluidn == 2) then
                    poromat%A(1,6,i) = -poromat%phi(i)*poromat%S2(i)*M1*poromat%b(i)
                    poromat%A(2,5,i) = -poromat%phi(i)*poromat%S2(i)*poromat%invT(i)/rhoast
                    poromat%A(3,6,i) =  poromat%phi(i)*poromat%S2(i)*Mti
                    poromat%A(4,5,i) =  poromat%phi(i)*poromat%S2(i)*(poromat%invT(i)**2-poromat%invT(i))/rhoast
                    poromat%A(5,2,i) =  M1*poromat%b(i) - poromat%phi(i)*(poromat%S1(i)*Mti + poromat%S2(i)*Mti1)
                    poromat%A(5,4,i) =  poromat%phi(i)*poromat%S1(i)*Mti
                    poromat%A(5,6,i) =  poromat%phi(i)*poromat%S2(i)*Mti1
                    poromat%A(6,1,i) =  (poromat%invT(i)-1)/rhoast
                    poromat%A(6,3,i) =  poromat%phi(i)*poromat%S1(i)*(poromat%invT(i)**2-poromat%invT(i))/rhoast
                    poromat%A(6,5,i) =  poromat%invT(i)/poromat%rho2(i) &
                        + poromat%phi(i)*poromat%S2(i)*(poromat%invT(i)**2-poromat%invT(i))/rhoast
                endif
                poromat%E(2,4,i) =  (poromat%phi(i)*poromat%S1(i)*poromat%invT(i)*lambdaast1_inv) / rhoast
                poromat%E(4,4,i) =  - lambdati21_inv
                if (par%fluidn == 2) then
                    poromat%E(2,6,i) = (poromat%phi(i)*poromat%S2(i)*poromat%invT(i)*lambdaast2_inv) / rhoast
                    poromat%E(4,6,i) = - lambdati12_inv
                    poromat%E(6,4,i) = - lambdati11_inv
                    poromat%E(6,6,i) = - lambdati22_inv
                    poromat%E(2,2,i) = -poromat%E(2,4,i)-poromat%E(2,6,i)
                    poromat%E(4,2,i) = -poromat%E(4,4,i)-poromat%E(4,6,i)
                    poromat%E(6,2,i) = -poromat%E(6,4,i)-poromat%E(6,6,i)
                else
                    poromat%E(2,2,i) = -poromat%E(2,4,i)
                    poromat%E(4,2,i) = -poromat%E(4,4,i)
                endif
            else
                poromat%A(4,1,i) = -1/rhoast
                if (par%fluidn == 2) then
                    poromat%A(6,1,i) = -1/rhoast
                endif
            endif

            !write (*,*) "A",i
            !do iii = 1,2+2*par%fluidn
            !    write (*,"(100g15.5)") (poromat%A(iii,iij,i),iij=1,2+2*par%fluidn)
            !enddo
            !write (*,*) "E",i
            !do iii = 1,2+2*par%fluidn
            !    write (*,"(100g15.5)") (poromat%E(iii,iij,i),iij=1,2+2*par%fluidn)
            !enddo

            call calcAPAM(poromat%A(:,:,i),poromat%AP(:,:,i),poromat%AM(:,:,i),poromat%vmax(i),poromat%vmin(i))
        enddo
        mat%types = mattypes
        close(lu)

        if (par%fluidn == 1) then
            deallocate(poromat%S1)
            deallocate(poromat%rho2)
            deallocate(poromat%S2)
            deallocate(poromat%K2)
            deallocate(poromat%ny2)
            deallocate(poromat%krel1)
            deallocate(poromat%krel2)
            deallocate(poromat%dpcdS1)
            deallocate(poromat%fitting_m)
            deallocate(poromat%fitting_n)
            deallocate(poromat%fitting_chi)
            deallocate(poromat%Sr1)
            deallocate(poromat%Sr2)
            deallocate(poromat%S1eff)
            deallocate(poromat%S2eff)
        endif
    end subroutine readPorousMaterialProperties

    subroutine calcAPAM(A,AP,AM,vmax,vmin)
        !input
        real(kind=custom_real), dimension(:,:), intent(in) :: A
        !output
        real(kind=custom_real), dimension(:,:), intent(out) :: AP, AM
        real(kind=custom_real), intent(out) :: vmax, vmin
        !local variables
        type(error_message) :: errmsg
        character(len=8) :: myname = 'calcAPAM'
        character(len=100) :: errstr
        real(kind=custom_real), dimension(:,:), allocatable :: Awork
        real(kind=custom_real), dimension(:), allocatable :: WR, WI    !contain real and imaginary part of eigenvalues, respectively
        real(kind=custom_real), dimension(:,:), allocatable :: VR      !matrix containing the eigenvectors as columns
        real(kind=custom_real), dimension(:,:), allocatable :: LambdaP, LambdaM  !matrix containing the eigenvalues on the diagonal elements
        real(kind=custom_real), dimension(:,:), allocatable :: VL,invVR
        real(kind=custom_real), dimension(:), allocatable :: work
        integer :: N, lwork, info, i
        integer, dimension(:), allocatable :: ipiv
        integer, dimension(2) :: shapeA

        shapeA=shape(A)
        N=shapeA(1)
        allocate(Awork(N,N))
        allocate(ipiv(N))
        allocate(WR(N))
        allocate(WI(N))
        allocate(VR(N,N))
        allocate(LambdaP(N,N))
        allocate(LambdaM(N,N))
        allocate(invVR(N,N))
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

        LambdaP = 0.
        LambdaM = 0.
        do i=1,N
            if (WR(i) > 0) then
                LambdaP(i,i)=WR(i)
            else
                LambdaM(i,i)=WR(i)
            endif
        enddo

        invVR = VR
        if (custom_real==size_double) then
            call dgetrf(N, N, invVR, N, ipiv, info)
            !do workspace query
            allocate(work(1)); lwork = -1
            call dgetri(N, invVR, N, ipiv, work, lwork, info)
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            call dgetri(N, invVR, N, ipiv, work, lwork, info)
            deallocate(work)
        elseif (custom_real==size_real) then
            call sgetrf(N, N, invVR, N, ipiv, info)
            !do workspace query
            allocate(work(1)); lwork = -1
            call sgetri(N, invVR, N, ipiv, work, lwork, info)
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            call sgetri(N, invVR, N, ipiv, work, lwork, info)
            deallocate(work)
        else
            call add(errmsg,2,'CUSTOM_REAL is neither SIZE_REAL nor SIZE_DOUBLE',myname)
            stop
        endif

        AP = matmul(VR,matmul(LambdaP,invVR))
        AM = matmul(VR,matmul(LambdaM,invVR))
        vmax = maxval(WR)
        vmin = vmax
        do i=1,size(WR)
            if (WR(i) > 0) vmin = min(vmin,WR(i))
        enddo

        !print *,WR

        if (allocated(Awork)) deallocate(Awork)
        if (allocated(ipiv)) deallocate(ipiv)
        if (allocated(WR)) deallocate(WR)
        if (allocated(WI)) deallocate(WI)
        if (allocated(VR)) deallocate(VR)
        if (allocated(LambdaP)) deallocate(LambdaP)
        if (allocated(LambdaM)) deallocate(LambdaM)
        if (allocated(invVR)) deallocate(invVR)
    end subroutine calcAPAM

    subroutine materials1D(par, mat, mesh, lu, filename, errmsg)
        !in/outut
        type(materialVar), intent(inout) :: mat
        type(meshVar), intent(inout) :: mesh
        !input
        type (parameterVar) :: par
        type (error_message) :: errmsg
        integer :: lu
        character(len=*) :: filename
        !locals
        integer :: i, j, l
        integer :: ios
        integer, dimension(par%matn) :: locOnGrid
        character(len=11) :: myname = 'materials1D'
        character(len=255) :: text

        call addTrace(errmsg,myname)
        open(lu,file = filename, status = 'old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            return
        endif

        text ='none'
        do while(text /= 'BEGIN')
            read(lu,'(a)') text
        enddo

        allocate(mat%b(par%matn))
        allocate(mat%i(par%matn))
        do i = 1, par%matn
            read(lu,*) mat%b(i), mat%i(i)
            if (par%physical_coordinates) then
                if (mat%b(i) > mesh%xmax) then
                    call add(errmsg,2,'Boundary not on grid. Abort...',myname)
                    return
                endif
                locOnGrid(i) = nint(((mat%b(i)+abs(mesh%xmin))*mesh%ncell/(mesh%xmax+abs(mesh%xmin)))+mesh%nghost)
                if (i == par%matn) locOnGrid(i) = nint(((mat%b(i)+abs(mesh%xmin))&
                *mesh%ncell/(mesh%xmax+abs(mesh%xmin)))+2*mesh%nghost)
            else
                locOnGrid(i) = nint(mat%b(i)) + mesh%nghost
                if (i == par%matn) locOnGrid(i) = nint(mat%b(i)) + 2*mesh%nghost
            endif
        enddo
        if (mat%b(par%matn) < mesh%xmax) then
            call add(errmsg,2,'Materials are not specified for the whole grid. Abort...',myname)
            return
        endif

        allocate(mesh%ntom(mesh%Np,mesh%K))
        allocate(mesh%etom(mesh%K))
        !Set up indexmatrix for the material properties
        l = 1
        do i = 1, par%matn
            do j = l, locOnGrid(i)
                mesh%ntom(:,j) = mat%i(i)
                mesh%etom(j) = mat%i(i)
            enddo
            l = locOnGrid(i)+1
        enddo
        close(lu)
    end subroutine

    subroutine deallocMatArrays(mat)
        type(materialVar) :: mat
        if (allocated(mat%b)) deallocate(mat%b)
        if (allocated(mat%i)) deallocate(mat%i)
    end subroutine

    subroutine deallocRoVsZimpMaterial(this)
        !Function to deallocate material object
        type (rovszimp_material) :: this

        if (associated(this%rho)) deallocate(this%rho)
        if (associated(this%vs)) deallocate(this%vs)
        if (associated(this%zimp)) deallocate(this%zimp)
    end subroutine deallocRoVsZimpMaterial

    subroutine deallocPorousMaterial(this)
        !Function to deallocate porous material object
        type (porous_material) :: this

        if (associated(this%rhos)) deallocate(this%rhos)
        if (associated(this%lambda)) deallocate(this%lambda)
        if (associated(this%my)) deallocate(this%my)
        if (associated(this%phi)) deallocate(this%phi)
        if (associated(this%kappa)) deallocate(this%kappa)
        if (associated(this%b)) deallocate(this%b)
        if (associated(this%invT)) deallocate(this%invT)
        if (associated(this%invN)) deallocate(this%invN)
        if (associated(this%rho1)) deallocate(this%rho1)
        if (associated(this%S1)) deallocate(this%S1)
        if (associated(this%K1)) deallocate(this%K1)
        if (associated(this%ny1)) deallocate(this%ny1)
        if (associated(this%rho2)) deallocate(this%rho2)
        if (associated(this%S2)) deallocate(this%S2)
        if (associated(this%K2)) deallocate(this%K2)
        if (associated(this%ny2)) deallocate(this%ny2)
        if (associated(this%krel1)) deallocate(this%krel1)
        if (associated(this%krel2)) deallocate(this%krel2)
        if (associated(this%dpcdS1)) deallocate(this%dpcdS1)
        if (associated(this%M)) deallocate(this%M)
        if (associated(this%A)) deallocate(this%A)
        if (associated(this%AP)) deallocate(this%AP)
        if (associated(this%AM)) deallocate(this%AM)
        if (associated(this%E)) deallocate(this%E)
        if (associated(this%vmax)) deallocate(this%vmax)
        if (associated(this%vmin)) deallocate(this%vmin)
        if (associated(this%Sr1)) deallocate(this%Sr1)
        if (associated(this%Sr2)) deallocate(this%Sr2)
        if (associated(this%S1eff)) deallocate(this%S1eff)
        if (associated(this%S2eff)) deallocate(this%S2eff)
    end subroutine deallocPorousMaterial
!-----------------------------------------------------------------
!> \brief return selected density value
!
    function getDensityRoVsZimpMaterial(this,j) result(res)
    type (rovszimp_material), intent(in) :: this
    integer, intent(in) :: j
    real(kind=custom_real) :: res
    res = this%rho(j)
    end function getDensityRoVsZimpMaterial
!-----------------------------------------------------------------
!> \brief return selected vs value
!
    function getSVelocityRoVsZimpMaterial(this,j) result(res)
    type (rovszimp_material), intent(in) :: this
    integer, intent(in) :: j
    real(kind=custom_real) :: res
    res = this%vs(j)
    end function getSVelocityRoVsZimpMaterial
!-----------------------------------------------------------------
!> \brief return selected impedance value
!
    function getImpedanceRoVsZimpMaterial(this,j) result(res)
    type (rovszimp_material), intent(in) :: this
    integer, intent(in) :: j
    real(kind=custom_real) :: res
    res = this%zimp(j)
    end function getImpedanceRoVsZimpMaterial
!-----------------------------------------------------------------
!> \brief return selected Matrix A
!
    function getMatrixAElasticMaterial(this,j) result(res)
        type (rovszimp_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(2,2) :: res
        res = 0.
        res(1,2) = -this%vs(j)*this%zimp(j)
        res(2,1) = -this%vs(j)/this%zimp(j)
    end function getMatrixAElasticMaterial

    function getMatrixAPElasticMaterial(this,j) result(res)
        type (rovszimp_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(2,2) :: res
        res(1,1) = this%vs(j)
        res(1,2) = -this%vs(j)*this%zimp(j)
        res(2,1) = -this%vs(j)/this%zimp(j)
        res(2,2) = this%vs(j)
        res = res/2.
    end function getMatrixAPElasticMaterial

    function getMatrixAMElasticMaterial(this,j) result(res)
        type (rovszimp_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(2,2) :: res
        res(1,1) = -this%vs(j)
        res(1,2) = -this%vs(j)*this%zimp(j)
        res(2,1) = -this%vs(j)/this%zimp(j)
        res(2,2) = -this%vs(j)
        res = res/2.
    end function getMatrixAMElasticMaterial

    function getMatrixAPorousMaterial(this,j) result(res)
        type (porous_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res
        if (associated(this%S2)) then
            allocate(res(6,6))
        else
            allocate(res(4,4))
        endif
        res = this%A(:,:,j)
    end function getMatrixAPorousMaterial

    function getMatrixAPPorousMaterial(this,j) result(res)
        type (porous_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res
        if (associated(this%S2)) then
            allocate(res(6,6))
        else
            allocate(res(4,4))
        endif
        res = this%AP(:,:,j)
    end function getMatrixAPPorousMaterial

    function getMatrixAMPorousMaterial(this,j) result(res)
        type (porous_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res
        if (associated(this%S2)) then
            allocate(res(6,6))
        else
            allocate(res(4,4))
        endif
        res = this%AM(:,:,j)
    end function getMatrixAMPorousMaterial

!-----------------------------------------------------------------
!> \brief return selected Matrix E
!
    function getMatrixEPorousMaterial(this,j) result(res)
        type (porous_material), intent(in) :: this
        integer, intent(in) :: j
        real, dimension(:,:), allocatable :: res
        if (associated(this%S2)) then
            allocate(res(6,6))
        else
            allocate(res(4,4))
        endif
        res = this%E(:,:,j)
    end function getMatrixEPorousMaterial
!------------------------------------------------------------------
!> \brief return maximum velocity
!
    function maximumVelocityRoVsZimpMaterial(this) result(res)
        type (rovszimp_material), intent(in) :: this
        real(kind=custom_real) :: res
        res = maxval(this%vs)
    end function maximumVelocityRoVsZimpMaterial
end module
