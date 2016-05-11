
!requires lapack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Modules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module inputData

        integer deltaStep
        integer nThreads

endmodule inputData

module atomData

        integer nAtoms
        real (kind=8), allocatable :: atomPos(:,:)
        real (kind=8), allocatable :: atomCharges(:)
        character*80 atomXyzFile
        character*80 atomPsfFile
        character*80 outFile
        real (kind=8), allocatable :: atomEspMat(:,:)

endmodule atomData

module cgData

        integer nCgs
        real (kind=8), allocatable :: cgPos(:,:)
        character*80 cgXyzFile
        real (kind=8), allocatable :: cgEspMat(:,:)

endmodule cgData


module integralData

        real (kind=8), allocatable :: A(:,:)
        real (kind=8), allocatable :: B(:,:)
        real (kind=8), allocatable :: C(:,:)
        real (kind=8), allocatable :: cgCharges(:)
        real (kind=8) gridRss
        real (kind=8) intRss

endmodule integralData

module minData

        integer nMinSteps
        real (kind=8) thresh
        real (kind=8), allocatable :: chargeForce(:,:)
        real (kind=8), allocatable :: chargeGradient(:)
        real (kind=8) maxGradient
!        real (kind=8), save :: lamda = 0.01
        real (kind=8) lambda
        real (kind=8), save :: delta = 0.1

endmodule minData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program esp_grid
        use atomData
        use cgData
        use inputData
        implicit none
        character*80 cfgFile
        real (kind=8) omp_get_wtime
        real (kind=8) ti,tf

        ti = omp_get_wtime()
        ! read config file, number of openMP threads and stride from command line
        call parse_command_line(cfgFile,nThreads,deltaStep)
        
        ! read the rest of the config parameters from the config file
        call parse_config_file(cfgFile)

        ! read the trajectories and perform the fits
        call read_xyz_fit_esp()

        ! compute pair coulomb energies
        call compute_dimer_elec_energy()

        tf = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') tf-ti

endprogram esp_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_dimer_elec_energy
        use atomData
        use cgData
        use integralData
        implicit none
        real (kind=8), parameter :: chargeConvert = 332.063697438
        integer atom1, atom2
        real (kind=8) diff(3)
        real (kind=8) dist2
        real (kind=8) coulombEnergyAtom
        real (kind=8) coulombEnergyCg

        ! compute coulomb energy
        coulombEnergyAtom = 0.0
        do atom1 = 1, nAtoms/2
                do atom2 = nAtoms/2 + 1, nAtoms
                        diff = atomPos(atom1,:) - atomPos(atom2,:)
                        dist2 = dot_product(diff,diff)
                        coulombEnergyAtom = coulombEnergyAtom + atomCharges(atom1)*atomCharges(atom2)/sqrt(dist2)
                enddo
        enddo
        ! compute coulomb energy
        coulombEnergyCg = 0.0
        do atom1 = 1, nCgs/2
                do atom2 = nCgs/2 + 1, nCgs
                        diff = cgPos(atom1,:) - cgPos(atom2,:)
                        dist2 = dot_product(diff,diff)
                        coulombEnergyCg = coulombEnergyCg + cgCharges(atom1)*cgCharges(atom2)/sqrt(dist2)
                enddo
        enddo

        open(72,file="atom_cg_elec_energy.dat")
        write(72,'(2f30.10)') coulombEnergyAtom*chargeConvert, coulombEnergyCg*chargeConvert
        close(72)

endsubroutine compute_dimer_elec_energy


! allocate arrays
subroutine allocate_arrays()
        use atomData
        use minData
        use cgData
        use integralData
        implicit none

        ! integral based matrices
        allocate(A(nCgs,nCgs),B(nCgs,nAtoms),C(nAtoms,nAtoms),cgCharges(nCgs))
        A = 0.0
        B = 0.0
        C = 0.0

        ! min arrays
        allocate(chargeForce(nCgs,3),chargeGradient(nCgs))

endsubroutine allocate_arrays

subroutine read_xyz_fit_esp()
        use inputData
        use atomData
        use cgData
        use minData
        use integralData
        implicit none
        integer atom1
        integer k
        integer nSteps
        integer step
        integer filenum
        integer minStep
        real (kind=8) atomAvg
        real (kind=8) residual
        real (kind=8) previousRss

                
        ! read the atom psf file to obtain number of atoms and charges.  atomPos and charge array are allocated in this routine
        call read_psf_file()
        call read_dcd_header(atomXyzFile,nAtoms,nSteps,20)
        call read_dcd_header(cgXyzFile,nCgs,nSteps,30)
        allocate(cgPos(nCgs,3))
        call read_dcd_step(atomPos,nAtoms,20)
        call read_dcd_step(cgPos,nCgs,30)
        ! allocate some integral arrays
        call allocate_arrays()
        ! place molecule at origin
        do k=1,3
                atomAvg = sum(atomPos(:,k))/ dble(nAtoms)
                atomPos(:,k) = atomPos(:,k) - atomAvg
                cgPos(:,k) = cgPos(:,k) - atomAvg
        enddo  
        ! update A, B and C matrices for integral method
        call compute_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
        ! fit using integral approach
        call integral_fit_charges(A, B, C, atomCharges,cgCharges, nAtoms, nCgs, intRss)

        ! perform minimzation
        open(35,file=outFile)
        lambda = 0.1
        previousRss = intRss
        do minStep = 1, nMinSteps


                call compute_charge_gradient(A,B,C,atomCharges,cgCharges,nAtoms,nCgs,intRss,atomPos,cgPos)
        
                call move_cg_sites(cgPos,nCgs,chargeForce)

                call recompute_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B)
                call integral_fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCgs, intRss)
                write(*,'("Step:", i10, " maximum gradient: ", f20.10, " residual:", f20.10)') minStep, maxval(chargeGradient), intRss
                ! print CG charges
                write(35,'(i10)') nCgs
                write(35,'(a10,a20,a20,a20,a20)') "CG Site", "X", "Y", "Z", "Integral Charges"
                do k=1, nCgs
                        write(35,'("C ",f20.10,f20.10,f20.10,f20.10)')  cgPos(k,1), cgPos(k,2), cgPos(k,3), cgCharges(k)
                enddo  
                flush(35)  

                if (maxval(chargeGradient) < thresh) then
                        exit
                endif       
                if (abs(intRss-previousRss) < 0.2) then 
                        lambda = 0.01
                endif
                if (abs(intRss-previousRss) < 0.01) then 
                        print*, "The residual has converged."
                        exit
                endif
                previousRss = intRss

        enddo
        close(35)

endsubroutine read_xyz_fit_esp

subroutine move_cg_sites(cgPos,nCgs,chargeForce)
        use minData, only : lambda
        implicit none
        integer nCgs
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) chargeForce(nCgs,3)

        ! move cg position vector
        cgPos = cgPos - lambda*chargeForce

endsubroutine move_cg_sites

! we will compute the gradient of moving the charges using numerical differentiation
subroutine compute_charge_gradient(A,B,C,atomCharges,cgCharges,nAtoms,nCgs,intRss,atomPos,cgPos)
        use inputData
        use minData
        implicit none
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nCgs)
        real (kind=8) B(nCgs,nAtoms)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) cgCharges(nCgs,1)
        real (kind=8) atomCharges(nAtoms,1)
        real (kind=8) intRss
        real (kind=8) plusRss
        real (kind=8) minusRss
        real (kind=8) temp
        real (kind=8) gradient
        integer cg1
        integer cg2
        integer k
      
        do cg1 = 1, nCgs
                gradient = 0
                do k=1,3   
                        cgPos(cg1,k) = cgPos(cg1,k) + delta
                        ! update A and B matrices for new CG pos
                        call update_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,cg1)
                        ! fit using integral approach
                        call integral_fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCgs, plusRss)
                        ! now shift negative
                        cgPos(cg1,k) = cgPos(cg1,k) - 2.0*delta
                        ! update A and B matrices for new CG pos
                        call update_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,cg1)
                        ! fit using integral approach
                        call integral_fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCgs, minusRss)
                        chargeForce(cg1,k) = (plusRss-minusRss) / (2.0*delta)
                        gradient = gradient + chargeForce(cg1,k)*chargeForce(cg1,k)
                        ! move it back
                        cgPos(cg1,k) = cgPos(cg1,k) + delta
                        call update_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,cg1)
                enddo
                chargeGradient(cg1) = sqrt(gradient)

        enddo 
        
endsubroutine compute_charge_gradient

subroutine parse_command_line(cfgFile,nThreads,deltaStep)
        implicit none
        integer i
        character*30 arg
        character*80 cfgFile
        integer nThreads
        integer deltaStep
        logical deltaStepFlag
        logical cfgFileFlag
        logical nThreadsFlag

        deltaStepFlag = .false.
        nThreadsFlag = .false.
        cfgFileFlag = .false.

        i=1
        do 
   
                call get_command_argument(i, arg) 
   
                select case (arg) 
   
                case ('-cfg')
                        i = i+1
                        call get_command_argument(i,cfgFile)
                        cfgFileFlag=.true.
                        print*, "CFG File: ", cfgFile
                case ('-stride') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') deltaStep
                        print*, "Stride: ", deltaStep
                        deltaStepFlag = .true.
                case ('-np') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') nThreads
                        print*, "Number of threads: ", nThreads
                        nThreadsFlag = .true.
                case default 
                        print '(a,a,/)', 'Unrecognized command-line option: ', arg 
                        print*, 'Usage: analytic_charge_min.x -cfg [cfg file] -stride [delta step size] -np [number of threads]'
                        stop 
                end select 
                i = i+1
                if (i.ge.command_argument_count()) exit
        enddo

        if (cfgFileFlag.eqv..false.) then
                write(*,'("Must provide a cfg file using command line argument -cfg [cfg file name]")')
                stop
        endif
        if (deltaStepFlag.eqv..false.) then
                write(*,'("Using default step size of 1.  Change this with command line argument -stride [delta step size]")')
                deltaStep = 1
        endif
        if (nThreadsFlag.eqv..false.) then
                write(*,'("Using default of 1 thread.  Change this with command line argument -np [number of threads]")')
                nThreads = 1
        endif

endsubroutine parse_command_line

!read config information from file
subroutine parse_config_file(cfgFile)
        use atomData
        use cgData
        use minData
        implicit none
        character*80 cfgFile
        character*200 line
        character*80,firstWord
        integer i, k
        character*30 arg
        character*30 sep
        real (kind=8) cutoff
        integer ios
        logical atomXyzFlag
        logical atomPsfFlag
        logical cgXyzFlag
        logical outFileFlag
        logical threshFlag
        logical nMinStepsFlag

        atomXyzFlag = .false.
        atomPsfFlag = .false.
        cgXyzFlag = .false.
        outFileFlag = .false.
        nMinStepsFlag = .false.
        threshFlag = .false.

        open(12,file=cfgFile)
        do while(ios>=0)
                read(12,'(a600)',IOSTAT=ios) line
                call split(line,'=',firstWord, sep)
                if (line .ne. "") then
                        if (firstWord .eq. "atomxyzfile") then
                                atomXyzFile = line
                                write(*,*) "Atom xyz file:", atomXyzFile
                                atomXyzFlag = .true.
                        else if (firstWord .eq. "atompsffile") then
                                atomPsfFile = line
                                write(*,*) "Atom PSF file:", atomPsfFile
                                atomPsfFlag = .true.
                        else if (firstWord .eq. "cgxyzfile") then
                                cgXyzFile = line
                                write(*,*) "CG xyz file:", cgXyzFile
                                cgXyzFlag = .true.
                        else if (firstWord .eq. "outfile") then
                                outFile = line
                                write(*,*) "out file:", outFile
                                outFileFlag = .true.
                        else if (firstWord .eq. "minsteps") then
                                read(line,'(i10)') nMinSteps
                                write(*,*) "Number of minimization steps:", nMinSteps
                                nMinStepsFlag = .true.
                        else if (firstWord .eq. "thresh") then
                                read(line,*) thresh
                                write(*,*) "Max gradient threshhold:", thresh
                                threshFlag = .true.
                        endif
                endif
        enddo

        close(12)


        if (atomXyzFlag.eqv..false.) then
                write(*,'("Must provide a atom xyz file using config file")')
                stop
        endif
        if (atomPsfFlag.eqv..false.) then
                write(*,'("Must provide a atom PSF file using config file")')
                stop
        endif
        if (cgXyzFlag.eqv..false.) then
                write(*,'("Must provide a CG xyz file using config file")')
                stop
        endif
        if (outFileFlag.eqv..false.) then
                write(*,'("Must provide an output file name  using config file")')
                stop
        endif
        if (threshFlag.eqv..false.) then
                write(*,'("Using default max gradient of 1e-4")')
                thresh = 1e-4
        endif
        if (nMinStepsFlag.eqv..false.) then
                write(*,'("Using default number of minimzation steps: 1000")')
                nMinSteps=1000
        endif

endsubroutine parse_config_file


!subroutine to compute CG charges from input matrices A and B.  The residual sum of squares is calculated with aid of all-atom matrix C
subroutine integral_fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg, rss)
        implicit none
        integer nAtoms
        integer nCg
        real (kind=8) A(nCg,nCg)
        real (kind=8) ATemp(nCg,nCg)
        real (kind=8) D(nCg-1,nCg)
        real (kind=8) B(nCg,nAtoms)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) BTemp(nCg,nAtoms)
        real (kind=8) cgCharges(nCg,1)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) atomChargesM(nAtoms,1)
        real (kind=8) newB(nCg)
        real (kind=8) temp(1,1)
        real (kind=8) rss
        integer j, i, k
        !lapack routine variables
        integer ipiv(nCg)
        integer info

        !First we need to modify A and B to have the correct matrix properties 
        !create D matrices
        D=0.0
        do i=1,nCg-1
                D(i,i)=1.0
                D(i,i+1)=-1.0
        enddo 
        !multiply A by D0 and B by D1 giving new matrices A and B the correct behavior
        ATemp(1:(nCg-1),:) = matmul(D,A)
        BTemp(1:(nCg-1),:) = matmul(D,B)
        !generate new matrices with last line having 1.0s forcing charge conservation
        ATemp(nCg,:) = 1.0
        BTemp(nCg,:) = 1.0

        ! multiple right hand side of equation by atomic charges
        newB = matmul(BTemp,atomCharges)
        ! determine the solution to system of linear equations ATemp*X = newB
        call dgesv(nCg,1, ATemp,nCg,ipiv,newB,nCg,info)
        
        cgCharges(:,1) = real(newB(1:nCg))
        atomChargesM(:,1) = atomCharges

        ! compute residual sum of squares
        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))+matmul(transpose(cgCharges),matmul(A,cgCharges))-2*matmul(transpose(cgCharges),matmul(B,atomChargesM))
        rss = temp(1,1)

endsubroutine integral_fit_charges

!requires lapack
subroutine compute_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
        implicit none
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nCgs)
        real (kind=8) B(nCgs,nAtoms)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2

        ! populate A matrix with negative distances between CG sites
        do cgSite1 = 1, nCgs-1
                do cgSite2 = cgSite1+1,nCgs
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-cgPos(cgSite2,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        A(cgSite1,cgSite2) = A(cgSite1,cgSite2)-dist
                        !symmetrize the matrix
                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                enddo
        enddo

        ! populate B matrix with negative distance between CG sites and atoms
        do cgSite1 = 1, nCgs
                do atom1 = 1,nAtoms
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-atomPos(atom1,j)
                                dist = dist + temp*temp
                        enddo
                        B(cgSite1,atom1) = B(cgSite1,atom1)-sqrt(dist)
                enddo
        enddo

        ! populate C matrix with negative distance between atoms
        do atom1 = 1, nAtoms-1
                do atom2 = atom1+1, nAtoms
                        dist = 0
                        do j=1,3
                                temp = atomPos(atom1,j)-atomPos(atom2,j)
                                dist = dist + temp*temp
                        enddo
                        C(atom1,atom2) = C(atom1,atom2)-sqrt(dist)
                        ! symmetrize the matrix
                        C(atom2,atom1) = C(atom1,atom2)
                enddo
        enddo

endsubroutine compute_A_B_C_matrices

!
subroutine update_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,cgSite)
        implicit none
        integer nAtoms
        integer nCgs
        integer cgSite
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nCgs)
        real (kind=8) B(nCgs,nAtoms)
        real (kind=8) dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1

        ! populate A matrix with negative distances between CG sites
        cgSite1 = cgSite
        do cgSite2 = 1,nCgs
                if (cgSite2 .ne. cgSite1) then
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-cgPos(cgSite2,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        A(cgSite1,cgSite2) = -dist
                        !symmetrize the matrix
                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                endif
        enddo

        ! populate B matrix with negative distance between CG sites and atoms
        cgSite1 = cgSite
        do atom1 = 1,nAtoms
                dist = 0
                do j=1,3
                        temp = cgPos(cgSite1,j)-atomPos(atom1,j)
                        dist = dist + temp*temp
                enddo
                dist = sqrt(dist)
                B(cgSite1,atom1) = -dist
        enddo

endsubroutine update_A_B_matrices

!
subroutine recompute_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B)
        implicit none
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nCgs)
        real (kind=8) B(nCgs,nAtoms)
        real (kind=8) dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1

        ! populate A matrix with negative distances between CG sites
        do cgSite1 = 1, nCgs-1 
                do cgSite2 = cgSite1+1,nCgs
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-cgPos(cgSite2,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        A(cgSite1,cgSite2) = -dist
                        !symmetrize the matrix
                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                 enddo
        enddo

        ! populate B matrix with negative distance between CG sites and atoms
        do cgSite1 = 1, nCgs
                do atom1 = 1,nAtoms
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-atomPos(atom1,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        B(cgSite1,atom1) = -dist
                enddo
        enddo

endsubroutine recompute_A_B_matrices

!Read atomic charges from some file
subroutine read_atom_xyz
        use atomData
        implicit none
!        real (kind=8), parameter :: chargeConvert = 553.43949573 ! to convert to kT/e with dielectric of 1
        real (kind=8), parameter :: chargeConvert = 1.0 ! 
        integer atom, j
        character*600 line
        character*80,firstWord
        character*30 sep
        integer ios

        !open the xyz file
        open(20,file=atomXyzFile)
        ! first line contains number of atoms 
        read(20,'(a600)') line
        call split(line,' ',firstWord, sep)
        read(firstWord,'(i10)') nAtoms
        print*, "Number of atoms:", nAtoms
        ! allocate position and charge arrays
        allocate(atomCharges(nAtoms))
        allocate(atomPos(nAtoms,3))
        ! skip second line
        read(20,*)
        ! read position and charge of all atoms
        do atom=1, nAtoms
                read(20,'(a600)') line
                call split(line,' ',firstWord, sep)
                do j = 1,3
                        call split(line,' ',firstWord, sep)
                        read(firstWord,'(f20.5)') atomPos(atom,j)
                enddo
                call split(line,' ',firstWord, sep)
                read(firstWord,'(f20.5)') atomCharges(atom)
        enddo
        close(20)

endsubroutine read_atom_xyz


!Read CG positions
subroutine read_cg_xyz
        use cgData
        implicit none
        integer atom, j
        character*600 line
        character*80,firstWord
        character*30 sep
        integer ios

        !open the xyz file
        open(30,file=cgXyzFile)
        ! first line contains number of atoms 
        read(30,'(a600)') line
        call split(line,' ',firstWord, sep)
        read(firstWord,'(i10)') nCgs
        print*, "Number of CG sites:", nCgs
        ! allocate position and charge arrays
        allocate(cgPos(nCgs,3))
        ! skip second line
        read(30,*)
        ! read position and charge of all atoms
        do atom=1, nCgs
                read(30,'(a600)') line
                call split(line,' ',firstWord, sep)
                do j = 1,3
                        call split(line,' ',firstWord, sep)
                        read(firstWord,'(f20.5)') cgPos(atom,j)
                enddo
        enddo
        close(30)

endsubroutine read_cg_xyz


!Read atomic charges from some file
subroutine read_psf_file
        use atomData
        implicit none
!        real (kind=8), parameter :: chargeConvert = 553.43949573 ! to convert to kT/e with dielectric of 1
        real (kind=8), parameter :: chargeConvert = 1.0 ! 
        integer atom, j
        character*6 check           !character to check if NATOM is in the line
        character*8 numChar         !character to read number of atoms.  must be converted to an integer
        character*4 atomCheck
        character*24 posChar
        integer ios

        !open the psf file
        print*, "psf file name here:", atomPsfFile
        open(10,file=atomPsfFile)

        !run through the header of the file and look for the number of atoms
        do 

                read(10,'(a8,2x,a6)') numChar, check

                !if we read the number of atoms exit this do loop
                if (check.eq.'NATOM ') then
                        !this converts the character natoms_char to the integer natoms
                        read(numChar,*) nAtoms
                        write(*,*) "Number of atoms=", nAtoms
                        !Now that we know the number of atoms we must allocate the arrays
                        allocate(atomCharges(nAtoms))
                        allocate(atomPos(nAtoms,3))
                        !Now we loop through the number of atoms and read the pertinent information
                        do atom=1,nAtoms

                                read(10,'(34x,f10.6)') atomCharges(atom) 
                                atomCharges(atom) = atomCharges(atom) * chargeConvert
                                !add to total charge

                        enddo
                        write(*,*) "Total charge in atom psf=", sum(atomCharges)
                elseif (check.eq.'NBOND:') then
                        exit
                endif

        enddo

        close(10)

endsubroutine read_psf_file

