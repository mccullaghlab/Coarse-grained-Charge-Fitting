
! USAGE analytic_charge_min.x -cfg [CONFIG FILE NAME]
!
! CONFIG FILE FORMAT:
! 
!atompsffile = [atomic PSF file name]
!atomdcdfile = [atomic DCD file name]
!cgdcdfile = [CG DCD file name] 
!outfile = [output file - will contain positions and charges]
!thresh = [minimization convergence criteria]
!minsteps = [number of minimization steps]
!
!
! OUTPUT: Minimization information will be written to standard out.  The CG positions and 
! charges will be written in XYZ format to the outfile
!
! REQUIREMENTS: lapack
!
! CITATION:
! P. McCullagh, P. T. Lake, M. McCullagh. Deriving Coarse-grained Charges from All-atom Systems: An Analytic Solution. J. Chem. Theory Comp., 2016.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Modules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module atomData

        integer nAtoms
        real (kind=8), allocatable :: atomPos(:,:)
        real (kind=8), allocatable :: atomCharges(:)
        character*80 atomDcdFile
        character*80 atomPsfFile
        character*80 outFile
        real (kind=8), allocatable :: atomEspMat(:,:)

endmodule atomData

module cgData

        integer nCgs
        real (kind=8), allocatable :: cgPos(:,:)
        character*80 cgDcdFile
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
        real (kind=8) lambda
        real (kind=8), save :: delta = 0.1

endmodule minData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program esp_grid
        use atomData
        use cgData
        implicit none
        character*80 cfgFile
        real (kind=8) omp_get_wtime
        real (kind=8) ti,tf

        ti = omp_get_wtime()
        ! read config file name from command line using option -cfg [config file name]
        call parse_command_line(cfgFile)
        
        ! read the rest of the config parameters from the config file
        call parse_config_file(cfgFile)

        ! read the trajectories and perform the fits
        call read_dcd_fit_esp()

        tf = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') tf-ti


        write(*,'("Please cite: P. McCullagh, P. T. Lake, M. McCullagh. Deriving Coarse-grained Charges from All-atom Systems: An Analytic Solution. J. Chem. Theory Comp., 2016.")')

endprogram esp_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate arrays
subroutine allocate_arrays()
        use atomData
        use minData
        use cgData
        use integralData
        implicit none

        ! integral based matrices
        allocate(A(nCgs,nAtoms),B(nCgs,nCgs),C(nAtoms,nAtoms),cgCharges(nCgs))
        A = 0.0
        B = 0.0
        C = 0.0

        ! min arrays
        allocate(chargeForce(nCgs,3),chargeGradient(nCgs))

endsubroutine allocate_arrays

! routine to perform minimization.  Will first read coordinates and charges of atom system and initial CG coordinates
subroutine read_dcd_fit_esp()
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
        ! read atom DCD header
        call read_dcd_header(atomDcdFile,nAtoms,nSteps,20)
        ! read CG DCD header
        call read_dcd_header(cgDcdFile,nCgs,nSteps,30)
        ! allocate CG position array
        allocate(cgPos(nCgs,3))
        ! read coordinates from DCD files
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
        lambda = 0.01
        previousRss = intRss
        do minStep = 1, nMinSteps


                call compute_charge_gradient(atomCharges,cgCharges,nAtoms,nCgs,intRss,atomPos,cgPos)
        
                call move_cg_sites(cgPos,nCgs,chargeForce)
                call recompute_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B)

                ! recompute charges and residual
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
                        print*, "The maximum gradient has converged."
                        exit
                endif 
                if (abs(intRss-previousRss) < 0.001) then 
                        lambda = 0.001
                endif
                if (abs(intRss-previousRss) < 0.00001) then 
                        lambda = 0.00001
                endif
                if (abs(intRss-previousRss) < thresh) then 
                        print*, "The residual has converged."
                        exit
                endif
                previousRss = intRss

        enddo
        close(35)

endsubroutine read_dcd_fit_esp

! move the CG sites along the spacial gradient
subroutine move_cg_sites(cgPos,nCgs,chargeForce)
        use minData, only : lambda
        implicit none
        integer nCgs
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) chargeForce(nCgs,3)

        ! move cg position vector
        cgPos = cgPos + lambda*chargeForce

endsubroutine move_cg_sites

! Compute the spacial gradient of the charge residual
subroutine compute_charge_gradient(atomCharges,cgCharges,nAtoms,nCgs,intRss,atomPos,cgPos)
        use minData
        implicit none
        real (kind=8), parameter :: pi = 3.1415926535
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) cgCharges(nCgs)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) intRss
        real (kind=8) plusRss
        real (kind=8) minusRss
        real (kind=8) temp(3)
        real (kind=8) gradient
        integer cg1
        integer cg2
        integer atom
        integer k
      
        chargeForce = 0.0
        do cg1 = 1, nCgs
                do atom = 1, nAtoms
                        temp = cgPos(cg1,:) - atomPos(atom,:)
                        temp = -4*pi*cgCharges(cg1)*atomCharges(atom)* temp/norm2(temp)
                        chargeForce(cg1,:) = chargeForce(cg1,:) + temp
                enddo
                do cg2 = 1, nCgs
                        if (cg1 .ne. cg2) then
                                temp = cgPos(cg1,:) - cgPos(cg2,:)
                                temp = 4*pi*cgCharges(cg1)*cgCharges(cg2)* temp/norm2(temp)
                                chargeForce(cg1,:) = chargeForce(cg1,:) + temp
                        endif
                enddo
                chargeGradient(cg1) = norm2(chargeForce(cg1,:))
        enddo 
        
endsubroutine compute_charge_gradient

! read config file name from command line
subroutine parse_command_line(cfgFile)
        implicit none
        integer i
        character*30 arg
        character*80 cfgFile
        logical cfgFileFlag

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
                case default 
                        print '(a,a,/)', 'Unrecognized command-line option: ', arg 
                        print*, 'Usage: analytic_charge_min.x -cfg [cfg file]'
                        stop 
                end select 
                i = i+1
                if (i.ge.command_argument_count()) exit
        enddo

        if (cfgFileFlag.eqv..false.) then
                write(*,'("Must provide a cfg file using command line argument -cfg [cfg file name]")')
                stop
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
        logical atomDcdFlag
        logical atomPsfFlag
        logical cgDcdFlag
        logical outFileFlag
        logical threshFlag
        logical nMinStepsFlag

        atomDcdFlag = .false.
        atomPsfFlag = .false.
        cgDcdFlag = .false.
        outFileFlag = .false.
        nMinStepsFlag = .false.
        threshFlag = .false.

        open(12,file=cfgFile)
        do while(ios>=0)
                read(12,'(a600)',IOSTAT=ios) line
                call split(line,'=',firstWord, sep)
                if (line .ne. "") then
                        if (firstWord .eq. "atomdcdfile") then
                                atomDcdFile = line
                                write(*,*) "Atom DCD file:", atomDcdFile
                                atomDcdFlag = .true.
                        else if (firstWord .eq. "atompsffile") then
                                atomPsfFile = line
                                write(*,*) "Atom PSF file:", atomPsfFile
                                atomPsfFlag = .true.
                        else if (firstWord .eq. "cgdcdfile") then
                                cgDcdFile = line
                                write(*,*) "CG DCD file:", cgDcdFile
                                cgDcdFlag = .true.
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


        if (atomDcdFlag.eqv..false.) then
                write(*,'("Must provide a atom DCD file using config file")')
                stop
        endif
        if (atomPsfFlag.eqv..false.) then
                write(*,'("Must provide a atom PSF file using config file")')
                stop
        endif
        if (cgDcdFlag.eqv..false.) then
                write(*,'("Must provide a CG DCD file using config file")')
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
        real (kind=8), parameter :: pi = 3.1415926535
        integer nAtoms
        integer nCg
        real (kind=8) A(nCg,nAtoms)
        real (kind=8) ATemp(nCg,nAtoms)
        real (kind=8) B(nCg,nCg)
        real (kind=8) BTemp(nCg,nCg)
        real (kind=8) D(nCg-1,nCg)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) cgCharges(nCg,1)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) atomChargesM(nAtoms,1)
        real (kind=8) solution(nCg)
        real (kind=8) temp(1,1)
        real (kind=8) rss
        real (kind=8) rss2
        real (kind=8) rss3
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
        solution = matmul(ATemp,atomCharges)
        ! determine the solution to system of linear equations ATemp*X = newB
        call dgesv(nCg,1, BTemp,nCg,ipiv,solution,nCg,info)
        
        cgCharges(:,1) = real(solution(1:nCg))
        atomChargesM(:,1) = atomCharges

        ! compute residual sum of squares
        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))+matmul(transpose(cgCharges),matmul(B,cgCharges))-2*matmul(transpose(cgCharges),matmul(A,atomChargesM))
        rss = 2*pi*temp(1,1)
        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))-matmul(transpose(cgCharges),matmul(B,cgCharges))
        rss2 = 2*pi*temp(1,1)
        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))-matmul(transpose(cgCharges),matmul(A,atomChargesM))
        rss3 = 2*pi*temp(1,1)
        print*, rss, rss3, rss2

endsubroutine integral_fit_charges

!requires lapack
subroutine compute_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
        implicit none
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nAtoms)
        real (kind=8) B(nCgs,nCgs)
        real (kind=8) C(nAtoms,nAtoms)
        real (kind=8) dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2

        ! populate B matrix with negative distances between CG sites
        do cgSite1 = 1, nCgs-1
                do cgSite2 = cgSite1+1,nCgs
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-cgPos(cgSite2,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        B(cgSite1,cgSite2) = B(cgSite1,cgSite2)-dist
                        !symmetrize the matrix
                        B(cgSite2,cgSite1) = B(cgSite1,cgSite2)
                enddo
        enddo

        ! populate A matrix with negative distance between CG sites and atoms
        do cgSite1 = 1, nCgs
                do atom1 = 1,nAtoms
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-atomPos(atom1,j)
                                dist = dist + temp*temp
                        enddo
                        A(cgSite1,atom1) = A(cgSite1,atom1)-sqrt(dist)
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
subroutine recompute_A_B_matrices(atomPos,nAtoms,cgPos,nCgs,A,B)
        implicit none
        integer nAtoms
        integer nCgs
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) A(nCgs,nAtoms)
        real (kind=8) B(nCgs,nCgs)
        real (kind=8) dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1

        ! populate B matrix with negative distances between CG sites
        do cgSite1 = 1, nCgs-1 
                do cgSite2 = cgSite1+1,nCgs
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-cgPos(cgSite2,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        B(cgSite1,cgSite2) = -dist
                        !symmetrize the matrix
                        B(cgSite2,cgSite1) = B(cgSite1,cgSite2)
                 enddo
        enddo

        ! populate A matrix with negative distance between CG sites and atoms
        do cgSite1 = 1, nCgs
                do atom1 = 1,nAtoms
                        dist = 0
                        do j=1,3
                                temp = cgPos(cgSite1,j)-atomPos(atom1,j)
                                dist = dist + temp*temp
                        enddo
                        dist = sqrt(dist)
                        A(cgSite1,atom1) = -dist
                enddo
        enddo

endsubroutine recompute_A_B_matrices

!Read atomic charges from psf file
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

