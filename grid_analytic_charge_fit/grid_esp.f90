
!requires lapack
! CITATION:
! P. McCullagh, P. T. Lake, M. McCullagh. Deriving Coarse-grained Charges from All-atom Systems: An Analytic Solution. J. Chem. Theory Comp., 2016.

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
        real (kind=8), allocatable :: cgCharges(:)
        character*80 cgXyzFile
        real (kind=8), allocatable :: cgEspMat(:,:)

endmodule cgData

module gridData

        integer nGrids(3)
        integer totalGrids
        real (kind=8) gridMin(3)
        real (kind=8) gridMax(3)
        real (kind=8) deltaGrid
        real (kind=8) cutoff2
        integer gridCut
        integer maxCutInt
        real (kind=8), save :: minThresh2 = 0.01

endmodule gridData

module integralData

        real (kind=8), allocatable :: A(:,:)
        real (kind=8), allocatable :: B(:,:)
        real (kind=8), allocatable :: C(:,:)
        real (kind=8), allocatable :: intCgCharges(:)
        real (kind=8) gridRss
        real (kind=8) intRss

endmodule integralData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program esp_grid
        use atomData
        use cgData
        use inputData
        use gridData, only : totalGrids
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

        tf = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') tf-ti

endprogram esp_grid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate arrays
subroutine allocate_arrays()
        use atomData
        use gridData
        use cgData
        use integralData
        implicit none

        ! grid based array containing the ESP multiplied by the atomic charges
        allocate(atomEspMat(nCgs,nAtoms))
        ! grid based array containing 1/dist between grids and CG sites
        allocate(cgEspMat(nCgs,nCgs))
        ! zero out each array as they will be added to
        atomEspMat = 0.0
        cgEspMat = 0.0

        ! integral based matrices
        allocate(A(nCgs,nCgs),B(nCgs,nAtoms),C(nAtoms,nAtoms),intCgCharges(nCgs))

        A = 0.0
        B = 0.0
        C = 0.0

endsubroutine allocate_arrays

subroutine read_xyz_fit_esp()
        use inputData
        use atomData
        use cgData
        use gridData
        use integralData
        implicit none
        integer atom1
        integer k
        integer nSteps
        integer step
        real (kind=8) centerVec(3)
        real (kind=8), allocatable :: lagrangeCgCharges(:)
        real (kind=8) gridCenter(3)
        real (kind=8) temp
        real (kind=8) tempLagrange
        integer filenum
        integer totalSteps
        integer cutoff_int
        real (kind=8) cutoff
        real (kind=8) atomMin(3)
        real (kind=8) atomMax(3)
        real (kind=8) atomAvg

        gridCenter = (gridMax+gridMin) / 2.0
        print*, "Grid center:", gridCenter(1),gridCenter(2),gridCenter(3)
                
        ! read the atom psf file to obtain number of atoms and charges.  atomPos and charge array are allocated in this routine
        call read_psf_file()
        call read_dcd_header(atomXyzFile,nAtoms,nSteps,20)
        call read_dcd_header(cgXyzFile,nCgs,nSteps,30)
        allocate(cgPos(nCgs,3))
        allocate(cgCharges(nCgs))
        call read_dcd_step(atomPos,nAtoms,20)
        call read_dcd_step(cgPos,nCgs,30)
!        call read_atom_xyz()
!        call read_cg_xyz()
        allocate(lagrangeCgCharges(nCgs))
        ! allocate some grid and integral arrays
        call allocate_arrays()
        ! place molecule at origin
        do k=1,3
                atomAvg = sum(atomPos(:,k))/ dble(nAtoms)
                atomPos(:,k) = atomPos(:,k) - atomAvg
                cgPos(:,k) = cgPos(:,k) - atomAvg
                atomMin(k) = minval(atomPos(:,k))
                atomMax(k) = maxval(atomPos(:,k))
        enddo  
        ! update A, B and C matrices for integral method
        call update_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
        ! fit using integral approach
        call integral_fit_charges(A, B, C, atomCharges, intCgCharges, nAtoms, nCgs, intRss)
        open(77,file = "ssc_v_cutoff.dat")
        write(77,'(a20,a20,a20,a25,a20)') "cutoff", "Sum Grid Charges", "Grid SSD", "Sum Lagrange Charges", "Lagrange Grid SSD"
        do cutoff_int = 1, maxCutInt
                cutoff = cutoff_int*1.0
                cutoff2 = cutoff*cutoff
                gridMin = atomMin - cutoff - deltaGrid
                gridMax = atomMax + cutoff + deltaGrid
                nGrids =  int( (gridMax-gridMin) / deltaGrid ) + 1
                totalGrids = nGrids(1)*nGrids(2)*nGrids(3)
                write(*,'("Iteration ", i5, " cutoff value:", f10.5, " nGrids", i10,i10,i10)') cutoff_int, cutoff, nGrids(1),nGrids(2),nGrids(3)
                atomEspMat = 0.0
                cgEspMat = 0.0
                ! ??
                call compute_grid_esp_mats(atomPos,nAtoms,cgPos,nCgs,atomEspMat,cgEspMat)
                call fit_cg_charges_lagrange(atomEspMat,nAtoms,cgEspMat,nCgs,atomCharges,lagrangeCgCharges)
                call fit_cg_charges(atomEspMat,nAtoms,cgEspMat,nCgs,atomCharges,cgCharges)
                temp = 0.0
                tempLagrange = 0.0
                do k = 1, nCgs
                        temp = temp + (cgCharges(k)-intCgCharges(k))**2
                        tempLagrange = tempLagrange + (lagrangeCgCharges(k)-intCgCharges(k))**2
                enddo
                write(77,'(f20.10,f20.10,f20.10,f25.10,f20.10)') cutoff,sum(cgCharges), temp, sum(lagrangeCgCharges),tempLagrange
        enddo
        close(77)

        ! print CG charges
        open(35,file=outFile)
        write(35,'(a10,a20,a20,a20,a20,a20)') "CG Site", "Grid Charges", "Integral Charges", "Squared difference", "Lagrange Charges", "SD"
        write(*,'(a10,a20,a20,a20,a20,a20)') "CG Site", "Grid Charges", "Integral Charges", "Squared difference", "Lagrange Charges", "SD"
        temp = 0
        do k=1, nCgs
                temp = temp + (cgCharges(k)-intCgCharges(k))**2
                write(35,'(i10,f20.10,f20.10,f20.10,f20.10,f20.10)') k , cgCharges(k), intCgCharges(k), (cgCharges(k)-intCgCharges(k))**2, lagrangeCgCharges(k), (lagrangeCgCharges(k)-intCgCharges(k))**2
                write(*,'(i10,f20.10,f20.10,f20.10,f20.10,f20.10)') k , cgCharges(k), intCgCharges(k), (cgCharges(k)-intCgCharges(k))**2, lagrangeCgCharges(k), (lagrangeCgCharges(k)-intCgCharges(k))**2
        enddo    
        write(35,'("Sums:", a20,a20,a20,a20)') "Atomic", "Grid Charges", "Integral Charges", "Sum of sq diff" 
        write(35,'(4x,f20.10,f20.10,f20.10,f20.10)') sum(atomCharges), sum(cgCharges), sum(intCgCharges), temp
        close(35)
        write(*,'("Sums:", a20,a20,a20,a20)') "Atomic", "Grid Charges", "Integral Charges", "Sum of sq diff"
        write(*,'(4x,f20.10,f20.10,f20.10,f20.10)') sum(atomCharges), sum(cgCharges), sum(intCgCharges), temp
        close(35)
!        write(*,'("CG grid charge:", f10.5, "CG int charges:", f10.5, " total atom charge:", f10.5)') sum(cgCharges), sum(intCgCharges), sum(atomCharges)

endsubroutine read_xyz_fit_esp


! solve set of linear equations to solve for CG charges 
subroutine fit_cg_charges(atomEspMat,nAtoms,cgEspMat,nCgs,atomCharges,cgCharges)
        implicit none
        integer nCgs
        integer nAtoms
        real (kind=8) atomEspMat(nCgs,nAtoms)
        real (kind=8) cgEspMat(nCgs,nCgs)
        real (kind=8) cgCharges(nCgs)
        real (kind=8) atomCharges(nAtoms,1)
        integer info
        integer ipiv(nCgs)
        real(kind=8) B(nCgs,1)
        integer cg1, cg2

!        ! first symmetrize cgEspMat matrix
!        do cg1 = 1, nCgs-1
!                do cg2 = cg1+1,nCgs
!                        cgEspMat(cg2,cg1) = cgEspMat(cg1,cg2)
!                enddo
!        enddo

        ! now compute right hand side of equation as A*q
        B = matmul(atomEspMat,atomCharges)

        call dgesv(nCgs, 1, cgEspMat, nCgs, ipiv, B, nCgs, info)

        cgCharges = B(1:nCgs,1)

endsubroutine fit_cg_charges

! solve set of linear equations to solve for CG charges with Lagrange multipliers to enforce charge conservation
subroutine fit_cg_charges_lagrange(atomEspMat,nAtoms,cgEspMat,nCgs,atomCharges,cgCharges)
        implicit none
        integer nCgs
        integer nAtoms
        real (kind=8) atomEspMat(nCgs,nAtoms)
        real (kind=8) cgEspMat(nCgs,nCgs)
        real (kind=8) cgCharges(nCgs)
        real (kind=8) atomCharges(nAtoms,1)
        integer info
        integer ipiv(nCgs+1)
        real(kind=8) B(nCgs+1,1)
        real(kind=8) A(nCgs+1,nCgs+1)
        integer cg1, cg2

        ! first symmetric cgEspMat matrix
        do cg1 = 1, nCgs-1
                do cg2 = cg1+1,nCgs
                        cgEspMat(cg2,cg1) = cgEspMat(cg1,cg2)
                enddo
        enddo

        ! now compute right hand side of equation as A*q
        B(1:nCgs,:) = matmul(atomEspMat,atomCharges)
        B(nCgs+1,1) = sum(atomCharges)

        A(1:nCgs,1:nCgs) = cgEspMat
        A(nCgs+1,1:nCgs) = 1.0
        A(1:nCgs,nCgs+1) = 1.0
        A(nCgs+1,nCgs+1) = 0.0

        call dgesv(nCgs+1, 1, A, nCgs+1, ipiv, B, nCgs+1, info)

        cgCharges = B(1:nCgs,1)

endsubroutine fit_cg_charges_lagrange

!
subroutine compute_grid_esp_mats(atomPos,nAtoms,cgPos,nCgs,atomEspMat,cgEspMat)
        use inputData, only : nThreads
        use gridData
        implicit none
        integer nAtoms
        real (kind=8) atomPos(nAtoms,3)
        integer nCgs
        real (kind=8) cgPos(nCgs,3)
        real (kind=8) atomEspMat(nCgs,nAtoms)
        real (kind=8) cgEspMat(nCgs,nCgs)
        real (kind=8) atomEspMat_temp(nCgs,nAtoms)
        real (kind=8) cgEspMat_temp(nCgs,nCgs)
        integer atom
        integer cg1, cg2
        real (kind=8) diff1(3), diff2(3)
        real (kind=8) dist1_2, dist2_2
        real (kind=8) dist1, dist2
        integer grid, gridCount(3)
        real (kind=8) gridPos(3)
        logical gridPass 

        !$omp parallel num_threads(nThreads) private(gridPass,atomEspMat_temp,cgEspMat_temp,grid,gridCount,gridPos,diff1,diff2,dist1_2,dist1,dist2_2,dist2,cg1,cg2,atom) shared(nGrids,deltaGrid,gridMin,cutoff2,cgEspMat,atomEspMat,nAtoms,nCgs)
        !$omp do
        do grid = 1, totalGrids
        
                gridCount(1) = int(grid/dble(nGrids(3)*nGrids(2))) 
                gridCount(2) = int( (grid - gridCount(1)*nGrids(3)*nGrids(2)) / dble(nGrids(3)))
                gridCount(3) = mod(grid-1,nGrids(3))

                gridPos = (gridCount+0.5)*deltaGrid + gridMin
                
                atomEspMat_temp = atomEspMat
                cgEspMat_temp = cgEspMat
                
                gridPass = .true.

                do cg1 = 1, nCgs

                        diff1 = cgPos(cg1,:) - gridPos
                        dist1_2 = dot_product(diff1,diff1)
                        if (dist1_2 < cutoff2 .and. dist1_2 > minThresh2) then
                                dist1 = sqrt(dist1_2)
                                do cg2 = cg1, nCgs
                                        diff2 = cgPos(cg2,:) - gridPos
                                        dist2_2 = dot_product(diff2,diff2)
                                        if (dist2_2 < cutoff2 .and. dist2_2 > minThresh2) then
                                                dist2 = sqrt(dist2_2)
                                                !$omp atomic
                                                cgEspMat_temp(cg1,cg2) = cgEspMat_temp(cg1,cg2) + 1/(dist1*dist2)
                                        elseif (dist2_2 < minThresh2) then
                                                gridPass = .false.
                                                exit
                                        endif

                                enddo
                
                                if (gridPass.eqv..true.) then
                                        do atom=1, nAtoms
                                                diff2 = atomPos(atom,:) - gridPos
                                                dist2_2 = dot_product(diff2,diff2)
                                                if (dist2_2 < cutoff2 .and. dist2_2 > minThresh2) then
                                                        dist2 = sqrt(dist2_2)
                                                        !$omp atomic
                                                        atomEspMat_temp(cg1,atom) = atomEspMat_temp(cg1,atom) + 1/(dist1*dist2)
                                                elseif (dist2_2 < minThresh2) then
                                                        gridPass = .false.
                                                        exit
                                                endif
                                        enddo
                                endif
                        elseif (dist1_2 < minThresh2) then
                                gridPass = .false.
                                exit
                        endif
                        if (gridPass .eqv. .false.) then
                                exit
                        endif

                enddo

                if (gridPass .eqv. .true.) then
                        atomEspMat = atomEspMat_temp
                        cgEspMat = cgEspMat_temp
                endif

        enddo
        !$omp end do
        !$omp end parallel

endsubroutine compute_grid_esp_mats


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
                        print*, 'Usage: esp_grid.x -cfg [cfg file] -stride [delta step size] -np [number of threads]'
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
        use gridData
        use atomData
        use cgData
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
        logical gridMaxFlag
        logical gridMinFlag
        logical deltaGridFlag
        logical cutoffFlag
        logical maxCutFlag

        atomXyzFlag = .false.
        atomPsfFlag = .false.
        cgXyzFlag = .false.
        outFileFlag = .false.
        gridMaxFlag = .false.
        gridMinFlag = .false.
        deltaGridFlag = .false.
        cutoffFlag = .false.
        maxCutFlag = .false.

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
                        else if (firstWord .eq. "deltagrid") then
                                read(line,'(f10.5)') deltaGrid
                                write(*,*) "deltaGrid:", deltaGrid
                                deltaGridFlag = .true.
                        else if (firstWord .eq. "cutoff") then
                                read(line,'(f10.5)') cutoff
                                cutoff2 = cutoff * cutoff
                                write(*,*) "cutoff:", cutoff
                                cutoffFlag = .true.
                        else if (firstWord .eq. "maxcut") then
                                read(line,'(i10)') maxCutInt
                                write(*,*) "maxCut:", maxCutInt
                                maxCutFlag = .true.
                        else if (firstWord .eq. "maxgrid") then
                                do k = 1,3
                                        call split(line,' ',firstWord, sep)
                                        read(firstWord,'(f10.5)') gridMax(k)
                                enddo
                                write(*,*) "max grid:", gridMax(1), gridMax(2), gridMax(3)
                                gridMaxFlag = .true.
                        else if (firstWord .eq. "mingrid") then
                                do k = 1,3
                                        call split(line,' ',firstWord, sep)
                                        read(firstWord,'(f10.5)') gridMin(k)
                                enddo
                                write(*,*) "min grid:", gridMin(1),gridMin(2),gridMin(3)
                                gridMinFlag = .true.
                        endif
                endif
        enddo

        close(12)


        if (cutoffFlag.eqv..false.) then
                write(*,'("cutoff not defined.  using default of 20.0.  change with cutoff = [cutoff] in config file")')
                cutoff = 20.0
                cutoff2 = cutoff * cutoff
        endif
        if (maxCutFlag.eqv..false.) then
                write(*,'("maxCut not defined.  using default of 50.0.  change with maxcut = [maxcut] in config file")')
                maxCutInt = 50
        endif
        if (atomXyzFlag.eqv..false.) then
                write(*,'("Must provide a atom xyz file using command line argument -adcd [atom dcd file name]")')
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
                write(*,'("Must provide an output file name  using command line argument -o [output file name]")')
                stop
        endif
        if (gridMaxFlag.eqv..false.) then
                write(*,'("Using default max radius of 400. Change this with command line argument -max [max radius]")')
                gridMax = 100.0
        endif
        if (gridMinFlag.eqv..false.) then
                write(*,'("Using default min radius of 0. Change this with command line argument -min [min radius]")')
                gridMin = 0.0
        endif
        if (deltaGridFlag.eqv..false.) then
                write(*,'("Using default of deltaGrid=1.0.")')
                deltaGrid=1.0
        endif
        ! compute total number of grid points
        totalGrids = 1.0
        do k=1,3
                nGrids(k) = int( (gridMax(k)-gridMin(k))/deltaGrid) + 1
                totalGrids = totalGrids*nGrids(k)
        enddo
        gridCut = int(cutoff/deltaGrid) + 1

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
subroutine update_A_B_C_matrices(atomPos,nAtoms,cgPos,nCgs,A,B,C)
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

endsubroutine update_A_B_C_matrices


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
        allocate(cgCharges(nCgs))
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

