
module endianess

        logical nativeEndianess

endmodule endianess


!Open dcd binary trajectory file and copy the header to the new traj file
subroutine read_dcd_header(inputTraj, nAtomTotal, nSteps,filePointer)
   use endianess
   implicit none
   character*80 inputTraj
   integer nAtomTotal
   integer junk, i, nSteps
   integer filePointer
   character*4 hdr 
   integer startTimeStep
   integer delta_step
   integer lastTimeStep
   integer nFreeAtoms
   real    delta_t
   integer nTitle
   character*80 title
   character turnaround*4,tmp*1
   integer int_swap_endianess
   real real_swap_endianess

   open(filePointer,file=inputTraj,access='stream')
   !Read the initial integer (should be 84)
   read(filePointer) junk
   if (junk == 84) then
           nativeEndianess=.true.
           print*, 'DCD file has native endianess'
   else
           nativeEndianess=.false.
           print*, 'DCD file has not native endianess'
   endif
   if (nativeEndianess.eqv..true.) then
         !Read the header (should be 'CORD' char*4)
         read(filePointer) hdr
         !Read the number of steps (integer)
         read(filePointer) nSteps
         print*, "Number of steps in trajectory file:", nSteps
         !Read the first time step (integer)
         read(filePointer) startTimeStep
         !Read the steps between trajectory writes
         read(filePointer) delta_step
         !Read the last time step (integer)
         read(filePointer) lastTimeStep
         !skip 4 blank integers
         do i=1,4
            read(filePointer) junk
         enddo
         !read the number of free atoms (integer)
         read(filePointer) nFreeAtoms
         !read the time between steps (real single precision)
         read(filePointer) delta_t
         !skip some blank or uniteresting integers
         do i=1,12
            read(filePointer) junk
         enddo
         !read the number of lines for the title (each line is 80 chars long)
         read(filePointer) nTitle
         !read the title lines
         do i=1,nTitle
            read(filePointer) title
         enddo
         !read two more integers
         read(filePointer) junk
         read(filePointer) junk
         !read the number of atoms
         read(filePointer) nAtomTotal
         print*, 'Number of atoms in trajectory file:', nAtomTotal
         !read another integer
         read(filePointer) junk
   else 
         !Read the header (should be 'CORD' char*4)
         read(filePointer) hdr
         !Read the number of steps (integer)
         read(filePointer) nSteps
         nSteps = int_swap_endianess(nSteps)
         print*, "Number of steps in trajectory file:", nSteps
         !Read the first time step (integer)
         read(filePointer) startTimeStep
         startTimeStep = int_swap_endianess(StartTimeStep)
         !Read the steps between trajectory writes
         read(filePointer) delta_step
         delta_step = int_swap_endianess(delta_step)
         !Read the last time step (integer)
         read(filePointer) lastTimeStep
         lastTimeStep = int_swap_endianess(lastTimeStep)
         !skip 4 blank integers
         do i=1,4
            read(filePointer) junk
            junk = int_swap_endianess(junk)
         enddo
         !read the number of free atoms (integer)
         read(filePointer) nFreeAtoms
         nFreeAtoms = int_swap_endianess(nFreeAtoms)
         !read the time between steps (real single precision)
         read(filePointer) delta_t
         delta_t = real_swap_endianess(delta_t)
         !skip some blank or uniteresting integers
         do i=1,12
            read(filePointer) junk
            junk = int_swap_endianess(junk)
         enddo
         !read the number of lines for the title (each line is 80 chars long)
         read(filePointer) nTitle
         nTitle = int_swap_endianess(nTitle)
         !read the title lines
         do i=1,nTitle
            read(filePointer) title
         enddo
         !read two uninteresting integers
         read(filePointer) junk
         junk=int_swap_endianess(junk)
         read(filePointer) junk
         junk=int_swap_endianess(junk)
         !read the number of atoms
         read(filePointer) nAtomTotal
         nAtomTotal = int_swap_endianess(nAtomTotal)
         print*, 'Number of atoms in trajectory file:', nAtomTotal
         !read another integer
         read(filePointer) junk
         junk=int_swap_endianess(junk)
   endif
   

endsubroutine read_dcd_header

!subroutine to read the positions of one step in a NAMD dcd file
subroutine get_selected_dcd_coord(coord,nSelectedAtoms,selectedAtoms,nAtomsTotal,axis, filePointer)
   use endianess
   implicit none
   integer nAtomsTotal
   integer nSelectedAtoms
   integer selectedAtoms(nSelectedAtoms)
   real (kind=8) coord(3,nSelectedAtoms)
   integer filePointer
   integer i, j
   real tempPos
   integer junk
   integer counter
   real (kind=8) junkd
   real (kind=8) axis(3)
   real real_swap_endianess
   integer int_swap_endianess
   real (kind=8) dble_swap_endianess

   if (nativeEndianess.eqv..true.) then
        read(filePointer) junk
        read(filePointer) axis(1)
        read(filePointer) junkd
        read(filePointer) axis(2)
        read(filePointer) junkd
        read(filePointer) junkd
        read(filePointer) axis(3)
        read(filePointer) junk
        do j=1,3
           !nAtoms*4
           counter = 1
           read(filePointer) junk
           do i=1, nAtomsTotal
              read(filePointer) tempPos
              if(i==selectedAtoms(counter)) then
                     coord(j,counter)=dble(tempPos)
                     counter = counter + 1
              endif
           enddo
           !nAtoms*4
           read(filePointer) junk
        enddo
   else
        read(filePointer) junk
        junk = int_swap_endianess(junk)
        read(filePointer) axis(1)
        axis(1) = dble_swap_endianess(axis(1))
        read(filePointer) junkd
        junkd = dble_swap_endianess(junkd)
        read(filePointer) axis(2)
        axis(2) = dble_swap_endianess(axis(2))
        read(filePointer) junkd
        junkd = dble_swap_endianess(junkd)
        read(filePointer) junkd
        junkd = dble_swap_endianess(junkd)
        read(filePointer) axis(3)
        axis(3) = dble_swap_endianess(axis(3))
        read(filePointer) junk
        junk = int_swap_endianess(junk)
        do j=1,3
           !nAtoms*4
           counter = 1
           read(filePointer) junk
           junk = int_swap_endianess(junk)
           do i=1, nAtomsTotal
              read(filePointer) tempPos
              if(i==selectedAtoms(counter)) then
                     tempPos = real_swap_endianess(tempPos)
                     coord(j,counter)=dble(tempPos)
                     counter = counter + 1
              endif
           enddo
           !nAtoms*4
           read(filePointer) junk
           junk = int_swap_endianess(junk)
        enddo
   endif

endsubroutine get_selected_dcd_coord


!subroutine to read the positions of one step in a NAMD dcd file
subroutine read_dcd_step(coord,nAtoms, filePointer)
   use endianess
   implicit none
   integer nAtoms
   real (kind=8) coord(nAtoms,3)
   integer filePointer
   integer i, j
   real tempPos
   integer junk
   integer counter
   real (kind=8) junkd
   real (kind=8) axis(3)
   real real_swap_endianess
   integer int_swap_endianess
   real (kind=8) dble_swap_endianess

   if (nativeEndianess.eqv..true.) then
        read(filePointer) junk
        read(filePointer) axis(1)
        read(filePointer) junkd
        read(filePointer) axis(2)
        read(filePointer) junkd
        read(filePointer) junkd
        read(filePointer) axis(3)
        read(filePointer) junk
        do j=1,3
           !nAtoms*4
           counter = 1
           read(filePointer) junk
           do i=1, nAtoms
              read(filePointer) tempPos
              coord(i,j)=dble(tempPos)
           enddo
           !nAtoms*4
           read(filePointer) junk
        enddo
   else
        read(filePointer) junk
        junk = int_swap_endianess(junk)
        read(filePointer) axis(1)
        axis(1) = dble_swap_endianess(axis(1))
        read(filePointer) junkd
        junkd = dble_swap_endianess(junkd)
        read(filePointer) axis(2)
        axis(2) = dble_swap_endianess(axis(2))
        read(filePointer) junkd
        junkd = dble_swap_endianess(junkd)
        read(filePointer) junkd
        junkd = dble_swap_endianess(junkd)
        read(filePointer) axis(3)
        axis(3) = dble_swap_endianess(axis(3))
        read(filePointer) junk
        junk = int_swap_endianess(junk)
        do j=1,3
           !nAtoms*4
           counter = 1
           read(filePointer) junk
           junk = int_swap_endianess(junk)
           do i=1, nAtoms
              read(filePointer) tempPos
              tempPos = real_swap_endianess(tempPos)
              coord(i,j)=dble(tempPos)
           enddo
           !nAtoms*4
           read(filePointer) junk
           junk = int_swap_endianess(junk)
        enddo
   endif

endsubroutine read_dcd_step




!!!!  Endianess Routines

integer function int_swap_endianess(intIn)
        implicit none
        integer, parameter :: nBytes=4
        integer intIn
        integer (selected_int_kind(1)),dimension(nBytes)      :: intArray
        integer (selected_int_kind(1)),dimension(nBytes)      :: intReversed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
        intArray = TRANSFER( intIn, intArray )
        ! Transfer reversed order bytes to 64 bit REAL space 
        intReversed = intArray(nBytes:1:-1)
        int_swap_endianess = TRANSFER( intReversed, int_swap_endianess )

endfunction int_swap_endianess

real function real_swap_endianess(realIn)
        implicit none
        integer, parameter :: nBytes=4
        real (kind=4) realIn
        integer (selected_int_kind(1)),dimension(nBytes)      :: intIn
        integer (selected_int_kind(1)),dimension(nBytes)      :: intReversed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
        intIn = TRANSFER( realIn, intIn )
        ! Transfer reversed order bytes to 64 bit REAL space 
        intReversed = intIn(nBytes:1:-1)
        real_swap_endianess = TRANSFER( intReversed, real_swap_endianess )

endfunction real_swap_endianess

real (kind=8) function dble_swap_endianess(dbleIn)
        implicit none
        integer, parameter :: nBytes=8
        real (kind=8) dbleIn
        integer (selected_int_kind(1)),dimension(nBytes)      :: intIn
        integer (selected_int_kind(1)),dimension(nBytes)      :: intReversed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
        intIn = TRANSFER( dbleIn, intIn )
        ! Transfer reversed order bytes to 64 bit REAL space 
        intReversed = intIn(nBytes:1:-1)
        dble_swap_endianess = TRANSFER( intReversed, dble_swap_endianess )

endfunction dble_swap_endianess


