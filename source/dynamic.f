c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular or stochastic dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program dynamic
      use sizes
      use atoms
      use bath
      use bndstr
      use bound
      use inform
      use iounit
      use keys
      use mdstuf
      use potent
      use solute
      use stodyn
      use usage
      implicit none
      integer i,istep,nstep
      integer mode,next
      real*8 dt,dtdump
      logical exist
      character*20 keyword
      character*120 record
      character*120 string

      ! Variables for key file management
      character*120 pert_k(25000), fpert_k(25000), bpert_k(25000)
      integer pert_nk, fpert_nk, bpert_nk
      logical pert_exist, fpert_exist, bpert_exist

c     Load key files into memory if they exist
      inquire(file='pert.key', exist=pert_exist)
      if (pert_exist) then
         call loadk('pert.key', pert_k, pert_nk)
      endif

      inquire(file='fpert.key', exist=fpert_exist)
      if (fpert_exist) then
         call loadk('fpert.key', fpert_k, fpert_nk)
      endif

      inquire(file='bpert.key', exist=bpert_exist)
      if (bpert_exist) then
         call loadk('bpert.key', bpert_k, bpert_nk)
      endif

c     Set up the structure and mechanics calculation
      call initial
      call getxyz
      call mechanic

c     Initialize temperature, pressure and coupling baths
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.

c     Check for keywords with altered parameters
      integrate = 'BEEMAN'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword(record,integrate,next)
            call upcase(integrate)
         endif
      enddo

c     Initialize simulation length (time steps)
      nstep = -1
      call nextarg(string,exist)
      if (exist) read(string,*,err=10,end=10) nstep
   10 continue
      dowhile (nstep .lt. 0)
         write(iout,20)
   20    format(/,' Enter Number of Dynamics Steps :  ',$)
         read(input,30,err=40) nstep
   30    format(i10)
         if (nstep .lt. 0) nstep = 0
   40    continue
      enddo

c     Get time step length in picoseconds
      dt = -1.0d0
      call nextarg(string,exist)
      if (exist) read(string,*,err=50,end=50) dt
   50 continue
      do while (dt .lt. 0.0d0)
         write(iout,60)
   60    format(/,' Enter Time Step Length (fs) [1.0] :  ',$)
         read(input,70,err=80) dt
   70    format(f20.0)
         if (dt .le. 0.0d0) dt = 1.0d0
   80    continue
      enddo
      dt = 0.001d0 * dt

c     Enforce bounds on coupling times
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)

c     Set time between coordinate dumps
      dtdump = -1.0d0
      call nextarg(string,exist)
      if (exist) read(string,*,err=90,end=90) dtdump
   90 continue
      do while (dtdump .lt. 0.0d0)
         write(iout,100)
  100    format(/,' Time between Dumps (ps) [0.1] :  ',$)
         read(input,110,err=120) dtdump
  110    format(f20.0)
         if (dtdump .le. 0.0d0) dtdump = 0.1d0
  120    continue
      enddo
      iwrite = nint(dtdump/dt)

c     Get ensemble choice for periodic system
      if (use_bounds) then
         mode = -1
         call nextarg(string,exist)
         if (exist) read(string,*,err=130,end=130) mode
  130    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write(iout,140)
  140       format(/,' Available Ensembles:',//,4x,'(1) NVE',/,
     &              4x,'(2) NVT',/,4x,'(3) NPH',/,4x,'(4) NPT',//,
     &              ' Enter Choice [1] :  ',$)
            read(input,150,err=160) mode
  150       format(i10)
            if (mode .le. 0) mode = 1
  160       continue
         enddo
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write(iout,170)
  170          format(/,' Switching to NPT Ensemble')
            endif
         endif
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg(string,exist)
            if (exist) read(string,*,err=180,end=180) kelvin
  180       continue
            do while (kelvin .lt. 0.0d0)
               write(iout,190)
  190          format(/,' Desired Temperature (K) [298] :  ',$)
               read(input,200,err=210) kelvin
  200          format(f20.0)
               if (kelvin .le. 0.0d0) kelvin = 298.0d0
  210          continue
            enddo
         endif
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg(string,exist)
            if (exist) read(string,*,err=220,end=220) atmsph
  220       continue
            do while (atmsph .lt. 0.0d0)
               write(iout,230)
  230          format(/,' Desired Pressure (Atm) [1.0] :  ',$)
               read(input,240,err=250) atmsph
  240          format(f20.0)
               if (atmsph .le. 0.0d0) atmsph = 1.0d0
  250          continue
            enddo
         endif
      endif

c     For nonperiodic systems
      if (.not. use_bounds) then
         mode = -1
         call nextarg(string,exist)
         if (exist) read(string,*,err=260,end=260) mode
  260    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write(iout,270)
  270       format(/,' Available Modes:',//,4x,'(1) Constant E',/,
     &              4x,'(2) Constant T',//,' Enter Choice [1] :  ',$)
            read(input,280,err=290) mode
  280       format(i10)
            if (mode .le. 0) mode = 1
  290       continue
         enddo
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg(string,exist)
            if (exist) read(string,*,err=300,end=300) kelvin
  300       continue
            do while (kelvin .lt. 0.0d0)
               write(iout,310)
  310          format(/,' Desired Temperature (K) [298] :  ',$)
               read(input,320,err=330) kelvin
  320          format(f20.0)
               if (kelvin .le. 0.0d0) kelvin = 298.0d0
  330          continue
            enddo
         endif
      endif

c     Initialize constraints and dynamics
      call shakeup
      call mdinit

c     Print dynamics header
      if (integrate .eq. 'VERLET') then
         write(iout,340)
  340    format(/,' MD via Velocity Verlet')
      else if (integrate .eq. 'STOCHASTIC') then
         write(iout,350)
  350    format(/,' Stochastic Dynamics via Verlet')
      else if (integrate .eq. 'BUSSI') then
         write(iout,360)
  360    format(/,' MD via Bussi-Parrinello NPT')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write(iout,370)
  370    format(/,' MD via Nose-Hoover NPT')
      else if (integrate .eq. 'GHMC') then
         write(iout,380)
  380    format(/,' Stochastic via GHMC')
      else if (integrate .eq. 'RIGIDBODY') then
         write(iout,390)
  390    format(/,' MD via Rigid Body')
      else if (integrate .eq. 'RESPA') then
         write(iout,400)
  400    format(/,' MD via r-RESPA MTS')
      else
         write(iout,410)
  410    format(/,' MD via Modified Beeman')
      endif

c     Integrate equations of motion
      do istep = 1, nstep
         if (integrate .eq. 'VERLET') then
            call verlet(istep,dt)
         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep(istep,dt)
         else if (integrate .eq. 'BUSSI') then
            call bussi(istep,dt)
         else if (integrate .eq. 'NOSE-HOOVER') then
            call nose(istep,dt)
         else if (integrate .eq. 'GHMC') then
            call ghmcstep(istep,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep(istep,dt)
         else if (integrate .eq. 'RESPA') then
            call respa(istep,dt)
         else
            call beeman(istep,dt)
         endif

c        Calculate perturbed energies
         if (pert_exist) then
            call calcpe(istep, pert_k, pert_nk)
         endif
         if (fpert_exist) then
            call calcfpe(istep, fpert_k, fpert_nk)
         endif
         if (bpert_exist) then
            call calcbpe(istep, bpert_k, bpert_nk)
         endif
      enddo

c     Final tasks
      call final
      end

c     Subroutine to calculate PE with pert.key
      subroutine calcpe(istep, keylines, nkeys)
      use sizes
      use atoms
      use analyz
      use potent
      use keys
      implicit none
      integer istep
      character*120 keylines(25000)
      integer nkeys
      real*8 pe_perturbed
      integer iunit
      logical file_exists
      real*8 energy
      character*120 orig_key(25000)
      integer orig_nk
      logical first_call
      save orig_key, orig_nk, first_call
      data first_call /.true./

      if (first_call) then
         call savekeys(orig_key, orig_nk)
         first_call = .false.
      endif

      inquire(file='pe.log', exist=file_exists)
      if (file_exists) then
         iunit = 30
         open(unit=iunit, file='pe.log', status='old',
     &     position='append')
      else
         iunit = 30
         open(unit=iunit, file='pe.log', status='new')
         write(iunit, '(A)') "step  PE_perturbed"
      endif

      call loadkeys(keylines, nkeys)
      call mechanic 
      pe_perturbed = energy()
      call restorekeys(orig_key, orig_nk)

      write(iunit, '(I5, F16.8)') istep, pe_perturbed
      close(iunit)
      end

c     Subroutine to calculate FPE with fpert.key
      subroutine calcfpe(istep, keylines, nkeys)
      use sizes
      use atoms
      use analyz
      use potent
      use keys
      implicit none
      integer istep
      character*120 keylines(25000)
      integer nkeys
      real*8 pe_perturbed
      integer iunit
      logical file_exists
      real*8 energy
      character*120 orig_key(25000)
      integer orig_nk
      logical first_call
      save orig_key, orig_nk, first_call
      data first_call /.true./

      if (first_call) then
         call savekeys(orig_key, orig_nk)
         first_call = .false.
      endif

      inquire(file='fpe.log', exist=file_exists)
      if (file_exists) then
         iunit = 30
         open(unit=iunit, file='fpe.log', status='old',
     &     position='append')
      else
         iunit = 30
         open(unit=iunit, file='fpe.log', status='new')
         write(iunit, '(A)') "step  PE_perturbed"
      endif

      call loadkeys(keylines, nkeys)
      call mechanic 
      pe_perturbed = energy()
      call restorekeys(orig_key, orig_nk)

      write(iunit, '(I5, F16.8)') istep, pe_perturbed
      close(iunit)
      end

c     Subroutine to calculate BPE with bpert.key
      subroutine calcbpe(istep, keylines, nkeys)
      use sizes
      use atoms
      use analyz
      use potent
      use keys
      implicit none
      integer istep
      character*120 keylines(25000)
      integer nkeys
      real*8 pe_perturbed
      integer iunit
      logical file_exists
      real*8 energy
      character*120 orig_key(25000)
      integer orig_nk
      logical first_call
      save orig_key, orig_nk, first_call
      data first_call /.true./

      if (first_call) then
         call savekeys(orig_key, orig_nk)
         first_call = .false.
      endif

      inquire(file='bpe.log', exist=file_exists)
      if (file_exists) then
         iunit = 30
         open(unit=iunit, file='bpe.log', status='old',
     &     position='append')
      else
         iunit = 30
         open(unit=iunit, file='bpe.log', status='new')
         write(iunit, '(A)') "step  PE_perturbed"
      endif

      call loadkeys(keylines, nkeys)
      call mechanic 
      pe_perturbed = energy()
      call restorekeys(orig_key, orig_nk)

      write(iunit, '(I5, F16.8)') istep, pe_perturbed
      close(iunit)
      end

c     Subroutine to save current key state
      subroutine savekeys(saved_key, saved_nk)
      use keys
      implicit none
      character*120 saved_key(25000)
      integer saved_nk
      integer i

      saved_nk = nkey
      do i = 1, nkey
         saved_key(i) = keyline(i)
      enddo
      end

c     Subroutine to load key file into memory
      subroutine loadk(filename, memory, nkeys)
      implicit none
      character*(*) filename
      character*120 memory(25000)
      integer nkeys
      integer i, ios

      open(unit=20, file=filename, status='old', iostat=ios)
      if (ios .ne. 0) then
         print *, "Error opening ", filename
         stop
      endif

      nkeys = 0
      do i = 1, 25000
         read(20, '(A)', iostat=ios) memory(i)
         if (ios .ne. 0) exit
         nkeys = nkeys + 1
      enddo
      close(20)
      end

c     Subroutine to load keys from memory
      subroutine loadkeys(memory, nkeys)
      use keys
      implicit none
      character*120 memory(25000)
      integer nkeys
      integer i

      nkey = nkeys
      do i = 1, nkey
         keyline(i) = memory(i)
      enddo
      end

c     Subroutine to restore original key state
      subroutine restorekeys(saved_key, saved_nk)
      use keys
      implicit none
      character*120 saved_key(25000)
      integer saved_nk
      integer i

      nkey = saved_nk
      do i = 1, nkey
         keyline(i) = saved_key(i)
      enddo
      call mechanic
      end
