! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine read_xc_info()
C
C     Reads the exchange-correlation functional information
C     and calls setXC to store it
C
      use precision, only : dp
      use gridXC,  only : setXC=>gridxc_setXC
      use fdf
      use parallel,  only : Node
      use sys,       only : die

      implicit none

C     LOCAL variables

      integer, parameter      :: MaxFunc = 20
      integer                 :: nXCfunc
      character(len=50)       :: XCauth(MaxFunc)
      character(len=20)       :: XCfunc(MaxFunc)
      real(dp)                :: XCweightX(MaxFunc)
      real(dp)                :: XCweightC(MaxFunc)

        integer            :: n
        integer            :: ni
        integer            :: nn
        integer            :: nr

        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline

        logical  :: found_block
        
!      Read XC functionals

        found_block = .false.
        
        if (fdf_block('xc.mix',bfdf)) then
           found_block = .true.
        else if (fdf_block('xc.cocktail',bfdf)) then
           call message(' (info)',
     $          'XC.cocktail block is deprecated. Use XC.mix')
           found_block = .true.
        else if (fdf_block('xc.hybrid',bfdf)) then
           call message(' (info)',
     $          'XC.hybrid block is deprecated. Use XC.mix')
           found_block = .true.
        endif


        if (found_block) then
           
          if (.not. fdf_bline(bfdf,pline)) then
            call die('read_xc: ERROR no data in xc.mix block')
          endif
          ni = fdf_bnintegers(pline)

          if (ni .eq. 0) then
            call die('read_xc: Number of functionals missing in '
     $            // 'xc.mix')
          endif
          nXCfunc = abs(fdf_bintegers(pline,1))
          if (nXCfunc .gt. MaxFunc) then
            call die('read_xc: Too many (>20) functionals in xc.mix')
          endif
          do n= 1, nXCfunc
            if (.not. fdf_bline(bfdf,pline)) then
              call die('read_xc: Too few functional lines in xc.mix')
            endif
            nn = fdf_bnnames(pline)
            nr = fdf_bnreals(pline)

            if (nn .gt. 0) then
              XCfunc(n) = fdf_bnames(pline,1)
            else
              XCfunc(n) = 'LDA'
            endif
            if (nn .gt. 1) then
              XCauth(n) = fdf_bnames(pline,2)
            else
              XCauth(n) = 'PZ'
            endif
            if (nr .gt. 1) then
              XCweightX(n) = fdf_breals(pline,1)
              XCweightC(n) = fdf_breals(pline,2)
            elseif (nr .eq. 1) then
              XCweightX(n) = fdf_breals(pline,1)
              XCweightC(n) = fdf_breals(pline,1)
            else
              XCweightX(n) = 1.0_dp
              XCweightC(n) = 1.0_dp
            endif
          enddo
          call fdf_bclose(bfdf)
        else
          nXCfunc = 1 
          XCfunc(1) = fdf_string('xc.functional','LDA')
          XCauth(1) = fdf_string('xc.authors','PZ')
          XCweightX(1) = 1.0_dp
          XCweightC(1) = 1.0_dp
        endif

        ! Announce the deprecated flags (if used)...
        ! ... after closing the block
        call fdf_deprecated('XC.hybrid','XC.mix')
        call fdf_deprecated('XC.cocktail','XC.mix')

C Output data for mix of functionals
        if ((nXCfunc .gt. 1) .and. (Node .eq. 0)) then
          write(6,'(/,''xc:'')')
          write(6,'(''xc: Exchange-correlation functional mix:'')
     .      ')
          write(6,'(''xc:'')')
          write(6,'(''xc: Number     Functional     Authors  '',
     .      ''   Weight(Ex)   Weight(Ec)'')')
          do n = 1,nXCfunc
          write(6,'(''xc: '',i4,3x,a20,2x,a20,3x,f5.3,8x,f5.3)')
     .        n,XCfunc(n),XCauth(n),XCweightX(n),XCweightC(n)
          enddo
          write(6,'(''xc:'')')
        endif

C     Store information in module

      call setXC (n=nXCfunc, func=XCfunc, auth=XCauth,
     $            wx=XCweightX, wc=XCweightC)

      end subroutine read_xc_info

