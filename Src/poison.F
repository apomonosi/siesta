! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine poison( CELL, N1, N2, N3, Mesh, RHO, U, V, STRESS, 
     &                   NSM )

C *********************************************************************
C SOLVES POISSON'S EQUATION.
C ENERGY AND POTENTIAL RETURNED IN RYDBERG UNITS.
C WRITTEN BY J.M.SOLER, JUNE 1995.
C **************** INPUT **********************************************
C REAL*8  CELL(3,3)     : UNIT CELL VECTORS
C INTEGER N1,N2,N3      : NUMBER OF MESH DIVISIONS IN EACH CELL VECTOR
C INTEGER Mesh(3)       : Number of global mesh divisions
C REAL*4  RHO(N1,N2,N3) : DENSITIY AT MESH POINTS
C **************** OUTPUT *********************************************
C REAL*8  U             : ELECTROSTATIC ENERGY (IN RY)
C REAL*4  V(N1,N2,N3)   : ELECTROSTATIC POTENTIAL (IN RY)
C                         V AND RHO MAY BE THE SAME PHYSICAL ARRAY
C REAL*8  STRESS(3,3) : Electrostatic-energy contribution to stress
C                       tensor (in Ry/Bohr**3) assuming constant density
C                       (not charge), i.e. r->r' => rho'(r') = rho(r)
C                       For plane-wave and grid (finite difference)
C                       basis sets, density rescaling gives an extra
C                       term (not included) equal to -2*U/cell_volume
C                       for the diagonal elements of stress. For other
C                       basis sets, the extra term is, in general:
C                       IntegralOf( V * d_rho/d_strain ) / cell_volume
C INTEGER NSM           : Number of sub-mesh points per mesh point
C                       : along each axis
C *********************************************************************

C     Modules
      use precision,   only : dp, grid_p
      use parallel,    only : Node, Nodes
      use moreMeshSubs, only : UNIFORM, getMeshBox
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc
      use m_fft,       only : fft, fft_init, getFFTLims
      use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
      use cellsubs,    only : volcel  ! Finds unit cell volume
      use m_chkgmx,    only : chkgmx  ! Checks planewave cutoff of a mesh
#ifdef MPI
      use mpi_siesta
#endif
      implicit          none
C     Input/output variables
      integer               :: N1, N2, N3, Mesh(3), NSM
      real(grid_p)          :: RHO(N1*N2*N3), V(N1*N2*N3)
      real(dp)              :: CELL(3,3), STRESS(3,3), U
C     Local variables
      integer               :: I, I1, I2, I3, IX, J, J1, J2, J3, JX,
     &                         NP, NG, NG1, NG2, NG3, j1l, J2L, J3L, JJ
      integer               :: lbox(2,3)
      integer,      pointer :: box(:,:,:)
      real(dp)              :: C, B(3,3), DU, G(3), G2, G2MAX, 
     &                         PI, VG, VOLUME, PI8
      real(grid_p), pointer :: CG1(:,:), CG2(:,:)
      real(dp),   parameter :: K0(3)= (/0.d0,0.d0,0.d0/), TINY= 1.d-15
      integer :: sz1, sz2
#ifdef MPI
      integer               :: MPIerror
#endif

#ifdef DEBUG
      call write_debug( '    PRE POISON' )
#endif
#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 2 )
#endif

C     Start time counter
      call timer( 'POISON', 1 )

C     Find fourier transform of density

      ! For this use case we can avoid going back to the x-pencil
      ! distribution immediately, reducing the number of communications

      call fft_init( Mesh, sz1, sz2, pin_x_distribution=.false. )

C     Allocate local memory
      nullify( CG1, CG2 )
      call re_alloc( CG1, 1, 2, 1, sz1, 'CG1', 'poison' )
      call re_alloc( CG2, 1, 2, 1, sz2, 'CG2', 'poison' )

C     Find unit cell volume
      VOLUME = VOLCEL( CELL )

C     Find reciprocal lattice vectors
      call reclat(CELL, B, 1 )

C     Find maximun planewave cutoff
      NP = N1 * N2 * N3
      G2MAX = 1.0e30_dp
      call CHKGMX( K0, B, Mesh, G2MAX )

C     Copy density to complex array
      do I = 1, NP
        CG1(1,I) = RHO(I)
        CG1(2,I) = 0.0_grid_p
      enddo

C     Find fourier transform of density
      call fft( CG1, CG2, -1 )

C     Initialize stress contribution
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = 0.0_dp
        enddo
      enddo

C     Work out processor grid dimensions
      call getFFTLims( lbox )
      lbox = lbox -1

C     Multiply by 8*PI/G2 to get the potential
      PI  = 4.0_dp * atan(1.0_dp)
      PI8 = PI * 8._dp
      U = 0.0_dp

      NG1 = Mesh(1)
      NG2 = Mesh(2)
      NG3 = Mesh(3)

C     getFFTLims returns the limits of current distribution
      J = 1
      do J3 = lbox(1,3), lbox(2,3)           ! Z direction
        J3L = J3 - lbox(1,3)
        I3 = merge( J3 - NG3, J3, J3.gt.NG3/2 )
        do J2 = lbox(1,2), lbox(2,2)               ! Y directiron
          J2L = J2 - lbox(1,2)
          I2 = merge( J2 - NG2, J2, J2.gt.NG2/2 )
          do J1 = lbox(1,1), lbox(2,1)             ! X directiron
            J1L = J1 - lbox(1,1)
            I1 = merge( J1 - NG1, J1, J1.gt.NG1/2 )

            G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
            G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
            G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3
            G2 = G(1)**2 + G(2)**2 + G(3)**2
            if (G2.LT.G2MAX .AND. G2.GT.TINY) then
              VG = PI8 / G2
              DU = VG * ( CG1(1,J)**2 + CG1(2,J)**2 )
              U = U + DU
              C = 2.0_dp * DU / G2
              DO IX = 1,3
                DO JX = 1,3
                  STRESS(JX,IX) = STRESS(JX,IX) + C * G(IX) * G(JX)
                ENDDO
              ENDDO
              CG1(1,J) = VG * CG1(1,J)
              CG1(2,J) = VG * CG1(2,J)
            else
              CG1(1,J) = 0.0_dp
              CG1(2,J) = 0.0_dp
            endif
            J = J+1
          enddo
        enddo
      enddo

      NG = Mesh(1)*Mesh(2)*Mesh(3)
      U = 0.5_dp * U * VOLUME / DBLE(NG)**2
      C = 0.5_dp / DBLE(NG)**2
      do IX = 1,3
        do JX = 1,3
          STRESS(JX,IX) = C * STRESS(JX,IX)
        enddo
        STRESS(IX,IX) = STRESS(IX,IX) + U / VOLUME
      enddo
 
C     Go back to real space
      call fft( CG1, CG2, +1 )

C     Copy potential to array V
      do I = 1, NP
        V(I) = CG1(1,I)
      enddo
 
C     Free local memory
      call de_alloc( CG1, 'CG1', 'poison' )
      call de_alloc( CG2, 'CG2', 'poison' )

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif

C     Stop time counter
      call timer( 'POISON', 2 )

#ifdef DEBUG
      call write_debug( '    POS POISON' )
#endif
      end subroutine poison
