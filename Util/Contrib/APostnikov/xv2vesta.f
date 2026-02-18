C
C    xf2xsf,  a quick script to reformat the XV file 
C             with lattice/coordinates into a XCrysden xsf file.
C
C             Written by Andrei Postnikov, Mar 2006   Vers_0.3
C             apostnik@uos.de
C


C
C    xv2vesta,  a quick script to reformat the XV file 
C             with lattice/coordinates into a VESTA vesta file.
C
C             Written by Fawzi Kaoui, June 2012   Vers_0.3
C             fawzi.kaoui@gmail.com
C
      program xv2vesta
      implicit none
      integer ii1,io1,jj,ialloc
      parameter (ii1=11,io1=14)
      integer ii,iat,nat,nz
      integer, allocatable :: iz(:),ityp(:)
      character inpfil*60,outfil*60,syslab*30,suffix*6
      character(len=2),allocatable :: label(:)
      double precision b2ang,cc_bohr(3,3),cc_ang(3,3),cc_angi(3,3),
     .                 vector_aa(3),vector_bb(3),vector_cc(3),coordr(3),
     .                 rad,alpha,beta,jamma,xx,yy,zz,aa,bb,cc,occ
      double precision, allocatable :: mass(:),coord(:,:)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
C      parameter (rad  =57.295780)  !  radian to degree
      external opnout,inver3,test_xv,read_xv
C     external inver3(cc_ang,cc_angi)              ! inverse matrix
C
C     string manipulation functions in Fortran used below:
C     len_trim(string): returns the length of string 
C                       without trailing blank characters,
C     trim(string)    : returns the string with railing blanks removed
      
      write (6,701,advance="no")
  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
      read (5,*) syslab
C     inpfil = syslab(1:len_trim(syslab))//'.XV'
      inpfil = trim(syslab)//'.XV'
      open (ii1,file=inpfil,form='formatted',status='old',err=801)
      write (6,*) 'Found and opened: ',inpfil
       call test_xv(ii1,nat)
      allocate (ityp(nat))
      allocate (iz(nat))
      allocate (mass(nat))
      allocate (label(nat))
      allocate (coord(1:3,1:nat),STAT=ialloc)
      if (ialloc.ne.0) then
         write (6,*) ' Fails to allocate space for ',nat,' atoms.'
         stop
      endif
      call read_xv(ii1,nat,ityp,iz,cc_ang,mass,label,coord)
C     bohr to Angstroem for cc_ang (translation vectors coordinates)
C     and coord : Cartesian coordinates of atoms in Angstroem
      call inver3(cc_ang,cc_angi)                              
      close (ii1)

C --- open output file:
C     outfil = syslab(1:len_trim(syslab))//'.vesta'
      outfil = trim(syslab)//'.vesta'
      call opnout(io1,outfil)
C --  extract crystal cell parameters
C --  read in translation vectors, convert into Ang:
C     Angstroem
      aa = sqrt(cc_ang(1,1)**2+cc_ang(2,1)**2+cc_ang(3,1)**2)     
      bb = sqrt(cc_ang(1,2)**2+cc_ang(2,2)**2+cc_ang(3,2)**2)     
      cc = sqrt(cc_ang(1,3)**2+cc_ang(2,3)**2+cc_ang(3,3)**2)     
      vector_aa=(/cc_ang(1,1),cc_ang(2,1),cc_ang(3,1)/)
      vector_bb=(/cc_ang(1,2),cc_ang(2,2),cc_ang(3,2)/)
      vector_cc=(/cc_ang(1,3),cc_ang(2,3),cc_ang(3,3)/)
C     dimensionless
      xx = dot_product(vector_bb,vector_cc)/(bb*cc)                
      yy = dot_product(vector_aa,vector_cc)/(aa*cc)                
      zz = dot_product(vector_aa,vector_bb)/(aa*bb)                
      rad = 45./atan(1.0)
C     degree
      alpha = acos(xx)*rad     ! degree
      beta  = acos(yy)*rad     ! degree
      jamma = acos(zz)*rad     ! degree
C     write (io1,*) '# crystal structure from ',
C    .              syslab(1:len_trim(syslab)),'.XV'
      write(io1,202)
C     does not leave space after # if the second cell
  202 format ('#VESTA_FORMAT_VERSION 3.0.0')            
      write (io1,*) '# crystal structure from ',trim(syslab),'.XV'
      write (io1,*) trim(syslab),'.XV'
      write (io1,*) 'CRYSTAL'
      write (io1,*) '# Cell vectors in Angstroem and angles in degrees:'
      write (io1,*) 'CELLP'
      write  (io1,210) aa,bb,cc,alpha,beta,jamma
  210 format(6f12.7,/, 6('   0.0000000'))
      call inver3(cc_ang,cc_angi)
      write (io1,*) '# relative coordinates of the atoms dimensionless:'
      write (io1,*) 'STRUC'
      occ=1.0000

      do iat=1,nat                                            
        coordr=matmul (cc_angi, coord(:,iat))                    
        write (io1,201) iat,label(iat),label(iat),occ,coordr 
C    coord: Cartesian coordinates of atoms in Angstroem
C    cc_angi: Cartesian coordinates of reverse translation vectors 
C             in 1/Angstroem
C    coordr: relative coordinates of the atoms dimensionless 
  201  format(i9,2a6,f8.4, 3f14.7,/25x, 3('     0.0000000'))  
C 201 format(i5,2a7,f8.4,3f12.6)
      enddo
      write (io1,*) '  0 0 0 0 0 0 0 0 0'

      close (ii1)
      close (io1)
      stop

  801 write (6,*) ' Error opening file ',
     .              trim(inpfil),' as old formatted'
      stop
  802 write (6,*) ' Error opening file ',
     .              trim(outfil),' as new formatted'
      stop
  803 write (6,*) ' End/Error reading XV for cell vector ',ii
      stop
  804 write (6,*) ' End/Error reading XV for number of atoms line'
      stop
  805 write (6,*) ' End/Error reading XV for atom number ',iat
      stop
      end

