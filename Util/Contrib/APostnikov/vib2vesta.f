C
C   vib2xsf,  a script to transform phonon eigenvectors 
C             provided by "vibrator"
C             either in XBS file with arrows,
C             or into animated AXBS file vizualizing vibration modes,
C             to be further worked on with XCrySDen
C
C             Written by Andrei Postnikov, Mar 2006   Vers_0.3
C             apostnik@uos.de
C


C
C   vib2vesta,  a script to transform phonon eigenvectors 
C             provided by "vibrator"
C             either in XBS file with arrows,
C             or into static AXBS file vizualizing vibration modes,
C             to be further worked on with VESTA
C
C             Written by Fawzi Kaoui June 2012, June 2012   Vers_0.3
C             fawzi.kaoui@gmail.com
C
      program vib2vesta
      implicit none
      integer ii1,ii2,io1,io2,is1,pen
      parameter (ii1=11,ii2=12,io1=13,io2=14,is1=15)
      integer ialloc,iat,jat,nat,ivmin,ivmax,iev,nsteps,nbox,
     .        color(3)
      integer, allocatable :: ityp(:),iz(:) 
      double precision cc_ang(3,3),dist,dstmin,
     .                 obox(3),rbox(3,3),rboxi(3,3),fact,diam
      double precision, allocatable :: mass(:),coord(:,:),freq(:),
     .                                 evr(:,:,:),disp(:,:,:)
      character inpfil*60,outfil*60,syslab*30,alab*1,
     .          itochar*10,modlab*10
      character(len=2), allocatable :: label(:)
      external test_xv,read_xv,read_ev,itochar,inver3,opnout,fillbox,displa
      write (6,701,advance="no")
  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
      read (5,*) syslab
C     inpfil = syslab(1:len_trim(syslab))//'.XV'
      inpfil = trim(syslab)//'.XV'
C     open ii1
      open (ii1,file=inpfil,form='formatted',status='old',err=801)    
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
      close (ii1)                                                     
C ---   find smallest interatomic distance, 
C       for nice scaling of displacements:
      dstmin = 1.0d6  !     distance inter atomiques
      do iat=1,nat
      do jat=iat+1,nat
        dist = (coord(1,iat)-coord(1,jat))**2 +
     +         (coord(2,iat)-coord(2,jat))**2 +
     +         (coord(3,iat)-coord(3,jat))**2 
        if (dist.lt.dstmin) dstmin = dist
      enddo
      enddo
      dstmin = sqrt(dstmin)
      write (6,*)'  dstmin=',dstmin   
C --- set up and fill output box:
      call makebox_vesta(obox,rbox)
C     makebox : (boite)
C     look inside a subroutine makebox ago axes
C     regarder makebox a l'interieur il ya les axes                        
C --- invert the box vectors; will need it in a minute...  
      call inver3(rbox,rboxi)                              
C --- write selected atoms first into a scratch file (is1), for the case
C     there are zero. Then the label 'ATOMS' with no  atoms following
C     will crash vesta.
C     open (is1, file='tmpfil',form='formatted',status='scratch')
      open (is1, file='IS1',   form='formatted',status='unknown')
      call fillbox(is1,obox,rbox,rboxi,cc_ang,nat,coord,nbox)     
      
      write (6,702,advance="no")
  702 format(' Specify SystemLabel of vibrator calculation ',
     .       "(or 'siesta' if none): ")
      read (5,*) syslab
C     inpfil = syslab(1:len_trim(syslab))//'.vectors'
      inpfil = trim(syslab)//'.vectors'
      open (ii2,file=inpfil,form='formatted',status='old',err=801)  
      write (6,703) nat*3
  703 format(' select first and last modes (out of ',i9,
     .       ' ) for analysis.'/' A separate vesta/AVESTA file',
     .       ' will be created for each mode.')
  101 write (6,704,advance="no")           
  704 format(' First mode : ')
C     modes de vibration (ivmax,ivmin)
      read (5,*) ivmin                     
      if (ivmin.le.0.or.ivmin.gt.nat*3) goto 101 
  102 write (6,705,advance="no") 
  705 format('  Last mode : ')
      read (5,*) ivmax
      if (ivmax.lt.ivmin.or.ivmax.gt.nat*3) goto 102 
      allocate (freq(ivmin:ivmax))
      allocate (evr(1:3,nat,ivmin:ivmax))
      allocate (disp(1:3,nat,ivmin:ivmax),STAT=ialloc)
      if (ialloc.ne.0) then
         write (6,*) ' Fails to allocate space for ',
     .                 ivmax-ivmin+1,' modes'
         stop
      endif

      call read_ev(ii2,nat,ivmin,ivmax,evr,freq)
      close (ii2)                                                  
C --- recover atom displacements within each mode with convenient scaling:
C     displacement en fonction de distance inter atomique dstmin et masse
      call displa(nat,ivmin,ivmax,mass,dstmin,evr,disp)            
c      write (6,*) ' Enter number of animation steps:'    
c      read  (5,*) nsteps
C --- write each pattern into a separate file:
      write (6,*) ' Parameters of the ARROWS :'
      write (6,*) ' Enter first ,second and third ',
     .            ' numbers color of ARROWS between 0 and 255 :'
      read  (5,*) color
      write (6,*) ' Enter value of invoice length to ARROWS :'
      read  (5,*) fact
      write (6,*) ' Enter value of diameter to ARROWS :'
      read  (5,*) diam
      write (6,*) ' Enter 0 or 1 , 1 to penetrating atoms ',
     .            ' 0 to not penetrating atoms:'
      read  (5,*) pen

      do iev=ivmin,ivmax
        modlab = itochar(iev)  
C       outfil = syslab(1:len_trim(syslab))//'.Mode_'//
C    .           modlab(1:len_trim(modlab))//'.vesta'
        outfil = trim(syslab)//'.Mode_'//trim(modlab)//'.vesta'   
        write (6,*) ' outfil=',outfil
        call opnout(io1,outfil)
C     io1 ouvert
        call write_arrow(io1,is1,nbox,ivmin,ivmax,iev,             
     .                   cc_ang,nat,iz,freq,disp,label,
     .                   syslab,color,fact,diam,pen)
        close (io1)
C     ajouter syslab
C       outfil = syslab(1:len_trim(syslab))//'.Mode_'//      
C     .           modlab(1:len_trim(modlab))//'a.vesta'
C       outfil = trim(syslab)//'.Mode_'//trim(modlab)//'a.vesta'  
C       call opnout(io2,outfil)
C       call write_movie(io2,is1,nbox,ivmin,ivmax,iev,            
C     .                   cc_ang,nat,iz,freq,disp,nsteps)
C       close (io2)                                          
      enddo
      close (is1)
      stop

  801 write (6,*) ' Error opening file ',
     .              trim(inpfil),' as old formatted'
      stop

      end
