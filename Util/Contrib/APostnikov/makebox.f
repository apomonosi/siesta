C
      subroutine makebox(obox,rbox)
      use units, only : Ang
C
C     asks for origin ans spanning vectors,
C     constructs the output box and its inversion.
C
C     Input:   none (interactive in/out)
C     Output:  obox - origin of output box,
C              rbox - its three spanning vectors
C
      implicit none
      integer ii
      double precision obox(3),rbox(3,3)
      character unitlab*1,labunit*4
      logical unitb

      write (6,702)
  101 write (6,703,advance="no")
      read (5,*) unitlab
      if (unitlab.eq.'B'.or.unitlab.eq.'b') then
        unitb = .true.
        labunit = 'Bohr'
      elseif (unitlab.eq.'A'.or.unitlab.eq.'a') then
        unitb = .false.
        labunit = 'Ang '
      else
        write (6,*) ' Sorry, no third choice.'
        goto 101
      endif
      write (6,704,advance="no") labunit
      read  (5,*) (obox(ii),ii=1,3)
      write (6,705,advance="no") '1st',labunit
      read  (5,*) (rbox(ii,1),ii=1,3)
      write (6,705,advance="no") '2nd',labunit
      read  (5,*) (rbox(ii,2),ii=1,3)
      write (6,705,advance="no") '3rd',labunit
      read  (5,*) (rbox(ii,3),ii=1,3)
      if (unitb) then
C       transform anything into Ang, as is standard in XCrysden:
        obox=obox/Ang
        rbox=rbox/Ang
      endif
      return

  702 format(" Now define the grid cell for your XCrysDen plot.",/
     .       " Note that it can be arbitrarily chosen",
     .       " with respect to the Siesta simulation cell,",
     .       " and it needs not to be orthogonal.",
     .       " We'll define it by the origin point and three",
     .       " spanning vectors. They can be given in Bohr or Ang.")
  703 format (" Would you use Bohr (B) or Ang (A) ? ")
  704 format(' Enter origin point in ',a4,' : ')
  705 format(' Enter ',a3,' spanning vector in ',a4,' : ')

      end
!AAT adds a little different subroutine to make box to VESTA
      subroutine makebox_vesta(obox,rbox)
C
C     asks for reference point (someplace inside or outside the box),
C     three direction vectors passing through it,
C     min and max limits along each vector, measured from the reference point.
C     Reformulates it into the box origin and its spanning vectors. 
C
C     Input:   none (interactive in/out)
C     Output:  obox - origin of output box,
C              rbox - its three spanning vectors
C
      implicit none
      integer ii
      double precision obox(3),rbox(3,3),b2ang,small,
     .                 bref(3),dvec(3,3),dmin(3),dmax(3),dmod,vol
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      parameter (small=1.0d-8)
      character unitlab*1,labunit*4
      logical unitb

      write (6,702)
  101 write (6,703,advance="no")
      read (5,*) unitlab 
      if (unitlab.eq.'B'.or.unitlab.eq.'b') then
        unitb = .true.
        labunit = 'Bohr'
      elseif (unitlab.eq.'A'.or.unitlab.eq.'a') then
        unitb = .false.
        labunit = 'Ang '
      else 
        write (6,*) ' Sorry, no third choice.' 
        goto 101
      endif
C ---- old scheme:
C     write (6,704,advance="no") labunit
C     read  (5,*) (obox(ii),ii=1,3)
C     write (6,705,advance="no") '1st',labunit
C     read  (5,*) (rbox(ii,1),ii=1,3)
C     write (6,705,advance="no") '2nd',labunit
C     read  (5,*) (rbox(ii,2),ii=1,3)
C     write (6,705,advance="no") '3rd',labunit
C     read  (5,*) (rbox(ii,3),ii=1,3)
C ---- new scheme:
      write (6,706,advance="no") labunit
      read  (5,*) bref
C
  201 write (6,707,advance="no") '1st'
      read  (5,*) dvec(:,1)
      dmod=sqrt(dvec(1,1)**2+dvec(2,1)**2+dvec(3,1)**2)
      if (dmod.lt.small) then
        write (6,*) ' Zero vector! Try again...'
        goto 201
      else 	
        dvec(:,1)=dvec(:,1)/dmod
      endif	
      write (6,708,advance="no") labunit
      read  (5,*) dmin(1),dmax(1)
C
  202 write (6,707,advance="no") '2nd'
      read  (5,*) dvec(:,2)
      dmod=sqrt(dvec(1,2)**2+dvec(2,2)**2+dvec(3,2)**2)
      if (dmod.lt.small) then
        write (6,*) ' Zero vector! Try again...'
        goto 202
      else 	
        dvec(:,2)=dvec(:,2)/dmod
      endif	
      write (6,708,advance="no") labunit
      read  (5,*) dmin(2),dmax(2)
C
  203 write (6,707,advance="no") '3rd'
      read  (5,*) dvec(:,3)
      dmod=sqrt(dvec(1,3)**2+dvec(2,3)**2+dvec(3,3)**2)
      if (dmod.lt.small) then
        write (6,*) ' Zero vector! Try again...'
        goto 203
      else 	
        dvec(:,3)=dvec(:,3)/dmod
      endif	
      write (6,708,advance="no") labunit
      read  (5,*) dmin(3),dmax(3)
C
C --- Origin of the box:
      obox(:)=bref(:)+dvec(:,1)*dmin(1) +	
     +                dvec(:,2)*dmin(2) +	
     +                dvec(:,3)*dmin(3) 	
C --- Spanning vectors:
      rbox(:,1) = dvec(:,1)*(dmax(1)-dmin(1))
      rbox(:,2) = dvec(:,2)*(dmax(2)-dmin(2))
      rbox(:,3) = dvec(:,3)*(dmax(3)-dmin(3))
C --- Check the volume:
      vol = rbox(1,1)*rbox(2,2)*rbox(3,3) +
     +      rbox(1,2)*rbox(2,3)*rbox(3,1) +
     +      rbox(1,3)*rbox(2,1)*rbox(3,2) -
     -      rbox(1,3)*rbox(2,2)*rbox(3,1) -
     -      rbox(1,2)*rbox(2,1)*rbox(3,3) -
     -      rbox(1,1)*rbox(2,3)*rbox(3,2) 
      write (6,709) labunit,vol
      if (abs(vol).lt.small) stop
      if (unitb) then
C       transform anything into Ang, as is standard in XCrysden:
        obox=obox*b2ang
        rbox=rbox*b2ang
      endif
      return

  702 format(" Now define the grid cell for your XCrysDen plot.",/
     .       " It is a parallelepiped, which can be arbitrarily",
     .       " chosen",/" with respect to the Siesta simulation cell,",
     .       " and needs not to be orthogonal.",/" We'll define it",
     .       " by the reference point,",/
     .       " three Cartesian direction vectors passing through it,",/
     .       " and limiting distances, measured along each",
     .       " vector from the reference point.",/" The reference",
     .       " point and distances can be given in Bohr or Ang;",/
     .       " the units of direction vectors are irrelevant.")
  703 format(' Would you use Bohr (B) or Ang (A) ? ')
  704 format(' Enter origin point in ',a4,' : ')
  705 format(' Enter ',a3,' spanning vector in ',a4,' : ')
  706 format(' Enter reference point in ',a4,' : ')
  707 format(' Enter ',a3,' direction vector in any units : ')
  708 format(' Enter min./max. limits along this vector,',/
     .       ' measured from reference point in ',a4,' : ')
  709 format(' Box volume in ',a4,' cube : ',f15.8)

      end
