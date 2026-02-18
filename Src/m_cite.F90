! ---
! Copyright (C) 1996-2025	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_cite
  implicit none

  private

  character(len=*), parameter :: STR_NULL = 'NULL'
  integer, parameter :: LEN_COMMENT = 256
  integer, parameter :: LEN_TYPE = 32
  integer, parameter :: LEN_DOI = 64
  integer, parameter :: LEN_CITEKEY = 64
  integer, parameter :: LEN_TITLE = 256
  integer, parameter :: LEN_AUTHOR = 512
  integer, parameter :: LEN_JOURNAL = 128
  integer, parameter :: LEN_VOLUME = 32
  integer, parameter :: LEN_ISSUE = 32
  integer, parameter :: LEN_PAGE = 32
  integer, parameter :: LEN_URL = 256

  type citation
    character(len=LEN_COMMENT) :: comment = STR_NULL
    character(len=LEN_TYPE) :: type = 'article'
    character(len=LEN_AUTHOR) :: author = STR_NULL
    character(len=LEN_TITLE) :: title = STR_NULL
    character(len=LEN_JOURNAL) :: journal = STR_NULL
    integer :: year = 0
    character(len=LEN_VOLUME) :: volume = STR_NULL
    character(len=LEN_ISSUE) :: issue = STR_NULL
    character(len=LEN_PAGE) :: page = STR_NULL
    character(len=LEN_URL) :: url = STR_NULL
    character(len=LEN_CITEKEY) :: cite_key = STR_NULL
    character(len=LEN_DOI) :: doi = STR_NULL
  end type citation

  ! The master citation database
  ! Split in two (or more) sections to avoid continuation-line limits
  ! This idiom needs F2008 support in the compiler
  type(citation), parameter :: citations1(*) = [ &
    citation( &
      comment = "Primary SIESTA paper", &
      author = "J. M. Soler and E. Artacho and J. D. Gale and A. Garcia and J. Junquera and P. Ordejon and D. Sanchez-Portal", &
      title = "The SIESTA method for ab initio order-N materials simulations", &
      journal = "Journal of Physics: Condensed Matter", &
      year = 2002, &
      volume = "14", &
      issue = "11", &
      page = "2745", &
      cite_key = "Soler2002", &
      doi = "10.1088/0953-8984/14/11/302"), &
    citation( &
      comment = "Primary TranSiesta paper", &
      author = "Brandbyge, M. and Mozos, J.-L. and Ordejón, P. and Taylor, J. and Stokbro, K.", &
      journal = "Physical Review B", &
      year = 2002, &
      volume = "65", &
      issue = "16", &
      cite_key = "Brandbyge2002", &
      doi = "10.1103/PhysRevB.65.165401"), &
    citation( &
      comment = "Slab-dipole correction", &
      journal = "Physical Review B", &
      year = 1999, &
      volume = "59", &
      issue = "16", &
      cite_key = "Bengtsson1999", &
      doi = "10.1103/PhysRevB.59.12301"), &
    citation( &
      comment = "LDA+U implementation", &
      journal = "Physical Review B", &
      year = 1998, &
      volume = "57", &
      issue = "3", &
      cite_key = "Dudarev1998", &
      doi = "10.1103/PhysRevB.57.1505"), &
    citation( &
      comment = "Charge/Hartree gate model", &
      journal = "Phys. Chem. Chem. Phys.", &
      year = 2016, &
      volume = "18", &
      issue = "2", &
      cite_key = "Papior2016A", &
      doi = "10.1039/C5CP04613K"), &
    citation( &
      comment = "TranSiesta N-electrode", &
      title = "Improvements on non-equilibrium and transport Green function techniques: The next-generation TranSiesta", &
      author = "Papior, N. and Lorente, N. and Frederiksen, T. and Garcia, A. and Brandbyge, M", &
      journal = "Computer Physics Communications", &
      year = 2017, &
      volume = "212", &
      cite_key = "Papior2017", &
      doi = "10.1016/j.cpc.2016.09.022"), &
    citation( &
      comment = "SIESTA-PEXSI", &
      title = "SIESTA-PEXSI: Massively parallel method for efficient and accurate &
      &ab initio materials simulation without matrix diagonalization", &
      journal = "Journal of Physics: Condensed Matter", &
      year = 2014, &
      volume = "26", &
      issue = "30", &
      cite_key = "Lin2014", &
      doi = "10.1088/0953-8984/26/30/305503"), &
    citation( &
      comment = "Real-space self-energies for TranSiesta", &
      author = "Papior, N. and Calogero, G. and Leitherer, S and Brandbyge, M", &
      journal = "Physical Review B", &
      year = 2019, &
      volume = "100", &
      issue = "19", &
      cite_key = "Papior2019A", &
      doi = "10.1103/PhysRevB.100.195417"), &
    citation( &
      comment = "Spin-orbit coupling (on-site approximation)", &
      title = "On-site approximation for spin–orbit coupling inlinear &
      &combination of atomic orbitals densityfunctional methods", &
      author = "Fernandez-Seivane, L. and Oliveira, M. A. and Sanvito, S. and Ferrer, J.", &
      journal = "Journal of Physics: Condensed Matter", &
      year = 2006, &
      volume = "18", &
      issue = "7999", &
      cite_key = "FernandezSeivane2006", &
      doi = "10.1088/0953-8984/19/19/489001"), &
    citation( &
      comment = "Spin-orbit coupling (off-site approximation)", &
      title = "Fully relativistic pseudopotential formalism under an atomic orbital basis: &
      &spin–orbit splittings and magnetic anisotropies", &
      author = "Cuadrado, R. and Cerda, J. I.", &
      journal = "Journal of Physics: Condensed Matter", &
      year = 2012, &
      volume = "24", &
      issue = "086005", &
      cite_key = "Cuadrado2012", &
      doi = "10.1088/0953-8984/24/8/086005"), &
    citation( &
      comment = "Maximally Localized Wannier Function paper", &
      author = "N. Marzari and D. Vanderbilt", &
      title = "Maximally localized generalized Wannier functions for composite energy bands", &
      journal = "Physical Review B", &
      year = 1997, &
      volume = "56", &
      issue = "20", &
      page = "12847", &
      cite_key = "Marzari1997", &
      doi = "10.1103/PhysRevB.56.12847"), &
    citation( &
      comment = "Wannier90 paper", &
      author = "A. A. Mostofi and J. R. Yates and G. Pizzi and Y.-S. Lee and I. Souza and D. Vanderbilt and N. Marzari", &
      title = "An updated version of wannier90: A tool for obtaining maximally-localised Wannier functions", &
      journal = "Comput. Phys. Commun.", &
      year = 2014, &
      volume = "185", &
      issue = "8", &
      page = "2309", &
      cite_key = "Wannier90", &
      doi = "10.1016/j.cpc.2014.05.003"), &
    citation( &
      comment = "Lowdin orthonormalization", &
      author = "Jorge Iniguez and Taner Yildirim", &
      title = "Unusual structural tuning of magnetism in cuprate perovskites", &
      journal = "Physical Review B", &
      year = 2005, &
      volume = "71", &
      issue = "18", &
      page = "180415", &
      cite_key = "Iniguez2005", &
      doi = "10.1103/PhysRevB.71.180415"), &
    citation( &
      comment = "Wannier90 paper, version 3.X", &
      author = "G. Pizzi et al.", &
      title = "Wannier90 as a community code: new features and applications", &
      journal = "J. Phys.: Condens. Matter", &
      year = 2020, &
      volume = "32", &
      issue = "16", &
      page = "165902", &
      cite_key = "Wannier90-v3.0.0", &
      doi = "10.1088/1361-648x/ab51ff"), &
    citation( &
      comment = "Review paper on Maximally Localized Wannier functions", &
      author = "N. Marzari et al.", &
      title = "Maximally localized Wannier functions: Theory and applications", &
      journal = "Rev. Mod. Phys.", &
      year = 2012, &
      volume = "84", &
      issue = "4", &
      page = "1419", &
      cite_key = "MLWF-review", &
      doi = "10.1103/RevModPhys.84.1419"), &
    citation( &
      comment = "ELSI library interface", &
      title = "ELSI: A unified software interface for Kohn–Sham electronic structure solvers", &
      author = "Victor Yu et al", &
      journal = "Computer Physics Communications", &
      year = 2018, &
      volume = "222", &
      issue = "267", &
      cite_key = "Yu2018", &
      doi = "10.1016/j.cpc.2017.09.007"), &
    citation( &
      comment = "SelInv algorithm", &
      title = "Fast algorithm for extracting the diagonal of the inverse matrix &
      &with application to the electronic structure analysis of metallic systems", &
      author = "Lin, L. and Lu, J. and Ying, L. and Car, R. and E, W.", &
      journal = "Comm. Math. Sci.", &
      year = 2009, &
      volume = "7", &
      issue = "755", &
      cite_key = "CMS2009", &
      doi = "projecteuclid.org/euclid.cms/12565628222"), &
    citation( &
      comment = "Bulk-bias calculations", &
      title = "Simple approach to current-induced bond weakening in ballistic conductors", &
      author = "Papior, N. and Leitherer, S. and Brandbyge, M.", &
      journal = "Physical Review B", &
      year = 2022, &
      volume = "106", &
      issue = "155401", &
      cite_key = "Papior22BulkBias", &
      doi = "10.1103/PhysRevB.106.155401") &
      ]

  ! second section
  type(citation), parameter :: citations2(*) = [ &
    citation( &
      comment = "ELPA parallel eigenvalue solutions", &
      author = "T. Auckenthaler and V. Blum and H.-J. Bungartz and &
      &T. Huckle and R. Johanni and L. Krämer and B. Lang and H. Lederer and P. R. Willems", &
      title = "Parallel solution of partial symmetric eigenvalue problems from electronic structure calculations", &
      journal = "Parallel Computing", &
      year = 2011, &
      volume = "37", &
      page = "783", &
      cite_key = "Auckenthaler2011", &
      doi = "10.1016/j.parco.2011.05.002"), &
    citation( &
      comment = "ELPA library", &
      author = "Andreas Marek and Volker Blum and Rainer Johanni and &
      &Ville Havu and Bruno Lang and Thomas Auckenthaler &
      &and Alexander Heinecke and Hans-Joachim Bungartz and Hermann Lederer", &
      title = "The ELPA library: scalable parallel eigenvalue solutions &
      &for electronic structure theory and computational science", &
      journal = "Journal of Physics: Condensed Matter", &
      year = 2014, &
      volume = "26", &
      page = "213201", &
      cite_key = "Marek2014", &
      doi = "10.1088/0953-8984/26/21/213201"), &
    citation( &
      comment = "ELPA initial GPU optimization", &
      type = "inproceedings", &
      author = "Pavel Kus and Andreas Marek and Hermann Lederer", &
      title = "GPU Optimization of Large-Scale Eigenvalue Solver", &
      journal = "Numerical Mathematics and Advanced Applications ENUMATH 2017", &
      year = 2017, &
      volume = "126", &
      cite_key = "Kus2017", &
      doi = "10.1007/978-3-319-96415-7"), &
    citation( &
      comment = "ELPA GPU acceleration", &
      author = "Victor When-zhe Yu and Jonathan Moussa and Pavel Kus &
      &and Andreas Marek and Peter Messmer and Mina Yoon and Hermann Lederer and Volker Blum", &
      title = "GPU-acceleration of the ELPA2 distributed eigensolver &
      &for dense symmetric and hermitian eigenpromlems", &
      journal = "Computer Physics Communications", &
      year = 2021, &
      volume = "262", &
      page = "107808", &
      cite_key = "Yu2021", &
      doi = "10.1016/j.cpc.2020.107808"), &
    citation( &
      comment = "ELPA eigenvalue-solver optimizations", &
      author = "P. Kus and A. Marek and S. S. Koecher and H.-H. Kowalski &
      &and Ch. Carbogno and Ch. Scheurer and K. Reuter and M. Scheffler and H. Lederer", &
      title = "Optimizations of the Eigenvaluesolvers in the ELPA Library", &
      journal = "Parallel Computing", &
      year = 2019, &
      volume = "85", &
      page = "167", &
      cite_key = "Kus2019", &
      doi = "10.1016/j.parco.2019.04.003"), &
    citation( &
      comment = "NTPoly library", &
      author = "William Dawson and Takahito Nakajima", &
      title = "Massively parallel sparse matrix function calculations with NTPoly", &
      journal = "Computer Physics Communications", &
      year = 2018, &
      volume = "225", &
      page = "154-165", &
      cite_key = "Dawson2018", &
      doi = "10.1016/j.cpc.2017.12.010"), &
    citation( &
      comment = "Simple-DFT-D3 library", &
      author = "Sebastian Ehlert", &
      title = "Simple DFT-D3: Library first implementation of the D3 dispersion correction", &
      journal = "Journal of Open Source Software", &
      year = 2024, &
      volume = "9", &
      issue = "103", &
      page = "7169", &
      cite_key = "Ehlert2024", &
      doi = "10.21105/joss.07169"), &
    citation( &
      comment = "DFTD3 BJ damping", &
      title = "Effect of the damping function in dispersion corrected density functional theory", &
      author = "Grimme, S. and  Ehrlich, S. and Goerigk, L.", &
      journal = "J. Comput. Chem.", &
      year = 2011, &
      volume = "32", &
      issue = "7", &
      cite_key = "DFTD3BJ", &
      doi = "10.1002/jcc.21759"), &
    citation( &
      comment = "DFTD3 Zero damping", &
      title = "A consistent and accurate ab initio parametrization &
      &of density functional dispersion correction (DFT-D) for the 94 elements H-Pu", &
      author = "Grimme, S. and Antony, J. and  Ehrlich, S. and Krieg, H.", &
      journal = "J. Chem. Phys.", &
      year = 2010, &
      volume = "132", &
      issue = "15", &
      cite_key = "DFTD3Zero", &
      doi = "10.1063/1.3382344") , &
    citation( &
      comment = "PSML pseudopotential format", &
      author = "Alberto Garcia and Matthieu J. Verstraete and Yann Pouillon and Javier Junquera", &
      title = "The psml format and library for norm-conserving pseudopotential data curation &
      &and interoperability", &
      journal = "Computer Physics Communications", &
      year = 2018, &
      volume = "227", &
      page = "51-71", &
      cite_key = "Garcia2018", &
      doi = "10.1016/j.cpc.2018.02.011"), &
    citation( &
      comment = "LIBXC library", &
      author = "Susi Lehtola and Conrad Steigemann and Micael J.T. Oliveira and Miguel A.L. Marques", &
      title = "Recent developments in libxc — A comprehensive library of functionals &
      &for density functional theory", &
      journal = "SoftwareX", &
      year = 2018, &
      volume = "7", &
      issue = "", &
      page = "1-5", &
      cite_key = "Lehtola2018", &
      doi = "10.1016/j.softx.2017.11.002") &
  ]

  ! Concatenation of the sections
  type(citation), parameter :: citations(*) = [citations1, citations2]

  ! To track which citations have been used in this run
  logical, allocatable :: used_citations(:)

  ! Default file-name
  character(len=64) :: cite_file = "CITATIONS.bib"

  public :: init_citation, add_citation, announce_citations
  public :: reset_citations

contains

  subroutine init_citation(prefix, delete)
    character(len=*), intent(in)  :: prefix
    logical, intent(in), optional :: delete
    logical :: ldelete
    integer :: iu

    ldelete = .true.
    if (present(delete)) ldelete = delete

    ! Initialize the used_citations array
    if (.not. allocated(used_citations)) then
      allocate(used_citations(size(citations)))
      used_citations = .false.
    end if

    cite_file = trim(prefix) // ".bib"
    if (.not. ldelete) return

    ! Be sure to have an "empty" citation file

    open(newunit=iu, file=cite_file, form='formatted', position='APPEND')
    close(iu, status='DELETE')

  end subroutine init_citation

  ! This is the main callable routine. It adds the citation to the 'used' list
  ! It can be called by any processor, not just the IONode, since it does not
  ! do any I/O
  subroutine add_citation(doi)
    character(len=*), intent(in) :: doi
    integer :: i

    ! Find the citation with matching DOI
    do i = 1, size(citations)
      if (citations(i)%doi == doi) then
        if (.not. used_citations(i)) then
          used_citations(i) = .true.
        endif
        RETURN
      end if
    end do

  end subroutine add_citation

  ! Routine that can be called at the end of the run (or, actually,
  ! several times along the run, if needed, since the .bib file is
  ! rewound).

  subroutine announce_citations(version_str)
    character(len=*), intent(in) :: version_str
    integer :: i, iu

    write(*,'(/,2a)') 'cite: Please indicate the Siesta version in published work: ',trim(version_str)
    write(*,'(a)') 'cite: This calculation has made use of features in the following articles'
    write(*,'(a)') 'cite: which we encourage you to cite in published work.'
    write(*,'(3a)') 'cite: (Please see "',trim(cite_file),'" for a BibTeX file.)'

    open(newunit=iu, file=cite_file, form='formatted', position='REWIND', &
         action='WRITE')

    write(iu,'(a)') '# This file contains articles that we encourage you to cite'
    write(iu,'(a)') '# in published work.'
    write(iu,'(a)') '# Each entry corresponds to a feature that has'
    write(iu,'(a)') '# been enabled via FDF-flags and which is described in'
    write(iu,'(a)') '# a paper or other reference, as indicated.'
    write(iu,*)

    do i = 1, size(citations)
       if (used_citations(i)) then
          ! Write summary in the output file
          write(*,'(tr8,a,/,tr10,2a)') trim(citations(i)%comment), &
               'DOI: www.doi.org/', trim(citations(i)%doi)

          ! Write full information in the .bib file
          call write_single_citation(citations(i), iu)
       end if
    end do
    close(iu)
    write(*,*)

  end subroutine announce_citations

  ! Helper routine
  subroutine write_single_citation(cit,iu)
    type(citation), intent(in) :: cit
    integer, intent(in)        :: iu

    if (cit%comment /= STR_NULL) then
      write(iu,'(2a)') '# ', trim(cit%comment)
    end if

    write(iu,'(5a)') '@',trim(cit%type),'{',trim(cit%cite_key),','

    if (cit%author /= STR_NULL) then
      write(iu,'(t3,3a)') 'author = {{',trim(cit%author),'}},'
    end if
    if (cit%title /= STR_NULL) then
      write(iu,'(t3,3a)') 'title = {{',trim(cit%title),'}},'
    end if
    if (cit%journal /= STR_NULL) then
      write(iu,'(t3,3a)') 'journal = {{',trim(cit%journal),'}},'
    end if
    if (cit%year /= 0) then
      write(iu,'(t3,a,i0,a)') 'year = {',cit%year,'},'
    end if
    if (cit%volume /= STR_NULL) then
      write(iu,'(t3,3a)') 'volume = {',trim(cit%volume),'},'
    end if
    if (cit%issue /= STR_NULL) then
      write(iu,'(t3,3a)') 'issue = {',trim(cit%issue),'},'
    end if
    if (cit%page /= STR_NULL) then
      write(iu,'(t3,3a)') 'page = {',trim(cit%page),'},'
    end if
    if (cit%url /= STR_NULL) then
      write(iu,'(t3,3a)') 'url = {',trim(cit%url),'},'
    end if
    if (cit%doi /= STR_NULL) then
      write(iu,'(t3,3a)') 'doi = {',trim(cit%doi),'},'
    end if

    write(iu,'(a)') '}'
    write(iu,*)

  end subroutine write_single_citation

  subroutine reset_citations()
    implicit none
    if (allocated(used_citations)) deallocate(used_citations)
  end subroutine reset_citations
end module m_cite
