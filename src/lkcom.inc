c LKCOM.INC
c
      integer nlkcolsize,nlkrowsize,nlkbufsize,nlkiter,
     &        ilkmin,ilkmax,ilkiniter,nlkterms
      logical lkleakage,lkincludesubgrid,lkcheck,lksolving,
     &        lkaveragehead,lkupdateresistances,lkexclude,
     &        lknosubgridfile
      real(8) rlkx0,rlky0,rlkrelax,rQ_leakage,rQ_target_leakage
      complex(8) ci
c
      common/lkcom/ nlkcolsize,nlkrowsize,nlkbufsize,nlkiter,
     &              lkleakage,lkincludesubgrid,lkcheck,lksolving,
     &              lkaveragehead,lkupdateresistances,lkexclude,
     &              lknosubgridfile,
     &              rlkx0,rlky0,rlkrelax,rQ_leakage,rQ_target_leakage,
     &              ci,ilkmin,ilkmax,ilkiniter,nlkterms
c
c
c     nlkcolsize              is the number of columns in the MODFLOW grid
c     nlkrowsize              is the number of rows in the MODFLOW grid
c     nlkbufsize              is the number of elements in the subgrid arrays
c     lkleakage               true if leakage grid is to be used
c     lkcheck                 true if a leakage comparison is to be calculated and written to .xtr file
c     lksolving               when true the change in leakage is used for omega, W, and normal flow functions.
c     lkaveragehead           true when upper heads are to be averaged over subgrid centers
c     lkupdateresistances     true if upper resistances need to be obtained from GFLOW and updated
c     lkexclude               true if MODFLOW cells with subgrids are to be excluded from the matrix
c     lknosubgridfile         true if normal *.grs file is found
c     nlkiter                 number of iterations on the leakage grid
c     lkincludesubgrid        true if subgrid is to be included in calculations
c     rlkx0,rlky0             coordinates of the upper left corner of the grid; row#=1,column#=1
c     rlkrelax                relaxation factor for leakage grid calculations
c     rQ_leakage              total leakage summed over a sub-set of MODFLOW cells
c     rQ_target_leakage       total leakage integrated (numerically) as calculated from upper head
c                                                      minus lower head divided by the resistance
c     ci                      (0.0,1.0)
c     ilkmin                  first equation for the leakage grid (to be used in lkmatcorrect_actual)
c     ilkmax                  last equation for the leakage grid (to be used in lkmatcorrect_actual)
c     ilkiniter               iteration at which to start the subgrid solution (if specified).
c     nlkterms                number of terms used in cell farfield expansion
c
c     allocatable arrays in LKIN:
c
c     rlkdeltax(j)            grid spacing in x-direction
c     rlkdeltay(i)            grid spacing in y-direction
c     rlkheadlower(i,j)       head in the lower aquifer (MODFLOW) in MODFLOW cell i,j
c     rlkleakage(i,j)         leakage into lower aquifer (MODFLOW) in MODFLOW cell i,j
c     rlkrecharge(i,j)        recharge into GFLOW aquifer (used to match MODFLOW recharge)
c     rlkdeltaleakage(i,j)    leakage increment between iterations
c     rlkresist(i,j)          resistance for MODFLOW cell i,j
c     ilkresolution(i,j)      number of sub-cells in MODFLOW cell i,j
c     rlksubheadlower(n)      head in lower aquifer (MODFLOW) for sub-cell n 1)
c     rlksubpotupper(n)       potential in upper aquifer (GFLOW) for sub-cell n
c     rlksubleakage(n)        leakage in sub-cell n minus the leakage in the associated MODFLOW cell 2)
c
c     notes:
c            1)  The sub-cells are numbered sequentially from MODFLOW cell to MODFLOW cell as follows:
c                first all sub-cells in MODFLOW cell (1,1), next all sub-cells in MODFLOW cell (2,1), etc.
c                The sub-cells within a MODFLOW cell are numbered differently than the MODFLOW cells in
c                the MODFLOW grid. The subcells are numbered from the bottom left to the right top of the MODFLOW cell.
c            2)  The leakage in sub-cells represent the deviation from the leakage in the associated MODFLOW
c                cell. Hence, the actual leakage is the leakage in the MODFLOW cell + the leakage in the
c                sub-cells.
c            3)  ilkresolution(i,j)=0  cell is to be excluded from the matrix (initial leakage value to be kept)
c                ilkresolution(i,j)=1  cell has NO subgrid, only one collocation point at the center.
c                ilkresolution(i,j)=2  cell has a 2*2 subgrid, hence 4 sub-cells.
c                ilkresolution(i,j)=3  cell has a 3*3 subgrid, hence 9 sub-cells.
c                       etc.
c

