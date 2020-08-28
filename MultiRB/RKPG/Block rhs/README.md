#### Block right-hand side figures 

## Independedent tests
- b_Ab_n500_bones: n=500, rhs is given by b for 1 column and [b, Ab] for 2 columns, where b is a vector of ones, 8 Zolotarev poles
- dif_size_8roots_denom_2colalt: n varies for tests with 8 Zolotarev poles, rhs is 2 cols of alternate 1 & -1 values (one where top half is alternating 1&-1 and bottom is zeros, and the other where top half is zeros and bottom half is alternating 1&-1)
- different_nr_poles: n=100, rhs is 2 columns, 1st ones, 2nd alternate 1&-1, number of Zolotarev poles vary
- different_nr_poles_3colrhs: n=100, rhs is 3 columns, 1st ones, 2nd alternate 1&-1, 3rd top half ones, bottom half zeros, number of Zolotarev poles vary
- n100_6poles_differentrhs: n=100, 3 tests where rhs varies (1 col: ones, 2 cols: [ones, alternate 1&-1], 3 cols:[ones, alternate 1&-1, top half ones bottom zero]), 6 Zolotarev poles
- n500-identity: n=500, rhs is 1, 2, 3 columns respectively (1 col: first col of identity, 2 col: [1st col I, 2nd col I], 3 col: [1st col I, 2nd col I, 3rd col I]), 8 Zolotarev poles
- randorthrhs_n500: n=500, rhs is 1, 2, 3 columns respectively, all random and built on top of eachother to be orthogonal, 8 Zolotarev poles


## Build on random rhs -- summary of below figs in comp_buildonrandrhs_nrpoles
- buildrand_6poles: n=500, rhs is in turn 1, 2, 3 columns of random values (the first column is kept the same for the 2 col test, and the 1st & 2nd columns are kept the same for the 3 col test), 6 Zolotarev poles
- buildrand_8poles: n=500, rhs is in turn 1, 2, 3 columns of random values (the first column is kept the same for the 2 col test, and the 1st & 2nd columns are kept the same for the 3 col test), 8 Zolotarev poles (same rhs as in prev fig)
- buildrand_10poles: n=500, rhs is in turn 1, 2, 3 columns of random values (the first column is kept the same for the 2 col test, and the 1st & 2nd columns are kept the same for the 3 col test), 10 Zolotarev poles (same rhs as in prev figs)
- buildrand_12poles: n=500, rhs is in turn 1, 2, 3 columns of random values (the first column is kept the same for the 2 col test, and the 1st & 2nd columns are kept the same for the 3 col test), 12 Zolotarev poles (same rhs as in prev figs)

## Alternate 1&-1 2 columns tests -- summary of below figs in comp_alternate_nrpoles
- 6poles_alternate: n=500, rhs = alternate 1&-1, residual plot with 1 and 2 columns, 6 Zolotarev poles
- 8poles_alternate: n=500, rhs = alternate 1&-1, residual plot with 1 and 2 columns, 8 Zolotarev poles
- 10poles_alternate: n=500, rhs = alternate 1&-1, residual plot with 1 and 2 columns, 10 Zolotarev poles
- 12poles_alternate: n=500, rhs = alternate 1&-1, residual plot with 1 and 2 columns, 12 Zolotarev poles

## Including Beckermann bounds:
- counts iterations as we add a new column (regardless of the number of columns in the rhs)
    - beckermann_n100_randrhs3_2poles: n=100, 3 column random rhs, 2 Zolotarev poles
    - beckermann_n100_randrhs3_8poles: n=100, 3 column random rhs, 8 Zolotarev poles

- counts iterations as we add (the number of columns in the rhs) at each step
    - beckermann_rkpg2_n100_rhs1-1_6poles: n=100, rhs is given by one column of ones and one column of alternating 1&-1, 6 Zolotarev poles
    - beckermann_rkpg2_n500_rhs_half10-01_8poles: n=500, rhs has 2 columns, one where top half is alternating 1&-1 and bottom is zeros, and the other where top half is zeros and bottom half is alternating 1&-1, 8 Zolotarev poles
    - n500_alt0-0alt_8rootsdenom_beckermann: n=500, rhs 1 and 2 columns (1 col: alternate 1&-1, 2 col: one where top half is alternating 1&-1 and bottom is zeros, and the other where top half is zeros and bottom half is alternating 1&-1), 8 Zolotarev poles
    - n500_rhsI_8rootsden: n=500, rhs is 1, 2, 3 columns respectively (1 col: first column of identity matrix, 2 col: [1 and zeros, 1 1 and zeros], 3 col: [1 and zeros, 1 1 and zeros, 1 1 1 and zeros]), 8 Zolotarev poles