#include "private/matimpl.h"
#include "private/vecimpl.h"  
/*
  This functions to be added to petsc's next release
*/

#undef __FUNCT__  
#define __FUNCT__ "MatZeroRowsStencil"
/*@C
   MatZeroRowsStencil - Zeros all entries (except possibly the main diagonal)
   of a set of rows of a matrix.

   Collective on Mat

   Input Parameters:
+  mat - the matrix
.  numRows - the number of rows to remove
.  rows - the grid coordinates (and component number when dof > 1) for matrix rows
-  diag - value put in all diagonals of eliminated rows (0.0 will even eliminate diagonal entry)

   Notes:
   For the AIJ and BAIJ matrix formats this removes the old nonzero structure,
   but does not release memory.  For the dense and block diagonal
   formats this does not alter the nonzero structure.

   If the option MatSetOption(mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE) the nonzero structure
   of the matrix is not changed (even for AIJ and BAIJ matrices) the values are
   merely zeroed.

   The user can set a value in the diagonal entry (or for the AIJ and
   row formats can optionally remove the main diagonal entry from the
   nonzero structure as well, by passing 0.0 as the final argument).

   The grid coordinates are across the entire grid, not just the local portion

   In Fortran idxm and idxn should be declared as
$     MatStencil idxm(4,m)
   and the values inserted using
$    idxm(MatStencil_i,1) = i
$    idxm(MatStencil_j,1) = j
$    idxm(MatStencil_k,1) = k
$    idxm(MatStencil_c,1) = c
   etc

   For periodic boundary conditions use negative indices for values to the left (below 0; that are to be 
   obtained by wrapping values from right edge). For values to the right of the last entry using that index plus one
   etc to obtain values that obtained by wrapping the values from the left edge. This does not work for the DA_NONPERIODIC
   wrap.

   For indices that don't mean anything for your case (like the k index when working in 2d) or the c index when you have
   a single value per point) you can skip filling those indices.

   Level: intermediate

   Concepts: matrices^zeroing rows

.seealso: MatZeroRows(), MatZeroRowsIS(), MatZeroEntries(), MatZeroRowsLocal(), MatSetOption()
@*/
PetscErrorCode PETSCMAT_DLLEXPORT MatZeroRowsStencil(Mat mat,PetscInt numRows,const MatStencil rows[],PetscScalar diag)
{
  PetscInt       dim    = mat->stencil.dim;
  PetscInt       sdim   = dim - (1 - (PetscInt) mat->stencil.noc);
  PetscInt      *dims   = mat->stencil.dims+1;
  PetscInt      *starts = mat->stencil.starts;
  PetscInt      *dxm    = (PetscInt *) rows;
  PetscInt      *jdxm,i,j,tmp,numNewRows = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mat,MAT_COOKIE,1);
  PetscValidType(mat,1);
  if (numRows) PetscValidIntPointer(rows,3);

  ierr = PetscMalloc(numRows * sizeof(PetscInt),&jdxm);CHKERRQ(ierr);
  for(i = 0; i < numRows; ++i) {
    /* Skip unused dimensions (they are ordered k, j, i, c) */
    for(j = 0; j < 3-sdim; ++j) dxm++;
    /* Local index in X dir */
    tmp = *dxm++ - starts[0];
    /* Loop over remaining dimensions */
    for(j = 0; j < dim-1; ++j) {
      /* If nonlocal, set index to be negative */
      if ((*dxm++ - starts[j+1]) < 0 || tmp < 0) tmp = PETSC_MIN_INT;
      /* Update local index */
      else                                       tmp = tmp*dims[j] + *(dxm-1) - starts[j+1];
    }
    /* Skip component slot if necessary */
    if (mat->stencil.noc) dxm++;
    /* Local row number */
    if (tmp >= 0) {
      jdxm[numNewRows++] = tmp;
    }
  }
  ierr = MatZeroRowsLocal(mat,numNewRows,jdxm,diag);CHKERRQ(ierr);
  ierr = PetscFree(jdxm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

