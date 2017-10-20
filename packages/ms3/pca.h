/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/11/2016
*******************************************************/

#ifndef PCA_H
#define PCA_H

#include "mda.h"
#include "mda32.h"

// see info below
void pca(Mda& components, Mda& features, Mda& sigma, const Mda& X, bigint num_features, bool subtract_mean);
void pca(Mda32& components, Mda32& features, Mda32& sigma, const Mda32& X, bigint num_features, bool subtract_mean);
void pca_subsampled(Mda32& components, Mda32& features, Mda32& sigma, const Mda32& X, bigint num_features, bool subtract_mean, bigint max_samples);

// same as pca, except input it X*X', and features are not computed (because how could they be?)
void pca_from_XXt(Mda& components, Mda& sigma, const Mda& XXt, bigint num_features);
void pca_from_XXt(Mda32& components, Mda32& sigma, const Mda32& XXt, bigint num_features);

// get the whitening matrix as described below
void whitening_matrix_from_XXt(Mda& W, const Mda& XXt);
void whitening_matrix_from_XXt(Mda32& W, const Mda32& XXt);

void pca_unit_test();

/*
  pca: compute the top K principal components along with the corresponding feature vectors and
  singular values

  Input:
  X -- MxN matrix

  Output:
  C = components -- MxK, where K=num_features
  F = features -- KxN
  sigma = singular values -- Kx1

  components are normalized so that
  C'*C = eye(K,K)

  F = C'*X
  X is approximated by C*F = C*C'*X

  The components are eigenvectors of X*X':
  X*X'*C = C * diag(sigma)

  where sigma are the singular values associated with the top K components

  sigma(1)>=...>=sigma(K)

  if (K=M) then
  C'*C=C*C'=eye(M,M)
  X*X' = C * diag(sigma) * C' is the svd of X*X'

  The whitening matrix is:
  W = C * diag(sigma^(-1/2)) *C' = CDC'
  because
  (WX)*(WX)' = WXX'W' = CDC'XX'CDC' = C*D*diag(sigma)*D*C' = C*C' = eye(M,M)
*/

#endif // PCA_H
