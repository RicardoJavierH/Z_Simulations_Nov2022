#ifndef IDMATERIALSTOGIDFILESHH
#define IDMATERIALSTOGIDFILESHH

#include "pzvec.h"

#include "TPZReadGIDGrid.h"

#include "pzgmesh.h"

// Utilitaries

// Utilitaries
void SetMaterialIdForElements(TPZGeoMesh* gmesh, int matid);
void SetMatBCIdForElements(TPZGeoMesh* gmesh, int matid, TPZVec<REAL>& P0, TPZVec<REAL>& P1);

// Inserting material id to elements from GID mesh
// Taptype: t - Top, b - Bottom, l - Left, r - Right, f - Front, v - Verso
void SetMatBCIdForElementsSameByHeigh(TPZGeoMesh* gmesh, int matid, char taptype);
// inserting material id to elements 3D without heigh
void SetMatBCIdForElementsOnOneCoordinate(TPZGeoMesh* gmesh, int matid, int coord, REAL coordfixed, TPZVec<REAL>& P0, TPZVec<REAL>& P1);

#endif
