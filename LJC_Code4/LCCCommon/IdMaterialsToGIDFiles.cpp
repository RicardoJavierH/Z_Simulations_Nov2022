#include "IdMaterialsToGIDFiles.h"

// Inserting material id to elements from GID mesh
void SetMaterialIdForElements(TPZGeoMesh* gmesh, int matid) {
	int DimProblem = gmesh->Dimension();
	if (DimProblem < 0 || DimProblem > 3)
		return;
	int nel = gmesh->NElements();
	for (int iel = 0; iel < nel; iel++) {
		TPZGeoEl* gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem) continue;
		gel->SetMaterialId(matid);
	}///iel
}
// Inserting material id to elements from GID mesh
void SetMatBCIdForElements(TPZGeoMesh * gmesh, int matid, TPZVec<REAL> & P0, TPZVec<REAL> & P1) {
	int DimProblem = gmesh->Dimension();
	if (DimProblem < 0 || DimProblem > 3)
		return;
	int count = 0, nel = gmesh->NElements();
	TPZManVector<REAL> Center(3, 0.), Point(3, 0.);
	REAL t = -1., r = -1., p = -1.;
	for (int iel = 0; iel < nel; iel++) {
		TPZGeoEl* gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem - 1) continue;
		gel->CenterPoint(gel->NSides() - 1, Center);
		gel->X(Center, Point);
		if (IsZero(P1[0] - P0[0])) {
			if (!IsZero(Point[0] - P0[0])) continue;
			t = 0.;
		}
		else
			t = (Point[0] - P0[0]) / (P1[0] - P0[0]);
		if (IsZero(P1[1] - P0[1])) {
			if (!IsZero(Point[1] - P0[1])) continue;
			r = 0.;
		}
		else
			r = (Point[1] - P0[1]) / (P1[1] - P0[1]);
		if (IsZero(P1[2] - P0[2])) {
			if (!IsZero(Point[2] - P0[2])) continue;
			p = 0.;
		}
		else
			p = (Point[2] - P0[2]) / (P1[2] - P0[2]);

		if (r < 0. || t < 0. || p < 0.) continue;
		if (r > 1. || t > 1. || p > 1.) continue;
		if ((IsZero(t) && IsZero(r)) || (IsZero(t) && IsZero(p)) || (IsZero(r) && IsZero(p))) {
			gel->SetMaterialId(matid); count++;
			continue;
		}
		if ((IsZero(t) && IsZero(r - p)) || (IsZero(r) && IsZero(t - p)) || (IsZero(p) && IsZero(r - t))) {
			gel->SetMaterialId(matid);
			count++;
		}
	}///iel
	std::cout << std::endl << "BC Condition " << matid << " com " << count << " elements." << std::endl;
}
// Inserting material id to elements from GID mesh
// Taptype: t - Top, b - Bottom, l - Left, r - Right, f - Front, v - Verso
void SetMatBCIdForElementsSameByHeigh(TPZGeoMesh * gmesh, int matid, char taptype) {
	int DimProblem = gmesh->Dimension();
	if (DimProblem < 0 || DimProblem > 3)
		return;
	int64_t iel, count = 0, nel = gmesh->NElements();
	TPZManVector<REAL> Center(DimProblem, 0.), Point(3, 0.);
	REAL Tap;
	switch (taptype) {
	case 'b':
	case 'l':
	case 'v':
		Tap = 1000000.;
		break;
	case 't':
	case 'r':
	case 'f':
		Tap = -1000000.;
		break;
	default:
		Tap = 0.;
	}
	for (iel = 0; iel < nel; iel++) {
		TPZGeoEl* gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		// Determining bc IdMat into elements with codimension 1
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem) continue;
		gel->CenterPoint(gel->NSides() - 1, Center);
		gel->X(Center, Point);
		switch (taptype) {
		case 'b':
			if (Tap > Point[2])
				Tap = Point[2];
			break;
		case 't':
			if (Tap < Point[2])
				Tap = Point[2];
			break;
		case 'l':
			if (Tap > Point[1])
				Tap = Point[1];
			break;
		case 'r':
			if (Tap < Point[1])
				Tap = Point[1];
			break;
		case 'f':
			if (Tap < Point[0])
				Tap = Point[0];
			break;
		case 'v':
			if (Tap > Point[0])
				Tap = Point[0];
			break;
		default:
			Tap = 0.;
		}
	}
	Center.Resize(DimProblem - 1);
	for (iel = 0; iel < nel; iel++) {
		TPZGeoEl* gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		// Determining bc IdMat into elements with codimension 1
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem - 1) continue;
		gel->CenterPoint(gel->NSides() - 1, Center);
		gel->X(Center, Point);
		switch (taptype) {
		case 'b':
		case 't':
			if (IsZero(Tap - Point[2])) {
				gel->SetMaterialId(matid); count++;
			};
			break;
		case 'v':
		case 'f':
			if (IsZero(Tap - Point[0])) {
				gel->SetMaterialId(matid); count++;
			};
			break;
		case 'r':
		case 'l':
			if (IsZero(Tap - Point[1])) {
				gel->SetMaterialId(matid); count++;
			};
			break;
		default:
			gel->SetMaterialId(matid); count++;
			break;
		}
	}///iel
	std::cout << std::endl << "BC Condition " << matid << " com " << count << " elements (" << taptype << ")." << std::endl;
}
// Inserting material id to elements from GID mesh
void SetMatBCIdForElementsOnOneCoordinate(TPZGeoMesh * gmesh, int matid, int coord, REAL coordfixed, TPZVec<REAL> & P0, TPZVec<REAL> & P1) {
	int DimProblem = gmesh->Dimension();
	if (DimProblem < 0 || DimProblem > 3)
		return;
	if (coord < 0 || coord > 2) return;
	int coordv = (coord + 1) % 3;
	if (!IsZero(P0[coordv] - coordfixed))
		coordv = (coordv + 1) % 3;

	int count = 0, nel = gmesh->NElements();
	TPZManVector<REAL> Center(3, 0.), Point(3, 0.);
	REAL t;
	for (int iel = 0; iel < nel; iel++) {
		t = -1.;
		TPZGeoEl* gel = gmesh->ElementVec()[iel];
		if (!gel) continue;
		if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
		if (gel->Dimension() != DimProblem - 1) continue;
		if (gel->MaterialId() > -1) continue;
		gel->CenterPoint(gel->NSides() - 1, Center);
		gel->X(Center, Point);
		if (!IsZero(Point[coordv] - coordfixed))
			continue;
		if (IsZero(P1[coord] - P0[coord])) {
			t = 0.;
		}
		else
			t = (Point[coord] - P0[coord]) / (P1[coord] - P0[coord]);

		if (t < 0. || t > 1.) continue;
		gel->SetMaterialId(matid);
		count++;
		continue;
	}///iel
	std::cout << std::endl << "BC Condition " << matid << " com " << count << " elements." << std::endl;
}
