/**
 * @file
 * @brief Contains the implementation of the TPZReadGIDGrid methods. 
 */

#include "TPZReadGIDGrid.h"
#include "pzgmesh.h"
#include "pzgeoelside.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "tpzcube.h"
#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include <tpzarc3d.h>

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "TPZVTKGeoMesh.h"
#include <pzgengrid.h>

#include <iostream>
#include <fstream>

TPZReadGIDGrid::TPZReadGIDGrid() {
	MatNumber = 0;
	BCNumber = 0;
	fProblemDimension = 0;
	fDimensionlessL = 1.0;
}//method

TPZReadGIDGrid::~TPZReadGIDGrid() {
	
}//method

TPZGeoMesh * TPZReadGIDGrid::GeometricGIDMesh(std::string FiletoRead)
{
	// File to read Dump file generated by GID
	
	std::string FileName;
	std::string stringTemp;
	FileName = FiletoRead;
	
	// Definitions
	int nMats = 0;
	int64_t numnodes=0;
	int64_t numelements=0;
	int64_t elements3DT=0;
	int64_t elements3DH=0;	
	int64_t elements2DT=0;
	int64_t elements2DQ=0;	
	int64_t elements1D=0;
	int64_t elements0D=0;
	
	//	Scanning for total Number of Nodes and differents Dimension Elements
	int64_t NumEntitiestoRead;
	TPZStack <std::string> SentinelString;
	{
		
		// reading a general mesh information by filter
		std::ifstream read (FileName.c_str());
		if (!read.is_open()) return 0;
		std::string FlagString;
		int flag = -1;		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			flag = str.find("---");
			
			if (flag >= 0)
			{
				if(str != "--- ELEMENTS ---" && str != "--- ELEMENTS ---\r") 
				{
					SentinelString.push_back(str);
				}
			}
			if(str == "LINEAR" || str == "LINEAR\r" ) 
			{
				SentinelString.push_back(str);
			}
			if(str == "TRIANGLE" || str == "TRIANGLE\r") 
			{
				SentinelString.push_back(str);
			}
			if(str == "QUADRILATERAL" || str == "QUADRILATERAL\r") 
			{
				SentinelString.push_back(str);
			}
			if(str == "TETRAHEDRA" || str == "TETRAHEDRA\r") 
			{
				SentinelString.push_back(str);
			}
			if(str == "HEXAHEDRA" || str == "HEXAHEDRA\r") 
			{
				SentinelString.push_back(str);
			}			
			
		}
		
		FlagString = "EndReading";
		SentinelString.push_back(FlagString);
	}
	
	NumEntitiestoRead = SentinelString.size();
	TPZStack <int> GeneralData(NumEntitiestoRead,0);
	TPZStack <int> DataToProcess(NumEntitiestoRead,-1);		
	
	{		
		// reading a general mesh information by filter
		std::ifstream read (FileName.c_str());
		std::string FlagString;
		int64_t cont = 0;

		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			
			// Reading General Data
			if(str == SentinelString[cont]) 
			{
				FlagString = str;	
			}
			if(SentinelString[cont] == "" || SentinelString[cont] == "\r") 
			{
				cont++;	
			}			
			
			if(SentinelString[cont] == "EndReading") 
			{
				break;	
			}			
			
			if( (str != "" || str != "\r") && FlagString == SentinelString[cont]) 
			{
				// Data scaning
				while (read) {
					char buftemp[1024];
					read.getline(buftemp, 1024);
					std::string strtemp(buftemp);
					GeneralData[cont]++;
					if(strtemp == "" || strtemp == "\r") 
					{
						FlagString = "";
						GeneralData[cont]--;
						cont++;
//						std::cout << "Scanning General Data -> done!" << std::endl;
						break;
					}
				}
				
			}	
			
			
		}
	}
	
	for (int64_t i = 0 ; i < NumEntitiestoRead; i++ )
	{
		if(SentinelString[i] == "--- USED MATERIALS ---" || SentinelString[i] == "--- USED MATERIALS ---\r") 
		{
			nMats=GeneralData[i];
			DataToProcess[i]=0;
		}
		if(SentinelString[i] == "--- CONDITIONS OVER NODES ---" || SentinelString[i] == "--- CONDITIONS OVER NODES ---\r") 
		{
			if(GeneralData[i] !=0)
			{
				GeneralData[i]--;	
				elements0D=GeneralData[i];
			}
			else 
			{
				
			}
			
			DataToProcess[i]=1;			
		}	
		if(SentinelString[i] == "--- NODES ---" || SentinelString[i] == "--- NODES ---\r") 
		{
			numnodes=GeneralData[i];
			DataToProcess[i]=2;			
		}		
		if(SentinelString[i] == "LINEAR" || SentinelString[i] == "LINEAR\r") 
		{
			elements1D=GeneralData[i];
			DataToProcess[i]=3;		
		}
		if(SentinelString[i] == "TRIANGLE" || SentinelString[i] == "TRIANGLE\r") 
		{
			elements2DT=GeneralData[i];
			DataToProcess[i]=4;		
		}
		if(SentinelString[i] == "QUADRILATERAL" || SentinelString[i] == "QUADRILATERAL\r") 
		{
			elements2DQ=GeneralData[i];
			DataToProcess[i]=5;			
		}
		if(SentinelString[i] == "TETRAHEDRA" || SentinelString[i] == "TETRAHEDRA\r") 
		{
			elements3DT=GeneralData[i];
			DataToProcess[i]=6;			
		}
		if(SentinelString[i] == "HEXAHEDRA" || SentinelString[i] == "HEXAHEDRA\r") 
		{
			elements3DH=GeneralData[i];
			DataToProcess[i]=7;			
		}		
		
	}

	numelements=elements3DT+elements3DH+elements2DT+elements2DQ+elements1D+elements0D;

	//  Mesh Creation
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;	
	gmesh -> NodeVec().Resize(numnodes);
	// needed for node insertion
	const int64_t Tnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Tnodes);
	int64_t nodeId;
	double nodecoordX , nodecoordY , nodecoordZ ;
	int64_t elementId = 0;
	int matElId = 0;
	int Layerid = 0;
	int ContMats = 0;
	int64_t ContNode = 0;
	int64_t ContPoint = 0;
	int64_t ContLine = 0;
	int64_t ContTrian = 0;
	int64_t ContQuad = 0;
	int64_t ContTet = 0;
	int64_t ContHex = 0;	
	
	TPZManVector <int64_t> TopolPoint(1);	
	TPZManVector <int64_t> TopolLine(2);
	TPZManVector <int64_t> TopolTriangle(3);
	TPZManVector <int64_t> TopolQuad(4);
	TPZManVector <int64_t> TopolTet(4);
	TPZManVector <int64_t> TopolHex(8);	
	
	{
		
		// reading a general mesh information by filter
		std::ifstream read (FileName.c_str());
		std::string FlagString;
		int64_t cont = 0;
		int dim = 0;
		int flag = 0;
		while(read)
		{
			char buf[1024];	
			read.getline(buf, 1024);
			std::string str(buf);
			std::string strtemp="InitialState";			
			
			// Reading General Data
			if(str == SentinelString[cont]) 
			{
				FlagString = str;	
			}
			
			if(SentinelString[cont] == "" || SentinelString[cont] == "\r") 
			{
				cont++;	
			}			
			
			if(SentinelString[cont] == "EndReading") 
			{
				break;	
			}			
			
			if( (str != "" || str != "\r" )&& FlagString == SentinelString[cont]) 
			{
				// Data scaning
				while (read) {
					
					switch (DataToProcess[cont]) {
						case 0:
						{
							//"--- USED MATERIALS ---"
							if (GeneralData[cont] != 0)
							{
								MaterialDataV MatTemp;
								read.getline(buf, 1024);
								int64_t spacecont = 0;
								char *p=buf, *q;
								while (p) 
								{
									q = strchr(p,' ');
									if(!q)
									{
										if (p) 
										{
											MatTemp.fProperties.push_back(atof(p));
										}
										break;
									}
									*q = 0;
									if(spacecont==0)
									{ 
										MatTemp.fMatID = atoi(p);
									}									
									else if(spacecont==1)
									{ 
										MatTemp.fMaterial = p;
										if(*p=='D' || *p=='d') 
										{
											MatNumber++;
										}
										else
										{
											BCNumber++;
										}

									}
									else 
									{		
										MatTemp.fProperties.push_back(atof(p));
										
									}
									
									p = q+1;							
									while(q && *p==' ')
										p++;
									spacecont++;
								}
								
								
								if (MatTemp.fMaterial[0]=='D' || MatTemp.fMaterial[0]=='d') 
								{
									fMaterialDataVec.Push(MatTemp);
								}
								else
								{
									fBCMaterialDataVec.Push(MatTemp);
								}

								ContMats++;
							}
							if(ContMats == nMats)
							{
								strtemp = "";	
							}					
						}
							break;
						case 1:
						{
							//"--- CONDITIONS OVER NODES ---"
							// 0D Elements
							if (GeneralData[cont] != 0)
							{
								if(flag == 0)
								{
									read.getline(buf, 1024);
									flag++;
								}
								int MatID;
								read >> TopolPoint[0]; //node 1	
								read >> MatID; //Material identity									
								read.getline(buf, 1024);
								TopolPoint[0]--;
								ContPoint++;						
								new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (numelements - elements0D + ContPoint, TopolPoint, MatID,*gmesh); 								
								
							}
							if(ContPoint == elements0D)
							{
								strtemp = "";	
							}					
						}
							break;
						case 2:
						{
							if (GeneralData[cont] != 0)
							{
								//"--- NODES ---"
								
								if (!dim)
								{							
									read.getline(buf, 1024);
									int64_t spacecont = 0;
									char *p=buf, *q;
									while (p) {
										q = strchr(p,' ');
										if(!q)
										{
											break;
										}
										*q = 0;
										if(spacecont==0) nodeId = atoi(p);
										else if(spacecont==1) Layerid = atoi(p);
										else if(spacecont==2) nodecoordX = atof(p);
										else if(spacecont==3) nodecoordY = atof(p);
										p = q+1;							
										while(q && *p==' ')
											p++;
										spacecont++;
									}
									if (spacecont == 3)
									{
										dim = 2;
										fProblemDimension = dim;
										nodecoordY = atof(p);							
										nodecoordZ = 0.0;							
									}
									else 
									{
										dim = 3;
										fProblemDimension = dim;
										nodecoordZ = atof(p);				
									}
									
								}
								else 
								{
									read >> nodeId;
									read >> Layerid;						
									read >> nodecoordX;
									read >> nodecoordY;
									if (dim == 2) {
										nodecoordZ = 0.0;
									}		
									else 
									{
										read >> nodecoordZ;
									}
								}
								int64_t nodeid = gmesh->CreateUniqueNodeId();
								Node[nodeId-1].SetNodeId(nodeid);
								Node[nodeId-1].SetCoord(0,nodecoordX/fDimensionlessL);
								Node[nodeId-1].SetCoord(1,nodecoordY/fDimensionlessL);
								Node[nodeId-1].SetCoord(2,nodecoordZ/fDimensionlessL);
								gmesh->NodeVec()[nodeid] = Node[nodeId-1];
								ContNode++;
							}
							if(ContNode == numnodes)
							{
								strtemp = "";	
							}					
						}
							break;
						case 3:
						{
							//"LINEAR"
							if (GeneralData[cont] != 0)
							{
								read >> elementId;
								read >> matElId;  // Material ID 
								read >> Layerid;						
								read >> TopolLine[0]; //node 2
								read >> TopolLine[1]; //node 3
								elementId--;
								TopolLine[0]--;
								TopolLine[1]--;
								new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementId,TopolLine,matElId,*gmesh);
								ContLine++;					
							}
							if(ContLine == elements1D)
							{
								strtemp = "";	
							}						
						}
							break;
						case 4:
						{
							if (GeneralData[cont] != 0)
							{
								//"TRIANGLE"
								read >> elementId;
								read >> matElId;  // Material ID
								read >> Layerid;						
								read >> TopolTriangle[0]; //node 1
								read >> TopolTriangle[1]; //node 2
								read >> TopolTriangle[2]; //node 3
								elementId--;						
								TopolTriangle[0]--;
								TopolTriangle[1]--;
								TopolTriangle[2]--;					
								new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementId, TopolTriangle, matElId, *gmesh);
								ContTrian++;							
							}
							if(ContTrian == elements2DT)
							{
								strtemp = "";	
							}					
						}
							break;
						case 5:
						{
							//"QUADRILATERAL"
							if (GeneralData[cont] != 0)
							{
								read >> elementId;
								read >> matElId;  // Material ID
								read >> Layerid;						
								read >> TopolQuad[0]; //node 1
								read >> TopolQuad[1]; //node 2
								read >> TopolQuad[2]; //node 3
								read >> TopolQuad[3]; //node 4						
								elementId--;						
								TopolQuad[0]--;
								TopolQuad[1]--;
								TopolQuad[2]--;	
								TopolQuad[3]--;							
								new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementId, TopolQuad, matElId, *gmesh);
								ContQuad++;						
							}
							if(ContQuad == elements2DQ)
							{
								strtemp = "";	
							}						
						}
							break;
						case 6:
						{
							//"TETRAHEDRA"
							if (GeneralData[cont] != 0)
							{
								read >> elementId;
								read >> matElId;  // Material ID
								read >> Layerid;						
								read >> TopolTet[0]; //node 1
								read >> TopolTet[1]; //node 2
								read >> TopolTet[2]; //node 3
								read >> TopolTet[3]; //node 4						
								elementId--;						
								TopolTet[0]--;
								TopolTet[1]--;
								TopolTet[2]--;	
								TopolTet[3]--;	
								new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (elementId, TopolTet, matElId, *gmesh);
								ContTet++;
							}
							if(ContTet == elements3DT)
							{
								strtemp = "";	
							}
						}
							break;
						case 7:
						{
							//"HEXAHEDRA"
							if (GeneralData[cont] != 0)
							{
								read >> elementId;
								read >> matElId;  // Material ID
								read >> Layerid;						
								read >> TopolHex[0]; //node 1
								read >> TopolHex[1]; //node 2
								read >> TopolHex[2]; //node 3
								read >> TopolHex[3]; //node 4
								read >> TopolHex[4]; //node 5
								read >> TopolHex[5]; //node 6
								read >> TopolHex[6]; //node 7
								read >> TopolHex[7]; //node 8						
								elementId--;						
								TopolHex[0]--;
								TopolHex[1]--;
								TopolHex[2]--;	
								TopolHex[3]--;
								TopolHex[4]--;
								TopolHex[5]--;
								TopolHex[6]--;	
								TopolHex[7]--;						
								new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elementId, TopolHex, matElId, *gmesh);					
								ContHex++;
							}
							if(ContHex == elements3DH)
							{
								strtemp = "";	
							}
						}
							break;					
						default:
						{
							strtemp = "";
						}
							break;
					}		
					
					if(strtemp == "" || strtemp == "\r") 
					{
						FlagString = "";
						cont++;
						break;
					}
				}
				
			}	
			
			
		}
	}	
	
	std::cout << "Read General Mesh Data -> done!" << std::endl;	
	gmesh->BuildConnectivity();
	std::cout << "Geometric Mesh Connectivity -> done!" << std::endl;
	return gmesh;
	
}// End Method

void TPZReadGIDGrid::SetfDimensionlessL(REAL fDimensionlessLValue)
{
	fDimensionlessL = fDimensionlessLValue;
}
