#ifndef LAPLACESOLVER_H
#define LAPLACESOLVER_H

#include <iostream>
using namespace std;
#include "fvmesh.h"
#include "matrix.h"
#include "math.h"
#include "iostream"
#include "exportensight.h"
#include "GlobalHeader/glbheader.h"



class LaplaceSolver
{
public:
    Fvmesh *fvmesh;
    matrix *Matrix;
    ExportEnSight *exportensight;

    REAL Gamma;
    REAL ro;
    REAL alphaU;
    REAL alphaP;

    LaplaceSolver();
    LaplaceSolver(Fvmesh *_fvmesh , INT _celltype);

    void Run(INT _celltype ,INT _equationtype , INT _gridtype);  // _celltype 0 : triangle , _celltype 1: quad
                                                                // equationtype 0 : triangle diffusion , equationtype 1 : quad diffusion , equtiontyeq 2: quad diffusion - convection (laplace transport equaiton)
                                                               // gridtype 0 : triangle , gridtype 1: quad
    INT FaceN();
    INT CellN();
    void TriCalCorrectionValue (REAL *_X, vector<REAL> &_corr, vector<REAL> &_facevalue , vector<REAL> &_flux , vector<INT> &_bctype);
    void QuadCalCorrectionValue(REAL *_X, vector<REAL> &_corr, vector<REAL> &_facevalue , vector<REAL> &_flux , vector<INT> &_bctype);
    void VectorAdd(vector<REAL> &_corr,REAL *_Bvector,REAL *_totalB);


    void TriMatrixA(INT _CellN ,vector<INT> &_bctype);
    void TrivectorB(INT _CellN ,vector<REAL> &_facevalue ,vector<REAL> &_flux,vector<INT> &_bctype,REAL *__Bvector);

    void QuadMatrixA(INT _CellN ,vector<INT> &_bctype);
    void QuadMatrixA(INT _CellN ,vector<INT> &_bctype,vector<vector<REAL>> &_velocity);
    void CQuadMatrixA(INT _CellN , vector<INT> &_bctype ,vector<vector<REAL>> &_velocity);
    void QuadvectorB(INT _CellN ,vector<REAL> &_facevalue ,vector<REAL> &_flux,vector<INT> &_bctype,REAL *__Bvector);
    void CQuadvectorB(INT _CellN ,vector<REAL> &_facevalue, vector<REAL> &_flux ,vector<INT> &_bctype,vector<vector<REAL>> &_velocity,REAL *__Bvector);

    void QuadvectorB(INT _CellN ,vector<REAL> &_facevalue ,vector<REAL> &_flux,vector<INT> &_bctype,vector<vector<REAL>> &_velocity,REAL *__Bvector);

    void FacePI(REAL *_X , vector<REAL> &__facePI);
    void TriCellgradient(INT _CellN, vector<REAL> &_facevalue,vector<REAL> &_flux,REAL *_X, vector<INT> &_bctype, vector<vector<REAL>> &__CellGrad);
    void QuadCellgradient(INT _CellN, vector<REAL> &_facevalue,vector<REAL> &_flux,REAL *_X,vector<INT> &_bctype, vector<vector<REAL>> &__CellGrad);
    void Facegradient(INT _FaceN, vector<vector<REAL>> &_cellGrad0, vector<vector<REAL>> &_cellGrad1, vector<vector<REAL>> &__FaceGrad,REAL *_X);
    void NormalVectorandZita0(INT _neicell0, INT _neicell1,vector<REAL> &_NormalVector,INT _faceID, vector<vector<REAL>> &__Value );
    void NormalVectorandZita1(INT _neicell0, INT _neicell1,vector<REAL> &_NormalVector,INT _faceID, vector<vector<REAL>> &__Value );
    void NormalVectorandZita2(INT _faceID, INT _CellID,vector<REAL> &_NormalVector, vector<vector<REAL>> &__Value );// unit(Normalvector - ZitaVector)
    void cellvelocity();
    void write0(REAL *_X);
    void write1(REAL *_X);
    void write2(REAL (*_V)[3]);


    void PQuadMatrixA(INT _CellN ,vector<INT> &_bctype,vector<vector<REAL>> &_velocity);
};

#endif // LAPLACESOLVER_H
