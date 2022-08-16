#include "laplacesolver.h"

LaplaceSolver::LaplaceSolver()
{
    Gamma = 0.0;
    ro = 0.0;

}
LaplaceSolver::LaplaceSolver(Fvmesh *_fvmesh ,INT _celltype)
{
    fvmesh = _fvmesh;
    Gamma = 0.01;
    ro  = 1;
    alphaU = 0.5;
    alphaP = 0.8;

    switch(_celltype)
    {
    case 0:
    {
        INT nodeXn0 = fvmesh->getXN();
        INT nodeYn0 = fvmesh->getYN();

        INT totalN0 = ( (nodeXn0-1) * ((nodeYn0 -1) * 2) );

        Matrix = new matrix(totalN0,totalN0);
        break;
    }
    case 1:
        INT nodeXn1 = fvmesh->getXN();
        INT nodeYn1 = fvmesh->getYN();

        INT totalN1 = ( (nodeXn1-1) * (nodeYn1-1) );

        Matrix = new matrix(totalN1 , totalN1);

        break;
    }

}

INT LaplaceSolver::FaceN()
{
    return fvmesh->getfaecN();
}

INT LaplaceSolver::CellN()
{
    return fvmesh->getcellN();
}

void LaplaceSolver::Run(INT _celltype ,INT _equationtype , INT _gridtype)
{
    INT cn , fn;
    INT xn , yn;
    cn = this->CellN();
    fn = this->FaceN();
    xn = fvmesh->getXN();
    yn = fvmesh->getYN();

    INT PatchN = 4;
    vector<vector<REAL>> bcvalue(PatchN,vector<REAL>(2,0.0));
    vector<REAL> facepi(yn,0);
    REAL Bvector[cn];
    REAL X[cn];
//    REAL totalB[cn];
    vector<INT> pif(fn,0);
    vector<REAL> corr(cn,0);
    vector<vector<REAL>> V(fn,vector<REAL>(2,0));

    REAL bcvalueX[PatchN][2];
    REAL bcvalueY[PatchN][2];


    for(INT i=0; i<cn ; i++)
    {
        Bvector[i] = X[i] = 0.0;
    }

    fvmesh->boundry(facepi);


    switch(_celltype)
    {
    case 0:
        fvmesh->geomesh->TriPatchidofFace(pif,xn,yn);
        break;
    case 1:
        fvmesh->geomesh->QuadPatchidofFace(pif,xn,yn);
        break;
    }

    bcvalue[0][0] = 0;  // bcvalue[i][0]  = Dirichlet , bcvalue[i][1] = Neumann
    bcvalue[0][1] = 0;

    bcvalue[1][0] = 0;
    bcvalue[1][1] = 0;

    bcvalue[2][0] = 0;
    bcvalue[2][1] = 0;

    bcvalue[3][0] = 0;
    bcvalue[3][1] = 0;
//  ㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡ//
    bcvalueX[0][0] = 0;  // bcvalue[i][0]  = Dirichlet , bcvalue[i][1] = Neumann
    bcvalueX[0][1] = 0;

    bcvalueX[1][0] = 0;
    bcvalueX[1][1] = 0;

    bcvalueX[2][0] = 1.e-4;
    bcvalueX[2][1] = 0;

    bcvalueX[3][0] = 0;
    bcvalueX[3][1] = 0;
// ㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡ//
    bcvalueY[0][0] = 0;  // bcvalue[i][0]  = Dirichlet , bcvalue[i][1] = Neumann
    bcvalueY[0][1] = 0;

    bcvalueY[1][0] = 0;
    bcvalueY[1][1] = 0;

    bcvalueY[2][0] = 0;
    bcvalueY[2][1] = 0;

    bcvalueY[3][0] = 0;
    bcvalueY[3][1] = 0;


//    double facevalueX[fn];
//    double facevalueY[fn];

    INT Neicell;
    INT Pid;
    vector<INT> bctype(fn,-1);
    vector<REAL> facevalue(fn,0);
    vector<REAL> flux(fn,0);
    INT A;
    A=0;


    for(INT i=0 ; i<fn; i++)
    {
        facevalue[i] = 0.0; // Dirichlet
        flux[i] = 0.0; // Neumann
        bctype[i] = -1;
        Neicell = fvmesh->getNeicell(i);
        if(Neicell == -1)
        {
            Pid = pif[i];

            if(Pid == 0)
            {
                bctype[i] = 1 ;
                facevalue[i] = bcvalue[Pid][1];
//                facevalueX[i] = bcvalueX[Pid][0];
//                facevalueY[i] = bcvalueY[Pid][0];
//                flux[i] = bcvalue[Pid][1];
//                fluxX[i] = bcvalueX[Pid][1];
//                fluxY[i] = bcvalueY[Pid][1];
            }
            else if(Pid == 1)
            {
                bctype[i] = 1;
                facevalue[i] = bcvalue[Pid][1];
//                facevalueX[i] = bcvalueX[Pid][0];
//                facevalueY[i] = bcvalueY[Pid][0];
//                flux[i] = bcvalue[Pid][1];
//                fluxX[i] = bcvalueX[Pid][1];
//                fluxY[i] = bcvalueY[Pid][1];
            }
            else if(Pid == 2)
            {
                bctype[i] = 0;
                facevalue[i] = bcvalue[Pid][0];
//                facevalueX[i] = bcvalueX[Pid][0];
//                facevalueY[i] = bcvalueY[Pid][0];
//                flux[i] = bcvalue[Pid][1];
//                fluxX[i] = bcvalueX[Pid][1];
//                fluxY[i] = bcvalueY[Pid][1];
            }
            else if(Pid == 3)
            {
                bctype[i] = 0;
//                facevalue[i] = bcvalue[Pid][1];
//                facevalueX[i] = bcvalueX[Pid][0];
//                facevalueY[i] = bcvalueY[Pid][0];
//                flux[i] = bcvalue[Pid][1];
//                fluxX[i] = bcvalueX[Pid][1];
//                fluxY[i] = bcvalueY[Pid][1];
                facevalue[i] = facepi[A];

                A++;
            }
            else{};

        }

//        cout << "id: "<<i << endl;
//        cout << facevalue[i] << endl;
    }

fvmesh->facevelocity0(V);

    switch(_equationtype)
    {
    case 0:
//        this->TriMatrixA(cn,bctype);
//        this->TrivectorB(cn,facevalue,flux,bctype,Bvector);
        break;
    case 1:
//        this->QuadMatrixA(cn,bctype);
//        this->QuadvectorB(cn,facevalue,flux,bctype,Bvector);
        break;
    case 2:
//        this->QuadMatrixA(cn,bctype,V);
        this->PQuadMatrixA(cn,bctype,V);
        this->QuadvectorB(cn,facevalue,flux,bctype,V,Bvector);
//        this->CQuadMatrixA(cn,bctype,V);
//        this->CQuadvectorB(cn,facevalue,flux,bctype,V,Bvector);
//         Matrix->show_element();
        break;
    }



Matrix->GaussSidel(10,1.e-5,Bvector,X);


    this->cellvelocity();

    switch (_gridtype)
    {
    case 0:
        this->write0(X);
        break;
    case 1:
//        this->write1(P);
        this->write1(X);
//        this->write2(NewUV);
        break;
    }

}

void LaplaceSolver::TriMatrixA(INT _CellN ,vector<INT> &_bctype)
{
    REAL D, S;
    INT Pcid,Neicid;
    INT Fid;
    INT temp;

    REAL value;
    for(int i=0 ; i<_CellN; i++)
    {
        Pcid = i;
        for(int j=0; j<3; j++)
        {
            value =0.0;
            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);

            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->area(Fid);

                    value = (Gamma * S ) /D ; //change
                    Matrix->add_each_element(Pcid,Pcid,value);

                }

                else if(_bctype[Fid] == 1)
                {

                }
            }

            else
            {
                D = fvmesh->CellCenterDistance(fvmesh->geomesh->facelist.getneicells0(fvmesh->geomesh->cellist.getcellfid(i,j)),fvmesh->geomesh->facelist.getneicells1(fvmesh->geomesh->cellist.getcellfid(i,j)));
                S = fvmesh->area(Fid);
                if(Neicid == i)
                {
                    temp = fvmesh->geomesh->facelist.getneicells0(Fid);
                    Neicid = temp;

                    value = (-1.0) * (Gamma * S) /D  ;  //change
                    Matrix->set_each_element(Pcid,Neicid,value);
                }
                else
                {
                    value = (-1.0) * (Gamma * S) /D ; //change
                    Matrix->set_each_element(Pcid,Neicid,value);
                }

                value= (Gamma * S) /D ; //change
                Matrix->add_each_element(Pcid,Pcid,value);

            }
        }
     }
}

void LaplaceSolver::QuadMatrixA(INT _CellN ,vector<INT> &_bctype,vector<vector<REAL>> &_velocity)
{
    INT fn;
    fn = this->FaceN();
    REAL D, S ,d;
    INT Pcid,Neicid;
    INT Fid;
    INT temp;
    vector<REAL> Vector(3,0.0);
    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0.0));
    vector<REAL> dotvalue(fn,0.0);
    REAL C ;
    REAL dvalue,cvalue , ccvalue;
    vector<REAL> Peclet(fn,0.0);
    REAL m;

    for(int i=0 ; i<_CellN; i++)
    {
        Pcid = i;
        for(int j=0; j<4; j++)
        {

            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];


            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->geomesh->facelist.getarea(Fid);

//                    cout << D << "   " << S << "   " << Gamma << endl;
                    dvalue = (Gamma * S ) /D ; //change
                    Matrix->add_each_element(Pcid,Pcid,dvalue);
                }

                else if(_bctype[Fid] == 1)
                {
                    S = fvmesh->area(Fid);
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);

                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    cvalue = ( ro * S * dotvalue[Fid] );
                    Matrix->add_each_element(Pcid,Pcid,cvalue);

                }
            }

            else
            {
                D = fvmesh->CellCenterDistance(fvmesh->getcellID(Fid),fvmesh->getNeicell(Fid));
                S = fvmesh->area(Fid);
                d = fvmesh->CellandFaceCenterDistance(i,Fid);

//                cout << D << "  " << S << endl;

                if(Neicid == i)
                {
                    fvmesh->peclet1(ro,Gamma,d,Fid,_velocity,Peclet);

                    temp = fvmesh->getcellID(Fid);
                    Neicid = temp;
                    dvalue = (-1.0) * (Gamma * S) /D  ;  //change
                    Matrix->set_each_element(Pcid,Neicid,dvalue);

                    C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                    m = ro * S * dotvalue[Fid];

                    if(Peclet[Fid] >2 )
                    {
                        if(m > 0)
                        {
                            cvalue = (ro * S * dotvalue[Fid] );

                            Matrix->add_each_element(Pcid,Pcid,cvalue);
                        }

                        else
                        {
                            cvalue = (ro * S * dotvalue[Fid] );
                            Matrix->add_each_element(Pcid,Neicid,cvalue);
                        }
                    }

                    else
                    {
                        cvalue = (( ro * S * C * dotvalue[Fid]));
                        Matrix->add_each_element(Pcid,Neicid,cvalue);


                        ccvalue = ( ro * S * (1-C) * dotvalue[Fid]) ; //change
                        Matrix->add_each_element(Pcid,Pcid,ccvalue);
                    }
                }
                else
                {
                    fvmesh->peclet1(ro,Gamma,d,Fid,_velocity,Peclet);

                    dvalue = (-1.0) * (Gamma * S) /D ; //change
                    Matrix->set_each_element(Pcid,Neicid,dvalue);

                    C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                    m = ro * S * dotvalue[Fid];

                    if(Peclet[Fid] >2 )
                    {
                        if(m > 0)
                        {
                            cvalue = (ro * S * dotvalue[Fid] );

                            Matrix->add_each_element(Pcid,Pcid,cvalue);
                        }

                        else
                        {
                            cvalue = (ro * S * dotvalue[Fid] );
                            Matrix->add_each_element(Pcid,Neicid,cvalue);
                        }
                    }

                    else
                    {
                        cvalue = (( ro * S * C * dotvalue[Fid]));
                        Matrix->add_each_element(Pcid,Neicid,cvalue);


                        ccvalue = ( ro * S * (1-C) * dotvalue[Fid]) ; //change
                        Matrix->add_each_element(Pcid,Pcid,ccvalue);
                    }

                }

                C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                dvalue = (Gamma * S) /D ; //change
                Matrix->add_each_element(Pcid,Pcid,dvalue);

//                cvalue =(-1.0) * ( ro * S * (1-C) * dotvalue[Fid]) ; //change
//                Matrix->add_each_element(Pcid,Pcid,cvalue);
            }
        }

     }
}


void LaplaceSolver::PQuadMatrixA(INT _CellN ,vector<INT> &_bctype,vector<vector<REAL>> &_velocity)
{
    INT fn;
    fn = this->FaceN();
    REAL D, S,d;
    INT Pcid,Neicid;
    INT Fid;
    INT temp;
    vector<REAL> Vector(3,0.0);
    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0.0));
    vector<REAL> dotvalue(fn,0.0);
    REAL C ;
    REAL dvalue,cvalue , ccvalue;
    vector<REAL> Peclet(fn,0.0);
    REAL m;

    for(int i=0 ; i<_CellN; i++)
    {
        Pcid = i;
        for(int j=0; j<4; j++)
        {

            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];


            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->geomesh->facelist.getarea(Fid);

//                    cout << D << "   " << S << "   " << Gamma << endl;
                    dvalue = (Gamma * S ) /D ; //change
                    Matrix->add_each_element(Pcid,Pcid,dvalue);
                }

                else if(_bctype[Fid] == 1)
                {
                    S = fvmesh->area(Fid);
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);

                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    cvalue = ( ro * S * dotvalue[Fid] );
                    Matrix->add_each_element(Pcid,Pcid,cvalue);

                }
            }

            else
            {
                D = fvmesh->CellCenterDistance(fvmesh->getcellID(Fid),fvmesh->getNeicell(Fid));
                d = fvmesh->CellandFaceCenterDistance(i,Fid);
                S = fvmesh->area(Fid);

//                cout << d << "  " << S << endl;

                if(Neicid == i)
                {
                    temp = fvmesh->getcellID(Fid);
                    Neicid = temp;
                    dvalue = (-1.0) * (Gamma * S) /D  ;  //change
                    Matrix->set_each_element(Pcid,Neicid,dvalue);

                    C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                    cvalue = (( ro * S * C * dotvalue[Fid]));
                    Matrix->add_each_element(Pcid,Neicid,cvalue);
                }
                else
                {
                    dvalue = (-1.0) * (Gamma * S) /D ; //change
                    Matrix->set_each_element(Pcid,Neicid,dvalue);

                    C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                    cvalue = (( ro * S * C * dotvalue[Fid]));
                    Matrix->add_each_element(Pcid,Neicid,cvalue);
                }

                C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                dvalue = (Gamma * S) /D ; //change
                Matrix->add_each_element(Pcid,Pcid,dvalue);

                cvalue = ( ro * S * (1-C) * dotvalue[Fid]) ; //change
                Matrix->add_each_element(Pcid,Pcid,cvalue);
            }
        }

     }
}
void LaplaceSolver::CQuadMatrixA(INT _CellN , vector<INT> &_bctype ,vector<vector<REAL>> &_velocity)
{
    INT fn;
    fn = this->FaceN();
    REAL S;
    INT Pcid,Neicid;
    INT Fid;
    INT temp;
    vector<REAL> Vector(3,0);
    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0));
    vector<REAL> dotvalue(fn,0);
    REAL C ;
    REAL value;

    for(int i=0 ; i<_CellN; i++)
    {
        Pcid = i;
        for(int j=0; j<4; j++)
        {
            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];

            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {

                }

                else if(_bctype[Fid] == 1)
                {
                    S = fvmesh->area(Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    value = ( ro * S * dotvalue[Fid] );
                    Matrix->add_each_element(Pcid,Pcid,value);
                }
            }

            else
            {
                S = fvmesh->geomesh->facelist.getarea(Fid);
                if(Neicid == i)
                {
                    C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    temp = fvmesh->getcellID(Fid);
                    Neicid = temp;
                    value = ( ro * S * C * dotvalue[Fid]);
                    Matrix->add_each_element(Pcid,Neicid,value);

                }
                else
                {
                    C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    value = ( ro * S * C * dotvalue[Fid]);
                    Matrix->add_each_element(Pcid,Neicid,value);
                }

                C = fvmesh->DistanceRatioCoefficient(i,Neicid,Fid);
                fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

//                cout << C << endl;

                value = ( ro * S * (1-C) * dotvalue[Fid]) ; //change
                Matrix->add_each_element(Pcid,Pcid,value);
            }
        }
     }
}


void LaplaceSolver::TrivectorB(INT _CellN ,vector<REAL> &_facevalue ,vector<REAL> &_flux,vector<INT> &_bctype,REAL *__Bvector)
{

    REAL D, S;
    INT fn =this->FaceN();
    INT Pcid,Neicid;
    INT Fid;
    vector<REAL> Vector(3,0);
    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0));
    REAL flux;

    for(int i=0 ; i<_CellN; i++)
    {
        __Bvector[i] = 0.0;
        Pcid = i;
        for(int j=0; j<3; j++)
        {
            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];

            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->area(Fid);

                    __Bvector[Pcid] += ((Gamma * S *_facevalue[Fid]) /D) ; //change
//                    cout << __Bvector[Pcid] <<endl;
                }

                else if(_bctype[Fid] == 1)
                {
                    flux =( (-1) * ((norvector[Fid][0] * _flux[Fid]) + (norvector[Fid][1] * _flux[Fid])) );
                    S = fvmesh->area(Fid);
                    __Bvector[Pcid] += (flux) * S ;


                }

                else{};
            }

            else
            {

            }
        }

//        cout << __Bvector[Pcid] << endl;
     }
}

void LaplaceSolver::QuadvectorB(INT _CellN ,vector<REAL> &_facevalue ,vector<REAL> &_flux,vector<INT> &_bctype,REAL *__Bvector)
{
    REAL D, S;
    INT fn =this->FaceN();
    INT Pcid,Neicid;
    INT Fid;
    vector<REAL> Vector(3,0);
//    double norvector[fn][2];
    vector<REAL> dotvalue(fn,0);
    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0));
    REAL flux;
    vector<REAL> facePI(fn,0);

    for(int i=0 ; i<_CellN; i++)
    {
        __Bvector[i] = 0.0;
        Pcid = i;
        for(int j=0; j<4; j++)
        {
            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];

            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->area(Fid);

                    __Bvector[Pcid] += ((Gamma * S *_facevalue[Fid]) /D) ; //change
//                    cout << __Bvector[Pcid] <<endl;
                }

                else if(_bctype[Fid] == 1)
                {
                    flux =( (-1) * ((norvector[Fid][0] * _flux[Fid]) + (norvector[Fid][1] * _flux[Fid])) );
                    S = fvmesh->area(Fid);
                    __Bvector[Pcid] += (flux) * S ;

//                    cout << __Bvector[Pcid] << endl;
                }

                else{};
            }

            else
            {

            }
        }

//        cout << __Bvector[Pcid] << endl;
     }
}

void LaplaceSolver::QuadvectorB(INT _CellN ,vector<REAL> &_facevalue ,vector<REAL> &_flux,vector<INT> &_bctype,vector<vector<REAL>> &_velocity,REAL *__Bvector)
{
    REAL D, S;
    INT fn =this->FaceN();
    INT Pcid,Neicid;
    INT Fid;
    vector<REAL> Vector(3,0);
//    double norvector[fn][2];
    vector<REAL> dotvalue(fn,0);
    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0));
    REAL flux;
    vector<REAL> facePI(fn,0);

    for(int i=0 ; i<_CellN; i++)
    {
        __Bvector[i] = 0.0;
        Pcid = i;
        for(int j=0; j<4; j++)
        {
            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];

            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->area(Fid);
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);

                    __Bvector[Pcid] += ((Gamma * S *_facevalue[Fid]) /D) ; //change
                    __Bvector[Pcid] += ( (-1) * (ro * S * dotvalue[Fid] * _facevalue[Fid]) ) ;
//                    cout << __Bvector[Pcid] <<endl;
                }

                else if(_bctype[Fid] == 1)
                {
                    flux =( (-1) * ((norvector[Fid][0] * _flux[Fid]) + (norvector[Fid][1] * _flux[Fid])) );
                    S = fvmesh->area(Fid);
                    __Bvector[Pcid] += ((flux) * S)  ;

                }

                else{};
            }

            else
            {

            }
        }

//        cout << __Bvector[Pcid] << endl;
     }
}


void LaplaceSolver::CQuadvectorB(INT _CellN ,vector<REAL> &_facevalue, vector<REAL> &_flux ,vector<INT> &_bctype,vector<vector<REAL>> &_velocity,REAL *__Bvector)
{
    INT fn = this->FaceN();
    REAL D, S;
    INT Pcid,Neicid;
    INT Fid;
    vector<REAL> Vector(3,0);

    vector<vector<REAL>> norvector(fn,vector<REAL>(3,0));
    vector<REAL> dotvalue(fn,0);
    vector<REAL> facePI(fn,0);

    for(int i=0 ; i<_CellN; i++)
    {
        Pcid = i;
        for(int j=0; j<4; j++)
        {
            Fid = fvmesh->getFid(Pcid,j);
            Neicid = fvmesh->getNeicell(Fid);
            fvmesh->getNormalVector(i,Fid,Vector);

            norvector[Fid][0] = Vector[0];
            norvector[Fid][1] = Vector[1];
            norvector[Fid][2] = Vector[2];

            if(Neicid == -1)
            {
                if(_bctype[Fid] == 0)
                {
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    S = fvmesh->area(Fid);

                    __Bvector[Pcid] += ( (-1) * (ro * S * dotvalue[Fid] * _facevalue[Fid]) ) ; //change
//                    cout << __Bvector[Pcid] <<endl;
                }

                else if(_bctype[Fid] == 1)
                {
                    fvmesh->dotproduct(norvector,_velocity,Fid,dotvalue);
                    D = fvmesh->CellandFaceCenterDistance(i,Fid);
                    S = fvmesh->area(Fid);

                    facePI[Fid] = ( _flux[Fid] * D );

                    __Bvector[Pcid] += ( (-1) * (ro * S * dotvalue[Fid] * facePI[Fid]) ) ; // cg

//                    cout << __Bvector[Pcid] << endl;
                }

                else{};
            }

            else
            {

            }
        }

//        cout << __Bvector[Pcid] << endl;
     }
}


void LaplaceSolver::FacePI(REAL *_X , vector<REAL> &__facePI)
{
    INT fn;
    fn = this->FaceN();
    INT neicell0;
    INT neicell1;

    for(int i=0; i<fn; i++)
    {
        neicell0=fvmesh->getcellID(i);
        neicell1=fvmesh->getNeicell(i);

        if(neicell1 != -1 )
        {
            __facePI[i] =fvmesh->CellcenterScalarDistanceCoe(_X[neicell0],_X[neicell1],i,neicell0,neicell1);
        }

        else{};
//        cout << __facePI[i] << endl;
    }
}

void LaplaceSolver::TriCellgradient(INT _CellN, vector<REAL> &_facevalue,vector<REAL> &_flux,REAL *_X, vector<INT> &_bctype, vector<vector<REAL>> &__CellGrad)
{
    for(int i=0 ; i<_CellN; i++)
    {
//        cout << i << endl;
        fvmesh->TriGradientPi(i,_facevalue,_flux,_X,_bctype,__CellGrad);
    }

//    for(int i = 0; i <_CellN; i++)
//    {
//       cout << "id:" << i << endl;
//        cout << __CellGrad[i][0] << "  " << __CellGrad[i][1] << endl;
//    }
}


void LaplaceSolver::QuadCellgradient(INT _CellN, vector<REAL> &_facevalue,vector<REAL> &_flux,REAL *_X,vector<INT> &_bctype, vector<vector<REAL>> &__CellGrad)
{
    for(int i=0 ; i<_CellN; i++)
    {
//        cout << i << endl;
        fvmesh->QuadGradientPi(i,_facevalue,_flux,_X,_bctype,__CellGrad);
    }
}
void LaplaceSolver::Facegradient(INT _FaceN,vector<vector<REAL>> &_cellGrad0, vector<vector<REAL>> &_cellGrad1, vector<vector<REAL>> &__FaceGrad,REAL *_X)
{
    INT neicell0;
    INT neicell1;
    INT cellID;
    INT faceID0;
    INT faceID1;
//    double D0[2];
//    double D1[2];
    REAL cell0,cell1;

    vector<REAL> D0(2,0);
    vector<REAL> D1(2,0);


    for(int i=0; i<_FaceN; i++)
    {
//        cout << i<< endl;

        neicell0= fvmesh->geomesh->facelist.getneicells0(i);
        neicell1= fvmesh->geomesh->facelist.getneicells1(i);

        cellID =  fvmesh->geomesh->facelist.getneicells0(i);
        faceID0 =  fvmesh->geomesh->facelist.getfacenids0(i);
        faceID1 = fvmesh->geomesh->facelist.getfacenids1(i);

        if(neicell1 == -1)
        {
            D0[0] = fvmesh->geomesh->nodelist.getcoordinateX(faceID0) - fvmesh->geomesh->cellist.GetCenterX(cellID);
            D0[1] = fvmesh->geomesh->nodelist.getcoordinateY(faceID0) - fvmesh->geomesh->cellist.GetCenterY(cellID);

            D1[0] = fvmesh->geomesh->nodelist.getcoordinateX(faceID1) - fvmesh->geomesh->cellist.GetCenterX(cellID);
            D1[1] = fvmesh->geomesh->nodelist.getcoordinateY(faceID1) - fvmesh->geomesh->cellist.GetCenterY(cellID);

            cell0 = _X[cellID] + ( (_cellGrad0[cellID][0] * D0[0]) + (_cellGrad0[cellID][1] * D0[1]) );
            cell1 = _X[cellID] + ( (_cellGrad0[cellID][0] * D1[0]) + (_cellGrad0[cellID][1] * D1[1]) );

            __FaceGrad[i][0] = ((cell0 - cell1) * (D0[0] - D1[0])) / ( (D0[0]-D1[0]) * (D0[0]-D1[0]) + ((D0[1]-D1[1]) * (D0[1]-D1[1])) );
            __FaceGrad[i][1] = ((cell0 - cell1) * (D0[1] - D1[1])) / ( (D0[0]-D1[0]) * (D0[0]-D1[0]) + ((D0[1]-D1[1]) * (D0[1]-D1[1])) );
            __FaceGrad[i][2] = 0.0;

//                    cout <<"id:" << i << endl;
//                    cout << __FaceGrad[i][0] << "  "  << __FaceGrad[i][1] << endl;
        }

        else
        {
            fvmesh->CellcenterVectorDistanceCoe(_cellGrad0,_cellGrad1,i,neicell0,neicell1,__FaceGrad);

//            cout << "id:" << i << endl;
        }
    }
//    for(int i=0; i<_FaceN; i++)
//    {
//        cout <<"id:" << i << endl;
//        cout << __FaceGrad[i][0] << "  "  << __FaceGrad[i][1] << endl;
//    }
}


void LaplaceSolver::NormalVectorandZita0(INT _neicell0, INT _neicell1,vector<REAL> &_NormalVector,INT _faceID, vector<vector<REAL>> &__Value )
{
    INT fn;
    fn = this->FaceN();
    REAL VectorX,VectorY,VectorZ;
    V3 Zita;

    for(int i=0; i<fn; i++)
    {
        __Value[i][0] = __Value[i][1] = __Value[i][2] = 0.0;
    }

    VectorX = fvmesh->geomesh->cellist.GetCenterX(_neicell1) - fvmesh->geomesh->cellist.GetCenterX(_neicell0);
    VectorY = fvmesh->geomesh->cellist.GetCenterY(_neicell1) - fvmesh->geomesh->cellist.GetCenterY(_neicell0);
    VectorZ = 0.0;

    V3 Vector(VectorX,VectorY,VectorZ);

    Zita = UNIT(Vector);

    __Value[_faceID][0] = ( _NormalVector[0] ) - Zita[0];
    __Value[_faceID][1] = ( _NormalVector[1] ) - Zita[1];
    __Value[_faceID][2] = ( _NormalVector[2] ) - Zita[2];

}



void LaplaceSolver::NormalVectorandZita1(INT _neicell0, INT _neicell1,vector<REAL> &_NormalVector,INT _faceID,vector<vector<REAL>> &__Value )
{
    INT fn;
    fn = this->FaceN();

    V3 Zita;
    REAL VectorX,VectorY,VectorZ;

    for(int i=0; i<fn; i++)
    {
        __Value[i][0] = __Value[i][1] = __Value[i][2] = 0.0;
    }

    VectorX = fvmesh->geomesh->cellist.GetCenterX(_neicell0) - fvmesh->geomesh->cellist.GetCenterX(_neicell1);
    VectorY = fvmesh->geomesh->cellist.GetCenterY(_neicell0) - fvmesh->geomesh->cellist.GetCenterY(_neicell1);
    VectorZ = 0.0;

    V3 Vector(VectorX,VectorY,VectorZ);

    __Value[_faceID][0] = ( _NormalVector[0] ) - Zita[0];
    __Value[_faceID][1] = ( _NormalVector[1] ) - Zita[1];
    __Value[_faceID][2] = ( _NormalVector[2] ) - Zita[2];

}

void LaplaceSolver::NormalVectorandZita2(INT _faceID, INT _CellID,vector<REAL> &_NormalVector,vector<vector<REAL>> &__Value )
{
    INT fn;
    fn = this->FaceN();
    V3 Zita;
    REAL VectorX,VectorY,VectorZ;

    for(int i=0; i<fn; i++)
    {
        __Value[i][0] = __Value[i][1] = __Value[i][2];
    }

    VectorX = fvmesh->geomesh->facelist.getCenterX(_faceID) - fvmesh->geomesh->cellist.GetCenterX(_CellID);
    VectorY = fvmesh->geomesh->facelist.getCenterY(_faceID) - fvmesh->geomesh->cellist.GetCenterY(_CellID);
    VectorZ = 0.0;

    V3 Vector(VectorX,VectorY,VectorZ);
    Zita = UNIT(Vector);
//    cout << Zita[_faceID][0] << "  " << Zita[_faceID][1] << endl;
//    cout << _NormalVector[0] << "  " << _NormalVector[1] << endl;

//    cout << endl;

    __Value[_faceID][0] = ( _NormalVector[0] ) - Zita[0];
    __Value[_faceID][1] = ( _NormalVector[1] ) - Zita[1];
    __Value[_faceID][2] = ( _NormalVector[2] ) - Zita[2];


//    cout << __Value[_faceID][0] << "  " << __Value[_faceID][1] <<endl;
}


void LaplaceSolver::TriCalCorrectionValue(REAL *_X, vector<REAL> &_corr, vector<REAL> &_facevalue , vector<REAL> &_flux , vector<INT> &_bctype)
{
    INT fn,cn;
    fn = this->FaceN();
    cn = this->CellN();

    INT neicell1;
    INT cellID;

// facePI

    this->FacePI(_X,_facevalue);

// cellgradient
//    double cellGrad0[cn][2];
//    double cellGrad1[cn][2];

    vector<vector<REAL>> cellGrad0(cn,vector<REAL>(3,0));
    vector<vector<REAL>> cellGrad1(cn,vector<REAL>(3,0));


    this->TriCellgradient(cn,_facevalue,_flux,_X,_bctype,cellGrad0);
    this->TriCellgradient(cn,_facevalue,_flux,_X,_bctype,cellGrad1);

//    for(int i=0 ; i<fn; i++)
//    {
//        cout << cellGrad0[i][0] << endl;
//    }

// facegradient
//    double faceGrad[fn][2];
    vector<vector<REAL>> faceGrad(fn,vector<REAL>(3,0));
    this->Facegradient(fn,cellGrad0,cellGrad1,faceGrad,_X);

//for(int i=0; i<fn; i++)
//{
//    cout << "id:" << i << endl;
//    cout << faceGrad[i][0] <<"  " << faceGrad[i][1] << endl;
//}
// correction
//    double value[fn][2];
    vector<vector<REAL>> value(fn,vector<REAL>(3,0));
    vector<REAL> Vector(3,0);
    vector<REAL> graddotvalue(fn,0);
    INT Pcid , faceID;
    REAL Gamma1 = 0.7;
    REAL X[cn];
    REAL _corr0[cn];
    REAL _corr1[cn];
    REAL S;
    S=0.0;


    for(int i=0; i<cn; i++)
    {
        Pcid = i ;
        X[i] = 0.0;
        _corr0[i] =0.0;
        _corr1[i] =0.0;
        for(int j=0; j<3; j++)
        {
            faceID = fvmesh->getFid(Pcid,j);
            cellID=fvmesh->getcellID(faceID);
            neicell1=fvmesh->getNeicell(faceID);
            fvmesh->getNormalVector(i,faceID,Vector);
            S = fvmesh->area(faceID);

            graddotvalue[faceID] =0.0;
           if(neicell1 == -1)
           {
               this->NormalVectorandZita2(faceID,cellID,Vector,value);
               fvmesh->dotproduct(faceGrad,value,faceID,graddotvalue);

//               cout << graddotvalue[faceID] << endl;
//               cout << S <<endl;
//               cout << Gamma << endl;
               _corr0[Pcid] += (-1) * graddotvalue[faceID] * S * Gamma * Gamma1;

//               cout << _corr0[Pcid] << endl;

           }

           else
           {
               if(cellID == i)
               {
                   this->NormalVectorandZita0(cellID,neicell1,Vector,faceID,value);
                   fvmesh->dotproduct(faceGrad,value,faceID,graddotvalue);

               }
               else
               {
                   this->NormalVectorandZita1(cellID,neicell1,Vector,faceID,value);
                   fvmesh->dotproduct(faceGrad,value,faceID,graddotvalue);
               }

               _corr1[Pcid] += graddotvalue[faceID] * S * Gamma * Gamma1; //change

//               cout << graddotvalue[faceID] << endl;
//               cout << S << endl;
//               cout << Gamma <<endl;
//               cout << _corr1[Pcid] << endl;
           }

         }

        _corr[Pcid] = (_corr0[Pcid] + _corr1[Pcid] );

    }

//        for(int i=0; i<cn; i++)
//        {
//            cout << i << endl;
//            cout << _corr[i] << endl;
//        }


}

void LaplaceSolver::QuadCalCorrectionValue(REAL *_X, vector<REAL> &_corr, vector<REAL> &_facevalue , vector<REAL> &_flux , vector<INT> &_bctype)
{
    INT fn,cn;
    fn = this->FaceN();
    cn = this->CellN();
    INT neicell0;
    INT neicell1;
    INT cellID;

// facePI
    this->FacePI(_X,_facevalue);

// cellgradient
//    double cellGrad0[cn][2];
//    double cellGrad1[cn][2];

    vector<vector<REAL>> cellGrad0(cn,vector<REAL>(3,0));
    vector<vector<REAL>> cellGrad1(cn,vector<REAL>(3,0));

    this->QuadCellgradient(cn,_facevalue,_flux,_X,_bctype,cellGrad0);
    this->QuadCellgradient(cn,_facevalue,_flux,_X,_bctype,cellGrad1);

// facegradient

    vector<vector<REAL>> faceGrad(fn,vector<REAL>(3,0));

    this->Facegradient(fn,cellGrad0,cellGrad1,faceGrad,_X);


// correction
//    double value[fn][2];
    vector<vector<REAL>> value(fn,vector<REAL>(3,0));
    vector<REAL> Vector(3,0);
//    double graddotvalue[fn];
    vector<REAL> graddotvalue(fn,0);
    INT Pcid , faceID;
    REAL Gamma1 = 0.7;
    REAL X[cn];
    REAL _corr0[cn];
    REAL _corr1[cn];
    REAL S;

    for(int i=0; i<cn; i++)
    {
        Pcid = i ;
        X[i] = 0.0;
        _corr0[i] =0.0;
        _corr1[i] =0.0;
        for(int j=0; j<4; j++)
        {
            faceID = fvmesh->getFid(Pcid,j);
            fvmesh->getNormalVector(i,faceID,Vector);
            neicell0=fvmesh->getcellID(faceID);
            neicell1=fvmesh->geomesh->facelist.getneicells1(faceID);
            S = fvmesh->geomesh->facelist.getarea(faceID);

            graddotvalue[faceID] =0.0;
           if(neicell1 == -1)
           {
               cellID = fvmesh->geomesh->facelist.getneicells0(faceID);

               this->NormalVectorandZita2(faceID,cellID,Vector,value);
               fvmesh->dotproduct(faceGrad,value,faceID,graddotvalue);

//               cout << graddotvalue[faceID] << endl;
//               cout << S <<endl;
//               cout << Gamma << endl;
               _corr0[Pcid] += (-1) * graddotvalue[faceID] * S * Gamma * Gamma1;

//               cout << _corr0[Pcid] << endl;

           }

           else
           {
               if(neicell0 == i)
               {
                   this->NormalVectorandZita0(neicell0,neicell1,Vector,faceID,value);
                   fvmesh->dotproduct(faceGrad,value,faceID,graddotvalue);

               }
               else
               {
                   this->NormalVectorandZita1(neicell0,neicell1,Vector,faceID,value);
                   fvmesh->dotproduct(faceGrad,value,faceID,graddotvalue);
               }

               _corr1[Pcid] +=  graddotvalue[faceID] * S * Gamma * Gamma1; //change

//               cout << graddotvalue[faceID] << endl;
//               cout << S << endl;
//               cout << Gamma <<endl;
//               cout << _corr1[Pcid] << endl;
           }
         }

        _corr[Pcid] = (_corr0[Pcid] + _corr1[Pcid] );

    }

//        for(int i=0; i<cn; i++)
//        {
//            cout << i << endl;
//            cout << _corr[i] << endl;
//        }


}

void LaplaceSolver::VectorAdd(vector<REAL> &_corr,REAL *_Bvector,REAL *_totalB)
{
    INT cn;
    cn = this->CellN();
    for(int i=0 ; i < cn ; i++)
    {
//        cout << i << endl;
        _totalB[i] =  (_corr[i] + _Bvector[i]);
//        cout << "id:" << i << endl;

//        cout << _totalB[i] << endl;

    }
//cout << endl ;
}

void LaplaceSolver::cellvelocity()
{
    int fn = fvmesh->geomesh->facelist.getfaceN();
    int cn = fvmesh->geomesh->cellist.getcellN();
//    double _faceV[fn][2];
    vector<vector<REAL>> _faceV(fn,vector<REAL>(2,0));
    double _cellV[cn][3];
    int Fid ;

    fvmesh->facevelocity0(_faceV);

    for(int i=0; i < cn; i++)
    {
        _cellV[i][0] = _cellV[i][1] = _cellV[i][2] = 0.0;
    }
    for(int i=0; i<cn; i++)
    {
        for(int j=0; j<4; j++)
        {
            Fid = fvmesh->geomesh->cellist.getcellfid(i,j);

//            cout << "faceid"<< Fid << endl;

//            cout << _faceV[Fid][0] << endl;
            _cellV[i][0] += (_faceV[Fid][0] / 4);
            _cellV[i][1] += (_faceV[Fid][1] / 4);
            _cellV[i][2] = 0;
        }
    }

    this->write2(_cellV);

    for(int i=0 ; i<cn; i++)
    {
//        cout <<"cellID : " << i << endl;
//       cout << _cellV[i][0] << " " <<_cellV[i][1] << "  " << _cellV[i][2] <<endl;
    }

}

void LaplaceSolver::write0(REAL *_X)
{
    INT NodeN;
    INT cellN;
    NodeN = fvmesh->geomesh->nodelist.getnodeN();
    cellN = fvmesh->geomesh->cellist.getcellN();

//    for(int i=0; i)
//    exportensight->EnSightGeoTriCellScalarResult("output1" , "scalar" ,NodeN,cn,, , X[] )


    FILE *in;
    in =fopen("laplace","w");

    for(int i=0; i<NodeN; i++)
    {
        fprintf(in,"10f" "10f" "10f" ,fvmesh->geomesh->nodelist.getcoordinateX(i) , fvmesh->geomesh->nodelist.getcoordinateY(i) , 0);
    }

    for(int i=0; i<cellN; i++)
    {
        fprintf(in, "10d" , "10d" , "10d", fvmesh->geomesh->cellist.getcellnids0(i), fvmesh->geomesh->cellist.getcellnids1(i) , fvmesh->geomesh->cellist.getcellnids2(i));
    }

    fclose(in);

    REAL coord[NodeN][3];

    for(int i=0; i<NodeN; i++)
    {
        coord[i][0] = fvmesh->geomesh->nodelist.getcoordinateX(i);
        coord[i][1] = fvmesh->geomesh->nodelist.getcoordinateY(i);
        coord[i][2] = 0;
    }

    INT tcoont[cellN][3];
    for(int i=0; i<cellN; i++)
    {
        tcoont[i][0] = fvmesh->geomesh->cellist.getcellnids0(i);
        tcoont[i][1] = fvmesh->geomesh->cellist.getcellnids1(i);
        tcoont[i][2] = fvmesh->geomesh->cellist.getcellnids2(i);
    }

    ExportEnSight exportensight;
//    exportensight.EnSightGeoTri(NodeN,cellN,coord,tcoont, "triangle");
    exportensight.EnSightGeoTriCellScalarResult("output2" ,"phi2",NodeN,cellN,coord,tcoont,_X);
}

void LaplaceSolver::write1(REAL *_X)
{
    INT NodeN;
    INT cellN;
    NodeN = fvmesh->geomesh->nodelist.getnodeN();
    cellN = fvmesh->geomesh->cellist.getcellN();

//    for(int i=0; i)
//    exportensight->EnSightGeoTriCellScalarResult("output1" , "scalar" ,NodeN,cn,, , X[] )


    FILE *in;
    in =fopen("laplace","w");

    for(int i=0; i<NodeN; i++)
    {
        fprintf(in,"10f" "10f" "10f" ,fvmesh->geomesh->nodelist.getcoordinateX(i) , fvmesh->geomesh->nodelist.getcoordinateY(i) , 0);
    }

    for(int i=0; i<cellN; i++)
    {
        fprintf(in, "10d" , "10d" , "10d", "10d", fvmesh->geomesh->cellist.getcellnids0(i), fvmesh->geomesh->cellist.getcellnids1(i) , fvmesh->geomesh->cellist.getcellnids2(i) , fvmesh->geomesh->cellist.getcellnids3(i));
    }

    fclose(in);

    REAL coord[NodeN][3];
    for(int i=0; i<NodeN; i++)
    {
        coord[i][0] = fvmesh->geomesh->nodelist.getcoordinateX(i);
        coord[i][1] = fvmesh->geomesh->nodelist.getcoordinateY(i);
        coord[i][2] = 0;

//        cout << fvmesh->geomesh->nodelist.getcoordinateX(i) << "  " << fvmesh->geomesh->nodelist.getcoordinateY(i) << endl;

    }

    INT tcoont[cellN][4];
    for(int i=0; i<cellN; i++)
    {
        tcoont[i][0] = fvmesh->geomesh->cellist.getcellnids0(i);
        tcoont[i][1] = fvmesh->geomesh->cellist.getcellnids1(i);
        tcoont[i][2] = fvmesh->geomesh->cellist.getcellnids2(i);
        tcoont[i][3] = fvmesh->geomesh->cellist.getcellnids3(i);

//        cout << fvmesh->geomesh->cellist.getcellnids0(i) << "  " << fvmesh->geomesh->cellist.getcellnids1(i) << "  " << fvmesh->geomesh->cellist.getcellnids2(i) << "  " << fvmesh->geomesh->cellist.getcellnids3(i) << endl;
    }

    ExportEnSight exportensight;
//    exportensight.EnSightGeoTri(NodeN,cellN,coord,tcoont, "triangle");
    exportensight.EnSightGeoQuadCellScalarResult("output1001" ,"Gamma1001",NodeN,cellN,coord,tcoont,_X);
}

void LaplaceSolver::write2(REAL (*_V)[3])
{
    int NodeN;
    int cellN;
    NodeN = fvmesh->geomesh->nodelist.getnodeN();
    cellN = fvmesh->geomesh->cellist.getcellN();

//    for(int i=0; i)
//    exportensight->EnSightGeoTriCellScalarResult("output1" , "scalar" ,NodeN,cn,, , X[] )


    FILE *in;
    in =fopen("laplace","w");

    for(int i=0; i<NodeN; i++)
    {
        fprintf(in,"10f" "10f" "10f" ,fvmesh->geomesh->nodelist.getcoordinateX(i) , fvmesh->geomesh->nodelist.getcoordinateY(i) , 0);
    }

    for(int i=0; i<cellN; i++)
    {
        fprintf(in, "10d" , "10d" , "10d", "10d", fvmesh->geomesh->cellist.getcellnids0(i), fvmesh->geomesh->cellist.getcellnids1(i) , fvmesh->geomesh->cellist.getcellnids2(i) , fvmesh->geomesh->cellist.getcellnids3(i));
    }

    fclose(in);

    REAL coord[NodeN][3];
    for(int i=0; i<NodeN; i++)
    {
        coord[i][0] = fvmesh->geomesh->nodelist.getcoordinateX(i);
        coord[i][1] = fvmesh->geomesh->nodelist.getcoordinateY(i);
        coord[i][2] = 0.0;

//        cout << fvmesh->geomesh->nodelist.getcoordinateX(i) << "  " << fvmesh->geomesh->nodelist.getcoordinateY(i) << endl;

    }

    INT tcoont[cellN][4];
    for(int i=0; i<cellN; i++)
    {
        tcoont[i][0] = fvmesh->geomesh->cellist.getcellnids0(i);
        tcoont[i][1] = fvmesh->geomesh->cellist.getcellnids1(i);
        tcoont[i][2] = fvmesh->geomesh->cellist.getcellnids2(i);
        tcoont[i][3] = fvmesh->geomesh->cellist.getcellnids3(i);

//        cout << fvmesh->geomesh->cellist.getcellnids0(i) << "  " << fvmesh->geomesh->cellist.getcellnids1(i) << "  " << fvmesh->geomesh->cellist.getcellnids2(i) << "  " << fvmesh->geomesh->cellist.getcellnids3(i) << endl;
    }

    ExportEnSight exportensight;
//    exportensight.EnSightGeoTri(NodeN,cellN,coord,tcoont, "triangle");
    exportensight.EnSightGeoQuadCellVectorResult("output4.0" ,"phi4.0",NodeN,cellN,coord,tcoont,_V);
}

