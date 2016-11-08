#include <ilcplex/ilocplex.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <string>
ILOSTLBEGIN

//-------------Global Variables--------------
int i,t,z,j,k,n,l;
const int tmax=4;
const int zmax=6;
const int imax=1;
const int kmax=2;
const int jmax=3;
const int capmax=100000;
const int capmin=0;
const int m=100000;
double E[imax][tmax];
double ae[imax][jmax][tmax];
double S[kmax][tmax];
double as[kmax][jmax][tmax];
double initial[zmax][jmax];

double epsilon=0.001;
double UpperBound=100;
double LowerBound=0;
double UpperBoundGlobal=100;
double fraction=0.90;
long double duration;  // tracks time
int start=0,BDFeasCuts=0,BDOptCuts=0;

double CiztValue[imax][zmax][tmax];
double DkztValue[kmax][zmax][tmax];
double FzjtValue[zmax][jmax][tmax];
double SCiztValue[imax][zmax][tmax];
double SDkztValue[kmax][zmax][tmax];
double XizjtValue[imax][zmax][jmax][tmax];
double YzkjtValue[zmax][kmax][jmax][tmax];
double IzjtValue[zmax][jmax][tmax];
double ThetaValue=0;

double OptimalCiztValue[imax][zmax][tmax];
double OptimalDkztValue[kmax][zmax][tmax];
double OptimalFzjtValue[zmax][jmax][tmax];
double OptimalSCiztValue[imax][zmax][tmax];
double OptimalSDkztValue[kmax][zmax][tmax];
double OptimalXizjtValue[imax][zmax][jmax][tmax];
double OptimalYzkjtValue[zmax][kmax][jmax][tmax];
double OptimalIzjtValue[zmax][jmax][tmax];
double OptimalThetaValue=0;

double OptimalOriginalObjFunction=0;
double OptimalMasterObjFunction=0;
double OptimalSlaveObjFunction=0;



//--------Declare the environment of CPLEX----------------
IloEnv env;
//--------Construct models----------------
IloModel modelSlave1 (env);
IloModel modelMaster (env);	
//--------Construct Matrices----------------
typedef IloArray<IloNumArray> IloNumMatrix2x2;
typedef IloArray<IloNumMatrix2x2> IloNumMatrix3x3; 
typedef IloArray<IloNumMatrix3x3> IloNumMatrix4x4;

typedef IloArray<IloNumVarArray> IloNumVarMatrix2x2;
typedef IloArray<IloNumVarMatrix2x2> IloNumVarMatrix3x3;
typedef IloArray<IloNumVarMatrix3x3> IloNumVarMatrix4x4;

typedef IloArray<IloRangeArray> IloRangeMatrix2x2;
typedef IloArray<IloRangeMatrix2x2> IloRangeMatrix3x3;
typedef IloArray<IloRangeMatrix3x3> IloRangeMatrix4x4;

//------Declare Decision Variables----------
IloNumVarMatrix3x3 Cizt(env,0);
IloNumVarMatrix3x3 Dkzt(env,0);
IloNumVarMatrix3x3 Fzjt(env,0);
IloNumVarArray Zn(env,0);
IloNumVarMatrix4x4 Xizjt(env,0);
IloNumVarMatrix4x4 Yzkjt(env,0);
IloNumVarMatrix3x3 Izjt(env,0);
IloNumVarMatrix3x3 SCizt(env,0);
IloNumVarMatrix3x3 SDkzt(env,0);

//--------Declare Slave constraints-------------
IloRangeMatrix3x3 SumXijt(env,0);
IloRangeMatrix3x3 DSumXijt(env,0);
IloRangeMatrix3x3 SumYkjt(env,0);
IloRangeMatrix3x3 DSumYkjt(env,0);
IloRangeMatrix3x3 CTIzjt(env,0);
IloRangeMatrix3x3 DCTIzjt(env,0); 
IloRangeMatrix2x2 Sum_Izt(env,0);
IloRangeMatrix3x3 CT1Fonctionement_Cizt(env,0);
IloRangeMatrix3x3 CT2Fonctionement_Cizt(env,0);
IloRangeMatrix3x3 CT1Fonctionement_Dkzt(env,0);
IloRangeMatrix3x3 CT2Fonctionement_Dkzt(env,0);
IloRangeMatrix3x3 CT1Melzjt(env,0);
IloRangeMatrix3x3 CT2Melzjt(env,0);
IloRangeMatrix3x3 SC2_Cizt(env,0);
IloRangeMatrix3x3 SD2_Dkzt(env,0);
IloRangeMatrix3x3 SC_Cizt(env,0);
IloRangeMatrix3x3 SD_Dkzt(env,0);

//--------Declare dual variables of each constraint----------------

//double valsDualCT3Melzt[zmax][tmax];

double S2valsDualSumXijt[imax][jmax][tmax];
double S2valsDualSumYkjt[kmax][jmax][tmax];
double S2valsDualCTIzjt[zmax][jmax][tmax];

double S22valsDualSumXijt[imax][jmax][tmax];
double S22valsDualSumYkjt[kmax][jmax][tmax];
double S22valsDualCTIzjt[zmax][jmax][tmax];

double S2valsDualSum_Izt[zmax][tmax];

double S2valsDualCT1Fonctionement_Cizt[imax][zmax][tmax];
double S2valsDualCT2Fonctionement_Cizt[imax][zmax][tmax];
double S2valsDualCT1Fonctionement_Dkzt[kmax][zmax][tmax];
double S2valsDualCT2Fonctionement_Dkzt[kmax][zmax][tmax];

double S2valsDualCT1Melzjt[zmax][jmax][tmax];
double S2valsDualCT2Melzjt[zmax][jmax][tmax];

double S2valsDualSC2_Cizt[imax][zmax][tmax];
double S2valsDualSD2_Dkzt[kmax][zmax][tmax];
double S2valsDualSC_Cizt[imax][zmax][tmax];
double S2valsDualSD_Dkzt[kmax][zmax][tmax];

//----------What does IloNum mean?---------------

IloNum valsDualRangeSumXijt;
IloNum valsDualRangeSumYkjt;
IloNum valsDualRangeCTIzjt;
IloNum valsDualRangeSum_Izt;
IloNum valsDualRangeCT1Fonctionement_Cizt;
IloNum valsDualRangeCT2Fonctionement_Cizt;
IloNum valsDualRangeCT1Fonctionement_Dkzt;
IloNum valsDualRangeCT2Fonctionement_Dkzt;

IloNum valsDualRangeCT1Melzjt;
IloNum valsDualRangeCT2Melzjt;
IloNum valsDualRangeSC2_Cizt;
IloNum valsDualRangeSD2_Dkzt;
IloNum valsDualRangeSC_Cizt;
IloNum valsDualRangeSD_Dkzt;

IloNumArray FeasvalsDualRangeSumXijt(env,0);

vector <double> LowerBoundArray;
vector <double> UpperBoundArray;
vector <double> UpperBoundGlobalArray;
vector <double> dTy;
vector <double> zCurrent;
vector <double> cTx;
vector <double> BestSlaveObjSoFar;
vector <double> Time;

typedef struct treenode_tag {
	double  lpbound;  // LP bound
	IloModel  lp;     // ptr to master
	IloModel  lp_cg;   // ptr to colgen
	treenode_tag  *nextnode;  // link to next node in tree
} treenode;

treenode_tag *BBTreeList;

void Found_Error(char *name)
{
	printf("%s failed, exiting...\n", name);
	printf("Press return to continue...\n");
    getchar();
}
int load_data(){
//-------------------Declare Data of the problem--------------------

for (i=0;i<imax;i++){
	for (t=0;t<tmax;t++){
		E[i][t]=0;
	}
}
for (i=0;i<imax;i++){
	for (j=0;j<jmax;j++){
		for (t=0;t<tmax;t++){
			ae[i][j][t]=0;
		}
	}
}

for (k=0;k<kmax;k++){
	for (j=0;j<jmax;j++){
		for(t=0;t<tmax;t++){
			as[k][j][t]=0;
		}
	}
}

for (k=0;k<kmax;k++){
	for(t=0;t<tmax;t++){
			S[k][t]=0;
	}
}

E[0][3]=120000;
//E[0][5]=110000;
//E[0][7]=90000;

ae[0][0][3]=1;
//ae[0][1][5]=1;
//ae[0][0][7]=1;

//DEMANDE
//CDU I
S[0][0]=0;
S[0][1]=30000;
S[0][2]=27750;
S[0][3]=9750;
/*
S[0][4]=27750;
S[0][5]=9750;
S[0][6]=35000;
S[0][7]=11250;
S[0][8]=31875;
S[0][9]=31875;
*/


//PourCentage Entree
//PourCentage Sorti CDU I 
as[0][0][0]=0;
as[0][1][0]=0;
as[0][2][0]=0;
//as[0][3][0]=0;

as[0][0][1]=0.50;
as[0][1][1]=0;
as[0][2][1]=0.50;
//as[0][3][1]=0;

as[0][0][2]=1;
as[0][1][2]=0;
as[0][2][2]=0;
//as[0][3][2]=0;

as[0][0][3]=1;
as[0][1][3]=0;
as[0][2][3]=0;
//as[0][3][3]=0;
/*
as[0][0][4]=1;
as[0][1][4]=0;
as[0][2][4]=0;
as[0][3][4]=0;

as[0][0][5]=1;
as[0][1][5]=0;
as[0][2][5]=0;
as[0][3][5]=0;

as[0][0][6]=0.28;
as[0][1][6]=0.28;
as[0][2][6]=0.44;
as[0][3][6]=0;

as[0][0][7]=0.53;
as[0][1][7]=0.47;
as[0][2][7]=0;
as[0][3][7]=0;

as[0][0][8]=0.53;
as[0][1][8]=0.47;
as[0][2][8]=0;
as[0][3][8]=0;

as[0][0][9]=0.53;
as[0][1][9]=0.47;
as[0][2][9]=0;
as[0][3][9]=0;
*/
//CDU II
S[1][0]=0;
S[1][1]=60000;
S[1][2]=30000;
S[1][3]=15600;
/*
S[1][4]=44400;
S[1][5]=60000;
S[1][6]=50000;
S[1][7]=40000;
S[1][8]=50000;
S[1][9]=60000;
*/
//PourCentage Sorti CDU II
as[1][0][0]=0;
as[1][1][0]=0;
as[1][2][0]=0;
//as[1][3][0]=0;

as[1][0][1]=0.41;
as[1][1][1]=0;
as[1][2][1]=0.59;
//as[1][3][1]=0;

as[1][0][2]=1;
as[1][1][2]=0;
as[1][2][2]=0;
//as[1][3][2]=0;

as[1][0][3]=0;
as[1][1][3]=0;
as[1][2][3]=1;
//as[1][3][3]=0;
/*
as[1][0][4]=0;
as[1][1][4]=0;
as[1][2][4]=1;
as[1][3][4]=0;

as[1][0][5]=1;
as[1][1][5]=0;
as[1][2][5]=0;
as[1][3][5]=0;

as[1][0][6]=0;
as[1][1][6]=0.60;
as[1][2][6]=0.40;
as[1][3][6]=0;

as[1][0][7]=0.38;
as[1][1][7]=0.62;
as[1][2][7]=0;
as[1][3][7]=0;

as[1][0][8]=0;
as[1][1][8]=0.20;
as[1][2][8]=0.80;
as[1][3][8]=0;

as[1][0][9]=1;
as[1][1][9]=0;
as[1][2][9]=0;
as[1][3][9]=0;
*/
//Initial Values


initial[0][0]=40000;
initial[0][1]=0;
initial[0][2]=0;
//initial[0][3]=0;

initial[1][0]=80000;
initial[1][1]=0;
initial[1][2]=0;
//initial[1][3]=0;

initial[2][0]=0;
initial[2][1]=0;
initial[2][2]=95000;
//initial[2][3]=0;

initial[3][0]=0;
initial[3][1]=0;
initial[3][2]=100000;
//initial[3][3]=0;

initial[4][0]=0;
initial[4][1]=0;
initial[4][2]=0;
//initial[4][3]=0;

initial[5][0]=0;
initial[5][1]=0;
initial[5][2]=0;
//initial[5][3]=0;

// DATA SOLUTION C,D 21//////////////////
for (i=0;i<imax;i++){
	for (z=0;z<zmax;z++){
		for (t=0;t<tmax;t++){
			CiztValue[i][z][t]=0;
			SCiztValue[i][z][t]=0;
			for (j=0;j<jmax;j++){
				XizjtValue[i][z][j][t];
				OptimalXizjtValue[i][z][j][t];
			}
		}
	}
}

for (k=0;k<kmax;k++){
	for (z=0;z<zmax;z++){
		for (t=0;t<tmax;t++){
			DkztValue[k][z][t]=0;
			SDkztValue[k][z][t]=0;
			for (j=0;j<jmax;j++){
				YzkjtValue[z][k][j][t];
				OptimalYzkjtValue[z][k][j][t];
			}
		}
	}
}

for (z=0;z<zmax;z++){
	for (j=0;j<jmax;j++){
		for (t=0;t<tmax;t++){
			FzjtValue[z][j][t]=0;
			IzjtValue[z][j][t]=0;
			OptimalFzjtValue[z][j][t]=0;
			OptimalIzjtValue[z][j][t]=0;
		}
	}
}

// FIN des DATA///////////////////////////
return 0;
}
int do_master(){

//double sum1[imax][tmax];
double sum2[kmax][tmax];
//double tot1[tmax];
//double tot2[tmax];			
double sum[jmax][tmax];	

//------------------------------------------------------------------------------
//---------------------------------- MASTER ------------------------------------
//------------------------------------------------------------------------------
//----------------------------- Master Variable --------------------------------
//-------------- Variable de Decision C ---------------------------------------

for (i=0;i<imax;i++){
	IloNumVarMatrix2x2 Czt(env,0);
	for (z=0;z<zmax;z++){
		IloNumVarArray Ct(env,0);
		for (t=0;t<tmax;t++){
			char Chargement[70];  
			sprintf(Chargement,"Cizt(i%d,z%d,t%d)",i,z,t); 
        		IloNumVar C(env,0,1,ILOINT,Chargement); 
        		Ct.add(C);
		}
		Czt.add(Ct);
	}
	Cizt.add(Czt);
}
//-------------- Variable de Decision D ---------------------------------------

for (k=0;k<kmax;k++){
	IloNumVarMatrix2x2 Dzt(env,0);
	for (z=0;z<zmax;z++){
		IloNumVarArray Dt(env,0);
		for (t=0;t<tmax;t++){
			char Dechargement[70];  
			sprintf(Dechargement,"Dkzt(k%d,z%d,t%d)",k,z,t); 
            		IloNumVar D(env,0,1,ILOINT,Dechargement); 
            		Dt.add(D);
		}
		Dzt.add(Dt);
	}
	Dkzt.add(Dzt);
}
//-------------- Variable de Decision f ---------------------------------------

for (z=0;z<zmax;z++){
	IloNumVarMatrix2x2 Fjt(env,0);
	for (j=0;j<jmax;j++){
		IloNumVarArray Ft(env,0);
		for (t=0;t<tmax;t++){
			char variableF[70];  
			sprintf(variableF,"Fzjt(z%d,j%d,t%d)",z,j,t); 
            		IloNumVar F(env,0,1,ILOINT,variableF); 
            		Ft.add(F);
	  	}
		Fjt.add(Ft);
	}
	Fzjt.add(Fjt);
}
//--------------------------- Variable de Decision Z ---------------------------


for (n=0;n<1;n++){
	char Theta[70];  
	sprintf(Theta,"Zn(n%d)",n); 
	IloNumVar Z(env,0,IloInfinity,ILOFLOAT,Theta); 
	Zn.add(Z);
}
			
//-----------------------------Finish of Master Variables --------------------------------

//-----------------------------------------------------------------------------
//-------------------------Start of Master Constraints-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------- Contrainte de melange CT3 -------------------------
IloRangeMatrix2x2 CT3Melzt(env,0);
for(z=0;z<zmax;z++){
	IloRangeArray CT3Melt(env,0);
	for(t=0;t<tmax;t++){
		IloExpr expr(env,0);
		for(j=0;j<jmax;j++){
			expr+=Fzjt[z][j][t];
		}
		char CT3Melange[60];
		sprintf(CT3Melange,"CT3Melzjt(z%d,t%d)",z,t);
		double LBCT3Melzt=0,UBCT3Melzt=1;
		IloRange CT3Mel(env,LBCT3Melzt,expr,UBCT3Melzt,CT3Melange);
		modelMaster.add(CT3Mel);
		CT3Melt.add(CT3Mel);
		expr.end();
	}
CT3Melzt.add(CT3Melt);
}

//------------------------------- Chargement CT3 ---------------------------------
IloRangeMatrix2x2 CT3C_ou_Dzt(env,0);
for (z=0;z<zmax;z++){
	IloRangeArray CT3C_ou_Dt(env,0);
	for (t=0;t<tmax;t++){
		for (i=0;i<imax;i++){
			for (k=0;k<kmax;k++){
				IloExpr expr(env,0);
				expr+=Cizt[i][z][t]+Dkzt[k][z][t];
				char Chargement_ou_Dechargement[60];
				sprintf(Chargement_ou_Dechargement,"CT3C_ou_Dzt(z%d,t%d)",z,t);
				double LBCT3C_ou_Dzt=0,UBCT3C_ou_Dzt=1;
				IloRange CT3C_ou_D(env,LBCT3C_ou_Dzt,expr,UBCT3C_ou_Dzt,Chargement_ou_Dechargement);
				modelMaster.add(CT3C_ou_D);
				CT3C_ou_Dt.add(CT3C_ou_D);
				expr.end();
			}
		}
	}
	CT3C_ou_Dzt.add(CT3C_ou_Dt);
}

//------------------------------- INITIAL VALUES (t=0) ---------------------------------
//------------------------------- z tank contains j type at t=0 ---------------------------------
IloRangeMatrix2x2 SupFzj0(env,0);
for (z=0;z<zmax;z++){
	IloRangeArray SupFj0(env,0);
	for (j=0;j<jmax;j++){
	
	if (z==0){
		if(j==0){//------ 0 tank contains 0 type at t=0 ------------
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=1,UB=1;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
		else{
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=0,UB=0;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
	}

	
	if (z==1){
		if(j==0){//------ 1 tank contains 0 type at t=0 ------------
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=1,UB=1;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
		else{
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=0,UB=0;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
	}
	
	
	if (z==2){
		if(j==2){//------ 2 tank contains 2 type at t=0 ------------
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=1,UB=1;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
		else{
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=0,UB=0;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
	}

	if (z==3){
		if(j==2){//------ 3 tank contains 2 type at t=0 ------------
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=1,UB=1;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}
	
		else{
			IloExpr expr(env,0);
			expr=Fzjt[z][j][0];
			char ConSupFzj0[60];
			sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
			double LB=0,UB=0;
			IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
			modelMaster.add(SupF0);
			SupFj0.add(SupF0);
			expr.end();
		}

	}



	if (z==4){//------ 4 tank is empty at t=0 ------------
		IloExpr expr(env,0);
		expr=Fzjt[z][j][0];
		char ConSupFzj0[60];
		sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
		double LB=0,UB=0;
		IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
		modelMaster.add(SupF0);
		SupFj0.add(SupF0);
		expr.end();
	}

	if (z==5){//------ 5 tank is empty at t=0 ------------
		IloExpr expr(env,0);
		expr=Fzjt[z][j][0];
		char ConSupFzj0[60];
		sprintf(ConSupFzj0,"ConSupFzj0(z%d,j%d)",z,j);
		double LB=0,UB=0;
		IloRange SupF0(env,LB,expr,UB,ConSupFzj0);
		modelMaster.add(SupF0);
		SupFj0.add(SupF0);
		expr.end();
	}

	}
	SupFzj0.add(SupFj0);
}

//--------------- No loading taking place at t=0 ------------
IloRangeMatrix2x2 SupCiz0(env,0);
for (i=0;i<imax;i++){
	IloRangeArray SupCz0(env,0);
	for (z=0;z<zmax;z++){
		IloExpr expr(env,0);
		expr=Cizt[i][z][0];
		char ConSupCiz0[60];
		sprintf(ConSupCiz0,"ConSupCiz0(i%d,z%d)",i,z);
		double LB=0,UB=0;
		IloRange SupC0(env,LB,expr,UB,ConSupCiz0);
		modelMaster.add(SupC0);
		SupCz0.add(SupC0);
		expr.end();
	}
SupCiz0.add(SupCz0);
}

//--------------- No unloading taking place at t=0 ------------
IloRangeMatrix2x2 SupDkz0(env,0);
for (k=0;k<kmax;k++){
	IloRangeArray SupDz0(env,0);
	for (z=0;z<zmax;z++){
		IloExpr expr(env,0);
		expr=Dkzt[k][z][0];
		char ConSupDkz0[60];
		sprintf(ConSupDkz0,"ConSupDkz0(k%d,z%d)",k,z);
		double LB=0,UB=0;
		IloRange SupD0(env,LB,expr,UB,ConSupDkz0);
		modelMaster.add(SupD0);
		SupDz0.add(SupD0);
		expr.end();
	}
SupDkz0.add(SupDz0);
}

//-------------------------------Finish of INITIAL VALUES (t=0) ---------------------------------

//Contrainte de feasabilitι////////////////////////////

IloRangeMatrix2x2 Con3W1it(env,0);
for (i=0;i<imax;i++){
	IloRangeArray Con3W1t(env,0);
	for (t=0;t<tmax;t++){
		IloExpr expr(env,0);
		for (z=0;z<zmax;z++){
			expr+=Cizt[i][z][t];
		}
		char CT3W1[60];
		sprintf(CT3W1,"CT3W1it(i%d,t%d)",i,t);
		double LB=E[i][t]/100000,UB=IloInfinity;
		IloRange Con3W1(env,LB,expr,UB,CT3W1);
		modelMaster.add(Con3W1);
		Con3W1t.add(Con3W1);
		expr.end();
	}
Con3W1it.add(Con3W1t);
}


IloRangeMatrix2x2 Con5W2kt(env,0);
for (k=0;k<kmax;k++){
	IloRangeArray Con5W2t(env,0);
	for (t=0;t<tmax;t++){
		IloExpr expr(env,0);
		for (z=0;z<zmax;z++){
			expr+=Dkzt[k][z][t];
		}
		char CT5W2[60];
		sprintf(CT5W2,"CT5W2kt(k%d,t%d)",k,t);
		double LB=S[k][t]/100000,UB=IloInfinity;
		IloRange Con5W2(env,LB,expr,UB,CT5W2);
		modelMaster.add(Con5W2);
		Con5W2t.add(Con5W2);
		expr.end();
	}
Con5W2kt.add(Con5W2t);
}



//---------------NO VALID INEQUALITIES-------------------

//-----------------------------------------------------------------------------
//-------------------------Finish of Master Constraints-----------------------------------------


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------------Fonction Objectif de master probleme--------------------------
//------------------------------------------------------------------------------
IloExpr expr1(env);

for (i=0;i<imax;i++){
	for (z=0;z<zmax;z++){
		for (t=0;t<tmax;t++){
			expr1+=Cizt[i][z][t];
		}
	}
}

for (k=0;k<kmax;k++){
	for (z=0;z<zmax;z++){
		for (t=0;t<tmax;t++){
			expr1+=Dkzt[k][z][t];
		}
	}
}
for (n=0;n<1;n++){
	expr1+=Zn[n];
}

modelMaster.add(IloMinimize(env, expr1));
expr1.end();		

return 0;
}
int do_slave(){

//--------------------------- PRIMAL SLAVE PROBLEM ----------------------------------
//------------------------------------------------------------------------------
//--------------------------- Slave Primal Variables ---------------------------
//--------------------------- Variable de Decision X ---------------------------

for (i=0;i<imax;i++){
	IloNumVarMatrix3x3 Xzjt(env,0);
	for (z=0;z<zmax;z++){
		IloNumVarMatrix2x2 Xjt(env,0);
		for (j=0;j<jmax;j++){
			IloNumVarArray Xt(env,0);
			for (t=0;t<tmax;t++){
				char Quantite_Charge[70];  
				sprintf(Quantite_Charge,"Xizjt(i%d,z%d,j%d,t%d)",i,z,j,t); 
                		IloNumVar X(env,0,IloInfinity,ILOFLOAT,Quantite_Charge); 
                		Xt.add(X);
			}
			Xjt.add(Xt);
		}
		Xzjt.add(Xjt);
	}
	Xizjt.add(Xzjt);
}
//--------------------------- Variable de Decision Y --------------------------

for (z=0;z<zmax;z++){
	IloNumVarMatrix3x3 Ykjt(env,0);
	for (k=0;k<kmax;k++){
		IloNumVarMatrix2x2 Yjt(env,0);
		for (j=0;j<jmax;j++){
			IloNumVarArray Yt(env,0);
			for (t=0;t<tmax;t++){
				char Quantite_Decharge[70];  
				sprintf(Quantite_Decharge,"Yzkjt(k%d,z%d,j%d,t%d)",k,z,j,t); 
				IloNumVar Y(env,0,IloInfinity,ILOFLOAT,Quantite_Decharge); 
                		Yt.add(Y);
			}
			Yjt.add(Yt);
		}
		Ykjt.add(Yjt);
	}
	Yzkjt.add(Ykjt);
}
//--------------------------- Variable de Decision I --------------------------

for (z=0;z<zmax;z++){
	IloNumVarMatrix2x2 Ijt(env,0);
	for (j=0;j<jmax;j++){
		IloNumVarArray It(env,0);
		for (t=0;t<tmax;t++){
			char Inventory[70];  
			sprintf(Inventory,"Izjt(z%d,j%d,t%d)",z,j,t); 
          		IloNumVar I(env,0,IloInfinity,ILOFLOAT,Inventory); 
        		It.add(I);
	}
	Ijt.add(It);
	}
Izjt.add(Ijt);
}

//-------------- Variable de Decision SC ---------------------------------------

for (i=0;i<imax;i++){
	IloNumVarMatrix2x2 SCzt(env,0);
	for (z=0;z<zmax;z++){
		IloNumVarArray SCt(env,0);
		for (t=0;t<tmax;t++){
			char SChargement[70];  
			sprintf(SChargement,"SCizt(i%d,z%d,t%d)",i,z,t); 
			IloNumVar SC(env,0,IloInfinity,ILOFLOAT,SChargement); 
            		SCt.add(SC);
		}
		SCzt.add(SCt);
	}
	SCizt.add(SCzt);
}
//-------------- Variable de Decision SD ---------------------------------------

for (k=0;k<kmax;k++){
	IloNumVarMatrix2x2 SDzt(env,0);
	for (z=0;z<zmax;z++){
		IloNumVarArray SDt(env,0);
		for (t=0;t<tmax;t++){
			char SDechargement[70];  
			sprintf(SDechargement,"SDkzt(k%d,z%d,t%d)",k,z,t); 
            		IloNumVar SD(env,0,IloInfinity,ILOFLOAT,SDechargement); 
            		SDt.add(SD);
		}
		SDzt.add(SDt);
	}
	SDkzt.add(SDzt);
}

//-----------------  FIN VARIABLE DE DECISION PRIMAL SLAVE 1 ------------------
//-----------------------------------------------------------------------------
//------------------------- Slave Primal  Constraintes ------------------------
//-----------------------------------------------------------------------------
//------------------------------ Quantite Charge ------------------------------

for (i=0;i<imax;i++){
	IloRangeMatrix2x2 SumXjt(env,0);
	for (j=0;j<jmax;j++){
		IloRangeArray SumXt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (z=0;z<zmax;z++){
				expr+=Xizjt[i][z][j][t];
			}
			char Quantite_Charge[60];
			sprintf(Quantite_Charge,"SumXizjt(i%d,j%d,t%d)",i,j,t);
			double LBSumXijt=E[i][t]*ae[i][j][t],UBSumXijt=IloInfinity;
			IloRange SumX(env,LBSumXijt,expr,UBSumXijt,Quantite_Charge);
			modelSlave1.add(SumX);
			SumXt.add(SumX);
			expr.end();
		}
	
		SumXjt.add(SumXt);
	}
	
	SumXijt.add(SumXjt);
}

for (i=0;i<imax;i++){
	IloRangeMatrix2x2 DSumXjt(env,0);
	for (j=0;j<jmax;j++){
		IloRangeArray DSumXt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (z=0;z<zmax;z++){
				expr+=-Xizjt[i][z][j][t];
			}
			char DQuantite_Charge[60];
			sprintf(DQuantite_Charge,"DSumXizjt(i%d,j%d,t%d)",i,j,t);
			double DLBSumXijt=-E[i][t]*ae[i][j][t],DUBSumXijt=IloInfinity;
			IloRange DSumX(env,DLBSumXijt,expr,DUBSumXijt,DQuantite_Charge);
			modelSlave1.add(DSumX);
			DSumXt.add(DSumX);
			expr.end();
		}
	
		DSumXjt.add(DSumXt);
	}
	
	DSumXijt.add(DSumXjt);
}

//-----------------------------------------------------------------------------
//-------------------------- Quantite Decharge --------------------------------

for (k=0;k<kmax;k++){
	IloRangeMatrix2x2 SumYjt(env,0);
	for (j=0;j<jmax;j++){
		IloRangeArray SumYt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env);
			for (z=0;z<zmax;z++){
				expr+=Yzkjt[z][k][j][t];
			}
			char Quantite_Decharge[60];
			sprintf(Quantite_Decharge,"SumYzkjt(k%d,j%d,t%d)",k,j,t);
			double LBSumYkjt=S[k][t]*as[k][j][t],UBSumYkjt=IloInfinity;
			IloRange SumY(env,LBSumYkjt,expr,UBSumYkjt,Quantite_Decharge);
			modelSlave1.add(SumY);
			SumYt.add(SumY);
			expr.end();
		}
	SumYjt.add(SumYt);
	}
SumYkjt.add(SumYjt);
}



for (k=0;k<kmax;k++){
	IloRangeMatrix2x2 DSumYjt(env,0);
	for (j=0;j<jmax;j++){
		IloRangeArray DSumYt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env);
			for (z=0;z<zmax;z++){
				expr+=-Yzkjt[z][k][j][t];
			}
			char DQuantite_Decharge[60];
			sprintf(DQuantite_Decharge,"DSumYzkjt(k%d,j%d,t%d)",k,j,t);
			double DLBSumYkjt=-S[k][t]*as[k][j][t],DUBSumYkjt=IloInfinity;
			IloRange DSumY(env,DLBSumYkjt,expr,DUBSumYkjt,DQuantite_Decharge);
			modelSlave1.add(DSumY);
			DSumYt.add(DSumY);
			expr.end();
		}
	DSumYjt.add(DSumYt);
	}
DSumYkjt.add(DSumYjt);
}



//-----------------------------------------------------------------------------
//---------------------------- Equation d'equilibre ---------------------------
 
for (z=0;z<zmax;z++){
	IloRangeMatrix2x2 CTIjt(env,0);
	for (j=0;j<jmax;j++){
		IloRangeArray CTIt(env,0);
		for (t=0;t<tmax;t++){
			if(t==0){
				IloExpr expr(env);
				for (i=0;i<imax;i++){
					expr-=Xizjt[i][z][j][t];
				}
				for (k=0;k<kmax;k++){
					expr+=Yzkjt[z][k][j][t];
				}
				expr+=Izjt[z][j][t];
				char Equoition_equilibre[60];
				sprintf(Equoition_equilibre,"CTIzjt(z%d,j%d,t%d)",z,j,t);
				double LBCTIzjt=initial[z][j],UBCTIzjt=IloInfinity;
				IloRange CTI(env,LBCTIzjt,expr,UBCTIzjt,Equoition_equilibre);
				modelSlave1.add(CTI);
				CTIt.add(CTI);
				expr.end();
		
			}
	
			else{
			
				IloExpr expr(env);
			
				for (i=0;i<imax;i++){
				expr-=Xizjt[i][z][j][t];
				}
				for (k=0;k<kmax;k++){
				expr+=Yzkjt[z][k][j][t];
				}
				expr+=Izjt[z][j][t]-Izjt[z][j][t-1];
				char Equoition_equilibre[60];
				sprintf(Equoition_equilibre,"CTIzjt(z%d,j%d,t%d)",z,j,t);
				double LBCTIzjt=0,UBCTIzjt=IloInfinity;
				IloRange CTI(env,LBCTIzjt,expr,UBCTIzjt,Equoition_equilibre);
				modelSlave1.add(CTI);
				CTIt.add(CTI);
				expr.end();
			}
		
		}
	CTIjt.add(CTIt);
	}
CTIzjt.add(CTIjt);
}


for (z=0;z<zmax;z++){
	IloRangeMatrix2x2 DCTIjt(env,0);
	for (j=0;j<jmax;j++){
		IloRangeArray DCTIt(env,0);
		for (t=0;t<tmax;t++){
			if(t==0){
				IloExpr expr(env);
				for (i=0;i<imax;i++){
					expr+=Xizjt[i][z][j][t];
				}
				for (k=0;k<kmax;k++){
					expr+=-Yzkjt[z][k][j][t];
				}
				expr+=-Izjt[z][j][t];
				char DEquoition_equilibre[60];
				sprintf(DEquoition_equilibre,"DCTIzjt(z%d,j%d,t%d)",z,j,t);
				double DLBCTIzjt=-initial[z][j],DUBCTIzjt=IloInfinity;
				IloRange DCTI(env,DLBCTIzjt,expr,DUBCTIzjt,DEquoition_equilibre);
				modelSlave1.add(DCTI);
				DCTIt.add(DCTI);
				expr.end();
		
			}
	
			else{
			
				IloExpr expr(env);
			
				for (i=0;i<imax;i++){
				expr+=Xizjt[i][z][j][t];
				}
				for (k=0;k<kmax;k++){
				expr+=-Yzkjt[z][k][j][t];
				}
				expr+=-Izjt[z][j][t]+Izjt[z][j][t-1];
				char DEquoition_equilibre[60];
				sprintf(DEquoition_equilibre,"DCTIzjt(z%d,j%d,t%d)",z,j,t);
				double DLBCTIzjt=0,DUBCTIzjt=IloInfinity;
				IloRange DCTI(env,DLBCTIzjt,expr,DUBCTIzjt,DEquoition_equilibre);
				modelSlave1.add(DCTI);
				DCTIt.add(DCTI);
				expr.end();
			}
		
		}
	DCTIjt.add(DCTIt);
	}
DCTIzjt.add(DCTIjt);
}

//-----------------------------------------------------------------------------
//----------------------------- Capacite de Stockage --------------------------

for (z=0;z<zmax;z++){
	IloRangeArray Sum_It(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (j=0;j<jmax;j++){
				expr+=-Izjt[z][j][t];
			}
			char Capacite_Stockage[60];
			sprintf(Capacite_Stockage,"Sum_Izt(z%d,t%d)",z,t);
			double LBSum_Izt=-1*capmax,UBSum_Izt=IloInfinity;
			IloRange Sum_I(env,LBSum_Izt,expr,UBSum_Izt,Capacite_Stockage);
			modelSlave1.add(Sum_I);
			Sum_It.add(Sum_I);
			expr.end();
		}
	Sum_Izt.add(Sum_It);
}

//-----------------------------------------------------------------------------
//------------------------------ Contrainte de Fonctionement ------------------
//-----------------------------------------------------------------------------
//------------------------------- Chargement CT1 ------------------------------

for(i=0;i<imax;i++){
	IloRangeMatrix2x2 CT1Fonctionement_Czt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray CT1Fonctionement_Ct(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (j=0;j<jmax;j++){
				expr+=-Xizjt[i][z][j][t];
			}
			char ChargementCT1[60];
			sprintf(ChargementCT1,"CT1Fonctionement_Cizt(i%d,z%d,t%d)",i,z,t);
			double LBCaCT1=-m*CiztValue[i][z][t],UBCaCT1=IloInfinity;
			IloRange CT1Fonctionement_C(env,LBCaCT1,expr,UBCaCT1,ChargementCT1);
			modelSlave1.add(CT1Fonctionement_C);
			CT1Fonctionement_Ct.add(CT1Fonctionement_C);
			expr.end();
		}
		CT1Fonctionement_Czt.add(CT1Fonctionement_Ct);
	}
	CT1Fonctionement_Cizt.add(CT1Fonctionement_Czt);
}
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//------------------------------- Chargement CT2 ---------------------------------
//--------------------------------------------------------------------------------

for(i=0;i<imax;i++){
	IloRangeMatrix2x2 CT2Fonctionement_Czt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray CT2Fonctionement_Ct(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (j=0;j<jmax;j++){
				expr+=Xizjt[i][z][j][t];
			}
			char ChargementCT2[60];
			sprintf(ChargementCT2,"CT2Fonctionement_Cizt(i%d,z%d,t%d)",i,z,t);
			double LBCaCT2=CiztValue[i][z][t],UBCaCT2=IloInfinity;
			IloRange CT2Fonctionement_C(env,LBCaCT2,expr,UBCaCT2,ChargementCT2);
			modelSlave1.add(CT2Fonctionement_C);
			CT2Fonctionement_Ct.add(CT2Fonctionement_C);
			expr.end();
		}
		CT2Fonctionement_Czt.add(CT2Fonctionement_Ct);
	}
	CT2Fonctionement_Cizt.add(CT2Fonctionement_Czt);
}

//------------------------------- Dechargement CT1 ---------------------------------

for(k=0;k<kmax;k++){
	IloRangeMatrix2x2 CT1Fonctionement_Dzt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray CT1Fonctionement_Dt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (j=0;j<jmax;j++){
				expr+=-Yzkjt[z][k][j][t];
			}
			char DechargementCT1[60];
			sprintf(DechargementCT1,"CT1Fonctionement_Dkzt(k%d,z%d,t%d)",k,z,t);
			double LBDeCT1=-m*DkztValue[k][z][t],UBDeCT1=IloInfinity;
			IloRange CT1Fonctionement_D(env,LBDeCT1,expr,UBDeCT1,DechargementCT1);
			modelSlave1.add(CT1Fonctionement_D);
			CT1Fonctionement_Dt.add(CT1Fonctionement_D);
			expr.end();
		}
		CT1Fonctionement_Dzt.add(CT1Fonctionement_Dt);
	}
	CT1Fonctionement_Dkzt.add(CT1Fonctionement_Dzt);
}
//------------------------------- Dechargement CT2 ---------------------------------

for(k=0;k<kmax;k++){
	IloRangeMatrix2x2 CT2Fonctionement_Dzt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray CT2Fonctionement_Dt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			for (j=0;j<jmax;j++){
				expr+=Yzkjt[z][k][j][t];
			}
			char DechargementCT2[60];
			sprintf(DechargementCT2,"CT2Fonctionement_Dkzt(k%d,z%d,t%d)",k,z,t);
			double LBDeCT2=DkztValue[k][z][t],UBDeCT2=IloInfinity;
			IloRange CT2Fonctionement_D(env,LBDeCT2,expr,UBDeCT2,DechargementCT2);
			modelSlave1.add(CT2Fonctionement_D);
			CT2Fonctionement_Dt.add(CT2Fonctionement_D);
			expr.end();
		}
		CT2Fonctionement_Dzt.add(CT2Fonctionement_Dt);
	}
	CT2Fonctionement_Dkzt.add(CT2Fonctionement_Dzt);
}

//------------------------------------------------------------------------------
//--------------------------- Contrainte de melange ----------------------------
//-------------------------- Contrainte de melange CT1 -------------------------

for(z=0;z<zmax;z++){
	IloRangeMatrix2x2 CT1Meljt(env,0);
	for(j=0;j<jmax;j++){
		IloRangeArray CT1Melt(env,0);
		for(t=0;t<tmax;t++){
			IloExpr expr(env,0);
			expr+=-Izjt[z][j][t];
			char CT1Melange[60];
			sprintf(CT1Melange,"CT1Melzjt(z%d,j%d,t%d)",z,j,t);
			double LBCT1Melzjt=-m*FzjtValue[z][j][t],UBCT1Melzjt=IloInfinity;
			IloRange CT1Mel(env,LBCT1Melzjt,expr,UBCT1Melzjt,CT1Melange);
			modelSlave1.add(CT1Mel);
			CT1Melt.add(CT1Mel);
			expr.end();
		}
		CT1Meljt.add(CT1Melt);
	}
	CT1Melzjt.add(CT1Meljt);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------- Contrainte de melange CT2 -------------------------

for(z=0;z<zmax;z++){
	IloRangeMatrix2x2 CT2Meljt(env,0);
	for(j=0;j<jmax;j++){
		IloRangeArray CT2Melt(env,0);
		for(t=0;t<tmax;t++){
			IloExpr expr(env,0);			
			expr+=Izjt[z][j][t];
			char CT2Melange[60];
			sprintf(CT2Melange,"CT2Melzjt(z%d,j%d,t%d)",z,j,t);
			double LBCT2Melzjt=FzjtValue[z][j][t],UBCT2Melzjt=IloInfinity;
			IloRange CT2Mel(env,LBCT2Melzjt,expr,UBCT2Melzjt,CT2Melange);
			modelSlave1.add(CT2Mel);
			CT2Melt.add(CT2Mel);
			expr.end();
		}
		CT2Meljt.add(CT2Melt);
	}
	CT2Melzjt.add(CT2Meljt);
}

//------------------------------------------------------------------------------
//--------------------------- On paye le Setup au debut ------------------------
//----------------------------------- Chargement -----------------------------

for (i=0;i<imax;i++){
	IloRangeMatrix2x2 SC2_Czt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray SC2_Ct(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			expr+=-SCizt[i][z][t];
			char CTSC2_Cizt[60];
			sprintf(CTSC2_Cizt,"SC2_Cizt(i%d,z%d,t%d)",i,z,t);
			double LBSC2_Cizt2=-1*CiztValue[i][z][t],UBSC2_Cizt2=IloInfinity;
			IloRange SC2_C(env,LBSC2_Cizt2,expr,UBSC2_Cizt2,CTSC2_Cizt);
			modelSlave1.add(SC2_C);
			SC2_Ct.add(SC2_C);
			expr.end();
		}
		SC2_Czt.add(SC2_Ct);
	}
	SC2_Cizt.add(SC2_Czt);
}

//------------------------------------------------------------------------------
//--------------------------- On paye le Setup au debut ------------------------
//----------------------------------- Dechargement -----------------------------

for (k=0;k<kmax;k++){
	IloRangeMatrix2x2 SD2_Dzt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray SD2_Dt(env,0);
		for (t=0;t<tmax;t++){
			IloExpr expr(env,0);
			expr+=-SDkzt[k][z][t];
			char CTSD2_Dkzt[60];
			sprintf(CTSD2_Dkzt,"SD2_Dkzt(k%d,z%d,t%d)",k,z,t);
			double LBSD2_Dkzt2=-1*DkztValue[k][z][t],UBSD2_Dkzt2=IloInfinity;
			IloRange SD2_D(env,LBSD2_Dkzt2,expr,UBSD2_Dkzt2,CTSD2_Dkzt);
			modelSlave1.add(SD2_D);
			SD2_Dt.add(SD2_D);
			expr.end();
		}
		SD2_Dzt.add(SD2_Dt);
	}
	SD2_Dkzt.add(SD2_Dzt);
}


//--------------------------------------------------------------------------------
//--------------------------- On paye le Setup au debut ------------------------
//----------------------------------- Chargement ---------------------------------

for (i=0;i<imax;i++){
	IloRangeMatrix2x2 SC_Czt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray SC_Ct(env,0);
		for (t=0;t<tmax;t++){
			if(t==0){
				IloExpr expr(env,0);
				expr+=SCizt[i][z][t];
				char CTSC_Cizt[60];
				sprintf(CTSC_Cizt,"SC_Cizt(i%d,z%d,t%d)",i,z,t);
				double LBSC_Cizt1=CiztValue[i][z][t],UBSC_Cizt1=IloInfinity;
				IloRange SC_C(env,LBSC_Cizt1,expr,UBSC_Cizt1,CTSC_Cizt);
				modelSlave1.add(SC_C);
				SC_Ct.add(SC_C);
				expr.end();
			}
		
			else {
				IloExpr expr(env,0);
				expr+=SCizt[i][z][t];
				char CTSC_Cizt[60];
				sprintf(CTSC_Cizt,"SC_Cizt(i%d,z%d,t%d)",i,z,t);
				double LBSC_Cizt2=CiztValue[i][z][t]-CiztValue[i][z][t-1],UBSC_Cizt2=IloInfinity;
				IloRange SC_C(env,LBSC_Cizt2,expr,UBSC_Cizt2,CTSC_Cizt);
				modelSlave1.add(SC_C);
				SC_Ct.add(SC_C);
				expr.end();
			}
		}
		SC_Czt.add(SC_Ct);
	}
	SC_Cizt.add(SC_Czt);
}
//------------------------------------------------------------------------------
//--------------------------- On paye le Setup au debut ------------------------
//----------------------------------- Dechargement -----------------------------

for (k=0;k<kmax;k++){
	IloRangeMatrix2x2 SD_Dzt(env,0);
	for (z=0;z<zmax;z++){
		IloRangeArray SD_Dt(env,0);
		for (t=0;t<tmax;t++){
			if(t==0){
				IloExpr expr(env,0);
				expr+=SDkzt[k][z][t];
				char CTSD_Dkzt[60];
				sprintf(CTSD_Dkzt,"SD_Dkzt(k%d,z%d,t%d)",k,z,t);
				double LBSD_Dkzt1=DkztValue[k][z][t],UBSD_Dkzt1=IloInfinity;
				IloRange SD_D(env,LBSD_Dkzt1,expr,UBSD_Dkzt1,CTSD_Dkzt);
				modelSlave1.add(SD_D);
				SD_Dt.add(SD_D);
				expr.end();
			}

			else {
				IloExpr expr(env,0);
				expr+=SDkzt[k][z][t];
				char CTSD_Dkzt[60];
				sprintf(CTSD_Dkzt,"SD_Dkzt(k%d,z%d,t%d)",k,z,t);
				double LBSD_Dkzt2=DkztValue[k][z][t]-DkztValue[k][z][t-1],UBSD_Dkzt2=IloInfinity;
				IloRange SD_D(env,LBSD_Dkzt2,expr,UBSD_Dkzt2,CTSD_Dkzt);
				modelSlave1.add(SD_D);
				SD_Dt.add(SD_D);
				expr.end();
			}
		}
		SD_Dzt.add(SD_Dt);
	}
	SD_Dkzt.add(SD_Dzt);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------------Fonction Objectif de slave probleme--------------------------
//------------------------------------------------------------------------------
IloExpr expr_slave(env);

for (i=0;i<imax;i++){
	for (z=0;z<zmax;z++){
		for (t=0;t<tmax;t++){
			expr_slave+=SCizt[i][z][t];
		}
	}
}

for (k=0;k<kmax;k++){
	for (z=0;z<zmax;z++){
		for (t=0;t<tmax;t++){
			expr_slave+=SDkzt[k][z][t];
		}
	}
}

modelSlave1.add(IloMinimize(env, expr_slave));
expr_slave.end();

return 0;
}

int BendersIteration(IloModel modelMaster_ptr, IloModel modelSlave_ptr){

bool InfeasibleMaster=false; 
int loop=0;
IloCplex cplexSlave_ptr(env);
IloCplex cplexMaster_ptr(env);	

cplexSlave_ptr.extract(modelSlave_ptr);
cplexSlave_ptr.exportModel("modelSlave1.lp");

cplexMaster_ptr.extract(modelMaster_ptr);
cplexMaster_ptr.exportModel("modelMaster1.lp");

double DTransposeY=0, SlaveObjFunction=0;
double BestSlaveObj=100;

LowerBoundArray.clear();
UpperBoundArray.clear();
UpperBoundGlobalArray.clear();
dTy.clear();
zCurrent.clear();
cTx.clear();
BestSlaveObjSoFar.clear();
Time.clear();
BDFeasCuts=0;
BDOptCuts=0;


while (UpperBoundGlobal>LowerBound && loop<10000){	
	loop++;
	cout<<"-----------------"<<endl;
	cout<<"Iteration ="<<loop<<endl;
	cout<<"-----------------"<<endl;
	DTransposeY=0;
	cplexMaster_ptr.extract(modelMaster_ptr);
	//--------------SOLVE MASTER PROBLEM----------------	
	try {
		cplexMaster_ptr.exportModel("CurrentMaster.lp");
		cplexMaster_ptr.solve();
			
		if (!cplexMaster_ptr.solve ()){ 
			env.error()<<"Failed to optimize Master 1 LP PAME ALLI MIA THA TA KATAFERIS."<<endl;
			env.out()<<"----------------------------------------------------------------"<<endl;
			InfeasibleMaster=true;
			break;
		}

		env.out()<<"Solution status Master1 = "<<cplexMaster_ptr.getStatus()<<endl;
		env.out()<<"Solution value Master1= "<<cplexMaster_ptr.getObjValue()<<endl;

		//--------LOWER BOUND------------
		LowerBound=cplexMaster_ptr.getObjValue();
		for (t=0;t<tmax;t++){
			for (z=0;z<zmax;z++){
				for (i=0;i<imax;i++){
					CiztValue[i][z][t]=cplexMaster_ptr.getValue(Cizt[i][z][t]);
				}
				for (k=0;k<kmax;k++){
					DkztValue[k][z][t]=cplexMaster_ptr.getValue(Dkzt[k][z][t]);
				}
				for(j=0;j<jmax;j++){
					FzjtValue[z][j][t]=cplexMaster_ptr.getValue(Fzjt[z][j][t]);
				}
			}
		}

		for (n=0;n<1;n++){
			ThetaValue=cplexMaster_ptr.getValue(Zn[n]);
		}	

		for (t=0;t<tmax;t++){
			for (z=0;z<zmax;z++){
				for (i=0;i<imax;i++){
					DTransposeY+=CiztValue[i][z][t];
				}
				for (k=0;k<kmax;k++){
					DTransposeY+=DkztValue[k][z][t];
				}
			}
		}
		dTy.push_back(DTransposeY);
		zCurrent.push_back(ThetaValue);
		
		OptimalThetaValue=ThetaValue;
/*
		for (t=0;t<tmax;t++){
			for (z=0;z<zmax;z++){
				for (i=0;i<imax;i++){
					if (CiztValue[i][z][t]<fraction){
						CiztValue[i][z][t]=0;
					}
					else if (CiztValue[i][z][t]>=fraction){
						CiztValue[i][z][t]=1;
					}
				}
				for (k=0;k<kmax;k++){
					if (DkztValue[k][z][t]<fraction){
						DkztValue[k][z][t]=0;
					}
					else if (DkztValue[k][z][t]>=fraction){
						DkztValue[k][z][t]=1;
					}
				}
				for(j=0;j<jmax;j++){
					if (FzjtValue[z][j][t]<fraction){
						FzjtValue[z][j][t]=0;
					}
					else if (FzjtValue[z][j][t]>=fraction){
						FzjtValue[z][j][t]=1;
					}
				}
			}
		}
		*/
		
	}	
	catch ( IloException& e){
		cerr << "concert exception caught Master:"<<e<<endl;
	}
	catch (...){
		cerr<<"Unknown exception caught Master " <<endl;
	}


//---------------Update the LB of the constraints of SLAVE problem----------------
	for(i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				CT1Fonctionement_Cizt[i][z][t].setLB(-m*CiztValue[i][z][t]);
				CT2Fonctionement_Cizt[i][z][t].setLB(CiztValue[i][z][t]);
			}
		}
	}

	for(k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				CT1Fonctionement_Dkzt[k][z][t].setLB(-m*DkztValue[k][z][t]);
				CT2Fonctionement_Dkzt[k][z][t].setLB(DkztValue[k][z][t]);
			}
		}
	}
	
	for(z=0;z<zmax;z++){
		for (j=0;j<jmax;j++){
			for (t=0;t<tmax;t++){
				CT1Melzjt[z][j][t].setLB(-m*FzjtValue[z][j][t]);
				CT2Melzjt[z][j][t].setLB(FzjtValue[z][j][t]);
			}
		}
	}

	for(i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				SC2_Cizt[i][z][t].setLB(-1*CiztValue[i][z][t]);
			}
		}
	}

	for(k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				SD2_Dkzt[k][z][t].setLB(-1*DkztValue[k][z][t]);
			}
		}
	}

	for(i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(t==0){
					SC_Cizt[i][z][t].setLB(CiztValue[i][z][t]);
				}else{
					SC_Cizt[i][z][t].setLB(CiztValue[i][z][t]-CiztValue[i][z][t-1]);
				}
			}
		}
	}
	for(k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(t==0){
					SD_Dkzt[k][z][t].setLB(DkztValue[k][z][t]);
				}else{
					SD_Dkzt[k][z][t].setLB(DkztValue[k][z][t]-DkztValue[k][z][t-1]);
				}
			}
		}
	}
	cplexSlave_ptr.extract(modelSlave_ptr);
//-----------Solve Slave Problem--------------
	try {
		
		cplexSlave_ptr.setParam(IloCplex::PreInd,0); 
		cplexSlave_ptr.setParam(IloCplex::ScaInd,-1); 
		cplexSlave_ptr.setParam(IloCplex::RootAlg, IloCplex::Dual);
		cplexSlave_ptr.exportModel("CurrentSlave.lp");
		cplexSlave_ptr.solve();
		env.out()<<"Solution status of SLAVE problem = " <<cplexSlave_ptr.getStatus()<<endl;
		env.out()<<"Solution value of SLAVE problem = " <<cplexSlave_ptr.getObjValue()<<endl;

		if (!cplexSlave_ptr.solve ()){ //---------IF SLAVE1 IS INFEASIBLE-----
			env.error()<<"Failed to optimize Slave 1 LP PAME ALLI MIA THA TA KATAFERIS."<<endl;
			env.out()<<"----------------------------------------------------------------"<<endl;
			//------Upper Bound Global remains the same--------
			UpperBound=100;
			UpperBoundGlobal=UpperBoundGlobal;

			cTx.push_back(0);
			

			//----------------Get an extreme ray of the DUAL SLAVE problem-------------
			//cout<<"size of Array ="<<FeasvalsDualRangeSumXijt.getSize()<<endl;
			
			cplexSlave_ptr.dualFarkas(SumXijt[0][0],FeasvalsDualRangeSumXijt);
			//cout<<"size of Array ="<<FeasvalsDualRangeSumXijt.getSize()<<endl;
			/*
			for (l=0;l<FeasvalsDualRangeSumXijt.getSize();l++){
				if(FeasvalsDualRangeSumXijt[l]!=0){
					cout<<"FeasvalsDualRangeSumXijt["<<l<<"]="<<FeasvalsDualRangeSumXijt[l]<<endl;
				}
			}
				*/
			l=0;
			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S2valsDualSumXijt[i][j][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}
			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S22valsDualSumXijt[i][j][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S2valsDualSumYkjt[k][j][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}
			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S22valsDualSumYkjt[k][j][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}
				
			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S2valsDualCTIzjt[z][j][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){						
						S22valsDualCTIzjt[z][j][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (z=0;z<zmax;z++){
				for (t=0;t<tmax;t++){
					S2valsDualSum_Izt[z][t]=FeasvalsDualRangeSumXijt[l];
					l++;
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualCT1Fonctionement_Cizt[i][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualCT2Fonctionement_Cizt[i][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualCT1Fonctionement_Dkzt[k][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualCT2Fonctionement_Dkzt[k][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}


			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S2valsDualCT1Melzjt[z][j][t]=FeasvalsDualRangeSumXijt[l];		
						l++;
					}
				}
			}

			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						S2valsDualCT2Melzjt[z][j][t]=FeasvalsDualRangeSumXijt[l];		
						l++;
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualSC2_Cizt[i][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}
			
			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualSD2_Dkzt[k][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualSC_Cizt[i][z][t]=FeasvalsDualRangeSumXijt[l];		
						l++;
					}
				}
			}
			
			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						S2valsDualSD_Dkzt[k][z][t]=FeasvalsDualRangeSumXijt[l];
						l++;
					}
				}
			}

//---------CREATION OF BENDERS FEASIBILITY CUT--------------- 

			IloExpr expr101(env);

			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualSumXijt[i][j][t]*(E[i][t]*ae[i][j][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						expr101+=S22valsDualSumXijt[i][j][t]*(-E[i][t]*ae[i][j][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualSumYkjt[k][j][t]*(S[k][t]*as[k][j][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						expr101+=S22valsDualSumYkjt[k][j][t]*(-S[k][t]*as[k][j][t]);
					}
				}
			}

			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					expr101+=S2valsDualCTIzjt[z][j][0]*(initial[z][j]);
				}
			}

			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					expr101+=S22valsDualCTIzjt[z][j][0]*(-initial[z][j]);
				}
			}

			for (z=0;z<zmax;z++){
				for (t=0;t<tmax;t++){
					expr101+=S2valsDualSum_Izt[z][t]*(-1*capmax);
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualCT1Fonctionement_Cizt[i][z][t]*(-m*Cizt[i][z][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualCT2Fonctionement_Cizt[i][z][t]*(Cizt[i][z][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualCT1Fonctionement_Dkzt[k][z][t]*(-m*Dkzt[k][z][t]);
					}
				}
			}


			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualCT2Fonctionement_Dkzt[k][z][t]*(Dkzt[k][z][t]);
					}
				}
			}

				
			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualCT1Melzjt[z][j][t]*(-m*Fzjt[z][j][t]);
					}
				}
			}


			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualCT2Melzjt[z][j][t]*(Fzjt[z][j][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualSC2_Cizt[i][z][t]*(-1*Cizt[i][z][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						expr101+=S2valsDualSD2_Dkzt[k][z][t]*(-1*Dkzt[k][z][t]);
					}
				}
			}


			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						if(t==0){
							expr101+=S2valsDualSC_Cizt[i][z][t]*(Cizt[i][z][t]);
						}else{
							expr101+=S2valsDualSC_Cizt[i][z][t]*(Cizt[i][z][t]-Cizt[i][z][t-1]);
						}
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						if(t==0){
							expr101+=S2valsDualSD_Dkzt[k][z][t]*(Dkzt[k][z][t]);
						}else{
							expr101+=S2valsDualSD_Dkzt[k][z][t]*(Dkzt[k][z][t]-Dkzt[k][z][t-1]);
						}
					}
				}
			}


//--------------ADD BENDERS FEASIBILITY CUT TO MASTER----------------

			char MasterFeasibilityCut[90];
			sprintf(MasterFeasibilityCut,"FeasibilityCut(iter%d)",loop);
			double LB101=-IloInfinity,UB101=0;
			IloRange CTMaster(env,LB101,expr101,UB101,MasterFeasibilityCut);
			modelMaster_ptr.add(CTMaster);
			expr101.end();
			BDFeasCuts++;
						
		}//Fin de IF QUI A TROUVE QUE SLAVE 1 EST INFEASIBLE
			
		else { //------------- IF SLAVE PROBLEM IS FEASIBLE------------
			
			//cplexMaster_ptr.exportModel("TELOS OPTIMISATION TIMER.lp");		

			env.error()<<"Found a FEASIBLE solution of Slave LP"<<endl;
			env.out()<<"----------------------------------------------------------------"<<endl;
			SlaveObjFunction=cplexSlave_ptr.getObjValue();
			UpperBound= DTransposeY + SlaveObjFunction;
		
			for (t=0;t<tmax;t++){
				for (z=0;z<zmax;z++){
					for (i=0;i<imax;i++){
						SCiztValue[i][z][t]=cplexSlave_ptr.getValue(SCizt[i][z][t]);
					}
					for (k=0;k<kmax;k++){
						SDkztValue[k][z][t]=cplexSlave_ptr.getValue(SDkzt[k][z][t]);
					}
				}
			}
			for (t=0;t<tmax;t++){
				for (z=0;z<zmax;z++){
					for (j=0;j<jmax;j++){
						for (i=0;i<imax;i++){
							XizjtValue[i][z][j][t]=cplexSlave_ptr.getValue(Xizjt[i][z][j][t]);
						}
						for (k=0;k<kmax;k++){
							YzkjtValue[z][k][j][t]=cplexSlave_ptr.getValue(Yzkjt[z][k][j][t]);
						}
						IzjtValue[z][j][t]=cplexSlave_ptr.getValue(Izjt[z][j][t]);
					}
				}
			}

			if( UpperBound >= UpperBoundGlobal){//-----We found a worse feasible solution---
				UpperBoundGlobal=UpperBoundGlobal;
			}else{//-----------We found a better feasible solution-------
				UpperBoundGlobal= UpperBound;//Update Upper Bound
				OptimalOriginalObjFunction=UpperBoundGlobal;
				OptimalMasterObjFunction=DTransposeY;
				OptimalSlaveObjFunction=SlaveObjFunction;
				for (t=0;t<tmax;t++){
					for (z=0;z<zmax;z++){
						for (i=0;i<imax;i++){
							OptimalCiztValue[i][z][t]=CiztValue[i][z][t];
						}
						for (k=0;k<kmax;k++){
							OptimalDkztValue[k][z][t]=DkztValue[k][z][t];
						}
						for(j=0;j<jmax;j++){
							OptimalFzjtValue[z][j][t]=FzjtValue[z][j][t];
						}
					}
				}

				for (t=0;t<tmax;t++){
					for (z=0;z<zmax;z++){
						for (i=0;i<imax;i++){
							OptimalSCiztValue[i][z][t]=SCiztValue[i][z][t];
						}
						for (k=0;k<kmax;k++){
							OptimalSDkztValue[k][z][t]=SDkztValue[k][z][t];
						}
					}
				}
				for (t=0;t<tmax;t++){
					for (z=0;z<zmax;z++){
						for (j=0;j<jmax;j++){
							for (i=0;i<imax;i++){
								OptimalXizjtValue[i][z][j][t]=XizjtValue[i][z][j][t];
							}
							for (k=0;k<kmax;k++){
								OptimalYzkjtValue[z][k][j][t]=YzkjtValue[z][k][j][t];
							}
							OptimalIzjtValue[z][j][t]=IzjtValue[z][j][t];
						}
					}
				}
				OptimalThetaValue=ThetaValue;
				
			}//end of else (better feasible solution found)

			cTx.push_back(SlaveObjFunction);
			


			
			//---------------------Get an extreme point of DUAL SLAVE problem--------------------

			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						valsDualRangeSumXijt=cplexSlave_ptr.getDual(SumXijt[i][j][t]);
						S2valsDualSumXijt[i][j][t]=valsDualRangeSumXijt;
						
						valsDualRangeSumXijt=cplexSlave_ptr.getDual(DSumXijt[i][j][t]);
						S22valsDualSumXijt[i][j][t]=valsDualRangeSumXijt;
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						valsDualRangeSumYkjt=cplexSlave_ptr.getDual(SumYkjt[k][j][t]);
						S2valsDualSumYkjt[k][j][t]=valsDualRangeSumYkjt;

						valsDualRangeSumYkjt=cplexSlave_ptr.getDual(DSumYkjt[k][j][t]);
						S22valsDualSumYkjt[k][j][t]=valsDualRangeSumYkjt;
					}
				}
			}
				
			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						valsDualRangeCTIzjt=cplexSlave_ptr.getDual(CTIzjt[z][j][t]);
						S2valsDualCTIzjt[z][j][t]=valsDualRangeCTIzjt;
						
						valsDualRangeCTIzjt=cplexSlave_ptr.getDual(DCTIzjt[z][j][t]);
						S22valsDualCTIzjt[z][j][t]=valsDualRangeCTIzjt;
					}
				}
			}

			for (z=0;z<zmax;z++){
				for (t=0;t<tmax;t++){
					valsDualRangeSum_Izt=cplexSlave_ptr.getDual(Sum_Izt[z][t]);
					S2valsDualSum_Izt[z][t]=valsDualRangeSum_Izt;

				}
			}


			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeCT1Fonctionement_Cizt=cplexSlave_ptr.getDual(CT1Fonctionement_Cizt[i][z][t]);
						S2valsDualCT1Fonctionement_Cizt[i][z][t]=valsDualRangeCT1Fonctionement_Cizt;
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeCT2Fonctionement_Cizt=cplexSlave_ptr.getDual(CT2Fonctionement_Cizt[i][z][t]);
						S2valsDualCT2Fonctionement_Cizt[i][z][t]=valsDualRangeCT2Fonctionement_Cizt;
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeCT1Fonctionement_Dkzt=cplexSlave_ptr.getDual(CT1Fonctionement_Dkzt[k][z][t]);
						S2valsDualCT1Fonctionement_Dkzt[k][z][t]=valsDualRangeCT1Fonctionement_Dkzt;
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeCT2Fonctionement_Dkzt=cplexSlave_ptr.getDual(CT2Fonctionement_Dkzt[k][z][t]);
						S2valsDualCT2Fonctionement_Dkzt[k][z][t]=valsDualRangeCT2Fonctionement_Dkzt;
					}
				}
			}


			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						valsDualRangeCT1Melzjt=cplexSlave_ptr.getDual(CT1Melzjt[z][j][t]);
						S2valsDualCT1Melzjt[z][j][t]=valsDualRangeCT1Melzjt;						
					}
				}
			}

			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						valsDualRangeCT2Melzjt=cplexSlave_ptr.getDual(CT2Melzjt[z][j][t]);
						S2valsDualCT2Melzjt[z][j][t]=valsDualRangeCT2Melzjt;						
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeSC2_Cizt=cplexSlave_ptr.getDual(SC2_Cizt[i][z][t]);
						S2valsDualSC2_Cizt[i][z][t]=valsDualRangeSC2_Cizt;							
					}
				}
			}
			
			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeSD2_Dkzt=cplexSlave_ptr.getDual(SD2_Dkzt[k][z][t]);
						S2valsDualSD2_Dkzt[k][z][t]=valsDualRangeSD2_Dkzt;						
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeSC_Cizt=cplexSlave_ptr.getDual(SC_Cizt[i][z][t]);
						S2valsDualSC_Cizt[i][z][t]=valsDualRangeSC_Cizt;						
					}
				}
			}
			
			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						valsDualRangeSD_Dkzt=cplexSlave_ptr.getDual(SD_Dkzt[k][z][t]);
						S2valsDualSD_Dkzt[k][z][t]=valsDualRangeSD_Dkzt;						
					}
				}
			}

//---------CREATION OF BENDERS OPTIMALITY CUT--------------- 

			IloExpr exprOptimality(env);

			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualSumXijt[i][j][t]*(E[i][t]*ae[i][j][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S22valsDualSumXijt[i][j][t]*(-E[i][t]*ae[i][j][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualSumYkjt[k][j][t]*(S[k][t]*as[k][j][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S22valsDualSumYkjt[k][j][t]*(-S[k][t]*as[k][j][t]);
					}
				}
			}

			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					exprOptimality+=S2valsDualCTIzjt[z][j][0]*(initial[z][j]);
				}
			}

			for (z=0;z<zmax;z++){
				for (j=0;j<jmax;j++){
					exprOptimality+=S22valsDualCTIzjt[z][j][0]*(-initial[z][j]);
				}
			}

			for (z=0;z<zmax;z++){
				for (t=0;t<tmax;t++){
					exprOptimality+=S2valsDualSum_Izt[z][t]*(-1*capmax);
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualCT1Fonctionement_Cizt[i][z][t]*(-m*Cizt[i][z][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualCT2Fonctionement_Cizt[i][z][t]*(Cizt[i][z][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualCT1Fonctionement_Dkzt[k][z][t]*(-m*Dkzt[k][z][t]);
					}
				}
			}


			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualCT2Fonctionement_Dkzt[k][z][t]*(Dkzt[k][z][t]);
					}
				}
			}

			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualCT1Melzjt[z][j][t]*(-m*Fzjt[z][j][t]);
					}
				}
			}

			for (z=0;z<zmax;z++){
				for(j=0;j<jmax;j++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualCT2Melzjt[z][j][t]*(Fzjt[z][j][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualSC2_Cizt[i][z][t]*(-1*Cizt[i][z][t]);
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						exprOptimality+=S2valsDualSD2_Dkzt[k][z][t]*(-1*Dkzt[k][z][t]);
					}
				}
			}

			for (i=0;i<imax;i++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						if(t==0){
							exprOptimality+=S2valsDualSC_Cizt[i][z][t]*(Cizt[i][z][t]);
						}else{
							exprOptimality+=S2valsDualSC_Cizt[i][z][t]*(Cizt[i][z][t]-Cizt[i][z][t-1]);
						}
					}
				}
			}

			for (k=0;k<kmax;k++){
				for (z=0;z<zmax;z++){
					for (t=0;t<tmax;t++){
						if(t==0){
							exprOptimality+=S2valsDualSD_Dkzt[k][z][t]*(Dkzt[k][z][t]);
						}else{
							exprOptimality+=S2valsDualSD_Dkzt[k][z][t]*(Dkzt[k][z][t]-Dkzt[k][z][t-1]);
						}
					}
				}
			}

			for (n=0;n<1;n++){
				exprOptimality-=Zn[n];
			}

//--------------ADD BENDERS OPTIMALITY CUT TO MASTER----------------

			char MasterOptimalityCut[90];
			sprintf(MasterOptimalityCut,"OptimalityCut(iter%d)",loop);
			double LB101=-IloInfinity,UB101=0;
			IloRange CTMaster(env,LB101,exprOptimality,UB101,MasterOptimalityCut);
			modelMaster_ptr.add(CTMaster);
			exprOptimality.end();
			BDOptCuts++;
		}//end of else

		BestSlaveObjSoFar.push_back(OptimalSlaveObjFunction);


	}//end of try(try of primal slave 1)
	
	catch ( IloException& e){
		cerr << "concert exception caught:"<<e<<endl;
	}
	catch (...){
		cerr<<"Unknown exception caught" <<endl;
	}

	LowerBoundArray.push_back(LowerBound);
	UpperBoundArray.push_back(UpperBound);
	UpperBoundGlobalArray.push_back(UpperBoundGlobal);
	Time.push_back((long double)(clock()-start)/CLOCKS_PER_SEC);

	
	//------------ Before we got to loop t+1 and solve the new Master Problem we update the LB and UB---- 
/*
	SlaveObjFunction=0;
	DTransposeY=0;
	

	//--------UPPER BOUND------------
	for (i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){	
				DTransposeY+=CiztValue[i][z][t];
			}
		}
	}
	for (k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				DTransposeY+=DkztValue[k][z][t];
			}
		}
	}

	if (!cplexSlave_ptr.solve ()){ //---------IF SLAVE1 IS INFEASIBLE-----
		UpperBound=100;
		UpperBoundGlobal=UpperBoundGlobal;
	}else{
		SlaveObjFunction=cplexSlave_ptr.getObjValue();
		
		UpperBound= DTransposeY + SlaveObjFunction;
		
		if( SlaveObjFunction < BestSlaveObj){
			BestSlaveObj=SlaveObjFunction;
		}
		if( DTransposeY + SlaveObjFunction > UpperBoundGlobal){
			UpperBoundGlobal=UpperBoundGlobal;
		}else{
			UpperBoundGlobal= DTransposeY + SlaveObjFunction;
		}
	}
	
	for (i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(CiztValue[i][z][t]!=0){
					cout<<"CiztValue"<<"["<<i<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<CiztValue[i][z][t]<<endl;
				}
			}
		}
	}
	cout<<"----------------------------------"<<endl;

	for (k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(DkztValue[k][z][t]!=0){
					cout<<"DkztValue"<<"["<<k<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<DkztValue[k][z][t]<<endl;
				}
			}
		}
	}
	cout<<"----------------------------------"<<endl;

	for (z=0;z<zmax;z++){
		for (j=0;j<jmax;j++){
			for (t=0;t<tmax;t++){
				if(FzjtValue[z][j][t]!=0){
					cout<<"FzjtValue"<<"["<<z<<"]"<<"["<<j<<"]"<<"["<<t<<"]"<<"="<<FzjtValue[z][j][t]<<endl;
				}
			}
		}
	}
	cout<<"----------------------------------"<<endl;

	for (i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(SCiztValue[i][z][t]!=0){
					cout<<"SCiztValue"<<"["<<i<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<SCiztValue[i][z][t]<<endl;
				}
			}
		}
	}
	cout<<"----------------------------------"<<endl;

	for (k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(SDkztValue[k][z][t]!=0){
					cout<<"SDkztValue"<<"["<<k<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<SDkztValue[k][z][t]<<endl;
				}
			}
		}
	}
	cout<<"----------------------------------"<<endl;
*/
	if(ThetaValue!=0){
		cout<<"OptimalThetaValue="<<OptimalThetaValue<<endl;
	}
	if(DTransposeY!=0){
		cout<<"DTransposeY="<<DTransposeY<<endl;
	}
	if(SlaveObjFunction!=0){
		cout<<"SlaveObjFunction="<<SlaveObjFunction<<endl;
	}
	if(OptimalSlaveObjFunction!=0){
		cout<<"OptimalSlaveObjFunction="<<OptimalSlaveObjFunction<<endl;
	}
	cout<<"LowerBound="<<LowerBound<<endl;

	cout<<"UpperBoundGlobal="<<UpperBoundGlobal<<endl;

	cout<<"UpperBound="<<UpperBound<<endl;
/*
	std::ofstream fsLowerBound; 
	fsLowerBound.open("C:\\Apotelesmata_Kwdika\\LowerBound.txt",std::ios::app );
	fsLowerBound<<LowerBound<< std::endl;
	fsLowerBound.close();

	std::ofstream fsUpperBound; 
	fsUpperBound.open("C:\\Apotelesmata_Kwdika\\UpperBound.txt",std::ios::app );
	fsUpperBound<<UpperBound<< std::endl;
	fsUpperBound.close();

	std::ofstream fsUpperBoundGlobalArray; 
	fsUpperBoundGlobalArray.open("C:\\Apotelesmata_Kwdika\\UpperBoundGlobalArray.txt",std::ios::app );
	fsUpperBoundGlobalArray<<UpperBoundGlobal<< std::endl;
	fsUpperBoundGlobalArray.close();
	*/
	cout<<"-----------------"<<endl;
	cout<<"------fIN--------"<<endl;

}//end of loop

return 0;
}//end of BendersIteration

int PrintOptimalSolution(){
	std::ostringstream os;
	os << "C:\\Results_CrudeOil\\WithoutCut - OptimalSolution.txt";
	std::string FileName = os.str();
	
	std::ofstream fsOptimalSolution; 
	fsOptimalSolution.open(FileName.c_str(),std::ios::app );
	fsOptimalSolution<<"TotalSolutionTime= "<<duration<<" seconds "<< std::endl;
	if((UpperBoundGlobal-LowerBound)/UpperBoundGlobal>0){
		fsOptimalSolution<<"OptimalityGap= "<<(UpperBoundGlobal-LowerBound)/UpperBoundGlobal<< std::endl;
	}else{
		fsOptimalSolution<<"OptimalityGap= 0"<< std::endl;
	}
	fsOptimalSolution<<"OptimalObjFunction= "<<OptimalOriginalObjFunction<< std::endl;
	fsOptimalSolution<<"OptimalMasterObjFunction= "<<OptimalMasterObjFunction<< std::endl;
	fsOptimalSolution<<"OptimalSlaveObjFunction= "<<OptimalSlaveObjFunction<< std::endl;
	fsOptimalSolution<<"----------------------------------"<< std::endl;
	if(OptimalThetaValue>0.01){
		fsOptimalSolution<<"OptimalThetaValue= "<<OptimalThetaValue<<std::endl;
	}	
	fsOptimalSolution<<"----------------------------------"<< std::endl;
	fsOptimalSolution<<"TotalNumberOfFeasibilityCuts= "<<BDFeasCuts<<std::endl;
	fsOptimalSolution<<"TotalNumberOfOptimalityCuts= "<<BDOptCuts<<std::endl;
	fsOptimalSolution<<"----------------------------------"<< std::endl;
	for (i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(OptimalCiztValue[i][z][t]>=0.01){
					fsOptimalSolution<<"OptimalCiztValue"<<"["<<i<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<OptimalCiztValue[i][z][t]<<std::endl;
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;

	for (k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(OptimalDkztValue[k][z][t]>=0.01){
					fsOptimalSolution<<"OptimalDkztValue"<<"["<<k<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<OptimalDkztValue[k][z][t]<<std::endl;
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;

	for (z=0;z<zmax;z++){
		for (j=0;j<jmax;j++){
			for (t=0;t<tmax;t++){
				if(OptimalFzjtValue[z][j][t]>=0.01){
					fsOptimalSolution<<"OptimalFzjtValue"<<"["<<z<<"]"<<"["<<j<<"]"<<"["<<t<<"]"<<"="<<OptimalFzjtValue[z][j][t]<<std::endl;
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;

	for (i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(OptimalSCiztValue[i][z][t]>=0.01){
					fsOptimalSolution<<"OptimalSCiztValue"<<"["<<i<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<OptimalSCiztValue[i][z][t]<<std::endl;
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;

	for (k=0;k<kmax;k++){
		for (z=0;z<zmax;z++){
			for (t=0;t<tmax;t++){
				if(OptimalSDkztValue[k][z][t]>=0.01){
					fsOptimalSolution<<"OptimalSDkztValue"<<"["<<k<<"]"<<"["<<z<<"]"<<"["<<t<<"]"<<"="<<OptimalSDkztValue[k][z][t]<<std::endl;
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;
	
	for (i=0;i<imax;i++){
		for (z=0;z<zmax;z++){
			for (j=0;j<jmax;j++){
				for (t=0;t<tmax;t++){
					if(OptimalXizjtValue[i][z][j][t]>=0.01){
						fsOptimalSolution<<"OptimalXizjtValue"<<"["<<i<<"]"<<"["<<z<<"]"<<"["<<j<<"]"<<"["<<t<<"]"<<"="<<OptimalXizjtValue[i][z][j][t]<<std::endl;
					}
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;

	for (z=0;z<zmax;z++){
		for (k=0;k<kmax;k++){
			for (j=0;j<jmax;j++){
				for (t=0;t<tmax;t++){
					if(OptimalYzkjtValue[z][k][j][t]>=0.01){
						fsOptimalSolution<<"OptimalYzkjtValue"<<"["<<z<<"]"<<"["<<k<<"]"<<"["<<j<<"]"<<"["<<t<<"]"<<"="<<OptimalYzkjtValue[z][k][j][t]<<std::endl;
					}
				}
			}
		}
	}
	fsOptimalSolution<<"----------------------------------"<<std::endl;

	
	for (z=0;z<zmax;z++){
		for (j=0;j<jmax;j++){
			for (t=0;t<tmax;t++){
				if(OptimalIzjtValue[z][j][t]>=0.01){
					fsOptimalSolution<<"OptimalIzjtValue"<<"["<<z<<"]"<<"["<<j<<"]"<<"["<<t<<"]"<<"="<<OptimalIzjtValue[z][j][t]<<std::endl;
				}
			}
		}
	}

	fsOptimalSolution.close();


	std::ostringstream LowerBound;
	LowerBound << "C:\\Results_CrudeOil\\WithoutCut - LowerBound.txt";
	std::string FileNameLB = LowerBound.str();
	std::ofstream fsLowerBound; 
	fsLowerBound.open(FileNameLB.c_str(),std::ios::app );
	for (i=0;i<LowerBoundArray.size();i++){
		fsLowerBound<<LowerBoundArray.at(i)<< std::endl;
	}
	fsLowerBound.close();

	std::ostringstream UpperBound;
	UpperBound << "C:\\Results_CrudeOil\\WithoutCut - UpperBound.txt";
	std::string FileNameUB = UpperBound.str();
	std::ofstream fsUpperBound; 
	fsUpperBound.open(FileNameUB.c_str(),std::ios::app );
	for (i=0;i<UpperBoundArray.size();i++){
		fsUpperBound<<UpperBoundArray.at(i)<< std::endl;
	}
	fsUpperBound.close();

	std::ostringstream UpperBoundGlobal;
	UpperBoundGlobal << "C:\\Results_CrudeOil\\WithoutCut - UpperBoundGlobal.txt";
	std::string FileNameUBG = UpperBoundGlobal.str();
	std::ofstream fsUpperBoundGlobal; 
	fsUpperBoundGlobal.open(FileNameUBG.c_str(),std::ios::app );
	for (i=0;i<UpperBoundGlobalArray.size();i++){
		fsUpperBoundGlobal<<UpperBoundGlobalArray.at(i)<< std::endl;
	}
	fsUpperBoundGlobal.close();


	std::ostringstream dTransY;
	dTransY << "C:\\Results_CrudeOil\\WithoutCut - DTrasnposeY.txt";
	std::string FileNameDTY = dTransY.str();
	std::ofstream fsdTransY; 
	fsdTransY.open(FileNameDTY.c_str(),std::ios::app );
	for (i=0;i<dTy.size();i++){
		fsdTransY<<dTy.at(i)<< std::endl;
	}
	fsdTransY.close();

	std::ostringstream cTransX;
	cTransX << "C:\\Results_CrudeOil\\WithoutCut - CTrasnposeX.txt";
	std::string FileNameCTX = cTransX.str();
	std::ofstream fscTransX; 
	fscTransX.open(FileNameCTX.c_str(),std::ios::app );
	for (i=0;i<cTx.size();i++){
		fscTransX<<cTx.at(i)<< std::endl;
	}
	fscTransX.close();

	std::ostringstream CurrentTheta;
	CurrentTheta << "C:\\Results_CrudeOil\\WithoutCut - CurrentTheta.txt";
	std::string FileNameCurrentTheta = CurrentTheta.str();
	std::ofstream fsCurrentTheta; 
	fsCurrentTheta.open(FileNameCurrentTheta.c_str(),std::ios::app );
	for (i=0;i<zCurrent.size();i++){
		fsCurrentTheta<<zCurrent.at(i)<< std::endl;
	}
	fsCurrentTheta.close();

	std::ostringstream BestSlaveObj;
	BestSlaveObj << "C:\\Results_CrudeOil\\WithoutCut - BestSlaveObjSoFar.txt";
	std::string FileNameBSO = BestSlaveObj.str();
	std::ofstream fsBestSlaveObj; 
	fsBestSlaveObj.open(FileNameBSO.c_str(),std::ios::app );
	for (i=0;i<BestSlaveObjSoFar.size();i++){
		fsBestSlaveObj<<BestSlaveObjSoFar.at(i)<< std::endl;
	}
	fsBestSlaveObj.close();


	std::ostringstream TimePath;
	TimePath << "C:\\Results_CrudeOil\\WithoutCut - Time.txt";
	std::string FileNameTime = TimePath.str();
	std::ofstream fsTime; 
	fsTime.open(FileNameTime.c_str(),std::ios::app );
	for (i=0;i<Time.size();i++){
		fsTime<<Time.at(i)<< std::endl;
	}
	fsTime.close();

return 0;
}

int main (int argc, char **argv)
{
int  stop, status;


start = clock();

status=load_data();
if (status != 0){
	Found_Error("load_data");
	return -1;
}

status=do_master();
if (status != 0){
	Found_Error("do_master");
	return -1;
}

status=do_slave();
if (status != 0){
	Found_Error("do_slave");
	return -1;
}

status=BendersIteration(modelMaster, modelSlave1);
if (status != 0){
	Found_Error("BendersIteration");
	return -1;
}
stop = clock();
duration = (long double) (stop-start)/CLOCKS_PER_SEC;

status=PrintOptimalSolution();
if (status != 0){
	Found_Error("PrintOptimalSolution");
	return -1;
}

env.end();


printf("Code terminated successfully \n");
printf("Execution time = %Lf seconds\n", duration);

return 0;

} //End main




