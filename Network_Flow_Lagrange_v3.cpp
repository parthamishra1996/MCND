/*
Version 1: Multiplier values were input manually 
Version 2: Cutting plane used but problem became unbounded
Version3: Solving simpler MCND using explicit upper bound
*/
//STANDARD C++ LIBRARIES
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <queue>		//For implementing branching
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <typeinfo>
#include <limits.h>		//For using INT_MAX

//CONCERT TECHNOLOGY LIBRARIES
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>
#include <ilconcert/ilocsvreader.h>

ILOSTLBEGIN


typedef IloArray<IloBoolArray> Bool2DMatrix;		//2D array of Bool
typedef IloArray<Bool2DMatrix> Bool3DMatrix;		//3D array of Bool
typedef IloArray<IloIntArray> Int2DMatrix;			//2D array of Int
typedef IloArray<Int2DMatrix> Int3DMatrix;			//3D array of Int
typedef IloArray<IloNumArray> Num2DMatrix;			//2D array of Num
typedef IloArray<Num2DMatrix> Num3DMatrix;			//3D array of Num
typedef IloArray<Num3DMatrix> Num4DMatrix;			//4D array of Num
typedef IloArray<IloBoolVarArray> BoolVar2DMatrix;	// 2D array of Bool Var
typedef IloArray<BoolVar2DMatrix> BoolVar3DMatrix;	//3D array of Var
typedef IloArray<IloIntVarArray> IntVar2DMatrix;	// 2D array of Int Var
typedef IloArray<IntVar2DMatrix> IntVar3DMatrix;	//3D array of Int Var
typedef IloArray<IloNumVarArray> NumVar2DMatrix;	// 2D array of Num Var
typedef IloArray<NumVar2DMatrix> NumVar3DMatrix;	//3D array of Num Var
typedef IloArray<NumVar3DMatrix> NumVar4DMatrix;	//4D array of Num Var
typedef IloArray<IloRangeArray> Range2DMatrix;		//2D Arrays of Ranges

using namespace std;

IloEnv env;

int getX(int n, int WS_Y)
{
	return n/WS_Y;
}

int getY(int n, int WS_Y)
{
	return n%WS_Y;
}

int main(int argc, char **argv)
{
	clock_t t1,t2;
	t1=clock();
	/*FILE READING*/
	IloNumArray r_val(env);	//Rh values for a priori cuts

	const char* data_filename = "Herman_congestion_test.txt";
	const char* filename = "C:/Users/Partha Mishra/Documents/Visual Studio 2012/Projects/MontrealProject/x64/Debug/r_val_0.01.dat";
	if (argc > 1)
	{
		data_filename = argv[1];
	}
	fstream datafile;
	fstream file;
	file.open(filename,ios::in);
	datafile.open(data_filename, ios::in);
	
	if (!datafile)
	{
		cerr << "ERROR: could not open file " << data_filename << " for reading" << endl;
		cerr << "usage:   " << argv[0] << " <datafile>" << endl;
		return -1;
	}
	//std::cout << "openned" << endl;
	//system("pause");

	/*PARAMETER READING*/

	int K;			//Total number of connections
	int L;			//Number of capacity levels
	int WS_X;		//Workspace size
	int WS_Y;		//Workspace size
	//int F;			//Fixed cost of each arc
	//int B;			//Capacity of ach arc
	int C;			//Variable cost for each arc
	float penalty;	//Unit penalty at each arc for congestion

	datafile>>WS_X>>WS_Y>>K>>L>>C>>penalty;

	IloNumArray fix(env, L);	//Fixed cost of all arcs at different capacity levels
	IloNumArray cap(env, L);	//Capacity of all arcs at different levels
	IloNumArray f(env, K);	//Flow/demand of each connection
	IloIntArray O(env, K);	//Origin of each connection
	IloIntArray D(env, K);	//Destination of each connection

	Num3DMatrix Fx(env,WS_X);
	Num3DMatrix Fy(env,WS_X-1);
	Num3DMatrix Bx(env,WS_X);
	Num3DMatrix By(env,WS_X-1);

	datafile>>fix>>cap>>f>>O>>D;

	for(int i = 0; i < WS_X; i++)
	{
		Fx[i] = Num2DMatrix(env,WS_Y-1);
		Bx[i] = Num2DMatrix(env,WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Fx[i][j] = IloNumArray(env,L);
			Bx[i][j] = IloNumArray(env,L);
			for(int l = 0; l < L; l++)
			{
				Fx[i][j][l] = fix[l];
				Bx[i][j][l] = cap[l];
			}
		}
	}

	for(int i = 0; i < WS_X-1; i++)
	{
		Fy[i] = Num2DMatrix(env,WS_Y);
		By[i] = Num2DMatrix(env,WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			Fy[i][j] = IloNumArray(env,L);
			By[i][j] = IloNumArray(env,L);
			for(int l = 0; l < L; l++)
			{
				Fy[i][j][l] = fix[l];
				By[i][j][l] = cap[l];
			}
		}
	}

	std::cout<<K<<":"<<fix<<":"<<cap<<":"<<penalty<<endl;
	std::cout<<f<<endl;
	std::cout<<O<<endl;
	std::cout<<D<<endl;
	
	//system("pause");
	
	Num3DMatrix alpha_x(env,WS_X);
	Num3DMatrix alpha_y(env,WS_X-1);

	for(int i = 0; i < WS_X; i++)
	{
		alpha_x[i] = Num2DMatrix(env,WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			alpha_x[i][j] = IloNumArray(env,K);
			for(int k = 0; k < K; k++)
			{
				alpha_x[i][j][k] = 0;
			}
		}
	}
	//cout<<"Fine!"<<endl;
	for(int i = 0; i < WS_X-1; i++)
	{
		alpha_y[i] = Num2DMatrix(env,WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			alpha_y[i][j] = IloNumArray(env,K);
			for(int k = 0; k < K; k++)
			{
				alpha_y[i][j][k] = 0;
			}
		}
	}
	
	Num3DMatrix P_val(env,WS_X);
	for(int i = 0; i < WS_X; i++)
	{
		//P[i] = NumVar2DMatrix(env, WS_Y-1);
		P_val[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			//P[i][j] = IloNumVarArray(env, K);
			P_val[i][j] = IloNumArray(env, K);
		}
	}

	Num3DMatrix Q_val(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		Q_val[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Q_val[i][j] = IloNumArray(env, K);
		}
	}

	Num3DMatrix U_val(env,WS_X-1);
	for(int i = 0; i < WS_X-1; i++)
	{
		U_val[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			U_val[i][j] = IloNumArray(env, K);
		}
	}
	//cout<<"Fine! U "<<endl;
	Num3DMatrix V_val(env,WS_X-1);
	for(int i = 0; i < WS_X-1; i++)
	{
		V_val[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			V_val[i][j] = IloNumArray(env, K);
		}
	}
	
	Num3DMatrix Wh_final(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		Wh_final[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Wh_final[i][j] = IloNumArray(env, K);
		}
	}

	Num3DMatrix Wv_final(env,WS_X-1);

	for(int i = 0; i < WS_X-1; i++)
	{
		Wv_final[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			Wv_final[i][j] = IloNumArray(env, K);
		}
	}
	
	Int3DMatrix P1_val(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		P1_val[i] = Int2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			P1_val[i][j] = IloIntArray(env, L);
		}
	}
	
	Int3DMatrix U1_val(env,WS_X-1);

	for(int i = 0; i < WS_X-1; i++)
	{
		U1_val[i] = Int2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			U1_val[i][j] = IloIntArray(env, L);
		}
	}

	IloNumVar theta(env,0, 20000);//IloInfinity);
	NumVar3DMatrix omega_x(env,WS_X);
	NumVar3DMatrix omega_y(env,WS_X-1);
	IloModel model_LR(env);
	IloCplex cplex_LR(env);
	model_LR.add(IloMaximize(env,theta));
	model_LR.add(theta <= 20000 );	
	for(int i = 0; i < WS_X; i++)
	{
		omega_x[i] = NumVar2DMatrix(env,WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			omega_x[i][j] = IloNumVarArray(env,K);
			for(int k = 0; k < K; k++)
			{
				omega_x[i][j][k] = IloNumVar(env,0,20000);
				model_LR.add(omega_x[i][j][k] <= 20000);
			}
		}
	}
	for(int i = 0; i < WS_X-1; i++)
	{
		omega_y[i] = NumVar2DMatrix(env,WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			omega_y[i][j] = IloNumVarArray(env,K);
			for(int k = 0; k < K; k++)
			{
				omega_y[i][j][k] = IloNumVar(env,0,20000);
				model_LR.add(omega_y[i][j][k] <= 20000);
			}
		}
	}

	IloNum tolerance = 0.01;
	IloNum LR_LB = -IloInfinity;
	IloNum LR_UB = IloInfinity;
	int iter = 0;
	while(LR_UB - LR_LB > tolerance*LR_LB  && iter < 5)
	{

		//cout<<"Hi"<<endl;
		IloNum sub_1 = 0.0;
	
		//Num2DMatrix r_val(env);	//Rh values for a priori cuts

		/*SUBPROBLEMS OF TYPE 1: SINGLE COMMODITY NETWORK FLOW PROBLEM*/
		for(int k = 0; k < K; k++)
		{
			//cout<<"Fine!"<<endl;
			/*VARIABLE DECLARATIONS*/
			NumVar3DMatrix P(env,WS_X);
			//Num3DMatrix P_val(env,WS_X);
			for(int i = 0; i < WS_X; i++)
			{
				P[i] = NumVar2DMatrix(env, WS_Y-1);
				//P_val[i] = Num2DMatrix(env, WS_Y-1);
				for(int j = 0; j < WS_Y-1; j++)
				{
					P[i][j] = IloNumVarArray(env, K);
					//P_val[i][j] = IloNumArray(env, K);
					for(int k = 0; k < K; k++)
					{
						P[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
					}
				}
			}

			NumVar3DMatrix Q(env,WS_X);
			//Num3DMatrix Q_val(env,WS_X);

			for(int i = 0; i < WS_X; i++)
			{
				Q[i] = NumVar2DMatrix(env, WS_Y-1);
				//Q_val[i] = Num2DMatrix(env, WS_Y-1);
				for(int j = 0; j < WS_Y-1; j++)
				{
					Q[i][j] = IloNumVarArray(env, K);
					//Q_val[i][j] = IloNumArray(env, K);
					for(int k = 0; k < K; k++)
					{
						Q[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
						//Q_val[i][j][k] = IloNum(env, 0.0, IloInfinity);
					}
				}
			}
			//cout<<"Fine P!"<<endl;
			NumVar3DMatrix U(env,WS_X-1);
			//Num3DMatrix U_val(env,WS_X-1);
			for(int i = 0; i < WS_X-1; i++)
			{
				U[i] = NumVar2DMatrix(env, WS_Y);
				//U_val[i] = Num2DMatrix(env, WS_Y);
				for(int j = 0; j < WS_Y; j++)
				{
					U[i][j] = IloNumVarArray(env, K);
					//U_val[i][j] = IloNumArray(env, K);
					for(int k = 0; k < K; k++)
					{
						U[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
						//U_val[i][j][k] = IloNum(env, 0.0, IloInfinity);
					}
				}
			}
			//cout<<"Fine! U "<<endl;
			NumVar3DMatrix V(env,WS_X-1);
			//Num3DMatrix V_val(env,WS_X-1);
			for(int i = 0; i < WS_X-1; i++)
			{
				V[i] = NumVar2DMatrix(env, WS_Y);
				//V_val[i] = Num2DMatrix(env, WS_Y);
				for(int j = 0; j < WS_Y; j++)
				{
					V[i][j] = IloNumVarArray(env, K);
					//V_val[i][j] = IloNumArray(env, K);
					for(int k = 0; k < K; k++)
					{
						V[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
					}
				}
			}
			std::cout<<"Variables declared!"<<std::endl;

			/*MODEL DESCRIPTION*/

			IloModel model(env);
		
			/*OBJECTIVE FUNCTION*/

			IloExpr obj(env);
			for(int i = 0; i < WS_X; i++)
			{
				for(int j = 0; j < WS_Y-1; j++)
				{
					//std::cout<<"Looping "<<k<<i<<j<<endl;
					obj += (f[k]*C + alpha_x[i][j][k])*(P[i][j][k] + Q[i][j][k]);// + U[i][j][k] + V[i][j][k]);
				}
			}
			for(int i = 0; i < WS_X-1; i++)
			{
				for(int j = 0; j < WS_Y; j++)
				{
					//std::cout<<"Looping "<<k<<i<<j<<endl;
					obj += (f[k]*C + alpha_y[i][j][k])*(U[i][j][k] + V[i][j][k]);
				}
			}
			model.add(IloMinimize(env,obj));
			obj.end();
			std::cout<<"Obj func declared!"<<std::endl;

			/*CONSTRAINT 1: ORIGIN, DESTINATION OR NORMAL NODE*/

			for(int i = 0; i < WS_X; i++)
			{
				for(int j = 0; j < WS_Y; j++)
				{
					IloExpr expr(env);
					//For input into node (i,j) from left side
					if((j-1) >=0 )
						expr += P[i][j-1][k];

					//For input into node (i,j) from right side
					if(j != WS_Y-1 )
						expr += Q[i][j][k];

					//For input into node (i,j) from top side
					if((i-1) >= 0 )
						expr += V[i-1][j][k];

					//For input into node (i,j) from bottom side
					if(i!= WS_X-1)
						expr += U[i][j][k];
					//std::std::cout<<"Pass"<<std::endl;
					//For output from node (i,j) to left side
					if((j-1) >=0)//  && j!=WS_Y-1)
						expr -= Q[i][j-1][k];
					//std::std::cout<<"Pass1"<<std::endl;
					//For output from node (i,j) to right side
					if(j!= WS_Y-1 )
						expr -= P[i][j][k];
					//std::std::cout<<"Pass2"<<std::endl;
					//For output from node (i,j) to top side
					if((i-1) >= 0)// && i!=WS_X-1)
						expr -= U[i-1][j][k];
					//std::std::cout<<"Pass3"<<std::endl;
					//For output from node (i,j) to bottom side
					if(i!= WS_X-1)
						expr -= V[i][j][k];
					//std::std::cout<<"Pass4"<<std::endl;
					//Check if (i,j) origin
					if(i==getX(O[k],WS_Y) && j==getY(O[k],WS_Y))
						model.add(expr == -1);
					else if(i==getX(D[k],WS_Y) && j==getY(D[k],WS_Y))
						model.add(expr == 1);
					else 
						model.add(expr == 0);
					//std::cout<<"At node: "<<k<<i<<j<<endl;
					expr.end();
				}
			}
			std::cout<<"Constraint 1 declared!"<<std::endl;

			IloCplex cplex(model);
		
			int res = cplex.solve();

			if(res)
			{
				float out = cplex.getObjValue();
				cout<<"Obj["<<k<<"]: "<<out<<endl;
				sub_1 += out;
				for(int i = 0; i < WS_X; i++)
				{
					for(int j = 0; j < WS_Y-1; j++)
					{
						P_val[i][j][k] = cplex.getValue(P[i][j][k]);
						Q_val[i][j][k] = cplex.getValue(Q[i][j][k]);
					}
				}
				for(int i = 0; i < WS_X-1; i++)
				{
					for(int j = 0; j < WS_Y; j++)
					{
						U_val[i][j][k] = cplex.getValue(U[i][j][k]);
						V_val[i][j][k] = cplex.getValue(V[i][j][k]);
					}
				}
			}
			else
			{
				cout<<"Couldn't solve for k: "<<k<<". Infeasible!"<<endl;
				return -1;
			}

			model.end();
		}

		cout<<"Sub_1: "<<sub_1<<endl;
	
		//system("pause");
		/*SUBPROBLEMS OF TYPE 2: NONLINEAR KNAPSACK PROBLEM*/
		IloNum sub_2x = 0.0;
		//IloNum sub_2x_LB = 0.0;
		IloNum sub_2y = 0.0;
		//IloNum sub_2y_LB = 0.0;
	
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				IloNum LB_l = IloInfinity;
				for(int l = 0; l < L; l++)
				{
					//cout<<"Fine"<<endl;
				
					//VARIABLE DECLARATIONS
				
					IloNumVar Rh(env,0,IloInfinity);
					IloNum Rh_val = 0;
				
					IloIntVarArray Wh(env,K);
					IloIntArray Wh_val(env, K);
					for(int k=0; k<K; k++)
					{
						Wh[k] = IloIntVar(env, 0, 1);
					}
						//}
					//}
					//cout<<"In"<<endl;
					//MODEL DESCRIPTION
					IloModel model(env);
					//cout<<"model"<<endl;
					//OBJECTIVE FUNCTION
					IloExpr obj(env);
					obj += Fx[i][j][l];
					obj += penalty*Rh;
					//cout<<"pen"<<endl;
					for(int k=0; k<K; k++)
					{
						//cout<<k<<endl;
						obj -= alpha_x[i][j][k]*Wh[k];
					}
					model.add(IloMinimize(env,obj));
					obj.end();
					//cout<<"Fine"<<endl;
					IloExpr lhs(env);
					for(int k=0; k<K; k++)
					{
						lhs += f[k]*Wh[k];
					}
					model.add(lhs <= Bx[i][j][l]);
					lhs.end();
					//cout<<"Fine"<<endl;
					IloNum UB_best = IloInfinity;	//Upper Bound
					IloNum UB = 0;		//Best Upper Bound
					IloNum LB = 0;					//Lower Bound
					IloNum tolerance = 0.01;
					IloInt iteration = 0;

					IloCplex cplex(model);
					cplex.setOut(env.getNullStream());
					//cout<<"Fine"<<endl;
					//Cutting Planes
					while ((UB_best - LB) > tolerance*LB && iteration < 100)//convex constraint not satisfied, ehnce continue
					{
						cout<<"========================================================================================="<<endl;
						cout<<"Iteration "<<iteration<<endl;
						cout<<"========================================================================================="<<endl;
						cplex.solve();//solving the MODEL
						if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
						{
							env.out() << "Sub Problem 2x infeasible!" << endl;
							return -1;
						}
						cout<<"Fine: "<<iteration<<endl;
						//cplex.getValues(X, X_val);
						Rh_val = cplex.getValue(Rh);
						LB = cplex.getObjValue();
						//IloNum LB = 0;//Lower Bound
						IloNum SumFlow_val = 0;
						UB = Fx[i][j][l];
						cout<<"UB: "<<UB<<endl;
						IloExpr cut_lhs(env);
						for(int k=0; k<K; k++)
						{
							Wh_val[k] = cplex.getValue(Wh[k]);
							UB -= alpha_x[i][j][k]*Wh_val[k];
							SumFlow_val+= f[k]*Wh_val[k];
							cut_lhs += f[k]*Wh[k];
						}
						cout<<"UB: "<<UB<<endl;
						cout<<"SumFlow: "<<SumFlow_val<<" Bx[i][j][l]: "<<Bx[i][j][l]<<endl;
						UB += penalty*SumFlow_val/(Bx[i][j][l] - SumFlow_val);
						UB_best = IloMin(UB, UB_best);
						cout<<"UB: "<<UB<<endl;
						cout<<"UB_best: "<<UB_best<<endl;
						cout<<"LB: "<<LB<<endl;
						//system("pause");
						
						IloNum cut_rhs = 0;
						IloNum R_new = 0;
						//R_new = SumFlow_val/(Bx[i][j][l] - SumFlow_val);
						cut_lhs -= (Bx[i][j][l]*Rh)/(IloPower(1.0 + Rh_val, 2));
						//cut_lhs-= (mu*R)/(IloPower(1+R_new, 2));
						cut_rhs = Bx[i][j][l]*IloPower((Rh_val)/(1+Rh_val), 2);
						//cut_rhs += (Bx[i][j][l]*(Rh[i][j] + R_new*R_new));
						//cout<<"Cut Added: Sum_i (\lambda_i*x_i) - "<<mu/(IloPower(1+R_val, 2))<<" R <= "<<mu*IloPower((R_val)/(1+R_val), 2)<<endl;
						//cout<<"New R Generated: "<<R_new<<endl;
						//cout<<"Cut Added: Sum_i (\lambda_i*x_i) - "<<mu/(IloPower(1+R_new, 2))<<" R <= "<<mu*IloPower((R_new)/(1+R_new), 2)<<endl;
						model.add(cut_lhs <= cut_rhs);
						cut_lhs.end();
						iteration+= 1;
					}
					//sub_2x_UB += UB_best;
					if(LB < LB_l)
					{
						LB_l = LB;
						P1_val[i][j][l] = 1;
						for(int k=0;k<K;k++)
						{
							Wh_final[i][j][k] = Wh_val[k];
						}
					}

					Rh.end();
					Wh.end();
					model.end();
					cplex.end();
				}
				sub_2x += LB_l;
			}
		}
		//cout<<"Sub_2x_UB: "<<sub_2x_UB<<endl;
		//cout<<"Sub_2x_LB: "<<sub_2x_LB<<endl;
		//system("pause");
	
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				IloNum LB_l = IloInfinity;
				for(int l = 0; l < L; l++)
				{				
					//VARIABLE DECLARATIONS
					cout<<"i: "<<i<<" j: "<<j<<" l: "<<l<<endl;
					IloNumVar Rv(env,0,IloInfinity);
					IloNum Rv_val = 0;
				
					IloIntVarArray Wv(env,K);
					IloIntArray Wv_val(env, K);
					for(int k=0; k<K; k++)
					{
						Wv[k] = IloIntVar(env, 0, 1);
					}
						//}
					//}
					//cout<<"In"<<endl;
					//MODEL DESCRIPTION
					IloModel model(env);
					//cout<<"model"<<endl;
					//OBJECTIVE FUNCTION
					IloExpr obj(env);
					obj += Fy[i][j][l];
					obj += penalty*Rv;
					//cout<<"pen"<<endl;
					for(int k=0; k<K; k++)
					{
						//cout<<k<<endl;
						obj -= alpha_y[i][j][k]*Wv[k];
					}
					model.add(IloMinimize(env,obj));
					obj.end();
					//cout<<"Fine"<<endl;
					IloExpr lhs(env);
					for(int k=0; k<K; k++)
					{
						lhs += f[k]*Wv[k];
					}
					model.add(lhs <= By[i][j][l]);
					lhs.end();
					//cout<<"Fine"<<endl;
					IloNum UB_best = IloInfinity;	//Upper Bound
					IloNum UB = 0;		//Best Upper Bound
					IloNum LB = 0;					//Lower Bound
					IloNum tolerance = 0.01;
					IloInt iteration = 0;

					IloCplex cplex(model);
					cplex.setOut(env.getNullStream());
					//cout<<"Fine"<<endl;
					//Cutting Planes
					while ((UB_best - LB) > tolerance*LB && iteration < 100)//convex constraint not satisfied, ehnce continue
					{
						cout<<"========================================================================================="<<endl;
						cout<<"Iteration "<<iteration<<endl;
						cout<<"========================================================================================="<<endl;
						cplex.solve();//solving the MODEL
						if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
						{
							env.out() << "Sub Problem 2y infeasible!" << endl;
							return -1;
						}
						//cout<<"Fine: "<<iteration<<endl;
						//cplex.getValues(X, X_val);
						Rv_val = cplex.getValue(Rv);
						LB = cplex.getObjValue();
						//IloNum LB = 0;//Lower Bound
						IloNum SumFlow_val = 0;
						UB = Fy[i][j][l];
						cout<<"UB: "<<UB<<endl;
						IloExpr cut_lhs(env);
						for(int k=0; k<K; k++)
						{
							Wv_val[k] = cplex.getValue(Wv[k]);
							UB -= alpha_y[i][j][k]*Wv_val[k];
							SumFlow_val+= f[k]*Wv_val[k];
							cut_lhs += f[k]*Wv[k];
						} 
						cout<<"UB: "<<UB<<endl;
						UB += penalty*SumFlow_val/(By[i][j][l] - SumFlow_val);
						UB_best = IloMin(UB, UB_best);
						cout<<"UB: "<<UB<<endl;
						cout<<"UB_best: "<<UB_best<<endl;
						cout<<"LB: "<<LB<<endl;
						//system("pause");
						
						IloNum cut_rhs = 0;
						IloNum R_new = 0;
						//R_new = SumFlow_val/(Bx[i][j][l] - SumFlow_val);
						cut_lhs -= (By[i][j][l]*Rv)/(IloPower(1.0 + Rv_val, 2));
						//cut_lhs-= (mu*R)/(IloPower(1+R_new, 2));
						cut_rhs = By[i][j][l]*IloPower((Rv_val)/(1+Rv_val), 2);
						//cut_rhs += (Bx[i][j][l]*(Rh[i][j] + R_new*R_new));
						//cout<<"Cut Added: Sum_i (\lambda_i*x_i) - "<<mu/(IloPower(1+R_val, 2))<<" R <= "<<mu*IloPower((R_val)/(1+R_val), 2)<<endl;
						//cout<<"New R Generated: "<<R_new<<endl;
						//cout<<"Cut Added: Sum_i (\lambda_i*x_i) - "<<mu/(IloPower(1+R_new, 2))<<" R <= "<<mu*IloPower((R_new)/(1+R_new), 2)<<endl;
						model.add(cut_lhs <= cut_rhs);
						cut_lhs.end();
						iteration+= 1;
						//cout<<"While loop! "<<iteration<<endl;
					}
					//sub_2y_UB += UB_best;
					
					cout<<"LB_l: "<<LB_l<<endl;
					if(LB < LB_l)
					{
						LB_l = LB;
						//cout<<"Here I was"<<endl;
						U1_val[i][j][l] = 1;
						//cout<<"Here I was"<<endl;
						for(int k=0;k<K;k++)
						{
							Wv_final[i][j][k] = Wv_val[k];
							//cout<<"Wv_final[i][j][k]: "<<Wv_final[i][j][k]<<endl;
						}
					}
					cout<<"Not dead!-2"<<endl;
					Rv.end();
					Wv.end();
					model.end();
					cplex.end();
				}
				sub_2y += LB_l;
			}
		}
		//cout<<"Not dead!-3"<<endl;
		cout<<"Sub1: "<<sub_1<<endl;
		cout<<"Sub_2x: "<<sub_2x<<endl;
		//cout<<"Sub_2x_LB: "<<sub_2x_LB<<endl;
		//cout<<"Sub_2y_UB: "<<sub_2y_UB<<endl;
		cout<<"Sub_2y: "<<sub_2y<<endl;
		float total = sub_1 + sub_2x + sub_2y;
		LR_LB = IloMax(LR_LB,total);
		cout<<"LR_LB: "<<LR_LB<<endl;
		
		IloExpr rhs_cons(env);
		IloExpr lhs_cons(env);
		lhs_cons += theta;
		//cout<<"OK"<<endl;
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				for(int k = 0; k < K; k++)
				{
					//cout<<"OK1"<<endl;
					rhs_cons += (f[k]*C )*(P_val[i][j][k] + Q_val[i][j][k]);// + U[i][j][k] + V[i][j][k]);
					lhs_cons -= omega_x[i][j][k]*(P_val[i][j][k] + Q_val[i][j][k]);
				}
			}
		}
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				for(int k = 0; k < K; k++)
				{
					cout<<"OK2"<<endl;
					rhs_cons += (f[k]*C)*(U_val[i][j][k] + V_val[i][j][k]);
					lhs_cons -= omega_y[i][j][k]*(U_val[i][j][k] + V_val[i][j][k]);
				}
			}
		}
		
		//std::cout<<"Total variable cost: "<<UB<<endl;
		//std::cout<<"out"<<endl;
		//For all (i,j) in A'
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				for(int l = 0;l<L;l++)
				{
					rhs_cons += Fx[i][j][l]*(P1_val[i][j][l]);// + penalty*(Rh_val[i][j]);// + U1[i][j]);
					//std::cout<<"UB_X: "<<UB<<endl;
				}
			}
		}
		
		//std::cout<<"U1: "<<U1_val[0][0]<<endl;
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				for(int l = 0;l<L;l++)
				{
					rhs_cons += Fy[i][j][l]*(U1_val[i][j][l]);// + penalty*(Rv_val[i][j]);
					//std::cout<<"UB_Y: "<<UB<<endl;
					
				}
			}
		}

		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				double flow1 = 0;
				for(int k = 0; k < K; k++)
				{
					flow1+= f[k]*(P_val[i][j][k] + Q_val[i][j][k]);
				}
				double arc_cap1 = 0;
				for(int l = 0;l<L;l++)
				{
					arc_cap1 += Bx[i][j][l]*P1_val[i][j][l];
				}
				
				if(abs(flow1 - arc_cap1) <= tolerance && flow1)
				{
					rhs_cons += 9999999;
					//cout<<"I was here!"<<endl;
				}
				else if(flow1)
				{	
					//td::cout<<"Flow: "<<flow1<<" Cap: "<<arc_cap1<<endl;
					//std::cout<<"Penalty[i][j]: "<<abs(penalty*(flow1/(arc_cap1-flow1)))<<endl;
					//std::cout<<"UB_x b4 = "<<UB<<endl;
					rhs_cons += abs(penalty*(flow1/(arc_cap1-flow1)));
					//std::cout<<"UB_x = "<<UB<<endl;
				}
			}
		}

		//std::cout<<"UB_x out= "<<UB<<endl;
		//UB = 0;
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				//std::cout<<"UB_y in= "<<UB<<endl;
				double flow2=0;
				for(int k = 0; k < K; k++)
				{
					flow2+= f[k]*(U_val[i][j][k] + V_val[i][j][k]);
				}
				double arc_cap2 = 0;
				for(int l = 0;l<L;l++)
				{
					arc_cap2 += By[i][j][l]*U1_val[i][j][l];
				}
				//std::cout<<"Flow: "<<flow2<<" Cap: "<<arc_cap2<<endl;
				if(abs(flow2 - arc_cap2) < tolerance && flow2)
				{
					rhs_cons += 9999999;
					//cout<<"I was here!"<<endl;
				}
				else if(flow2)
				{	
					//std::cout<<"Flow: "<<flow2<<" Cap: "<<arc_cap2<<endl;
					//std::cout<<"Penalty[i][j]: "<<abs(penalty*(flow2/(arc_cap2-flow2)))<<endl;
					//std::cout<<"UB_y b4= "<<UB<<endl;
					rhs_cons += abs(penalty*(flow2/(arc_cap2-flow2)));
					//std::cout<<"UB_y = "<<UB<<endl;
				}
			}
		}
		

		model_LR.add(lhs_cons <= rhs_cons);
		lhs_cons.end();
		rhs_cons.end();
		//IloNum bound = getUB();
		//IloNum rhs_val = 0;
		//rhs_val = getRhs();
		//model_LR.add(theta <= rhs_val);
		cplex_LR.extract(model_LR);
		//solving the MODEL
		if (!cplex_LR.solve()) // if the problem is infeasible
		{
			env.out() << "Lagrangian Infeasible!" << endl;
			return -1;
		}
		cplex_LR.exportModel("abc.lp");
		LR_UB = IloMin(LR_UB,cplex_LR.getObjValue());
		cout<<"LR_UB: "<<LR_UB<<endl;
		cout<<"Theta: "<<cplex_LR.getValue(theta);
		for(int i = 0; i < WS_X; i++)
		{
			//alpha_x[i] = Num2DMatrix(env,WS_Y-1);
			for(int j = 0; j < WS_Y-1; j++)
			{
				//alpha_x[i][j] = IloNumArray(env,K);
				for(int k = 0; k < K; k++)
				{
					IloNum val = cplex_LR.getValue(omega_x[i][j][k]);
					//cout<<"Fine: "<<val<<endl;
					alpha_x[i][j][k] = val;
				}
			}
		}
		
		for(int i = 0; i < WS_X-1; i++)
		{
			//alpha_y[i] = Num2DMatrix(env,WS_Y);
			for(int j = 0; j < WS_Y; j++)
			{
				//alpha_y[i][j] = IloNumArray(env,K);
				for(int k = 0; k < K; k++)
				{
					IloNum val = cplex_LR.getValue(omega_y[i][j][k]);
					//cout<<"Fine: "<<val<<endl;
					alpha_y[i][j][k] = val;
				}
			}
		}
		cout<<"-----------------------------"<<endl;
		cout<<"LR_UB: "<<LR_UB<<endl;
		cout<<"LR_LB: "<<LR_LB<<endl;
		iter ++;

		for(int i = 0; i < WS_X; i++)
		{
			//alpha_x[i] = Num2DMatrix(env,WS_Y-1);
			for(int j = 0; j < WS_Y-1; j++)
			{
				//alpha_x[i][j] = IloNumArray(env,K);
				//for(int k = 0; k < K; k++)
				//{
					cout<<alpha_x[i][j][0]<<" ";
				//}
			}
			cout<<endl;
		}
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				//for(int k = 0; k < K; k++)
				//{
					cout<<alpha_y[i][j][0]<<" ";
				//}
			}
			cout<<endl;
		}
		system("pause");
	}
	model_LR.end();
	cplex_LR.end();
	theta.end();
	omega_x.end();
	omega_y.end();
	alpha_x.end();
	alpha_y.end();
	env.end();

	return 0;
}
