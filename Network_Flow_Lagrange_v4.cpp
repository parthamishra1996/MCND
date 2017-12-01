/*
Version 1: Multiplier values were input manually 
Version 2: Cutting plane used but problem became unbounded
Version3: Solving simpler MCND using explicit upper bound
Version4: Use of sub gradient method
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

int WS_X;		//Workspace size
int WS_Y;		//Workspace size
int K;			//Total number of connections
int L;			//Number of capacity levels
int C;			//Variable cost for each arc
float penalty;	//Unit penalty at each arc for congestion
IloNumArray fix(env, L);	//Fixed cost of all arcs at different capacity levels
IloNumArray cap(env, L);	//Capacity of all arcs at different levels
IloNumArray f(env, K);	//Flow/demand of each connection
IloIntArray O(env, K);	//Origin of each connection
IloIntArray D(env, K);	//Destination of each connection

IloNumArray r_val(env);	//Rh values for a priori cuts

IloCplex cplex(env);

IloNum eps = 0.01;

int getX(int n, int WS_Y)
{
	return n/WS_Y;
}

int getY(int n, int WS_Y)
{
	return n%WS_Y;
}

IloNum solve(Num3DMatrix alpha_x, Num3DMatrix alpha_y, Num3DMatrix P_val, Num3DMatrix Q_val, Num3DMatrix U_val, Num3DMatrix V_val, Int3DMatrix Wh_vals, Int3DMatrix Wv_vals)
{
	clock_t t1,t2;
	t1=clock();

	std::cout<<K<<":"<<fix<<":"<<cap<<":"<<penalty<<endl;
	std::cout<<f<<endl;
	std::cout<<O<<endl;
	std::cout<<D<<endl;
	
	Num3DMatrix Fx(env,WS_X);
	Num3DMatrix Fy(env,WS_X-1);
	Num3DMatrix Bx(env,WS_X);
	Num3DMatrix By(env,WS_X-1);
	
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
				//cout<<i<<j<<l<<endl;
				Fx[i][j][l] = fix[l];
				Bx[i][j][l] = cap[l];
				//cout<<Fx[i][j][l]<<":"<<Bx[i][j][l]<<endl;
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
				//cout<<i<<j<<l<<endl;
				Fy[i][j][l] = fix[l];
				By[i][j][l] = cap[l];
			}
		}
	}
	

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
				//for(int k = 0; k < K; k++)
				//{					
				P[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
				//}
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
				//for(int k = 0; k < K; k++)
				//{					
				Q[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
					//Q_val[i][j][k] = IloNum(env, 0.0, IloInfinity);
				//}
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
				//for(int k = 0; k < K; k++)
				//{
				U[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
					//U_val[i][j][k] = IloNum(env, 0.0, IloInfinity);
				//}
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
				//for(int k = 0; k < K; k++)
				//{
				V[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
				//}
			}
		}
		std::cout<<"Variables declared!"<<std::endl;

		/*MODEL DESCRIPTION*/

		IloModel model(env);
		
		/*OBJECTIVE FUNCTION*/
		//cout<<"check1"<<endl;
		IloExpr obj(env);
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				//std::cout<<"Looping1: "<<k<<i<<j<<endl;
				obj += (f[k]*C + alpha_x[i][j][k])*(P[i][j][k] + Q[i][j][k]);// + U[i][j][k] + V[i][j][k]);
			}
		}
		//cout<<"check2"<<endl;
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				//std::cout<<"Looping2: "<<k<<i<<j<<endl;
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
		//IloCplex cplex(model);
		//cout<<"check0"<<endl;
		cplex.extract(model);
		//cout<<"check1"<<endl;
		int res = cplex.solve();
		//cout<<"check2"<<endl;
		if(res)
		{
			//cout<<"check3"<<endl;
			float out = cplex.getObjValue();
			//cout<<"Obj["<<k<<"]: "<<out<<endl;
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
				//cout<<i<<":"<<j<<":"<<l<<endl;
				
				//VARIABLE DECLARATIONS
				
				IloNumVar Rh(env,0,IloInfinity);
				IloNum Rh_val = 0;
				//cout<<"check1"<<endl;
				IloIntVarArray Wh(env,K);
				IloIntArray Wh_val(env, K);
				for(int k=0; k<K; k++)
				{
					Wh[k] = IloIntVar(env, 0, 1);
				}
				//cout<<"check2"<<endl;
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
				//cout<<"check3"<<endl;
				//cout<<"Fine"<<endl;
				IloExpr lhs(env);
				//cout<<"Fine1"<<endl;
				for(int k=0; k<K; k++)
				{
					//cout<<"Fine2: "<<endl;
					lhs += f[k]*Wh[k];
				}
				//cout<<"Fine: "<<endl;
				model.add(lhs <= Bx[i][j][l]);
				lhs.end();
				
				IloNum UB_best = IloInfinity;	//Upper Bound
				IloNum UB = 0;		//Best Upper Bound
				IloNum LB = 0;					//Lower Bound
				IloNum tolerance = 0.01;
				IloInt iteration = 0;
				//cout<<"check4"<<endl;
				cplex.extract(model);
				//cout<<"check5"<<endl;
				//IloCplex cplex(model);
				cplex.setOut(env.getNullStream());
				//cout<<"Fine"<<endl;
				//Cutting Planes
				
				while (IloAbs(UB_best - LB) > tolerance*IloAbs(LB) && iteration < 5)//convex constraint not satisfied, ehnce continue
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
					//cout<<"Fine: "<<iteration<<endl;
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
						Wh_vals[i][j][k] = Wh_val[k];
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
				LB_l = IloMin(LB_l,LB);

				Rh.end();
				Wh.end();
				model.end();
				//cplex.end();
			}
			sub_2x += LB_l;
		}
	}
	
	cout<<"Sub_2x: "<<sub_2x<<endl;
	
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
				while ((UB_best - LB)/LB > tolerance && iteration < 5)//convex constraint not satisfied, ehnce continue
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
					//cout<<"Fine1"<<endl;
					for(int k=0; k<K; k++)
					{
						//cout<<"Fine2"<<endl;
						Wv_val[k] = cplex.getValue(Wv[k]);
						//cout<<"Fine3"<<endl;
						Wv_vals[i][j][k] = Wv_val[k];
						//cout<<"Fine4"<<endl;
						UB -= alpha_y[i][j][k]*Wv_val[k];
						//cout<<"Fine5"<<endl;
						SumFlow_val+= f[k]*Wv_val[k];
						//cout<<"Fine6"<<endl;
						cut_lhs += f[k]*Wv[k];
						//cout<<"Fine7"<<endl;
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
					//cout<<"Not dead!-1"<<endl;
				}
				//sub_2y_UB += UB_best;
				LB_l = IloMin(LB_l,LB);
				//cout<<"Not dead!-2"<<endl;
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
	double total = sub_1 + sub_2x + sub_2y;
	cout<<"Total: "<<total<<endl;
	//alpha_x.end();
	//alpha_y.end();
	return total;
}

IloInt checkConstraints(Num3DMatrix P_val, Num3DMatrix Q_val, Num3DMatrix U_val, Num3DMatrix V_val, Int3DMatrix Wh_vals, Int3DMatrix Wv_vals, Num3DMatrix alpha_x, Num3DMatrix alpha_y, IloNum &sum_product_expr, IloNum &sum_squared_expr)
{
	cout<<"Constraints getting checked!"<<endl;
	IloInt status = 1;
	IloNum expr = 0;
	for(int i=0;i<WS_X;i++)
	{
		for(int j=0;j<WS_Y-1;j++)
		{
			for(int k=0;k<K;k++)
			{
				expr = P_val[i][j][k] + Q_val[i][j][k] - Wh_vals[i][j][k];
				if(expr > eps)
					status = 0;
				sum_product_expr += alpha_x[i][j][k]*expr;
				sum_squared_expr += expr*expr;
			}
		}
	}
	for(int i=0;i<WS_X-1;i++)
	{
		for(int j=0;j<WS_Y;j++)
		{
			for(int k=0;k<K;k++)
			{
				expr = U_val[i][j][k] + V_val[i][j][k] - Wv_vals[i][j][k];
				if(expr > eps)
					status = 0;
				sum_product_expr += alpha_y[i][j][k]*expr;
				sum_squared_expr += expr*expr;
			}
		}
	}
	cout<<"Constraints checked!"<<endl;
	cout<<"Status: "<<status<<endl;
	return status;
}


int main(int argc, char **argv)
{
	/*FILE READING*/

	const char* data_filename = "Test_small.txt";
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
	
	//int F;			//Fixed cost of each arc
	//int B;			//Capacity of ach arc
	
	datafile>>WS_X>>WS_Y>>K>>L>>C>>penalty;

	datafile>>fix>>cap;
	datafile>>f>>O>>D;

	std::cout<<WS_X<<":"<<WS_Y<<":"<<K<<":"<<L<<":"<<penalty<<endl;
	std::cout<<f<<endl;
	std::cout<<O<<endl;
	std::cout<<D<<endl;
	//cout<<Bx<<endl;
	//cout<<By<<endl;
	cout<<fix<<endl;
	cout<<cap<<endl;
	//IloIntArray X_vals(env, 4);

	

	//cout<<Bx[0][0][0]<<endl;
	//cout<<By[0][1][2]<<endl;
	


	//Multipliers
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
				alpha_x[i][j][k] = 1000;
				//cout<<alpha_x[i][j][k]<<endl;
			}
		}
	}
	
	for(int i = 0; i < WS_X-1; i++)
	{
		alpha_y[i] = Num2DMatrix(env,WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			alpha_y[i][j] = IloNumArray(env,K);
			for(int k = 0; k < K; k++)
			{
				alpha_y[i][j][k] = 100;
			}
		}
	}

	
	//Variable values

	Num3DMatrix P_val(env,WS_X);
	for(int i = 0; i < WS_X; i++)
	{
		P_val[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
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
	
	Num3DMatrix V_val(env,WS_X-1);
	for(int i = 0; i < WS_X-1; i++)
	{
		V_val[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			V_val[i][j] = IloNumArray(env, K);
		}
	}

	Int3DMatrix Wh_vals(env,WS_X);
	
	for(int i = 0; i < WS_X; i++)
	{
		Wh_vals[i] = Int2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Wh_vals[i][j] = IloIntArray(env, K);
		}
	}

	Int3DMatrix Wv_vals(env,WS_X-1);
	for(int i = 0; i < WS_X-1; i++)
	{
		Wv_vals[i] = Int2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			Wv_vals[i][j] = IloIntArray(env, K);
		}
	}


	IloNum t = 1, L_curr = 0, L_opt = IloInfinity, Z_opt = 0, del = 0.60;
	IloInt iter = 0, iterLimit = 20, imp = 0, impLimit = 3;
	IloNum sum_product_expr = 0, sum_squared_expr = 0, expr = 0;
	do
	{
		L_curr = solve(alpha_x, alpha_y, P_val, Q_val, U_val, V_val, Wh_vals, Wv_vals);
		imp++;
		//cout<<"Iter: "<<iter<<"		u: "<<u<<"		L_curr: "<<L_curr<<"		Z_opt: "<<Z_opt<<endl;
		//expr = 10 - 8*X_vals[0] - 2*X_vals[1] - X_vals[2]  - 4*X_vals[3];
		//cout<<"Expr: "<<expr<<endl;
		if(checkConstraints(P_val,Q_val,U_val,V_val,Wh_vals,Wv_vals,alpha_x,alpha_y,sum_product_expr,sum_squared_expr) && (L_curr - sum_product_expr) > Z_opt)	//if expr<eps and Z < Zopt then change Zopt
			Z_opt = L_curr - sum_product_expr;					//if all constraints are satisfied and Z<Zopt then update
		if(L_curr < L_opt)							//if Lcurr < Lopt then change Lopt
		{
			L_opt = L_curr;
			imp = 0;
		}
		if(imp == impLimit)
		{
			del/=2.0;
			imp = 0;
		}

		//u_prev = u;
		//u = IloMax(0,(u_prev - t*(sum_expr)));		//remains same
		for(int i=0;i<WS_X;i++)
		{
			for(int j=0;j<WS_Y-1;j++)
			{
				for(int k=0;k<K;k++)
				{
					expr = P_val[i][j][k] + Q_val[i][j][k] - Wh_vals[i][j][k];
					alpha_x[i][j][k] = alpha_x[i][j][k] - t*expr;
				}
			}
		}
		for(int i=0;i<WS_X-1;i++)
		{
			for(int j=0;j<WS_Y;j++)
			{
				for(int k=0;k<K;k++)
				{
					expr = U_val[i][j][k] + V_val[i][j][k] - Wv_vals[i][j][k];
					alpha_y[i][j][k] = alpha_y[i][j][k] - t*expr;						
				}
			}
		}
		//cout<<"alphax:"<<endl;
		//cout<<alpha_x<<endl;
		//cout<<"alphay:"<<endl;
		//cout<<alpha_y<<endl;
		t = (L_curr - Z_opt)*del/(sum_squared_expr);	//sum of squared expression for each constraint
		iter ++;
		cout<<"L_curr: "<<L_curr<<endl;
		cout<<"L_opt: "<<L_opt<<endl;
		cout<<"Z_opt: "<<Z_opt<<endl;
		cout<<"Iter: "<<iter<<" Del: "<<del<<" t: "<<t<<endl;
	}while( iter < iterLimit && del > eps && IloAbs((Z_opt - L_curr)) > eps*Z_opt); //expr*expr > eps && t > eps && 

	cout<<"L_curr: "<<L_curr<<endl;

	return 0;
}