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
	IloEnv env;
	clock_t t1,t2,t3;
	t1=clock();
	/*FILE READING*/
	IloNumArray r_val(env);	//Rh values for a priori cuts

	/*const char* data_filename = "Test1.txt";
	const char* filename = "C:/Users/Partha Mishra/Documents/Visual Studio 2012/Projects/MontrealProject/x64/Debug/r_val_0.001.dat";
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
	*/
	/*PARAMETER READING*/

	int K;			//Total number of connections
	int L = 3;			//Number of capacity levels
	int WS_X;		//Workspace size
	int WS_Y;		//Workspace size
	//int F;			//Fixed cost of each arc
	//int B;			//Capacity of ach arc
	int C = 1;			//Variable cost for each arc
	float penalty = 1000;	//Unit penalty at each arc for congestion
	cout<<"WS_X: ";
	cin>>WS_X;
	cout<<"WS_Y: ";
	cin>>WS_Y;
	cout<<"Number of commodity flows(K): ";
	cin>>K;
	//datafile>>WS_X>>WS_Y>>K>>L>>C>>penalty;

	IloNumArray fix(env, L);	//Fixed cost of all arcs at different capacity levels
	IloNumArray cap(env, L);	//Capacity of all arcs at different levels
	IloNumArray f(env, K);	//Flow/demand of each connection
	IloIntArray O(env, K);	//Origin of each connection
	IloIntArray D(env, K);	//Destination of each connection

	Num3DMatrix Fx(env,WS_X);
	Num3DMatrix Fy(env,WS_X-1);
	Num3DMatrix Bx(env,WS_X);
	Num3DMatrix By(env,WS_X-1);

	//datafile>>fix>>cap>>f>>O>>D;
	fix[0] = 100;
	fix[1] = 500;
	fix[2] = 1000;
	
	cap[0] = 5;
	cap[1] = 10;
	cap[2] = 15;



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

	//std::cout<<K<<":"<<fix<<":"<<cap<<":"<<penalty<<endl;
	//std::cout<<f<<endl;
	//std::cout<<O<<endl;
	//std::cout<<D<<endl;
	//system("pause");
	//system("pause");
	

	IloNum UB_best = IloInfinity;//Upper Bound
	IloNum LB = 0;//Best Lower Bound
	//Num2DMatrix r_val(env);	//Rh values for a priori cuts
	
	/*VARIABLE DECLARATIONS*/
	NumVar2DMatrix Rh(env,WS_X);	//For horizontal arcs
	Num2DMatrix Rh_val(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		Rh[i] = IloNumVarArray(env, WS_Y-1);
		Rh_val[i] = IloNumArray(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Rh[i][j] = IloNumVar(env, 0.0, IloInfinity);
		}
	}
	
	NumVar2DMatrix Rv(env,WS_X-1);	//For vertical arcs
	Num2DMatrix Rv_val(env,WS_X-1);

	for(int i = 0; i < WS_X-1; i++)
	{
		Rv[i] = IloNumVarArray(env, WS_Y);
		Rv_val[i] = IloNumArray(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			Rv[i][j] = IloNumVar(env, 0.0, IloInfinity);
		}
	}
	//std::cout<<"R declared"<<endl;
	NumVar3DMatrix Zh(env,WS_X);	//For linearization of horizontal arcs
	Num3DMatrix Zh_val(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		Zh[i] = NumVar2DMatrix(env, WS_Y-1);
		Zh_val[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Zh[i][j] = IloNumVarArray(env, L);
			Zh_val[i][j] = IloNumArray(env, L);
			for(int k = 0; k < L; k++)
			{
				Zh[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
			}
		}
	}

	NumVar3DMatrix Zv(env,WS_X-1);	//For linearization of vertical arcs
	Num3DMatrix Zv_val(env,WS_X-1);

	for(int i = 0; i < WS_X-1; i++)
	{
		Zv[i] = NumVar2DMatrix(env, WS_Y);
		Zv_val[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			Zv[i][j] = IloNumVarArray(env, L);
			Zv_val[i][j] = IloNumArray(env, L);
			for(int k = 0; k < L; k++)
			{
				Zv[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
			}
		}
	}

	NumVar3DMatrix P(env,WS_X);
	Num3DMatrix P_val(env,WS_X);
	for(int i = 0; i < WS_X; i++)
	{
		P[i] = NumVar2DMatrix(env, WS_Y-1);
		P_val[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			P[i][j] = IloNumVarArray(env, K);
			P_val[i][j] = IloNumArray(env, K);
			for(int k = 0; k < K; k++)
			{
				P[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
			}
		}
	}
	//std::cout<<"flag!"<<endl;

	IntVar3DMatrix P1(env,WS_X);
	Int3DMatrix P1_val(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		P1[i] = IntVar2DMatrix(env, WS_Y-1);
		P1_val[i] = Int2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			P1[i][j] = IloIntVarArray(env, L);
			P1_val[i][j] = IloIntArray(env, L);
			for(int k = 0; k < L; k++)
			{
				P1[i][j][k] = IloIntVar(env, 0, 1);
				//P1_val[i][j] = IloIntArray(env, K);
			}
		}
	}

	NumVar3DMatrix Q(env,WS_X);
	Num3DMatrix Q_val(env,WS_X);

	for(int i = 0; i < WS_X; i++)
	{
		Q[i] = NumVar2DMatrix(env, WS_Y-1);
		Q_val[i] = Num2DMatrix(env, WS_Y-1);
		for(int j = 0; j < WS_Y-1; j++)
		{
			Q[i][j] = IloNumVarArray(env, K);
			Q_val[i][j] = IloNumArray(env, K);
			for(int k = 0; k < K; k++)
			{
				Q[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
				//Q_val[i][j][k] = IloNum(env, 0.0, IloInfinity);
			}
		}
	}

	NumVar3DMatrix U(env,WS_X-1);
	Num3DMatrix U_val(env,WS_X-1);
	for(int i = 0; i < WS_X-1; i++)
	{
		U[i] = NumVar2DMatrix(env, WS_Y);
		U_val[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			U[i][j] = IloNumVarArray(env, K);
			U_val[i][j] = IloNumArray(env, K);
			for(int k = 0; k < K; k++)
			{
				U[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
				//U_val[i][j][k] = IloNum(env, 0.0, IloInfinity);
			}
		}
	}

	IntVar3DMatrix U1(env,WS_X-1);
	Int3DMatrix U1_val(env,WS_X-1);

	for(int i = 0; i < WS_X-1; i++)
	{
		U1[i] = IntVar2DMatrix(env, WS_Y);
		U1_val[i] = Int2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			U1[i][j] = IloIntVarArray(env, L);
			U1_val[i][j] = IloIntArray(env, L);
			for(int k = 0; k < L; k++)
			{
				U1[i][j][k] = IloIntVar(env, 0, 1);
			}
		}
	}

	NumVar3DMatrix V(env,WS_X-1);
	Num3DMatrix V_val(env,WS_X-1);
	for(int i = 0; i < WS_X-1; i++)
	{
		V[i] = NumVar2DMatrix(env, WS_Y);
		V_val[i] = Num2DMatrix(env, WS_Y);
		for(int j = 0; j < WS_Y; j++)
		{
			V[i][j] = IloNumVarArray(env, K);
			V_val[i][j] = IloNumArray(env, K);
			for(int k = 0; k < K; k++)
			{
				V[i][j][k] = IloNumVar(env, 0.0, IloInfinity);
			}
		}
	}
	std::cout<<"Variables declared!"<<std::endl;
	/* initialize random seed: */
	srand (time(NULL));
	int R = 0;
	vector <int> resource;
	int rand_num;
	for(int i=0;i<WS_X-3;i+=3)
	{
		for(int j=0;j<WS_Y-3;j+=2)
		{
			R++;
			
			rand_num = rand()%4;
			//cout<<"Rand: "<<rand_num<<endl;
			if(rand_num == 0)
				resource.push_back(WS_Y*(i+1)+j+1);
			if(rand_num == 1)
				resource.push_back(WS_Y*(i+1)+j+2);
			if(rand_num == 2)
				resource.push_back(WS_Y*(i+2)+j+1);
			if(rand_num == 3)
				resource.push_back(WS_Y*(i+2)+j+2);
		}
	}
	
	int rand1,rand2;
	for(int i=0;i<K;i++)
	{
		//srand (time(NULL));
		rand1 = rand()%R;
		//srand (time(NULL));
		rand2 = rand()%R;
		while(rand2==rand1)
			rand2 = rand()%R;
		O[i] = resource[rand1];
		D[i] = resource[rand2];
		f[i] = 1;
	}
	//cout<<f<<endl;
	cout<<"Origin/Destination"<<endl;
	cout<<O<<endl;
	cout<<D<<endl;
	system("pause");
	
	
	//system("pause");
	/*MODEL DESCRIPTION*/

	IloModel model(env);

	/*OBJECTIVE FUNCTION*/

	IloExpr obj(env);
	//For all combinations of (i,j) in A
	for(int k = 0; k < K; k++)
	{
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				//std::cout<<"Looping "<<k<<i<<j<<endl;
				obj += f[k]*C*(P[i][j][k] + Q[i][j][k]);// + U[i][j][k] + V[i][j][k]);
			}
		}
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				//std::cout<<"Looping "<<k<<i<<j<<endl;
				obj += f[k]*C*(U[i][j][k] + V[i][j][k]);
			}
		}
	}
	//std::cout<<"out"<<endl;
	//For all (i,j) in A'
	for(int i = 0; i < WS_X; i++)
	{
		for(int j = 0; j < WS_Y-1; j++)
		{
			for(int k = 0; k < L; k++)
			{
				obj += Fx[i][j][k]*(P1[i][j][k]);// + U1[i][j]);
			}
			obj += + penalty*(Rh[i][j]);
		}
	}
	for(int i = 0; i < WS_X-1; i++)
	{
		for(int j = 0; j < WS_Y; j++)
		{
			for(int k = 0; k < L; k++)
			{
				obj += Fy[i][j][k]*(U1[i][j][k]);
			}
			obj += + penalty*(Rv[i][j]);
		}
	}
	model.add(IloMinimize(env,obj));
	std::cout<<"Obj func declared!"<<std::endl;
	//system("pause");
	/*CONSTRAINT 1: ORIGIN, DESTINATION OR NORMAL NODE*/
	//For all nodes
	for(int k = 0; k < K; k++)
	{
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
			}
		}
	}
	std::cout<<"Constraint 1 declared!"<<std::endl;
	//system("pause");
	/*CONSTRAINT 2: ARC CAPACITY*/
	//For (i,j) in A'
	//FOR HORIZONTAL ARCS
	for(int i = 0; i < WS_X; i++)
	{
		for(int j = 0; j < WS_Y-1; j++)
		{
			IloExpr expr1(env);
			//IloExpr expr2(env);
			for(int k = 0; k < K; k++)
			{
				expr1 += f[k]*(P[i][j][k] + Q[i][j][k]);
				//std::std::cout<<"At node: "<<i<<j<<k<<std::endl;
				//expr2 += f[k]*(U[i][j][k] + V[j][i][k]);
			}
			IloExpr expr2(env);
			for(int k = 0; k < L; k++)
			{
				expr2 += Bx[i][j][k]*Zh[i][j][k];
				//std::std::cout<<"At node: "<<i<<j<<k<<std::endl;
				//expr2 += f[k]*(U[i][j][k] + V[j][i][k]);
			}
			model.add(expr1 == expr2);
			//model.add(expr2 <= B);
			expr1.end();
			expr2.end();
		}
	}
	//std::std::cout<<"Pass_HA"<<std::endl;
	//FOR VERTICAL ARCS
	for(int i = 0; i < WS_X-1; i++)
	{
		for(int j = 0; j < WS_Y; j++)
		{
			IloExpr expr2(env);
			for(int k = 0; k < K; k++)
			{
				//expr1 += f[k]*(P[i][j][k] + Q[j][i][k]);
				expr2 += f[k]*(U[i][j][k] + V[i][j][k]);
			}
			IloExpr expr1(env);
			for(int k = 0; k < L; k++)
			{
				expr1 += By[i][j][k]*Zv[i][j][k];
				//std::std::cout<<"At node: "<<i<<j<<k<<std::endl;
				//expr2 += f[k]*(U[i][j][k] + V[j][i][k]);
			}
			//model.add(expr1 <= B);
			model.add(expr2 == expr1);
			//expr1.end();
			expr2.end();
			expr1.end();
		}
	}
	std::cout<<"Constraint 2 declared!"<<std::endl;
	//system("pause");

	/*CONSTRAINT 3: LINK CONSTRAINT*/
	//For (i,j) in A'
	for(int i = 0; i < WS_X; i++)
	{
		for(int j = 0; j < WS_Y-1; j++)
		{
			for(int k = 0; k < L; k++)
			{
				model.add(Zh[i][j][k] <= P1[i][j][k]);
			}
			/*for(int k = 0; k < K; k++)
			{
				//model.add(P[i][j][k] <= P1[i][j]);
				//model.add(Q[i][j][k] <= P1[i][j]);
			}*/
		}
	}

	for(int i = 0; i < WS_X-1; i++)
	{
		for(int j = 0; j < WS_Y; j++)
		{
			for(int k = 0; k < L; k++)
			{
				model.add(Zv[i][j][k] <= U1[i][j][k]);
			}
			/*for(int k = 0; k < K; k++)
			{
				model.add(U[i][j][k] <= U1[i][j]);
				model.add(V[i][j][k] <= U1[i][j]);
			}*/
		}
	}
	std::cout<<"Constraint 3 declared!"<<std::endl;
	//system("pause");

	//CONSTRAINT 4: ARC CONSTRAINTS
	
	for(int i=1;i<WS_X-2;i+=3)
	{
		for(int j=1;j<WS_Y-2;j+=2)
		{
			for(int k = 0; k < L; k++)
			{
				model.add(P1[i][j][k] == 0);
				model.add(P1[i+1][j][k] == 0);
				model.add(U1[i][j][k] == 0);
				model.add(U1[i][j+1][k] == 0);
			}
		}
	}

	std::cout<<"Constraint 4 declared!"<<std::endl;

	//CONSTRAINT 5: SINGLE LEVEL ALLOWED ON EACH ARC
	
	for(int i = 0; i < WS_X; i++)
	{
		for(int j = 0; j < WS_Y-1; j++)
		{
			IloExpr expr1(env);
			for(int k = 0; k < L; k++)
			{
				expr1 += P1[i][j][k];
			}
			model.add(expr1<=1);
			expr1.end();
		}
	}

	for(int i = 0; i < WS_X-1; i++)
	{
		for(int j = 0; j < WS_Y; j++)
		{
			IloExpr expr1(env);
			for(int k = 0; k < L; k++)
			{
				expr1 += U1[i][j][k];
			}
			model.add(expr1<=1);
			expr1.end();
		}
	}
	
	std::cout<<"Constraint 5 declared!"<<std::endl;
	
	//CONSTRAINT 6: A PRIORI CUTS
	//For horizontal arcs
	for (IloInt k=0; k<r_val.getSize(); k++)
	{
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				for(int l = 0; l < L; l++)
				{
					model.add((1+r_val[k])*(1+r_val[k])*Zh[i][j][l] - Rh[i][j] <= r_val[k]*r_val[k]);
				}
			}
		}
	}
	//For vertical arcs
	for (IloInt k=0; k<r_val.getSize(); k++)
	{
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				for(int l = 0; l < L; l++)
				{
					model.add((1+r_val[k])*(1+r_val[k])*Zv[i][j][l] - Rv[i][j] <= r_val[k]*r_val[k]);
				}
			}
		}
	}
	
	std::cout<<"Constraint 6 declared!"<<std::endl;
	
	/*SOLVE THE MODEL*/
	IloCplex cplex(env);
	//cplex.exportModel("cap.lp");
	//system("pause");
	IloNum tolerance = 0.1;
	IloInt iteration = 0;
	double UB;
	cplex.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
	cplex.setParam(IloCplex::EpInt, 0);
	
	float fixed_cost = 0;
	float variable_cost = 0;
	float penalty_cost = 0;

	while ((UB_best - LB) > tolerance && iteration < 100 && (float(t1) / CLOCKS_PER_SEC) < 10800)//convex constraint not satisfied, hence continue
	{
		std::cout<<"========================================================================================="<<endl;
		std::cout<<"Iteration "<<iteration<<endl;
		std::cout<<"=========================================================================================="<<endl;
		cplex.extract(model);
		//cplex.exportModel("cap.lp");
		float elapsed = (float(t1) / CLOCKS_PER_SEC);
		cplex.setParam(IloCplex::TiLim,10800 - elapsed);
		cplex.solve();//solving the MODEL
		if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
		{
			env.out() << "Problem Infeasible" << endl; 
			return -1;
		}
		std::cout<<"Solved!"<<endl;
		//cplex.getValues(X, X_val);
		
		//READING VALUES OF VARIABLES
		/*for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				cplex.getValue(P1[i][j],P1_val[i][j]);
				cplex.getValue(U1[i][j],U1_val[i][j]);
				cplex.getValue(Zh[i][j],Zh_val[i][j]);
				cplex.getValue(Zv[i][j],Zv_val[i][j]);
				//std::cout<<cplex.getValue(P1[i][j])<<endl;
				//std::cout<<cplex.getValue(U1[i][j])<<endl;
			
				cplex.getValues(P[i][j],P_val[i][j]);
				cplex.getValues(Q[i][j],Q_val[i][j]);
				cplex.getValues(U[i][j],U_val[i][j]);
				cplex.getValues(V[i][j],V_val[i][j]);			
			}
		}*/
		
		for(int i=0;i<WS_X;i++)
		{
			for(int j=0;j<WS_Y-1;j++)
			{
				for(int k = 0;k<L;k++)
				{
					P1_val[i][j][k] = cplex.getValue(P1[i][j][k]);
					Zh_val[i][j][k] = cplex.getValue(Zh[i][j][k]);
				}
				
				cplex.getValues(P[i][j],P_val[i][j]);
				
				cplex.getValues(Q[i][j],Q_val[i][j]);
				//std::cout<<"Got this!"<<endl;
				Rh_val[i][j] = cplex.getValue(Rh[i][j]);
				//std::cout<<"P: "<<P_val[i][j]<<endl<<"Q: "<<Q_val[i][j]<<endl<<"Rh: "<<Rh_val[i][j]<<endl;
			}
		}
		for(int i=0;i<WS_X-1;i++)
		{
			for(int j=0;j<WS_Y;j++)
			{
				for(int k = 0;k<L;k++)
				{
					U1_val[i][j][k] = cplex.getValue(U1[i][j][k]);
					Zv_val[i][j][k] = cplex.getValue(Zv[i][j][k]);
				}
				cplex.getValues(U[i][j],U_val[i][j]);
				cplex.getValues(V[i][j],V_val[i][j]);
				Rv_val[i][j] = cplex.getValue(Rv[i][j]);
				//std::cout<<"U: "<<U_val[i][j]<<endl<<"U1: "<<U1_val[i][j]<<endl<<"V: "<<V_val[i][j]<<endl<<"Rv: "<<Rv_val[i][j]<<endl;
			}
		}
		
		std::cout<<endl<<endl<<endl;
		std::cout<<"Horizontal arc matrix"<<endl<<endl;
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				int flag = 0;
				for(int k = 0;k<L;k++)
				{
					if(cplex.getValue(P1[i][j][k]))
					{
						std::cout<<"-----"<<k+1<<"|||||";
						flag = 1;
					}
				}
				if(flag == 0)
					std::cout<<"-----"<<0<<"|||||";
			}
			std::cout<<endl;
		}
		std::cout<<"\nVertical arc matrix"<<endl<<endl;
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				int flag = 0;
				for(int k = 0;k<L;k++)
				{
					if(cplex.getValue(U1[i][j][k]))
					{
						std::cout<<"-----"<<k+1<<"|||||";
						flag = 1;
					}
					//std::cout<<cplex.getValue(U1[i][j][k])
				}
				if(flag == 0)
					std::cout<<"-----"<<0<<"|||||";
			}
			std::cout<<endl;
		}
	
		LB = cplex.getObjValue();
		/*IloNum LB = 0;//Lower Bound
		IloNum SumLambda_val = 0;
		for (IloInt i=0; i<N; i++)
		{
			LB+= profit[i]*X_val[i];
			SumLambda_val+= lambda[i]*X_val[i];
		} 
		LB-= cost_wait*SumLambda_val/(mu - SumLambda_val);
		*/

		//COMPUTING UB
		UB = 0.0;//Upper Bound

		for(int k = 0; k < K; k++)
		{
			for(int i = 0; i < WS_X; i++)
			{
				for(int j = 0; j < WS_Y-1; j++)
				{
					//std::cout<<"Looping "<<k<<i<<j<<endl;
					UB += f[k]*C*(P_val[i][j][k] + Q_val[i][j][k]);// + U[i][j][k] + V[i][j][k]);
				}
			}
			//std::cout<<"UB1: "<<UB<<endl;
			for(int i = 0; i < WS_X-1; i++)
			{
				for(int j = 0; j < WS_Y; j++)
				{
					//std::cout<<"Looping "<<k<<i<<j<<endl;
					UB += f[k]*C*(U_val[i][j][k] + V_val[i][j][k]);
				}
			}
		}
		variable_cost = UB;
		//std::cout<<"Total variable cost: "<<UB<<endl;
		//std::cout<<"out"<<endl;
		//For all (i,j) in A'
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				for(int l = 0;l<L;l++)
				{
					UB += Fx[i][j][l]*(P1_val[i][j][l]);// + penalty*(Rh_val[i][j]);// + U1[i][j]);
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
					UB += Fy[i][j][l]*(U1_val[i][j][l]);// + penalty*(Rv_val[i][j]);
					//std::cout<<"UB_Y: "<<UB<<endl;
					
				}
			}
		}
		fixed_cost = UB - variable_cost;
		//std::cout<<"Total variable + fixed  cost: "<<UB<<endl;
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
					UB += 9999999;
					//cout<<"I was here!"<<endl;
				}
				else if(flow1)
				{	
					//td::cout<<"Flow: "<<flow1<<" Cap: "<<arc_cap1<<endl;
					//std::cout<<"Penalty[i][j]: "<<abs(penalty*(flow1/(arc_cap1-flow1)))<<endl;
					//std::cout<<"UB_x b4 = "<<UB<<endl;
					UB += abs(penalty*(flow1/(arc_cap1-flow1)));
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
					UB += 9999999;
					//cout<<"I was here!"<<endl;
				}
				else if(flow2)
				{	
					//std::cout<<"Flow: "<<flow2<<" Cap: "<<arc_cap2<<endl;
					//std::cout<<"Penalty[i][j]: "<<abs(penalty*(flow2/(arc_cap2-flow2)))<<endl;
					//std::cout<<"UB_y b4= "<<UB<<endl;
					UB += abs(penalty*(flow2/(arc_cap2-flow2)));
					//std::cout<<"UB_y = "<<UB<<endl;
				}
			}
		}
		penalty_cost = UB - fixed_cost - variable_cost;
		//std::cout<<"UB3: "<<UB<<endl;

		//UB = computeUB();
		UB_best = IloMin(UB, UB_best);
		// Print results
		//std::cout<<"X : "<<X_val<<endl;
		//std::cout<<"Number of Customer Locations Served = "<<IloSum(X_val)<<endl;
		//std::cout<<"R : "<<R_val<<endl;
		
		std::cout<<"UB = "<<UB<<endl;
		float min_util = 0;
		float max_util = 0;
		float avg_util = 0;
		float util = 0;
		int num_arcs = 0;
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				//cplex.getValue(P1[i][j],P1_val[i][j]);
				//cplex.getValue(U1[i][j],U1_val[i][j]);
				//std::cout<<"P1["<<i<<"]["<<j<<"]"<<endl;
				//std::cout<<P1_val[i][j]<<endl;
				IloNum flow = 0;
				for(int k = 0; k < K; k++)
				{
					flow += f[k]*(P_val[i][j][k] + Q_val[i][j][k]);
				}
				IloNum cap = 0;
				for(int l = 0; l < L; l++)
				{
					cap += Bx[i][j][l]*P1_val[i][j][l];
					if(P1_val[i][j][l])
						num_arcs ++;
				}

				util = flow / cap;
				if(cap != 0)
					avg_util += util;
				if(util < min_util)
					min_util = util;
				else if(util > max_util)
					max_util = util;
				//if(flow!=0)
				//	std::cout<<"Utilisation: "<<flow/cap<<endl;
				//else std::cout<<"Utilisation: "<<0<<endl;
			
				//cplex.getValues(P[i][j],P_val[i][j]);
				//cplex.getValues(Q[i][j],Q_val[i][j]);
				//cplex.getValues(U[i][j],U_val[i][j]);
				//cplex.getValues(V[i][j],V_val[i][j]);
			
				//std::cout<<"P["<<i<<"]["<<j<<"]"<<endl;
				//std::cout<<P_val[i][j]<<endl;
				//std::cout<<"Q["<<i<<"]["<<j<<"]"<<endl;
				//std::cout<<Q_val[i][j]<<endl;
				//std::cout<<"U["<<i<<"]["<<j<<"]"<<endl;
				//std::cout<<U_val[i][j]<<endl;
				//std::cout<<"V["<<i<<"]["<<j<<"]"<<endl;
				//std::cout<<V_val[i][j]<<endl;
			
			}
		}
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				//std::cout<<"U1["<<i<<"]["<<j<<"]"<<endl;
				//std::cout<<U1_val[i][j]<<endl;
				IloNum flow = 0;
				for(int k = 0; k < K; k++)
				{
					flow += f[k]*(U_val[i][j][k] + V_val[i][j][k]);
				}
				IloNum cap = 0;
				for(int l = 0; l < L; l++)
				{
					cap += By[i][j][l]*U1_val[i][j][l];
					if(U1_val[i][j][l])
						num_arcs ++;
				}
				util = flow / cap;
				if(cap != 0)
					avg_util += util;
				if(util < min_util)
					min_util = util;
				else if(util > max_util)
					max_util = util;
					//if(flow!=0)
					//	std::cout<<"Utilisation: "<<flow/cap<<endl;
					//else std::cout<<"Utilisation: "<<0<<endl;
			}
		}
	
		avg_util /= num_arcs;
		t3=clock();
		float diff ((float)t3-(float)t1);
		float seconds = diff / CLOCKS_PER_SEC;

		//std::cout<<"TOTAL CPUTIME = "<<seconds<<" secs"<<std::endl;


		std::cout<<endl<<"REPORT"<<std::endl;
		std::cout<<"GRID SIZE: "<<WS_X<<"*"<<WS_Y<<std::endl;
		std::cout<<"LEVELS: "<<L<<std::endl;
		std::cout<<"COMMODITIES: "<<K<<std::endl;
		std::cout<<"FIXED COST: "<<fixed_cost<<std::endl;
		std::cout<<"VARIABLE COST: "<<variable_cost<<std::endl;
		std::cout<<"PENALTY COST: "<<penalty_cost<<std::endl;
		std::cout<<"UPPER BOUND: "<<UB_best<<std::endl;
		std::cout<<"LOWER BOUND: "<<LB<<std::endl;
		std::cout<<"GAP: "<<UB_best - LB<<std::endl;
		std::cout<<"MIN ARC UTILISATION: "<<min_util<<std::endl;
		std::cout<<"MAX ARC UTILISATION: "<<max_util<<std::endl;
		std::cout<<"AVERAGE ARC UTILISATION: "<<avg_util<<std::endl;
		std::cout<<"NUMBER OF ITERATIONS: "<<iteration<<std::endl;
		std::cout<<"NUMBER OF ARCS ADDED: "<<num_arcs<<std::endl;
		std::cout<<"CPUTIME = "<<seconds<<" secs"<<std::endl;
		//system("pause");
		//ADDING CUTTING PLANES (CONSTRAINTS ON THE FLY)
		
		for(int i = 0; i < WS_X; i++)
		{
			for(int j = 0; j < WS_Y-1; j++)
			{
				for(int k = 0;k<L;k++)
				{
					if(P1_val[i][j][k])
					{
						IloExpr lhs(env);
						lhs+=(1+ Rh_val[i][j])*(1+ Rh_val[i][j])*Zh[i][j][k];
						lhs-=Rh[i][j];
						IloExpr rhs(env);
						rhs+=Rh_val[i][j]*Rh_val[i][j];
						model.add(lhs<=rhs);
						lhs.end();
						rhs.end();
					}
					//std::cout<<(1+ Rh_val[i][j])*(1+ Rh_val[i][j])<<"*Zh["<<i<<"]["<<j<<"] + Rh["<<i<<"]["<<j<<"] <= "<<Rh_val[i][j]*Rh_val[i][j]<<"		added"<<endl;
				}
			}
		}
		for(int i = 0; i < WS_X-1; i++)
		{
			for(int j = 0; j < WS_Y; j++)
			{
				for(int k = 0;k<L;k++)
				{
					if(U1_val[i][j][k])
					{
						IloExpr lhs(env);
						lhs+=(1+ Rv_val[i][j])*(1+ Rv_val[i][j])*Zv[i][j][k];
						lhs-=Rv[i][j];
						IloExpr rhs(env);
						rhs+=Rv_val[i][j]*Rv_val[i][j];
						model.add(lhs<=rhs);
						//cout<<"Adding constraint!";
						lhs.end();
						rhs.end();
					}
					//std::cout<<(1+ Rv_val[i][j])*(1+ Rv_val[i][j])<<"*Zv["<<i<<"]["<<j<<"] + Rv["<<i<<"]["<<j<<"] <= "<<Rv_val[i][j]*Rv_val[i][j]<<"		added"<<endl;
				}
			}
		}
		/*
		IloExpr cut_lhs(env);
		IloNum cut_rhs = 0;
		for (IloInt i=0; i<N; i++)
		{
			cut_lhs+= lambda[i]*X[i];
		}
		IloNumArray R_new(env);
		R_new = SumLambda_val/(mu - SumLambda_val);
		//cut_lhs-= (mu*R)/(IloPower(1+R_val, 2));
		cut_lhs-= (mu*R)/(IloPower(1+R_new, 2));
		//cut_rhs = mu*IloPower((R_val)/(1+R_val), 2);
		cut_rhs = mu*IloPower((R_new)/(1+R_new), 2);
		//std::cout<<"Cut Added: Sum_i (\lambda_i*x_i) - "<<mu/(IloPower(1+R_val, 2))<<" R <= "<<mu*IloPower((R_val)/(1+R_val), 2)<<endl;
		std::cout<<"New R Generated: "<<R_new<<endl;
		std::cout<<"Cut Added: Sum_i (\lambda_i*x_i) - "<<mu/(IloPower(1+R_new, 2))<<" R <= "<<mu*IloPower((R_new)/(1+R_new), 2)<<endl;
		model.add(cut_lhs <= cut_rhs);
		cut_lhs.end();
		//t2=clock();*/
		//float diff ((float)t2-(float)t1);
		//float seconds = diff / CLOCKS_PER_SEC;
		//std::cout<<"TOTAL CPUTIME = "<<seconds<<" secs"<<endl;
		iteration+= 1;
	}
	//cplex.exportModel("LastModel.lp");
	/*int res = cplex.solve();
	cplex.exportModel("Test.lp");
	if(!res)
	{
		std::cout<<"PROBLEM INFEASIBLE!"<<endl;
		return -1;
	}*/
	t2=clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;

	std::cout<<"TOTAL CPUTIME = "<<seconds<<" secs"<<std::endl;
	
	//DISPLAY RESULTS
	//std::cout<<"Objective function value: "<<cplex.getObjValue()<<endl;
	/*for(int i = 0; i < WS_X; i++)
	{
		for(int j = 0; j < WS_Y-1; j++)
		{
			//cplex.getValue(P1[i][j],P1_val[i][j]);
			//cplex.getValue(U1[i][j],U1_val[i][j]);
			std::cout<<"P1["<<i<<"]["<<j<<"]"<<endl;
			std::cout<<P1_val[i][j]<<endl;
			IloNum flow = 0;
			for(int k = 0; k < K; k++)
			{
				flow += f[k]*(P_val[i][j][k] + Q_val[i][j][k]);
			}
			IloNum cap = 0;
			for(int l = 0; l < L; l++)
			{
				cap += Bx[i][j][l]*P1_val[i][j][l];
			}
			if(flow!=0)
				std::cout<<"Utilisation: "<<flow/cap<<endl;
			else std::cout<<"Utilisation: "<<0<<endl;
			
			cplex.getValues(P[i][j],P_val[i][j]);
			cplex.getValues(Q[i][j],Q_val[i][j]);
			cplex.getValues(U[i][j],U_val[i][j]);
			cplex.getValues(V[i][j],V_val[i][j]);
			
			std::cout<<"P["<<i<<"]["<<j<<"]"<<endl;
			std::cout<<P_val[i][j]<<endl;
			std::cout<<"Q["<<i<<"]["<<j<<"]"<<endl;
			std::cout<<Q_val[i][j]<<endl;
			std::cout<<"U["<<i<<"]["<<j<<"]"<<endl;
			std::cout<<U_val[i][j]<<endl;
			std::cout<<"V["<<i<<"]["<<j<<"]"<<endl;
			std::cout<<V_val[i][j]<<endl;
			
		}
	}
	for(int i = 0; i < WS_X-1; i++)
	{
		for(int j = 0; j < WS_Y; j++)
		{
			std::cout<<"U1["<<i<<"]["<<j<<"]"<<endl;
			std::cout<<U1_val[i][j]<<endl;
			IloNum flow = 0;
			for(int k = 0; k < K; k++)
			{
				flow += f[k]*(U_val[i][j][k] + V_val[i][j][k]);
			}
			IloNum cap = 0;
			for(int l = 0; l < L; l++)
			{
				cap += By[i][j][l]*U1_val[i][j][l];
			}
			if(flow!=0)
				std::cout<<"Utilisation: "<<flow/cap<<endl;
			else std::cout<<"Utilisation: "<<0<<endl;
		}
	}*/
	return 0;
}