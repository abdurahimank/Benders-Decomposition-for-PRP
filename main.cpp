////2 Stage - Production Routing Problem (PRP)////
// - A single product.
// - A single production plant.
// - "n" number of customers.
// - "t" number of time periods.
// - "m" number of vehicles.


#pragma warning(disable : 4996) //For Visual Studio 2012
#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <iosfwd>
#include <string>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector> //for vectors
#include <math.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

using namespace std;

typedef IloArray<IloNumArray> TwoDMatrix;

ILOSTLBEGIN

int main(int argc, char** argv)
{
	IloEnv env;
	try
	{
		//////////Part 1 - DEFINE VARIABLES AND PARAMETERS//////////
		
		/////CONSTANTS/////
		int N = 5;  // Number of customers + Plant (# of Plant = 1).
		int T = 3;  // Number of time periods.
		int M = 1;  // Number of vehicles available.

		int C = 100;  // Production capacity of the plant in a time period.
		int Q = 50;  // Vehicle capacity (all vehicles are identical).
		double u = 10;  // Unit cost of production.
		double f = 5;  // Setup cost for production.

		
		/////PARAMETERS/////
		TwoDMatrix demand(env);  // Demand at a customer node in a time period.
		TwoDMatrix tranport_cost(env);  // Transportation cost between two nodes.
		IloNumArray holding_cost(env);  // Inventory holding cost at each node.
		IloNumArray penalty(env);  // Penalt at node i, if demand is unmet in a period.

		IloNumArray init_Inventory(env);  // Initial inventory at each node.
		IloNumArray invemtory_cap(env);  // Inventory capacity at each node (node 0 is the plant).






		/////BASIC DECISION VARIABLES/////
		IloNumVarArray Y(env, T, 0, IloInfinity, ILOBOOL);  // 1 if production takesplace in time period T, ow = 0.








		IloNumVarArray X_dual(env, 2, 0, IloInfinity, ILOFLOAT);
		IloNumArray Y_val(env, 3);
		IloNum theta_val;


		/////DECISION VARIABLES AND PARAMETERS FOR EXTREM RAY PROBLEM/////
		IloNumVarArray X_dual_er(env, 2, 0, IloInfinity, ILOFLOAT);
		//IloNumVarArray Y_val_er(env, 3, 0, 5, ILOINT);
		//IloNum theta_val;


		/////DECISION VARIABLES AND PARAMETERS FOR MASTER PROBLEM/////
		IloNumVarArray Y(env, 3, 0, 5, ILOINT);
		IloNumArray X_dual_val(env, 2);
		IloNumArray X_dual_val_er(env, 2);
		//IloNumVar theta_var(env, -IloInfinity, IloInfinity, ILOFLOAT);
		IloNumVar theta_var(env, 0, IloInfinity, ILOFLOAT);





		//////////Part 2 - DEVELOP GENERIC MODEL//////////

		/////SET MASTER PROBLEM/////
		IloModel model_master(env);
		IloExpr Objective_master(env);
		Objective_master = theta_var - 5 * Y[0] + 2 * Y[1] - 9 * Y[2];
		model_master.add(IloMinimize(env, Objective_master));
		IloCplex cplex_master(env);
		Objective_master.end();
		cplex_master.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_master.setWarning(env.getNullStream()); // This is to supress warning messages on screen

		/////SET SUBPROBLEM (DUAL FORMULATION)/////
		IloModel model_sub(env);
		IloObjective Objective_sub = IloMaximize(env);
		model_sub.add(Objective_sub);

		model_sub.add(-2 * X_dual[0] - 3 * X_dual[1] <= -2);
		model_sub.add(-3 * X_dual[0] + X_dual[1] <= 3);
		model_sub.add(-6 * X_dual[0] - 3 * X_dual[1] <= -4);


		IloCplex cplex_sub(model_sub);
		IloNum eps = cplex_sub.getParam(IloCplex::EpInt);//Integer tolerance for MIP models; 
		//default value of EpInt remains 1e-5 http://www.iro.umontreal.ca/~gendron/IFT6551/CPLEX/HTML/relnotescplex/relnotescplex12.html
		cplex_sub.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_sub.setWarning(env.getNullStream()); //This is to supress warning messages on screen


		/////EXTREME RAY PROBLEM (DUAL FORMULATION)/////
		IloModel model_sub_er(env);
		IloObjective Objective_sub_er = IloMaximize(env);
		model_sub_er.add(Objective_sub_er);

		model_sub_er.add(-2 * X_dual_er[0] - 3 * X_dual_er[1] <= 0);
		model_sub_er.add(-3 * X_dual_er[0] + X_dual_er[1] <= 0);
		model_sub_er.add(-6 * X_dual_er[0] - 3 * X_dual_er[1] <= 0);
		model_sub_er.add(X_dual_er[0] + X_dual_er[1] == 1);


		IloCplex cplex_sub_er(model_sub_er);
		cplex_sub_er.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_sub_er.setWarning(env.getNullStream()); //This is to supress warning messages on screen






		//////////Part 3 - BEGIN ITERATIONS//////////

		IloNum GAP = IloInfinity;
		theta_val = 0;
		Y_val[0] = 0;
		Y_val[1] = 0;
		Y_val[2] = 0;
		IloNum sub_obj_val = 0;
		IloNum er_obj_val = 0;
		IloNum Upper_bound = IloInfinity;
		IloNum Lower_bound = 0;
		GAP = Upper_bound - Lower_bound;
		cout << "Y = " << Y_val << endl;
		IloInt Iter = 0;

		//while( Iter < MaxCut )
		while (Upper_bound - Lower_bound > eps)
		{
			Iter++;
			cout << "=========================================" << endl;
			cout << "============ITERATION " << Iter << "==============" << endl;
			//Define Object Function for the Dual of the Sub problem
			IloExpr sub_obj(env);
			sub_obj = (2 + 5 * Y_val[0] - 3 * Y_val[1] + 7 * Y_val[2]) * X_dual[0] + (-10 + 4 * Y_val[0] + 2 * Y_val[1] + 4 * Y_val[2]) * X_dual[1];
			Objective_sub.setExpr(IloMaximize(env, sub_obj));

			//Define Object Function for the Dual of the Extreme ray problem
			IloExpr sub_obj_er(env);
			sub_obj_er = (2 + 5 * Y_val[0] - 3 * Y_val[1] + 7 * Y_val[2]) * X_dual_er[0] + (-10 + 4 * Y_val[0] + 2 * Y_val[1] + 4 * Y_val[2]) * X_dual_er[1];
			Objective_sub_er.setExpr(IloMaximize(env, sub_obj_er));


			cplex_sub.setParam(cplex_sub.PreInd, 0);   //Disable presolve, otherwise, if dual is infeasible, 
			//we don't know if prime is unbounded or infeasible
			cplex_sub.setParam(IloCplex::RootAlg, IloCplex::Primal);//Solve the SP Dual using Primal Simplex
			cout << "SOLVING SUB PROBLEM" << endl;
			cplex_sub.solve();


			cout << "Sub Problem Solution Status: " << cplex_sub.getCplexStatus() << endl;
			if (cplex_sub.getCplexStatus() == CPX_STAT_OPTIMAL)
			{// Dual subproblem is bounded; Add Optimality Cut to the Master Problem
				cplex_sub.getValues(X_dual_val, X_dual);  // taking values of X_dual from SP and saves to X_dual_val
				cout << "X_dual = " << X_dual_val << endl;
				sub_obj_val = cplex_sub.getObjValue();
				cout << "sub_obj_val = " << sub_obj_val << endl;
				Upper_bound = IloMin(Upper_bound, (-5 * Y_val[0] + 2 * Y_val[1] - 9 * Y_val[2] + sub_obj_val));
				cout << "Upper_bound = " << Upper_bound << endl;

				//Add Cut to the Master Problem
				//cout << "Optimality Cut Added to Master Problem: " << "theta + " << (-5 * X_dual_val[0] - 4 * X_dual_val[1]) << " Y1 + " << (3 * X_dual_val[0] - 2 * X_dual_val[2]) << " Y2 + " << (-7 * X_dual_val[0] - 4 * X_dual_val[1]) << " Y3 >= " << (2 * X_dual_val[0]) - (10 * X_dual_val[1]) << endl;
				model_master.add(theta_var + (-5 * X_dual_val[0] - 4 * X_dual_val[1]) * Y[0] + (3 * X_dual_val[0] - 2 * X_dual_val[1]) * Y[1] + (-7 * X_dual_val[0] - 4 * X_dual_val[1]) * Y[2] >=
					2 * X_dual_val[0] - 10 * X_dual_val[1]);
				cout << "Optimality Cut Added to Master Problem: " << "theta + " << (-5 * X_dual_val[0] - 4 * X_dual_val[1]) << " Y1 + "
					<< (3 * X_dual_val[0] - 2 * X_dual_val[1]) << " Y2 + " << (-7 * X_dual_val[0] - 4 * X_dual_val[1]) << " Y3 >= " << 2 * X_dual_val[0] - 10 * X_dual_val[1] << endl;

			}

			if (cplex_sub.getCplexStatus() == CPX_STAT_UNBOUNDED)
			{// Dual subproblem is unbounded; Add Optimality Cut to the Master Problem
				cplex_sub.getValues(X_dual_val, X_dual);  // taking values of X_dual from SP and saves to X_dual_val
				cout << "X_dual = " << X_dual_val << endl;
				sub_obj_val = cplex_sub.getObjValue();
				cout << "sub_obj_val = " << sub_obj_val << endl;
				Upper_bound = IloMin(Upper_bound, (-5 * Y_val[0] + 2 * Y_val[1] - 9 * Y_val[2] + sub_obj_val));
				cout << "Upper_bound = " << Upper_bound << endl;

				//Add Cut to the Master Problem
				cout << "Optimality Cut Added to Master Problem: " << "theta + " << (-5 * X_dual_val[0] - 4 * X_dual_val[1]) << " Y1 + "
					<< (3 * X_dual_val[0] - 2 * X_dual_val[1]) << " Y2 + " << (-7 * X_dual_val[0] - 4 * X_dual_val[1]) << " Y3 >= " << 2 * X_dual_val[0] - 10 * X_dual_val[1] << endl;
				model_master.add(theta_var + (-5 * X_dual_val[0] - 4 * X_dual_val[1]) * Y[0] + (3 * X_dual_val[0] - 2 * X_dual_val[1]) * Y[1] + (-7 * X_dual_val[0] - 4 * X_dual_val[1]) * Y[2] >=
					2 * X_dual_val[0] - 10 * X_dual_val[1]);



				//////SOLVING EXTREME RAYS PROBLEM////
				cplex_sub_er.solve();
				cout << "Extreme Ray Problem Solution Status: " << cplex_sub_er.getCplexStatus() << endl;

				cout << "SOLVING EXTREME RAY PROBLEM" << endl;
				// Dual subproblem is unbounded; Hence add feasibility Cut to the Master Problem
				cplex_sub_er.getValues(X_dual_val_er, X_dual_er);  // taking values of X_dual from SP and saves to X_dual_val
				cout << "X_dual of extreme rays = " << X_dual_val_er << endl;
				//sub_obj_val = cplex_sub.getObjValue();
				//cout << "extreme_ray_obj_val = " << sub_obj_val << endl;
				//Upper_bound = IloMin(Upper_bound, (-5 * Y_val[0] + 2 * Y_val[1] - 9 * Y_val[2] + sub_obj_val));
				//cout << "Upper_bound = " << Upper_bound << endl;

				//Add Cut to the Master Problem
				cout << "Feasibility Cut Added to Master Problem: " << "0 + " << (-5 * X_dual_val_er[0] - 4 * X_dual_val_er[1]) << " Y1 + "
					<< (3 * X_dual_val_er[0] - 2 * X_dual_val_er[1]) << " Y2 + " << (-7 * X_dual_val_er[0] - 4 * X_dual_val_er[1]) << " Y3 >= " << 2 * X_dual_val_er[0] - 10 * X_dual_val_er[1] << endl;
				model_master.add(0 + (-5 * X_dual_val_er[0] - 4 * X_dual_val_er[1]) * Y[0] + (3 * X_dual_val_er[0] - 2 * X_dual_val_er[1]) * Y[1] + (-7 * X_dual_val_er[0] - 4 * X_dual_val_er[1]) * Y[2] >=
					2 * X_dual_val_er[0] - 10 * X_dual_val_er[1]);
			}




			cout << "SOLVING MASTER PROBLEM" << endl;
			cout << "Master Problem Solution Status: " << cplex_master.getCplexStatus() << endl;
			cplex_master.extract(model_master);
			if (!cplex_master.solve())
			{
				cout << "Failed" << endl;
				throw(-1);
			}
			Y_val[0] = cplex_master.getValue(Y[0]);
			Y_val[1] = cplex_master.getValue(Y[1]);
			Y_val[2] = cplex_master.getValue(Y[2]);
			theta_val = cplex_master.getValue(theta_var);
			cout << "theta_var = " << theta_val << endl;
			cout << "Y1 = " << Y_val[0] << endl;
			cout << "Y2 = " << Y_val[1] << endl;
			cout << "Y3 = " << Y_val[2] << endl;
			Lower_bound = theta_val - 5 * Y_val[0] + 2 * Y_val[1] - 9 * Y_val[2];
			cout << "Lower_bound = " << Lower_bound << endl;
		}//while(Upper_bound - Lower_bound > eps)
		model_master.end();
		model_sub.end();
		cplex_master.end();
		cplex_sub.end();
	}//try
	catch (IloException& e)
	{
		env.out() << "ERROR: " << e << endl;
	}
	catch (...)
	{
		env.out() << "Unknown exception" << endl;
	}
	env.end();
	return 0;
}