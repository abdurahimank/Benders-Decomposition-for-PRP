////2 Stage - Production Routing Problem (PRP)////
// - A single product.
// - A single production plant.
// - "n" number of customers.
// - "t" number of time periods.
// - "m" number of vehicles.
//Refer Main Paper - 
//Refer Supplementary Material - 


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
#include <random>
#include <math.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

using namespace std;

typedef IloArray<IloNumArray> Array2D;  // Defining a normal 2D array.
typedef IloArray<Array2D> Array3D;  // Defining a normal 3D array.
typedef IloArray<Array3D> Array4D;  // Defining a normal 4D array.


typedef IloArray<IloNumVarArray> VarArray2D;  // Defining a 2D array of decision variables
typedef IloArray<VarArray2D> VarArray3D;  // Defining a 3D array of decision variables
typedef IloArray<VarArray3D> VarArray4D;  // Defining a 4D array of decision variables


ILOSTLBEGIN

int main(int argc, char** argv)
{
	IloEnv env;
	try
	{
		//////////Part 1 - DEFINE VARIABLES AND PARAMETERS//////////

#pragma region Defining Constants
		/////CONSTANTS/////
		int N = 10;  // Number of customers + Plant (# of Plant = 1).
		int T = 3;  // Number of time periods.
		int K = 1;  // Number of vehicles available.

		int C = 100;  // Production capacity of the plant in a time period.
		int Q = 50;  // Vehicle capacity (all vehicles are identical).

		double u = 10;  // Unit cost of production.
		double f = 500;  // Setup cost for production.

		vector<pair<int, int>> E;  // Set of edges, such that i < j.
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				E.emplace_back(i, j);  // store edge (i,j) only when i < j
			}
		}

#pragma endregion


#pragma region Defining Parameters and Reading Values


		/////PARAMETERS/////
		Array2D demand(env);  // Demand at a customer node in a time period.

		Array2D tranport_cost(env);  // Transportation cost between two nodes.
		IloNumArray holding_cost(env);  // Inventory holding cost at each node.
		IloNumArray penalty(env);  // Penalty at node i, if demand is unmet in a period.

		IloNumArray init_inventory(env);  // Initial inventory at each node.
		IloNumArray inventory_cap(env);  // Inventory capacity at each node (node 0 is the plant).


		///READING PARAMETER DATA///
		const char* data_filename = "Data.dat";
		if (argc > 1)
		{
			data_filename = argv[1];
		}
		fstream datafile;
		datafile.open(data_filename, ios::in);

		if (!datafile)
		{
			cerr << "ERROR: could not open file " << data_filename << " for reading" << endl;
			cerr << "usage:   " << argv[0] << " <datafile>" << endl;
			throw(-1);
		}

		//Importing data in the form of Matrix in the .dat file into the variables as per order.
		datafile >> demand >> tranport_cost >> holding_cost >> penalty >> init_inventory >> inventory_cap;

		
		IloBool consistentData = (demand.getSize() == tranport_cost.getSize());  // Ensuring nodes are same in different matrices.
		if (!consistentData)
		{
			cerr << "ERROR: data file '" << data_filename << "' contains inconsistent data" << endl;
			throw(-1);
		}
		datafile.close();

		// READING DONE......

		
		///DISPLAYING PARAMETER DATA///
		//Printing Demand
		cout << "Demand at each node in each time period: " << endl;
		cout << "[";
		for (int i = 0; i < N; ++i) {
			cout << "[";
			for (int t = 0; t < T; ++t) {
				if (t < T - 1) {
					cout << demand[i][t] << ", ";
				}
				else {
					cout << demand[i][t];
				}

			}
			if (i < N - 1) {
				cout << "]," << endl;
			}
			else {
				cout << "]";
			}

		}
		cout << "]" << endl << endl << endl;



		//Printing Transportation Cost
		cout << "Transportation cost between nodes: " << endl;
		cout << "[";
		for (int i = 0; i < N; ++i) {
			cout << "[";
			for (int t = 0; t < N; ++t) {
				if (t < N - 1) {
					cout << tranport_cost[i][t] << ", ";
				}
				else {
					cout << tranport_cost[i][t];
				}

			}
			if (i < N - 1) {
				cout << "]," << endl;
			}
			else {
				cout << "]";
			}

		}
		cout << "]" << endl << endl << endl;



		//Printing Inventory Holding Cost
		cout << "Inventory Holding Cost: " << endl;
		cout << "[";
		for (int i = 0; i < N; ++i) {
			if (i < N - 1) {
				cout << holding_cost[i] << ", ";
			}
			else {
				cout << holding_cost[i];
			}

		}
		cout << "]" << endl << endl << endl;



		//Printing Penalty for Unmet Demand
		cout << "Penalty for Unmet Demand: " << endl;
		cout << "[";
		for (int i = 0; i < N; ++i) {
			if (i < N - 1) {
				cout << penalty[i] << ", ";
			}
			else {
				cout << penalty[i];
			}

		}
		cout << "]" << endl << endl << endl;



		//Printing Initial Inventory
		cout << "Initial Inventory: " << endl;
		cout << "[";
		for (int i = 0; i < N; ++i) {
			if (i < N - 1) {
				cout << init_inventory[i] << ", ";
			}
			else {
				cout << init_inventory[i];
			}

		}
		cout << "]" << endl << endl << endl;



		//Printing Inventory Capacity
		cout << "Inventory Capacity: " << endl;
		cout << "[";
		for (int i = 0; i < N; ++i) {
			if (i < N - 1) {
				cout << inventory_cap[i] << ", ";
			}
			else {
				cout << inventory_cap[i];
			}

		}
		cout << "]" << endl << endl << endl;

#pragma endregion


#pragma region Defining Decision Variables for Master Problem

		/////DECISION VARIABLES for MASTER PROBLEM/////
		IloNumVarArray Y(env, T, 0, IloInfinity, ILOBOOL);  // 1 if production takesplace in time period T, OW = 0.
		IloNumArray	Y_val(env, T);  // To store "Y" values.


		// Z[i][k][t] is equal to 1 iff node i is visited by vehicle k in period t, OW = 0.
		VarArray3D Z(env, T);  // Defining 3D array of decision variables corresponding to each time period.
		Array3D Z_val(env, T);  // To store Z values in a 3D Array.

		for (int i = 0; i < T; i++) {
			Z[i] = VarArray2D(env, K);  // Decision variables corresponding to each vehicle.
			Z_val[i] = Array2D(env, K);  // To store Z[i] values.

			for (int j=0; j < K; j++) {
				Z[i][j] = IloNumVarArray(env, N, 0, 1, ILOBOOL);  // Decision variables corresponding to each locations.
				Z_val[i][j] = IloNumArray(env, N);
			}
		}


		// X[i][j][k][t] represents the number of times vehicle k travels directly between node i and node j in period t.
		// Xijkt = {0, 1} and X0jkt = {0, 1, 2}.
		VarArray4D X(env, T);  // Defining 4D array of decision variables corresponding to each time period.
		Array4D X_val(env, T);  // To store X values in a 4D Array.

		for (int i = 0; i < T; i++) {
			X[i] = VarArray3D(env, K);  // Decision variables corresponding to each vehicle.
			X_val[i] = Array3D(env, K);  // To store X[i] values.

			for (int j=0; j < K; j++) {
				X[i][j] = VarArray2D(env, N);  // Decision variables corresponding to each locations.
				X_val[i][j] = Array2D(env, N);  // To store X[i][j] values.

				for (int k=0; k < N; k++) {
					// k==0 is plant, which can have value 0, 1, and 2. All other locations can have only values 0 and 1.
					if (k == 0) {
						X[i][j][k] = IloNumVarArray(env, N, 0, 2, ILOINT);  // Decision variables corresponding to each locations.
						X_val[i][j][k] = IloNumArray(env, N);  // To store X[i][j][k] values.
					}
					else {
						X[i][j][k] = IloNumVarArray(env, N, 0, 1, ILOBOOL);  // Decision variables corresponding to each locations.
						X_val[i][j][k] = IloNumArray(env, N);  // To store X[i][j][k] values.
					}
				}
			}
		}


		IloNumVar eta(env, 0, IloInfinity, ILOFLOAT);  // DSP gives the value of Theta.
		IloNum eta_val = 0;  // To store "Theta" values.
#pragma endregion


#pragma region Defining Decision Variables for Dual Sub Problem
		/////DUAL DECISION VARIABLES for DUAL PROBLEM/////
		
		//Alpha[t] - Defined for each 't' in T.
		IloNumVarArray Alpha(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray Alpha_val(env, T);

		
		//Beta[t][i] - Defined for each 'i' in Nc and each 't' in T.
		VarArray2D Beta(env, T); 
		Array2D Beta_val(env, T);  

		for (int i=0; i < T; i++) {
			Beta[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);  // For values from 1 to N, ie only customer locations.
			Beta_val[i] = IloNumArray(env, N);
		}

		
		//Gama[t] - Defined for 't' in T.
		IloNumVarArray Gamma(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray Gamma_val(env, T);


		//Theta - Defined for each 'i' in Nc and each 't' in T.
		VarArray2D Theta(env, T);
		Array2D Theta_val(env, T);

		for (int i = 0; i < T; i++) {
			Theta[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);  // For values from 1 to N, ie only customer locations.
			Theta_val[i] = IloNumArray(env, N);
		}

		//Delta - Defined for 't' in T.
		IloNumVarArray Delta(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray Delta_val(env, T);


		//Kappa - Defined for each 'k' in K and each 't' in T.
		VarArray2D Kappa(env, T);
		Array2D Kappa_val(env, T);

		for (int i = 0; i < T; i++) {
			Kappa[i] = IloNumVarArray(env, K, 0, IloInfinity, ILOFLOAT);  
			Kappa_val[i] = IloNumArray(env, K);
		}


		//Zeta - Defined for each 'i' in Nc, for each 'k' in K and each 't' in T.
		VarArray3D Zeta(env, T);  // Defining 3D array of decision variables corresponding to each time period.
		Array3D Zeta_val(env, T);  // To store Zeta values in a 3D Array.

		for (int i = 0; i < T; i++) {
			Zeta[i] = VarArray2D(env, K);  // Decision variables corresponding to each vehicle.
			Zeta_val[i] = Array2D(env, K);  // To store Z[i] values.

			for (int j = 0; j < K; j++) {
				Zeta[i][j] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);  // Decision variables corresponding to each customer locations, hence N = 1 to n.
				Zeta_val[i][j] = IloNumArray(env, N);
			}
		}

#pragma endregion



		//////////Part 2 - DEVELOP GENERIC MODEL//////////

#pragma region Defining Master Problem
		/////SET MASTER PROBLEM/////
		// Objective Function: Minimize: sum {t in 1..T} ((setup cost * Y[t]) + (sum {(i, j) in E}{k in 1..K}transport cost[i][j]*X[i][j][k][t])) + eta
		IloModel model_master(env);

		IloExpr obj_master(env);

		for (int t = 0; t < T; t++) {
			obj_master += f * Y[t];  // Total setup cost.

			for (int k = 0; k < K; k++) {
				for (const auto& edge : E) {
					int i = edge.first;
					int j = edge.second;
					obj_master += tranport_cost[i][j] * X[t][k][i][j];  // Total transportation cost.
				}
			}
		}
		obj_master += eta;  // Adding DSP objective value, ie total flow cost.

		model_master.add(IloMinimize(env, obj_master));
		obj_master.end();
		IloCplex cplex_master(env);
		cplex_master.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_master.setWarning(env.getNullStream()); // This is to supress warning messages on screen
#pragma endregion


#pragma region Defining Dual Sub Problem
		/////SET SUBPROBLEM (DUAL FORMULATION)/////
		IloModel model_sub(env);
		IloObjective Objective_sub = IloMaximize(env);
		model_sub.add(Objective_sub);


		// Refer Inequality 43 - Supplementary material
		//Constraint 1 - for {t in 1..T}: Alpha[t] - Delta[t] <= u;
		for (int t = 0; t < T; ++t) {
			model_sub.add(Alpha[t] - Delta[t] <= u);
		}


		// Refer Inequality 44 - Supplementary material
		//Constraint 2 - for {t in 1..T}: -Alpha[t] + Alpha[t+1] - Gamma[t+1] <= holding cost of plant;
		// Loop till T-1, because t+1 at the end will be undefined - To check why they have taken like this?
		for (int t = 0; t < T - 1; t++) {
			model_sub.add(-Alpha[t] + Alpha[t + 1] - Gamma[t + 1] <= holding_cost[0]);
		}


		// Refer Inequality 45 - Supplementary material
		//Constraint 3 - for {i in Nc} for {t in T}: -Beta[i][t] + Beta[i][t+1] - Theta[i][t] <= holding cost[i].
		// Loop till T-1, because t+1 at the end will be undefined - To check why they have taken like this?
		for (int t = 0; t < T - 1; t++) {
			for (int i = 1; i < N; i++) {
				model_sub.add(-Beta[t][i] + Beta[t + 1][i] - Theta[t][1] <= holding_cost[i]);
			}
		}


		// Refer Inequality 46 - Supplementary material
		//Constraint 4 - for {i in Nc} for {k in K} for {t in T}: -Alpha[t] + Beta[i][t] - Kappa[k][t] - Zeta[i][k][t] <= 0.
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (int i = 1; i < N; i++) {
					model_sub.add(-Alpha[t] + Beta[t][i] - Kappa[t][k] - Zeta[t][k][i] <= 0);
				}
			}
		}



		// Refer Inequality 47 - Supplementary material
		//Constraint 5 - for {i in Nc} for {t in T}: Beta[i][t] <= penalty[i].
		for (int t = 0; t < T; t++) {
			for (int i = 1; i < N; i++) {
				model_sub.add(Beta[t][i] <= penalty[i]);
			}
		}


		IloCplex cplex_sub(model_sub);
		cplex_sub.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_sub.setWarning(env.getNullStream()); //This is to supress warning messages on screen
#pragma endregion



		//////////Part 3 - BEGIN ITERATIONS//////////

		IloNum eps = cplex_sub.getParam(IloCplex::EpInt);//Integer tolerance for MIP models;
		//default value of EpInt remains 1e-5 http://www.iro.umontreal.ca/~gendron/IFT6551/CPLEX/HTML/relnotescplex/relnotescplex12.html




		/*
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

		*/
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