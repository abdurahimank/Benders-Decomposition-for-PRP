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

#pragma region Defining High Dimensional Arrays
typedef IloArray<IloNumArray> Array2D;  // Defining a normal 2D array.
typedef IloArray<Array2D> Array3D;  // Defining a normal 3D array.
typedef IloArray<Array3D> Array4D;  // Defining a normal 4D array.


typedef IloArray<IloNumVarArray> VarArray2D;  // Defining a 2D array of decision variables
typedef IloArray<VarArray2D> VarArray3D;  // Defining a 3D array of decision variables
typedef IloArray<VarArray3D> VarArray4D;  // Defining a 4D array of decision variables
#pragma endregion


ILOSTLBEGIN

int main(int argc, char** argv)
{
	IloEnv env;
	try
	{
		//////////Part 1 - DEFINE VARIABLES AND PARAMETERS//////////

#pragma region Defining Constants and Sets
		/////CONSTANTS/////
		int N = 5;  // Number of customers + Plant (# of Plant = 1).
		int T = 3;  // Number of time periods.
		int K = 1;  // Number of vehicles available.

		double C = 100;  // Production capacity of the plant in a time period.
		double Q = 50;  // Vehicle capacity (all vehicles are identical).

		double u = 10;  // Unit cost of production.
		double f = 500;  // Setup cost for production.

		///Defining set E.
		vector<pair<int, int>> E;  // Set of edges, such that i < j.
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				E.emplace_back(i, j);  // store edge (i,j) only when i < j
			}
		}

		/*
		//Print set E.
		cout << "Edges (i, j) where i < j:" << endl;
		for (const auto& edge : E) {
			cout << "(" << edge.first << ", " << edge.second << ")" << endl;
		}
		*/


		///Defining set D - 
		vector<vector<pair<int, int>>> D_i(N);  // delta[i] = list of (i,j) or (j,i) edges

		for (auto& e : E) {
			int u = e.first, v = e.second;
			D_i[u].emplace_back(u, v);
			D_i[v].emplace_back(u, v);  // because undirected
		}

		/*
		//Print set Del.
		for (int i = 0; i < N; i++) {
			std::cout << "D[" << i << "] = { ";
			for (const auto& edge : D_i[i]) {
				std::cout << "(" << edge.first << "," << edge.second << ") ";
			}
			std::cout << "}" << std::endl;
		}
		*/




		


#pragma endregion


#pragma region Defining Parameters and Reading Values


		/////PARAMETERS/////
		Array2D demand(env);  // Demand at a customer node in a time period.

		Array2D tranport_cost(env);  // Transportation cost between two nodes.
		IloNumArray holding_cost(env);  // Inventory holding cost at each node.
		IloNumArray penalty(env);  // Penalty at node i, if demand is unmet in a period.

		//Array2D init_inventory(env);  // Initial inventory at each node and at each time period.
		IloNumArray init_inventory(env);  // Initial inventory at each node and at each time period.
		IloNumArray inventory_cap(env);  // Inventory capacity at each node (node 0 is the plant).


		///READING PARAMETER DATA///
		const char* data_filename = "Data_T3_K1_N5.dat";  // Data with 3 time periods, 1 vehicle and 10 locations including plant.
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


		IloBool consistentData = (demand[0].getSize() == tranport_cost.getSize());  // Ensuring nodes are same in different matrices.
		if (!consistentData)
		{
			cerr << "ERROR: data file '" << data_filename << "' contains inconsistent data" << endl;
			throw(-1);
		}
		datafile.close();

		// READING DONE......

		//Defining set SetM[t]
		vector<double> SetM(T);  // M in each T time periods.

		for (int t = 0; t < T; t++) {
			double sum_d = 0.0;
			for (int j = t; j < T; j++) {
				for (int i = 1; i < N; i++) {
					sum_d += demand[j][i];
				}
			}
			SetM[t] = min(C, sum_d);
		}

		//Defining set SetN[t][i]
		vector<vector<double>> SetN(T, vector<double>(N));

		for (int t = 0; t < T; t++) {
			for (int i = 1; i < N; i++) {
				double sum_d = 0.0;
				for (int j = t; j < T; j++) {
					sum_d += demand[j][i];
				}
				SetN[t][i] = min({ inventory_cap[i], Q, sum_d });
			}
		}



		///DISPLAYING PARAMETER DATA///
		//Printing Demand
		cout << "Demand at each node in each time period: " << endl;
		cout << "[";
		for (int t = 0; t < T; ++t) {
			cout << "[";
			for (int i = 0; i < N; ++i) {
				if (i < N - 1) {
					cout << demand[t][i] << ", ";
				}
				else {
					cout << demand[t][i];
				}

			}
			if (t < T - 1) {
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

			for (int j = 0; j < K; j++) {
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

			for (int j = 0; j < K; j++) {
				X[i][j] = VarArray2D(env, N);  // Decision variables corresponding to each locations.
				X_val[i][j] = Array2D(env, N);  // To store X[i][j] values.

				for (int k = 0; k < N; k++) {
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

		/*

		for (int t = 0; t < T; t++) {
			X[t] = VarArray3D(env, K);  // Decision variables corresponding to each vehicle.
			X_val[t] = Array3D(env, K);  // To store X[i] values.

			for (int k = 0; k < K; k++) {
				X[t][k] = VarArray2D(env, N);  // Decision variables corresponding to each locations.
				X_val[t][k] = Array2D(env, N);
				for (const auto& edge : E) {
					int i = edge.first;
					int j = edge.second;
					if (i == 0) {
						X[t][k][i] = IloNumVarArray(env, j, 0, 2, ILOINT);  // Decision variables corresponding to each locations.
						X_val[t][k][i] = IloNumArray(env, j);  // To store X[i][j][k] values.
					}
					else {
						X[t][k][i] = IloNumVarArray(env, j, 0, 1, ILOBOOL);  // Decision variables corresponding to each locations.
						X_val[t][k][i] = IloNumArray(env, N);  // To store X[i][j][k] values.
					}
				}
			}
		}
		*/


		IloNumVar eta(env, 0, IloInfinity, ILOFLOAT);  // DSP gives the value of Theta.
		IloNum eta_val = 0;  // To store "Theta" values.
#pragma endregion


#pragma region Defining Decision Variables for Dual Sub Problem
		/////DUAL DECISION VARIABLES for DUAL SUB PROBLEM/////

		//Alpha[t] - Defined for each 't' in T.
		IloNumVarArray Alpha(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray Alpha_val(env, T);


		//typedef IloArray<IloNumVarArray> VarArray2D;  // Defining a 2D array of decision variables
		//typedef IloArray<IloNumArray> Array2D;  // Defining a normal 2D array.
		//Beta[t][i] - Defined for each 'i' in Nc and each 't' in T.
		VarArray2D Beta(env, T);
		Array2D Beta_val(env, T);

		for (int i = 0; i < T; i++) {
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


#pragma region Defining Decision Variables for Extreme Ray Problem
		/////DUAL DECISION VARIABLES for EXTREME RAY PROBLEM/////

		//Alpha[t] - Defined for each 't' in T.
		IloNumVarArray AlphaEr(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray AlphaEr_val(env, T);


		//Beta[t][i] - Defined for each 'i' in Nc and each 't' in T.
		VarArray2D BetaEr(env, T);
		Array2D BetaEr_val(env, T);

		for (int i = 0; i < T; i++) {
			BetaEr[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);  // For values from 1 to N, ie only customer locations.
			BetaEr_val[i] = IloNumArray(env, N);
		}


		//Gama[t] - Defined for 't' in T.
		IloNumVarArray GammaEr(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray GammaEr_val(env, T);


		//Theta - Defined for each 'i' in Nc and each 't' in T.
		VarArray2D ThetaEr(env, T);
		Array2D ThetaEr_val(env, T);

		for (int i = 0; i < T; i++) {
			ThetaEr[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);  // For values from 1 to N, ie only customer locations.
			ThetaEr_val[i] = IloNumArray(env, N);
		}

		//Delta - Defined for 't' in T.
		IloNumVarArray DeltaEr(env, T, 0, IloInfinity, ILOFLOAT);
		IloNumArray DeltaEr_val(env, T);


		//Kappa - Defined for each 'k' in K and each 't' in T.
		VarArray2D KappaEr(env, T);
		Array2D KappaEr_val(env, T);

		for (int i = 0; i < T; i++) {
			KappaEr[i] = IloNumVarArray(env, K, 0, IloInfinity, ILOFLOAT);
			KappaEr_val[i] = IloNumArray(env, K);
		}


		//Zeta - Defined for each 'i' in Nc, for each 'k' in K and each 't' in T.
		VarArray3D ZetaEr(env, T);  // Defining 3D array of decision variables corresponding to each time period.
		Array3D ZetaEr_val(env, T);  // To store Zeta values in a 3D Array.

		for (int i = 0; i < T; i++) {
			ZetaEr[i] = VarArray2D(env, K);  // Decision variables corresponding to each vehicle.
			ZetaEr_val[i] = Array2D(env, K);  // To store Z[i] values.

			for (int j = 0; j < K; j++) {
				ZetaEr[i][j] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);  // Decision variables corresponding to each customer locations, hence N = 1 to n.
				ZetaEr_val[i][j] = IloNumArray(env, N);
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

		///ADDING CONSTRAINTS///
		//Refer constrain 9 in paper
		/// for {i in Nc} for {t in T} sum {k in K} (Z[t][k][i] <= 1)
		for (int i = 1; i < N; i++) {
			for (int t = 0; t < T; t++) {
				IloExpr sum_Ztk(env);
				for (int k = 0; k < K; k++) {
					sum_Ztk += Z[t][k][i];
				}
				model_master.add(sum_Ztk <= 1);
				sum_Ztk.end();
			}
		}


		//Refer constrain 10 in paper - If a node is visited, there should be incoming and outgoing.
		/// for {t in T} for {k in K} for {i in N} sum {(j,j') in D_i} X[t][k][j][j'] = 2*Z[t][k][i]
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (int i = 0; i < N; i++) {
					IloExpr exp10(env);

					for (int j = 0; j < N; j++) {
						if (j == i) continue;

						if (i < j) {
							exp10 += X[t][k][i][j];  // edge (i,j)
						}
						else {
							exp10 += X[t][k][j][i];  // edge (j,i)
						}
					}

					model_master.add(exp10 == 2 * Z[t][k][i]);
					exp10.end();  // cleanup
				}
			}
		}


		//Refer constrain 11 in paper - Subtour elimination





		//Refer constrain 13 in paper
		/// Y[t] => {0, 1} and Z[t][k][i] => {0, 1}
		// Implemented when decision variable is defined.


		//Refer constrain 14 in paper
		/// X[t][k][i][j] where i != 0 => {0, 1}
		// Implemented when decision variable is defined.


		//Refer constrain 15 in paper
		/// X[t][k][0][j] => {0, 1, 2}
		// Implemented when decision variable is defined.


		//Refer constrain 16 in paper
		/// for {t in T} for {k in K} for {i in Nc} (Z[t][k][i] <= Z[t][k][0])
		for (int t = 0; t < T; t++) {           
			for (int k = 0; k < K; k++) {       
				for (int i = 1; i < N; i++) {  // i in N_c (excluding depot i = 0)
					model_master.add(Z[t][k][i] <= Z[t][k][0]);
				}
			}
		}


		//Refer constrain 17 in paper
		/// for {t in T} for {k in K} for {(i, j) in E(Nc)} X[t][k][i][j] <= Z[t][k][i] and X[t][k][i][j] <= Z[t][k][j]
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (auto& edge : E) {
					int i = edge.first;
					int j = edge.second;

					// Skip edges connected to depot (node 0)
					if (i == 0 || j == 0) continue;

					model_master.add(X[t][k][i][j] <= Z[t][k][i]);
					model_master.add(X[t][k][i][j] <= Z[t][k][j]);
				}
			}
		}


		//Refer constrain 18 in paper - Symmetry breaking constraint for Multi vehicle case


		//Refer constrain 19 in paper - Symmetry breaking constraint for Multi vehicle case



		IloCplex cplex_master(env);
		cplex_master.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_master.setWarning(env.getNullStream()); // This is to supress warning messages on screen
#pragma endregion


#pragma region Defining Dual Sub Problem
		/////DEFINING DUAL SUBPROBLEM/////
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
				model_sub.add(-Beta[t][i] + Beta[t + 1][i] - Theta[t][i] <= holding_cost[i]);
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


#pragma region Defining Extreme Ray Problem
		/////EXTREME RAY PROBLEM (DUAL FORMULATION)/////
		IloModel model_sub_er(env);
		IloObjective Objective_sub_er = IloMaximize(env);
		model_sub_er.add(Objective_sub_er);


		// Refer Inequality 43 - Supplementary material
		///Constraint 1 - for {t in 1..T}: Alpha[t] - Delta[t] <= u;
		for (int t = 0; t < T; t++) {
			model_sub_er.add(AlphaEr[t] - DeltaEr[t] <= 0);
		}


		// Refer Inequality 44 - Supplementary material
		///Constraint 2 - for {t in 1..T}: -Alpha[t] + Alpha[t+1] - Gamma[t+1] <= holding cost of plant;
		// Loop till T-1, because t+1 at the end will be undefined - To check why they have taken like this?
		for (int t = 0; t < T - 1; t++) {
			model_sub_er.add(-AlphaEr[t] + AlphaEr[t + 1] - GammaEr[t + 1] <= 0);
		}


		// Refer Inequality 45 - Supplementary material
		///Constraint 3 - for {i in Nc} for {t in T}: -Beta[i][t] + Beta[i][t+1] - Theta[i][t] <= holding cost[i].
		// Loop till T-1, because t+1 at the end will be undefined - To check why they have taken like this?
		for (int t = 0; t < T - 1; t++) {
			for (int i = 1; i < N; i++) {
				model_sub_er.add(-BetaEr[t][i] + BetaEr[t + 1][i] - ThetaEr[t][1] <= 0);
			}
		}


		// Refer Inequality 46 - Supplementary material
		///Constraint 4 - for {i in Nc} for {k in K} for {t in T}: -Alpha[t] + Beta[i][t] - Kappa[k][t] - Zeta[i][k][t] <= 0.
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (int i = 1; i < N; i++) {
					model_sub_er.add(-AlphaEr[t] + BetaEr[t][i] - KappaEr[t][k] - ZetaEr[t][k][i] <= 0);
				}
			}
		}


		// Refer Inequality 47 - Supplementary material
		///Constraint 5 - for {i in Nc} for {t in T}: Beta[i][t] <= penalty[i].
		for (int t = 0; t < T; t++) {
			for (int i = 1; i < N; i++) {
				model_sub_er.add(BetaEr[t][i] <= 0);
			}
		}


		///Constraint 6 - Sum of all duals equal to 1
		IloExpr sumDuals(env);

		// Add AlphaEr[t] 
		for (int t = 0; t < T; t++) {
			sumDuals += AlphaEr[t];
		}

		// Add BetaEr[t][i]
		for (int t = 0; t < T; t++) {
			for (int i = 1; i < N; i++) {
				sumDuals += BetaEr[t][i];
			}
		}

		// Add GammaEr[t]
		for (int t = 1; t < T; t++) {  // starts at 1 if GammaEr[0] isn't used
			sumDuals += GammaEr[t];
		}

		// Add ThetaEr[t][i]
		for (int t = 0; t < T - 1; t++) {
			for (int i = 1; i < N; i++) {
				sumDuals += ThetaEr[t][i];
			}
		}

		// Add DeltaEr[t]
		for (int t = 0; t < T; t++) {
			sumDuals += DeltaEr[t];
		}

		// Add KappaEr[t][k]
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				sumDuals += KappaEr[t][k];
			}
		}

		// Add ZetaEr[t][k][i]
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (int i = 1; i < N; i++) {
					sumDuals += ZetaEr[t][k][i];
				}
			}
		}

		// Add the normalization constraint
		model_sub_er.add(sumDuals == 1);
		sumDuals.end();  // To free the memory



		IloCplex cplex_sub_er(model_sub_er);
		cplex_sub_er.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_sub_er.setWarning(env.getNullStream()); //This is to supress warning messages on screen
#pragma endregion



		//////////Part 3 - ITERATIONS//////////

#pragma region Initiliazing Values
		/////INITILIAZING VALUES TO START ITERATION/////

		IloNum GAP = IloInfinity;  // Gap between LB and UB initially set to inifinity.
		IloNum eps = cplex_sub.getParam(IloCplex::EpInt);//Integer tolerance for MIP models;
		//default value of EpInt remains 1e-5 http://www.iro.umontreal.ca/~gendron/IFT6551/CPLEX/HTML/relnotescplex/relnotescplex12.html

		eta_val = 0;  // eta_val set to 0.

		// X[t][k][i][j] values set to 0.
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						X_val[t][k][i][j] = 0;
					}
				}
			}
		}

		// Z[t][k][i] values set to 0.
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				for (int i = 0; i < N; i++) {
					Z_val[t][k][i] = 0;
				}
			}
		}

		// Y[t] values set to 0.
		for (int t = 0; t < T; t++) {
			Y_val[t] = 0;
		}

		IloNum sub_obj_val = 0;
		IloNum Upper_bound = IloInfinity;
		IloNum Lower_bound = 0;
		// GAP = Upper_bound - Lower_bound;
		IloInt Iter = 0;
#pragma endregion


		/////ITERATION STARTING/////
		//while (Iter < 10)
		//while (GAP > eps)
		while (Upper_bound - Lower_bound > eps)
		{
			Iter++;
			cout << "=========================================" << endl;
			cout << "============ITERATION " << Iter << "==============" << endl;

#pragma region Solving Dual Sub Problem
			/////SOLVING SUB PROBLEM/////

			// Refer Equation 32 in paper.
			//Define Objective Function for the Dual Sub Problem.
			IloExpr sub_obj(env);
			// Term 1: -Initial Inventory[0][0] * alpha[1]
			sub_obj += -init_inventory[0] * Alpha[0];

			// Term 2: sum {i in Nc} (demand[t][i] - Initial Inventory[i][0]) * Beta[t][i]
			for (int i = 1; i < N; i++) {
				sub_obj += (demand[1][i] - init_inventory[i]) * Beta[1][i];
			}

			// Term 3: sum {i in Nc} sum {t = 2..T} (demand[t][i] * beta[t][i])
			for (int t = 1; t < T; t++) {
				for (int i = 1; i < N; i++) {
					sub_obj += demand[t][i] * Beta[t][i];
				}
			}

			// Term 4: sum {t in T} (Inventory Capacity[0] * gamma[t])
			for (int t = 0; t < T; t++) {
				sub_obj += -inventory_cap[0] * Gamma[t];
			}

			// Term 5: sum {t in T} sum {i in Nc) (Inventory Capacity[i] - demand[t][i]) * theta[t][i]
			for (int t = 0; t < T; t++) {
				for (int i = 1; i < N; i++) {
					sub_obj += -(inventory_cap[i] - demand[t][i]) * Theta[t][i];
				}
			}

			// Term 6: sum {t in T} (-SetM[t] * Y_val[t] * delta[t])
			for (int t = 0; t < T; t++) {
				sub_obj += -SetM[t] * Y_val[t] * Delta[t];
			}

			// Term 7: sum {t in T} sum {k in K} (-Q * Z[t][k][0] * kappa[t][k])
			for (int t = 0; t < T; t++) {
				for (int k = 0; k < K; k++) {
					sub_obj += -Q * Z_val[t][k][0] * Kappa[t][k];
				}
			}

			// Term 8: sum {t in T} sum {k in K} sum {i in Nc} (SetN[t][i] * Z[t][k][i] * zeta[t][k][i])
			for (int t = 0; t < T; t++) {
				for (int k = 0; k < K; k++) {
					for (int i = 1; i < N; i++) {
						sub_obj += -SetN[t][i] * Z_val[t][k][i] * Zeta[t][k][i];
					}
				}
			}

			Objective_sub.setExpr(IloMaximize(env, sub_obj));
			cout << "SOLVING SUB PROBLEM" << endl;

			cplex_sub.setParam(cplex_sub.PreInd, 0);   //Disable presolve, otherwise, if dual is infeasible,
			//we don't know if prime is unbounded or infeasible
			cplex_sub.setParam(IloCplex::RootAlg, IloCplex::Primal);//Solve the SP Dual using Primal Simplex


			cplex_sub.solve();
			cout << "Sub Problem Solution Status: " << cplex_sub.getCplexStatus() << endl;

			sub_obj_val = cplex_sub.getObjValue();
			cout << "sub_obj_val = " << sub_obj_val << endl << endl;
#pragma endregion


#pragma region Adding Benders Optimality Cut
			
			/////ADDING BENDERS OPTIMALITY CUT/////
			if (cplex_sub.getCplexStatus() == CPX_STAT_OPTIMAL || cplex_sub.getCplexStatus() == CPX_STAT_UNBOUNDED)
			{// Add Optimality Cut to the Master Problem for both optimal an unbounbed subproblem.
				
				// Storing values of dual variables to respective arrays
				cplex_sub.getValues(Alpha_val, Alpha);  
				
				// Since getValues() can be used for only 1D array.
				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						Beta_val[t][i] = cplex_sub.getValue(Beta[t][i]);
					}
				}
				
				cplex_sub.getValues(Gamma_val, Gamma);
				
				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						Theta_val[t][i] = cplex_sub.getValue(Theta[t][i]);
					}
				}
				
				cplex_sub.getValues(Delta_val, Delta);
				
				for (int t = 0; t < T; t++) {
					cplex_sub.getValues(Kappa_val[t], Kappa[t]);
				}
				
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						for (int i = 1; i < N; i++) {
							Zeta_val[t][k][i] = cplex_sub.getValue(Zeta[t][k][i]);
						}
					}
				}
				
				
				/// CALCULATING UPPER BOUND
				double new_UB = 0.0;
				for (int t = 0; t < T; t++) {
					new_UB += f * Y_val[t];  // Total setup cost.

					for (int k = 0; k < K; k++) {
						for (const auto& edge : E) {
							int i = edge.first;
							int j = edge.second;
							new_UB += tranport_cost[i][j] * X_val[t][k][i][j];  // Total transportation cost.
						}
					}
				}
				new_UB += sub_obj_val;

				Upper_bound = IloMin(Upper_bound, new_UB);
				cout << "UPPER BOUND = " << Upper_bound << endl << endl;

				
				///ADDING OPTIMALITY CUT TO THE MASTER PROBLEM
				IloExpr lhs(env);
				IloNum rhs_const = 0.0;

				// ==== LHS ====
				// Add eta variable
				lhs += eta;

				// Add decision variable terms: Y[t]
				for (int t = 0; t < T; t++) {
					lhs += SetM[t] * Delta_val[t] * Y[t];
				}

				// Add decision variable terms: Z[t][k][0]
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						lhs += Q * Kappa_val[t][k] * Z[t][k][0];
					}
				}

				// Add decision variable terms: Z[t][k][i]
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						for (int i = 1; i < N; i++) {
							lhs += SetN[t][i] * Zeta_val[t][k][i] * Z[t][k][i];
						}
					}
				}

				// ==== RHS ====
				// All terms without master decision variables (purely numerical values)

				// Term 1: -Initial Inventory[0][0] * alpha[1]
				rhs_const += -init_inventory[0] * Alpha_val[0];

				// Term 2: sum {i in Nc} (demand[t][i] - Initial Inventory[i][0]) * Beta[t][i]
				for (int i = 1; i < N; i++) {
					rhs_const += (demand[1][i] - init_inventory[i])* Beta_val[1][i];
				}

				// Term 3: sum {i in Nc} sum {t = 2..T} (demand[t][i] * beta[t][i])
				for (int t = 1; t < T; t++) {
					for (int i = 1; i < N; i++) {
						rhs_const += demand[t][i] * Beta_val[t][i];
					}
				}

				// Term 4: sum {t in T} (Inventory Capacity[0] * gamma[t])
				for (int t = 0; t < T; t++) {
					rhs_const += -inventory_cap[0] * Gamma_val[t];
				}

				// Term 5: sum {t in T} sum {i in Nc) (Inventory Capacity[i] - demand[t][i]) * theta[t][i]
				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						rhs_const += -(inventory_cap[i] - demand[t][i]) * Theta_val[t][i];
					}
				}


				// Add the cut to master problem
				model_master.add(lhs >= rhs_const);
				cout << "Optimality Cut Added to Master Problem" << endl;
				// Cleanup
				lhs.end();
				
			}
			
#pragma endregion


#pragma region Solving Extreme Ray Problem
			
			/////SOLVING EXTREME RAY PROBLEM/////
			if (cplex_sub.getCplexStatus() == CPX_STAT_UNBOUNDED)
			{
				cout << endl << endl << "SOLVING EXTREME RAY PROBLEM" << endl;

				//Define Objective Function for the Extreme Ray Problem.
				IloExpr sub_obj_er(env);
				// Term 1: -Initial Inventory[0][0] * alpha[1]
				sub_obj_er += -init_inventory[0] * AlphaEr[0];

				// Term 2: sum {i in Nc} (demand[t][i] - Initial Inventory[i][0]) * Beta[t][i]
				for (int i = 1; i < N; i++) {
					sub_obj_er += (demand[1][i] - init_inventory[i]) * BetaEr[1][i];
				}

				// Term 3: sum {i in Nc} sum {t = 2..T} (demand[t][i] * beta[t][i])
				for (int t = 1; t < T; t++) {
					for (int i = 1; i < N; i++) {
						sub_obj_er += demand[t][i] * BetaEr[t][i];
					}
				}

				// Term 4: sum {t in T} (Inventory Capacity[0] * gamma[t])
				for (int t = 0; t < T; t++) {
					sub_obj_er += -inventory_cap[0] * GammaEr[t];
				}

				// Term 5: sum {t in T} sum {i in Nc) (Inventory Capacity[i] - demand[t][i]) * theta[t][i]
				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						sub_obj_er += -(inventory_cap[i] - demand[t][i]) * ThetaEr[t][i];
					}
				}

				// Term 6: sum {t in T} (-SetM[t] * Y_val[t] * delta[t])
				for (int t = 0; t < T; t++) {
					sub_obj_er += -SetM[t] * Y_val[t] * DeltaEr[t];
				}

				// Term 7: sum {t in T} sum {k in K} (-Q * Z[t][k][0] * kappa[t][k])
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						sub_obj_er += -Q * Z_val[t][k][0] * KappaEr[t][k];
					}
				}

				// Term 8: sum {t in T} sum {k in K} sum {i in Nc} (SetN[t][i] * Z[t][k][i] * zeta[t][k][i])
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						for (int i = 1; i < N; i++) {
							sub_obj_er += -SetN[t][i] * Z_val[t][k][i] * ZetaEr[t][k][i];
						}
					}
				}

				Objective_sub_er.setExpr(IloMaximize(env, sub_obj_er));
				
				cplex_sub_er.setParam(cplex_sub.PreInd, 0);   //Disable presolve, otherwise, if dual is infeasible,
				//we don't know if prime is unbounded or infeasible
				cplex_sub_er.setParam(IloCplex::RootAlg, IloCplex::Primal);//Solve the SP Dual using Primal Simplex


				cplex_sub_er.solve();
				cout << "Extreme Ray Problem Solution Status: " << cplex_sub_er.getCplexStatus() << endl;

				//cout << "SOLVING EXTREME RAY PROBLEM" << endl;
				// Dual subproblem is unbounded; Hence add feasibility Cut to the Master Problem
				//cplex_sub_er.getValues(X_dual_val_er, X_dual_er);  // taking values of X_dual from SP and saves to X_dual_val
				//cout << "X_dual of extreme rays = " << X_dual_val_er << endl;
				sub_obj_val = cplex_sub_er.getObjValue();
				cout << "Extreme ray obj val = " << sub_obj_val << endl;

			}
#pragma endregion


#pragma region Adding Benders Feasibility Cut
			
			/////ADDING BENDERS FEASIBILTY CUT/////
			if (cplex_sub.getCplexStatus() == CPX_STAT_UNBOUNDED)
			{// Dual subproblem is unbounded; Add Optimality Cut to the Master Problem
				
				// Storing values of dual variables to respective arrays
				cplex_sub_er.getValues(AlphaEr_val, AlphaEr);
				
				// Since getValues() can be used for only 1D array.
				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						BetaEr_val[t][i] = cplex_sub_er.getValue(BetaEr[t][i]);
					}
				}

				cplex_sub_er.getValues(GammaEr_val, GammaEr);

				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						ThetaEr_val[t][i] = cplex_sub_er.getValue(ThetaEr[t][i]);
					}
				}

				cplex_sub_er.getValues(DeltaEr_val, DeltaEr);

				for (int t = 0; t < T; t++) {
					cplex_sub_er.getValues(KappaEr_val[t], KappaEr[t]);
				}

				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						for (int i = 1; i < N; i++) {
							ZetaEr_val[t][k][i] = cplex_sub_er.getValue(ZetaEr[t][k][i]);
						}
					}
				}

				

				///ADDING FEASIBILITY CUT TO THE MASTER PROBLEM
				IloExpr Erlhs(env);
				IloNum Errhs_const = 0.0;

				// ==== LHS ====
				// Add eta variable
				//Erlhs += eta;

				// Add decision variable terms: Y[t]
				for (int t = 0; t < T; t++) {
					Erlhs += SetM[t] * DeltaEr_val[t] * Y[t];
				}

				// Add decision variable terms: Z[t][k][0]
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						Erlhs += Q * KappaEr_val[t][k] * Z[t][k][0];
					}
				}

				// Add decision variable terms: Z[t][k][i]
				for (int t = 0; t < T; t++) {
					for (int k = 0; k < K; k++) {
						for (int i = 1; i < N; i++) {
							Erlhs += SetN[t][i] * ZetaEr_val[t][k][i] * Z[t][k][i];
						}
					}
				}

				// ==== RHS ====
				// All terms without master decision variables (purely numerical values)

				// Term 1: -Initial Inventory[0][0] * alpha[1]
				Errhs_const += -init_inventory[0] * AlphaEr_val[0];

				// Term 2: sum {i in Nc} (demand[t][i] - Initial Inventory[i][0]) * Beta[t][i]
				for (int i = 1; i < N; i++) {
					Errhs_const += (demand[1][i] - init_inventory[i]) * BetaEr_val[1][i];
				}

				// Term 3: sum {i in Nc} sum {t = 2..T} (demand[t][i] * beta[t][i])
				for (int t = 1; t < T; t++) {
					for (int i = 1; i < N; i++) {
						Errhs_const += demand[t][i] * BetaEr_val[t][i];
					}
				}

				// Term 4: sum {t in T} (Inventory Capacity[0] * gamma[t])
				for (int t = 0; t < T; t++) {
					Errhs_const += -inventory_cap[0] * GammaEr_val[t];
				}

				// Term 5: sum {t in T} sum {i in Nc) (Inventory Capacity[i] - demand[t][i]) * theta[t][i]
				for (int t = 0; t < T; t++) {
					for (int i = 1; i < N; i++) {
						Errhs_const += -(inventory_cap[i] - demand[t][i]) * ThetaEr_val[t][i];
					}
				}

				cout << "Errhs const: " << Errhs_const << endl;
				// Add the cut to master problem
				//model_master.add(0 + Erlhs >= Errhs_const);
				cout << "Feasibility Cut Added to Master Problem" << endl;
				// Cleanup
				Erlhs.end();






				/*
				
					*/

			}
		
#pragma endregion


#pragma region Solving Master Problem
			
			/////SOLVING MASTER PROBLEM/////
			cout << endl << endl << "SOLVING MASTER PROBLEM" << endl;
			cout << "Master Problem Solution Status: " << cplex_master.getCplexStatus() << endl;
			cplex_master.extract(model_master);
			if (!cplex_master.solve())
			{
				cout << "Failed" << endl;
				throw(-1);
			}

			///STORING AND DISPLAYING DECISION VARIABLE VALUES OF MASTER PROBLEM
			// Y[t] values.
			cplex_master.getValues(Y_val, Y);
			cout << "Y = " << Y_val << endl;

			// X[t][k][i][j] values.
			for (int t = 0; t < T; t++) {
				for (int k = 0; k < K; k++) {
					for (const auto& edge : E) {
						int i = edge.first;
						int j = edge.second;
						X_val[t][k][i][j] = cplex_master.getValue(X[t][k][i][j]);  // Total transportation cost.
					}
				}
			}
			cout << "X Values = " << X_val << endl;

			// eta value.
			eta_val = cplex_master.getValue(eta);
			cout << "Eta = " << eta_val << endl;


			///CALCULATING LOWER BOUND
			double new_LB = 0.0;

			for (int t = 0; t < T; t++) {
				new_LB += f * Y_val[t];  // Total setup cost.

				for (int k = 0; k < K; k++) {
					for (const auto& edge : E) {
						int i = edge.first;
						int j = edge.second;
						
						new_LB += tranport_cost[i][j] * X_val[t][k][i][j];  // Total transportation cost.
					}
				}
			}
			new_LB += eta_val;  // Adding DSP objective value, ie total flow cost.

			cout << endl << "LOWER BOUND = " << new_LB << endl;

#pragma endregion


			//GAP = 0;  // Uncomment this to stop iteration at just one iteration.
			cout << endl << endl << endl;
		}  //while(Upper_bound - Lower_bound > eps)
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