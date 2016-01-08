# IGCSimulation
Gene Conversion Simulation project

This simulation study aims at discovering the influence of inter-locus gene conversion (IGC) tract length on 
the statistical analysis which assumed independence among site. Part of the simulator was finished as a course 
project of CSC530. 

In the simulation, IGC follows an poisson process indenpendent from point mutation process with two new parameters: IGC initiation rate per site (rIGC), Average tract length (lIGC). It is assumed that IGC tract length follows a geometric distribution. 

The two new parameters relates with the proposed IGC parameter tau by rIGC x lIGC = Tau. 

In this simulation study, estimated parameter values from Yeast paralog pair YDR418W_YEL054C are chosen. lIGC of value 1, 10, 50, 100, 500 is used to generate 100 simulated dataset each time to study the influence of tract length on the standard deviation of parameter estimations from the IGC expansion model. 