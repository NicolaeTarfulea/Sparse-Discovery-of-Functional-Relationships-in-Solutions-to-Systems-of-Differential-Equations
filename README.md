# Sparse-Discovery-of-Functional-Relationships-in-Solutions-to-Systems-of-Differential-Equations
This repository accompanies the work:  *Sparse Discovery of Functional Relationships in Solutions to Systems of Differential Equations*. 
It provides a MATLAB implementation for discovering functional relationships among solution components of 
first-order systems of ordinary differential equations. The method leverages sparse identification 
techniques applied to numerical solutions of initial-value problems. The core assumption is that only a 
small number of terms govern the relationships between solution components, allowing the underlying structure 
to be identified efficiently.

**Features:**
* Sparse identification of relationships among ODE solution components  
* General-purpose MATLAB implementation  
* Applicable to nonlinear dynamical systems  
* Visualization of identified relationships and constraint preservation 

**Repository Structure:**
1. `functional_relationships.m` — Main algorithm implementation  
2. `enzyme.m` — Enzyme dynamics model (Subsection 4.1)  
3. `glycolytic.m` — Glycolytic oscillator model (Subsection 4.2) 

**Code description:**
The MATLAB function that implements the algorithm presented in the paper is named *functional\_relationships*. 
To ensure general usability, the dependent variables in all models are uniformly labeled using the notation X. 
The input of the function consist of odesys, t0, X0, T, m, p, S, where:
odesys is the MATLAB function for the ODE system
t0 is the initial value of the independent variable
x0 is the initial vector value of the dependent variable at t=t0
T is the end value of the independent variable
m is the number of partition intervals of [t0, T]
p is the number of Step 3 repetitions
S is the number of particular xi solutions - if they exist- the user wants

The output consists of the relationship between the solution components and graphical evidence of constraint preservation.

Comprehensive usage instructions and explanatory comments are included within the function's body.

The enzyme dynamics model described by equations (9)-(12) in Subsection 4.1 of the paper is defined by the MATLAB function 
called *enzyme*. The input of this function consist of time, t, and vector function X, whose components X(1), X(2), X(3), and X(4) 
correspond to the concentrations C_1(t), E(t), S(t), and P(t), respectively.

The glycolytic oscillator model described by equations (16)-(22) in Subsection 4.2 of the paper is defined by the MATLAB function 
called *glycolytic*. The input of this function consist of time, t, and vector function X, whose components X(1),..., X(7) 
correspond to the concentrations of seven biochemical species S_1(t),..., S_7(t), respectively.

**Instructions on how to run the MATLAB code:**

Concrete examples of use:

1. For numerical experiments related to the enzyme dynamics first type at the prompt in the Command Window

*functional\_relationships(@enzyme,0,\[1 0 1 1],1,100,3,6)*

where the input quantities are defined in the code description above.

When prompted: *Do you want to use powers, i.e., X^k? Answer 1 for Yes and 0 for No.* Enter 1

Next, specify powers, e.g.:

*For monomials X^k the powers k =* [0 1]

The output will consist of the general solution for coefficients and six particular functional relationships, together with their corresponding graphs.


2. For numerical experiments related to the glycolytic oscillator model first type at the prompt in the Command Window

*functional_relationships(@glycolytic,0,[1 1 0.1 0.2 0.2 1 0.1],1,10000,3,4)*

where the input quantities are defined in the code description above.

When propted: *Do you want to use powers, i.e., X^k? Answer 1 for Yes and 0 for No.* Enter 1

Next, specify powers, e.g.:

*For monomials X^k the powers k =* [0 1 2]

The output will consist of the general solution for coefficients and four particular functional relationships, together with their corresponding graphs.
