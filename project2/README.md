<h1>Project 2: Solving the Heat Equation with Finite Element Methods</h1>

<h2>About</h2>
This project uses the finite element Galerkin method and Euler time discretization to solve the heat equation at a given time, given a forcing function.

<h2>Included Files</h2>
<ul>
<li>feFunctions.py: contains supplementary functions for the finite element method</li>
<li>fe.py: code that allows user to input a number of spatial and time nodes and to choose between forward or backward Euler method; it then computes and plots the solution for the heat equation at the final timestep.</li>
</ul>

<h2>Executing Scripts</h2>
First, download both feFunctions.py and fe.py to the same directory in an environment that can both run Python 3 and has numpy installed. (Running the command <code>pip install numpy</code> may be necessary.)<br>
Then, run <code>python3 fe.py</code>, and enter inputs for number of spatial nodes, timesteps, and foward or backward Euler method according to the input prompts. The output will include a plot of the computed solution to the heat equation, compared to the analytic solution.<br>

<h2>Sample Output</h2>
<pre>
Enter number of nodes (positive integer): 
11
Enter number of timesteps (positive integer): 
551
Forward or backward Euler ("f"/"b"): 
f
<img src="https://github.com/bethanygrimm/coe352/blob/b0490cd3e32d1151b6753a90f906c644aa6c7da7/project2/HeatEquationSolution.png" alt="Computed Solution and Analytic Solution Plot">
</pre>
