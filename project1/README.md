<h1>Project 1: Solving a Spring-Mass System</h1>

<h2>About</h2>
This project utilizes SVD-decomposition to solve for unknowns in a simple, user-defined spring-mass system.

<h2>Included Files</h2>
<ul>
<li>svd.py: contains supplementary functions to compute SVD-decompositions</li>
<li>springSystem.py: code that allows user to input values in a spring-mass system, and then returns the system's resultant displacements, elongations, internal forces, singular values, eigenvalues, and L2-condition number.</li>
</ul>

<h2>Executing Scripts</h2>
First, download both svd.py and springSystem.py to the same directory in an environment that can both run Python 3 and has numpy installed. (Running the command <code>pip install numpy</code> may be necessary.)<br>
Then, run <code>python3 springSystem.py</code>, and enter inputs for number of masses, fixed ends, masses, and spring constants according to the input prompts. The output will include a displacement vector of the masses, an elongation vector of the springs, an internal force vector, the K matrix's singular values, eigenvalues, and L2-condition number.<br>

<h2>Sample Output</h2>
Example - one fixed end:<br>
<pre>
Enter number of masses (positive integer):
3
Enter number of fixed ends (1 or 2):
1
Enter mass 1 (positive number):
0.5
Enter mass 2 (positive number):
0.2
Enter mass 3 (positive number):
0.1
Enter spring constant 1 (positive number):
10
Enter spring constant 2 (positive number):
12
Enter spring constant 3 (positive number):
8
Displacement vector of masses: [0.7848   1.03005  1.152675]
Elongation vector of springs: [0.7848   0.24525  0.122625]
Internal stresses in springs: [7.848 2.943 0.981]
Singular values of K: [34.22235539 13.73533225  2.04231237], and eigenvalues of K: [1171.16960828  188.65935192    4.1710398 ]
L2 condition number of K: 18.314240057640827
</pre><br>
Example - two fixed ends:<br>
<pre>
Enter number of masses (positive integer):
4
Enter number of fixed ends (1 or 2):
2
Enter mass 1 (positive number):
7
Enter mass 2 (positive number):
6
Enter mass 3 (positive number):
7
Enter mass 4 (positive number):
5
Enter spring constant 1 (positive number):
10
Enter spring constant 2 (positive number):
12
Enter spring constant 3 (positive number):
13
Enter spring constant 4 (positive number):
10
Enter spring constant 5 (positive number):
14
Displacement vector of masses: [12.20110946 16.64620068 16.2216695   8.80277896]
Elongation vector of springs: [12.20110946  4.44509122 -0.42453118 -7.41889054 -8.80277896]
Internal stresses in springs: [ 122.01109461   53.34109461   -5.51890539  -74.18890539 -123.23890539]
Singular values of K: [43.14102584 29.52614229  4.36110258 16.97172929], and eigenvalues of K: [1861.14811039  871.79307881   19.01921568  288.03959512]
L2 condition number of K: 13.248128854512569
</pre>
