##### Introduction

Single header file solver for linear systems ($Mx=y$) with tridiagonal matrices.

The square (row major) matrix is zero where no element is defined:

$M=\begin{bmatrix}d_1&u_1&&&\alpha\cr{l_2}&d_2&u_2&&\cr&...&...&...&\cr&&l_{n-1}&d_{n-1}&u_{n-1}\cr{\beta}&&&l_n&d_n\end{bmatrix}$

Lower/sub- ($l_i$), upper/super- ($u_i$) and main- ($d_i$) diagonals. And cyclic elements $\alpha$ and $\beta$

With $\alpha=\beta=0$, the system is solved using the Llewellyn Thomas algortihm.


If $\alpha\ne0$ and/or $\beta\ne0$, the matrix is cyclic, and is solved using the Sherman-Morrison formula.

The cyclic solver uses the property that the matrix can be wrttien as $M = A+\vec{u}\vec{v}^T$. With $A$ being a non-cyclic matrix, and $\gamma=d_1$, $\vec{u}=\{\gamma,0,...,0,\beta\}$, $\vec{v}=\{1,0,...,0,\alpha/\gamma\}$.

##### Usage

Solve a cyclic tridiagonal matrix of order 5.

$\begin{bmatrix}2&2&0&0&2\cr3&2&2&0&0\cr0&2&8&7&0\cr0&0&2&3&2\cr3&0&0&7&5\end{bmatrix}\vec{x}=\begin{pmatrix}6\cr-6\cr1\cr2\cr3\end{pmatrix}$

	std::vector<std::double_t> upper{ 2, 3, 7, 2, 0 };
	std::vector<std::double_t> main{ 2, 2, 8, 3, 5 };
	std::vector<std::double_t> lower{ 0, 3, 2, 2, 7 };
	std::double_t const alpha = 2;
	std::double_t const beta = 3;

	std::vector<std::double_t> y{ 6, -6, 1, 2, 3 };

	std::vector<std::double_t> x = Tridiagonal::Solve( upper, main, lower, alpha, beta, y );

	if ( !x.size() || !Tridiagonal::ValidSolution( upper, main, lower, alpha, beta, y, x ) )
		std::cout << "No valid solution found.\n";
	else {
		std::cout << "Solution found:\n";
		for ( std::uint32_t i{ 0 }; i < x.size(); ++i )
			std::cout << "x[" << i + 1 << "] = " << x[i] << '\n';
	}

Output:

	Solution found:
	x[1] = -3.48077
	x[2] = 4.86923
	x[3] = -1.76538
	x[4] = 0.769231
	x[5] = 1.61154

##### Requirement

C++11
