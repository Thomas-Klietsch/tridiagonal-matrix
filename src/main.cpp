#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "./tridiagonal-matrix.hpp"

void ZeroDeterminant() {
	std::cout << "Test zero determinant.\n";
	// Upper diagonal (last element is unused)
	std::vector<std::double_t> a{ 2, 3, 7, 1, 0 };
	// Diagonal
	std::vector<std::double_t> b{ 2, 2, 8, 1, 5 };
	// Lower diagonal (first element is unused)
	std::vector<std::double_t> c{ 0, 3, 2, 2, 7 };

	// Vector, M x = y
	std::vector<std::double_t> y{ 6, 6, 6, 6, 6 };

	// Solve M x = y, returns x (vector)
	auto x = Tridiagonal::Solve( a, b, c, 0, 0, y );

	if ( !x.size() || !Tridiagonal::ValidSolution( a, b, c, 0, 0, y, x ) )
		std::cout << "No valid solution found.\n";
	else {
		std::cout << "Solution found:\n";
		for ( std::uint32_t i{ 0 }; i < x.size(); ++i )
			std::cout << "x[" << i + 1 << "] = " << x[i] << '\n';
	}
	std::cout << '\n';
};

void SolveNonCyclic() {
	std::cout << "Solve a (non-cyclic) tridiagonal of order 5. \n";
	// Upper diagonal (last element is unused)
	std::vector<std::double_t> a{ 2, 3, 7, 2, 0 };
	// Diagonal
	std::vector<std::double_t> b{ 2, 2, 8, 3, 5 };
	// Lower diagonal (first element is unused)
	std::vector<std::double_t> c{ 0, 3, 2, 2, 7 };

	std::vector<std::double_t> y{ 6, 6, 6, 6, 6 };

	// Solve for x; M x = y
	auto x = Tridiagonal::Solve( a, b, c, 0, 0, y );

	if ( !x.size() || !Tridiagonal::ValidSolution( a, b, c, 0, 0, y, x ) )
		std::cout << "No valid solution found.\n";
	else {
		std::cout << "Solution found:\n";
		for ( std::uint32_t i{ 0 }; i < x.size(); ++i )
			std::cout << "x[" << i + 1 << "] = " << x[i] << '\n';
	}

	std::cout << '\n';
};

void SolveCyclic() {
	std::cout << "Solve a cyclic tridiagonal of order 5. \n";
	// Upper diagonal (last element is unused)
	std::vector<std::double_t> upper{ 2, 3, 7, 2, 0 };
	// Diagonal
	std::vector<std::double_t> main{ 2, 2, 8, 3, 5 };
	// Lower diagonal (first element is unused)
	std::vector<std::double_t> lower{ 0, 3, 2, 2, 7 };
	// Top right element, first row last column
	std::double_t const alpha = 2;
	// Lower left element, last row first column
	std::double_t const beta = 3;

	std::vector<std::double_t> y{ 6, -6, 1, 2, 3 };

	// Solve for x; M x = y
	std::vector<std::double_t> x = Tridiagonal::Solve( upper, main, lower, alpha, beta, y );

	if ( !x.size() || !Tridiagonal::ValidSolution( upper, main, lower, alpha, beta, y, x ) )
		std::cout << "No valid solution found.\n";
	else {
		std::cout << "Solution found:\n";
		for ( std::uint32_t i{ 0 }; i < x.size(); ++i )
			std::cout << "x[" << i + 1 << "] = " << x[i] << '\n';
	}
	std::cout << '\n';
};

int main(
	[[maybe_unused]] int argc,
	[[maybe_unused]] char* argv[]
) {
	ZeroDeterminant();

	SolveNonCyclic();

	SolveCyclic();

	return EXIT_SUCCESS;
};
