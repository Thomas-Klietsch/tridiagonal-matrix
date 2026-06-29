// Copyright (c) 2026 Thomas Klietsch, all rights reserved.
//
// Licensed under the GNU Lesser General Public License, version 3.0 or later
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or ( at your option ) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General
// Public License along with this program.If not, see < https://www.gnu.org/licenses/>. 

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

// Square matrices, of at least order 5,
// defined by three (3) diagonal vectors.
namespace Tridiagonal {

	// Values less than this are considered to be zero
	constexpr std::double_t EPSILON_ZERO = 1.e-14;

	// Solve cyclic, and non-cyclic, tridiagonal matrices
	std::vector<std::double_t> Solve(
		std::vector<std::double_t> upper,
		std::vector<std::double_t> diagonal,
		std::vector<std::double_t> lower,
		std::double_t const alpha,
		std::double_t const beta,
		std::vector<std::double_t> y
	) {
		if ( std::abs( diagonal[0] ) < EPSILON_ZERO ) {
			std::cout << "Error! First diagonal element is zero.\n";
			return {};
		}

		std::uint32_t const n = y.size();

		if ( n < 5 ) {
			std::cout << "Error! Matrix needs to be at least order 5.\n";
			return {};
		}

		if ( n != upper.size() || n != diagonal.size() || n != lower.size() ) {
			std::cout << "Error! All input data must be of the same size.\n";
			return {};
		}

		// Lambda function
		// Solve non-cyclic tridiagonal, using (Llewellyn) Thomas algorithm.
		auto solve_non_cyclic = [n](
			// Upper diagonal, last element is ignored
			std::vector<std::double_t> const& upper,
			// Diagonal
			std::vector<std::double_t> const& diagonal,
			// Lower diagonal, first element is ignored
			std::vector<std::double_t> const& lower,
			// Result of matrix multiplied by x
			std::vector<std::double_t> const& y,
			// Check determinant if true (default)
			bool const f_determinant = true
			) -> std::vector<std::double_t>
		{
			if ( std::abs( diagonal[0] ) < EPSILON_ZERO ) {
				std::cout << "Error! First diagonal element is zero.\n";
				return {};
			}

			// Calculate determinant, using tridiagonal optimisation.
			if ( f_determinant ) {
				std::vector<std::double_t> K( n, 0 );
				K[0] = 1;
				K[1] = diagonal[0];
				for ( std::uint32_t i{ 2 }; i < n; ++i )
					K[i] = diagonal[i - 1] * K[i - 1] - upper[i - 2] * lower[i - 1] * K[i - 2];
				if ( std::abs( K[n - 1] ) < EPSILON_ZERO ) {
					std::cout << "Error! Determinant is zero.\n";
					return {};
				}
			}

			// Solution vector
			std::vector<std::double_t> x( n, 0 );

			// New upper diagonal matrix, last element is unused
			std::vector<std::double_t> modified_upper( n, 0 );

			// Stage 1: Reduce tridiagonal to upper triangular matrix
			// Initial values
			x[0] = y[0] / diagonal[0];
			modified_upper[0] = upper[0] / diagonal[0];
			// 1 to n-1
			for ( std::uint32_t i{ 1 }; i < n; ++i ) {
				std::double_t const temp = ( diagonal[i] - lower[i] * modified_upper[i - 1] );
				modified_upper[i] = upper[i] / temp;
				x[i] = ( y[i] - lower[i] * x[i - 1] ) / temp;
			}

			// Stage 2: Back substitute
			// n-2 to 0
			for ( std::uint32_t i{ n - 1 }; i-- > 0; )
				x[i] -= modified_upper[i] * x[i + 1];

			return x;
		}; // end of lambda function

		// Test if system can be solve with non-cyclic solver (i.e. Thomas algorithm).
		if ( std::abs( alpha ) < EPSILON_ZERO && std::abs( beta ) < EPSILON_ZERO )
			return solve_non_cyclic( upper, diagonal, lower, y, true );

		// Solve using cyclic method, Sherman-Morrison algorithm
		// Jack Sherman and Winifred J. Morrison

		std::double_t const gamma = -diagonal[0]; // Top left

		// u vector
		std::vector<std::double_t> u( n, 0 );
		u[0] = gamma;
		u[n - 1] = beta;
		// v vector = { 1, 0 ... 0, alpha/gamma }
		// M = A + uv^T

		// New diagonal, matrix A
		std::vector<std::double_t> modified_main{ diagonal };
		modified_main[0] = diagonal[0] - gamma;
		modified_main[n - 1] = diagonal[n - 1] - alpha * beta / gamma;

		// Solve Ax=y
		auto x = solve_non_cyclic( upper, modified_main, lower, y, true );
		// Check for no solution
		if ( !x.size() )
			return {};

		// Solve Az=u
		auto z = solve_non_cyclic( upper, modified_main, lower, u, false ); // Det. is already checked above
		// Check for no solution
		if ( !z.size() )
			return {};

		// Sherman-Morrison, (x dot v)/(1 + z dot v)
		std::double_t const scalar = ( x[0] + x[n - 1] * alpha / gamma ) / ( 1. + z[0] + z[n - 1] * alpha / gamma );
		for ( std::uint32_t i{ 0 }; i < n; ++i )
			x[i] -= z[i] * scalar;

		return x;
	};

	// The solver might return values even if solution is invalid.
	// Checks if it is valid (Mx=y)
	bool ValidSolution(
		std::vector<std::double_t> const& upper,
		std::vector<std::double_t> const& diagonal,
		std::vector<std::double_t> const& lower,
		std::double_t const alpha,
		std::double_t const beta,
		std::vector<std::double_t> const& y,
		std::vector<std::double_t> const& x
	) {
		std::uint32_t const n = y.size();

		if ( n < 5 ) {
			std::cout << "Error! Matrix needs to be at least order 5.\n";
			return false;
		}

		if ( n != upper.size() || n != diagonal.size() || n != lower.size() || n != x.size() ) {
			std::cout << "Error! All input data must be of the same size.\n";
			return false;
		}

		// First row
		std::double_t value;

		value = diagonal[0] * x[0] + upper[0] * x[1] + alpha * x[n - 1];
		if ( std::abs( y[0] - value ) > EPSILON_ZERO )
			return false;

		// Intermediate rows
		for ( std::uint32_t i{ 1 }; i < n - 1; ++i ) {
			value = lower[i] * x[i - 1] + diagonal[i] * x[i] + upper[i] * x[i + 1];
			if ( std::abs( y[i] - value ) > EPSILON_ZERO )
				return false;
		}

		// Last row
		value = beta * x[0] + lower[n - 1] * x[n - 2] + diagonal[n - 1] * x[n - 1];
		if ( std::abs( y[n - 1] - value ) > EPSILON_ZERO )
			return false;

		return true;
	};

};
