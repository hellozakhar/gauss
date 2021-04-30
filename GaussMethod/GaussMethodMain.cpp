#include "GaussMethodTest.h"

int main() {
	// isAuto true  => coeffs and y - generated,
	// isAuto false => coeffs and y - enter manually.
	// const bool isAuto = true;

	// isForCopy true => print in console each avg time only
	// isForCopy false => print size of matrices, min/max/avg times
	const bool isForCopy = true;

	//const int n = 100;
	const int times = 5;

	//testGaussDef(isAuto);
	//testGaussDiscrepancy();
	//testCompareSolutionGaussAndLU();

	//testDecomposeSolveV1(n);
	//testDecomposeSolveV2(n);

	/*testDecomposeSolveV1(100, times, isForCopy);
	testDecomposeSolveV1(200, times, isForCopy);
	testDecomposeSolveV1(500, times, isForCopy);
	testDecomposeSolveV1(1000, times, isForCopy);
	testDecomposeSolveV1(1500, times, isForCopy);
	testDecomposeSolveV1(2000, times, isForCopy);
	testDecomposeSolveV1(3000, times, isForCopy);
	testDecomposeSolveV1(4000, times, isForCopy);
	testDecomposeSolveV1(5000, times, isForCopy);
	testDecomposeSolveV1(6000, times, isForCopy);
	testDecomposeSolveV1(7000, times, isForCopy);*/

	/*testDecomposeSolveV2(100, times, isForCopy);
	testDecomposeSolveV2(200, times, isForCopy);
	testDecomposeSolveV2(500, times, isForCopy);
	testDecomposeSolveV2(1000, times, isForCopy);
	testDecomposeSolveV2(1500, times, isForCopy);
	testDecomposeSolveV2(2000, times, isForCopy);
	testDecomposeSolveV2(3000, times, isForCopy);
	testDecomposeSolveV2(4000, times, isForCopy);
	testDecomposeSolveV2(5000, times, isForCopy);
	testDecomposeSolveV2(6000, times, isForCopy);
	testDecomposeSolveV2(7000, times, isForCopy);*/

	testDecomposeSolveV2Par(100, times, isForCopy);
	testDecomposeSolveV2Par(200, times, isForCopy);
	testDecomposeSolveV2Par(500, times, isForCopy);
	testDecomposeSolveV2Par(1000, times, isForCopy);
	testDecomposeSolveV2Par(1500, times, isForCopy);
	testDecomposeSolveV2Par(2000, times, isForCopy);
	testDecomposeSolveV2Par(3000, times, isForCopy);
	testDecomposeSolveV2Par(4000, times, isForCopy);
	/*testDecomposeSolveV2Par(5000, times, isForCopy);
	testDecomposeSolveV2Par(6000, times, isForCopy);
	testDecomposeSolveV2Par(7000, times, isForCopy);*/

	int zzz;
	std::cin >> zzz;
	return 0;
}