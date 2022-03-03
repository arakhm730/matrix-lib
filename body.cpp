#include "Source.cpp"


int main() {
	setlocale(LC_ALL, "RUS");


	std::ifstream ifs("scores.txt");
	Matrix<double> txt;
	ifs >> txt;
	ifs.close(); 
	
	std::ofstream ofs("out.txt");
	ofs.clear();
	ofs << txt;
	ofs.close();



	std::ofstream out("numbers", std::ios::out | std::ios::binary);
	out.clear(); 
	Matrix<double> ot({ {9.5, -3.4, 1.0, 15}, {9.5, -3.4, 1.0, 17}, {9.5, -3.4, 40, 0}, {9.5, -3.4, 1.0, 13} });
	ot.save_to(out);
	out.close();

	std::ifstream in("numbers", std::ios::in | std::ios::binary);
	Matrix<double> ii(4, 4);
	ii.load_from(in);
	std::cout << ii << std::endl;

	in.close();


	/*Matrix<double> k(1, 3), l(1, 3);
	std::cin >> k;
    std::cin >> l;
	std::cout << "Rank of vector is " << k.rank() << std::endl;
	std::cout << "scalar product of vectors is " << scalarProduct(k, l) << std::endl;
	std::cout << "angle between vectors is " << vectors_angle(k, l) << " rad" << std::endl;
	std::cout << " vector max norm of first is " << k.vecNormMax() << " and vector norm euclid of second is " << l.vecNormEuclid() << std::endl;
	Matrix<double> a(3, 3);
	std::cin >> a;
	std::cout << a << std::endl;
	std::cout << "Matrix norm is " << a.matrixNorm() << std::endl;
	std::cout << "And matrix rank is " << a.rank() << std::endl;*/
	/*Matrix<double> b(3, 3);
	std::cin >> b;
	std::cout << "Rang is " << b.rank() << std::endl;
	std::cout << "trace is " << b.trace() << " and determinant is ";
	std::cout << b.determinant() << std::endl;
	std::cout << b.transpose() << std::endl;
	std::cout << b << std::endl;*/
	Matrix<double> abc({ {1, 4, 2, 3}, {1, -5, 2, 1}, {2, 13, 1, 1} });
	std::cout << abc;
	std::cout << abc.transpose();
	std::cout << ot.determinant();
	//Matrix<double> inv = abc.inverse();
	//std::cout << inv * abc << std::endl; 
	/*Matrix<int> p(3, 3);
	std::cin >> p;
	std::cout << p << std::endl;
	IdentityMatrix e(4);
	std::cout << 2*e << std::endl;
	DiagonalMatrix<double> d(b);
	UpperTriangularMatrix<double> up(a);
	LowerTriangularMatrix<double> down(b);
	SymmetricalMatrix<int> sym(p);
	std::cout << sym << std::endl;
	std::cout << sym.transpose() << std::endl;
	std::cout << up << down << std::endl;
	std::cout << 3.0 * d << std::endl;
	std::cout << "Adamar mult:" << std::endl;
	std::cout << AdamarMult(a, b) << std::endl;*/
	return 0;
}