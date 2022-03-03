#include "Header.h"


template<class T>
Matrix<T>::Matrix(const int& rows, const int& cols) : Row(rows), Col(cols) {
	try
	{
		if (rows < 0 || cols < 0) throw(123);
		mtr.assign(rows, std::vector<T>(cols));
	}
	catch (const int) {
		std::cout << "Incorrect number of rows or columns, should be greater than 0!" << std::endl;
	}
}

template <class T>
Matrix<T>::Matrix(std::vector <std::vector <T>> data) : mtr(data), Row(data.size()), Col((data.at(0).size() != 0) ? data.at(0).size() : 0) {}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& sample) {
	int rows = sample.GetRow();
	int cols = sample.GetCol();

	try
	{
		if (rows < 0 || cols < 0) throw('e');
		mtr.assign(rows, std::vector<T>(cols));
	}
	catch (const char) {
		std::cout << "Incorrect number of rows or columns, should be greater than 0!" << std::endl;
	}

	this->Row = rows;
	this->Col = cols;

	

	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			mtr.at(i).at(j) = sample.mtr.at(i).at(j);
}

template<class T>
int Matrix<T>::GetCol() const {
	return Col;
}

template<class T>
void Matrix<T>::SetRow(const int & row)
{
	Row = row;
}

template<class T>
void Matrix<T>::SetCol(const int & col)
{
	Col = col;
}

template<class T>
void Matrix<T>::SetMtr(std::vector<std::vector<T>> data)
{
	mtr = data;
}

template<class T>
int Matrix<T>::GetRow() const {
	return Row;
}

template<class T>
std::vector<std::vector<T>> Matrix<T>::GetMtrData() const {
	return mtr;
}

template<class T>
double Matrix<T>::trace() const {
	try {
		if (this->GetCol() != this->GetRow()) throw (2.34);
		else {
			double tr = 0.0;
			for (int i = 0; i < this->GetRow(); i++) tr += this->mtr.at(i).at(i);
			return tr;
		}
	}
	catch (double) { std::cout << "Cannot calculate trace of non-square matrix!" << std::endl; }
}

template<class T>
double Matrix<T>::vecNormMax() const {
	try {
		if (this->Row != 1) throw(-5);
		else {
			double max = 0.0;
			for (int j = 0; j < this->Col; j++)
				if ((double)this->mtr.at(0).at(j) >= max) max = this->mtr.at(0).at(j);
			return max;
		}
	}
	catch (int) { std::cout << "Cannot calculate vector norm of non-vector object!" << std::endl; }
}

template<class T>
double Matrix<T>::vecNormEuclid() const {
	try {
		if (this->Row != 1) throw(-5);
		else {
			T squares = 0.0;
			for (int j = 0; j < this->Col; j++) squares += pow(this->mtr.at(0).at(j), 2);
			return (double)pow(squares, 0.5);
		}
	}
	catch (int) { std::cout << "Cannot calculate vector norm of non-vector object!" << std::endl; }
}

template<class T>
double scalarProduct(const Matrix<T>& first, const Matrix<T>& second) {
	try {
		if (((first.GetRow() == second.GetRow() && second.GetRow() == 1) && first.GetCol() == second.GetCol() )) {
			double product = 0.0;
				for (int i = 0; i < first.GetCol(); i++) product += first.GetMtrData().at(0).at(i) * second.GetMtrData().at(0).at(i);
			return product;
		}
		else { throw (1); }
	}
	catch (int) { std::cout << "Vectors have inequal sizes, cannot calculate scalar product!" << std::endl; }
}

template <class T>
double Matrix<T>::matrixNorm() const {
	T squares = 0;
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			squares += pow(this->mtr.at(i).at(j), 2);
	return pow((double)squares, 0.5);
}

template<class T>
double vectors_angle(const Matrix<T>& first, const Matrix<T>& second) {
	double scalar = scalarProduct(first, second);
	return acos(scalar / (first.vecNormEuclid() * second.vecNormEuclid()));
}

template<class T>
double det(const Matrix<T>& inst, const int n) {
	try {
		if (inst.GetCol() != inst.GetRow()) throw (2.34);
		else {
			double deter = 0;
			Matrix<T> submatrix(n, n);
			if (n == 2)
				return ((inst.GetMtrData().at(0).at(0) * inst.GetMtrData().at(1).at(1)) - (inst.GetMtrData().at(1).at(0) * inst.GetMtrData().at(0).at(1)));
			else {
				for (int x = 0; x < n; x++) {
					int subi = 0;
					for (int i = 1; i < n; i++) {
						int subj = 0;
						for (int j = 0; j < n; j++) {
							if (j == x)
								continue;
							submatrix.GetMtrData().at(subi).at(subj) = inst.GetMtrData().at(i).at(j);
							subj++;
						}
						subi++;
					}
					deter = deter + (pow(-1, x) * inst.GetMtrData().at(0).at(x) * det(submatrix, n - 1));
				}
			}
			return deter;
		}
	}
	catch (double) { std::cout << "Cannot calculate determinant of non-square matrix!" << std::endl; }
}

template <class T>
double Matrix<T>::determinant() {
	Matrix<double> copy(this->Row, this->Col);
	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			copy.mtr[i][j] = this->mtr[i][j];
	return getDeterminant(copy.mtr);
}



template<class T>
Matrix<T> Matrix<T>::transpose() {
	std::vector <std::vector<T>> copy(this->Col, std::vector<T>());
	for(int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++) {
			copy.at(j).push_back(this->mtr.at(i).at(j));
		}
	Matrix<T> c(copy);
	return c;
}

template <class T>
double Matrix<T>::getDeterminant(const std::vector<std::vector<double>> vect) {
	if (vect.size() != vect[0].size()) {
		throw std::runtime_error("Matrix is not quadratic");
	}
	int dimension = vect.size();

	if (dimension == 0) {
		return 1;
	}

	if (dimension == 1) {
		return vect[0][0];
	}

	//Formula for 2x2-matrix
	if (dimension == 2) {
		return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
	}

	double result = 0;
	int sign = 1;
	for (int i = 0; i < dimension; i++) {

		//Submatrix
		std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double>(dimension - 1));
		for (int m = 1; m < dimension; m++) {
			int z = 0;
			for (int n = 0; n < dimension; n++) {
				if (n != i) {
					subVect[m - 1][z] = vect[m][n];
					z++;
				}
			}
		}

		//recursive call
		result = result + sign * vect[0][i] * getDeterminant(subVect);
		sign = -sign;
	}

	return result;
}

template <class T>
std::vector<std::vector<double>> Matrix<T>::getTranspose(const std::vector<std::vector<double>> matrix1) {

	//Transpose-matrix: height = width(matrix), width = height(matrix)
	std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double>(matrix1.size()));

	//Filling solution-matrix
	for (size_t i = 0; i < matrix1.size(); i++) {
		for (size_t j = 0; j < matrix1[0].size(); j++) {
			solution[j][i] = matrix1[i][j];
		}
	}
	return solution;
}

template <class T>
std::vector<std::vector<double>> Matrix<T>::getCofactor(const std::vector<std::vector<double>> vect) {
	try {
		if (vect.size() != vect[0].size()) {
			throw std::runtime_error("Matrix is not quadratic");
		}

		std::vector<std::vector<double>> solution(vect.size(), std::vector<double>(vect.size()));
		std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double>(vect.size() - 1));

		for (std::size_t i = 0; i < vect.size(); i++) {
			for (std::size_t j = 0; j < vect[0].size(); j++) {

				int p = 0;
				for (size_t x = 0; x < vect.size(); x++) {
					if (x == i) {
						continue;
					}
					int q = 0;

					for (size_t y = 0; y < vect.size(); y++) {
						if (y == j) {
							continue;
						}

						subVect[p][q] = vect[x][y];
						q++;
					}
					p++;
				}
				solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
			}
		}
		return solution;
	}
	catch (std::runtime_error) { std::cout << "Matrix is not quadratic" << std::endl; }
	std::vector<std::vector<double>> A{ 0 };
	return A;
}

template <class T>
std::vector<std::vector<double>> Matrix<T>::getInverse(const std::vector<std::vector<double>> vect) {
		if (getDeterminant(vect) == 0) {
			throw std::runtime_error("Determinant is 0");
		}

		double d = 1.0 / getDeterminant(vect);
		std::vector<std::vector<double>> solution(vect.size(), std::vector<double>(vect.size()));

		for (size_t i = 0; i < vect.size(); i++) {
			for (size_t j = 0; j < vect.size(); j++) {
				solution[i][j] = vect[i][j];
			}
		}

		solution = getTranspose(getCofactor(solution));

		for (size_t i = 0; i < vect.size(); i++) {
			for (size_t j = 0; j < vect.size(); j++) {
				solution[i][j] *= d;
			}
		}
		return solution;
}

template <class T>
Matrix<double> Matrix<T>::inverse() {
	try {
		Matrix<double> copy(this->Row, this->Col);
		for (int i = 0; i < Row; i++)
			for (int j = 0; j < Col; j++)
				copy.mtr[i][j] = double(this->mtr[i][j]);
		Matrix<double> inversed(getInverse(copy.GetMtrData()));
		return inversed;
	}
	catch (std::runtime_error) { std::cout << "Determinant is 0, matrix is non-inversible!" << std::endl; }
	std::vector<std::vector<double>> A{ {0} };
	return Matrix<double>(A);
}

const double EPS = 1E-9;

template <class T>
int Matrix<T>::compute_rank(std::vector<std::vector<T>> A) {
	int n = A.size();
	int m = A[0].size();

	int rg = 0;
	std::vector<bool> row_selected(n, false);
	for (int i = 0; i < m; ++i) {
		int j;
		for (j = 0; j < n; ++j) {
			if (!row_selected[j] && abs(A[j][i]) > EPS)
				break;
		}

		if (j != n) {
			++rg;
			row_selected[j] = true;
			for (int p = i + 1; p < m; ++p)
				A[j][p] /= A[j][i];
			for (int k = 0; k < n; ++k) {
				if (k != j && abs(A[k][i]) > EPS) {
					for (int p = i + 1; p < m; ++p)
						A[k][p] -= A[j][p] * A[k][i];
				}
			}
		}
	}
	return rg;
}

template <class T>
int Matrix<T>::rank() {
	return compute_rank(this->mtr);
}

template<class T>
Matrix<T> AdamarMult(Matrix<T>& first, Matrix<T>& second) {
	try {
		if (first.Row != second.Row || first.Col != second.Col) throw("help");

			Matrix<T> c(first.Row, first.Col);
			for (int i = 0; i < first.Row; i++)
				for (int j = 0; j < first.Col; j++)
					c.mtr.at(i).at(j) = first.mtr.at(i).at(j) * second.mtr.at(i).at(j);
			return c;
		
	}
	catch (const char*) {
	std::cout << "Incompatible matrix sizes!" << std::endl;
	}
}

template<class T>
Matrix<T> operator+ (Matrix<T>& first, Matrix<T>& second) {
	if (first.Row == second.Row && first.Col == second.Col) {
		Matrix<T> c(first.Row, first.Col);
		for (int i = 0; i < first.Row; i++)
			for (int j = 0; j < first.Col; j++)
				c.mtr.at(i).at(j) = first.mtr.at(i).at(j) + second.mtr.at(i).at(j);
		return c;
	}
	else std::cout << "Incompatible matrix sizes!" << std::endl;

}

template<class T>
Matrix<T> operator- (Matrix<T>& first, Matrix<T>& second) {
	if (first.Row == second.Row && first.Col == second.Col) {
		Matrix<T> c(first.Row, first.Col);
		for (int i = 0; i < first.Row; i++)
			for (int j = 0; j < first.Col; j++)
				c.mtr.at(i).at(j) = first.mtr.at(i).at(j) - second.mtr.at(i).at(j);
		return c;
	}
	else std::cout << "Incompatible matrix sizes!" << std::endl;

}

template<class T>
Matrix<T> operator*(const T& arg, Matrix<T>& inst) {
	Matrix<T> multi(inst.GetRow(), inst.GetCol());
	for (int i = 0; i < inst.GetRow(); i++)
		for (int j = 0; j < inst.GetCol(); j++)
			multi.mtr.at(i).at(j) = inst.mtr.at(i).at(j) * arg;
	return multi;
}

template<class T>
Matrix<T> operator*(Matrix<T>& inst, const T& arg) {
	Matrix<T> multi(inst.GetRow(), inst.GetCol());
	for (int i = 0; i < inst.GetRow(); i++)
		for (int j = 0; j < inst.GetCol(); j++)
			multi.mtr.Matrix(i).at(j) = inst.mtr.at(i).at(j) * arg;
	return multi;
}

template<class T>
Matrix<T> operator*(Matrix<T>& inst1, Matrix<T>& inst2) {
	if (inst1.GetCol() == inst2.GetRow()) {
		
		Matrix<T> mult(inst1.GetRow(), inst2.GetCol());
		for (int i = 0; i < inst1.GetRow(); i++) {
			for (int j = 0; j < inst2.GetCol(); j++) {
				T elm = 0;
				for (int k = 0; k < inst1.GetRow(); k++)
					elm += inst1.mtr.at(i).at(k) * inst2.mtr.at(k).at(j);

				mult.mtr.at(i).at(j) = elm;
			}
		}
		return mult;
	}
	else std::cout << "Incompatible matrices size!" << std::endl;
	
}




template<class T>
std::ostream& operator<< (std::ostream& out, const Matrix<T>& instance) {
	out << "Current matrix data:\n";
	int Row = instance.GetRow();
	int Col = instance.GetCol();
	for (int i = 0; i < Row; i++)
	{
		if (i == 0) out << "/ ";
		else if (i == Row - 1) out << "\\ ";
		else out << "| ";
		for (int j = 0; j < Col; j++) {
			(j!=0) ? out << '\t' : out << "";
			out << instance.mtr.at(i).at(j) << std::setw(2);
			//if (j != Col - 1) out
		}


		if (i == 0) out << " \\";
		else if (i == Row - 1) out << " /";
		else out << " |";
		out << "\n";
	}
	return out;
}

template <class T>
std::istream& operator>> (std::istream& in, Matrix<T>& dest) {
	int Row = dest.GetRow();
	int Col = dest.GetCol();
	std::cout << "Please, type in " << Row << " rows of matrix elemets, each has to have" << Col << " elements:\n";
	for (int i = 0; i < Row; i++)
		for (int j = 0; j < Col; j++)
			in >> dest.mtr.at(i).at(j);
	return in;
}

template <class T>
std::ofstream& operator<<(std::ofstream& ofs, const Matrix<T>& data) {
	for (int i = 0; i < data.GetRow(); i++) {
		for (int j = 0; j < data.GetCol(); j++) {
			ofs << data.GetMtrData().at(i).at(j);
			j != data.GetCol() - 1 ? (ofs << '\t') : (ofs << '\n');
		}
	}
	return ofs;
}



template<class T>
DiagonalMatrix<T>::DiagonalMatrix(const Matrix<T>& inst) : Matrix<T>::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i == j ? this->mtr.at(i).at(j) = inst.mtr.at(i).at(j) : this->mtr.at(i).at(j) = 0;
}

template<class T>
UpperTriangularMatrix<T>::UpperTriangularMatrix(const Matrix<T>& inst) : Matrix<T>::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i >= j ? this->mtr.at(i).at(j) = inst.mtr.at(i).at(j) : this->mtr.at(i).at(j) = 0;
}

template<class T>
LowerTriangularMatrix<T>::LowerTriangularMatrix(const Matrix<T>& inst) : Matrix<T>::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i <= j ? this->mtr.at(i).at(j) = inst.mtr.at(i).at(j) : this->mtr.at(i).at(j) = 0;
}

template<class T>
SymmetricalMatrix<T>::SymmetricalMatrix(const Matrix<T>& inst) : Matrix<T>::Matrix(inst) {
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			i > j ? this->mtr.at(i).at(j) = inst.mtr.at(j).at(i) : this->mtr.at(i).at(j) = inst.mtr.at(i).at(j);
}

template <class T>
std::ifstream& Matrix<T>::load_from(std::ifstream& file) {
	if (!file) {
		std::cout << "Cannot open file.";
		return file;
	}
	for (int i = 0; i < this->Row; i++)
		for (int j = 0; j < this->Col; j++)
			file.read(reinterpret_cast<char*>(&mtr[i][j]), sizeof T);

	return file;
}

template <class T>
std::ofstream& Matrix<T>::save_to(std::ofstream& bin) {
	try {
		if (!bin) throw "error";

		for (int i = 0; i < this->Row; i++)
			for (int j = 0; j < this->Col; j++)
				bin.write(reinterpret_cast<char *>(&mtr[i][j]), sizeof T);
	}
	catch (std::string) { std::cout << "File couldn't be opened!" << std::endl; }
	return bin;

}

inline std::ifstream& operator>>(std::ifstream& ifs, Matrix<double>& fill) {
	try {
		if (ifs.is_open()) {
			std::vector<std::vector<double>> numbers;

			std::string temp;
			int size = 0; int size1 = 0;
			int k = 0;
			while (std::getline(ifs, temp)) {
				std::istringstream buffer(temp);
				std::vector<double> line;
				for (std::string s; buffer >> s;) line.push_back(std::stod(s));
				size = line.size();
				if (size != size1 && k != 0) throw(std::runtime_error("Cannot build matrix!"));
				size1 = line.size();


				numbers.push_back(line);
				k++;
			}
			fill.SetRow(numbers.size());
			fill.SetCol(numbers[0].size());
			fill.SetMtr(numbers);
		}

		else { std::cout << "File wasn't opened!" << std::endl; }
	}
	catch (std::runtime_error&) { std::cout << "\nFix the text file contents and make it look like 2D matrix.\n"; };
	return ifs;
}