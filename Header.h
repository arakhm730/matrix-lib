#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <iterator>
#include <sstream>

template<class T> class Matrix;
template <class T> class DiagonalMatrix;
template <class T> class UpperTriangularMatrix;
template <class T> class LowerTriangularMatrix;
template <class T> class SymmetricalMatrix;
template<class T> Matrix<T> operator+(Matrix<T> &, Matrix<T> &);
template<class T> Matrix<T> operator-(Matrix<T> &, Matrix<T> &);
template<class T> Matrix<T> AdamarMult(Matrix<T>&, Matrix<T>&);
template<class T> Matrix<T> operator*(const T&, Matrix<T> &);
template<class T> Matrix<T> operator*(Matrix<T> &, const T&);
template<class T> Matrix<T> operator*(Matrix<T> &, Matrix<T> &);
template<class T> std::ostream& operator<< (std::ostream&, const Matrix<T>&);
template<class T> std::istream& operator>> (std::istream&, Matrix<T>&);

template<class T> std::ofstream& operator<< (std::ofstream&, const Matrix<T>&);

inline std::ifstream& operator>>(std::ifstream&, Matrix<double>&);

template<class T>
class Matrix
{
private:
	int Row, Col;
	std::vector <std::vector <T>> mtr;



	double getDeterminant(const std::vector<std::vector<double>> vect);
	std::vector<std::vector<double>> getTranspose(const std::vector<std::vector<double>> matrix1);
	std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>> vect);
	std::vector<std::vector<double>> getInverse(const std::vector<std::vector<double>> vect);
	int compute_rank(std::vector<std::vector<T>> A);

public:
	Matrix() = default;
	Matrix(std::vector <std::vector <T>> data);
	Matrix(const int&, const int&);
	virtual ~Matrix() = default;
	Matrix(const Matrix<T> &);
	std::vector<std::vector<T>> GetMtrData() const;
	int GetRow() const;
	int GetCol() const;

	void SetRow(const int&);
	void SetCol(const int&);
	void SetMtr(std::vector<std::vector<T>>);


	double determinant();
	double trace() const;
	double vecNormMax() const;
	double vecNormEuclid() const;
	double matrixNorm() const;

	Matrix<T> transpose();
	int rank();
	Matrix<double> inverse();

	std::ifstream& load_from(std::ifstream&);
	std::ofstream& save_to(std::ofstream&);

	friend class IdentityMatrix;
	friend class DiagonalMatrix<T>;
	friend class UpperTriangularMatrix<T>;
	friend class LowerTriangularMatrix<T>;
	friend class SymmetricalMatrix<T>;

	friend Matrix<T> AdamarMult<>(Matrix&, Matrix&);

	friend Matrix<T> operator+ <>(Matrix &, Matrix &);
	friend Matrix<T> operator- <>(Matrix &, Matrix &);

	friend Matrix<T> operator* <>(const T&, Matrix &);
	friend Matrix<T> operator* <>(Matrix &, const T&);
	friend Matrix<T> operator* <>(Matrix &, Matrix &);

	friend std::ofstream& operator<< <>(std::ofstream&, const Matrix&);
	friend std::ostream& operator<< <>(std::ostream&, const Matrix&);
	friend std::istream& operator>> <>(std::istream&, Matrix&);

};


class IdentityMatrix: public Matrix<int> {
public:

	IdentityMatrix(const int& dim) : Matrix<int>::Matrix(dim, dim) {
		for (int i = 0; i < Row; i++) for (int j = 0; j < Col; j++)
			i == j ? this->mtr.at(i).at(j) = 1 : this->mtr.at(i).at(j) = 0;
	}
	friend Matrix<int> AdamarMult<>(Matrix&, Matrix&);

	friend Matrix<int> operator+ <>(Matrix &, Matrix &);
	friend Matrix<int> operator- <>(Matrix &, Matrix &);

	friend Matrix<int> operator* <>(const int&, Matrix &);
	friend Matrix<int> operator* <>(Matrix &, const int&);
	friend Matrix<int> operator* <>(Matrix &, Matrix &);
	friend std::ofstream& operator<< <>(std::ofstream&, const Matrix&);
	friend std::ostream& operator<< <>(std::ostream&, const Matrix&);
	friend std::istream& operator>> <>(std::istream&, Matrix&);
};

template<class T>
class DiagonalMatrix: public Matrix<T> {
public:
	
	DiagonalMatrix(const Matrix<T>&);

	friend Matrix<T> AdamarMult<>(Matrix<T>&, Matrix<T>&);

	friend Matrix<T> operator+ <>(Matrix<T> &, Matrix<T>&);
	friend Matrix<T> operator- <>(Matrix<T> &, Matrix<T> &);

	friend Matrix<T> operator* <>(const T&, Matrix<T> &);
	friend Matrix<T> operator* <>(Matrix<T> &, const T&);
	friend Matrix<T> operator* <>(Matrix<T> &, Matrix<T> &);

	friend std::ofstream& operator<< <>(std::ofstream&, const Matrix<T>&);
	friend std::ostream& operator<< <>(std::ostream&, const Matrix<T>&);
	friend std::istream& operator>> <>(std::istream&, Matrix<T>&);
};

template<class T>
class UpperTriangularMatrix : public Matrix<T> {
public:

	UpperTriangularMatrix(const Matrix<T>&);

	friend Matrix<T> AdamarMult<>(Matrix<T>&, Matrix<T>&);

	friend Matrix<T> operator+ <>(Matrix<T> &, Matrix<T>&);
	friend Matrix<T> operator- <>(Matrix<T> &, Matrix<T> &);

	friend Matrix<T> operator* <>(const T&, Matrix<T> &);
	friend Matrix<T> operator* <>(Matrix<T> &, const T&);
	friend Matrix<T> operator* <>(Matrix<T> &, Matrix<T> &);

	friend std::ofstream& operator<< <>(std::ofstream&, const Matrix<T>&);
	friend std::ostream& operator<< <>(std::ostream&, const Matrix<T>&);
	friend std::istream& operator>> <>(std::istream&, Matrix<T>&);
};

template<class T>
class LowerTriangularMatrix : public Matrix<T> {
public:

	LowerTriangularMatrix(const Matrix<T>&);

	friend Matrix<T> AdamarMult<>(Matrix<T>&, Matrix<T>&);

	friend Matrix<T> operator+ <>(Matrix<T> &, Matrix<T>&);
	friend Matrix<T> operator- <>(Matrix<T> &, Matrix<T> &);

	friend Matrix<T> operator* <>(const T&, Matrix<T> &);
	friend Matrix<T> operator* <>(Matrix<T> &, const T&);
	friend Matrix<T> operator* <>(Matrix<T> &, Matrix<T> &);

	friend std::ofstream& operator<< <>(std::ofstream&, const Matrix<T>&);
	friend std::ostream& operator<< <>(std::ostream&, const Matrix<T>&);
	friend std::istream& operator>> <>(std::istream&, Matrix<T>&);
};

template<class T>
class SymmetricalMatrix : public Matrix<T> {
public:
	SymmetricalMatrix(const Matrix<T>&);

	friend Matrix<T> AdamarMult<>(Matrix<T>&, Matrix<T>&);

	friend Matrix<T> operator+ <>(Matrix<T> &, Matrix<T>&);
	friend Matrix<T> operator- <>(Matrix<T> &, Matrix<T> &);

	friend Matrix<T> operator* <>(const T&, Matrix<T> &);
	friend Matrix<T> operator* <>(Matrix<T> &, const T&);
	friend Matrix<T> operator* <>(Matrix<T> &, Matrix<T> &);

	friend std::ofstream& operator<< <>(std::ofstream&, const Matrix<T>&);
	friend std::ostream& operator<< <>(std::ostream&, const Matrix<T>&);
	friend std::istream& operator>> <>(std::istream&, Matrix<T>&);
};