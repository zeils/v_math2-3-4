using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace matematica2
{
    public class Matrix
    {
        public int rows;
        public int cols;
        public double[,] mat;

        public Matrix L;
        public Matrix U;

        public Matrix Q;
        public Matrix R;

        public int[] pi;
        private double detType = 1;

        public Matrix(int iRows, int iCols)
        {
            rows = iRows;
            cols = iCols;
            mat = new double[rows, cols];
        }



        public double this[int st, int col]
        {
            get { return mat[st, col]; }
            set { mat[st, col] = value; }
        }

        public static Matrix IdentityMatrix(int iRows, int iCols)
        {
            Matrix matrix = ZeroMatrix(iRows, iCols);
            for (int i = 0; i < Math.Min(iRows, iCols); i++)
                matrix[i, i] = 1;
            return matrix;
        }



        public void MakeLU()
        {

            L = IdentityMatrix(rows, cols);
            U = Duplicate();

            pi = new int[rows];
            for (int i = 0; i < rows; i++) pi[i] = i;

            double p = 0;
            double pom2;
            int k0 = 0;
            int pom1 = 0;

            for (int k = 0; k < cols - 1; k++)
            {
                p = 0;
                for (int i = k; i < rows; i++)
                {
                    if (Math.Abs(U[i, k]) > p)
                    {
                        p = Math.Abs(U[i, k]);
                        k0 = i;
                    }
                }


                pom1 = pi[k]; pi[k] = pi[k0]; pi[k0] = pom1;

                for (int i = 0; i < k; i++)
                {
                    pom2 = L[k, i]; L[k, i] = L[k0, i]; L[k0, i] = pom2;
                }

                if (k != k0) detType *= -1;

                for (int i = 0; i < cols; i++)
                {
                    pom2 = U[k, i]; U[k, i] = U[k0, i]; U[k0, i] = pom2;
                }

                for (int i = k + 1; i < rows; i++)
                {
                    L[i, k] = U[i, k] / U[k, k];
                    for (int j = k; j < cols; j++)
                        U[i, j] = U[i, j] - L[i, k] * U[k, j];
                }
            }
        }


        public Matrix SolveWith(Matrix v)
        {

            if (L == null) MakeLU();

            Matrix b = new Matrix(rows, 1);
            for (int i = 0; i < rows; i++) b[i, 0] = v[pi[i], 0];

            Matrix z = SubsForth(L, b);
            Matrix x = SubsBack(U, z);

            return x;
        }

        public Matrix SolveWith(Matrix v, ref int iterator)
        {

            if (L == null) MakeLU();
            iterator = iterator = 3* this.cols * this.cols * this.cols;

            Matrix b = new Matrix(rows, 1);
            for (int i = 0; i < rows; i++)
            {
                iterator++;
                b[i, 0] = v[pi[i], 0];
            }

            Matrix z = SubsForth(L, b , ref iterator);
            Matrix x = SubsBack(U, z, ref iterator);

            return x;
        }

        public void SetCol(Matrix v, int k)
        {
            for (int i = 0; i < rows; i++) mat[i, k] = v[i, 0];
        }

        public Matrix Invert()
        {
            if (L == null) MakeLU();

            Matrix inv = new Matrix(rows, cols);

            for (int i = 0; i < rows; i++)
            {
                Matrix Ei = Matrix.ZeroMatrix(rows, 1);
                Ei[i, 0] = 1;
                Matrix col = SolveWith(Ei);
                inv.SetCol(col, i);
            }
            return inv;
        }


        public double Det()
        {
            if (L == null) MakeLU();
            double det = detType;
            for (int i = 0; i < rows; i++) det *= U[i, i];
            return det;
        }

        public Matrix GetP()
        {
            if (L == null) MakeLU();

            Matrix matrix = ZeroMatrix(rows, cols);
            for (int i = 0; i < rows; i++) matrix[pi[i], i] = 1;
            return matrix;
        }

        public Matrix Duplicate()
        {
            Matrix matrix = new Matrix(rows, cols);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    matrix[i, j] = this.mat[i, j];
            return matrix;
        }

        public static Matrix SubsForth(Matrix A, Matrix b)
        {
            if (A.L == null) A.MakeLU();
            int n = A.rows;
            Matrix x = new Matrix(n, 1);

            for (int i = 0; i < n; i++)
            {
                x[i, 0] = b[i, 0];
                for (int j = 0; j < i; j++) x[i, 0] -= A[i, j] * x[j, 0];
                x[i, 0] = x[i, 0] / A[i, i];
            }
            return x;
        }

        public static Matrix SubsForth(Matrix A, Matrix b, ref int iterator)
        {
            if (A.L == null) A.MakeLU();
            int n = A.rows;
            Matrix x = new Matrix(n, 1);

            for (int i = 0; i < n; i++)
            {
                x[i, 0] = b[i, 0];
                for (int j = 0; j < i; j++) x[i, 0] -= A[i, j] * x[j, 0];
                x[i, 0] = x[i, 0] / A[i, i];
                iterator = iterator + i;
            }
            return x;
        }

        public static Matrix SubsBack(Matrix A, Matrix b)
        {
            if (A.L == null) A.MakeLU();
            int n = A.rows;
            Matrix x = new Matrix(n, 1);

            for (int i = n - 1; i > -1; i--)
            {
                x[i, 0] = b[i, 0];
                for (int j = n - 1; j > i; j--) x[i, 0] -= A[i, j] * x[j, 0];
                x[i, 0] = x[i, 0] / A[i, i];
            }
            return x;
        }

        public static Matrix SubsBack(Matrix A, Matrix b, ref int iterator)
        {
            if (A.L == null) A.MakeLU();
            int n = A.rows;
            Matrix x = new Matrix(n, 1);

            for (int i = n - 1; i > -1; i--)
            {
                x[i, 0] = b[i, 0];
                for (int j = n - 1; j > i; j--) x[i, 0] -= A[i, j] * x[j, 0];
                x[i, 0] = x[i, 0] / A[i, i];
                iterator = iterator + i +1;
            }
            return x;
        }

        public static Matrix ZeroMatrix(int iRows, int iCols)
        {
            Matrix matrix = new Matrix(iRows, iCols);
            for (int i = 0; i < iRows; i++)
                for (int j = 0; j < iCols; j++)
                    matrix[i, j] = 0;
            return matrix;
        }



        public static Matrix RandomMatrix(int iRows, int iCols, int dispersion)
        {
            Random random = new Random();
            Matrix matrix = new Matrix(iRows, iCols);
            for (int i = 0; i < iRows; i++)
                for (int j = 0; j < iCols; j++)
                    matrix[i, j] = random.Next(-dispersion, dispersion);
            
            return matrix;
        }


        public static Matrix Multiply(Matrix m1, Matrix m2)
        {

            Matrix result = ZeroMatrix(m1.rows, m2.cols);
            for (int i = 0; i < result.rows; i++)
                for (int j = 0; j < result.cols; j++)
                    for (int k = 0; k < m1.cols; k++)
                        result[i, j] += m1[i, k] * m2[k, j];
            return result;
        }
        private static Matrix Multiply(double n, Matrix m)
        {
            Matrix r = new Matrix(m.rows, m.cols);
            for (int i = 0; i < m.rows; i++)
                for (int j = 0; j < m.cols; j++)
                    r[i, j] = m[i, j] * n;
            return r;
        }
        private static Matrix Add(Matrix m1, Matrix m2)
        {

            Matrix r = new Matrix(m1.rows, m1.cols);
            for (int i = 0; i < r.rows; i++)
                for (int j = 0; j < r.cols; j++)
                    r[i, j] = m1[i, j] + m2[i, j];
            return r;
        }




        public static Matrix Eye(int n)
        {
            Matrix t = new Matrix(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    if (i == j) t.mat[i, j] = 1;
                    else t.mat[i, j] = 0;

                }
            return t;
        }

        public Matrix Transp()
        {
            int n = rows;

            Matrix t = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    t[i, j] = this[j, i];
                }
            }
            return t;
        }

        public double Rank()
        {

            int n = rows;
            int i = 1;
            int j = 1;
            int f = 0;
            int rank = 0;
            while (i<n)
            {
                if (this.mat[i,j] != 0)
                {
                    i++;
                    rank = i - 1;
                }
                else
                {
                    if(f==0)
                    {
                        f = 1;
                        this.mat[j, n] = this.mat[n, j];
                    }
                    else
                    {
                        f = 0;
                        rank = i - 1;
                        j++;
                    }
                }

            }



            return rank;
        }

        public static void CheckRank()
        {
            Matrix nums = new Matrix(3, 3);
            
            nums[0, 0] = 1;
            nums[0, 1] = 2;
            nums[0, 2] = 3;

            nums[1, 0] = 3;
            nums[1, 1] = 2;
            nums[1, 2] = 1;

            nums[2, 0] = 4;
            nums[2, 1] = 4;
            nums[2, 2] = 4;

            Console.WriteLine(nums.Rank());






        }

        public static double norma(Matrix a)
        {
            double max = 0;
            for (int j = 0; j < a.rows; j++)
            {
                double sum = 0;
                for (int i = 0; i < a.cols; i++)
                {

                    sum = sum + Math.Abs(a.mat[j, i]);
                }
                if (sum > max) max = sum;
            }

            return max;
        }





        public static bool equal(Matrix m1, Matrix m2)
        {
            int c = m1.cols;
            int r = m1.rows;
            if ((r != m2.rows) || (c != m2.cols)) return false;
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    if (m1[i, j] != m2[i, j]) return false;
                }
            }



            return true;
        }












       


        





        public void QR_solution()
        {
            int n = rows;
            R = this;
            Q = Eye(n);

            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                {
                    Matrix q_temp = Eye(n);
                    double s = -R.mat[j, i] / Math.Sqrt(R.mat[i, i] * R.mat[i, i] + R.mat[j, i] * R.mat[j, i]);
                    double c = R.mat[i, i] / Math.Sqrt(R.mat[i, i] * R.mat[i, i] + R.mat[j, i] * R.mat[j, i]);
                    q_temp.mat[i, i] = c;
                    q_temp.mat[j, i] = s;
                    q_temp.mat[j, j] = c;
                    q_temp.mat[i, j] = -s;
                    R = Multiply(q_temp, R);
                    q_temp.mat[i, j] = s;
                    q_temp.mat[j, i] = -s;
                    Q = Multiply(Q, q_temp);

                }
        }





        public static Matrix operator -(Matrix m)
        { return Matrix.Multiply(-1, m); }

        public static Matrix operator +(Matrix m1, Matrix m2)
        { return Matrix.Add(m1, m2); }

        public static Matrix operator -(Matrix m1, Matrix m2)
        { return Matrix.Add(m1, -m2); }









    }
}