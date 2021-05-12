using System;
using System.Collections.Generic;
using System.Text;
using matematica2;

namespace matematica2
{
    class jacobi_method
    {
        public static void Solution()
        {
            int n = 4;
            Matrix A = Matrix.RandomMatrix(n, n, 10);
            Matrix b = Matrix.RandomMatrix(n, 1, 10);
            double[,] numsA = new double[,] { { 88, 5, 5, 7 }, { 1, 77, 2, 5 }, { 8, 1, 88, 2 }, { 1, 2, 3, 88 } };
            double[,] numsB = new double[,] { { 7 }, { 1 }, { 3 }, { 2 } };
            A.mat = numsA;
            b.mat = numsB;



            Program.Show(A.mat, "A:");
            Program.Show(b.mat, n, 1, "b:");


            Matrix d = Matrix.ZeroMatrix(n, n);


            for (int i = 0; i < n; i++)
            {
                d.mat[i, i] = A.mat[i, i];
            }
            Matrix c = Matrix.Multiply(d.Invert(), b);
            Matrix mat_d = new Matrix(n, n);
            Matrix mat_r = new Matrix(n, n);
            Matrix mat_l = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    mat_d[i, i] = A[i, i];
                    for (int j = i + 1; j < n; j++)
                    {
                        mat_r[i, j] = A[i, j];
                    }
                }
                else
                {
                    for (int j = 0; j < i; j++)
                    {
                        mat_l[i, j] = A[i, j];
                    }

                    mat_d[i, i] = A[i, i];

                    for (int j = i + 1; j < n; j++)
                    {
                        mat_r[i, j] = A[i, j];
                    }

                }
            }


            Matrix mat_b = -Matrix.Multiply(mat_d.Invert(), (mat_l + mat_r));



            Matrix x = jacobi(A, d, c, b, Matrix.ZeroMatrix(1, n), 0, mat_b);
            Program.Show(x.mat, x.rows, x.cols, "x:");


        }



        public static Matrix jacobi(Matrix a, Matrix d, Matrix c, Matrix b, Matrix x, int iterator, Matrix mat_b)
        {
            iterator++;
            Matrix x2 = new Matrix(x.rows, x.cols);
            if (Matrix.equal(x, Matrix.ZeroMatrix(x.rows, x.cols)))
            {

                x2 = Matrix.Multiply(d.Invert(), b);
                return jacobi(a, d, c, b, x2, iterator, mat_b);
            }
            else
            {

                x2 = Matrix.Multiply(mat_b, x) + c;
            }

            if ((Matrix.norma(x2 - x)) < Math.Pow(10, -6) * ((1 - Matrix.norma(mat_b)) / Matrix.norma(mat_b)))
                return x2;
            return jacobi(a, d, c, b, x2, iterator, mat_b);

        }





    }
}
