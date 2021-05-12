using System;
using System.Collections.Generic;
using System.Text;

namespace matematica2
{
    class seidel_method
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


            Matrix d = Matrix.ZeroMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                d.mat[i, i] = A.mat[i, i];
            }
            Matrix c = Matrix.Multiply(d.Invert(), b);





            Matrix mat_bl = new Matrix(n, n);
            Matrix mat_bdr = new Matrix(n, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    mat_bl.mat[i, j] = mat_b[i, j];
                }

                for (int j = i; j < n; j++)
                {
                    mat_bdr.mat[i, j] = mat_b[i, j];
                }

            }


            Matrix x = seidel(mat_bl, mat_bdr, c, Matrix.ZeroMatrix(n, 1), mat_b, 0, mat_l, mat_d, mat_r);

            Program.Show(x.mat, x.rows, x.cols, "x:");
        }


        public static Matrix seidel(Matrix bl, Matrix bdr, Matrix c, Matrix x, Matrix mat_b, int iterator, Matrix l, Matrix d, Matrix r)
        {
            iterator++;
            List<double> solution = new List<double>();
            if (Matrix.equal(x, Matrix.ZeroMatrix(x.rows, x.cols)))
            {
                x = c;
                Matrix ld = l + d;

            }
            for (int i = 0; i < mat_b.rows; i++)
            {
                double x2 = 0;
                for (int j = 0; j < mat_b.cols; j++)
                {
                    try
                    {
                        x2 = x2 + (solution[j] * mat_b[i, j]);
                    }
                    catch
                    {
                        x2 = x2 + x.mat[j, 0] * mat_b[i, j];
                    }
                }
                x2 = x2 + c.mat[i, 0];
                solution.Add(x2);

            }



            Matrix solution2 = new Matrix(solution.Count, 1);
            for (int i = 0; i < solution.Count; i++)
            {
                solution2.mat[i, 0] = solution[i];
            }

            if (Matrix.norma(solution2 - x) < Math.Pow(10, -6) * (1 - Matrix.norma(mat_b)) / Matrix.norma(bdr))
                return solution2;
            return seidel(bl, bdr, c, solution2, mat_b, iterator, l, d, r);

        }



    }
}
