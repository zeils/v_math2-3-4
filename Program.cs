using System;
namespace matematica2
{
    class Program
    {

        static int n = 4;

        public static void one()
        {
            Matrix A = Matrix.RandomMatrix(n, n, 10);
            Matrix b = Matrix.RandomMatrix(n, 1, 10);
            Show(A.mat, "А:");
            Console.WriteLine("DetA={0}", A.Det());
            A.MakeLU();
            Show(A.L.mat, "L:");
            Show(A.U.mat, "U:");
            Show(Matrix.StupidMultiply(A.L, A.U).mat, "LU=");
            Console.WriteLine("Det(LU)={0}", Matrix.StupidMultiply(A.L, A.U).Det());
            Show(Matrix.StupidMultiply(A.GetP(), A).mat, "AP=");
            Console.WriteLine("Det(PA)={0}", Matrix.StupidMultiply(A.GetP(), A).Det());


            Console.WriteLine("--------------------------------------------------------------------");
            Console.WriteLine("det(A) = {0}, det (LUP) = {1}", A.Det(), A.L.Det() * A.U.Det() * A.GetP().Det());
            Console.WriteLine("--------------------------------------------------------------------");
            Show(b.mat, n, 1, "b:");

            Show(A.SolveWith(b).mat, n, 1, "solved Ax=b:");
            Console.WriteLine("--------------------------------------------------------------------");
            Show(A.Invert().mat, "inverted A");
            Console.WriteLine("Cond(A) = det(inv(A))* det(A) = {0}", A.Invert().Det() * A.Det());
            Console.WriteLine("--------------------------------------------------------------------");






        }

        public static void two()
        {
            Matrix A = Matrix.RandomMatrix(n, n, 10);
            Show(A.mat, "А:");
            Console.WriteLine("Rank =  {0}", A.Rank());


        }

        public static void three()
        {
            Matrix A = Matrix.RandomMatrix(n, n, 10);
            Matrix b = Matrix.RandomMatrix(n, 1, 10);
            Show(A.mat, "A:");
            A.QR_solution();
            Show(Matrix.StupidMultiply(A.Q, A.R).mat, "QR=");
            Show(A.R.SolveWith(Matrix.StupidMultiply(A.Q.Transp(), b)).mat, n, 1, "x:");



        }

        public static void fourv1()
        {
            Matrix A = Matrix.RandomMatrix(n, n, 10);
            Matrix b = Matrix.RandomMatrix(n, 1, 10);
            double[,] numsA = new double[,] { { 88, 5, 5, 7 }, { 1, 77, 2, 5 }, { 8, 1, 88, 2 }, { 1, 2, 3, 88 } };
            double[,] numsB = new double[,] { { 7 }, { 1 }, { 3 }, { 2 } };
            A.mat = numsA;
            b.mat = numsB;



            Show(A.mat, "A:");
            Show(b.mat, n, 1, "b:");


            Matrix d = Matrix.ZeroMatrix(n, n);


            for (int i = 0; i < n; i++)
            {
                d.mat[i, i] = A.mat[i, i];
            }
            Matrix c = Matrix.StupidMultiply(d.Invert(), b);
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


            Matrix mat_b = -Matrix.StupidMultiply(mat_d.Invert(), (mat_l + mat_r));



            Matrix x = Matrix.jacobi(A, d, c, b, Matrix.ZeroMatrix(1, n), 0, mat_b);
            Show(x.mat, x.rows, x.cols, "x:");



        }


        public static void fourv2()
        {
            Matrix A = Matrix.RandomMatrix(n, n, 10);
            Matrix b = Matrix.RandomMatrix(n, 1, 10);
            double[,] numsA = new double[,] { { 88, 5, 5, 7 }, { 1, 77, 2, 5 }, { 8, 1, 88, 2 }, {1, 2, 3, 88 } };
            double[,] numsB = new double[,] { { 7 }, { 1 }, { 3 }, { 2 } };
            A.mat = numsA;
            b.mat = numsB;


            Show(A.mat, "A:");
            Show(b.mat, n, 1, "b:");

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


            Matrix mat_b = -Matrix.StupidMultiply(mat_d.Invert(), (mat_l + mat_r));


            Matrix d = Matrix.ZeroMatrix(n, n);
            for (int i = 0; i < n; i++)
            {
                d.mat[i, i] = A.mat[i, i];
            }
            Matrix c = Matrix.StupidMultiply(d.Invert(), b);





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


            Matrix x = Matrix.seidel(mat_bl, mat_bdr, c, Matrix.ZeroMatrix(n, 1), mat_b, 0, mat_l, mat_d, mat_r);

            Show(x.mat, x.rows, x.cols, "x:");










        }

        public static void Show(double[,] arr, string name)
        {
            Console.WriteLine(name);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(" {0} ", arr[i, j]);

                }
                Console.WriteLine(" ");
            }
            Console.WriteLine(" ");



        }
        public static void Show(double[,] arr, int col, int row, string name)
        {
            Console.WriteLine(name);
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    Console.Write(" {0} ", arr[j, i]);

                }
                Console.WriteLine(" ");
            }
            Console.WriteLine(" ");



        }

        public static void Show(double[] arr, int col, string name)
        {
            Console.WriteLine(name);
            for (int i = 0; i < col; i++)
            {
                Console.WriteLine(arr[i]);
            }
            Console.WriteLine(" ");

        }











        static void Main(string[] args)
        {

            //one();
            //two();
            //three();
            //fourv1();
            //fourv2();
            
            Show(newton_method.Solution(0).mat, 10, 1, "Решение методом Ньютона:");
            


        }
    }
}
