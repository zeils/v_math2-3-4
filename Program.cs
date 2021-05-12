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
            Show(Matrix.Multiply(A.L, A.U).mat, "LU=");
            Console.WriteLine("Det(LU)={0}", Matrix.Multiply(A.L, A.U).Det());
            Show(Matrix.Multiply(A.GetP(), A).mat, "AP=");
            Console.WriteLine("Det(PA)={0}", Matrix.Multiply(A.GetP(), A).Det());


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
            Show(Matrix.Multiply(A.Q, A.R).mat, "QR=");
            Show(A.R.SolveWith(Matrix.Multiply(A.Q.Transp(), b)).mat, n, 1, "x:");



        }

        public static void fourv1()
        {
            jacobi_method.Solution();
        }


        public static void fourv2()
        {
            seidel_method.Solution();
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

            one();
            two();
            three();
            fourv1();
            fourv2();

            Show(newton_method.Solution(0).mat, 10, 1, "Решение методом Ньютона:");
            


        }
    }
}
