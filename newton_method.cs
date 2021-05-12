using System;
using System.Collections.Generic;
using System.Text;
using matematica2;




namespace matematica2
{
     

    class Function
    {
        private double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10; // элементы матрицы x
        public Matrix F = new Matrix(10,1);
        

        private static double sin(double x)
        {
            return Math.Sin(x);
        }
        private static double cos(double x)
        {
            return Math.Cos(x);
        }
        private static double pow(double x, double n)
        {
            return Math.Pow(x, n);
        }
        private static double exp(double x)
        {
            return Math.Exp(x);
        }

        public Function(Matrix x)
        {
            F_invalidate(x);
        }

        public void F_invalidate(Matrix x)
        {
            x1 = x.mat[0, 0];
            x2 = x.mat[1, 0];
            x3 = x.mat[2, 0];
            x4 = x.mat[3, 0];
            x5 = x.mat[4, 0];
            x6 = x.mat[5, 0];
            x7 = x.mat[6, 0];
            x8 = x.mat[7, 0];
            x9 = x.mat[8, 0];
            x10 = x.mat[9, 0];


            F[0, 0] = cos(x2 * x1) - exp(-(3.0 * x3)) + x4 * x5 * x5 - x6 - Math.Sinh((2.0 * x8)) * x9 + (2.0 * x10) + 2.000433974165385440;
            F[1, 0] = sin(x2 * x1) + x3 * x9 * x7 - exp(-x10 + x6) + 3.0 * x5 * x5 - x6 * (x8 + 1.0) + 10.886272036407019994;
            F[2, 0] = x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8 + x9 - x10 - 3.1361904761904761904;
            F[3, 0] = 2.0 * cos(-x9 + x4) + x5 / (x3 + x1) - sin(x2 * x2) + pow(cos(x7 * x10), 2.0) - x8 - 0.1707472705022304757;
            F[4, 0] = sin(x5) + 2.0 * x8 * (x3 + x1) - exp(-x7 * (-x10 + x6)) + 2.0 * cos(x2) - 1.0 / (-x9 + x4) - 0.3685896273101277862;
            F[5, 0] = exp(x1 - x4 - x9) + x5 * x5 / x8 + cos(3.0 * x10 * x2) / 2.0 - x6 * x3 + 2.0491086016771875115;
            F[6, 0] = pow(x2, 3.0) * x7 - sin(x10 / x5 + x8) + (x1 - x6) * cos(x4) + x3 - 0.7380430076202798014;
            F[7, 0] = x5 * pow(x1 - 2.0 * x6, 2.0) - 2.0 * sin(-x9 + x3) + 1.5 * x4 - exp(x2 * x7 + x10) + 3.5668321989693809040;
            F[8, 0] = 7.0 / x6 + exp(x5 + x4) - 2.0 * x2 * x8 * x10 * x7 + 3.0 * x9 - 3.0 * x1 - 8.4394734508383257499;
            F[9, 0] = x10 * x1 + x9 * x2 - (x8 * x3) + sin(x4 + x5 + x6) * x7 - 0.78238095238095238096;




        }




    }

    class Jacobi_mat
    {
        private double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10; // элементы матрицы x
        public Matrix J = new Matrix(10, 10);

        private static double sin(double x)
        {
            return Math.Sin(x);
        }
        private static double cos(double x)
        {
            return Math.Cos(x);
        }
        private static double pow(double x, double n)
        {
            return Math.Pow(x,n);
        }
        private static double exp(double x)
        {
            return Math.Exp(x);
        }


        public Jacobi_mat(Matrix x)
        {
            this.jacobian_invalidate(x);
        }

        public void jacobian_invalidate(Matrix x)
        {
            x1 = x.mat[0, 0];
            x2 = x.mat[1, 0];
            x3 = x.mat[2, 0];
            x4 = x.mat[3, 0];
            x5 = x.mat[4, 0];
            x6 = x.mat[5, 0];
            x7 = x.mat[6, 0];
            x8 = x.mat[7, 0];
            x9 = x.mat[8, 0];
            x10 = x.mat[9, 0];
            //double[,] J = new double[10, 10];


            J[0, 0] = -x2 * sin(x2 * x1);
            J[0, 1] = -x1 * sin(x2 * x1);
            J[0, 2] = 3.0 * exp(-(3.0 * x3));
            J[0, 3] = x5 * x5;
            J[0, 4] = 2.0 * x4 * x5;
            J[0, 5] = -1.0;
            J[0, 6] = 0.0;
            J[0, 7] = -2.0 * Math.Cosh(2.0 * x8) * x9;
            J[0, 8] = -Math.Sinh(2.0 * x8);
            J[0, 9] = 2.0;
            J[1, 0] = x2 * cos(x2 * x1);
            J[1, 1] = x1 * cos(x2 * x1);
            J[1, 2] = x9 * x7;
            J[1, 3] = 0.0;
            J[1, 4] = 6.0 * x5;
            J[1, 5] = -exp(-x10 + x6) - x8 - 1.0;
            J[1, 6] = x3 * x9;
            J[1, 7] = -x6;
            J[1, 8] = x3 * x7;
            J[1, 9] = exp(-x10 + x6);
            J[2, 0] = 1;
            J[2, 1] = -1;
            J[2, 2] = 1;
            J[2, 3] = -1;
            J[2, 4] = 1;
            J[2, 5] = -1;
            J[2, 6] = 1;
            J[2, 7] = -1;
            J[2, 8] = 1;
            J[2, 9] = -1;
            J[3, 0] = -x5 * pow(x3 + x1, -2.0);
            J[3, 1] = -2.0 * x2 * cos(x2 * x2);
            J[3, 2] = -x5 * pow(x3 + x1, -2.0);
            J[3, 3] = -2.0 * sin(-x9 + x4);
            J[3, 4] = 1.0 / (x3 + x1);
            J[3, 5] = 0;
            J[3, 6] = -2.0 * cos(x7 * x10) * x10 * sin(x7 * x10);
            J[3, 7] = -1;
            J[3, 8] = 2.0 * sin(-x9 + x4);
            J[3, 9] = -2.0 * cos(x7 * x10) * x7 * sin(x7 * x10);
            J[4, 0] = 2 * x8;
            J[4, 1] = -2.0 * sin(x2);
            J[4, 2] = 2 * x8;
            J[4, 3] = pow(-x9 + x4, -2.0);
            J[4, 4] = cos(x5);
            J[4, 5] = x7 * exp(-x7 * (-x10 + x6));
            J[4, 6] = -(x10 - x6) * exp(-x7 * (-x10 + x6));
            J[4, 7] = (2 * x3) + 2.0 * x1;
            J[4, 8] = -pow(-x9 + x4, -2.0);
            J[4, 9] = -x7 * exp(-x7 * (-x10 + x6));
            J[5, 0] = exp(x1 - x4 - x9);
            J[5, 1] = -3.0 / 2.0 * x10 * sin(3.0 * x10 * x2);
            J[5, 2] = -x6;
            J[5, 3] = -exp(x1 - x4 - x9);
            J[5, 4] = 2 * x5 / x8;
            J[5, 5] = -x3;
            J[5, 6] = 0;
            J[5, 7] = -x5 * x5 * pow(x8, (-2));
            J[5, 8] = -exp(x1 - x4 - x9);
            J[5, 9] = -3.0 / 2.0 * x2 * sin(3.0 * x10 * x2);
            J[6, 0] = cos(x4);
            J[6, 1] = 3.0 * x2 * x2 * x7;
            J[6, 2] = 1;
            J[6, 3] = -(x1 - x6) * sin(x4);
            J[6, 4] = x10 * pow(x5, (-2)) * cos(x10 / x5 + x8);
            J[6, 5] = -cos(x4);
            J[6, 6] = pow(x2, 3.0);
            J[6, 7] = -cos(x10 / x5 + x8);
            J[6, 8] = 0;
            J[6, 9] = -1.0 / x5 * cos(x10 / x5 + x8);
            J[7, 0] = 2.0 * x5 * (x1 - 2.0 * x6);
            J[7, 1] = -x7 * exp(x2 * x7 + x10);
            J[7, 2] = -2.0 * cos(-x9 + x3);
            J[7, 3] = 1.5;
            J[7, 4] = pow(x1 - 2.0 * x6, 2.0);
            J[7, 5] = -4.0 * x5 * (x1 - 2.0 * x6);
            J[7, 6] = -x2 * exp(x2 * x7 + x10);
            J[7, 7] = 0;
            J[7, 8] = 2.0 * cos(-x9 + x3);
            J[7, 9] = -exp(x2 * x7 + x10);
            J[8, 0] = -3;
            J[8, 1] = -2.0 * x8 * x10 * x7;
            J[8, 2] = 0;
            J[8, 3] = exp((x5 + x4));
            J[8, 4] = exp((x5 + x4));
            J[8, 5] = -7.0 * pow(x6, -2.0);
            J[8, 6] = -2.0 * x2 * x8 * x10;
            J[8, 7] = -2.0 * x2 * x10 * x7;
            J[8, 8] = 3;
            J[8, 9] = -2.0 * x2 * x8 * x7;
            J[9, 0] = x10;
            J[9, 1] = x9;
            J[9, 2] = -x8;
            J[9, 3] = cos(x4 + x5 + x6) * x7;
            J[9, 4] = cos(x4 + x5 + x6) * x7;
            J[9, 5] = cos(x4 + x5 + x6) * x7;
            J[9, 6] = sin(x4 + x5 + x6);
            J[9, 7] = -x3;
            J[9, 8] = x2;
            J[9, 9] = x1;



        }


    }



    class newton_method
    {

        public static double eps = Math.Pow(10, -6);

        public static Matrix Solution(int k)
        {
            Matrix x = Matrix.ZeroMatrix(10, 1);
            int iterator = 0;
            x.mat[0, 0] = 0.5;
            x.mat[1, 0] = 0.5;
            x.mat[2, 0] = 1.5;
            x.mat[3, 0] = -1.0;
            x.mat[4, 0] = -0.5;
            x.mat[5, 0] = 1.5;
            x.mat[6, 0] = 0.5;
            x.mat[7, 0] = -0.5;
            x.mat[8, 0] = 1.5;
            x.mat[9, 0] = -1.5;
            return Calculate(x, ref iterator, x, k);
            
        }

        public static Matrix Calculate(Matrix x, ref int iterator, Matrix x0, int k)
        {
            iterator++;
            Matrix nextX = Matrix.ZeroMatrix(10, 1);
            Matrix deltaX = Matrix.ZeroMatrix(10, 1);
            if (iterator <= k)
            {
                x0 = x;
            }
                Jacobi_mat j = new Jacobi_mat(x0);
            // k = 1
            // iter 1 , 

            
             // v1
            Function f = new Function(x);
            deltaX = j.J.SolveWith(-f.F);
            nextX = x + deltaX; // x(k+1) = deltaX + x(k)



            // x(k+1) = x(k) + J^-1(x0)*F(x(k))
            // x(k+1) - x(k) = delta(x)
            // delta(x) = J^-1(x0)*F(x(k))
            // J(x0) *delta(x) = F(x(k))


            
            Console.WriteLine("--------------------------------------------------------------------------}");
            Console.WriteLine("iterator = {0}", iterator);
            Program.Show(x.mat, 10, 1, "x=");
            Program.Show(deltaX.mat, 10, 1, "Delta");
            Console.WriteLine("--------------------------------------------------------------------------}");


            if (Check(x, nextX))
            {
               
                return nextX;
            }
            return Calculate(nextX, ref iterator,x0,k);
















        }

        private static bool Check(Matrix currentX, Matrix nextX)
        {
           
            double max = 0;
            for (int i = 0; i<currentX.rows; i++)
            {
                if (max < Math.Abs(currentX.mat[i, 0] - nextX.mat[i, 0])) max = Math.Abs(currentX.mat[i, 0] - nextX.mat[i, 0]);

            }
            if (eps > max)
            {
                return true;
            }
            else
            {
                return false;
            }

        }

       


       
        

        
        

      
        
        
    





    }
}
