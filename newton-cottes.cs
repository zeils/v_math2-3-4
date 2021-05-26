using System;
using System.Collections.Generic;
using System.Text;

namespace vMath2
{

    static class z
    {
        private static double[] z_mass;

        private static int n = 3; // трехточечное

        
        private static double a = 1.5;
        private static double b = 2.3;
        private static double alpha = 0.2;
        private static double beta = 0;

       

        //z[1] 
        //z[1.5]
        //z[2]

        // *2
        //z[2] 
        //z[3]
        //z[4]
        public static void Calculate(int k)
        {
            z_mass = new double[2*(k+1)];
            double h = (b - a) / k;
            for (int i = 0; i <= k; i++)
            {
                // i = 2i
                z_mass[2 * i] = a + i * h;

            }

            for (int i = 0; i <= k; i++)
            {
                z_mass[2 * i + 1] = a + (i + 0.5) * h;
                
            }
        }

        // i 1
        // i 0.5


        public static double Get(double i)
        {
            int k = Convert.ToInt32(i * 2);
            return z_mass[k];

        }
                
    } // complete





    static class newton_cottes
    {


        private static int n = 3; // трехточечное

        private static double a = 1.5;
        private static double b = 2.3;
        private static double alpha = 0.2;
        private static double beta = 0;

        private static double answer0 = 32.21951452884234295708696008290380201405;

        //private static List<List<double>> A = new List<List<double>>() ;

        //private static List<List<double>> Nu = new List<List<double>>();

        private static double[,] A;
        private static double[,] Nu;

        public static void Solution()
        {
            double eps = Math.Pow(10, -6);
            int k = 2;
            double answer = CalculateI(k);
           
            while(CalculateI(k) - CalculateI(k-1) > eps)
            {
                k++;
                answer = CalculateI(k);
            }



            Console.WriteLine("answer = {0}", answer);
            Console.WriteLine("answer0 = {0}", answer0);
            Console.WriteLine("delta_answer = {0}", answer - answer0);



        }

         private static double CalculateI(int k)
        {
            z.Calculate (k);

            A = new double[k+1,4];
            Nu = new double[k+1, 4];

            CalculateNu(k);
            CalculateA(k);
            return Calculation(k);

            


            




        }

        private static double F (double x)
        {
            double f;
            f = 2 * Math.Cos(3.5 * x) * Math.Exp(5 * x / 3) + 3 * Math.Sin(1.5 * x) * Math.Exp(-4 * x) + 3;
            return f;

        }


        private static double Calculation( int k ) // complete
        {
            double I = A[1, 1]* F(z.Get(0)) + A[k,3]*F(z.Get(k)) ;
            for (int i = 1; i<= k; i++)
            {
                I = I + A[i, 2] * F(z.Get(i - (1 / 2)));
            }

            for (int i = 1; i <= k -1; i++)
            {
                I = I + (A[i, 3]+A[i+1, 1])*F(z.Get(i));
            }


            return I;
        }

        private static void CalculateNu(int k)
        {
            
            // n =3
            for(int i = 1; i <= k; i++)
            {
                
                double t = (Math.Pow(z.Get(i) - a, 1 - alpha) - Math.Pow(z.Get(i - 1) - a, 1 - alpha)) / (1 - alpha);
                Nu[i, 0] = (Math.Pow(z.Get(i)- a, 1 - alpha) - Math.Pow(z.Get(i - 1) - a, 1 - alpha)) / (1 - alpha);
                Nu[i, 1] = (Math.Pow(z.Get(i) - a, 2 - alpha) - Math.Pow(z.Get(i - 1) - a, 2 - alpha)) / (2 - alpha) + a * Nu[i, 0];
                Nu[i, 2] = (Math.Pow(z.Get(i) - a, 3 - alpha) - Math.Pow(z.Get(i - 1) - a, 3 - alpha)) / (3 - alpha) + 2 * a * Nu[i, 1] - a * a * Nu[i, 0];
                

            }    

           
           
            

        }

        private static void CalculateA(int k)
        {
            for (int i = 1; i <= k; i++)
            {

                A[i, 1] = (Nu[i, 2] - Nu[i, 1]*(z.Get(i- 0.5) +z.Get(i)) + Nu[i, 0]*z.Get(i- 0.5) *z.Get(i)) / ((z.Get(i- 0.5) -z.Get(i-1))*(z.Get(i)-z.Get(i-1)));

                A[i, 2] = -(Nu[i, 2] - Nu[i, 1] * (z.Get(i - 1) + z.Get(i)) + Nu[i, 0] * z.Get(i - 1) * z.Get(i)) / ((z.Get(i - 0.5) - z.Get(i - 1)) * (z.Get(i) - z.Get(i - 0.5)));

                A[i, 3] = (Nu[i, 2] - Nu[i, 1] * (z.Get(i - 0.5) + z.Get(i-1)) + Nu[i, 0] * z.Get(i - 0.5) * z.Get(i-1)) / ((z.Get(i) - z.Get(i - 0.5)) * (z.Get(i) - z.Get(i - 1)));
            }
        }

    }
}
