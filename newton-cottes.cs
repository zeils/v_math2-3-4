using System;
using System.Collections.Generic;
using System.Text;

namespace matematica2
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

        private static double eps = Math.Pow(10, -6);
        private static int n = 3; // трехточечное

        private static double a = 1.5;
        private static double b = 2.3;
        private static double alpha = 0.2;
        private static double beta = 0;

        //private static double a = 0.1;
        //private static double b = 2.3;
        //private static double alpha = 0.2;
        //private static double beta = 0;
        




        //private static double exact_value = 32.21951452884234295708696008290380201405;
        private static double exact_value = 3.578861536040539915439859609644293194417;


        private static double[,] A;
        private static double[,] Nu;

        public static void Solution_newton_cottes()
        {
            
            int k = 3;





            Console.WriteLine("Точное значение интеграла = {0}", exact_value);
            Formula_srednik_pg(k);
            Console.WriteLine("Интерполяционная формула Ньютона-Котса = {0}", newt_cots(a,b,k));
            Console.WriteLine("Погрешность = {0}", Math.Abs(newt_cots(a,b,k) - exact_value));
            Console.WriteLine("Методическая погрешность {0}", 6.17853);
            List<double> answ = skf_newton_kots(a,b,k);
            Console.WriteLine("Составная квадратурная формула Ньютона-Котса {0}", answ[0]);
            Console.WriteLine("Погрешность = {0}", Math.Abs(answ[0] - exact_value));

            Console.WriteLine("Составная квадратурная формула Ньютона-Котса с оптимальным шагом {0}" , skf_newton_kots(a, b, n,  Math.Round((b - a) / answ[1]))[0]);
            // Mn = 2262.14 // ограничение третьей производной
            // 2262.14/6*abs((t-1.5)*(t-2.3)*(t-(1.5+2.3)/2)/(t-1.5)^0.2) // 6.17853




        }

        

        private static double richardson(int r, List<double> s, double h, double L, double m)
        {
            r = r + 1;
            Matrix matrix = new Matrix(r, r);
            for (int i = 0; i< r; i++)
            {
                for (int j = 0; j < r-1; j++)
                {
                    matrix[i, j] = Math.Pow(h / Math.Pow(L, i), m + j);
                }
                matrix[i, r - 1] = -1;
            }
            Matrix b = new Matrix(r, 1);
            for (int i =0; i<r; i++)
            {
                b[i, 0] = -1 * s[i];
            }
            Matrix c = matrix.SolveWith(b);
            return c[r - 1,0];
        }

        private static double eitken(List<double> s, double L)
        {
            double m = -Math.Log(Math.Abs((s[s.Count - 1] - s[s.Count - 2]) / (s[s.Count - 2] - s[s.Count - 3]))) / Math.Log(L);
            Console.WriteLine("m= {0}", m);
            return m;
                
        }

        private static List<double> skf_newton_kots(double a, double b, int n, double k =1)
        {
            double eps2 = 1.5;
            double Rn = 1;
            int it = 1;
            double m = 0;
            int r = 0;
            double L = 2;
            double hopt = -1;
            List<double> s = new List<double>();
            s.Add(0);
            double hconst = (b - a) / k;
            double result = 0;

          
            while (Math.Abs(Rn) > eps)
            {
                double h = (b - a) / k;
                for (int i = 0; i < k; i++)
                    s[it - 1] = s[it - 1] + newt_cots(a + i * h, a + (i + 1) * h, n);
                if (it >= 3)
                {
                    m = eitken(s, L);
                }
                if (Math.Abs(m - 3) < eps2)
                {
                    if (r >= 2)
                    {
                        Rn = richardson(r, s, hconst, L, m) - s[s.Count - 1];
                        if ((r >= 2) && (hopt == -1))
                        {
                            hopt = h * (Math.Pow(1e-6 / Math.Abs(Rn), (1 / m))) * 1.95;
                        }
                    }
                }
                r++;
                k = k * L;
                it++;
                result = s[s.Count - 1];
                Console.WriteLine("k {0} - s {1}", k / L, s[s.Count - 1]);
                s.Add(0);
               
            }
            Console.WriteLine("kopt {0}", (b - a) / hopt);
            List<double> answer = new List<double>();
            answer.Add(s[s.Count - 2]);
            answer.Add(hopt);
            return answer;

        }













        private static void Formula_srednik_pg(int k)
        {
            double h = (b - a) / k;
            double s = 0;
            for (int i =1; i <=k;i++)
            {
                s = s + F(a + (i - 0.5) * h);
            }
            s = s * h;
            Console.WriteLine("Составная формула средник прямоугольников = {0} ", s);
            
        }



       

        private static double F (double x)
        {
            double f;
            f = 2 * Math.Cos(3.5 * x) * Math.Exp(5 * x / 3) + 3 * Math.Sin(1.5 * x) * Math.Exp(-4 * x) + 3;
            //f = 2.5 * Math.Cos(2 * x) * Math.Exp(2 * x / 3) + 4 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 3 * x;
            return f;

        }

        


       

        private static double newt_cots(double x0, double xn, int k)
        {
            // a,b,n



            Matrix matr = new Matrix(3, 3);
            matr[0, 0] = 1;
            matr[0, 1] = 1;
            matr[0, 2] = 1;

            matr[1, 0] = x0;
            matr[1, 1] = (x0 + xn) * 0.5;
            matr[1, 2] = xn;

            matr[2, 0] = Math.Pow(matr[1, 0], 2);
            matr[2, 1] = Math.Pow(matr[1, 1], 2);
            matr[2, 2] = Math.Pow(matr[1, 2], 2);

        //    private static double a = 0.1;
        //private static double b = 2.3;
        //private static double alpha = 0.2;
        //private static double beta = 0;


        Matrix m = new Matrix(3, 1);
            m[0, 0] = (Math.Pow(xn - a, 1 - alpha) - Math.Pow(x0 - a, 1 - alpha)) / (1 - alpha);
            m[1, 0] = (Math.Pow(xn - a, 2 - alpha) - Math.Pow(x0 - a, 2 - alpha)) / (2 - alpha) + a * m[0, 0];
            m[2, 0] = (Math.Pow(xn - a, 3 - alpha) - Math.Pow(x0 - a, 3 - alpha)) / (3 - alpha) + 2 * a * m[1, 0] - a * a * m[0, 0];
            Matrix t = matr.SolveWith(m);
            double sum = 0;
            for (int i = 0; i < 3; i++)
            {
                
                
                sum = sum + t[i, 0] * F(matr[1, i]);
            }
            //Console.WriteLine("m = [{0}, {1}, {2}]",m[0,0], m[1, 0], m[2, 0]);

           
            

          

            return sum;

        }

       

       

        

    }
}
