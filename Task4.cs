using System;
using System.Collections.Generic;
using System.Text;
using MathNet;
using MathNet.Numerics.LinearAlgebra.Complex;

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





    static class Task4
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





        private static double exact_value = 32.21951452884234295708696008290380201405;
        //private static double exact_value = 3.578861536040539915439859609644293194417;


        private static double[,] A;
        private static double[,] Nu;

        public static void Solution()
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
            // Mn = 2262.14 // ограничение третьей производной
            // 2262.14/6*abs((t-1.5)*(t-2.3)*(t-(1.5+2.3)/2)/(t-1.5)^0.2) // 6.17853
            Console.WriteLine("Составная квадратурная формула Ньютона-Котса с оптимальным шагом {0}" , skf_newton_kots(a, b, n,  Math.Ceiling((b - a) / answ[1]))[0]);
            


            Console.WriteLine("-------------------------------------------------------------------------------------------------");
            
            Console.WriteLine("Квадратурная формула Гаусса {0}", gaus(a, b, k));
            Console.WriteLine("Погрешность {0}", gaus(a, b, k)-exact_value);
            Console.WriteLine("Методическая погрешность {0}", 0);
            answ = skf_gaus(a, b, k);
            Console.WriteLine("Составная квадратурная формула Гаусса {0}", answ[0]);
            Console.WriteLine("Погрешность{0}", Math.Abs(answ[0]-exact_value));
            Console.WriteLine("Составная формула Гаусса с оптимальным шагом {0}", skf_gaus(a, b, n, Math.Ceiling((b - a) / answ[1]))[0]);





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
                Console.WriteLine("k {0} - s {1} - delta{2}", k / L, s[s.Count - 1], s[s.Count - 1]-exact_value);
                s.Add(0);
               
            }
            Console.WriteLine("kopt {0}", (b - a) / hopt);
            List<double> answer = new List<double>();
            answer.Add(s[s.Count - 2]);
            answer.Add(hopt);
            return answer;

        }

        private static List<double> skf_gaus(double a, double b, int n, double k = 1)
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

            int iterator = 0;
            while (Math.Abs(Rn) > eps)
            {
                iterator++;
                double h = (b - a) / k;
                for (int i = 0; i < k; i++)
                    s[it - 1] = s[it - 1] + gaus(a + i * h, a + (i + 1) * h, n);
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
               
                Console.WriteLine("k {0} - s {1} - delta{2}", k / L, s[s.Count - 1], s[s.Count - 1] - exact_value);
                if ((s.Count > 3))
                {
                    if ((Math.Abs(s[s.Count - 1] - exact_value) > Math.Abs(s[s.Count - 2] - exact_value)))
                    {
                        hopt = k / L;
                        //Console.WriteLine("Расхождение");
                        break;
                    }
                }

                s.Add(0);

               
            }
            Console.WriteLine("kopt {0}", (b - a) / hopt);
            List<double> answer = new List<double>();
            answer.Add(s[s.Count - 2]);
            answer.Add(Math.Abs(hopt));
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
        private static double gaus(double x0, double xn, int k)
        {
            //    private static double a = 0.1;
            //private static double b = 2.3;
            //private static double alpha = 0.2;
            //private static double beta = 0;
            double[] m = new double[2 * k];
            m[0] = (Math.Pow(xn - a, 1 - alpha) - Math.Pow(x0 - a, 1 - alpha)) / (1 - alpha);
            m[1] = (Math.Pow(xn - a, 2 - alpha) - Math.Pow(x0 - a, 2 - alpha)) / (2 - alpha) + a * m[0];
            m[2] = (Math.Pow(xn - a, 3 - alpha) - Math.Pow(x0 - a, 3 - alpha)) / (3 - alpha) + 2 * a * m[1] - a * a * m[0];
            //m[3] = (5 * Math.Pow(10 * xn - 1, 1 - alpha) * (672 * Math.Pow(xn, 3) + 72 * Math.Pow(xn, 2) + 8 * xn + 1)) / (12768 * Math.Pow(10, 1 - alpha)) - (5 * Math.Pow(10 * x0 - 1, 1 - alpha) * (672 * Math.Pow(x0, 3) + 72 * Math.Pow(x0, 2) + 8 * x0 + 1)) / (12768 * Math.Pow(10, 1 - alpha));
            //m[4] = (Math.Pow(10 * xn - 1, 1 - alpha) * (31920 * Math.Pow(xn, 4) + 3360 * Math.Pow(xn, 3) + 360 * Math.Pow(xn, 2) + 40 * xn + 5)) / (153216 * Math.Pow(10, 1 - alpha)) - (Math.Pow(10 * x0 - 1, 1 - alpha) * (31920 * Math.Pow(x0, 4) + 3360 * Math.Pow(x0, 3) + 360 * Math.Pow(x0, 2) + 40 * x0 + 5)) / (153216 * Math.Pow(10, 1 - alpha));
            //m[5] = (5 * Math.Pow(10 * xn - 1, 1 - alpha) * (306432 * Math.Pow(xn, 5) + 31920 * Math.Pow(xn, 4) + 3360 * Math.Pow(xn, 3) + 360 * Math.Pow(xn, 2) + 40 * xn + 5)) / (8886528 * Math.Pow(10, 1 - alpha)) - (5 * Math.Pow(10 * x0 - 1, 1 - alpha) * (306432 * Math.Pow(x0, 5) + 31920 * Math.Pow(x0, 4) + 3360 * Math.Pow(x0, 3) + 360 * Math.Pow(x0, 2) + 40 * x0 + 5)) / (8886528 * Math.Pow(10, 1 - alpha));

            //Console.WriteLine("m[3] - {0}", m[3]);
            //double test = 0;

            //double test1 = (Math.Pow(xn - a, 0.8) * (62.5 * Math.Pow(a, 3) + 50 * Math.Pow(a, 2) * xn + 45 * a * Math.Pow(xn, 2) + 42 * Math.Pow(xn, 3)) - Math.Pow(x0 - a, 0.8) * (62.5 * Math.Pow(a, 3) + 50 * Math.Pow(a, 2) * x0 + 45 * a * Math.Pow(x0, 2) + 42 * Math.Pow(x0, 3))) * 5 / 798;
            //test1 = (Math.Pow(xn - a, 0.8) * (62.5 * Math.Pow(a, 3) + 50 * Math.Pow(a, 2) * xn + 45 * a * Math.Pow(xn, 2) + 42 * Math.Pow(xn, 3)) - Math.Pow(x0 - a, 0.8) * (62.5 * Math.Pow(a, 3) + 50 * Math.Pow(a, 2) * x0 + 45 * a * Math.Pow(x0, 2) + 42 * Math.Pow(x0, 3))) * 5 / 798;

            //m[4] = (Math.Pow(xn - a, 0.8) * (625  * Math.Pow(a, 4) + 500 * Math.Pow(a, 3) * xn + 450 * Math.Pow(a, 2) * Math.Pow(xn, 2) + 420 * a * Math.Pow(xn, 3) + 399 * Math.Pow(xn, 4))- Math.Pow(x0 - a, 0.8) * (625 * Math.Pow(a, 4) + 500 * Math.Pow(a, 3) * x0 + 450 * Math.Pow(a, 2) * Math.Pow(x0, 2) + 420* a * Math.Pow(x0, 3) + 399 * Math.Pow(x0, 4)))*5/9576;

            m[3] = 5 * (Math.Pow(xn - 1.5, 0.8) * (224 * Math.Pow(xn, 3) + 360 * Math.Pow(xn, 2) + 600 * xn + 1125) - Math.Pow(x0 - 1.5, 0.8) * (224 * Math.Pow(x0, 3) + 360 * Math.Pow(x0, 2) + 600 * x0 + 1125)) / 4256;
            m[4] = 5*(Math.Pow(xn - 1.5, 0.8) * (2128 * Math.Pow(xn, 4) + 3360 * Math.Pow(xn, 3) + 5400 * Math.Pow(xn, 2) + 9000 * xn + 16875) - Math.Pow(x0 - 1.5, 0.8) * (2128 * Math.Pow(x0, 4) + 3360 * Math.Pow(x0, 3) + 5400 * Math.Pow(x0, 2) + 9000 * x0 + 16875))/ 51072;
            m[5] = 5*(Math.Pow(xn - 1.5, 0.8) * (34048 * Math.Pow(xn, 5) + 53200 * Math.Pow(xn, 4) + 84000 * Math.Pow(xn, 3) + 135000 * Math.Pow(xn, 2) + 225000 * xn + 421875) - Math.Pow(x0 - 1.5, 0.8) * (34048 * Math.Pow(x0, 5) + 53200 * Math.Pow(x0, 4) + 84000 * Math.Pow(x0, 3) + 135000 * Math.Pow(x0, 2) + 225000 * x0 + 421875)) / 987392;
            //Console.WriteLine("test1 - {0} - {1} // exact {2}", test1, test1 - m[3], m[3]);
            //Console.WriteLine("test2 - {0} - {1} // exact {2}", test2, test2- m[4], m[4]);
            //Console.WriteLine("test3 - {0} - {1} // exact {2}", test3, test3 - m[5], m[5]);
            Matrix moments = new Matrix(k, k);
            moments[0, 0] = m[0];
            moments[0, 1] = m[1];
            moments[1, 0] = m[1];
            moments[0, 2] = m[2];
            moments[1, 1] = m[2];
            moments[2, 0] = m[2];
            moments[2, 1] = m[3];
            moments[1, 2] = m[3];
            moments[2, 2] = m[4];
            Matrix b = new Matrix(k, 1);
            b[0, 0] = -m[3];
            b[1, 0] = -m[4];
            b[2, 0] = -m[5];

            Matrix aa = moments.SolveWith(b);
            //Matrix aaa = new Matrix ()
          
            double[] aaa_ = new double[4] { aa[0, 0], aa[1, 0], aa[2, 0] , 1 };
            System.Numerics.Complex[] x_sol_ = new System.Numerics.Complex[4];
            x_sol_ = MathNet.Numerics.FindRoots.Polynomial(aaa_);
            double[] x_sol = new double[x_sol_.Length];
            for (int i = 0; i < x_sol_.Length; i++)
            {
                if (x_sol_[i].Imaginary == 0)
                {
                    x_sol[i] = x_sol_[i].Real;
                }

            }
            Matrix matr = new Matrix(3, 3);
            Matrix mat_b = new Matrix(3, 1);
            Matrix s = new Matrix(3, 1);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    matr[i, j] = (Math.Pow(x_sol[j], i));
                    //matr[i, j] = MathNet.Numeric
                }
            }
            for (int i = 0; i < 3; i++)
            {
                mat_b[i, 0] = m[i];
            }


            //MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex> matr = DenseMatrix.OfArray(new System.Numerics.Complex[,] { { 1, 1, 1 }, { 0, 0, x_sol_[2] }, { 0, 0, x_sol_[2] * x_sol_[2].Real } });
            //MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex> mat_b = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>.Build.Dense(new System.Numerics.Complex[] { m[0], m[1], m[2] });

            //var coef = matr.Solve(mat_b);









            Matrix coef = matr.SolveWith(mat_b);
            double sum = 0;
            //System.Numerics.Complex sum = 0;
            for (int i =0; i<3; i++)
            {
                sum = sum + F(x_sol[i]) * coef[i,0];
            }

            //return 0;
            return sum;
            
        }









    }
}
