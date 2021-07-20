using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Cylinder_test
    {
        public Cylinder_test()
        {
            int size = 10;
            double Re = 100.0;
            double[] min = new double[] { -2.0, -2.0 };
            double[] max = new double[] { 20.0, 2.1 };
            int[] delta = new int[] { 22 * size, 41 * size / 10 };
            Grid grid = new Grid(min, max, delta);
            Dictionary<string, object> n = new Dictionary<string, object>();
            n.Add("mu", 1.0 / Re);
            NavierStokes2D ns = new NavierStokes2D(grid, n);
            Dictionary<string, object> bc1 = new Dictionary<string, object>();
            bc1.Add("u", 0.0);
            bc1.Add("v", 0.0);
            ns.Set_bc("south+north", bc1);
            double f(double x, double y, double t)
            {
                return 6.0 * (y + 2.0) * (2.1 - y) / (Math.Pow(4.1, 2));
            }
            Func<double,double, double, double> inl = f;
            Dictionary<string, object> bc2 = new Dictionary<string, object>();
            bc2.Add("u", inl);
            bc2.Add("v", 0.0);
            ns.Set_bc("west", bc2);
            Dictionary<string, object> bc3 = new Dictionary<string, object>();
            bc3.Add("dudn", 0.0);
            bc3.Add("dvdn", 0.0);
            ns.Set_bc("east", bc3);
            bool q(double x, double y)
            {
                return ((Math.Pow(x,2) + Math.Pow(y,2)) < 0.25);
            }
            Func<double, double, bool> w = q;
            ns.Set_obstacle(w);
            double[] e = new double[] { -3.0, -2.0 };
            double[] r = new double[] { -1.0, 0.0 };
            ns.Add_tracer(e, r);
            ns.Step(ns.Find_suitable_dt());
            ns.Step(ns.Find_suitable_dt());
            ns.Step(ns.Find_suitable_dt());
            ns.Step(ns.Find_suitable_dt());
            ns.Step(ns.Find_suitable_dt());
        }
    }
}
