using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Drop_test
    {
        public Drop_test()
        {
            int size = 100;
            double[] min = new double[] { 0.0, 0.0 };
            double[] max = new double[] { 6.0, 6.0 };
            int[] dim = new int[] { size, size };
            Grid grid = new Grid(min, max, dim);
            Dictionary<string, object> n = new Dictionary<string, object>();
            n.Add("mu", 0.1);
            n.Add("mu_gas", 0.0016);
            n.Add("rho", 1.0);
            n.Add("rho_gas", 0.0013);
            n.Add("surface_tension_coeff", 730.0);
            double[] g = new double[] { 0.0, -100.00 };
            n.Add("gravity", g);
            NavierStokes2D ns = new NavierStokes2D(grid, n);
            Dictionary<string, object> bc1 = new Dictionary<string, object>();
            bc1.Add("u", 0.0);
            bc1.Add("v", 0.0);
            ns.Set_bc("west+east", bc1);
            Dictionary<string, object> bc2 = new Dictionary<string, object>();
            bc2.Add("u", 0.0);
            bc2.Add("v", 0.0);
            ns.Set_bc("south+north", bc2);
            double[] g_m = new double[] { 0.0, 0.0 };
            double[] g_n = new double[] { 6.0, 6.0 };
            ns.Add_gas(g_m, g_n);
            double[] a_l = new double[] { 6.0, 2.0 };
            ns.Add_liquid(g_m, a_l);
            double[] a_d = new double[] { 3.0, 4.5 };
            ns.Add_drop(a_d, 0.3);
            ns.Step(ns.Find_suitable_dt());
        }
    }
}
