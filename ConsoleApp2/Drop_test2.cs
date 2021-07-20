using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Drop_test2
    {
        public Drop_test2()
        {
            int size = 96;
            double surf_coeff = 0.357;
            double radius = 0.25;
            double[] min = new double[] { 0.0, 0.0 };
            double[] max = new double[] { 1.0, 1.0 };
            int[] delta = new int[] { size, size };
            Grid grid = new Grid(min, max, delta);
            Dictionary<string, object> n = new Dictionary<string, object>();
            n.Add("mu", 1.0);
            n.Add("mu_gas", 1.0);
            n.Add("rho", 4.0);
            n.Add("rho_gas", 4.0);
            n.Add("surface_tension_coeff", surf_coeff);
            double[] g = new double[] { 0.0, 0.0 };
            n.Add("gravity", g);
            n.Add("curv_method", "csf");
            NavierStokes2D ns = new NavierStokes2D(grid, n);
            Dictionary<string, object> bc1 = new Dictionary<string, object>();
            bc1.Add("u", 0.0);
            bc1.Add("v", 0.0);
            ns.Set_bc("west+east", bc1);
            Dictionary<string, object> bc2 = new Dictionary<string, object>();
            bc2.Add("u", 0.0);
            bc2.Add("v", 0.0);
            ns.Set_bc("south+north", bc2);
            double[] q = new double[] { 0.0, 0.0 };
            double[] qq = new double[] { 1.0, 1.0 };
            ns.Add_gas(q, qq);
            double[] w = new double[] { 0.5, 0.5 };
            ns.Add_drop(w, radius);
            double dt = 1e-5;
            while (true)
            {
                ns.Step(dt);
                char a = Console.ReadKey().KeyChar;
                if (a == 'q')
                    break;
                else if (a == 'p')
                    continue;
            }
        }
    }
}
