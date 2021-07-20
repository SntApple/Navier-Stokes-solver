using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Bubble_test
    {
        public Bubble_test()
        {
            int size = 25;
            double[] min_coord = new double[] { 0, 0 };
            double[] max_corrd = new double[] { 2, 4 };
            int[] delta = new int[] { size, 2 * size };
            Grid grid = new Grid(min_coord, max_corrd, delta);
            Dictionary<string, object> q = new Dictionary<string, object>();
            q.Add("mu", 0.002);
            q.Add("mu_gas", 3.2e-5);
            q.Add("rho", 1.0);
            q.Add("rho_gas", 0.0013);
            q.Add("surface_tension_coeff", 1.46);
            double[] gravity = new double[] { 0, -5 };
            q.Add("gravity", gravity);
            NavierStokes2D ns = new NavierStokes2D(grid, q);
            Dictionary<string, object> bc1 = new Dictionary<string, object>();
            bc1.Add("u", 0.0);
            bc1.Add("v", 0.0);
            ns.Set_bc("west+east", bc1);
            Dictionary<string, object> bc2 = new Dictionary<string, object>();
            bc2.Add("u", 0.0);
            bc2.Add("v", 0.0);
            ns.Set_bc("south+north", bc2);
            double[] centre = new double[] { 1, 1 };
            ns.Add_bubble(centre, 0.5);
            double t = 0.0;
            double dt;
            while(t < 0.5)
            {
                dt = ns.Find_suitable_dt();
                ns.Step(dt);
                t += dt;
                while (t < 0.1)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 0.2)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 0.3)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 0.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
            }
        }
    }
}
