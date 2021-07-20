using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Cavity_test
    {
        public Cavity_test()
        {
            int power = 6;
            double Re = 1000.0;
            double[] min = new double[] { 0.0, 0.0 };
            double[] max = new double[] { 1.0, 1.0 };
            int[] delta = new int[] { (int)Math.Pow(2, power), (int)Math.Pow(2, power) };
            Grid grid = new Grid(min, max, delta);
            Dictionary<string, object> n = new Dictionary<string, object>();
            n.Add("nu", 1.0 / Re);
            NavierStokes2D ns = new NavierStokes2D(grid, n);
            Dictionary<string, object> bc1 = new Dictionary<string, object>();
            bc1.Add("u", 0.0);
            bc1.Add("v", 0.0);
            ns.Set_bc("west+east+south", bc1);
            Dictionary<string, object> bc2 = new Dictionary<string, object>();
            bc2.Add("u", -1.0);
            bc2.Add("v", 0.0);
            ns.Set_bc("north", bc2);
            double t = 0;
            double dt;
            while(t < 17)
            {
                dt = ns.Find_suitable_dt();
                ns.Step(dt);
                t += dt;
                while (t < 1.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 2.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 3.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 4.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 5.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 6.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 7.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 8.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 9.4)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 12)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 15)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
            }
        }
    }
}
