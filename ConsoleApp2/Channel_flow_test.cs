using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Channel_test
    {
        public Channel_test()
        {
            double Re = 1.0;
            int size = 25;
            double[] min_coord = new double[] { 0.0, 0.0 };
            double[] max_coord = new double[] { 3.0, 1.0 };
            int[] divs = new int[] { 3 * size, 2 * size + 1 };
            Grid grid = new Grid(min_coord, max_coord, divs);
            Dictionary<string, object> d = new Dictionary<string, object>();
            d.Add("nu", (double)1.0 / Re);
            NavierStokes2D ns = new NavierStokes2D(grid, d);
            Dictionary<string, object> bc1 = new Dictionary<string, object>();
            bc1.Add("u", 0.0);
            bc1.Add("v", 0.0);
            ns.Set_bc("north+south", bc1);
            Dictionary<string, object> bc2 = new Dictionary<string, object>();
            bc2.Add("u", 1.0);
            bc2.Add("v", 0.0);
            ns.Set_bc("west", bc2);
            Dictionary<string, object> bc3 = new Dictionary<string, object>();
            bc3.Add("dudn", 0.0);
            bc3.Add("dvdn", 0.0);
            ns.Set_bc("east", bc3);
            double t = 0;
            double dt;
            while(t < 0.25)
            {
                dt = ns.Find_suitable_dt();
                ns.Step(dt);
                t += dt;
                while (t < 0.05)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 0.1)
                {
                    dt = ns.Find_suitable_dt();
                    ns.Step(dt);
                    t += dt;
                }
                while (t < 0.15)
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
            }
        }
    }
}
