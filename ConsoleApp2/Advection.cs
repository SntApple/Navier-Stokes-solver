using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Advection
    {
        public Advection()
        {
            int size = 50;
            double[] min = new double[] { -1.0, -1.0 };
            double[] max = new double[] { 1.0, 1.0 };
            int[] delta = new int[] { size, size };
            Grid grid = new Grid(min, max, delta);
            Dictionary<string, object> n = new Dictionary<string, object>();
            n.Add("mu_liquid", 1.0);
            n.Add("mu_gas", 1.0);
            n.Add("rho_liquid", 1.0);   
            n.Add("rho_gas", 1.0);
            NavierStokes2D ns = new NavierStokes2D(grid, n);
        }
    }
}
