using System;
using System.Collections.Generic;
using System.Linq;
namespace ConsoleApp2
{
    class NavierStokes2D
    {
        readonly string EAST = "east";
        readonly string WEST = "west";
        readonly string NORTH = "north";
        readonly string SOUTH = "south";

        readonly string MASS_SCALE = "scale";
        readonly string MASS_ADD = "add";
        readonly string MASS_IGNORE = "ignore";

        readonly string CURV_CSF = "csf";
        readonly string CURV_DAC = "dac";
        readonly string CURV_MDAC = "mdac";

        Grid grid;
        Utils utils = new Utils();
        Volume vof = new Volume();
        double t = 0.0;
        double[,] p;
        double[,] u;
        double[,] v;
        long[,] mask;
        double[,] c;
        double[,] dye;
        Dictionary<string, object> bc = new Dictionary<string, object>();
        bool bc_finalized = false;
        List<double[]> particles = new List<double[]>();
        int it = 0;
        List<double> gas_cell_fraction = new List<double>();
        List<int> gas_fraction = new List<int>();
        List<double[]> gas_centre = new List<double[]>();
        List<double[]> liquid_centre = new List<double[]>();
        List<double> record_time = new List<double>();
        Dictionary<string, object> default_params = new Dictionary<string, object>(11);
        
        public NavierStokes2D(Grid grid, Dictionary<string, object> kwargs)
        {
            this.grid = grid;
            int[] shape = grid.Get_shape();
            this.p = new double[shape[0] + 1, shape[1] + 1];
            this.u = new double[p.GetLength(0) - 1, p.GetLength(1)];
            this.v = new double[p.GetLength(0), p.GetLength(1) - 1];
            this.mask = new long[p.GetLength(0), p.GetLength(1)];
            for (int i = 0; i < mask.GetLength(0); ++i)
            {
                for (int j = 0; j < mask.GetLength(1); ++j)
                    mask[i, j] = 1;
            }
            this.c = new double[p.GetLength(0), p.GetLength(1)];
            this.dye = new double[p.GetLength(0), p.GetLength(1)];
            Dictionary<string, object> ea = new Dictionary<string, object>();
            Dictionary<string, object> we = new Dictionary<string, object>();
            Dictionary<string, object> no = new Dictionary<string, object>();
            Dictionary<string, object> so = new Dictionary<string, object>();
            this.bc.Add(this.EAST, ea);
            this.bc.Add(this.WEST, we);
            this.bc.Add(this.NORTH, no);
            this.bc.Add(this.SOUTH, so);

            this.default_params.Add("rho_liquid", 1.0);
            this.default_params.Add("rho_gas", 0.0013);
            this.default_params.Add("mu_liquid", 1.0);
            this.default_params.Add("mu_gas", 0.016);
            this.default_params.Add("surface_tension_coeff", 0.0);
            double[] ar = new double[2] { 0.0, 0.0 };
            this.default_params.Add("gravity", ar);
            this.default_params.Add("two_phase", false);
            this.default_params.Add("use_dye", false);
            this.default_params.Add("curv_method", this.CURV_MDAC);
            this.default_params.Add("mass_conservation", this.MASS_ADD);
            this.default_params.Add("property_smoothing", false);
            for (int index = 0; index < kwargs.Count; ++index)
            {
                var item = kwargs.ElementAt(index);
                var itemKey = item.Key;
                if (itemKey == "mu")
                    kwargs["mu_liquid"] = kwargs["mu"];
                if (itemKey == "nu")
                    kwargs["nu_fluid"] = kwargs["nu"];
                if (itemKey == "rho")
                    kwargs["rho_liquid"] = kwargs["rho"];
            }

            foreach (KeyValuePair<string, object> en in kwargs)
                if (en.Key == "rho_liquid")
                    default_params["rho_liquid"] = (double)kwargs["mu"];
            foreach (KeyValuePair<string, object> en in kwargs)
                if (en.Key == "rho_gas")
                    default_params["rho_gas"] = (double)kwargs["rho_gas"];

            for (int index = 0; index < kwargs.Count; ++index)
            {
                var item = kwargs.ElementAt(index);
                var itemKey = item.Key;
                if (itemKey == "nu_fluid")
                    kwargs["mu_liquid"] = (double)kwargs["nu_fluid"] * (double)default_params["rho_liquid"];
                else if (itemKey == "nu_gas")
                    kwargs["mu_gas"] = (double)kwargs["nu_gas"] * (double)default_params["rho_gas"];
            }
            for (int index = 0; index < this.default_params.Count; ++index)
            {
                var item = this.default_params.ElementAt(index);
                var itemKey = item.Key;
                for (int i = 0; i < kwargs.Count; ++i)
                {
                    var it = kwargs.ElementAt(i);
                    var itKey = it.Key;
                    if (itemKey == itKey)
                        this.default_params[itemKey] = kwargs[itKey];
                }
            }

        }


        public void Set_ic(Dictionary<string, object> kwargs)
        {
            foreach (KeyValuePair<string, object> ic in kwargs)
            {
                if (ic.Key != "u" && ic.Key != "v" && ic.Key != "p" && ic.Key != "c")
                    continue;
                object f = ic.Value;

                double[,] a = (double[,])this.default_params[ic.Key];
                if (f.GetType() == typeof(Func<double, double, object>))
                {
                    double xOffset;
                    if (ic.Key == "p" || ic.Key == "v" || ic.Key == "c")
                        xOffset = 0.5;
                    else
                        xOffset = 0;
                    double yOffset;
                    if (ic.Key == "p" || ic.Key == "u" || ic.Key == "c")
                        yOffset = 0.5;
                    else
                        yOffset = 0;
                    for (int i = 0; i < a.GetLength(0); ++i)
                    {
                        for (int j = 0; j < a.GetLength(1); ++j)
                        {
                            double[] node = this.grid.Getitem(i - xOffset, j - yOffset);
                            a[i, j] = (double)((Func<double, double, object>)f)(node[0], node[1]);
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < a.GetLength(0); ++i)
                    {
                        for (int j = 0; j < a.GetLength(1); ++j)
                            a[i, j] = (double)f;
                    }
                }
            }
        }

        public void Set_bc(string boundaries, Dictionary<string, object> kwargs)
        {
            string[] b = boundaries.Split('+');
            for (int index = 0; index < this.bc.Count; ++index)
            {
                var item = this.bc.ElementAt(index);
                var itemKey = item.Key;
                foreach (string bs in b)
                {
                    if (bs == itemKey)
                    {
                        Dictionary<string, object> itemVal = (Dictionary<string, object>)item.Value;
                        foreach (KeyValuePair<string, object> kw in kwargs)
                        {
                            itemVal[kw.Key] = kw.Value;
                        }
                    }
                    foreach (KeyValuePair<string, object> c in kwargs)
                    {
                        if (c.Key == "p")
                        {
                            if (bs == this.NORTH)
                                for (int i = 0; i < mask.GetLength(0); ++i)
                                {
                                    this.mask[i, mask.GetLength(1) - 1] = 3;
                                }
                            else if (bs == this.SOUTH)
                                for (int i = 0; i < mask.GetLength(0); ++i)
                                {
                                    mask[i, 0] = 3;
                                }
                            else if (bs == this.WEST)
                                for (int j = 0; j < mask.GetLength(1); ++j)
                                {
                                    mask[0, j] = 3;
                                }
                            else if (bs == this.EAST)
                                for (int j = 0; j < mask.GetLength(1); ++j)
                                {
                                    mask[mask.GetLength(0) - 1, j] = 3;
                                }
                        }
                        if (c.Key == "v")
                        {
                            if (bs == this.NORTH || bs == this.SOUTH)
                            {
                                Dictionary<string, object> o = (Dictionary<string, object>)this.bc[bs];
                                o["dpdn"] = (double)0.0;
                            }
                        }
                        if (c.Key == "u")
                        {
                            if (bs == this.EAST || bs == this.WEST)
                            {
                                Dictionary<string, object> o = (Dictionary<string, object>)this.bc[bs];
                                o["dpdn"] = (double)0.0;
                            }
                        }
                    }
                }
            }
        }



        void Reset_time()
        {
            this.t = 0.0;
        }

        public void Add_obstacle(double[] min_coord, double[] max_coord)
        {
            double[] scale = new double[2];
            for (int i = 0; i < 2; ++i)
            {
                scale[i] = this.grid.divs[i] / (this.grid.max_coord[i] - this.grid.min_coord[i]);

            }
            double[] low = new double[2];
            double[] high = new double[2];
            for (int i = 0; i < 2; ++i)
            {
                low[i] = scale[i] * (min_coord[i] - this.grid.min_coord[i]) + 1;
                high[i] = scale[i] * (max_coord[i] - this.grid.min_coord[i]) + 1;
            }
            for (int i = (int)Math.Round(low[0]); i < (int)Math.Round(high[0]); ++i)
            {
                for (int j = (int)Math.Round(low[1]); j < (int)Math.Round(high[1]); ++j)
                    this.mask[i, j] = 0;
            }
        }

        public void Set_obstacle(Func<double, double, bool> f)
        {
            long[,] mask = this.mask;
            for (int i = 1; i < mask.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < mask.GetLength(1); ++j)
                {
                    double[] node = grid.Getitem(i - 0.5, j - 0.5);
                    if (((Func<double, double, bool>)f)(node[0], node[1]))
                    {
                        mask[i, j] = 0;
                    }
                    else
                    {
                        mask[i, j] = 1;
                    }
                }
            }
        }

        public void Add_particle(double[] pos)
        {
            this.particles.Add(pos);
        }

        public void Add_tracer(double[] min_coord, double[] max_coord)
        {
            Volume volf = new Volume();
            this.default_params["use_dye"] = true;
            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Get_delta();

            double[] x = utils.Linspace(low[0] - delta[0], high[0] + delta[0], this.grid.divs[0] + 3);
            double[] y = utils.Linspace(low[1] - delta[1], high[1] + delta[1], this.grid.divs[1] + 3);
            this.dye = utils.Maximum(this.dye, volf.Rectangle_fraction(min_coord[0], min_coord[1], max_coord[0], max_coord[1], x, y));
        }

        public void Add_gas(double[] min_coord, double[] max_coord)
        {
            Volume volf = new Volume();
            this.default_params["two_phase"] = true;
            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Delta();

            double[] x = utils.Linspace(low[0] - delta[0], high[0] + delta[0], this.grid.divs[0] + 3);
            double[] y = utils.Linspace(low[1] - delta[1], high[1] + delta[1], this.grid.divs[1] + 3);
            this.c = utils.Maximum(this.c, volf.Rectangle_fraction(min_coord[0], min_coord[1], max_coord[0], max_coord[1], x, y));
        }

        public void Add_bubble(double[] centre, double radius, double scale_x = 1.0, double scale_y = 1.0)
        {
            Volume volf = new Volume();
            this.default_params["two_phase"] = true;
            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Delta();

            double[] x = utils.Linspace(low[0] - delta[0], high[0] + delta[0], this.grid.divs[0] + 3);
            double[] y = utils.Linspace(low[1] - delta[1], high[1] + delta[1], this.grid.divs[1] + 3);
            for (int i = 0; i < x.Length; ++i)
                x[i] /= scale_x;
            for (int i = 0; i < y.Length; ++i)
                y[i] /= scale_y;
            double[,] d = volf.Circle_fraction(centre[0] / scale_x, centre[1] / scale_y, radius, x, y);
            this.c = utils.Maximum(this.c, volf.Circle_fraction(centre[0] / scale_x, centre[1] / scale_y, radius, x, y));
        }

        public void Add_liquid(double[] min_coord, double[] max_coord)
        {
            this.default_params["two_phase"] = true;
            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Delta();
            double[] x = utils.Linspace(low[0] - delta[0], high[0] + delta[0], this.grid.divs[0] + 3);
            double[] y = utils.Linspace(low[1] - delta[1], high[1] + delta[1], this.grid.divs[1] + 3);
            double[,] d = vof.Rectangle_fraction(min_coord[0], min_coord[1], max_coord[0], max_coord[1], x, y);
            for(int i = 0; i < d.GetLength(0); ++i)
            {
                for(int j = 0; j < d.GetLength(1); ++j)
                {
                    d[i, j] = 1.0 - d[i, j];
                }
            }
            this.c = utils.Minimum(this.c, d);
        }

        public void Add_drop(double[] centre, double radius, double scale_x = 1.0, double scale_y = 1.0)
        {
            this.default_params["two_phase"] = true;
            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Delta();
            double[] x = utils.Linspace(low[0] - delta[0], high[0] + delta[0], this.grid.divs[0] + 3);
            double[] y = utils.Linspace(low[1] - delta[1], high[1] + delta[1], this.grid.divs[1] + 3);
            for (int i = 0; i < x.Length; ++i)
                x[i] /= scale_x;
            for (int i = 0; i < y.Length; ++i)
                y[i] /= scale_y;
            double[,] d = vof.Circle_fraction(centre[0] / scale_x, centre[1] / scale_y, radius, x, y);
            for(int i = 0; i < d.GetLength(0); ++i)
            {
                for(int j = 0; j < d.GetLength(1); ++j)
                {
                    d[i, j] = 1.0 - d[i, j];
                }
            }
            this.c = utils.Minimum(this.c, d);
        }

        public void Set_interface(int orientation, Func<double, double> F)
        {
            this.default_params["two_phase"] = true;
            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Delta();
            double[,] c = this.c;

            double[] x = utils.Linspace(low[0] - delta[0], high[0] + delta[0], this.grid.divs[0] + 3);
            double[] y = utils.Linspace(low[1] - delta[1], high[1] + delta[1], this.grid.divs[1] + 3);
            int j = 0;
            if (orientation == 0 || orientation == 2)
            {
                for (int i = 0; i < x.Length - 1; ++i)
                {
                    double h = (F(x[i + 1]) - F(x[i])) / (delta[0] * delta[1]) + 1.0;
                    int n = (int)Math.Floor(h);
                    double r = h - n;
                    if (orientation == 0)
                    {
                        for (j = 0; j < n; ++j)
                        {
                            c[i, j] = 0.0;
                        }
                        for (j = n; j < c.GetLength(1); ++j)
                        {
                            c[i, j] = 1.0;
                        }
                        c[i, n] -= r;
                    }
                    else
                    {
                        for (j = c.GetLength(1) - n; j < c.GetLength(1); ++j)
                        {
                            c[i, j] = 0.0;
                        }
                        for (j = 0; j < c.GetLength(1) - n; ++j)
                        {
                            c[i, j] = 1.0;
                        }
                        c[i, c.GetLength(1) - n - 1] -= r;
                    }
                }
            }
            else
            {
                for (j = 0; j < y.Length - 1; ++j)
                {
                    double h = (F(y[j + 1]) - F(y[j])) / (delta[0] * delta[1]) + 1.0;
                    int n = (int)Math.Floor(h);
                    double r = h - Math.Floor(h);

                    if (orientation == 3)
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            c[i, j] = 0.0;
                        }
                        for (int i = n; i < c.GetLength(0); ++i)
                        {
                            c[i, j] = 1.0;
                        }
                        c[n, j] -= r;
                    }
                    else
                    {
                        for (int i = c.GetLength(0) - n; i < c.GetLength(0); ++i)
                        {
                            c[i, j] = 0.0;
                        }
                        for (int i = 0; i < c.GetLength(0) - n; ++i)
                        {
                            c[i, j] = 1.0;
                        }
                        c[c.GetLength(0) - n - 1, j] -= r;
                    }
                }
            }
        }

        public void Step(double dt, double beta = 1.0, double gamma = 1.0)
        {
            double[] delta = this.grid.Delta();
            long[,] mask = this.mask;
            double[,] c = this.c;
            double[,] u = this.u;
            double[,] v = this.v;
            double[,] p = this.p;
            Csf csf = new Csf();
            if (!this.bc_finalized)
            {
                utils.Remove_singularity(this.mask);
                this.bc_finalized = true;
            }
            Update_velocity_bc();
            Update_pressure_bc();
            Update_level_bc();
            for (int i = 0; i < this.particles.Count; ++i)
            {
                double[] interpolate_particles = Interpolate_velocity(this.particles[i]);
                this.particles[i][0] += dt * interpolate_particles[0];
                this.particles[i][1] += dt * interpolate_particles[1];
            }

            if ((bool)this.default_params["use_dye"])
            {
                vof.Advect(this.dye, u, v, mask, delta, dt, this.it);
            }
            if ((bool)this.default_params["two_phase"])
            {
                vof.Advect(c, u, v, mask, delta, dt, this.it);
                utils.Neumann_bc(c, mask);
                if ((bool)this.default_params["property_smoothing"])
                {
                    c = csf.Smooth(c, mask, delta);
                    utils.Neumann_bc(c, mask);
                }
            }
            double[,] rho = new double[c.GetLength(0), c.GetLength(1)];
            double[,] mu = new double[c.GetLength(0), c.GetLength(1)];
            for (int i = 0; i < rho.GetLength(0); ++i)
            {
                for (int j = 0; j < rho.GetLength(1); ++j)
                {
                    rho[i, j] = (double)this.default_params["rho_liquid"] * (1.0 - c[i, j]) + (double)this.default_params["rho_gas"] * c[i, j];
                }
            }
            for (int i = 0; i < mu.GetLength(0); ++i)
            {
                for (int j = 0; j < mu.GetLength(1); ++j)
                {
                    mu[i, j] = (double)this.default_params["mu_liquid"] * (1.0 - c[i, j]) + (double)this.default_params["mu_gas"] * c[i, j];
                }
            }
            double[,] rho_u = new double[rho.GetLength(0) - 1, rho.GetLength(1)];
            for (int i = 0; i < rho_u.GetLength(0); ++i)
            {
                for (int j = 0; j < rho_u.GetLength(1); ++j)
                {
                    rho_u[i, j] = 0.5 * (rho[i + 1, j] + rho[i, j]);
                }
            }
            double[,] rho_v = new double[rho.GetLength(0), rho.GetLength(1) - 1];
            for (int i = 0; i < rho_v.GetLength(0); ++i)
            {
                for (int j = 0; j < rho_v.GetLength(1); ++j)
                    rho_v[i, j] = 0.5 * (rho[i, j + 1] + rho[i, j]);
            }

            double[,] du;
            double[,] dv;
            Calc_velocity_change(dt, rho, rho_u, rho_v, mu, beta, gamma, out du, out dv);
            for (int i = 0; i < u.GetLength(0); ++i)
            {
                for (int j = 0; j < u.GetLength(1); ++j)
                {
                    u[i, j] += du[i, j];
                }
            }
            for (int i = 0; i < v.GetLength(0); ++i)
            {
                for (int j = 0; j < v.GetLength(1); ++j)
                {
                    v[i, j] += dv[i, j];
                }
            }
            this.t += dt;
            this.it += 1;
            Update_tentative_velocity_bc();
            double[,] f = Poisson_rhs(dt);
            double[,] phi = new double[p.GetLength(0), p.GetLength(1)];
            Update_phi_bc(phi, beta);
            for (int i = 0; i < p.GetLength(0); ++i)

            {
                for (int j = 0; j < p.GetLength(1); ++j)
                {
                    p[i, j] *= beta;
                }
            }
            PPE pois = new PPE();
            double [,]poi = pois.Poisson(phi, f, delta, mask, rho, this.bc, this.grid, this.t);
            for (int i = 0; i < p.GetLength(0); ++i)
            {
                for (int j = 0; j < p.GetLength(1); ++j)
                {
                    p[i, j] += poi[i, j];
                }
            }
            for (int i = 1; i < p.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < p.GetLength(1) - 1; ++j)
                {
                    p[i, j] *= (mask[i, j] & 1);
                }
            }
            for (int i = 1; i < u.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < u.GetLength(1) - 1; ++j)
                {
                    u[i, j] -= dt * (phi[i + 1, j] - phi[i, j]) /
                        (delta[0] * rho_u[i, j]) * (mask[i + 1, j] & mask[i, j] & 1);
                }
            }
            for(int i = 1; i < v.GetLength(0) - 1; ++i)
            {
                for(int j = 1; j < v.GetLength(1) - 1; ++j)
                {
                    v[i, j] -= dt * (phi[i, j + 1] - phi[i, j]) /
                        (delta[1] * rho_v[i, j]) * (mask[i, j + 1] & mask[i, j] & 1);
                }
            }
            for (int j = 1; j < u.GetLength(1) - 1; ++j)
            {
                u[0, j] -= dt * (phi[1, j] - phi[0, j]) /
                    (delta[0] * rho_u[0, j]) * (mask[1, j] & 1);
                u[u.GetLength(0) - 1, j] -= dt * (phi[phi.GetLength(0) - 1, j] - phi[phi.GetLength(0) - 2, j]) /
                   (delta[0] * rho_u[rho_u.GetLength(0) - 1, j]) * (mask[mask.GetLength(0) - 2, j] & 1);
            }

            for (int i = 1; i < v.GetLength(0) - 1; ++i)
            {
                v[i, 0] -= dt * (phi[i, 1] - phi[i, 0]) /
                    (delta[1] * rho_v[i, 0]) * (mask[i, 1] & 1);
                v[i, v.GetLength(1) - 1] -= dt * (phi[i, phi.GetLength(1) - 1] - phi[i, phi.GetLength(1) - 2]) /
                   (delta[1] * rho_v[i, rho_v.GetLength(1) - 1]) * (mask[i, mask.GetLength(1) - 2] & 1);
            }
            utils.Write_f(p, @"C:\Users\Grisha\Desktop\cs.csv");
            utils.Write_f(u, @"C:\Users\Grisha\Desktop\cs_u.csv");
            utils.Write_f(v, @"C:\Users\Grisha\Desktop\cs_v.csv");
            utils.Write_f(rho, @"C:\Users\Grisha\Desktop\cs_rho.csv");
            if ((bool)this.default_params["two_phase"])
                Update_records(rho);
        }

        double[,] Poisson_rhs(double dt)
        {
            double[,] u = this.u;
            double[,] v = this.v;
            double[,] p = this.p;
            double[] delta = this.grid.Delta();
            long[,] mask = this.mask;

            double net_outflow;
            double sum_u = 0.0;
            double sum_v = 0.0;
            for (int j = 1; j < u.GetLength(1) - 1; ++j)
            {
                sum_u += u[u.GetLength(0) - 1, j] - u[0, j];
            }
            for(int i = 1; i < v.GetLength(0) - 1; ++i)
            {
                sum_v += v[i, v.GetLength(1) - 1] - v[i, 0];
            }
            net_outflow = sum_u * delta[1] + sum_v * delta[0];
            double outflow_length = 0.0;
            double outflow = 0.0;
            bool dirichlet_used = false;
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[NORTH])
            {
                if (q.Key == "p")
                {
                    dirichlet_used = true;
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[SOUTH])
            {
                if (q.Key == "p")
                {
                    dirichlet_used = true;
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[EAST])
            {
                if (q.Key == "p")
                {
                    dirichlet_used = true;
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[WEST])
            {
                if (q.Key == "p")
                {
                    dirichlet_used = true;
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[NORTH])
            {
                if (q.Key == "dvdn")
                {
                    double mask_sum = 0;
                    double v_sum = 0;
                    for (int i = 1; i < mask.GetLength(0) - 1; ++i)
                    {
                        mask_sum += mask[i, mask.GetLength(1) - 2];
                        v_sum += v[i, v.GetLength(1) - 1];
                    }

                    outflow_length += mask_sum * delta[0];
                    outflow += v_sum * delta[0];
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[SOUTH])
            {
                if (q.Key == "dvdn")
                {
                    double mask_sum = 0;
                    double v_sum = 0;
                    for (int i = 1; i < mask.GetLength(0) - 1; ++i)
                    {
                        mask_sum += mask[i, 1];
                        v_sum += v[i, 0];
                    }
                    outflow_length += mask_sum * delta[0];
                    outflow -= v_sum * delta[0];
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[WEST])
            {
                if (q.Key == "dudn")
                {
                    double mask_sum = 0;
                    double v_sum = 0;
                    for (int j = 1; j < mask.GetLength(1) - 1; ++j)
                    {
                        mask_sum += mask[1, j];
                        v_sum += v[0, j];
                    }
                    outflow_length += mask_sum * delta[1];
                    outflow -= v_sum * delta[1];
                }
            }
            foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[EAST])
            {
                if (q.Key == "dudn")
                {
                    double mask_sum = 0;
                    double v_sum = 0;
                    for (int j = 1; j < mask.GetLength(1) - 1; ++j)
                    {
                        mask_sum += mask[mask.GetLength(1) - 2, j];
                        v_sum += v[v.GetLength(1) - 1, j];
                    }
                    outflow_length += mask_sum * delta[1];
                    outflow += v_sum * delta[1];
                }
            }
            if ((!dirichlet_used) && outflow_length > 0.0 && ((string)this.default_params["mass_conservation"] == MASS_ADD || (string)this.default_params["mass_conservation"] == MASS_SCALE))
            {
                if (outflow == 0.0 || (string)this.default_params["mass_conservation"] == this.MASS_ADD)
                {
                    double flow_corr = net_outflow / outflow_length;
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[NORTH])
                    {
                        if (q.Key == "dvdn")
                        {
                            for (int i = 1; i < v.GetLength(0) - 1; ++i)
                            {
                                v[i, v.GetLength(1) - 1] -= mask[i, mask.GetLength(1) - 2] * flow_corr;
                            }
                        }
                    }
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[SOUTH])
                    {
                        if (q.Key == "dvdn")
                        {
                            for (int i = 1; i < v.GetLength(0) - 1; ++i)
                            {
                                v[i, 0] += mask[i, 1] * flow_corr;
                            }
                        }
                    }
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[WEST])
                    {
                        if (q.Key == "dudn")
                        {
                            for (int j = 1; j < u.GetLength(1) - 1; ++j)
                            {
                                u[0, j] += mask[1, j] * flow_corr;
                            }
                        }
                    }
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[EAST])
                    {
                        if (q.Key == "dudn")
                        {
                            for (int j = 1; j < u.GetLength(1) - 1; ++j)
                            {
                                u[u.GetLength(0) - 1, j] -= mask[mask.GetLength(0) - 2, j] * flow_corr;
                            }
                        }
                    }
                }
                else
                {
                    double flow_corr = 1.0 - net_outflow / outflow;
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[NORTH])
                    {
                        if (q.Key == "dvdn")
                        {
                            for (int i = 1; i < v.GetLength(0) - 1; ++i)
                            {
                                v[i, v.GetLength(1) - 1] *= flow_corr;
                            }
                        }
                    }
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[SOUTH])
                    {
                        if (q.Key == "dvdn")
                        {
                            for (int i = 1; i < v.GetLength(0) - 1; ++i)
                            {
                                v[i, 0] *= flow_corr;
                            }
                        }
                    }
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[WEST])
                    {
                        if (q.Key == "dudn")
                        {
                            for (int j = 1; j < u.GetLength(1) - 1; ++j)
                            {
                                u[0, j] *= flow_corr;
                            }
                        }
                    }
                    foreach (KeyValuePair<string, object> q in (Dictionary<string, object>)this.bc[EAST])
                    {
                        if (q.Key == "dudn")
                        {
                            for (int j = 1; j < u.GetLength(1) - 1; ++j)
                            {
                                u[u.GetLength(0) - 1, j] *= flow_corr;
                            }
                        }
                    }
                }
            }
            double[,] f = new double[p.GetLength(0), p.GetLength(1)];
            for (int i = 1; i < f.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < f.GetLength(1) - 1; ++j)
                {
                    f[i, j] = ((u[i, j] - u[i - 1, j]) / delta[0] +
                        (v[i, j] - v[i, j - 1]) / delta[1]);
                }
            }
            if ((!dirichlet_used) && (outflow_length == 0.0 || ((string)this.default_params["mass_conservation"] != MASS_ADD && (string)this.default_params["mass_conservation"] != MASS_SCALE)))
            {
                double f_sum = 0.0;
                int f_size = 0;
                for (int i = 1; i < f.GetLength(0) - 1; ++i)
                {
                    for (int j = 1; j < f.GetLength(1) - 1; ++j)
                    {
                        f_sum += f[i, j];
                        f_size++;
                    }
                }
                for (int i = 1; i < f.GetLength(0) - 1; ++i)
                {
                    for (int j = 1; j < f.GetLength(1) - 1; ++j)
                    {
                        f[i, j] -= f_sum / f_size;
                    }
                }
            }
            for (int i = 0; i < f.GetLength(0); ++i)
            {
                for (int j = 0; j < f.GetLength(1); ++j)
                {
                    f[i, j] *= (-1.0 / dt);
                }
            }
            return f;
        }

        void Calc_velocity_change(double dt, double[,] rho, double[,] rho_u,
               double[,] rho_v, double[,] mu, double beta, double gamma, out double[,] du, out double[,] dv)
        {
            double[] delta = this.grid.Delta();
            long[,] mask = this.mask;
            double[,] u = this.u;
            double[,] v = this.v;
            double[,] p = this.p;

            du = new double[u.GetLength(0), u.GetLength(1)];
            dv = new double[v.GetLength(0), v.GetLength(1)];

            double[,] mu_off = new double[mu.GetLength(0) - 1, mu.GetLength(1) - 1];
            for (int i = 0; i < mu_off.GetLength(0); ++i)
            {
                for (int j = 0; j < mu_off.GetLength(1); ++j)
                {
                    mu_off[i, j] = 0.25 * (mu[i + 1, j + 1] + mu[i + 1, j] + mu[i, j + 1] + mu[i, j]);
                }
            }
            //Viscosity
            //d/dy(mu*(du/dy+dv/dx))/rho
            for (int i = 0; i < du.GetLength(0); ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] += (mu_off[i, j] * ((u[i, j + 1] - u[i, j]) /
                        (Math.Pow(delta[1], 2)) + (v[i + 1, j] - v[i, j]) / (delta[0] * delta[1])) -
                        mu_off[i, j - 1] * ((u[i, j] - u[i, j - 1]) / (Math.Pow(delta[1], 2)) + (v[i + 1, j - 1] - v[i, j - 1]) / (delta[0] * delta[1]))) / rho_u[i, j];
                }
            }
            //d/dx(2*mu*dx/dx)/rho
            for (int i = 1; i < du.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] += (mu[i + 1, j] * (u[i + 1, j] - u[i, j]) -
                        mu[i, j] * (u[i, j] - u[i - 1, j])) / (0.5 * (Math.Pow(delta[0], 2)) * rho_u[i, j]);
                }
            }

            //boundary (assuming dirichlet bc for pressure)
            //(if another bc is used, these values will be owerriten later in the algorithm
            for (int j = 1; j < du.GetLength(1) - 1; ++j)
            {
                du[0, j] -= (mu[1, j] * (v[1, j] - v[1, j - 1]) - mu[0, j] * (v[0, j] - v[0, j - 1])) /
                    (0.5 * delta[0] * delta[1] * rho_u[0, j]);
                du[du.GetLength(0) - 1, j] -= (mu[mu.GetLength(0) - 1, j] * (v[v.GetLength(0) - 1, j] - v[v.GetLength(0) - 1, j - 1]) - mu[mu.GetLength(0) - 2, j] * v[v.GetLength(0) - 2, j] - v[v.GetLength(0) - 2, j - 1]) /
                    (0.5 * delta[0] * delta[1] * rho_u[rho_u.GetLength(0) - 1, j]);
            }

            //d/dx(mu*(du/dy + dv/dx))/rho
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 0; j < dv.GetLength(1); ++j)
                {
                    dv[i, j] += (mu_off[i, j] * ((v[i + 1, j] - v[i, j]) / (Math.Pow(delta[0], 2)) +
                            (u[i, j + 1] - u[i, j]) / (delta[0] * delta[1])) -
                            mu_off[i - 1, j] * ((v[i, j] - v[i - 1, j]) / (Math.Pow(delta[0], 2)) +
                            (u[i - 1, j + 1] - u[i - 1, j]) / (delta[0] * delta[1]))) / rho_v[i, j];
                }
            }

            // d/dy(2*mu*dv/dy)/rho
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < dv.GetLength(1) - 1; ++j)
                {
                    dv[i, j] += (mu[i, j + 1] * (v[i, j + 1] - v[i, j]) - mu[i, j] *
                        (v[i, j] - v[i, j - 1])) / (0.5 * (Math.Pow(delta[1], 2)) * rho_v[i, j]);
                }
            }

            // boundary
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                dv[i, 0] -= (mu[i, 1] * (u[i, 1] - u[i - 1, 1]) -
                            mu[i, 0] * (u[i, 0] - u[i - 1, 0])) /
                            (0.5 * delta[0] * delta[1] * rho_v[i, 0]);
                dv[i, dv.GetLength(1) - 1] -= (mu[i, mu.GetLength(1) - 1] * (u[i, u.GetLength(1) - 1] - u[i - 1, u.GetLength(1) - 1]) -
                            mu[i, mu.GetLength(1) - 2] * (u[i, u.GetLength(1) - 2] - u[i - 1, u.GetLength(1) - 2])) /
                            (0.5 * delta[0] * delta[1] * rho_v[i, rho_v.GetLength(1) - 1]);
            }
            //Convection
            // d(uu)/dx
            for (int i = 1; i < du.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] -= 0.25 * (Math.Pow(u[i, j] + u[i + 1, j], 2) - Math.Pow(u[i, j] + u[i - 1, j], 2)) / delta[0];
                }
            }

            //boundary 2u*du/dx => -2u*dv/dy
            for (int j = 1; j < du.GetLength(1) - 1; ++j)
            {
                du[0, j] += u[0, j] * ((v[0, j] + v[1, j]) - (v[0, j - 1] + v[1, j - 1])) / delta[1];
                du[du.GetLength(0) - 1, j] += u[du.GetLength(0) - 1, j] * ((v[v.GetLength(0) - 2, j] + v[v.GetLength(0) - 1, j]) - (v[v.GetLength(0) - 2, j - 1] + v[v.GetLength(0) - 1, j - 1])) / delta[1];
            }
            //d(uv)/dy
            for (int i = 0; i < du.GetLength(0); ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] -= 0.25 * ((u[i, j] + u[i, j + 1]) * (v[i, j] + v[i + 1, j]) -
                        (u[i, j] + u[i, j - 1]) * (v[i, j - 1] + v[i + 1, j - 1])) / delta[1];
                }
            }
            //d(vv)/dy
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < dv.GetLength(1) - 1; ++j)
                {
                    dv[i, j] -= 0.25 * (Math.Pow(v[i, j] + v[i, j + 1], 2) - Math.Pow(v[i, j] + v[i, j - 1], 2)) / delta[1];
                }
            }
            //boundary
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                dv[i, 0] += v[i, 0] * ((u[i, 0] + u[i, 1]) - (u[i - 1, 0] + u[i - 1, 1])) / delta[0];
                dv[i, dv.GetLength(1) - 1] += v[i, v.GetLength(1) - 1] * ((u[i, v.GetLength(1) - 2] + u[i, u.GetLength(1) - 1]) -
                    (u[i - 1, u.GetLength(1) - 2] + u[i - 1, u.GetLength(1) - 1])) / delta[0];
            }
            //d(vu)/dx
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 0; j < dv.GetLength(1); ++j)
                {
                    dv[i, j] -= 0.25 * ((v[i, j] + v[i + 1, j]) * (u[i, j] + u[i, j + 1]) -
                        (v[i, j] + v[i - 1, j]) * (u[i - 1, j] + u[i - 1, j + 1])) / delta[0];
                }
            }
            //Upwind deifference for convection (ref.griebel p.29)
            // d(uu)/dx
            for (int i = 1; i < du.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] -= gamma * 0.25 * (Math.Abs(u[i, j] + u[i + 1, j]) *
                         (u[i, j] - u[i + 1, j]) - Math.Abs(u[i - 1, j] + u[i, j]) *
                         (u[i - 1, j] - u[i, j])) / delta[0];
                }
            }

            // d(uv)/dy
            for (int i = 1; i < du.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] -= gamma * 0.25 * ((u[i, j] - u[i, j + 1]) *
                            Math.Abs(v[i, j] + v[i + 1, j]) - (u[i, j - 1] - u[i, j]) * Math.Abs(v[i, j - 1] + v[i + 1, j - 1])) / delta[1];
                }
            }
            // d(vv)/dy
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < dv.GetLength(1) - 1; ++j)
                {
                    dv[i, j] -= gamma * 0.25 * (Math.Abs(v[i, j] + v[i, j + 1]) *
                        (v[i, j] - v[i, j + 1]) - Math.Abs(v[i, j] + v[i, j - 1]) *
                        (v[i, j - 1] - v[i, j])) / delta[1];
                }
            }

            // d(vu)/dx
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < dv.GetLength(1) - 1; ++j)
                {
                    dv[i, j] -= gamma * 0.25 * ((v[i, j] - v[i + 1, j]) *
                        Math.Abs(u[i, j] + u[i, j + 1]) - (v[i - 1, j] - v[i, j]) *
                        Math.Abs(u[i - 1, j] + u[i - 1, j + 1])) / delta[0];
                }
            }
            //gravity
            double[] duv = (double[])this.default_params["gravity"];
            for (int i = 0; i < du.GetLength(0); ++i)
            {
                for (int j = 0; j < du.GetLength(1); ++j)
                {
                    du[i, j] += duv[0];
                }
            }
            for (int i = 0; i < dv.GetLength(0); ++i)
            {
                for (int j = 0; j < dv.GetLength(1); ++j)
                {
                    dv[i, j] += duv[1];
                }
            }
            //surface tension
            if ((bool)this.default_params["two_phase"] && (double)this.default_params["surface_tension_coeff"] != 0.0)
            {
                double[,] F_x = new double[rho_u.GetLength(0) - 2, rho_u.GetLength(1) - 2];
                double[,] F_y = new double[rho_v.GetLength(0) - 2, rho_v.GetLength(1) - 2];
                Surface_tension(out F_x, out F_y);
                for (int i = 0; i < F_x.GetLength(0); ++i)
                {
                    for (int j = 0; j < F_x.GetLength(1); ++j)
                    {
                        F_x[i, j] /= rho_u[i + 1, j + 1];
                    }
                }
                for (int i = 0; i < F_y.GetLength(0); ++i)
                {
                    for (int j = 0; j < F_y.GetLength(1); ++j)
                    {
                        F_y[i, j] /= rho_v[i + 1, j + 1];
                    }
                }
                for (int i = 1; i < du.GetLength(0) - 1; ++i)
                {
                    for (int j = 1; j < du.GetLength(1) - 1; ++j)
                    {
                        du[i, j] += F_x[i - 1, j - 1];
                    }
                }
                for (int i = 1; i < dv.GetLength(0) - 1; ++i)
                {
                    for (int j = 1; j < dv.GetLength(1) - 1; ++j)
                    {
                        dv[i, j] += F_y[i - 1, j - 1];
                    }
                }
            }
            //pressure forces
            //dp/dx
            for (int i = 0; i < du.GetLength(0); ++i)
            {
                for (int j = 1; j < du.GetLength(1) - 1; ++j)
                {
                    du[i, j] -= (beta / delta[0]) * (p[i + 1, j] - p[i, j]) / rho_u[i, j];
                }
            }

            //dp/dy
            for (int i = 1; i < dv.GetLength(0) - 1; ++i)
            {
                for (int j = 0; j < dv.GetLength(1); ++j)
                {
                    dv[i, j] -= (beta / delta[1]) * (p[i, j + 1] - p[i, j]) / rho_v[i, j];
                }
            }
            for (int i = 0; i < du.GetLength(0); ++i)
            {
                for (int j = 0; j < du.GetLength(1); ++j)
                {
                    du[i, j] *= dt;
                }
            }
            for (int i = 0; i < dv.GetLength(0); ++i)
            {
                for (int j = 0; j < dv.GetLength(1); ++j)
                {
                    dv[i, j] *= dt;
                }
            }

        }

        void Surface_tension(out double[,] first, out double[,] second)
        {
            double[,] c = this.c;
            double[] delta = this.grid.Delta();
            long[,] mask = this.mask;
            double[,] u = this.u;
            double[,] v = this.v;

            double[,] F_x;
            double[,] F_y;
            if ((string)this.default_params["curv_method"] == this.CURV_CSF)
            {
                Csf csf = new Csf();
                csf.Surface_tension(c, mask, delta, (double)this.default_params["surface_tension_coeff"], out F_x, out F_y);
                first = new double[F_x.GetLength(0) - 2, F_x.GetLength(1) - 2];
                second = new double[F_y.GetLength(0) - 2, F_y.GetLength(1) - 2];
                for (int i = 0; i < F_x.GetLength(0) - 2; ++i)
                {
                    for (int j = 0; j < F_x.GetLength(1) - 2; ++j)
                    {
                        first[i, j] = F_x[i + 1, j + 1];
                    }
                }
                for (int i = 0; i < F_y.GetLength(0) - 2; ++i)
                {
                    for (int j = 0; j < F_y.GetLength(1) - 2; ++j)
                    {
                        second[i, j] = F_y[i + 1, j + 1];
                    }
                }
                return;
            }
            double[,] kappa_tmp = vof.Curvature(c, mask, delta, (string)this.default_params["curv_method"] == this.CURV_MDAC);
            double[,] kappa = new double[kappa_tmp.GetLength(0) - 2, kappa_tmp.GetLength(1) - 2];
            for (int i = 0; i < kappa.GetLength(0); ++i)
            {
                for (int j = 0; j < kappa.GetLength(1); ++j)
                {
                    kappa[i, j] = kappa_tmp[i + 1, j + 1];
                }
            }
            double[,] w = new double[c.GetLength(0) - 2, c.GetLength(1) - 2];
            for (int i = 0; i < w.GetLength(0); ++i)
            {
                for (int j = 0; j < w.GetLength(1); ++j)
                {
                    w[i, j] = c[i + 1, j + 1] * (1.0 - c[i + 1, j + 1]) + 1e-16;
                }
            }
            double[,] grad_x = new double[c.GetLength(0) - 1, c.GetLength(1) - 2];
            for (int i = 0; i < grad_x.GetLength(0); ++i)
            {
                for (int j = 0; j < grad_x.GetLength(1); ++j)
                {
                    grad_x[i, j] = (c[i + 1, j + 1] - c[i, j + 1]) / delta[0];
                }
            }
            double[,] grad_y = new double[c.GetLength(0) - 2, c.GetLength(1) - 1];
            for (int i = 0; i < grad_y.GetLength(0); ++i)
            {
                for (int j = 0; j < grad_y.GetLength(1); ++j)
                {
                    grad_y[i, j] = (c[i + 1, j + 1] - c[i + 1, j]) / delta[1];
                }
            }
            for (int i = 1; i < grad_x.GetLength(0) - 1; ++i)
            {
                for (int j = 0; j < grad_x.GetLength(1); ++j)
                {
                    grad_x[i, j] *= (mask[i + 1, j + 1] & mask[i, j + 1] & 1);
                }
            }
            for (int i = 0; i < grad_y.GetLength(0); ++i)
            {
                for (int j = 1; j < grad_y.GetLength(1) - 1; ++j)
                {
                    grad_y[i, j] *= (mask[i + 1, j + 1] & mask[i + 1, j] & 1);
                }
            }
            first = new double[kappa.GetLength(0) - 1, kappa.GetLength(1)];
            for (int i = 0; i < first.GetLength(0); ++i)
            {
                for (int j = 0; j < first.GetLength(1); ++j)
                {
                    first[i, j] = (double)this.default_params["surface_tension_coeff"] * (kappa[i + 1, j] * w[i + 1, j] + kappa[i, j] * w[i, j]) * grad_x[i + 1, j] / (w[i + 1, j] + w[i, j]);
                }
            }
            second = new double[kappa.GetLength(0), kappa.GetLength(1) - 1];
            for (int i = 0; i < second.GetLength(0); ++i)
            {
                for (int j = 0; j < second.GetLength(1); ++j)
                {
                    second[i, j] = (double)this.default_params["surface_tension_coeff"] * (kappa[i, j + 1] * w[i, j + 1] + kappa[i, j] * w[i, j]) * grad_y[i, j + 1] / (w[i, j + 1] + w[i, j]);
                }
            }
        }

        void Update_tentative_velocity_bc()
        {
            double[,] u = this.u;
            double[,] v = this.v;

            long[,] u_mask = new long[mask.GetLength(0) - 1, mask.GetLength(1)];
            long[,] v_mask = new long[mask.GetLength(0), mask.GetLength(1) - 1];

            for (int i = 0; i < u_mask.GetLength(0); ++i)
            {
                for (int j = 0; j < u_mask.GetLength(1); ++j)
                {
                    u_mask[i, j] = (this.mask[i, j] | this.mask[i + 1, j]) & 1;
                }
            }

            for (int i = 0; i < v_mask.GetLength(0); ++i)
            {
                for (int j = 0; j < v_mask.GetLength(1); ++j)
                {
                    v_mask[i, j] = (this.mask[i, j] | this.mask[i, j + 1]) & 1;
                }
            }

            for (int i = 0; i < u.GetLength(0); ++i)
            {
                for (int j = 0; j < u.GetLength(1); ++j)
                {
                    u[i, j] *= (this.mask[i, j] & this.mask[i + 1, j] & 1);
                }
            }

            for (int i = 0; i < u.GetLength(0); ++i)
            {
                for (int j = 1; j < u.GetLength(1) - 2; ++j)
                {
                    u[i, j] -= (1 - u_mask[i, j]) * u[i, j + 1];
                }
            }
            for (int i = 0; i < u.GetLength(0); ++i)
            {
                for (int j = 2; j < u.GetLength(1) - 1; ++j)
                {
                    u[i, j] -= (1 - u_mask[i, j]) * u[i, j - 1];
                }
            }

            //zero velocity inside and on the boundary of obstacles
            for (int i = 0; i < v.GetLength(0); ++i)
            {
                for (int j = 0; j < v.GetLength(1); ++j)
                {
                    v[i, j] *= (this.mask[i, j] * this.mask[i, j + 1] & 1);
                }
            }
            for (int i = 1; i < v.GetLength(0) - 2; ++i)
            {
                for (int j = 0; j < v.GetLength(1); ++j)
                {
                    v[i, j] -= (1 - v_mask[i, j]) * v[i + 1, j];
                }
            }
            for (int i = 2; i < v.GetLength(0) - 1; ++i)
            {
                for (int j = 0; j < v.GetLength(1); ++j)
                {
                    v[i, j] -= (1 - v_mask[i, j]) * v[i - 1, j];
                }
            }


            //north boundary
            Dictionary<string, object> b = (Dictionary<string, object>)this.bc[this.NORTH];
            foreach (KeyValuePair<string, object> q in b)
            {
                if (q.Key == "v")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < v.GetLength(0); ++i)
                        {
                            double[] node = grid.Getitem(i - 0.5, v.GetLength(1) - 1);
                            v[i, v.GetLength(1) - 1] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t) * (this.mask[i, mask.GetLength(1) - 2] & 1);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < v.GetLength(0); ++i)
                        {
                            v[i, v.GetLength(1) - 1] = (double)f;
                        }
                    }
                }
            }
            //south boundary
            Dictionary<string, object> s = (Dictionary<string, object>)this.bc[this.SOUTH];
            foreach (KeyValuePair<string, object> q in s)
            {
                if (q.Key == "v")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < v.GetLength(0); ++i)
                        {
                            double[] node = grid.Getitem(i - 0.5, 0);
                            v[i, 0] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t) * (this.mask[i, 1] & 1);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < v.GetLength(0); ++i)
                        {
                            v[i, 0] = (double)f;
                        }
                    }
                }
            }
            //east boundary
            Dictionary<string, object> e = (Dictionary<string, object>)this.bc[this.EAST];
            foreach (KeyValuePair<string, object> q in e)
            {
                if (q.Key == "u")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < u.GetLength(1); ++i)
                        {
                            double[] node = grid.Getitem(u.GetLength(0) - 1, i - 0.5);
                            u[u.GetLength(0) - 1, i] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t) * (this.mask[mask.GetLength(0) - 2, i] & 1);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < v.GetLength(1); ++i)
                        {
                            u[u.GetLength(0) - 1, i] = (double)f;
                        }
                    }
                }
            }
            //west boundary
            Dictionary<string, object> w = (Dictionary<string, object>)this.bc[this.WEST];
            foreach (KeyValuePair<string, object> q in w)
            {
                if (q.Key == "u")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < u.GetLength(1); ++i)
                        {
                            double[] node = grid.Getitem(0, i - 0.5);
                            u[0, i] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t) * (this.mask[1, i] & 1);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < u.GetLength(1); ++i)
                        {
                            u[0, i] = (double)f;
                        }
                    }
                }
            }
        }

        void Update_pressure_bc()
        {
            double[,] p = this.p;

            //north boundary
            Dictionary<string, object> n = (Dictionary<string, object>)this.bc[this.NORTH];
            foreach (KeyValuePair<string, object> q in n)
            {
                if (q.Key == "dpdn")
                {
                    for (int i = 0; i < p.GetLength(0); ++i)
                    {
                        p[i, p.GetLength(1) - 1] = p[i, p.GetLength(1) - 2];
                    }
                }
                else if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(0); ++i)
                        {
                            double[] node = grid.Getitem(i - 0.5, p.GetLength(1) - 2);
                            p[i, p.GetLength(1) - 1] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - p[i, p.GetLength(1) - 2];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < p.GetLength(0); ++i)
                        {
                            p[i, p.GetLength(1) - 1] = 2 * (double)f - p[i, p.GetLength(1) - 2];
                        }
                    }
                }
            }

            //south boundary
            Dictionary<string, object> s = (Dictionary<string, object>)this.bc[this.SOUTH];
            foreach (KeyValuePair<string, object> q in s)
            {
                if (q.Key == "dpdn")
                {
                    for (int i = 0; i < p.GetLength(0); ++i)
                    {
                        p[i, 0] = p[i, 1];
                    }
                }
                else if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(0); ++i)
                        {
                            double[] node = grid.Getitem(i - 0.5, 0);
                            p[i, 0] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - p[i, 1];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < p.GetLength(0); ++i)
                        {
                            p[i, 0] = 2 * (double)f - p[i, 1];
                        }
                    }
                }
            }

            //east boundary
            Dictionary<string, object> e = (Dictionary<string, object>)this.bc[this.EAST];
            foreach (KeyValuePair<string, object> q in e)
            {
                if (q.Key == "dpdn")
                {
                    for (int i = 0; i < p.GetLength(1); ++i)
                    {
                        p[p.GetLength(0) - 1, i] = p[p.GetLength(0) - 2, i];
                    }
                }
                else if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(1); ++i)
                        {
                            double[] node = grid.Getitem(p.GetLength(0) - 2, i - 0.5);
                            p[p.GetLength(0) - 1, i] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - p[p.GetLength(0) - 2, i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < p.GetLength(1); ++i)
                        {
                            p[p.GetLength(0) - 1, i] = 2 * (double)f - p[p.GetLength(0) - 2, i];
                        }
                    }
                }
            }

            //west boundary
            Dictionary<string, object> w = (Dictionary<string, object>)this.bc[this.WEST];
            foreach (KeyValuePair<string, object> q in w)
            {
                if (q.Key == "dpdn")
                {
                    for (int i = 0; i < p.GetLength(1); ++i)
                    {
                        p[0, i] = p[1, i];
                    }
                }
                else if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(1); ++i)
                        {
                            double[] node = grid.Getitem(0, i - 0.5);
                            p[0, i] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - p[1, i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < p.GetLength(1); ++i)
                        {
                            p[0, i] = 2 * (double)f - p[1, i];
                        }
                    }
                }
            }
        }

        void Update_records(double[,] rho)
        {
            if (!(bool)this.default_params["two_phase"])
                return;

            double[,] c = this.c;
            this.record_time.Add(this.t);
            int s = (c.GetLength(0) - 2) * (c.GetLength(1) - 2);
            double gas_sum = 0;
            int gas_cell_sum = 0;
            for (int i = 1; i < c.GetLength(0) - 1; ++i)
            {
                for (int l = 1; l < c.GetLength(1) - 1; ++l)
                {
                    gas_sum += c[i, l];
                    if (c[i, l] > 0.5)
                    {
                        gas_cell_sum++;
                    }
                }
            }
            double liquid_sum = (double)s - gas_sum;
            this.gas_cell_fraction.Add((double)gas_cell_sum / (double)s);
            this.gas_fraction.Add((int)gas_sum / s);

            double[] low = this.grid.min_coord;
            double[] high = this.grid.max_coord;
            double[] delta = this.grid.Delta();
            double[] x = utils.Linspace(low[0] * 0.5 + delta[0], high[0] - 0.5 * delta[0], this.grid.divs[0]);
            double[] y = utils.Linspace(low[1] * 0.5 + delta[1], high[1] - 0.5 * delta[1], this.grid.divs[1]);
            double c_sum_t = 0;
            double c_sum = 0;
            int m = 0;
            int j = 0;
            for (int i = 1; i < c.GetLength(0) - 1; ++i)
            {
                int k = 0;
                for (j = 1; j < c.GetLength(1) - 1; ++j)
                {
                    c_sum_t = c[i, j] * x[m];
                    c_sum += c[i, j] * y[k++];
                }
                m++;
            }
            double x_mc = c_sum_t / gas_sum;
            double y_mc = c_sum / gas_sum;
            double[] mc = new double[2];
            mc[0] = x_mc;
            mc[1] = y_mc;
            this.gas_centre.Add(mc);
            double[,] one_minus_c = new double[c.GetLength(0), c.GetLength(1)];
            for (int i = 0; i < c.GetLength(0); ++i)
            {
                for (j = 0; j < c.GetLength(1); ++j)
                {
                    one_minus_c[i, j] = 1 - c[i, j];
                }
            }
            c_sum_t = 0;
            c_sum = 0;
            m = 0;
            for (int i = 1; i < c.GetLength(0) - 1; ++i)
            {
                int k = 0;
                for (j = 1; j < c.GetLength(1) - 1; ++j)
                {
                    c_sum_t += one_minus_c[i, j] * x[m];
                    c_sum += one_minus_c[i, j] * y[k++];
                }
                m++;
            }
            double[] mc_l = new double[2];
            mc_l[0] = c_sum_t / liquid_sum;
            mc_l[1] = c_sum / liquid_sum;
            this.liquid_centre.Add(mc_l);
        }

        void Update_phi_bc(double[,] phi, double beta)
        {
            double[,] p = this.p;

            //north boundary
            Dictionary<string, object> n = (Dictionary<string, object>)this.bc[this.NORTH];
            foreach (KeyValuePair<string, object> q in n)
            {
                if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(0); ++i)
                        {
                            double[] node = this.grid.Getitem(i - 0.5, p.GetLength(1) - 2);
                            phi[i, phi.GetLength(1) - 1] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - beta * (p[i, p.GetLength(1) - 2] + p[i, p.GetLength(1) - 1]);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < phi.GetLength(0); ++i)
                        {
                            phi[i, phi.GetLength(1) - 1] = 2 * (double)f - beta * (p[i, p.GetLength(1) - 2] + p[i, p.GetLength(1) - 1]);
                        }
                    }
                }
            }

            //south boundary
            Dictionary<string, object> s = (Dictionary<string, object>)this.bc[this.SOUTH];
            foreach (KeyValuePair<string, object> q in s)
            {
                if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(0); ++i)
                        {
                            double[] node = this.grid.Getitem(i - 0.5, 0);
                            phi[i, 0] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - beta * (p[i, 1] + p[i, 0]);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < phi.GetLength(0); ++i)
                        {
                            phi[i, 0] = 2 * (double)f - beta * (p[i, 1] + p[i, 0]);
                        }
                    }
                }
            }

            //east boundary
            Dictionary<string, object> e = (Dictionary<string, object>)this.bc[this.EAST];
            foreach (KeyValuePair<string, object> q in e)
            {
                if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(1); ++i)
                        {
                            double[] node = this.grid.Getitem(p.GetLength(1) - 2, i - 0.5);
                            phi[phi.GetLength(0) - 1, i] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - beta * (p[phi.GetLength(0) - 2, i] + p[phi.GetLength(0) - 1, i]);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < phi.GetLength(1); ++i)
                        {
                            phi[phi.GetLength(0) - 1, i] = 2 * (double)f - beta * (p[phi.GetLength(0) - 2, i] + p[phi.GetLength(0) - 1, i]);
                        }
                    }
                }
            }

            //west boundary
            Dictionary<string, object> w = (Dictionary<string, object>)this.bc[this.WEST];
            foreach (KeyValuePair<string, object> q in w)
            {
                if (q.Key == "p")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < p.GetLength(1); ++i)
                        {
                            double[] node = this.grid.Getitem(0, i - 0.5);
                            phi[0, i] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - beta * (p[1, i] + p[0, i]);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < phi.GetLength(1); ++i)
                        {
                            phi[0, i] = 2 * (double)f - beta * (p[1, i] + p[0, i]);
                        }
                    }
                }
            }
        }



        void Update_velocity_bc()
        {
            double[,] u = this.u;
            double[,] v = this.v;
            double[,] c = this.c;

            Update_tentative_velocity_bc();
            //teh rest
            //north boundary
            Dictionary<string, object> n = (Dictionary<string, object>)this.bc[this.NORTH];
            foreach (KeyValuePair<string, object> q in n)
            {
                if (q.Key == "dvdn")
                {
                    for (int i = 0; i < v.GetLength(0); ++i)
                    {
                        v[i, v.GetLength(1) - 1] = v[i, v.GetLength(1) - 2];
                    }
                }
                if (q.Key == "dudn")
                {
                    for (int i = 0; i < u.GetLength(0); ++i)
                    {
                        u[i, v.GetLength(1) - 1] = u[i, v.GetLength(1) - 2];
                    }
                }
                else if (q.Key == "u")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < u.GetLength(0); ++i)
                        {
                            double[] node = this.grid.Getitem(i, u.GetLength(1) - 2);
                            u[i, u.GetLength(1) - 1] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - u[i, u.GetLength(1) - 2];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < u.GetLength(0); ++i)
                        {
                            u[i, u.GetLength(1) - 1] = 2 * (double)f - u[i, u.GetLength(1) - 2];
                        }
                    }
                }
            }

            //south boundary
            Dictionary<string, object> s = (Dictionary<string, object>)this.bc[this.SOUTH];
            foreach (KeyValuePair<string, object> q in s)
            {
                if (q.Key == "dvdn")
                {
                    for (int i = 0; i < v.GetLength(0); ++i)
                    {
                        v[i, 0] = v[i, 1];
                    }
                }
                if (q.Key == "dudn")
                {
                    for (int i = 0; i < u.GetLength(0); ++i)
                    {
                        u[i, 0] = u[i, 1];
                    }
                }
                else if (q.Key == "u")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < u.GetLength(0); ++i)
                        {
                            double[] node = this.grid.Getitem(i, 0);
                            u[i, 0] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - u[i, 1];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < u.GetLength(0); ++i)
                        {
                            u[i, 0] = 2 * (double)f - u[i, 1];
                        }
                    }
                }
            }

            //east boundary
            Dictionary<string, object> e = (Dictionary<string, object>)this.bc[this.EAST];
            foreach (KeyValuePair<string, object> q in e)
            {
                if (q.Key == "dudn")
                {
                    for (int i = 0; i < u.GetLength(1); ++i)
                    {
                        u[u.GetLength(0) - 1, i] = u[u.GetLength(0) - 2, i];
                    }
                }
                if (q.Key == "dvdn")
                {
                    for (int i = 0; i < v.GetLength(1); ++i)
                    {
                        v[v.GetLength(0) - 1, i] = v[v.GetLength(0) - 2, i]; ;
                    }
                }
                else if (q.Key == "v")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < v.GetLength(1); ++i)
                        {
                            double[] node = this.grid.Getitem(v.GetLength(0) - 2, i);
                            v[v.GetLength(0) - 1, i] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - v[v.GetLength(0) - 2, i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < v.GetLength(1); ++i)
                        {
                            v[v.GetLength(0) - 1, i] = 2 * (double)f - v[v.GetLength(0) - 2, i];
                        }
                    }
                }
            }

            //west boundary
            Dictionary<string, object> w = (Dictionary<string, object>)this.bc[this.WEST];
            foreach (KeyValuePair<string, object> q in w)
            {
                if (q.Key == "dudn")
                {
                    for (int i = 0; i < u.GetLength(1); ++i)
                    {
                        u[0, i] = u[1, i];
                    }
                }
                if (q.Key == "dvdn")
                {
                    for (int i = 0; i < v.GetLength(1); ++i)
                    {
                        v[0, i] = v[1, i]; ;
                    }
                }
                else if (q.Key == "v")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < v.GetLength(1); ++i)
                        {
                            double[] node = this.grid.Getitem(0, i);
                            v[0, i] = 2 * ((Func<double, double, double, double>)f)(node[0], node[1], this.t) - v[1, i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < v.GetLength(1); ++i)
                        {
                            v[0, i] = 2 * (double)f - v[1, i];
                        }
                    }
                }
            }
        }

        void Update_level_bc()
        {
            double[,] c = this.c;
            //north boundary
            Dictionary<string, object> n = (Dictionary<string, object>)this.bc[this.NORTH];
            foreach (KeyValuePair<string, object> q in n)
            {
                if (q.Key == "c")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < c.GetLength(0); ++i)
                        {
                            double[] node = this.grid.Getitem(i - 0.5, c.GetLength(1) - 1.5);
                            c[i, c.GetLength(1) - 1] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < c.GetLength(0); ++i)
                        {
                            c[i, c.GetLength(1) - 1] = (double)f;
                        }
                    }
                }
            }

            //south boundary
            Dictionary<string, object> s = (Dictionary<string, object>)this.bc[this.SOUTH];
            foreach (KeyValuePair<string, object> q in s)
            {
                if (q.Key == "c")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < c.GetLength(0); ++i)
                        {
                            double[] node = this.grid.Getitem(i - 0.5, -0.5);
                            c[i, 0] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < c.GetLength(0); ++i)
                        {
                            c[i, 0] = (double)f;
                        }
                    }
                }
            }
            //east boundary
            Dictionary<string, object> e = (Dictionary<string, object>)this.bc[this.EAST];
            foreach (KeyValuePair<string, object> q in e)
            {
                if (q.Key == "c")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < c.GetLength(1); ++i)
                        {
                            double[] node = this.grid.Getitem(c.GetLength(0) - 1.5, i - 0.5);
                            c[c.GetLength(0) - 1, i] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < c.GetLength(1); ++i)
                        {
                            c[c.GetLength(0) - 1, i] = (double)f;
                        }
                    }
                }
            }
            //west boundary
            Dictionary<string, object> w = (Dictionary<string, object>)this.bc[this.WEST];
            foreach (KeyValuePair<string, object> q in w)
            {
                if (q.Key == "c")
                {
                    object f = q.Value;
                    if (f.GetType() == typeof(Func<double, double, double, double>))
                    {
                        for (int i = 0; i < c.GetLength(1); ++i)
                        {
                            double[] node = this.grid.Getitem(-0.5, i - 0.5);
                            c[0, i] = ((Func<double, double, double, double>)f)(node[0], node[1], this.t);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < c.GetLength(1); ++i)
                        {
                            c[0, i] = (double)f;
                        }
                    }
                }
            }
        }

        //retrieving
        double Interpolate_u(double[] pos)
        {
            for (int i = 0; i < 2; ++i)
            {
                pos[i] = (pos[i] - this.grid.min_coord[i]) / (this.grid.max_coord[i] - this.grid.min_coord[i]);
            }
            pos[0] = pos[0] * this.grid.divs[0];
            pos[1] = pos[1] * this.grid.divs[1] + 0.5;
            int[] low = new int[pos.Length];
            low[0] = (int)pos[0];
            low[1] = (int)pos[1];
            if (low[0] < 0)
                low[0] = 0;
            if (low[0] > this.u.GetLength(0) - 2)
                low[0] = this.u.GetLength(0) - 2;
            if (low[1] < 0)
                low[1] = 0;
            if (low[1] > this.u.GetLength(1) - 2)
                low[1] = this.u.GetLength(1) - 2;
            double[] frac = new double[2];
            for (int i = 0; i < 2; ++i)
            {
                frac[i] = pos[i] - low[i];
            }
            double left = utils.Lerp(this.u[low[0], low[1]], this.u[low[0], low[1] + 1], frac[1]);
            double right = utils.Lerp(this.u[low[0] + 1, low[1]], this.u[low[0] + 1, low[1] + 1], frac[1]);
            return utils.Lerp(left, right, frac[0]);
        }

        double Interpolate_v(double[] pos)
        {
            for (int i = 0; i < 2; ++i)
            {
                pos[i] = (pos[i] - this.grid.min_coord[i]) / (this.grid.max_coord[i] - this.grid.min_coord[i]);
            }
            pos[0] = pos[0] * this.grid.divs[0] + 0.5;
            pos[1] = pos[1] * this.grid.divs[1];
            int[] low = new int[pos.Length];
            low[0] = (int)pos[0];
            low[1] = (int)pos[1];
            if (low[0] < 0)
                low[0] = 0;
            if (low[0] > this.v.GetLength(0) - 2)
                low[0] = this.v.GetLength(0) - 2;
            if (low[1] < 0)
                low[1] = 0;
            if (low[1] > this.v.GetLength(1) - 2)
                low[1] = this.v.GetLength(1) - 2;
            double[] frac = new double[2];
            for (int i = 0; i < 2; ++i)
            {
                frac[i] = pos[i] - low[i];
            }
            double left = utils.Lerp(this.v[low[0], low[1]], this.v[low[0], low[1] + 1], frac[1]);
            double right = utils.Lerp(this.v[low[0] + 1, low[1]], this.v[low[0] + 1, low[1] + 1], frac[1]);
            return utils.Lerp(left, right, frac[0]);
        }

        double[] Interpolate_velocity(double[] pos)
        {
            double[] int_pos = new double[2];
            int_pos[1] = Interpolate_u(pos);
            int_pos[2] = Interpolate_u(pos);
            return (int_pos);
        }

        void Avg_velocity(out double[,] u_avg, out double[,] v_avg, bool middle = true)
        {

            if (middle)
            {
                u_avg = new double[u.GetLength(0) - 1, u.GetLength(1) - 2];
                v_avg = new double[v.GetLength(0) - 2, v.GetLength(1) - 1];
                for (int i = 0; i < u_avg.GetLength(0); ++i)
                {
                    for (int j = 0; j < u_avg.GetLength(1); ++j)
                    {
                        u_avg[i, 0] = 0.5 * (this.u[i, j + 1] + this.u[i + 1, j + 1]);
                    }
                }
                for (int i = 0; i < v_avg.GetLength(0); ++i)
                {
                    for (int j = 0; j < v_avg.GetLength(1); ++j)
                    {
                        v_avg[0, i] = 0.5 * (this.v[i + 1, j] + this.v[i + 1, j + 1]);
                    }
                }

            }
            else
            {
                u_avg = new double[u.GetLength(0), u.GetLength(1) - 1];
                v_avg = new double[v.GetLength(0) - 1, v.GetLength(1)];
                for (int i = 0; i < u_avg.GetLength(0); ++i)
                {
                    for (int j = 0; j < u_avg.GetLength(1); ++j)
                    {
                        u_avg[i, j] = 0.5 * (this.u[i, j + 1] + this.u[i, j]);
                    }
                }
                for (int i = 0; i < v_avg.GetLength(0); ++i)
                {
                    for (int j = 0; j < v_avg.GetLength(1); ++j)
                    {
                        v_avg[j, i] = 0.5 * (this.v[i + 1, j] + this.v[i, j]);
                    }
                }
            }
        }

        List<double[]> Mass_centre()
        {
            List<double[]> cm = new List<double[]>();
            for (int i = 0; i < gas_centre.Count; ++i)
            {
                double[] G = new double[2];
                G = gas_centre[i];
                double[] L = new double[2];
                L = liquid_centre[i];
                double[] tcm = new double[2];
                tcm[0] = (double)default_params["rho_gas"] * G[0] + (double)this.default_params["rho_liquid"] * L[0];
                tcm[1] = (double)default_params["rho_gas"] * G[1] + (double)this.default_params["rho_liquid"] * L[1];
                cm.Add(tcm);
            }
            return (cm);
        }

        double[,] Diveregence()
        {
            double[] delta = this.grid.Delta();
            double[,] f = new double[this.p.GetLength(0), this.p.GetLength(1)];
            double[,] u = this.u;
            double[,] v = this.v;
            for (int i = 1; i < f.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < f.GetLength(1) - 1; ++j)
                {
                    f[i, j] = (u[i, j] - u[i - 1, j]) / delta[0] +
                        (v[i, j] - v[i, j - 1]) / delta[1];
                }
            }
            return (f);
        }

        double[,] Vorticity()
        {
            double[] delta = this.grid.Delta();
            double[,] u = this.u;
            double[,] v = this.v;
            double[,] f = new double[v.GetLength(0) - 1, v.GetLength(1)];
            for (int i = 0; i < f.GetLength(0); ++i)
            {
                for (int j = 0; j < f.GetLength(1); ++j)
                {
                    f[i, j] = (v[i + 1, j] - v[i, j]) / delta[0] -
                        (u[i, j + 1] - u[i, j]) / delta[1];
                }
            }
            return (f);
        }

        double[,] Stream_function()
        {
            double[,] u = this.u;
            double[,] v = this.v;
            double[,] phi = new double[this.grid.Shape()[0], this.grid.Shape()[1]];
            long[,] mask = new long[this.mask.GetLength(0) - 1, this.mask.GetLength(1) - 2];
            for (int i = 0; i < mask.GetLength(0); ++i)
            {
                for (int j = 0; j < mask.GetLength(1); ++j)
                {
                    mask[i, j] = (this.mask[i, j + 1] | this.mask[i + 1, j + 1]) & 1;
                }
            }
            for (int i = 1; i < v.GetLength(0) - 1; ++i)
            {
                v[i, 0] *= this.grid.Delta()[0];
            }
            for (int i = 2; i < v.GetLength(0) - 1; ++i)
            {
                v[i, 0] += v[i - 1, 0];
            }
            for (int i = 1; i < phi.GetLength(0); ++i)
            {
                phi[i, 0] = -v[i, 0];
            }
            for (int i = 0; i < phi.GetLength(0); ++i)
            {
                for (int j = 1; j < phi.GetLength(1); ++j)
                {
                    phi[i, j] = mask[i, j] * u[i, j] * this.grid.Delta()[1];
                }
            }
            for (int i = 0; i < phi.GetLength(0); ++i)
            {
                for (int j = 1; j < phi.GetLength(1); ++j)
                {
                    phi[i, j] += phi[i, j - 1];
                }
            }
            return (phi);
        }

        public double Find_suitable_dt(double safety_factor = 0.95)
        {
            double[] delta = this.grid.Delta();
            double[,] u = new double[this.u.GetLength(0), this.u.GetLength(1)];
            double[,] v = new double[this.v.GetLength(0), this.v.GetLength(1)];
            for (int i = 0; i < u.GetLength(0); ++i)
            {
                for (int j = 0; j < u.GetLength(1); ++j)
                {
                    int q = this.u[i, j] == 0.0 ? 1 : 0;
                    u[i, j] = Math.Abs(this.u[i, j] + 0.000001 * q);
                }
            }
            for (int i = 0; i < v.GetLength(0); ++i)
            {
                for (int j = 0; j < v.GetLength(1); ++j)
                {
                    int q = this.v[i, j] == 0.0 ? 1 : 0;
                    v[i, j] = Math.Abs(this.v[i, j] + 0.000001 * q);
                }
            }
            double u_max = u[0, 0];
            for (int i = 0; i < u.GetLength(0); ++i)
            {
                for (int j = 0; j < u.GetLength(1); ++j)
                {
                    if (u_max < u[i, j])
                        u_max = u[i, j];
                }
            }
            double v_max = v[0, 0];
            for (int i = 0; i < v.GetLength(0); ++i)
            {
                for (int j = 0; j < v.GetLength(1); ++j)
                {
                    if (v_max < v[i, j])
                        v_max = v[i, j];
                }
            }
            double st_delta = 0;
            for (int i = 0; i < 2; ++i)
            {
                st_delta += 1 / (delta[i] * delta[i]);
            }
            if ((bool)this.default_params["two_phase"])
            {
                double[] nus = new double[2];
                nus[0] = (double)this.default_params["mu_liquid"] / (double)this.default_params["rho_liquid"];
                nus[1] = (double)this.default_params["mu_gas"] / (double)this.default_params["rho_gas"];


                double dt1 = 0.5 * Math.Min(delta[0] / u_max, delta[1] / v_max);

                double nus_max = nus[0] > nus[1] ? nus[0] : nus[1];
                double dt2 = 0.5 / (nus_max * st_delta);
                double nus_min = nus[0] < nus[1] ? nus[0] : nus[1];
                double dt3 = 2.0 * nus_min / Math.Pow(Math.Max(u_max, v_max), 2);
                if ((double)this.default_params["surface_tension_coeff"] > 0.0)
                {
                    double delta_min = delta[0] < delta[1] ? delta[0] : delta[1];
                    double dt4 = Math.Sqrt((((double)this.default_params["rho_liquid"] + (double)this.default_params["rho_gas"]) * Math.Pow(delta_min, 3) / (4.0 * Math.PI * (double)this.default_params["surface_tension_coeff"])));
                    double[] dt = new double[] { dt1, dt2, dt3, dt4 };
                    double dt_min = dt[0];
                    for (int i = 1; i < dt.Length; ++i)
                    {
                        dt_min = dt_min < dt[i] ? dt_min : dt[i];
                    }
                    return (dt_min * safety_factor);
                }
                else
                {
                    double[] dt = new double[] { dt1, dt2, dt3 };
                    double dt_min = dt[0];
                    for (int i = 1; i < dt.Length; ++i)
                    {
                        if (dt_min > dt[i])
                            dt_min = dt[i];
                    }
                    return (dt_min * safety_factor);
                }
            }
            else
            {
                double nu = (double)this.default_params["mu_liquid"] / (double)this.default_params["rho_liquid"];
                double dt1 = Math.Min(delta[0] / u_max, delta[1] / v_max);
                if ((bool)this.default_params["use_dye"])
                {
                    dt1 *= 0.5;
                }
                double dt2 = 0.5 / (nu * st_delta);
                double dt3 = 2.0 * nu / Math.Pow(Math.Max(u_max, v_max), 2);
                double[] dt = new double[] { dt1, dt2, dt3 };
                double dt_min = dt[0];
                for (int i = 1; i < dt.Length; ++i)
                {
                    if (dt_min > dt[i])
                        dt_min = dt[i];
                }
                return (dt_min * safety_factor);
            }
        }
    }
}

