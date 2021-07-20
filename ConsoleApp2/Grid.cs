namespace ConsoleApp2
{
    class Grid
    {
        public double[] min_coord;
        public double[] max_coord;
        public int[] divs;
        public Grid(double[] min_coord, double[] max_coord, int[] divs)
        {
            this.min_coord = min_coord;
            this.max_coord = max_coord;
            this.divs = divs;
        }

        public double[] Getitem(double index1, double index2)
        {
            double[] index = new double[2] { index1, index2 };
            double[] i = new double[2];
            for (int j = 0; j < 2; ++j)
            {
                i[j] = this.min_coord[j] + index[j] / this.divs[j] *
                    (this.max_coord[j] - this.min_coord[j]);
            }
            return (i);
        }
        public int[] Get_shape()
        {
            int[] q = new int[2];
            for (int i = 0; i < 2; ++i)
            {
                q[i] = this.divs[i] + 1;
            }
            return (q);
        }

        public double[] Get_delta()
        {
            double[] l = new double[2];
            for (int i = 0; i < 2; ++i)
            {
                l[i] = (this.max_coord[i] - this.min_coord[i])
                    / this.divs[i];
            }
            return (l);
        }

        double[,] Get_coordinate(int i)
        {
            Utils u = new Utils();
            double[,] x = new double[this.divs[0] + 1, this.divs[1] + 1];
            double[] lin = u.Linspace(this.min_coord[i], this.max_coord[i], this.divs[i] + 1);
            if (i == 0)
            {
                for (int l = 0; l < divs[i] + 1; ++l)
                {
                    for (int m = 0; m < divs[i] + 1; ++m)
                    {
                        x[l, m] = lin[l];
                    }
                }
            }
            if (i == 1)
            {
                for (int l = 0; l < divs[i] + 1; ++l)
                {
                    for (int m = 0; m < divs[i] + 1; ++m)
                    {
                        x[l, m] = lin[m];
                    }
                }
            }
            return x;
        }


        public int[] Shape()
        {
            return (Get_shape());
        }
        public double[] Delta()
        {
            return (Get_delta());
        }
        public double[,] Get_y()
        {
            return (Get_coordinate(1));
        }
        public double[,] Get_z()
        {
            return (Get_coordinate(2));
        }
    }
}
