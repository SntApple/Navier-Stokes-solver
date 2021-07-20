using System;
using System.Diagnostics;

namespace ConsoleApp2
{
    public class Csf
    {
        static Utils u = new Utils();
        public double[,] Smooth(double[,] phi, long[,] mask, double[] delta)
        {
            double[,] stencil = new double[3, 3] { { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 } };

            long[,] m = new long[mask.GetLength(0), mask.GetLength(1)];
            for (int i = 0; i < m.GetLength(0); ++i)
            {
                for (int j = 0; j < m.GetLength(1); ++j)
                {
                    m[i, j] = mask[i, j] & 1;
                }
            }

            int[] mid = new int[2];

            mid[0] = stencil.GetLength(0) / 2;
            mid[1] = stencil.GetLength(1) / 2;
            double[,] sum_phi = new double[phi.GetLength(0), phi.GetLength(1)];
            double[,] sum_w = (double[,])sum_phi.Clone();
            for (int i = 0; i < stencil.GetLength(0); ++i)
            {
                int x_offs = i - mid[0];
                int x1 = Math.Max(0, x_offs) + 1;
                int x2 = phi.GetLength(0) + Math.Min(0, x_offs) - 1;
                for (int j = 0; j < stencil.GetLength(1); ++j)
                {
                    int y_offs = j - mid[1];
                    int y1 = Math.Max(0, y_offs) + 1;
                    int y2 = phi.GetLength(1) + Math.Min(0, y_offs) - 1;
                    for(int k = x1-x_offs; k < x2-x_offs; ++k)
                    {
                        for(int l = y1-y_offs; l < y2-y_offs; ++l)
                        {
                            sum_phi[k, l] += stencil[i, j] * phi[k + x_offs, l + y_offs] * m[k + x_offs, l + y_offs];
                            sum_w[k, l] += stencil[i, j] * m[k + x_offs, l + y_offs];
                        }
                    }
                }
            }
            double[,] answer = new double[sum_phi.GetLength(0), sum_phi.GetLength(1)];
            for (int i = 0; i < answer.GetLength(0); ++i)
            {
                for (int j = 0; j < answer.GetLength(1); ++j)
                {
                    int zero = sum_w[i, j] == 0.0 ? 1 : 0;
                    answer[i, j] = sum_phi[i, j] / (sum_w[i, j] + zero);
                }
            }
            return answer;
        }

        void Gradient(double[,] phi, long[,] mask, double[] delta, out double[,] gx, out double[,] gy)
        {
            double dx = delta[0];
            double dy = delta[1];
            gx = new double[phi.GetLength(0) - 1, phi.GetLength(1)];
            gy = new double[phi.GetLength(0), phi.GetLength(1) - 1];
            for (int i = 1; i < gx.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < gx.GetLength(1) - 1; ++j)
                {
                    gx[i, j] = (phi[i + 1, j] - phi[i, j]) / dx;
                    gx[i, j] *= mask[i, j] & mask[i + 1, j] & 1;
                }
            }
            for(int i = 1; i < gy.GetLength(0) - 1; ++i)
            {
                for(int j = 1; j < gy.GetLength(1) - 1; ++j)
                {
                    gy[i, j] = (phi[i, j + 1] - phi[i, j]) / dy;
                    gy[i, j] *= mask[i, j] & mask[i, j + 1] & 1;
                }
            }
        }

        double[,] Curvature(double[,] gx, double[,] gy, double[] delta)
        {
            double dx = delta[0];
            double dy = delta[1];
            double[,] div = new double[gy.GetLength(0), gx.GetLength(1)];

            double[,] vgx = new double[gx.GetLength(0), gx.GetLength(1) - 1];
            for (int i = 0; i < vgx.GetLength(0); ++i)
            {
                for (int j = 0; j < vgx.GetLength(1); ++j)
                {
                    vgx[i, j] = 0.5 * (gx[i, j + 1] + gx[i, j]);
                }
            }
            double[,] vgy = new double[gy.GetLength(0) - 1, gy.GetLength(1)];
            for (int i = 0; i < vgy.GetLength(0); ++i)
            {
                for (int j = 0; j < vgy.GetLength(1); ++j)
                {
                    vgy[i, j] = 0.5 * (gy[i + 1, j] + gy[i, j]);
                }
            }

            double[,] len_vg = new double[vgx.GetLength(0), vgx.GetLength(1)];
            for (int i = 0; i < len_vg.GetLength(0); ++i)
            {
                for (int j = 0; j < len_vg.GetLength(1); ++j)
                {
                    len_vg[i, j] = Math.Sqrt(vgx[i, j] * vgx[i, j] + vgy[i, j] * vgy[i, j]);
                }
            }
            double[,] div_g = new double[vgx.GetLength(0) - 1, vgx.GetLength(1) - 1];
            for (int i = 0; i < div_g.GetLength(0); ++i)
            {
                for (int j = 0; j < div_g.GetLength(1); ++j)
                {
                    div_g[i, j] = (vgx[i + 1, j + 1] + vgx[i + 1, j] - vgx[i, j + 1] - vgx[i, j]) / (2.0 * dx)
                        + (vgy[i + 1, j + 1] - vgy[i + 1, j] + vgy[i, j + 1] - vgy[i, j]) / (2.0 * dy);
                }
            }
            double[,] cgx = new double[vgx.GetLength(0) - 1, vgx.GetLength(1) - 1];
            for (int i = 0; i < cgx.GetLength(0); ++i)
            {
                for (int j = 0; j < cgx.GetLength(1); ++j)
                {
                    cgx[i, j] = 0.25 * (vgx[i + 1, j + 1] + vgx[i + 1, j] + vgx[i, j + 1] + vgx[i, j]);
                }
            }
            double[,] cgy = new double[vgy.GetLength(0) - 1, vgy.GetLength(1) - 1];
            for (int i = 0; i < cgy.GetLength(0); ++i)
            {
                for (int j = 0; j < cgy.GetLength(1); ++j)
                {
                    cgy[i, j] = 0.25 * (vgy[i + 1, j + 1] + vgy[i + 1, j] + vgy[i, j + 1] + vgy[i, j]);
                }
            }
            double[,] len_cg = new double[cgx.GetLength(0), cgx.GetLength(1)];
            for (int i = 0; i < len_cg.GetLength(0); ++i)
            {
                for (int j = 0; j < len_cg.GetLength(1); ++j)
                {
                    len_cg[i, j] = Math.Sqrt(cgx[i, j] * cgx[i, j] + cgy[i, j] * cgy[i, j]);
                    if (len_cg[i, j] == 0.0)
                        len_cg[i, j]++;
                }
            }
           double[,] n_dot_grad = new double[len_vg.GetLength(0) - 1, len_vg.GetLength(1) - 1];
            for (int i = 0; i < n_dot_grad.GetLength(0); ++i)
            {
                for (int j = 0; j < n_dot_grad.GetLength(1); ++j)
                {
                    n_dot_grad[i, j] = (len_vg[i + 1, j + 1] + len_vg[i + 1, j] - len_vg[i, j + 1] - len_vg[i, j]) / (2.0 * dx)
                        * (cgx[i, j] / len_cg[i, j]) + (len_vg[i + 1, j + 1] - len_vg[i + 1, j] + len_vg[i, j + 1] - len_vg[i, j])
                        / (2.0 * dy) * (cgy[i, j] / len_cg[i, j]);
                }
            }
            for (int i = 1; i < div.GetLength(0) - 1; ++i)
            {
                for (int j = 1; j < div.GetLength(1) - 1; ++j)
                {
                    div[i, j] = (n_dot_grad[i - 1, j - 1] - div_g[i - 1, j - 1]) / len_cg[i - 1, j - 1];
                }
            }
            return (div);
        }

        public void Surface_tension(double[,] phi, long[,] mask, double[] delta, double sigma, out double[,] Fx, out double[,] Fy)
        {
            double[,] phi_kappa = Smooth(phi, mask, delta);
            Gradient(phi_kappa, mask, delta, out double[,] gx, out double[,] gy);
            double[,] kappa = Curvature(gx, gy, delta);
            Gradient(phi, mask, delta, out gx, out gy);
            Fx = new double[kappa.GetLength(0) - 1, kappa.GetLength(1)];
            Fy = new double[kappa.GetLength(0), kappa.GetLength(1) - 1];
            for (int i = 0; i < Fx.GetLength(0); ++i)
            {
                for (int j = 0; j < Fx.GetLength(1); ++j)
                {
                    Fx[i, j] = 0.5 * sigma * gx[i, j] * (kappa[i + 1, j] + kappa[i, j]);
                }
            }

            for (int i = 0; i < Fy.GetLength(0); ++i)
            {
                for (int j = 0; j < Fy.GetLength(1); ++j)
                {
                    Fy[i, j] = 0.5 * sigma * gy[i, j] * (kappa[i, j + 1] + kappa[i, j]);
                }
            }
        }
    }
}

