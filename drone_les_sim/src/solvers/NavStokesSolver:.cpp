void NavStokesSolver::computePredictorVelocity() {
    const double dx = grid_.dx();
    const double dy = grid_.dy();
    const int nx = grid_.nx();
    const int ny = grid_.ny();
  

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {

            // =================================================================
            // 1. COMPOSANTE X (u_star)
            // =================================================================

            // A. Advection de u : d(u^2)/dx + d(uv)/dy
            double u_E = (u_(i, j) + u_(i+1, j)) / 2.0;
            double u_W = (u_(i-1, j) + u_(i, j)) / 2.0;
            double adv_EW = (u_E * u_E - u_W * u_W) / dx;

            double u_N = (u_(i, j) + u_(i, j+1)) / 2.0;
            double v_N = (v_(i, j) + v_(i, j+1)) / 2.0;
            double uv_N = u_N * v_N;

            double u_S = (u_(i, j-1) + u_(i, j)) / 2.0;
            double v_S = (v_(i, j-1) + v_(i, j)) / 2.0;
            double uv_S = u_S * v_S;

            double adv_NS = (uv_N - uv_S) / dy;
            double adv_x = adv_EW + adv_NS;

            // B. Viscosité pour u : d/dx(2*nu*du/dx) + d/dy(nu*(du/dy + dv/dx))
            double nu_E = (nu_t_(i+1, j) + nu_t_(i, j)) / 2.0 + nu_mol_;
            double nu_W = (nu_t_(i, j) + nu_t_(i-1, j)) / 2.0 + nu_mol_;
            double v_x1 = (2.0 * nu_E * (u_(i+1, j) - u_(i, j)) - 2.0 * nu_W * (u_(i, j) - u_(i-1, j))) / (dx * dx);

            double nu_N = (nu_t_(i, j+1) + nu_t_(i, j)) / 2.0 + nu_mol_;
            double nu_S = (nu_t_(i, j) + nu_t_(i, j-1)) / 2.0 + nu_mol_;

            double dvdx_N = (v_(i+1, j+1) - v_(i-1, j+1) + v_(i+1, j) - v_(i-1, j)) / (4.0 * dx);
            double dvdx_S = (v_(i+1, j) - v_(i-1, j) + v_(i+1, j-1) - v_(i-1, j-1)) / (4.0 * dx);

            double v_x2_N = nu_N * ((u_(i, j+1) - u_(i, j)) / dy + dvdx_N);
            double v_x2_S = nu_S * ((u_(i, j) - u_(i, j-1)) / dy + dvdx_S);
            double v_x2 = (v_x2_N - v_x2_S) / dy;

            double diff_x = v_x1 + v_x2;

            // C. Mise à jour u_star
            double Rx = -adv_x + diff_x;
            u_star_(i, j) = u_(i, j) + dt_ * Rx;


            // =================================================================
            // 2. COMPOSANTE Y (v_star)
            // =================================================================

            // A. Advection de v : d(uv)/dx + d(v^2)/dy
            double u_E_v = (u_(i, j) + u_(i, j+1)) / 2.0;
            double v_E_v = (v_(i, j) + v_(i+1, j)) / 2.0;
            double vu_E = u_E_v * v_E_v;

            double u_W_v = (u_(i-1, j) + u_(i-1, j+1)) / 2.0;
            double v_W_v = (v_(i-1, j) + v_(i, j)) / 2.0;
            double vu_W = u_W_v * v_W_v;

            double adv_y_EW = (vu_E - vu_W) / dx;

            double v_N_v = (v_(i, j) + v_(i, j+1)) / 2.0;
            double v_S_v = (v_(i, j-1) + v_(i, j)) / 2.0;
            double adv_y_NS = (v_N_v * v_N_v - v_S_v * v_S_v) / dy;

            double adv_y = adv_y_EW + adv_y_NS;

            // B. Viscosité pour v : d/dx(nu*(du/dy + dv/dx)) + d/dy(2*nu*dv/dy)
            double dudy_E = (u_(i+1, j+1) - u_(i+1, j-1) + u_(i, j+1) - u_(i, j-1)) / (4.0 * dy);
            double dudy_W = (u_(i, j+1) - u_(i, j-1) + u_(i-1, j+1) - u_(i-1, j-1)) / (4.0 * dy);

            double v_y1_E = nu_E * (dudy_E + (v_(i+1, j) - v_(i, j)) / dx);
            double v_y1_W = nu_W * (dudy_W + (v_(i, j) - v_(i-1, j)) / dx);
            double v_y1 = (v_y1_E - v_y1_W) / dx;

            double v_y2 = (2.0 * nu_N * (v_(i, j+1) - v_(i, j)) - 2.0 * nu_S * (v_(i, j) - v_(i, j-1))) / (dy * dy);

            double diff_y = v_y1 + v_y2;

            // C. Mise à jour v_star
            double Ry = -adv_y + diff_y;
            v_star_(i, j) = v_(i, j) + dt_ * Ry;
        }
    }
}