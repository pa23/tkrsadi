/*
  tkrsadi
  Turbocharger parameters calculation.

  File: calc.cpp

  Copyright (C) 2021 Artem Petrov <pa23666@yandex.ru>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "calc.hpp"
#include "const.hpp"
#include "prgid.hpp"
#include "conf.hpp"
#include "auxf.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <iomanip>

using std::cout;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::shared_ptr;
using std::istreambuf_iterator;
using std::setw;
using std::setprecision;
using std::fixed;
using std::setfill;

Calc::Calc(const shared_ptr<Conf> &conf) {
    m_conf = conf;
}

bool Calc::calculate() {

    const double Ne       = m_conf->val_Ne();
    const double ge       = m_conf->val_ge();
    const double alpha    = m_conf->val_alpha();
    const double phip     = m_conf->val_phip();
    const double pik      = m_conf->val_pik();
    const double itkr     = m_conf->val_itkr();
    const double p0_z     = m_conf->val_p0() / 1000.0;
    const double t0_z     = m_conf->val_t0() + 273;
    const double sigma_in = m_conf->val_sigma_in();
    const double ca       = m_conf->val_ca();
    const double etakn    = m_conf->val_etakn();
    const double phi      = m_conf->val_phi();
    const double zeta_in  = m_conf->val_zeta_in();
    const double D0D1     = m_conf->val_D0D1();
    const double D1D2     = m_conf->val_D1D2();
    const double ia       = m_conf->val_ia();
    const double tau1     = m_conf->val_tau1();
    const double zeta1    = m_conf->val_zeta1();
    const double zeta2    = m_conf->val_zeta2();
    const double kc2r     = m_conf->val_kc2r();
    const double alphap   = m_conf->val_alphap();
    const double zk       = m_conf->val_zk();
    const double tau2     = m_conf->val_tau2();
    const double kb3      = m_conf->val_kb3();
    const double kD3      = m_conf->val_kD3();
    const double zetad    = m_conf->val_zetad();
    const double kD4      = m_conf->val_kD4();
    const double kb4      = m_conf->val_kb4();
    const double kalpha4  = m_conf->val_kalpha4();
    const double zd       = m_conf->val_zd();
    const double n4       = m_conf->val_n4();
    const double tau3     = m_conf->val_tau3();
    const double tau4     = m_conf->val_tau4();
    const double nul      = m_conf->val_nul();
    const double kck      = m_conf->val_kck();

    const double pr_z     = m_conf->val_pr_z() / 1000.0;
    const double tr_z     = m_conf->val_tr_z() + 273.0;
    const double sigmay   = m_conf->val_sigmay();
    const double kD1t     = m_conf->val_kD1t();
    const double etat     = m_conf->val_etat();
    const double kD2t     = m_conf->val_kD2t();
    const double kDbt     = m_conf->val_kDbt();
    const double rhot     = m_conf->val_rhot();
    const double phit     = m_conf->val_phit();
    const double alpha1t  = m_conf->val_alpha1t();
    const double psi      = m_conf->val_psi();
    const double delta    = m_conf->val_delta() / 1000.0;
    const double kNrd     = m_conf->val_kNrd();
    const double etamt    = m_conf->val_etamt();

    // radial flow compressor calculation

    m_Gk = 4.028 * 0.001 * Ne * ge / 1000.0 * alpha * phip;
    m_ta_z = t0_z;
    m_pa_z = p0_z * sigma_in;
    m_ta = m_ta_z - pow(ca, 2) / 2.0 / CPair;
    m_pa = m_pa_z * pow(m_ta / m_ta_z, Kair / (Kair - 1.0));
    m_La = CPair * m_ta * (pow(pik, (Kair - 1.0) / Kair) - 1.0);
    m_u2 = sqrt(m_La / etakn);
    m_c1a = m_u2 * phi;
    m_c1 = m_c1a;
    m_t1 = m_ta + (pow(ca, 2.0) - pow(m_c1, 2.0)) / 2.0 / CPair;
    m_Lrin = zeta_in * pow(m_c1, 2.0) / 2.0;
    double tmp = Kair / (Kair - 1.0) - m_Lrin / Rair / (m_t1 - m_ta);
    m_n1 = tmp / (tmp - 1.0);
    m_p1 = m_pa * pow(m_t1 / m_ta, m_n1 / (m_n1 - 1.0));
    m_rho1 = m_p1 * 1000000.0 / Rair / m_t1;
    m_A1 = m_Gk / itkr / m_c1a / m_rho1;
    m_D1 = sqrt((4.0 * m_A1) / (PI * (1.0 - pow(D0D1, 2.0))));
    m_D0 = D0D1 * m_D1;
    m_D2 = m_D1 / D1D2;
    m_ntk = 60.0 * m_u2 / PI / m_D2;
    m_Dm = sqrt(0.5 * (pow(m_D1, 2.0) + pow(m_D0, 2.0)));
    m_um = PI * m_Dm * m_ntk / 60.0;
    m_betam = atan(m_c1a / m_um) * 180.0 / PI;
    m_betamb = m_betam + ia;
    m_c1m = m_c1a / tau1;
    m_wm = sqrt(pow(m_um, 2.0) + pow(m_c1m, 2.0));
    m_Mw1 = sqrt((pow(PI * m_D1 * m_ntk / 60.0, 2.0) + pow(m_c1m, 2.0)) / (Kair * Rair * m_t1));
    m_Lr1 = zeta1 * pow(m_wm, 2.0) / 2.0;
    m_c2r = kc2r * m_c1a;
    m_Lr2 = zeta2 * pow(m_c2r, 2.0) / 2.0;
    m_Lrg = alphap * pow(m_u2, 2.0);
    m_mu = 1.0 / (1 + (2.0 / 3.0) * (PI / zk) / (1.0 - pow(m_Dm, 2.0) / pow(m_D2, 2.0)));
    m_t2 = m_t1 + (m_mu + alphap - 0.5 * pow(m_mu, 2.0)) * pow(m_u2, 2.0) / CPair;
    tmp = Kair / (Kair - 1.0) - (m_Lr1 + m_Lr2 + m_Lrg) / (Rair * (m_t2 - m_t1));
    m_n2 = tmp / (tmp - 1.0);
    m_p2 = m_p1 * pow(m_t2 / m_t1, m_n2 / (m_n2 - 1.0));
    m_rho2 = m_p2 * 1000000.0 / Rair / m_t2;
    m_c2u = m_mu * m_u2;
    m_c2 = sqrt(pow(m_c2u, 2.0) + pow(m_c2r, 2.0));
    m_w2u = m_u2 - m_c2u;
    m_w2r = m_c2r;
    m_w2 = sqrt(pow(m_w2u, 2.0) + pow(m_w2r, 2.0));
    m_beta2 = atan(m_w2r / m_w2u) * 180.0 / PI;
    m_alpha2 = atan(m_c2r / m_c2u) * 180.0 / PI;
    m_b2 = m_Gk / (PI * m_D2 * m_c2r * tau2 * m_rho2) / itkr;
    m_t2z = m_t2 + pow(m_c2, 2.0) / 2.0 / CPair;
    m_b3 = kb3 * m_b2;
    m_D3 = kD3 * m_D2;
    tmp = Kair / (Kair - 1.0) * zetad;
    m_n3 = tmp / (tmp - 1.0);
    double sigma = ((Kair - 1.0) / 2.0) * pow(m_c2, 2.0) / (Kair * Rair * m_t2);
    double m = 2.0 / (m_n3 - 1.0);
    double q = pow(m_D2 * m_b2 / m_D3 / m_b3, 2.0);
    double beta = T3T2_MINVAL;
    double x = 0;
    while (1) {
        if (beta > T3T2_MAXVAL) {
            cout << ERRORMSGBLANK
                 << "No solution for T3/T2 ratio.\n";
            break;
        }
        x = 1.0 / pow(beta, m) + (beta - 1.0) / (sigma * q) - 1.0 / q;
        if (fabs(x) <= T3T2_ACCUR) {
            break;
        }
        beta += T3T2_STEP;
    }
    m_t3 = beta * m_t2;
    m_p3 = m_p2 * pow(beta, m_n3 / (m_n3 - 1.0));
    m_rho3 = m_p3 * 1000000.0 / Rair / m_t3;
    m_c3 = m_c2 * (m_D2 * m_b2 * m_rho2) / (m_D3 * m_b3 * m_rho3);
    m_D4 = kD4 * m_D2;
    m_b4 = kb4 * m_b3;
    m_alpha4 = m_alpha2 + kalpha4;
    sigma = 0.5 * (Kair - 1.0) * pow(m_c3, 2.0) / (Kair * Rair * m_t3);
    m = 2.0 / (n4 - 1.0);
    q = pow((m_D3 * m_b3 * tau3 * sin(m_alpha2 * PI / 180.0)) /
            (m_D4 * m_b4 * tau4 * sin(m_alpha4 * PI / 180.0)), 2.0);
    beta = T4T3_MINVAL;
    x = 0;
    while (1) {
        if (beta > T4T3_MAXVAL) {
            cout << ERRORMSGBLANK
                 << "No solution for T4/T3 ratio.\n";
            break;
        }
        x = 1.0 / pow(beta, m) + (beta - 1.0) / (sigma * q) - 1.0 / q;
        if (fabs(x) <= T4T3_ACCUR) {
            break;
        }
        beta += T4T3_STEP;
    }
    m_t4 = beta * m_t3;
    m_p4 = m_p3 * pow(beta, n4 / (n4 - 1.0));
    m_rho4 = m_p4 * 1000000.0 / Rair / m_t4;
    m_c4 = m_c3 * (m_D3 * m_b3 * m_rho3 * tau3 * sin(m_alpha2 * PI / 180.0)) /
        (m_D4 * m_b4 * m_rho4 * tau4 * sin(m_alpha4 * PI / 180.0));
    m_ck = m_c4 / kck;
    m_tk = m_t4 + (pow(m_c4, 2.0) - pow(m_ck, 2.0)) / 2.0 / CPair;
    m_pk = m_p4 * pow(m_tk / m_t4, nul / (nul - 1.0));
    m_pikd = m_pk / m_pa;
    m_Lad = CPair * m_ta * (pow(m_pikd, (Kair - 1.0) / Kair) - 1.0);
    m_etaka = (pow(m_pikd, (Kair - 1.0) / Kair) - 1.0) / (m_tk / m_ta - 1.0);
    m_tkz = m_tk + pow(m_ck, 2.0) / 2.0 / CPair;
    m_etakn = m_Lad / pow(m_u2, 2.0);
    m_Nk = (m_Gk * m_Lad) / (1000.0 * m_etaka * itkr);

    //

    if ((m_n1 < 1.32) || (m_n1 > 1.39)) {
        cout << WARNMSGBLANK << "Strange value of n1.\n";
    }

    if (m_Mw1 > 0.9) {
        cout << WARNMSGBLANK << "Too large value of Mach number Mw1.\n";
    }

    if ((m_n2 < 1.45) || (m_n2 > 1.6)) {
        cout << WARNMSGBLANK << "Strange value of n2.\n";
    }

    const double b2D2 = m_b2 / m_D2;

    if ((b2D2 < 0.04) || (b2D2 > 0.08)) {
        cout << WARNMSGBLANK << "Strange value of b2/D2 ratio.\n";
    }

    if (fabs(m_etakn - etakn) > ETAKN_ACCUR) {
        cout << WARNMSGBLANK
             << "Low convergence of the parameter etakn. "
             << "Please change value of etakn and restart caclulation.\n";
    }

    // centripetal turbine calculation

    m_Gt = m_Gk * (1.0 + 1.0 / (alpha * phip * L0)) * sigmay;
    m_D1t = kD1t * m_D2;
    m_u1t = PI * m_D1t * m_ntk / 60.0;
    m_Ht = m_Lad / m_etaka / etat * m_Gk / m_Gt;
    m_ca = sqrt(2.0 * m_Ht);
    m_chi = m_u1t / m_ca;
    m_pt_z = pr_z * pow(1.0 - m_Ht / CPexh / tr_z, -Kexh / (Kexh - 1.0));
    m_D2t = kD2t * m_D1t;
    m_Dbt = kDbt * m_D1t;
    m_Dmt = sqrt(0.5 * (pow(m_D2t, 2.0) + pow(m_Dbt, 2.0)));
    m_mut = m_Dmt / m_D1t;
    m_rhot_min = pow(m_chi, 2.0) * (1.0 - pow(m_mut, 2.0));
    m_Hc = (1.0 - rhot) * m_Ht;
    m_c1t = phit * sqrt(2.0 * m_Hc);
    m_c1rt = m_c1t * sin(alpha1t * PI / 180.0);
    m_c1ut = m_c1t * cos(alpha1t * PI / 180.0);
    m_beta1t = 90.0 + atan((m_u1t - m_c1ut) / m_c1rt) * 180.0 / PI;
    m_t1t = tr_z - pow(m_c1t, 2.0) / 2.0 / CPexh;
    m_p1t = m_pt_z * pow(1.0 - m_Hc / (CPexh * tr_z), Kexh / (Kexh - 1.0));
    m_rho1t = (m_p1t * 1000000.0) / (Rexh * m_t1t);
    m_b1t = m_Gt / (PI * m_D1t * m_rho1t * m_c1rt) / itkr;
    m_w1t = sqrt(pow(m_c1t, 2.0) + pow(m_u1t, 2.0) - 2.0 * m_c1t * m_u1t * cos(alpha1t * PI / 180.0));
    m_Hl = rhot * m_Ht;
    m_w2t = psi * sqrt(pow(m_w1t, 2.0) + 2.0 * m_Hl - pow(m_u1t, 2.0) * (1.0 - pow(m_mut, 1.0)));
    m_t2t = m_t1t - (pow(m_w2t, 2.0) - pow(m_w1t, 2.0) + pow(m_u1t, 2.0) * (1.0 - pow(m_mut, 2.0))) / (2.0 * CPexh);
    m_t2at = tr_z * pow(pr_z / m_pt_z, (Kexh - 1.0) / Kexh);
    m_rho2t = (pr_z * 1000000.0) / (Rexh * m_t2t);
    m_A2t = PI / 4.0 * (pow(m_D2t, 2.0) - pow(m_Dbt, 2.0));
    m_beta2st = asin(m_Gt / (m_w2t * m_A2t * m_rho2t) / 2.0) * 180.0 / PI;
    m_Gl = 0.45 * 2.0 * delta * (1.0 + (m_D2t - m_Dbt) / (2.0 * m_Dmt)) / (m_D2t - m_Dbt);
    m_beta2t = asin((m_Gt - m_Gl) / m_Gt * sin(m_beta2st * PI / 180.0)) * 180.0 / PI;
    m_umt = m_mut * m_u1t;
    m_c2ut = m_umt - m_w2t * cos(m_beta2t * PI / 180.0);
    m_c2rt = m_w2t * sin(m_beta2t * PI / 180.0);
    m_c2t = sqrt(pow(m_c2ut, 2.0) + pow(m_c2rt, 2.0));
    m_alpha2t = 90.0 + asin(m_c2ut / m_c2rt) * 180.0 / PI;
    m_Lu = m_u1t * m_c1ut - m_umt * m_c2ut;
    m_etau = m_Lu / m_Ht;
    m_zc = (1.0 / pow(phit, 2.0) - 1.0) * pow(m_c1t, 2.0) / 2.0;
    m_zl = (1.0 / pow(psi, 2.0) - 1.0) * pow(m_w2t, 2.0) / 2.0;
    m_etaat = 1.0 - (m_zc + m_zl) / m_Ht;
    m_zout = pow(m_c2t, 2.0) / 2.0;
    m_zleak = m_Lu * m_Gl / m_Gt;
    m_Nrd = kNrd * pow(m_D1t, 2.0) * pow(m_u1t / 100.0, 3.0) * ((m_rho1t + m_rho2t) / 2.0) * 0.735;
    m_zrd = m_Nrd * 1000.0 * itkr / m_Gt;
    m_etait = (m_Lu - m_zrd - m_zleak) / m_Ht;
    m_etat = m_etait * etamt;
    m_Nt = m_Ht * m_Gt * m_etat / itkr / 1000.0;

    //

    if ((m_chi < 0.64) || (m_chi > 0.7)) {
        cout << WARNMSGBLANK << "Strange value of chi ratio.\n";
    }

    if ((m_mut < 0.5) || (m_mut > 0.8)) {
        cout << WARNMSGBLANK << "Strange value of mut.\n";
    }

    if ((m_beta1t < 80.0) || (m_beta1t > 100.0)) {
        cout << WARNMSGBLANK << "Strange value of beta1t.\n";
    }

    if (fabs(m_etat - etat) > ETAT_ACCUR) {
        cout << WARNMSGBLANK
             << "Low convergence of the parameter etat. "
             << "Please change value of etat and restart caclulation.\n";
    }

    double N_diff = 0;

    if (m_Nk > m_Nt) {
        N_diff = (m_Nk - m_Nt) / m_Nk * 100.0;
    }
    else {
        N_diff = (m_Nt - m_Nk) / m_Nt * 100.0;
    }

    if (N_diff > 5.0) {
        cout << WARNMSGBLANK
             << "Low convergence between parameters Nk and Nt (< 5%). "
             << "Please correct source data and restart caclulation.\n";
    }

    //

    return true;
}

bool Calc::createReport() const {

    string programName(PRGNAME);
    string dateTime = currDateTime();
    string srcDataFilename = programName + "_source-data_" + dateTime + ".txt";
    string reportFilename = programName + "_results_" + dateTime + ".txt";

    //

    ifstream fin_srcdata(CONFIGFILE);
    string srcdata;

    if (!fin_srcdata) {
        cout << WARNMSGBLANK << "Can not open file \""
             << CONFIGFILE << "\" to read!\n";
        cout << WARNMSGBLANK << "Copying source data will be skipping.\n\n";
    }
    else {
        srcdata.assign((istreambuf_iterator<char>(fin_srcdata)),
                       (istreambuf_iterator<char>()));
    }

    fin_srcdata.close();

    ofstream fout_srccopy(srcDataFilename);

    if (!fout_srccopy) {
        cout << WARNMSGBLANK << "Can not open file \""
             << srcDataFilename << "\" to write!\n";
        cout << WARNMSGBLANK << "Copying source data will be skipping.\n\n";
    }
    else {
        fout_srccopy << srcdata;
    }

    fout_srccopy.close();

    cout << MSGBLANK << "Source data copied to \""
         << srcDataFilename << "\".\n";

    //

    ofstream fout(reportFilename);

    if (!fout) {
        cout << ERRORMSGBLANK << "Can not open file \""
             << reportFilename << "\" to write!\n";
        return false;
    }

    fout << PRGNAME << "\nv" << PRGVERSION << "\n\n";
    fout << "Calculation source data: " << srcDataFilename << ".\n\n";

    fout << "Results of radial flow compressor calculation\n\n";
    fout << "Gk     = " << m_Gk      << "\n";
    fout << "ta_z   = " << m_ta_z    << "\n";
    fout << "pa_z   = " << m_pa_z    << "\n";
    fout << "ta     = " << m_ta      << "\n";
    fout << "pa     = " << m_pa      << "\n";
    fout << "La     = " << m_La      << "\n";
    fout << "u2     = " << m_u2      << "\n";
    fout << "c1     = " << m_c1      << "\n";
    fout << "c1a    = " << m_c1a     << "\n";
    fout << "t1     = " << m_t1      << "\n";
    fout << "Lrin   = " << m_Lrin    << "\n";
    fout << "n1     = " << m_n1      << "\n";
    fout << "p1     = " << m_p1      << "\n";
    fout << "rho1   = " << m_rho1    << "\n";
    fout << "A1     = " << m_A1      << "\n";
    fout << "D1     = " << m_D1      << "\n";
    fout << "D0     = " << m_D0      << "\n";
    fout << "D2     = " << m_D2      << "\n";
    fout << "ntk    = " << m_ntk     << "\n";
    fout << "Dm     = " << m_Dm      << "\n";
    fout << "um     = " << m_um      << "\n";
    fout << "betam  = " << m_betam   << "\n";
    fout << "betamb = " << m_betamb  << "\n";
    fout << "c1m    = " << m_c1m     << "\n";
    fout << "wm     = " << m_wm      << "\n";
    fout << "Mw1    = " << m_Mw1     << "\n";
    fout << "Lr1    = " << m_Lr1     << "\n";
    fout << "c2r    = " << m_c2r     << "\n";
    fout << "Lr2    = " << m_Lr2     << "\n";
    fout << "Lrg    = " << m_Lrg     << "\n";
    fout << "mu     = " << m_mu      << "\n";
    fout << "t2     = " << m_t2      << "\n";
    fout << "n2     = " << m_n2      << "\n";
    fout << "p2     = " << m_p2      << "\n";
    fout << "rho2   = " << m_rho2    << "\n";
    fout << "c2u    = " << m_c2u     << "\n";
    fout << "c2     = " << m_c2      << "\n";
    fout << "w2u    = " << m_w2u     << "\n";
    fout << "w2r    = " << m_w2r     << "\n";
    fout << "w2     = " << m_w2      << "\n";
    fout << "beta2  = " << m_beta2   << "\n";
    fout << "alpha2 = " << m_alpha2  << "\n";
    fout << "b2     = " << m_b2      << "\n";
    fout << "t2z    = " << m_t2z     << "\n";
    fout << "b3     = " << m_b3      << "\n";
    fout << "D3     = " << m_D3      << "\n";
    fout << "n3     = " << m_n3      << "\n";
    fout << "t3     = " << m_t3      << "\n";
    fout << "p3     = " << m_p3      << "\n";
    fout << "rho3   = " << m_rho3    << "\n";
    fout << "c3     = " << m_c3      << "\n";
    fout << "D4     = " << m_D4      << "\n";
    fout << "b4     = " << m_b4      << "\n";
    fout << "alpha4 = " << m_alpha4  << "\n";
    fout << "t4     = " << m_t4      << "\n";
    fout << "p4     = " << m_p4      << "\n";
    fout << "rho4   = " << m_rho4    << "\n";
    fout << "c4     = " << m_c4      << "\n";
    fout << "ck     = " << m_ck      << "\n";
    fout << "tk     = " << m_tk      << "\n";
    fout << "pk     = " << m_pk      << "\n";
    fout << "pikd   = " << m_pikd    << "\n";
    fout << "Lad    = " << m_Lad     << "\n";
    fout << "etaka  = " << m_etaka   << "\n";
    fout << "tkz    = " << m_tkz     << "\n";
    fout << "etakn  = " << m_etakn   << "\n";
    fout << "Nk     = " << m_Nk      << "\n\n";

    fout << "Results of centripetal turbine calculation\n\n";
    fout << "m_Gt       = " << m_Gt       << "\n";
    fout << "m_D1t      = " << m_D1t      << "\n";
    fout << "m_u1t      = " << m_u1t      << "\n";
    fout << "m_Ht       = " << m_Ht       << "\n";
    fout << "m_ca       = " << m_ca       << "\n";
    fout << "m_chi      = " << m_chi      << "\n";
    fout << "m_pt_z     = " << m_pt_z     << "\n";
    fout << "m_D2t      = " << m_D2t      << "\n";
    fout << "m_Dbt      = " << m_Dbt      << "\n";
    fout << "m_Dmt      = " << m_Dmt      << "\n";
    fout << "m_mut      = " << m_mut      << "\n";
    fout << "m_rhot_min = " << m_rhot_min << "\n";
    fout << "m_Hc       = " << m_Hc       << "\n";
    fout << "m_c1t      = " << m_c1t      << "\n";
    fout << "m_c1rt     = " << m_c1rt     << "\n";
    fout << "m_c1ut     = " << m_c1ut     << "\n";
    fout << "m_beta1t   = " << m_beta1t   << "\n";
    fout << "m_t1t      = " << m_t1t      << "\n";
    fout << "m_p1t      = " << m_p1t      << "\n";
    fout << "m_rho1t    = " << m_rho1t    << "\n";
    fout << "m_b1t      = " << m_b1t      << "\n";
    fout << "m_w1t      = " << m_w1t      << "\n";
    fout << "m_Hl       = " << m_Hl       << "\n";
    fout << "m_w2t      = " << m_w2t      << "\n";
    fout << "m_t2t      = " << m_t2t      << "\n";
    fout << "m_t2at     = " << m_t2at     << "\n";
    fout << "m_rho2t    = " << m_rho2t    << "\n";
    fout << "m_A2t      = " << m_A2t      << "\n";
    fout << "m_beta2st  = " << m_beta2st  << "\n";
    fout << "m_Gl       = " << m_Gl       << "\n";
    fout << "m_beta2t   = " << m_beta2t   << "\n";
    fout << "m_umt      = " << m_umt      << "\n";
    fout << "m_c2ut     = " << m_c2ut     << "\n";
    fout << "m_c2rt     = " << m_c2rt     << "\n";
    fout << "m_c2t      = " << m_c2t      << "\n";
    fout << "m_alpha2t  = " << m_alpha2t  << "\n";
    fout << "m_Lu       = " << m_Lu       << "\n";
    fout << "m_etau     = " << m_etau     << "\n";
    fout << "m_zc       = " << m_zc       << "\n";
    fout << "m_zl       = " << m_zl       << "\n";
    fout << "m_etaat    = " << m_etaat    << "\n";
    fout << "m_zout     = " << m_zout     << "\n";
    fout << "m_zleak    = " << m_zleak    << "\n";
    fout << "m_Nrd      = " << m_Nrd      << "\n";
    fout << "m_zrd      = " << m_zrd      << "\n";
    fout << "m_etait    = " << m_etait    << "\n";
    fout << "m_etat     = " << m_etat     << "\n";
    fout << "m_Nt       = " << m_Nt       << "\n";

    fout.close();

    cout << MSGBLANK << "Report file \"" << reportFilename << "\"created.\n\n";

    return true;
}
