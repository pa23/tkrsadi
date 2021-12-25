/*
  tkrsadi
  Turbocharger parameters calculation.

  File: calc.hpp

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

#ifndef CALC_HPP
#define CALC_HPP

#include <vector>
#include <memory>

#include "conf.hpp"

class Calc {

public:

    Calc(const std::shared_ptr<Conf> &conf);

    bool calculate();
    bool createReport() const;

private:

    std::shared_ptr<Conf> m_conf;

    double m_Gk     = 0;
    double m_ta_z   = 0;
    double m_pa_z   = 0;
    double m_ta     = 0;
    double m_pa     = 0;
    double m_La     = 0;
    double m_u2     = 0;
    double m_c1     = 0;
    double m_c1a    = 0;
    double m_t1     = 0;
    double m_Lrin   = 0;
    double m_n1     = 0;
    double m_p1     = 0;
    double m_rho1   = 0;
    double m_A1     = 0;
    double m_D1     = 0;
    double m_D0     = 0;
    double m_D2     = 0;
    double m_ntk    = 0;
    double m_Dm     = 0;
    double m_um     = 0;
    double m_betam  = 0;
    double m_betamb = 0;
    double m_c1m    = 0;
    double m_wm     = 0;
    double m_Mw1    = 0;
    double m_Lr1    = 0;
    double m_c2r    = 0;
    double m_Lr2    = 0;
    double m_Lrg    = 0;
    double m_mu     = 0;
    double m_t2     = 0;
    double m_n2     = 0;
    double m_p2     = 0;
    double m_rho2   = 0;
    double m_c2u    = 0;
    double m_c2     = 0;
    double m_w2u    = 0;
    double m_w2r    = 0;
    double m_w2     = 0;
    double m_beta2  = 0;
    double m_alpha2 = 0;
    double m_b2     = 0;
    double m_t2z    = 0;
    double m_b3     = 0;
    double m_D3     = 0;
    double m_n3     = 0;
    double m_t3     = 0;
    double m_p3     = 0;
    double m_rho3   = 0;
    double m_c3     = 0;
    double m_D4     = 0;
    double m_b4     = 0;
    double m_alpha4 = 0;
    double m_t4     = 0;
    double m_p4     = 0;
    double m_rho4   = 0;
    double m_c4     = 0;
    double m_ck     = 0;
    double m_tk     = 0;
    double m_pk     = 0;
    double m_pikd   = 0;
    double m_Lad    = 0;
    double m_etaka  = 0;
    double m_tkz    = 0;
    double m_etakn  = 0;
    double m_Nk     = 0;

    double m_Gt       = 0;
    double m_D1t      = 0;
    double m_u1t      = 0;
    double m_Ht       = 0;
    double m_ca       = 0;
    double m_chi      = 0;
    double m_pt_z     = 0;
    double m_D2t      = 0;
    double m_Dbt      = 0;
    double m_Dmt      = 0;
    double m_mut      = 0;
    double m_rhot_min = 0;
    double m_Hc       = 0;
    double m_c1t      = 0;
    double m_c1rt     = 0;
    double m_c1ut     = 0;
    double m_beta1t   = 0;
    double m_t1t      = 0;
    double m_p1t      = 0;
    double m_rho1t    = 0;
    double m_b1t      = 0;
    double m_w1t      = 0;
    double m_Hl       = 0;
    double m_w2t      = 0;
    double m_t2t      = 0;
    double m_t2at     = 0;
    double m_rho2t    = 0;
    double m_A2t      = 0;
    double m_beta2st  = 0;
    double m_Gl       = 0;
    double m_beta2t   = 0;
    double m_umt      = 0;
    double m_c2ut     = 0;
    double m_c2rt     = 0;
    double m_c2t      = 0;
    double m_alpha2t  = 0;
    double m_Lu       = 0;
    double m_etau     = 0;
    double m_zc       = 0;
    double m_zl       = 0;
    double m_etaat    = 0;
    double m_zout     = 0;
    double m_zleak    = 0;
    double m_Nrd      = 0;
    double m_zrd      = 0;
    double m_etait    = 0;
    double m_etat     = 0;
    double m_Nt       = 0;

};

#endif // CALC_HPP
