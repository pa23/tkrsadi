/*
  tkrsadi
  Turbocharger parameters calculation.

  File: conf.cpp

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

#include "conf.hpp"
#include "const.hpp"
#include "prgid.hpp"
#include "auxf.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>

using std::cout;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::regex;
using std::regex_match;

bool Conf::readConfigFile() {

    ifstream fin(CONFIGFILE);

    if (!fin) {

        cout << WARNMSGBLANK << "Can not open file \""
             << CONFIGFILE << "\" to read!\n";
        cout << WARNMSGBLANK << "I will try to create blank.\n\n";

        if (createBlank()) {
            cout << MSGBLANK << "Blank created. "
                 << "Please edit it and restart calculation.\n\n";
        }
        else {
            cout << ERRORMSGBLANK << "Can not create blank!\n\n";
        }

        return false;
    }

    string s;
    vector<string> elem;

    while (!fin.eof()) {

        getline(fin, s);

        if (regex_match(s, regex(COMMENTREGEX)) || s.empty()) {
            continue;
        }

        splitString(s, elem, PARAMDELIMITER);

        if (elem.size() != 2) {

            s.clear();
            elem.clear();

            continue;
        }

        if      ( elem[0] == "Ne"       ) { m_Ne       = stringToDouble(elem[1]); }
        else if ( elem[0] == "ge"       ) { m_ge       = stringToDouble(elem[1]); }
        else if ( elem[0] == "alpha"    ) { m_alpha    = stringToDouble(elem[1]); }
        else if ( elem[0] == "phip"     ) { m_phip     = stringToDouble(elem[1]); }
        else if ( elem[0] == "pik"      ) { m_pik      = stringToDouble(elem[1]); }
        else if ( elem[0] == "itkr"     ) { m_itkr     = stringToDouble(elem[1]); }
        else if ( elem[0] == "p0"       ) { m_p0       = stringToDouble(elem[1]); }
        else if ( elem[0] == "t0"       ) { m_t0       = stringToDouble(elem[1]); }
        else if ( elem[0] == "sigma_in" ) { m_sigma_in = stringToDouble(elem[1]); }
        else if ( elem[0] == "ca"       ) { m_ca       = stringToDouble(elem[1]); }
        else if ( elem[0] == "etakn"    ) { m_etakn    = stringToDouble(elem[1]); }
        else if ( elem[0] == "phi"      ) { m_phi      = stringToDouble(elem[1]); }
        else if ( elem[0] == "zeta_in"  ) { m_zeta_in  = stringToDouble(elem[1]); }
        else if ( elem[0] == "D0D1"     ) { m_D0D1     = stringToDouble(elem[1]); }
        else if ( elem[0] == "D1D2"     ) { m_D1D2     = stringToDouble(elem[1]); }
        else if ( elem[0] == "ia"       ) { m_ia       = stringToDouble(elem[1]); }
        else if ( elem[0] == "tau1"     ) { m_tau1     = stringToDouble(elem[1]); }
        else if ( elem[0] == "zeta1"    ) { m_zeta1    = stringToDouble(elem[1]); }
        else if ( elem[0] == "zeta2"    ) { m_zeta2    = stringToDouble(elem[1]); }
        else if ( elem[0] == "kc2r"     ) { m_kc2r     = stringToDouble(elem[1]); }
        else if ( elem[0] == "alphap"   ) { m_alphap   = stringToDouble(elem[1]); }
        else if ( elem[0] == "zk"       ) { m_zk       = stringToDouble(elem[1]); }
        else if ( elem[0] == "tau2"     ) { m_tau2     = stringToDouble(elem[1]); }
        else if ( elem[0] == "kb3"      ) { m_kb3      = stringToDouble(elem[1]); }
        else if ( elem[0] == "kD3"      ) { m_kD3      = stringToDouble(elem[1]); }
        else if ( elem[0] == "zetad"    ) { m_zetad    = stringToDouble(elem[1]); }
        else if ( elem[0] == "kD4"      ) { m_kD4      = stringToDouble(elem[1]); }
        else if ( elem[0] == "kb4"      ) { m_kb4      = stringToDouble(elem[1]); }
        else if ( elem[0] == "kalpha4"  ) { m_kalpha4  = stringToDouble(elem[1]); }
        else if ( elem[0] == "zd"       ) { m_zd       = stringToDouble(elem[1]); }
        else if ( elem[0] == "n4"       ) { m_n4       = stringToDouble(elem[1]); }
        else if ( elem[0] == "tau3"     ) { m_tau3     = stringToDouble(elem[1]); }
        else if ( elem[0] == "tau4"     ) { m_tau4     = stringToDouble(elem[1]); }
        else if ( elem[0] == "nul"      ) { m_nul      = stringToDouble(elem[1]); }
        else if ( elem[0] == "kck"      ) { m_kck      = stringToDouble(elem[1]); }

        else if ( elem[0] == "pr"       ) { m_pr_z     = stringToDouble(elem[1]); }
        else if ( elem[0] == "tr"       ) { m_tr_z     = stringToDouble(elem[1]); }
        else if ( elem[0] == "sigmay"   ) { m_sigmay   = stringToDouble(elem[1]); }
        else if ( elem[0] == "kD1t"     ) { m_kD1t     = stringToDouble(elem[1]); }
        else if ( elem[0] == "etat"     ) { m_etat     = stringToDouble(elem[1]); }
        else if ( elem[0] == "kD2t"     ) { m_kD2t     = stringToDouble(elem[1]); }
        else if ( elem[0] == "kDbt"     ) { m_kDbt     = stringToDouble(elem[1]); }
        else if ( elem[0] == "rhot"     ) { m_rhot     = stringToDouble(elem[1]); }
        else if ( elem[0] == "phit"     ) { m_phit     = stringToDouble(elem[1]); }
        else if ( elem[0] == "alpha1t"  ) { m_alpha1t  = stringToDouble(elem[1]); }
        else if ( elem[0] == "psi"      ) { m_psi      = stringToDouble(elem[1]); }
        else if ( elem[0] == "delta"    ) { m_delta    = stringToDouble(elem[1]); }
        else if ( elem[0] == "kNrd"     ) { m_kNrd     = stringToDouble(elem[1]); }
        else if ( elem[0] == "etamt"    ) { m_etamt    = stringToDouble(elem[1]); }

        s.clear();
        elem.clear();
    }

    fin.close();

    return true;
}

bool Conf::createBlank() const {

    ofstream fout(CONFIGFILE);

    if (!fout) {
        cout << ERRORMSGBLANK << "Can not open file \""
             << CONFIGFILE << "\" to write!\n";
        return false;
    }

    fout << "//\n"
         << "// This is " << PRGNAME << " configuration file.\n"
         << "// Parameter-Value delimiter is symbol \"" << PARAMDELIMITER << "\".\n"
         << "// Vector elements delimiter is symbol \"" << ELEMDELIMITER << "\".\n"
         << "// Table data is entered line by line.\n"
         << "// Text after \"//\" is comment.\n" << "//\n\n";

    fout << "// Мощность двигателя, кВт\n"
         << "Ne=185\n\n"
         << "// Удельный эффективный расход топлива, г/кВтч\n"
         << "ge=230\n\n"
         << "// Коэффициент избытка воздуха\n"
         << "alpha=1.8\n\n"
         << "// Коэффициент продувки\n"
         << "phip=1.15\n\n"
         << "// Степень повышения давления в компрессоре\n"
         << "pik=1.4143\n\n"
         << "// Количество турбокомпрессоров\n"
         << "itkr=2\n\n"
         << "// Барометрическое давление, кПа\n"
         << "p0=101.3\n\n"
         << "// Температура окружающего воздуха, грЦ\n"
         << "t0=20\n\n"
         << "// Коэффициент потерь давления в системе очистки воздуха\n"
         << "sigma_in=0.98\n\n"
         << "// Скорость воздуха на входе в компрессор (30..80), м/с\n"
         << "ca=50\n\n"
         << "// Напорный адиабатный КПД (0.58..0.73)\n"
         << "etakn=0.7\n\n"
         << "// Коэффициент расхода компрессора (0.25..0.35)\n"
         << "phi=0.3\n\n"
         << "// Коэффициент потерь энергии во входном устройстве компрессора (0.03..0.06)\n"
         << "zeta_in=0.035\n\n"
         << "// Отношение диаметра втулки колеса компрессора к диаметру колеса на входе (0.28..0.55)\n"
         << "D0D1=0.4\n\n"
         << "// Отношение диаметра колеса компрессора на входе к диаметру колеса (0.55..0.7)\n"
         << "D1D2=0.6\n\n"
         << "// Угол атаки лопаток колеса компрессора на среднем диаметре (2..4), градус\n"
         << "ia=2.7\n\n"
         << "// Коэффициент загромождения на входе в колесо компрессора, учитывающий толщину лопаток (0.8..0.9)\n"
         << "tau1=0.8\n\n"
         << "// Коэффициент потерь на ВНА (0.1..0.3)\n"
         << "zeta1=0.15\n\n"
         << "// Коэффициент потерь на трение по поворот потока в межлопаточных каналах колеса компрессора (0.1..0.2)\n"
         << "zeta2=0.15\n\n"
         << "// Коэффициент радиальной составляющей абс. скорости на выходе из колеса компрессора (1..1.1)\n"
         << "kc2r=1\n\n"
         << "// Коэффициент дисковых потерь колеса компрессора (0.04..0.08)\n"
         << "alphap=0.06\n\n"
         << "// Количество лопаток колеса компрессора (9..34)\n"
         << "zk=15\n\n"
         << "// Коэффициент загромождения на выходе из колеса компрессора (0.92..0.96)\n"
         << "tau2=0.93\n\n"
         << "// Коэффициент ширины безлопаточной части диффузора на выходе (0.95..1)\n"
         << "kb3=1\n\n"
         << "// Коэффициент наружного диаметра безлопаточной части диффузора (1.05..1.2)\n"
         << "kD3=1.1\n\n"
         << "// Политропный КПД диффузора компрессора (0.5..0.75)\n"
         << "zetad=0.7\n\n"
         << "// Коэффициент наружного диаметра лопаточной части диффузора (1.35..1.7)\n"
         << "kD4=1.4\n\n"
         << "// Коэффициент ширины лопаточной части диффузора (1..1.3)\n"
         << "kb4=1\n\n"
         << "// Коэффициент угла наклона лопаток диффузора (12..18), гр\n"
         << "kalpha4=16\n\n"
         << "// Количество лопаток диффузора (9..36)\n"
         << "zd=17\n\n"
         << "// Показатель политропы сжатия в лопаточном диффузоре (1.5..1.8)\n"
         << "n4=1.6\n\n"
         << "// Коэффициент загромождения на входе в лопаточный диффузор (0.95..0.97)\n"
         << "tau3=0.96\n\n"
         << "// Коэффициент загромождения на выходе из лопаточного диффузора (0.95..0.97)\n"
         << "tau4=0.96\n\n"
         << "// Показатель политропы сжатия улитки (1.9..2.2)\n"
         << "nul=2.2\n\n"
         << "// Коэффициент скорости на выходе из компрессора (1.3..1.4)\n"
         << "kck=1.3\n\n"
         << "// Давление газов за турбиной, кПа\n"
         << "pr=106\n\n"
         << "// Температура газов за турбиной, грЦ\n"
         << "tr=650\n\n"
         << "// Коэффициент, учитывающий утечку газа (0.98..0.99)\n"
         << "sigmay=0.99\n\n"
         << "// Коэффициент диаметра колеса турбины\n"
         << "kD1t=1\n\n"
         << "// Эффективный КПД турбины\n"
         << "etat=0.8\n\n"
         << "// Коэффициент наружного диаметра колеса турбины (0.7..0.85)\n"
         << "kD2t=0.85\n\n"
         << "// Коэффициент диаметра втулки колеса турбины на выходе (0.25..0.32)\n"
         << "kDbt=0.32\n\n"
         << "// Степень реактивности турбины (0.45..0.55)\n"
         << "rhot=0.5\n\n"
         << "// Коэффициент скорости, учитывающий потери в сопловом аппарате турбины (0.93..0.96)\n"
         << "phit=0.96\n\n"
         << "// Угол выхода из соплового аппарата турбины (15..25), градус\n"
         << "alpha1t=22\n\n"
         << "// Коэффициент скорости учитывающий потери в рабочем колесе турбины (0.85..0.94)\n"
         << "psi=0.9\n\n"
         << "// Зазор между корпусом и колесом турбины (0.5..1.5), мм\n"
         << "delta=0.5\n\n"
         << "// Коэффициент, зависящий от типа рабочего колеса турбины (3.5..5)\n"
         << "kNrd=4\n\n"
         << "// Механический КПД турбины (0.9..0,96)\n"
         << "etamt=0.948\n\n";

    fout.close();

    return true;
}
