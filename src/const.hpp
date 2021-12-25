/*
  tkrsadi
  Turbocharger parameters calculation.

  File: const.hpp

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

#ifndef CONST_HPP
#define CONST_HPP

#define CONFIGFILE     "tkrsadi_conf.txt"
#define COMMENTREGEX   "^[ ]*//.*"
#define PARAMDELIMITER "="
#define ELEMDELIMITER  ","
#define CSVDELIMITER   ";"
#define MSGBLANK       "tkrsadi =>\t"
#define ERRORMSGBLANK  "tkrsadi ERROR =>\t"
#define WARNMSGBLANK   "tkrsadi WARNING =>\t"

#define PI 3.14159265
#define E 2.71828182

#define Kair 1.4
#define Rair 287.2
#define CPair 1005.0

#define Kexh 1.34
#define Rexh 286.4
#define CPexh 1128.7

#define L0 14.5

#define T3T2_MINVAL 1
#define T3T2_MAXVAL 1.5
#define T3T2_STEP 0.0001
#define T3T2_ACCUR 0.001

#define T4T3_MINVAL 1
#define T4T3_MAXVAL 1.5
#define T4T3_STEP 0.0001
#define T4T3_ACCUR 0.001

#define ETAKN_ACCUR 0.01
#define ETAT_ACCUR 0.01

#endif // CONST_HPP
