/*
  tkrsadi
  Turbocharger parameters calculation.

  File: auxf.hpp

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

#ifndef AUXF_HPP
#define AUXF_HPP

#include <string>
#include <vector>

std::string uintToString(size_t);
double stringToDouble(const std::string &);
bool stringToBool(const std::string &);
void splitString(
    const std::string &,        // source string
    std::vector<std::string> &, // vector for elements
    const std::string &         // delimeter
    );

std::string currDateTime();
std::string trimDate(const std::string &);

double maxValue(const std::vector<double> &);

#endif // AUXF_HPP
