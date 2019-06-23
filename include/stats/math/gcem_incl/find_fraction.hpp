/*################################################################################
  ##
  ##   Copyright (C) 2016-2018 Keith O'Hara
  ##
  ##   This file is part of the GCE-Math C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * find the fraction part of x = n + r, where -0.5 < r < 0.5
 */

#ifndef _gcem_find_fraction_HPP
#define _gcem_find_fraction_HPP

template<typename T>
constexpr
T
find_fraction(const T x)
{
    return( abs(x - T(int(x))) > T(0.5) ? \
            // if 
                x - T(int(x)) - sgn(x) : 
            //else 
                x - T(int(x)) );
}

#endif
