## Copyright (C) 2018 Windi
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} leastSquareWithHouseHolder (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Windi <Windi@WIN8>
## Created: 2018-03-13

function [x] = leastSquareWithHouseholder (A, b)
  [Q, R, p] = householder(A);
  [m, n] = size(A);
  x = zeros(n,1);
  
  %permutate b
  tmp = b;
  for i = 1:m
    b(i) = tmp(p(i));
  endfor
  
  %back subtitution
  Qb = Q' * b;
  Qb
  x(n) = Qb(n) / R(n,n);
  for i = n-1:-1:1
    x(i) = (Qb(i) - R(i, i+1:n) * x(i+1:n)) / R(i,i);
  endfor

endfunction
