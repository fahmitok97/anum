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
## @deftypefn {} {@var{retval} =} householder (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Windi <Windi@WIN8>
## Created: 2018-03-13

function [Q, R, p] = householder (A)
  %input matrix A
  %output matrix Q & R and pivot vector p
  
  [m, n] = size(A);
  Q = eye(m);
  p = 1:m;
  
  %loop columns
  for i = 1:n
    %pivoting
    
    %find max element
    [maxElem, maxElemRow] = max(abs(A(i:m,i)));
    %swap A
    tmp = A(i,:); A(i,:) = A(maxElemRow,:); A(maxElemRow,:) = tmp;
    %swap Q
    tmp = Q(i,:); Q(i,:) = Q(maxElemRow,:); Q(maxElemRow,:) = tmp;
    %swap p
    tmp = p(i); p(i) = p(maxElemRow); p(maxElemRow) = tmp;
    
    %init v
    v = A(i:m, i);
    v(1) = v(1) + sign(v(1))*norm(v);
    
    tmp = 2 / (v' * v); %avoid compute this constant twice
    A(i:m, i:n) = A(i:m, i:n) - (tmp * v) * (v' * A(i:m, i:n));
    
    Q(:, i:m) = Q(:, i:m) - (Q(:, i:m) * v) * (tmp * v');
  endfor
  R = A;
endfunction
