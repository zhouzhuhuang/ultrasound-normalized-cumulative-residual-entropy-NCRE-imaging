% This callback function is for computation of a local NCRE estimate
%
% Input:
%    block_struct: the Matlab-defined structure for the blockproc() function
% Output:
%     cre: the local NCRE estimate
%
% Please see the following paper for the details of the algorithm:
% [1] Ruiyang Gao, Po-Hsiang Tsui, Sinan Li, Guangyu Bin, Dar-In Tai, Shuicai Wu, Zhuhuang Zhou,
% Ultrasound normalized cumulative residual entropy imaging: Theory, methodology, and application,
% Computer Methods and Programs in Biomedicine,
% 2024,
% 108374,
% ISSN 0169-2607,
% https://doi.org/10.1016/j.cmpb.2024.108374.
% (https://www.sciencedirect.com/science/article/pii/S0169260724003675)
%
%
% ABOUT:
%     author               - Zhuhuang Zhou
%     date                 - 14th August 2024
%     last update          - 14th August 2024
%
% Copyright (C) 2024 Zhuhuang Zhou
% zhouzhuhuang@126.com
%
%
% This code is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this code. If not, see <http://www.gnu.org/licenses/>.

function cre = blockCRE(block_struct)

blocks = block_struct.data;
blocks = blocks(:);

% data normalization
x1 = blocks/rms(blocks);

% CDF estimation
[yCDF1, xCDF1] = ecdf(x1);
yCDF1 = yCDF1(2:end);
xCDF1 = xCDF1(2:end);

xCDF1(isnan(yCDF1)) = [];
yCDF1(isnan(yCDF1)) = [];

% normalized cumulative residual entropy estimation
cre = -trapz(xCDF1, (1-yCDF1).*log((1-yCDF1)+eps));

if isnan(cre)
  cre = 0;
end



