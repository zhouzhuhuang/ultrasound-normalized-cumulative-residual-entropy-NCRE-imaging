% This function is for parallel computation of local NCRE estimates
%
% Input:
%    fs: sampling frequency [Hz]
%    v: sound speed [m/s]
%    scan_step: distance between two adjacent scan lines [m]
%    pulselength: pulse length of the transducer [m]
%    cc: for the side length of the block, cc time(s) of the pulselength
%    envelope: the envelopes of rf signals
% Output:
%     ncre: NCRE map
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

function ncre = ent_map_block_cre(fs, v, scan_step, pulselength, cc, envelope)

%obtain image size
image_range = size(envelope);
line_number = image_range(2);
datalength = image_range(1);

% determine the block size
y = 0.5 * v / fs;
inte = round(cc * pulselength / y);
inte1 = round(cc * pulselength / scan_step);
bs1 = inte;
bs2 = inte1;

% overlap ratio between adjacent blocks
overlap1 = 0.9;
overlap2 = 0.9;

% parameters for the blockproc() function
% refer to the Matlab doc for the blockproc() function
v = floor(inte*overlap1/2);
h = floor(inte1*overlap2/2);
m = round((1-overlap1)*bs1);
n = round((1-overlap2)*bs2);

% parallel computation of local NCRE estimates,
% where blockCRE is the callback function
ncre = blockproc(envelope, [m,n], @(block_struct) blockCRE(block_struct), 'BorderSize', [v,h], 'UseParallel', true, 'TrimBorder', false);

