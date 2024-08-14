% This code demonstrates ultrasound normalized cumulative residual entropy (NCRE) imaging
%
%
% Please see the following paper for the algorithm of the ultrasound NCRE
% imaging
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

clear all;
close all;

% constant parameters
v = 1540;   % sound speed [m/s]
MHz = 10^6; % [Hz]
center_f = 5 * MHz; % central frequency [Hz]
fs = 10 * center_f; % sampling frequency [Hz]
pulselength = 0.89 * 10^(-3) % pulse length of the transducer [m]

% rfdata is simulated using the method described in [1], with 8
% scatterers/mm^2
load rfdata;

% envelope detection
envelope = abs(hilbert(rfdata));

% obtain image size
realh = size(rfdata, 1);
realw = size(rfdata, 2);
image_range = size(rfdata);
line_number = image_range(2);
datalength = image_range(1);

xxx = (v * (realw/fs) * 10^3)/line_number : (v * (realw/fs) * 10^3)/line_number : (v * (realw/fs) * 10^2);% image width [cm]
yyy = (v * (realh/fs) * 10^3)/datalength : (v * (realh/fs) * 10^3)/datalength : (v * (realh/fs) * 10^2); % image height [cm]

% dynamic range
dBrange = [-40, 0];

% B-mode imaging
figure('visible','on');
imagesc(xxx, yyy, db(envelope/max(max(envelope))), dBrange);
colormap gray;
h = colorbar;
h.Title.String = 'dB';
h.TickLabelInterpreter = 'latex';
h.TickLabels = {'$-$40', '$-$30', '$-$20', '$-$10', '0'};
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
xtickformat('$%g$')
set(gca, 'fontsize', 20);
set(gca, 'fontname', 'Euclid');
set(gca, 'linewidth', 2);
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);
title('B-mode image');
axis image;

scan_step = ( xxx(end) - xxx(1) ) * 10^(-2) / realw; % distance between two adjacent scan lines [m]

% parallel computation of local NCRE estimates to obtain the NCRE map
% (ncre)
ncre = ent_map_block_cre(fs, v, scan_step, pulselength, 1, envelope);

% resize the NCRE map (ncre) to the original size of rf/envelope data to
% obtain ncre1
ncre1 = imresize(ncre, [realh realw], 'bicubic');

% NCRE parametric imaging
figure('visible','on');
imagesc(xxx, yyy, ncre1, [0.2 0.9]);
h = colorbar;
set(get(h, 'Title'), 'string', 'a.u.');
colormap(jet)
set(gca, 'TickLabelInterpreter', 'latex')
ytickformat('$%g$')
xtickformat('$%g$')
set(gca, 'fontsize', 20);
set(gca, 'fontname', 'Euclid');
set(gca, 'linewidth', 2);
xlabel('Lateral distance [cm]','fontsize',20);
ylabel('Axial distance [cm]','fontsize',20);
title('NCRE image');
axis image;