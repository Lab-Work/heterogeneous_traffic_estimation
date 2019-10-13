# heterogeneous_tse_IEEE_ITS_2019

Yanbing Wang
October, 2019

### Overview
This repository contains the source code for heterogeneous traffic estimation. This resulting journal paper is submitted to IEEE Transactions on Intelligent Transportation Systems (T-ITS). Suggested citation: 

Yanbing Wang and Daniel B. Work. "Enhanced estimation for heterogeneous traffic". *IEEE Transactions on Intelligent Transportation Systems*, 2019. (In review).

### Usage
#### Running variations of particle filters:
1. Standard PF: run `main_pf_creeping.m`, set `spatial_correlation` as "false"
2. PF+SCNM: run `main_pf_creeping.m`, set `spatial_correlation` as "true"
3. PAPF: run `main_pf_dual.m`, set `spatial_correlation` as "false"
4. PAPF+SCNM: run `main_pf_dual.m`, set `spatial_correlation` as "true"

#### Change traffic scenarios:
In `main_pf_creeping.m` or `main_pf_dual.m`
1. Set `test = 1` for freeflow traffic with overtaking
2. Set `test = 2` for congested traffic
3. Set `test = 3` for queue clearance
4. Set `test = 4` for creeping traffic

#### Toggle image display:
1. Set `show_sim = true` to turn on simulation plots, for example:
<img src="https://github.com/yanb514/heterogeneous_tse_IEEE_ITS_2019/blob/master/figures_ieee/test1_sim_snapshots.png" width="350" height="200" />
2. Set `show_est = true` to turn on estimation plots, for example:
<img src="https://github.com/yanb514/heterogeneous_tse_IEEE_ITS_2019/blob/master/figures_ieee/test1_filter_snapshots.png" width="350" height="200" />
